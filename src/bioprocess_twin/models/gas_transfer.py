"""
SI.7 gas-liquid mass transfer (Table SI.7.1) for O2, CO2 (as carbon), and NH3 (as nitrogen).

Temperature-reference conventions (Casagli et al. SI; see docs/HYDROCHEMISTRY.md):

- Kinetic prefactor: theta_kla raised to (T_C - T_REF_CELSIUS_KLA) with T_REF_CELSIUS_KLA = 20 C.
  The tabulated oxygen k_L a is the value at 20 C; theta scales it to the operating T_C (C).
  This is not the same reference as aqueous equilibrium.

- Thermodynamic SSOT at 25 C (298.15 K): Henry correlations SI.7.3-SI.7.5 use 1/298.15 K^-1
  in the exponential anchor; dissociation constants for volatile fractions follow chemistry.T_REF_K
  and scale_dissociation_constants_at_t (K_a,ref at 298.15 K, van't Hoff to T_C).

Do not merge these references into a single T_ref symbol in calling code.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from bioprocess_twin.core.state import StateVector
from bioprocess_twin.models.chemistry import (
    M_C,
    M_N,
    default_dissociation_constants_ref_molar,
    default_dissociation_enthalpy_j_per_mol,
    mol_c_per_m3_carbon_from_s_ic,
    mol_n_per_m3_nitrogen_from_s_nh,
    scale_dissociation_constants_at_t,
    speciate_ammonia,
    speciate_carbonate,
)
from bioprocess_twin.models.kinetics import EnvConditions

# SI.7 / Table SI.7.1: k_L a for O2 is defined at this reference temperature [C].
T_REF_CELSIUS_KLA = 20.0

# MATH_MODEL.md section 1.2.6 [m^2 s^-1] (Perry, 2007)
D_O2_M2_S = 2.5e-9
D_CO2_M2_S = 2.1e-9
D_NH3_M2_S = 2.4e-9

# MATH_MODEL.md section 1.2.6 central values (defaults for GasTransferConditions)
DEFAULT_KLA_O2_D = 34.0
DEFAULT_P_O2_ATM = 0.21
DEFAULT_P_CO2_ATM = 4.0e-4
DEFAULT_P_NH3_ATM = 1.5e-6


@dataclass(frozen=True, slots=True)
class HenryConstantsGas:
    """
    Henry coefficients at t_celsius from SI.7.3-SI.7.5.

    The exponential uses 1/298.15 K^-1 as reference (25 C anchor in the SI), distinct
    from T_REF_CELSIUS_KLA used only for k_L a temperature scaling.
    """

    h_o2: float
    """[g O2 m^-3 atm^-1]"""

    h_co2: float
    """[g C as CO2 m^-3 atm^-1] (SI.7.4 includes 12/44)."""

    h_nh3: float
    """[g N as NH3 m^-3 atm^-1] (SI.7.5 includes 14/17)."""


def henry_o2_co2_nh3(t_celsius: float) -> HenryConstantsGas:
    """
    Evaluate SI.7.3-SI.7.5 at t_celsius [C].

    Reference inside the exponentials is 298.15 K (25 C), not T_REF_CELSIUS_KLA.
    """
    inv_t_k = 1.0 / (273.15 + t_celsius)
    inv_298 = 1.0 / 298.15
    d_inv = inv_t_k - inv_298
    h_o2 = 42.15 * math.exp(1700.0 * d_inv)
    h_co2 = (1511.13 * math.exp(2400.0 * d_inv)) * (12.0 / 44.0)
    h_nh3 = (4.63e5 * math.exp(2100.0 * d_inv)) * (14.0 / 17.0)
    return HenryConstantsGas(h_o2=h_o2, h_co2=h_co2, h_nh3=h_nh3)


def kinetic_kla_factor(theta_kla: float, t_celsius: float) -> float:
    """
    Return theta_kla raised to (t_celsius - T_REF_CELSIUS_KLA).

    SI.7 / Table SI.7.1 Arrhenius-style correction for k_L a relative to 20 C, independent
    of the 298.15 K reference used for Henry constants and K_a.
    """
    return float(theta_kla ** (t_celsius - T_REF_CELSIUS_KLA))


def effective_kla_o2_per_day(
    kla_o2_ref_per_day: float,
    theta_kla: float,
    t_celsius: float,
) -> float:
    """
    Effective k_L a for oxygen [d^-1] at t_celsius.

    kla_o2_ref_per_day is the value at 20 C (T_REF_CELSIUS_KLA); multiply by
    kinetic_kla_factor for other temperatures.
    """
    return kla_o2_ref_per_day * kinetic_kla_factor(theta_kla, t_celsius)


def diffusivity_ratio_sqrt_co2() -> float:
    """Return sqrt(D_CO2 / D_O2) for SI.7.1 / Table SI.7.1 rho_21."""
    return math.sqrt(D_CO2_M2_S / D_O2_M2_S)


def diffusivity_ratio_sqrt_nh3() -> float:
    """Return sqrt(D_NH3 / D_O2) for SI.7.1 / Table SI.7.1 rho_22."""
    return math.sqrt(D_NH3_M2_S / D_O2_M2_S)


def liquid_volatile_co2_gC_per_m3(
    s_ic_g_per_m3: float,
    h_plus_mol_per_m3: float,
    ka_co2_molar: float,
    ka_hco3_molar: float,
) -> float:
    """
    Dissolved volatile CO2 as g C m^-3 from SI.6.1 rows 8-9 (full carbonate speciation).

    Table SI.7.1 uses a shorthand in H_ion+; the full speciate_carbonate path matches
    HYDROCHEMISTRY.md section 7.4 and speciate_aqueous in chemistry.py.
    """
    c_tot_ic = mol_c_per_m3_carbon_from_s_ic(s_ic_g_per_m3)
    co2_mol, _, _ = speciate_carbonate(h_plus_mol_per_m3, ka_co2_molar, ka_hco3_molar, c_tot_ic)
    return co2_mol * M_C


def liquid_volatile_nh3_gN_per_m3(
    s_nh_g_per_m3: float,
    h_plus_mol_per_m3: float,
    ka_nh4_molar: float,
) -> float:
    """
    Unionized ammonia nitrogen as g N m^-3 from SI.6.1 row 2 (speciate_ammonia).

    Same K_a and [H+] unit convention as the rest of the chemistry module.
    """
    c_tot_nh = mol_n_per_m3_nitrogen_from_s_nh(s_nh_g_per_m3)
    nh3_mol, _ = speciate_ammonia(h_plus_mol_per_m3, ka_nh4_molar, c_tot_nh)
    return nh3_mol * M_N


@dataclass(frozen=True, slots=True)
class GasTransferConditions:
    """Atmosphere and oxygen k_L a reference (SI.7.2 / MATH_MODEL section 1.2.6)."""

    kla_o2_per_day: float
    p_o2_atm: float
    p_co2_atm: float
    p_nh3_atm: float

    @staticmethod
    def default_math_model() -> GasTransferConditions:
        """Central values from docs/MATH_MODEL.md section 1.2.6."""
        return GasTransferConditions(
            kla_o2_per_day=DEFAULT_KLA_O2_D,
            p_o2_atm=DEFAULT_P_O2_ATM,
            p_co2_atm=DEFAULT_P_CO2_ATM,
            p_nh3_atm=DEFAULT_P_NH3_ATM,
        )


@dataclass(frozen=True, slots=True)
class GasTransferRates:
    """Volumetric gas-liquid rates from Table SI.7.1 (positive sign = net dissolution)."""

    rho_o2_g_m3_d: float
    rho_co2_gC_m3_d: float
    rho_nh3_gN_m3_d: float


def calculate_gas_transfer(
    state: StateVector,
    env_conditions: EnvConditions,
    h_plus_mol_m3: float,
    *,
    theta_kla: float,
    gas_conditions: GasTransferConditions | None = None,
) -> GasTransferRates:
    """
    Compute rho_20, rho_21, rho_22 [g component m^-3 d^-1] per Table SI.7.1.

    Inputs: StateVector pools S_O2, S_IC, S_NH; env_conditions.temperature_C [C] for Henry,
    K_a(T), and theta^(T - 20); h_plus_mol_m3 must match solve_pH / charge balance (mol m^-3).

    Temperature: theta_kla scales k_L a from its reference at 20 C (T_REF_CELSIUS_KLA).
    Henry SI.7.3-7.5 and scale_dissociation_constants_at_t use the 298.15 K thermodynamic anchor
    described in chemistry.py and the module docstring above.
    """
    # Default atmosphere + kLa at 20 C reference if caller did not pass a bundle.
    gc = gas_conditions or GasTransferConditions.default_math_model()
    # Operating temperature [C] for Henry, van't Hoff K_a, and theta^(T-20).
    t_c = env_conditions.temperature_C
    # Dissociation constants at t_c (mol/L); needed for volatile CO2 and NH3 fractions.
    k_at_t = scale_dissociation_constants_at_t(
        default_dissociation_constants_ref_molar(),
        default_dissociation_enthalpy_j_per_mol(),
        t_c,
    )
    # Henry coefficients H_O2, H_CO2, H_NH3 at t_c (SI.7.3-7.5).
    henry = henry_o2_co2_nh3(t_c)
    # Oxygen kLa at t_c: kLa_ref * theta_kla ** (t_c - 20 C).
    kla_eff = effective_kla_o2_per_day(gc.kla_o2_per_day, theta_kla, t_c)

    # rho_20: linear driving force (g O2 m^-3 d^-1); no diffusivity ratio vs O2 in SI.7.1.
    delta_o2 = henry.h_o2 * gc.p_o2_atm - float(state.S_O2)
    rho_o2 = kla_eff * delta_o2

    # rho_21: volatile dissolved CO2 as carbon (g C m^-3) from full carbonate speciation.
    s_co2_vol = liquid_volatile_co2_gC_per_m3(float(state.S_IC), h_plus_mol_m3, k_at_t.ka_co2, k_at_t.ka_hco3)
    delta_co2 = henry.h_co2 * gc.p_co2_atm - s_co2_vol
    rho_co2 = kla_eff * diffusivity_ratio_sqrt_co2() * delta_co2

    # rho_22: unionized NH3 as nitrogen (g N m^-3); sqrt(D_NH3/D_O2) per SI.7.1.
    s_nh3_vol = liquid_volatile_nh3_gN_per_m3(float(state.S_NH), h_plus_mol_m3, k_at_t.ka_nh4)
    delta_nh3 = henry.h_nh3 * gc.p_nh3_atm - s_nh3_vol
    rho_nh3 = kla_eff * diffusivity_ratio_sqrt_nh3() * delta_nh3

    return GasTransferRates(
        rho_o2_g_m3_d=rho_o2,
        rho_co2_gC_m3_d=rho_co2,
        rho_nh3_gN_m3_d=rho_nh3,
    )
