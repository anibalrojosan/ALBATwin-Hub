"""
ALBA aqueous chemistry: reference dissociation constants and van't Hoff scaling.

Reference K_a and K_w in mol/L are taken from Table SI.6.1 in
docs/supporting_informations/SI.6 (third column, header 293 K).

Equation SI.6.4 in the same file uses the van't Hoff form with T in Celsius and T_ref in the
text given as 298.15 K; the tabulated K values are at 293 K. This module uses T_REF_K = 293.15
so K_ref in ka_at_T matches the table temperature (same equation shape as SI.6.4).

The numeric Delta H values in default_dissociation_enthalpy_j_per_mol are placeholders (not in
the SI.6 excerpt in this repo). Replace when project has SI enthalpy data or calibration.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# SI.6.4: gas constant [J K^-1 mol^-1]
R_GAS = 8.314462618

# Temperature of K_A in Table SI.6.1 (293 K) [K]; used as T_ref in ka_at_T with tabulated K
T_REF_K = 293.15

# Molar masses for ALBA state totals [g mol^-1] (SI.6 uses 12, 14, 31)
M_C = 12.0
M_N = 14.0
M_P = 31.0


@dataclass(frozen=True, slots=True)
class AlbaDissociationConstantsRef:
    """Dissociation constants at T_REF_K (298.15 K), in mol/L. Field names match SI.6.1."""

    ka_co2: float
    ka_hco3: float
    ka_nh4: float
    ka_hno2: float
    ka_hno3: float
    ka_h3po4: float
    ka_h2po4: float
    ka_hpo4: float
    kw: float


@dataclass(frozen=True, slots=True)
class AlbaDissociationEnthalpy:
    """
    Standard reaction enthalpy Delta H in J/mol per equilibrium for SI.6.4.

    Same keys as AlbaDissociationConstantsRef; dh_kw is for water ionization.
    """

    dh_co2: float
    dh_hco3: float
    dh_nh4: float
    dh_hno2: float
    dh_hno3: float
    dh_h3po4: float
    dh_h2po4: float
    dh_hpo4: float
    dh_kw: float


def default_dissociation_constants_ref_molar() -> AlbaDissociationConstantsRef:
    """
    Build K_a,ref from docs/MATH_MODEL.md section 1.2.7 (pK_a at 298.15 K).

    Returns all values in mol/L. kw is 1e-14 at 25 C (common convention).
    """
    return AlbaDissociationConstantsRef(
        ka_co2=10.0**-6.37,
        ka_hco3=10.0**-10.33,
        ka_nh4=10.0**-9.25,
        ka_hno2=10.0**-3.35,
        ka_hno3=10.0**1.64,
        ka_h3po4=10.0**-2.14,
        ka_h2po4=10.0**-7.21,
        ka_hpo4=10.0**-12.67,
        kw=1.0e-14,
    )


def default_dissociation_enthalpy_j_per_mol() -> AlbaDissociationEnthalpy:
    """
    Placeholder Delta H in J/mol for van't Hoff (SI.6.4). Not sourced from repo SI files.

    Override when project adopts final values.
    """
    return AlbaDissociationEnthalpy(
        dh_co2=0,
        dh_hco3=0,
        dh_nh4=0,
        dh_hno2=0,
        dh_hno3=0,
        dh_h3po4=0,
        dh_h2po4=0,
        dh_hpo4=0,
        dh_kw=0,
    )


def ka_at_T(ka_ref_molar: float, delta_h_j_per_mol: float, t_celsius: float) -> float:
    """
    Scale equilibrium constant K from T_REF_K to t_celsius using SI.6.4.

    ka_ref_molar: K at T_REF_K in mol/L (for K_w use the same formula; units stay consistent with the reference).
    delta_h_j_per_mol: Delta H in J/mol for that reaction.
    t_celsius: temperature in degrees Celsius.

    Returns K at t_celsius in the same concentration units as ka_ref_molar.
    """
    t_k = t_celsius + 273.15
    return ka_ref_molar * math.exp(
        (delta_h_j_per_mol / R_GAS) * (1.0 / T_REF_K - 1.0 / t_k),
    )


def kw_at_T(kw_ref_molar: float, delta_h_kw_j_per_mol: float, t_celsius: float) -> float:
    """Water ionization constant at t_celsius; same math as ka_at_T. kw_ref_molar is often 1e-14 (mol/L)^2 at 25 C."""
    return ka_at_T(kw_ref_molar, delta_h_kw_j_per_mol, t_celsius)


def mol_c_per_m3_carbon_from_s_ic(s_ic_g_per_m3: float) -> float:
    """Moles of inorganic carbon per m3 from S_IC in g C/m3 (SI.6 row 7)."""
    return s_ic_g_per_m3 / M_C


def mol_n_per_m3_nitrogen_from_s_nh(s_nh_g_per_m3: float) -> float:
    """Moles of ammoniacal nitrogen per m3 from S_NH in g N/m3."""
    return s_nh_g_per_m3 / M_N


def mol_n_per_m3_nitrogen_from_s_no2(s_no2_g_per_m3: float) -> float:
    """Moles of nitrite nitrogen per m3 from S_NO2 in g N/m3."""
    return s_no2_g_per_m3 / M_N


def mol_n_per_m3_nitrogen_from_s_no3(s_no3_g_per_m3: float) -> float:
    """Moles of nitrate nitrogen per m3 from S_NO3 in g N/m3."""
    return s_no3_g_per_m3 / M_N


def mol_p_per_m3_phosphorus_from_s_po4(s_po4_g_per_m3: float) -> float:
    """Moles of inorganic phosphorus per m3 from S_PO4 in g P/m3."""
    return s_po4_g_per_m3 / M_P


def scale_dissociation_constants_at_t(
    ref: AlbaDissociationConstantsRef,
    dh: AlbaDissociationEnthalpy,
    t_celsius: float,
) -> AlbaDissociationConstantsRef:
    """Apply ka_at_T / kw_at_T to each field of ref using dh."""
    return AlbaDissociationConstantsRef(
        ka_co2=ka_at_T(ref.ka_co2, dh.dh_co2, t_celsius),
        ka_hco3=ka_at_T(ref.ka_hco3, dh.dh_hco3, t_celsius),
        ka_nh4=ka_at_T(ref.ka_nh4, dh.dh_nh4, t_celsius),
        ka_hno2=ka_at_T(ref.ka_hno2, dh.dh_hno2, t_celsius),
        ka_hno3=ka_at_T(ref.ka_hno3, dh.dh_hno3, t_celsius),
        ka_h3po4=ka_at_T(ref.ka_h3po4, dh.dh_h3po4, t_celsius),
        ka_h2po4=ka_at_T(ref.ka_h2po4, dh.dh_h2po4, t_celsius),
        ka_hpo4=ka_at_T(ref.ka_hpo4, dh.dh_hpo4, t_celsius),
        kw=kw_at_T(ref.kw, dh.dh_kw, t_celsius),
    )
