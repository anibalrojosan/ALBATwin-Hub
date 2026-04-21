"""Unit tests for SI.7 gas–liquid transfer (Table SI.7.1)."""

from __future__ import annotations

import pytest

from bioprocess_twin.core.state import StateVector
from bioprocess_twin.models.chemistry import (
    default_dissociation_constants_ref_molar,
    default_dissociation_enthalpy_j_per_mol,
    scale_dissociation_constants_at_t,
)
from bioprocess_twin.models.gas_transfer import (
    T_REF_CELSIUS_KLA,
    calculate_gas_transfer,
    effective_kla_o2_per_day,
    henry_o2_co2_nh3,
    kinetic_kla_factor,
    liquid_volatile_co2_gC_per_m3,
    liquid_volatile_nh3_gN_per_m3,
)
from bioprocess_twin.models.kinetics import EnvConditions


def _minimal_state(*, s_o2: float, s_ic: float, s_nh: float) -> StateVector:
    return StateVector(
        X_ALG=0.0,
        X_AOB=0.0,
        X_NOB=0.0,
        X_H=0.0,
        X_S=0.0,
        X_I=0.0,
        S_S=0.0,
        S_I=0.0,
        S_IC=s_ic,
        S_ND=0.0,
        S_NH=s_nh,
        S_NO2=0.0,
        S_NO3=0.0,
        S_N2=0.0,
        S_PO4=0.0,
        S_O2=s_o2,
        S_H2O=0.0,
    )


def test_rho_o2_matches_table_si_7_1_closed_form() -> None:
    """rho_20 = theta^(T-20) * kLa * (H_O2*p_O2 - S_O2); no diffusivity ratio for oxygen."""
    t_c = 20.0
    theta = 1.024
    kla = 34.0
    p_o2 = 0.21
    s_o2 = 5.0
    henry = henry_o2_co2_nh3(t_c)
    expected = kinetic_kla_factor(theta, t_c) * kla * (henry.h_o2 * p_o2 - s_o2)

    st = _minimal_state(s_o2=s_o2, s_ic=30.0, s_nh=1.0)
    env = EnvConditions(temperature_C=t_c, pH=7.0, irradiance_umol_m2_s=0.0)
    rates = calculate_gas_transfer(st, env, h_plus_mol_m3=1.0e-5 * 1000.0, theta_kla=theta)
    assert rates.rho_o2_g_m3_d == pytest.approx(expected, rel=0.0, abs=1e-9)


def test_effective_kla_at_reference_temperature_is_identity() -> None:
    assert effective_kla_o2_per_day(34.0, 1.024, T_REF_CELSIUS_KLA) == pytest.approx(34.0)


def test_co2_driving_force_increases_with_lower_ph() -> None:
    """More unionized CO2 at low pH → larger volatile liquid term → smaller (H*p - S_vol) if S fixed."""
    t_c = 25.0
    k_at_t = scale_dissociation_constants_at_t(
        default_dissociation_constants_ref_molar(),
        default_dissociation_enthalpy_j_per_mol(),
        t_c,
    )
    s_ic = 50.0
    h_high = 1.0e-4 * 1000.0  # ~pH 4
    h_low = 1.0e-9 * 1000.0  # ~pH 8
    s_vol_acid = liquid_volatile_co2_gC_per_m3(s_ic, h_high, k_at_t.ka_co2, k_at_t.ka_hco3)
    s_vol_base = liquid_volatile_co2_gC_per_m3(s_ic, h_low, k_at_t.ka_co2, k_at_t.ka_hco3)
    assert s_vol_acid > s_vol_base
    henry = henry_o2_co2_nh3(t_c)
    p_co2 = 4.0e-4
    df_acid = henry.h_co2 * p_co2 - s_vol_acid
    df_base = henry.h_co2 * p_co2 - s_vol_base
    assert df_acid < df_base


def test_nh3_volatile_increases_with_ph() -> None:
    """Unionized NH3 rises with pH (lower [H+]) for fixed S_NH."""
    t_c = 25.0
    k_at_t = scale_dissociation_constants_at_t(
        default_dissociation_constants_ref_molar(),
        default_dissociation_enthalpy_j_per_mol(),
        t_c,
    )
    s_nh = 20.0
    h_acid = 1.0e-4 * 1000.0
    h_basic = 1.0e-10 * 1000.0
    v_acid = liquid_volatile_nh3_gN_per_m3(s_nh, h_acid, k_at_t.ka_nh4)
    v_basic = liquid_volatile_nh3_gN_per_m3(s_nh, h_basic, k_at_t.ka_nh4)
    assert v_basic > v_acid


def test_henry_o2_increases_as_temperature_decreases() -> None:
    """SI.7 narrative: lower T → higher solubility (larger H) for O2 and CO2."""
    h_cold = henry_o2_co2_nh3(10.0)
    h_warm = henry_o2_co2_nh3(30.0)
    assert h_cold.h_o2 > h_warm.h_o2
    assert h_cold.h_co2 > h_warm.h_co2


def test_henry_nh3_same_trend() -> None:
    """SI.7 narrative: lower T → higher solubility (larger H) for NH3."""
    h_cold = henry_o2_co2_nh3(5.0)
    h_warm = henry_o2_co2_nh3(35.0)
    assert h_cold.h_nh3 > h_warm.h_nh3


def test_calculate_gas_transfer_nh3_rate_depends_on_h_plus() -> None:
    """End-to-end: rho_22 changes with [H+] because the unionized NH3 pool (SI.6) does."""
    t_c = 25.0
    theta = 1.024
    st = _minimal_state(s_o2=8.0, s_ic=20.0, s_nh=15.0)
    env = EnvConditions(temperature_C=t_c, pH=7.0, irradiance_umol_m2_s=0.0)
    h_acidic = 1.0e-4 * 1000.0  # ~pH 4 in mol m^-3 convention
    h_basic = 1.0e-10 * 1000.0  # ~pH 10
    r_acid = calculate_gas_transfer(st, env, h_acidic, theta_kla=theta)
    r_base = calculate_gas_transfer(st, env, h_basic, theta_kla=theta)
    assert abs(r_acid.rho_nh3_gN_m3_d - r_base.rho_nh3_gN_m3_d) > 1.0
