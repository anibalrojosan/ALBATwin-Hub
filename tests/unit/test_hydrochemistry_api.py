"""
Stage 5: tests for hydrochemistry_api (StateVector + EnvConditions entry points).

Includes day vs night EnvConditions snapshots (temperature and irradiance) with
kinetics after aligning EnvConditions.pH to SI.6 solved pH. This exercises the
pattern for time-varying forcing and dynamic pH before Sprint 4 Orchestration and
Numerical Solver (ODE assembly, LSODA, transport). Here all checks are instantaneous
at two fixed env snapshots, not a time integration.
"""

from __future__ import annotations

import math

from bioprocess_twin.core.state import StateVector
from bioprocess_twin.models.chemistry import (
    ChargeBalanceInputs,
    SpeciationTotals,
    solve_pH,
)
from bioprocess_twin.models.gas_transfer import calculate_gas_transfer
from bioprocess_twin.models.hydrochemistry_api import (
    hydrochemistry_step,
    solve_pH_from_state,
    speciation_totals_from_state,
)
from bioprocess_twin.models.kinetic_parameters import default_alba
from bioprocess_twin.models.kinetics import EnvConditions, calculate_rates


def _state_vector_stage3_internal_regression() -> StateVector:
    """State totals matching test_solve_ph_internal_regression_case (mol pools via g m^-3)."""
    return StateVector(
        X_ALG=0.0,
        X_AOB=0.0,
        X_NOB=0.0,
        X_H=0.0,
        X_S=0.0,
        X_I=0.0,
        S_S=0.0,
        S_I=0.0,
        S_IC=0.025 * 12.0,
        S_ND=0.0,
        S_NH=0.008 * 14.0,
        S_NO2=0.0005 * 14.0,
        S_NO3=0.0025 * 14.0,
        S_N2=0.0,
        S_PO4=0.0008 * 31.0,
        S_O2=8.0,
        S_H2O=0.0,
    )


def test_speciation_totals_from_state_matches_stage3_mol_totals() -> None:
    """Totals built from StateVector match explicit SpeciationTotals used in stage 3 tests."""
    st = _state_vector_stage3_internal_regression()
    got = speciation_totals_from_state(st)
    assert math.isclose(got.c_tot_nh3, 0.008, rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(got.c_tot_no2, 0.0005, rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(got.c_tot_no3, 0.0025, rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(got.c_tot_ic, 0.025, rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(got.c_tot_po4, 0.0008, rel_tol=0.0, abs_tol=1e-15)


def test_solve_pH_from_state_matches_direct_charge_balance_inputs() -> None:
    """solve_pH_from_state agrees with solve_pH on equivalent ChargeBalanceInputs (stage 3 parity)."""
    st = _state_vector_stage3_internal_regression()
    env = EnvConditions(temperature_C=25.0, pH=6.5, irradiance_umol_m2_s=100.0)
    direct = ChargeBalanceInputs(
        totals=SpeciationTotals(
            c_tot_nh3=0.008,
            c_tot_no2=0.0005,
            c_tot_no3=0.0025,
            c_tot_ic=0.025,
            c_tot_po4=0.0008,
        ),
        t_celsius=25.0,
        delta_cat_an_mol_per_m3=0.006,
    )
    r_direct = solve_pH(direct, initial_ph=7.2)
    r_api = solve_pH_from_state(st, env, delta_cat_an_mol_per_m3=0.006, initial_ph=7.2)
    assert r_api.converged == r_direct.converged
    assert math.isclose(r_api.ph, r_direct.ph, rel_tol=0.0, abs_tol=1e-10)
    assert math.isclose(r_api.h_plus_mol_per_m3, r_direct.h_plus_mol_per_m3, rel_tol=0.0, abs_tol=1e-10)


def test_hydrochemistry_step_gas_matches_manual_chain() -> None:
    """hydrochemistry_step gas_rates match calculate_gas_transfer after solve_pH."""
    st = _state_vector_stage3_internal_regression()
    env = EnvConditions(temperature_C=25.0, pH=7.0, irradiance_umol_m2_s=0.0)
    theta = default_alba().theta_kla
    step = hydrochemistry_step(st, env, theta_kla=theta, delta_cat_an_mol_per_m3=0.006, initial_ph=7.2)
    manual = calculate_gas_transfer(st, env, step.ph_result.h_plus_mol_per_m3, theta_kla=theta)
    assert math.isclose(step.gas_rates.rho_o2_g_m3_d, manual.rho_o2_g_m3_d, rel_tol=0.0, abs_tol=1e-9)
    assert math.isclose(step.gas_rates.rho_co2_gC_m3_d, manual.rho_co2_gC_m3_d, rel_tol=0.0, abs_tol=1e-9)
    assert math.isclose(step.gas_rates.rho_nh3_gN_m3_d, manual.rho_nh3_gN_m3_d, rel_tol=0.0, abs_tol=1e-9)


def test_diel_env_snapshots_change_kinetics_when_ph_aligned() -> None:
    """
    Day vs night: same StateVector, same temperature for identical SI.6 pH, but
    irradiance high vs zero. After replacing EnvConditions.pH with solved pH, rho1
    (algal growth on NH4) should be higher in daytime. Prepares Sprint 4 where env
    and pH will vary along a time axis inside the ODE RHS.
    """
    st = StateVector(
        X_ALG=80.0,
        X_AOB=0.0,
        X_NOB=0.0,
        X_H=0.0,
        X_S=0.0,
        X_I=0.0,
        S_S=0.0,
        S_I=0.0,
        S_IC=2.0,
        S_ND=0.0,
        S_NH=0.15,
        S_NO2=0.01,
        S_NO3=0.05,
        S_N2=0.0,
        S_PO4=0.08,
        S_O2=10.0,
        S_H2O=0.0,
    )
    t_c = 25.0
    env_for_ph = EnvConditions(temperature_C=t_c, pH=7.0, irradiance_umol_m2_s=400.0)

    ph_res = solve_pH_from_state(st, env_for_ph, delta_cat_an_mol_per_m3=0.0, initial_ph=7.0)
    assert ph_res.converged

    params = default_alba()
    env_day = EnvConditions(
        temperature_C=t_c,
        pH=ph_res.ph,
        irradiance_umol_m2_s=400.0,
    )
    env_night = EnvConditions(
        temperature_C=t_c,
        pH=ph_res.ph,
        irradiance_umol_m2_s=0.0,
    )
    rho_day = calculate_rates(st, env_day, params)
    rho_night = calculate_rates(st, env_night, params)
    assert rho_day[0] > rho_night[0]
    assert rho_night[0] >= 0.0
