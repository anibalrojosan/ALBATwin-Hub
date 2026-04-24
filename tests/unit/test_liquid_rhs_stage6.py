"""Stage 6 tests for 22-process liquid RHS assembly."""

from __future__ import annotations

import math

import numpy as np

from bioprocess_twin.core.state import StateVector
from bioprocess_twin.models.kinetic_parameters import default_alba
from bioprocess_twin.models.kinetics import EnvConditions, calculate_rates
from bioprocess_twin.models.stoichiometry import (
    N_PROCESSES,
    N_PROCESSES_WITH_GAS_TRANSFER,
    N_STATE,
    S_IC,
    S_NH,
    S_O2,
    get_gas_transfer_matrix,
    get_petersen_matrix,
    get_petersen_matrix_with_gas_transfer,
)
from bioprocess_twin.simulator.liquid_rhs import evaluate_liquid_rhs


def _stage6_state() -> StateVector:
    """Representative non-zero liquid state for Stage 6 checks."""
    return StateVector(
        X_ALG=80.0,
        X_AOB=25.0,
        X_NOB=22.0,
        X_H=120.0,
        X_S=50.0,
        X_I=35.0,
        S_S=40.0,
        S_I=20.0,
        S_IC=35.0,
        S_ND=5.0,
        S_NH=12.0,
        S_NO2=1.2,
        S_NO3=5.5,
        S_N2=8.0,
        S_PO4=4.0,
        S_O2=7.0,
        S_H2O=0.0,
    )


def test_get_petersen_matrix_with_gas_transfer_shape_and_mapping() -> None:
    """22-row matrix appends SI.3.1 rows 20-22 with one-hot O2/IC/NH mapping."""
    s_bio = get_petersen_matrix()
    s_gas = get_gas_transfer_matrix()
    s_full = get_petersen_matrix_with_gas_transfer()

    assert s_bio.shape == (N_PROCESSES, N_STATE)
    assert s_gas.shape == (3, N_STATE)
    assert s_full.shape == (N_PROCESSES_WITH_GAS_TRANSFER, N_STATE)

    assert np.allclose(s_full[:N_PROCESSES, :], s_bio)
    assert np.allclose(s_full[N_PROCESSES:, :], s_gas)

    expected_20 = np.zeros(N_STATE)
    expected_20[S_O2] = 1.0
    expected_21 = np.zeros(N_STATE)
    expected_21[S_IC] = 1.0
    expected_22 = np.zeros(N_STATE)
    expected_22[S_NH] = 1.0

    assert np.allclose(s_gas[0], expected_20)
    assert np.allclose(s_gas[1], expected_21)
    assert np.allclose(s_gas[2], expected_22)


def test_liquid_rhs_equals_bio_plus_sparse_gas_sum() -> None:
    """RHS equals bio block plus sparse gas block in SI columns."""
    st = _stage6_state()
    env = EnvConditions(temperature_C=25.0, pH=8.8, irradiance_umol_m2_s=350.0)
    out = evaluate_liquid_rhs(st, env)

    s_bio = get_petersen_matrix()
    s_gas = get_gas_transfer_matrix()
    s_full = get_petersen_matrix_with_gas_transfer()

    rhs_bio = s_bio.T @ out.rho_bio_g_m3_d
    rhs_gas = s_gas.T @ out.rho_gas_g_m3_d
    rhs_sum = rhs_bio + rhs_gas
    rhs_full = s_full.T @ out.rho_full_g_m3_d

    assert np.allclose(rhs_sum, out.dcdt_g_m3_d)
    assert np.allclose(rhs_full, out.dcdt_g_m3_d)

    gas_nonzero_idx = np.where(np.abs(rhs_gas) > 0.0)[0].tolist()
    assert gas_nonzero_idx == [S_IC, S_NH, S_O2]
    assert math.isclose(rhs_gas[S_O2], out.rho_gas_g_m3_d[0], rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(rhs_gas[S_IC], out.rho_gas_g_m3_d[1], rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(rhs_gas[S_NH], out.rho_gas_g_m3_d[2], rel_tol=0.0, abs_tol=1e-12)


def test_liquid_rhs_aligns_kinetics_with_solved_ph() -> None:
    """Orchestrator uses solved pH for biological rates in Stage 6 path."""
    st = _stage6_state()
    params = default_alba()
    env = EnvConditions(temperature_C=22.0, pH=3.5, irradiance_umol_m2_s=150.0)

    out = evaluate_liquid_rhs(st, env, kinetic_parameters=params)

    env_aligned = EnvConditions(
        temperature_C=env.temperature_C,
        pH=out.ph_result.ph,
        irradiance_umol_m2_s=env.irradiance_umol_m2_s,
    )
    rho_aligned = calculate_rates(st, env_aligned, params)
    rho_input_ph = calculate_rates(st, env, params)

    assert np.allclose(out.rho_bio_g_m3_d, rho_aligned)
    assert not np.allclose(rho_aligned, rho_input_ph)
    assert not math.isclose(out.ph_result.ph, env.pH, rel_tol=0.0, abs_tol=1e-6)


def test_liquid_rhs_o2_gas_flux_changes_sign_with_o2_state() -> None:
    """Gas O2 contribution flips sign between undersaturated and supersaturated states."""
    env = EnvConditions(temperature_C=25.0, pH=7.5, irradiance_umol_m2_s=0.0)

    st_low = _stage6_state().model_copy(update={"S_O2": 2.0})
    st_high = _stage6_state().model_copy(update={"S_O2": 20.0})
    out_low = evaluate_liquid_rhs(st_low, env)
    out_high = evaluate_liquid_rhs(st_high, env)

    rhs_gas_low = get_gas_transfer_matrix().T @ out_low.rho_gas_g_m3_d
    rhs_gas_high = get_gas_transfer_matrix().T @ out_high.rho_gas_g_m3_d

    assert out_low.gas_rates.rho_o2_g_m3_d > 0.0
    assert out_high.gas_rates.rho_o2_g_m3_d < 0.0
    assert rhs_gas_low[S_O2] > 0.0
    assert rhs_gas_high[S_O2] < 0.0


def test_liquid_rhs_stage6_regression_snapshot() -> None:
    """Pinned Stage 6 snapshot for pH and rho vector lengths."""
    st = _stage6_state()
    env = EnvConditions(temperature_C=25.0, pH=7.0, irradiance_umol_m2_s=300.0)
    out = evaluate_liquid_rhs(st, env)

    assert out.rho_bio_g_m3_d.shape == (N_PROCESSES,)
    assert out.rho_gas_g_m3_d.shape == (3,)
    assert out.rho_full_g_m3_d.shape == (N_PROCESSES_WITH_GAS_TRANSFER,)
    assert out.dcdt_g_m3_d.shape == (N_STATE,)
    assert math.isclose(out.ph_result.ph, 4.3158655499793985, rel_tol=0.0, abs_tol=1e-12)
