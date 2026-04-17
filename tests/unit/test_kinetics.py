"""Unit tests for ALBA process rates (``calculate_rates``, §4)."""

import numpy as np
import pytest
from pydantic import ValidationError

from bioprocess_twin.core.state import StateVector
from bioprocess_twin.models import EnvConditions, calculate_rates, default_alba
from bioprocess_twin.models.stoichiometry import N_PROCESSES, N_STATE


def _state_kwargs_17(**overrides: float) -> dict[str, float]:
    base = {
        "X_ALG": 100.0,
        "X_AOB": 10.0,
        "X_NOB": 5.0,
        "X_H": 50.0,
        "X_S": 20.0,
        "X_I": 15.0,
        "S_S": 30.0,
        "S_I": 10.0,
        "S_IC": 40.0,
        "S_ND": 2.0,
        "S_NH": 15.0,
        "S_NO2": 1.0,
        "S_NO3": 20.0,
        "S_N2": 0.0,
        "S_PO4": 5.0,
        "S_O2": 8.0,
        "S_H2O": 1000.0,
    }
    base.update(overrides)
    return base


def _env_plateau_algae() -> EnvConditions:
    """T_opt, pH_opt for algae; I_opt for Haldane (f_i_dim = 1)."""
    p = default_alba()
    return EnvConditions(
        temperature_C=p.cardinal_temp_alg.t_opt,
        pH=p.cardinal_ph_alg.ph_opt,
        irradiance_umol_m2_s=p.i_opt,
    )


def test_calculate_rates_shape_and_dtype() -> None:
    state = StateVector(**_state_kwargs_17())
    rho = calculate_rates(state, _env_plateau_algae())

    assert isinstance(rho, np.ndarray)
    assert rho.shape == (N_PROCESSES,)
    assert N_PROCESSES == 19
    assert rho.dtype == np.float64
    assert not np.any(np.isnan(rho))
    # rho[k] is process ρ_{k+1} (row k of Petersen); see calculate_rates docstring.
    assert len(rho) == 19


def test_calculate_rates_ndarray_17_and_18_slice() -> None:
    v = np.arange(1.0, 18.0, dtype=np.float64)
    env = EnvConditions(temperature_C=20.0, pH=7.0, irradiance_umol_m2_s=100.0)
    rho = calculate_rates(v, env)
    assert rho.shape == (N_PROCESSES,)

    v18 = np.concatenate([v, [99.0]])
    assert np.array_equal(calculate_rates(v18, env), rho)


def test_calculate_rates_rejects_wrong_length() -> None:
    env = EnvConditions(temperature_C=20.0, pH=7.0, irradiance_umol_m2_s=100.0)
    with pytest.raises(ValueError, match="must have length"):
        calculate_rates(np.zeros(16, dtype=np.float64), env)
    with pytest.raises(ValueError, match="must have length"):
        calculate_rates(np.zeros(20, dtype=np.float64), env)


def test_env_conditions_frozen() -> None:
    env = EnvConditions(temperature_C=20.0, pH=7.0, irradiance_umol_m2_s=100.0)
    with pytest.raises(ValidationError):
        env.temperature_C = 0.0  # type: ignore[misc]


def test_plateau_rho1_scales_with_mu_max_alg_times_X_alg() -> None:
    """Near-1 Monod, f_T=f_ph=f_i_dim=1; S_O2=0 ⇒ f_DO,g=1 (§3.4)."""
    p = default_alba()
    env = _env_plateau_algae()
    kwargs = _state_kwargs_17(
        S_IC=1e6,
        S_NH=1e6,
        S_PO4=1e6,
        S_O2=0.0,
    )
    state = StateVector(**kwargs)
    rho = calculate_rates(state, env, kinetic_parameters=p)

    expected_factor = p.mu_max_alg * kwargs["X_ALG"]
    # CTMI/CPM at cardinals can differ from 1.0 at ~1e-7 rel near t_opt/ph_opt.
    assert rho[0] == pytest.approx(expected_factor, rel=1e-6, abs=1e-6)


def test_zero_biomass_all_rates_zero() -> None:
    env = EnvConditions(temperature_C=20.0, pH=7.0, irradiance_umol_m2_s=200.0)
    kwargs = _state_kwargs_17(X_ALG=0.0, X_AOB=0.0, X_NOB=0.0, X_H=0.0)
    state = StateVector(**kwargs)
    rho = calculate_rates(state, env)
    assert np.all(rho == 0.0)


def test_rho2_nh4_block_goes_to_zero_when_s_nh_large() -> None:
    p = default_alba()
    env = _env_plateau_algae()
    kwargs = _state_kwargs_17(S_NH=1e9, S_IC=1e6, S_NO3=1e6, S_PO4=1e6, S_O2=0.0)
    state = StateVector(**kwargs)
    rho = calculate_rates(state, env, kinetic_parameters=p)

    assert rho[0] > 0.0
    assert abs(rho[1]) <= 1e-6 * max(abs(rho[0]), 1.0)


def test_rho8_anox_switch_zero_when_aerobic() -> None:
    p = default_alba()
    env = EnvConditions(
        temperature_C=p.cardinal_temp_h.t_opt,
        pH=p.cardinal_ph_h.ph_opt,
        irradiance_umol_m2_s=0.0,
    )
    # K_O/(K_O+S_O2) must be negligible vs aerobic ρ5–ρ7 pathways.
    kwargs = _state_kwargs_17(S_O2=1e6, S_S=50.0, S_NO2=5.0, S_PO4=5.0)
    state = StateVector(**kwargs)
    rho = calculate_rates(state, env, kinetic_parameters=p)
    assert rho[7] < 1e-4


def test_rho10_limited_by_no2_when_no2_starved() -> None:
    p = default_alba()
    env = EnvConditions(
        temperature_C=p.cardinal_temp_h.t_opt,
        pH=p.cardinal_ph_h.ph_opt,
        irradiance_umol_m2_s=0.0,
    )
    s_no2 = 1e-6
    s_no3 = 100.0
    kwargs = _state_kwargs_17(S_O2=0.01, S_NO2=s_no2, S_NO3=s_no3)
    state = StateVector(**kwargs)
    rho = calculate_rates(state, env, kinetic_parameters=p)

    m_no2 = s_no2 / (p.k_no2_h + s_no2)
    m_no3 = s_no3 / (p.k_no3_h + s_no3)
    assert m_no2 < m_no3
    anox = p.k_o_h / (p.k_o_h + 0.01)
    expected = p.b_max_r_h * p.eta_anox * 1.0 * 1.0 * anox * m_no2 * kwargs["X_H"]
    assert rho[9] == pytest.approx(expected, rel=1e-6)


def test_rho11_increases_when_x_s_doubles_at_fixed_x_h() -> None:
    p = default_alba()
    env = EnvConditions(
        temperature_C=p.cardinal_temp_h.t_opt,
        pH=p.cardinal_ph_h.ph_opt,
        irradiance_umol_m2_s=0.0,
    )
    base = _state_kwargs_17(X_S=10.0, X_H=50.0)
    r_lo = calculate_rates(StateVector(**base), env, kinetic_parameters=p)[10]
    r_hi = calculate_rates(StateVector(**{**base, "X_S": 20.0}), env, kinetic_parameters=p)[10]
    assert r_hi > r_lo


def test_ndarray_clips_negative_s_o2_for_hill() -> None:
    """Raw ndarray: negative S_O2 is clipped to 0 so f_DO Hill does not raise."""
    env = EnvConditions(temperature_C=20.0, pH=7.0, irradiance_umol_m2_s=100.0)
    v = np.array([1.0] * N_STATE, dtype=np.float64)
    v[15] = -1.0  # S_O2 column index in SI layout
    rho = calculate_rates(v, env)
    assert rho.shape == (N_PROCESSES,)
    assert not np.any(np.isnan(rho))


def test_rho1_zero_when_temperature_below_algal_cardinal_min() -> None:
    p = default_alba()
    env = EnvConditions(
        temperature_C=p.cardinal_temp_alg.t_min - 5.0,
        pH=p.cardinal_ph_alg.ph_opt,
        irradiance_umol_m2_s=p.i_opt,
    )
    kwargs = _state_kwargs_17(S_IC=1e6, S_NH=1e6, S_PO4=1e6, S_O2=0.0)
    rho = calculate_rates(StateVector(**kwargs), env, kinetic_parameters=p)
    assert rho[0] == pytest.approx(0.0, abs=1e-12)


def test_rho3_algal_respiration_hand_calc_half_monod() -> None:
    """ρ₃ at T_opt, pH_opt: only S_O2/(K_O+S_O2) varies; pick S_O2 = K_O for 0.5 Monod."""
    p = default_alba()
    env = _env_plateau_algae()
    s_o2 = p.k_o_alg
    kwargs = _state_kwargs_17(S_O2=s_o2)
    rho = calculate_rates(StateVector(**kwargs), env, kinetic_parameters=p)
    m_o2 = s_o2 / (p.k_o_alg + s_o2)
    expected = p.b_max_r_alg * m_o2 * kwargs["X_ALG"]
    assert rho[2] == pytest.approx(expected, rel=1e-6, abs=1e-9)


def test_rho14_aob_growth_hand_calc_ammonia_limited() -> None:
    p = default_alba()
    env = EnvConditions(
        temperature_C=p.cardinal_temp_aob.t_opt,
        pH=p.cardinal_ph_aob.ph_opt,
        irradiance_umol_m2_s=0.0,
    )
    s_nh = 0.05
    kwargs = _state_kwargs_17(
        S_NH=s_nh,
        S_O2=1e3,
        S_IC=1e3,
        S_PO4=1e3,
    )
    rho = calculate_rates(StateVector(**kwargs), env, kinetic_parameters=p)
    m_nh = s_nh / (p.k_n_aob + s_nh)
    expected = p.mu_max_aob * m_nh * kwargs["X_AOB"]
    assert rho[13] == pytest.approx(expected, rel=1e-6, abs=1e-9)


def test_mu_max_h_scales_rho5_linearly() -> None:
    p0 = default_alba()
    p1 = p0.model_copy(update={"mu_max_h": p0.mu_max_h * 2.0})
    env = EnvConditions(
        temperature_C=p0.cardinal_temp_h.t_opt,
        pH=p0.cardinal_ph_h.ph_opt,
        irradiance_umol_m2_s=0.0,
    )
    kwargs = _state_kwargs_17(S_S=1e3, S_O2=1e3, S_NH=1e3, S_PO4=1e3)
    state = StateVector(**kwargs)
    r0 = calculate_rates(state, env, kinetic_parameters=p0)[4]
    r1 = calculate_rates(state, env, kinetic_parameters=p1)[4]
    assert r1 == pytest.approx(2.0 * r0, rel=1e-9)
