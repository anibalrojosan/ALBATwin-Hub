import numpy as np
import pytest
from pydantic import ValidationError

from bioprocess_twin.core.state import StateVector
from bioprocess_twin.models import EnvConditions, calculate_rates
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


def _env() -> EnvConditions:
    return EnvConditions(temperature_C=20.0, pH=7.5, irradiance_umol_m2_s=200.0)


def test_calculate_rates_state_vector_returns_zeros_stub() -> None:
    state = StateVector(**_state_kwargs_17())
    rho = calculate_rates(state, _env())

    assert isinstance(rho, np.ndarray)
    assert rho.shape == (N_PROCESSES,)
    assert rho.dtype == np.float64
    assert np.all(rho == 0.0)
    assert not np.any(np.isnan(rho))


def test_calculate_rates_ndarray_17() -> None:
    v = np.arange(1.0, 18.0, dtype=np.float64)
    rho = calculate_rates(v, _env())

    assert rho.shape == (N_PROCESSES,)
    assert rho.dtype == np.float64
    assert np.all(rho == 0.0)


def test_calculate_rates_ndarray_18_uses_first_17_for_stub() -> None:
    v18 = np.concatenate([np.arange(1.0, 18.0, dtype=np.float64), [99.0]])
    v17 = v18[:N_STATE]

    r18 = calculate_rates(v18, _env())
    r17 = calculate_rates(v17, _env())

    assert np.array_equal(r18, r17)


def test_calculate_rates_rejects_wrong_length() -> None:
    env = _env()
    with pytest.raises(ValueError, match="must have length"):
        calculate_rates(np.zeros(16, dtype=np.float64), env)
    with pytest.raises(ValueError, match="must have length"):
        calculate_rates(np.zeros(20, dtype=np.float64), env)


def test_env_conditions_frozen() -> None:
    env = _env()
    with pytest.raises(ValidationError):
        env.temperature_C = 0.0  # type: ignore[misc]
