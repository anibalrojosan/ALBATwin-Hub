import numpy as np
import pytest
from pydantic import ValidationError

from bioprocess_twin.core.state import StateVector, StateVectorVariant, state_array_len


def _state_kwargs_17(
    *,
    x_alg: float = 100.0,
    s_h2o: float = 1000.0,
) -> dict[str, float]:
    return {
        "X_ALG": x_alg,
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
        "S_H2O": s_h2o,
    }


def test_state_vector_instantiation():
    """Verify that StateVector initializes with 17 required + optional S_H_PROTON."""
    state = StateVector(**_state_kwargs_17())

    assert state.X_ALG == 100.0
    assert state.S_O2 == 8.0
    assert state.S_H_PROTON == 0.0
    assert len(StateVector.model_fields) == 18


def test_state_array_len_variants() -> None:
    assert state_array_len(StateVectorVariant.SI) == 17
    assert state_array_len(StateVectorVariant.OXYGEN_CLOSURE) == 17
    assert state_array_len(StateVectorVariant.OXYGEN_AND_PROTON_CLOSURE) == 18


def test_state_vector_to_array_default_si_17():
    """to_array() defaults to SI (17 components)."""
    state = StateVector(**_state_kwargs_17(x_alg=1.0, s_h2o=17.0))

    arr = state.to_array()

    assert isinstance(arr, np.ndarray)
    assert arr.shape == (17,)
    assert arr[0] == 1.0
    assert arr[16] == 17.0


def test_state_vector_to_array_oxygen_closure_17():
    state = StateVector(**_state_kwargs_17())
    arr = state.to_array(variant=StateVectorVariant.OXYGEN_CLOSURE)
    assert arr.shape == (17,)


def test_state_vector_to_array_proton_18():
    state = StateVector(**_state_kwargs_17(), S_H_PROTON=1.23e-6)
    arr = state.to_array(variant=StateVectorVariant.OXYGEN_AND_PROTON_CLOSURE)
    assert arr.shape == (18,)
    assert arr[17] == 1.23e-6


def test_state_vector_from_array_17_resets_proton():
    arr = np.arange(1.0, 18.0, dtype=float)

    state = StateVector.from_array(arr[:17])

    assert state.S_H2O == 17.0
    assert state.S_H_PROTON == 0.0


def test_state_vector_from_array_18():
    arr = np.concatenate([np.arange(1.0, 18.0, dtype=float), [0.42]])

    state = StateVector.from_array(arr)

    assert state.S_H2O == 17.0
    assert state.S_H_PROTON == 0.42


def test_state_vector_from_array_explicit_variant_mismatch():
    arr = np.arange(1.0, 18.0, dtype=float)
    with pytest.raises(ValueError, match="Expected array length 18"):
        StateVector.from_array(arr, variant=StateVectorVariant.OXYGEN_AND_PROTON_CLOSURE)


def test_state_vector_from_array_infer_rejects_other_lengths():
    with pytest.raises(ValueError, match="Expected array length 17 or 18"):
        StateVector.from_array(np.array([1.0, 2.0, 3.0]))


def test_state_vector_physical_validation():
    """Verify that negative concentrations raise a ValidationError."""
    with pytest.raises(ValidationError):
        StateVector(**{**_state_kwargs_17(), "X_ALG": -10.0})


def test_state_vector_from_array_validation():
    """Verify that from_array() also enforces non-negative constraints."""
    arr = np.array([-1.0] + [2.0] * 16)
    with pytest.raises(ValidationError):
        StateVector.from_array(arr)


def test_state_vector_immutability():
    """Ensure the state object is frozen."""
    state = StateVector(**_state_kwargs_17(x_alg=1.0, s_h2o=17.0))

    with pytest.raises(ValidationError):
        state.X_ALG = 100.0
