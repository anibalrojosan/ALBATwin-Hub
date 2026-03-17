import numpy as np
import pytest
from pydantic import ValidationError

from bioprocess_twin.core.state import StateVector


def test_state_vector_instantiation():
    """Verify that StateVector initializes exactly 17 variables correctly."""
    state = StateVector(
        X_ALG=100.0,
        X_AOB=10.0,
        X_NOB=5.0,
        X_H=50.0,
        X_S=20.0,
        X_I=15.0,
        S_S=30.0,
        S_I=10.0,
        S_IC=40.0,
        S_ND=2.0,
        S_NH=15.0,
        S_NO2=1.0,
        S_NO3=20.0,
        S_N2=0.0,
        S_PO4=5.0,
        S_O2=8.0,
        S_H2O=1000.0,
    )

    assert state.X_ALG == 100.0
    assert state.S_O2 == 8.0
    # Also check if exactly 17 fields are defined in the model
    assert len(StateVector.model_fields) == 17


def test_state_vector_to_array():
    """Ensure to_array() returns a numpy array with the correct order."""
    state = StateVector(
        X_ALG=1.0,
        X_AOB=2.0,
        X_NOB=3.0,
        X_H=4.0,
        X_S=5.0,
        X_I=6.0,
        S_S=7.0,
        S_I=8.0,
        S_IC=9.0,
        S_ND=10.0,
        S_NH=11.0,
        S_NO2=12.0,
        S_NO3=13.0,
        S_N2=14.0,
        S_PO4=15.0,
        S_O2=16.0,
        S_H2O=17.0,
    )

    arr = state.to_array()

    assert isinstance(arr, np.ndarray)
    assert arr.shape == (17,)
    assert arr[0] == 1.0  # X_ALG
    assert arr[10] == 11.0  # S_NH
    assert arr[15] == 16.0  # S_O2
    assert arr[16] == 17.0  # S_H2O


def test_state_vector_from_array():
    """Ensure from_array() correctly populates the named attributes."""
    arr = np.array(
        [
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
            8.0,
            9.0,
            10.0,
            11.0,
            12.0,
            13.0,
            14.0,
            15.0,
            16.0,
            17.0,
        ]
    )

    state = StateVector.from_array(arr)

    assert state.X_ALG == 1.0
    assert state.S_NH == 11.0
    assert state.S_O2 == 16.0
    assert state.S_H2O == 17.0


def test_state_vector_physical_validation():
    """Verify that negative concentrations raise a ValidationError."""
    with pytest.raises(ValidationError):
        StateVector(
            X_ALG=-10.0,  # Invalid
            X_AOB=10.0,
            X_NOB=5.0,
            X_H=50.0,
            X_S=20.0,
            X_I=15.0,
            S_S=30.0,
            S_I=10.0,
            S_IC=40.0,
            S_ND=2.0,
            S_NH=15.0,
            S_NO2=1.0,
            S_NO3=20.0,
            S_N2=0.0,
            S_PO4=5.0,
            S_O2=8.0,
            S_H2O=1000.0,
        )


def test_state_vector_from_array_validation():
    """Verify that from_array() also enforces non-negative constraints."""
    arr = np.array(
        [
            -1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
            8.0,
            9.0,
            10.0,
            11.0,
            12.0,
            13.0,
            14.0,
            15.0,
            16.0,
            17.0,
        ]
    )
    with pytest.raises(ValidationError):
        StateVector.from_array(arr)


def test_state_vector_immutability():
    """Ensure the state object is frozen."""
    state = StateVector(
        X_ALG=1.0,
        X_AOB=2.0,
        X_NOB=3.0,
        X_H=4.0,
        X_S=5.0,
        X_I=6.0,
        S_S=7.0,
        S_I=8.0,
        S_IC=9.0,
        S_ND=10.0,
        S_NH=11.0,
        S_NO2=12.0,
        S_NO3=13.0,
        S_N2=14.0,
        S_PO4=15.0,
        S_O2=16.0,
        S_H2O=17.0,
    )

    with pytest.raises(ValidationError):
        state.X_ALG = 100.0
