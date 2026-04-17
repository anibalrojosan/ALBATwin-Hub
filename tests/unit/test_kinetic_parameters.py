import pytest
from pydantic import ValidationError

from bioprocess_twin.models.kinetic_parameters import (
    CardinalPH,
    CardinalTemperature,
    KineticParameters,
    default_alba,
)


def test_default_alba_matches_math_model_nominals() -> None:
    p = default_alba()
    assert p.mu_max_alg == 2.5
    assert p.k_n_alg == 0.1
    assert p.k_a == 1.0
    assert p.k_c_nob == 0.5
    assert p.ec50_o2 == 20.0
    assert p.hill_n_o2 == 15.0
    assert p.i_opt == 300.0
    assert p.alpha_light == 0.01
    assert p.theta_amm == 1.12
    assert p.cardinal_temp_alg.t_min == -10.0
    assert p.cardinal_temp_alg.t_opt == 20.0
    assert p.cardinal_temp_alg.t_max == 42.0
    assert p.cardinal_ph_h.ph_opt == 7.0


def test_kinetic_parameters_frozen() -> None:
    p = default_alba()
    with pytest.raises(ValidationError):
        p.mu_max_alg = 1.0  # type: ignore[misc]


def test_cardinal_temperature_frozen() -> None:
    c = CardinalTemperature(t_min=0.0, t_opt=20.0, t_max=40.0)
    with pytest.raises(ValidationError):
        c.t_opt = 25.0  # type: ignore[misc]


def test_default_alba_is_kinetic_parameters_instance() -> None:
    assert isinstance(default_alba(), KineticParameters)


def test_cardinal_ph_bounds() -> None:
    ph = CardinalPH(ph_min=2.0, ph_opt=8.4, ph_max=12.0)
    assert ph.ph_min < ph.ph_opt < ph.ph_max
