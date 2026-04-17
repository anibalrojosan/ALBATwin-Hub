import pytest

from bioprocess_twin.models.kinetic_modifiers import (
    f_do_decay_hill,
    f_do_growth_hill,
    f_i_haldane,
    f_ph_cpm,
    f_t_decay_arrhenius,
    f_t_growth_ctmi,
)
from bioprocess_twin.models.kinetic_parameters import default_alba


def test_f_t_growth_ctmi_is_one_at_t_opt() -> None:
    card = default_alba().cardinal_temp_alg
    assert f_t_growth_ctmi(card.t_opt, card) == pytest.approx(1.0, rel=0, abs=1e-9)


def test_f_t_growth_ctmi_zero_outside_window() -> None:
    card = default_alba().cardinal_temp_alg
    assert f_t_growth_ctmi(card.t_min - 5.0, card) == 0.0
    assert f_t_growth_ctmi(card.t_max + 5.0, card) == 0.0


def test_f_t_growth_ctmi_zero_at_endpoints() -> None:
    card = default_alba().cardinal_temp_h
    assert f_t_growth_ctmi(card.t_min, card) == pytest.approx(0.0, abs=1e-12)
    assert f_t_growth_ctmi(card.t_max, card) == pytest.approx(0.0, abs=1e-12)


def test_f_t_decay_arrhenius_reference_temperature() -> None:
    assert f_t_decay_arrhenius(20.0, 1.07) == pytest.approx(1.0)
    assert f_t_decay_arrhenius(25.0, 1.07) == pytest.approx(1.07**5.0)


def test_f_t_decay_arrhenius_rejects_nonpositive_theta() -> None:
    with pytest.raises(ValueError, match="theta"):
        f_t_decay_arrhenius(20.0, 0.0)


def test_f_ph_cpm_is_one_at_ph_opt() -> None:
    card = default_alba().cardinal_ph_alg
    assert f_ph_cpm(card.ph_opt, card) == pytest.approx(1.0, rel=0, abs=1e-9)


def test_f_ph_cpm_zero_outside_range() -> None:
    card = default_alba().cardinal_ph_aob
    assert f_ph_cpm(card.ph_min - 0.5, card) == 0.0
    assert f_ph_cpm(card.ph_max + 0.5, card) == 0.0


def test_f_i_haldane_at_i_opt_equals_mu_max() -> None:
    p = default_alba()
    mu = 2.5
    f = f_i_haldane(p.i_opt, mu, p.alpha_light, p.i_opt)
    assert f == pytest.approx(mu)


def test_f_i_haldane_requires_positive_parameters() -> None:
    with pytest.raises(ValueError):
        f_i_haldane(100.0, 1.0, 0.0, 300.0)


def test_f_do_growth_hill_limits() -> None:
    """High DO inhibits growth (§3.4); at S_O2=0 the Hill factor is 1 (no high-O2 penalty)."""
    p = default_alba()
    assert f_do_growth_hill(0.0, p.ec50_o2, p.hill_n_o2) == pytest.approx(1.0, abs=1e-12)
    assert f_do_growth_hill(1e6, p.ec50_o2, p.hill_n_o2) == pytest.approx(0.0, abs=1e-6)


def test_f_do_decay_hill_limits() -> None:
    p = default_alba()
    assert f_do_decay_hill(0.0, p.ec50_o2, p.hill_n_o2) == pytest.approx(0.0, abs=1e-12)
    assert f_do_decay_hill(1e6, p.ec50_o2, p.hill_n_o2) == pytest.approx(1.0, rel=1e-6)


def test_f_do_rejects_invalid_args() -> None:
    with pytest.raises(ValueError):
        f_do_growth_hill(-1.0, 20.0, 15.0)
