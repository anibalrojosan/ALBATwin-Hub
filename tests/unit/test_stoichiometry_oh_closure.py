"""Tests for experimental O/H Petersen closure layer (SI matrix unchanged)."""

from __future__ import annotations

import numpy as np
import pytest

from bioprocess_twin.models.stoichiometry import S_H2O, get_composition_matrix, get_petersen_matrix
from bioprocess_twin.models.stoichiometry_closure import (
    build_petersen_matrix_with_oh_closure,
    compute_mass_balance_matrix,
    compute_stoichiometric_s_h2o_total_for_row,
    get_petersen_matrix_for_simulation,
    list_oh_mass_balance_violations,
)
from bioprocess_twin.stoichiometry_validation import MASS_BALANCE_ATOL


def test_list_oh_mass_balance_violations_matches_si_failures() -> None:
    S_si = get_petersen_matrix()
    comp = get_composition_matrix()
    B = compute_mass_balance_matrix(S_si, comp)
    viol = list_oh_mass_balance_violations(B)
    for i, b_o, b_h in viol:
        assert abs(b_o) > MASS_BALANCE_ATOL or abs(b_h) > MASS_BALANCE_ATOL
        assert abs(B[i, 1] - b_o) < 1e-15
        assert abs(B[i, 5] - b_h) < 1e-15
    for i in range(19):
        if i in {v[0] for v in viol}:
            continue
        assert abs(B[i, 1]) <= MASS_BALANCE_ATOL
        assert abs(B[i, 5]) <= MASS_BALANCE_ATOL


def test_inventory_non_empty_until_full_oh_closure() -> None:
    viol = list_oh_mass_balance_violations()
    assert len(viol) > 0


def test_closure_preserves_non_h2o_columns() -> None:
    S_si = get_petersen_matrix()
    S_cl, _, detail = build_petersen_matrix_with_oh_closure()
    mask = np.ones(S_si.shape[1], dtype=bool)
    mask[S_H2O] = False
    assert np.allclose(S_si[:, mask], S_cl[:, mask])


def test_closure_cod_c_n_p_balances_unchanged() -> None:
    S_si = get_petersen_matrix()
    comp = get_composition_matrix()
    B_si = compute_mass_balance_matrix(S_si, comp)
    _, B_cl, _ = build_petersen_matrix_with_oh_closure()
    for k in (0, 2, 3, 4):
        np.testing.assert_allclose(B_cl[:, k], B_si[:, k], rtol=0, atol=1e-12)


def test_closure_oxygen_near_zero_for_adjusted_rows() -> None:
    _, B_cl, detail = build_petersen_matrix_with_oh_closure()
    for i in detail.rows_adjusted:
        assert abs(B_cl[i, 1]) < 1e-6


def test_closure_max_oxygen_residual_small() -> None:
    _, B_cl, _ = build_petersen_matrix_with_oh_closure()
    assert np.max(np.abs(B_cl[:, 1])) < 2e-3


def test_stoichiometric_s_h2o_total_matches_delta_from_si_oxygen_balance() -> None:
    """alpha_i = -L_i/w_o equals s0 - B_O^SI/w_o for adjusted rows."""
    S_si = get_petersen_matrix()
    comp = get_composition_matrix()
    B_si = compute_mass_balance_matrix(S_si, comp)
    w_o = float(comp[1, S_H2O])
    for i in range(S_si.shape[0]):
        b_o = float(B_si[i, 1])
        s0 = float(S_si[i, S_H2O])
        alpha = compute_stoichiometric_s_h2o_total_for_row(i, S_si, comp)
        np.testing.assert_allclose(alpha, s0 - b_o / w_o, rtol=0, atol=1e-12)


def test_stoichiometric_water_and_oxygen_exact_identical() -> None:
    S_a, _, _ = build_petersen_matrix_with_oh_closure(strategy="stoichiometric_water")
    S_b, _, _ = build_petersen_matrix_with_oh_closure(strategy="oxygen_exact")
    np.testing.assert_allclose(S_a, S_b, rtol=0, atol=0.0)


def test_get_petersen_matrix_for_simulation_toggle() -> None:
    S_si = get_petersen_matrix_for_simulation(closure_mode="si")
    S_cl = get_petersen_matrix_for_simulation(closure_mode="oxygen")
    assert S_si.shape == S_cl.shape
    assert not np.allclose(S_si[:, S_H2O], S_cl[:, S_H2O])


@pytest.mark.parametrize(
    "strategy",
    ["stoichiometric_water", "oxygen_exact", "least_squares"],
)
def test_closure_strategies_return_finite_matrices(strategy: str) -> None:
    S, B, d = build_petersen_matrix_with_oh_closure(strategy=strategy)  # type: ignore[arg-type]
    assert np.all(np.isfinite(S))
    assert np.all(np.isfinite(B))
    assert d.strategy == strategy
