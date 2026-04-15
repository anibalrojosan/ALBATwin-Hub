"""Tests for O + proton (H) extended Petersen closure (19×18 audit layer)."""

from __future__ import annotations

import numpy as np

from bioprocess_twin.models.stoichiometry import (
    N_STATE,
    S_H2O,
    get_composition_matrix,
    get_petersen_matrix,
)
from bioprocess_twin.models.stoichiometry_closure import (
    N_STATE_EXTENDED,
    S_H_PROTON,
    build_petersen_matrix_with_oh_closure,
    build_petersen_matrix_with_oxygen_and_proton_closure,
    compute_mass_balance_matrix,
    compute_stoichiometric_s_h_proton_total_for_row,
    get_composition_matrix_proton_closure,
    get_petersen_matrix_for_simulation,
)
from bioprocess_twin.stoichiometry_validation import MASS_BALANCE_ATOL


def test_composition_proton_closure_shape_and_h_unit() -> None:
    comp18 = get_composition_matrix_proton_closure()
    comp17 = get_composition_matrix()
    assert comp18.shape == (6, N_STATE_EXTENDED)
    np.testing.assert_allclose(comp18[:, :N_STATE], comp17, rtol=0, atol=0.0)
    assert comp18[5, S_H_PROTON] == 1.0
    assert np.all(comp18[:5, S_H_PROTON] == 0.0)


def test_extended_matrix_first_17_columns_match_oxygen_closure() -> None:
    S_o, _, _ = build_petersen_matrix_with_oh_closure()
    S_ext, _, _ = build_petersen_matrix_with_oxygen_and_proton_closure()
    assert S_ext.shape == (19, N_STATE_EXTENDED)
    np.testing.assert_allclose(S_ext[:, :N_STATE], S_o, rtol=0, atol=0.0)


def test_proton_coef_minus_residual_h_after_oxygen_closure() -> None:
    S_o, _, _ = build_petersen_matrix_with_oh_closure()
    comp17 = get_composition_matrix()
    S_ext, _, detail = build_petersen_matrix_with_oxygen_and_proton_closure()
    for i in range(19):
        b_h = float(np.dot(S_o[i, :], comp17[5, :]))
        beta = compute_stoichiometric_s_h_proton_total_for_row(i, S_o, comp17)
        np.testing.assert_allclose(beta, -b_h, rtol=0, atol=1e-12)
        if abs(b_h) <= MASS_BALANCE_ATOL:
            assert S_ext[i, S_H_PROTON] == 0.0
        else:
            np.testing.assert_allclose(S_ext[i, S_H_PROTON], beta, rtol=0, atol=1e-12)
            assert i in detail.s_h_proton_by_row


def test_proton_closure_h_and_oxygen_residuals_small() -> None:
    _, B_ext, detail = build_petersen_matrix_with_oxygen_and_proton_closure()
    for i in detail.rows_adjusted_proton:
        assert abs(B_ext[i, 5]) < 1e-9
    assert np.max(np.abs(B_ext[:, 1])) < 2e-3
    assert np.max(np.abs(B_ext[:, 5])) <= MASS_BALANCE_ATOL + 1e-12


def test_proton_closure_cod_c_n_p_unchanged_vs_si() -> None:
    S_si = get_petersen_matrix()
    comp17 = get_composition_matrix()
    B_si = compute_mass_balance_matrix(S_si, comp17)
    _, B_ext, _ = build_petersen_matrix_with_oxygen_and_proton_closure()
    for k in (0, 2, 3, 4):
        np.testing.assert_allclose(B_ext[:, k], B_si[:, k], rtol=0, atol=1e-12)


def test_get_petersen_matrix_for_simulation_shapes() -> None:
    assert get_petersen_matrix_for_simulation(closure_mode="si").shape == (19, N_STATE)
    assert get_petersen_matrix_for_simulation(closure_mode="oxygen").shape == (19, N_STATE)
    assert get_petersen_matrix_for_simulation(closure_mode="oxygen_and_protons").shape == (
        19,
        N_STATE_EXTENDED,
    )


def test_legacy_use_oh_closure_overrides_mode() -> None:
    S_a = get_petersen_matrix_for_simulation(use_oh_closure=True)
    S_b = get_petersen_matrix_for_simulation(closure_mode="oxygen")
    np.testing.assert_allclose(S_a, S_b, rtol=0, atol=0.0)
    S_si = get_petersen_matrix_for_simulation(use_oh_closure=False)
    np.testing.assert_allclose(S_si, get_petersen_matrix(), rtol=0, atol=0.0)


def test_proton_column_does_not_change_oxygen_balance() -> None:
    S_o, B_o, _ = build_petersen_matrix_with_oh_closure()
    _, B_ext, _ = build_petersen_matrix_with_oxygen_and_proton_closure()
    np.testing.assert_allclose(B_ext[:, 1], B_o[:, 1], rtol=0, atol=1e-12)
