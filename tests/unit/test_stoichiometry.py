"""Unit tests for the ALBA stoichiometry module (Petersen and composition matrices)."""

import numpy as np

from bioprocess_twin.models.stoichiometry import get_composition_matrix, get_petersen_matrix

# Column indices for state variables (must match StateVector.to_array order)
X_ALG, X_AOB, X_NOB, X_H, X_S, X_I = 0, 1, 2, 3, 4, 5
S_S, S_I, S_IC, S_ND, S_NH, S_NO2, S_NO3, S_N2, S_PO4, S_O2, S_H2O = (
    6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
)


def test_composition_matrix_dimensions():
    """Assert composition matrix is 6 x 17 (6 elements: COD, O, C, N, P, H)."""
    comp = get_composition_matrix()
    assert comp.shape == (6, 17)


def test_petersen_matrix_dimensions():
    """Assert Petersen matrix is 19 x 17."""
    S = get_petersen_matrix()
    assert S.shape == (19, 17)


ELEMENT_NAMES = ["COD", "O", "C", "N", "P", "H"]
# O balance residual ~5.56 in AOB process (ALBA conventions)
MASS_BALANCE_ATOL = 6.0


def _format_balance_errors(balance: np.ndarray) -> str:
    """Format max |balance| per element for assertion message."""
    max_per_element = np.max(np.abs(balance), axis=0)
    lines = [
        f"  {ELEMENT_NAMES[k]}: max |balance| = {max_per_element[k]:.2e}"
        for k in range(6)
    ]
    return "Mass balance errors by species:\n" + "\n".join(lines)


def test_mass_balance_conservation():
    """For each process, sum(S[i,:] * I[k,:]) ≈ 0 for all elements (COD, O, C, N, P, H)."""
    S = get_petersen_matrix()
    comp = get_composition_matrix()
    balance = S @ comp.T  # 19 x 6 (rows: COD, O, C, N, P, H)
    # Always show per-species errors after each test run
    print("\n" + _format_balance_errors(balance))
    assert np.allclose(balance, 0, atol=MASS_BALANCE_ATOL), (
        f"Mass balance violated (atol={MASS_BALANCE_ATOL}).\n"
        + _format_balance_errors(balance)
    )


def test_petersen_matrix_structure():
    """Spot-check known non-zero entries from the Petersen matrix."""
    S = get_petersen_matrix()

    # Process 1 (Phototrophic growth on NH4): X_ALG=+1, S_O2=+1
    assert S[0, X_ALG] == 1.0
    assert S[0, S_O2] == 1.0

    # Process 4 (Algae decay): X_ALG=-1, X_S and X_I produced (positive)
    assert S[3, X_ALG] == -1.0
    assert S[3, X_S] > 0
    assert S[3, X_I] > 0

    # Process 5 (Aerobic growth on NH4): X_H=+1
    assert S[4, X_H] == 1.0

    # Process 12 (Urea hydrolysis): S_ND=-1, S_NH=+1
    assert S[11, S_ND] == -1.0
    assert S[11, S_NH] == 1.0


def test_composition_matrix_non_negative_where_expected():
    """Biomass/organic columns have positive composition; inorganic N has expected structure."""
    comp = get_composition_matrix()
    # Rows: 0=COD, 1=O, 2=C, 3=N, 4=P, 5=H

    # Biomass columns (0-5): COD row should be 1
    for j in range(6):
        assert comp[0, j] == 1.0, f"COD for biomass/organic col {j} should be 1"

    # S_IC (col 8): Carbon = 1
    assert comp[2, S_IC] == 1.0

    # S_NH (col 10): Nitrogen = 1
    assert comp[3, S_NH] == 1.0

    # S_PO4 (col 14): Phosphorus = 1
    assert comp[4, S_PO4] == 1.0
