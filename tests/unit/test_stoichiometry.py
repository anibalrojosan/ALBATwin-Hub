"""Unit tests for the ALBA stoichiometry module (Petersen and composition matrices)."""

import numpy as np
import pytest
from stoichiometry_mass_balance_shared import (
    MASS_BALANCE_ATOL,
    compute_mass_balance_matrix,
    format_mass_balance_all_cells,
    format_mass_balance_by_element_summary,
)

from bioprocess_twin.models.stoichiometry import (
    I_H_ND,
    _alpha_alg_o2_growth_nh4,
    _alpha_alg_o2_growth_no3,
    get_composition_matrix,
    get_petersen_matrix,
)

# Column indices for state variables (must match StateVector.to_array order)
X_ALG, X_AOB, X_NOB, X_H, X_S, X_I = 0, 1, 2, 3, 4, 5
S_S, S_I, S_IC, S_ND, S_NH, S_NO2, S_NO3, S_N2, S_PO4, S_O2, S_H2O = (
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
)


def test_composition_matrix_dimensions():
    """Assert composition matrix is 6 x 17 (6 elements: COD, O, C, N, P, H)."""
    comp = get_composition_matrix()
    assert comp.shape == (6, 17)


def test_petersen_matrix_dimensions():
    """Assert Petersen matrix is 19 x 17."""
    S = get_petersen_matrix()
    assert S.shape == (19, 17)


@pytest.mark.strict_si_mass_balance
def test_mass_balance_conservation():
    """For each process, sum(S[i,:] * I[k,:]) ≈ 0 for all elements (COD, O, C, N, P, H)."""
    balance = compute_mass_balance_matrix()
    print("\n" + format_mass_balance_by_element_summary(balance))
    print("\n" + format_mass_balance_all_cells(balance, MASS_BALANCE_ATOL))
    assert np.allclose(balance, 0, atol=MASS_BALANCE_ATOL), (
        f"Mass balance violated (atol={MASS_BALANCE_ATOL}).\n"
        + format_mass_balance_by_element_summary(balance)
        + "\n"
        + format_mass_balance_all_cells(balance, MASS_BALANCE_ATOL)
    )


def test_petersen_matrix_structure():
    """Spot-check known non-zero entries from the Petersen matrix."""
    S = get_petersen_matrix()

    # rho1/rho2: S_O2 from Table SI.3.3 (not +1)
    assert S[0, X_ALG] == 1.0
    assert S[0, S_O2] == _alpha_alg_o2_growth_nh4()
    assert S[1, X_ALG] == 1.0
    assert S[1, S_O2] == _alpha_alg_o2_growth_no3()

    # Process 4 (Algae decay): X_ALG=-1, X_S and X_I produced (positive)
    assert S[3, X_ALG] == -1.0
    assert S[3, X_S] > 0
    assert S[3, X_I] > 0

    # Process 5 (Aerobic growth on NH4): X_H=+1
    assert S[4, X_H] == 1.0

    # rho12 (Urea hydrolysis): Table SI.3.3
    assert S[11, S_ND] == -1.0
    assert S[11, S_NH] == 1.0
    assert S[11, S_H2O] == I_H_ND


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
