"""Per-cell mass balance: 19 processes × 6 elements (114 checks).

Run with ``-s`` to see the full 114-line residual table from
``test_mass_balance_print_all_cell_residuals`` or from
``test_stoichiometry.py::test_mass_balance_conservation``.
"""

import numpy as np
import pytest
from stoichiometry_mass_balance_shared import (
    ELEMENT_NAMES,
    MASS_BALANCE_ATOL,
    PETERSEN_PROCESS_LABELS,
    compute_mass_balance_matrix,
    format_mass_balance_all_cells,
)


@pytest.fixture(scope="module")
def mass_balance_matrix() -> np.ndarray:
    """S @ I.T once per module."""
    return compute_mass_balance_matrix()


def test_mass_balance_print_all_cell_residuals(mass_balance_matrix: np.ndarray) -> None:
    """Print every process×element residual (use ``pytest -s``). Always passes."""
    print("\n" + format_mass_balance_all_cells(mass_balance_matrix, MASS_BALANCE_ATOL))


@pytest.mark.strict_si_mass_balance
@pytest.mark.parametrize("element_idx", range(6), ids=ELEMENT_NAMES)
@pytest.mark.parametrize("process_idx", range(19), ids=[f"rho{r + 1}" for r in range(19)])
def test_mass_balance_each_process_element(
    mass_balance_matrix: np.ndarray,
    process_idx: int,
    element_idx: int,
) -> None:
    """Each cell: balance[i,k] = sum_j S[i,j]*I[k,j] must be ≈ 0."""
    residual = float(mass_balance_matrix[process_idx, element_idx])
    assert abs(residual) <= MASS_BALANCE_ATOL, (
        f"S[{process_idx}] {PETERSEN_PROCESS_LABELS[process_idx]} / "
        f"{ELEMENT_NAMES[element_idx]}: balance = {residual:+.6e} "
        f"(|.| > atol={MASS_BALANCE_ATOL})"
    )
