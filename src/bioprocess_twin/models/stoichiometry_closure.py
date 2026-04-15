"""
Experimental O/H closure layer on top of the SI-faithful Petersen matrix.

`get_petersen_matrix()` remains the literal ALBA SI transcription. This module
builds a **copy** with a **stoichiometrically determined** coefficient on column
``S_H2O`` for rows that fail the elemental O/H audit.

**Stoichiometric principle (default strategy)**

For process row :math:`i`, let :math:`L_i=\\sum_{j\\neq S_{\\mathrm{H2O}}} S_{i,j} I_{\\mathrm{O},j}`.
That is the elemental oxygen carried by every **explicit** SI species in the row.
If all oxygen that is still missing from the ASM-style formulation is assigned to
**water oxygen** (via the composition row for :math:`S_{\\mathrm{H2O}}` in g H), then
the **unique** total Petersen coefficient :math:`\\alpha_i` on :math:`S_{\\mathrm{H2O}}`
that satisfies :math:`B_{i,\\mathrm{O}}=0` is

.. math::

    \\alpha_i = -\\frac{L_i}{I_{\\mathrm{O},S_{\\mathrm{H2O}}}}.

This is algebraically the same as adjusting the SI row by
:math:`\\Delta_i=-B_{i,\\mathrm{O}}^{\\mathrm{SI}}/I_{\\mathrm{O},S_{\\mathrm{H2O}}}` when
only the water column moves, but the **interpretation** is stoichiometric: :math:`\\alpha_i`
is the net **g H** (per process normalization) tied to water formation implied by the SI
reaction once water is made explicit. See ``docs/mass_balances/stoichiometric_water_rationale.md``.

Hydrogen and charge may still not close with a single water column; see ADR 007.

Simulation note: use ``get_petersen_matrix_for_simulation`` to select SI vs closure
matrix; the ODE stack is not wired here yet.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from bioprocess_twin.stoichiometry_validation import MASS_BALANCE_ATOL

from .stoichiometry import (
    N_STATE,
    S_H2O,
    get_composition_matrix,
    get_petersen_matrix,
)

ClosureStrategy = Literal["stoichiometric_water", "oxygen_exact", "least_squares"]

# One-line chemistry anchors for the auxiliary water term (interpretation only).
# Coefficients still follow L_i / I_O_S_H2O from the SI row + composition matrix.
AUXILIARY_WATER_PROCESS_NOTES: tuple[str, ...] = (
    "rho1: phototrophic growth on NH4+ (SI alpha_H2O; no change unless audit fails).",
    "rho2: phototrophic growth on NO3- (SI alpha_H2O).",
    "rho3: algal aerobic respiration (SI alpha_H2O).",
    "rho4: algal decay — implicit solvent/water in ASM-style decay products.",
    "rho5: heterotrophic aerobic growth on NH4+ — O in H2O from redox + assimilation.",
    "rho6: heterotrophic aerobic growth on NO3- — same with nitrate as N source.",
    "rho7: heterotrophic aerobic respiration.",
    "rho8: anoxic growth on NO3- (N2 gas).",
    "rho9: anoxic growth on NO2- (N2 gas).",
    "rho10: anoxic respiration (NO2/NO3 to N2).",
    "rho11: hydrolysis X_S — organic fragmentation; auxiliary water from O bookkeeping.",
    "rho12: urea hydrolysis — SI gives alpha_H2O=I_H_ND; O closure refines total water-H.",
    "rho13: heterotrophic decay.",
    "rho14: AOB growth on NH4+ — net H2O from NH4 oxidation + autotrophic synthesis (see balance worksheet).",
    "rho15: AOB respiration.",
    "rho16: AOB decay.",
    "rho17: NOB growth on NO2-.",
    "rho18: NOB respiration.",
    "rho19: NOB decay.",
)


@dataclass
class OHClosureDetail:
    """Diagnostics for ``build_petersen_matrix_with_oh_closure``."""

    strategy: ClosureStrategy
    atol: float
    rows_adjusted: tuple[int, ...] = field(default_factory=tuple)
    delta_h2o_by_row: dict[int, float] = field(default_factory=dict)
    s_h2o_total_by_row: dict[int, float] = field(default_factory=dict)
    balance_si: np.ndarray | None = None
    balance_closed: np.ndarray | None = None


def compute_mass_balance_matrix(
    S: np.ndarray,
    comp: np.ndarray | None = None,
) -> np.ndarray:
    """Return ``S @ comp.T`` with shape (19, 6)."""
    if comp is None:
        comp = get_composition_matrix()
    return S @ comp.T


def compute_stoichiometric_s_h2o_total_for_row(
    i: int,
    S_si: np.ndarray,
    comp: np.ndarray,
) -> float:
    """
    Total :math:`S_{i,S_{\\mathrm{H2O}}}` that zeros elemental **oxygen** for row ``i``.

    Uses :math:`L_i=\\sum_{j\\neq S_{\\mathrm{H2O}}} S_{i,j} I_{\\mathrm{O},j}` and
    :math:`\\alpha_i=-L_i/I_{\\mathrm{O},S_{\\mathrm{H2O}}}` (see module docstring).
    """
    w_o = float(comp[1, S_H2O])
    if w_o == 0.0:
        return float(S_si[i, S_H2O])
    net_o = 0.0
    for j in range(N_STATE):
        if j == S_H2O:
            continue
        net_o += float(S_si[i, j]) * float(comp[1, j])
    return -net_o / w_o


def list_oh_mass_balance_violations(
    balance: np.ndarray | None = None,
    *,
    atol: float | None = None,
) -> list[tuple[int, float, float]]:
    """
    Inventory: process rows (0-based) where oxygen or hydrogen balance fails.

    Returns tuples ``(i, B_{i,O}, B_{i,H})`` for rows with
    ``|B_{i,O}| > atol`` or ``|B_{i,H}| > atol``.
    """
    if atol is None:
        atol = MASS_BALANCE_ATOL
    if balance is None:
        balance = compute_mass_balance_matrix(get_petersen_matrix())
    out: list[tuple[int, float, float]] = []
    for i in range(balance.shape[0]):
        b_o = float(balance[i, 1])
        b_h = float(balance[i, 5])
        if abs(b_o) > atol or abs(b_h) > atol:
            out.append((i, b_o, b_h))
    return out


def _delta_h2o_least_squares(
    b_o: float,
    b_h: float,
    w_o: float,
    w_h: float,
) -> float:
    """Minimize (b_o + d*w_o)^2 + (b_h + d*w_h)^2 in d (increment to SI S_H2O)."""
    den = w_o * w_o + w_h * w_h
    if den == 0.0:
        return 0.0
    return -(b_o * w_o + b_h * w_h) / den


def build_petersen_matrix_with_oh_closure(
    *,
    atol: float | None = None,
    strategy: ClosureStrategy = "stoichiometric_water",
) -> tuple[np.ndarray, np.ndarray, OHClosureDetail]:
    """
    Copy SI Petersen matrix and set ``S[i, S_H2O]`` from closure strategy.

    * ``stoichiometric_water`` / ``oxygen_exact``: total water-H coefficient
      :math:`\\alpha_i=-L_i/I_{\\mathrm{O},S_{\\mathrm{H2O}}}` (same numeric result;
      ``oxygen_exact`` kept as an alias).
    * ``least_squares``: add increment ``d`` to SI ``S_H2O`` minimizing O/H sum of squares.

    Non-``S_{H2O}`` entries match the SI matrix.
    """
    if atol is None:
        atol = MASS_BALANCE_ATOL

    strat: ClosureStrategy = "stoichiometric_water" if strategy == "oxygen_exact" else strategy

    S_si = get_petersen_matrix()
    comp = get_composition_matrix()
    B_si = compute_mass_balance_matrix(S_si, comp)

    w_o = float(comp[1, S_H2O])
    w_h = float(comp[5, S_H2O])

    S_closed = S_si.copy()
    deltas: dict[int, float] = {}
    totals: dict[int, float] = {}
    rows_adj: list[int] = []

    for i in range(S_si.shape[0]):
        b_o = float(B_si[i, 1])
        b_h = float(B_si[i, 5])
        if abs(b_o) <= atol and abs(b_h) <= atol:
            continue

        s0 = float(S_si[i, S_H2O])
        if strat == "least_squares":
            delta = _delta_h2o_least_squares(b_o, b_h, w_o, w_h)
            s1 = s0 + delta
        else:
            s1 = compute_stoichiometric_s_h2o_total_for_row(i, S_si, comp)
            delta = s1 - s0

        S_closed[i, S_H2O] = s1
        deltas[i] = delta
        totals[i] = s1
        rows_adj.append(i)

    B_closed = compute_mass_balance_matrix(S_closed, comp)

    detail = OHClosureDetail(
        strategy=strategy,
        atol=atol,
        rows_adjusted=tuple(rows_adj),
        delta_h2o_by_row=deltas,
        s_h2o_total_by_row=totals,
        balance_si=B_si,
        balance_closed=B_closed,
    )
    return S_closed, B_closed, detail


def get_petersen_matrix_for_simulation(*, use_oh_closure: bool = False) -> np.ndarray:
    """
    Return Petersen matrix for ``dC/dt = S^T rho``.

    Parameters
    ----------
    use_oh_closure:
        If True, use ``build_petersen_matrix_with_oh_closure`` (stoichiometric water
        on ``S_H2O``). If False, use literal SI ``get_petersen_matrix()``.

    Note
    ----
    When ``use_oh_closure`` is True, validate impacts on ``S_H2O`` dynamics and
    pH / hydrochemistry before production runs (not yet integrated in solver).
    """
    if use_oh_closure:
        S, _, _ = build_petersen_matrix_with_oh_closure()
        return S
    return get_petersen_matrix()


def format_closure_inventory_markdown() -> str:
    """Short Markdown table of SI O/H violations at ``MASS_BALANCE_ATOL``."""
    viol = list_oh_mass_balance_violations()
    lines = [
        "# O/H mass-balance inventory (SI Petersen)",
        "",
        f"Tolerance: `MASS_BALANCE_ATOL = {MASS_BALANCE_ATOL:g}` (`stoichiometry_validation.py`).",
        "",
        "| $\\rho$ | $i$ | $B_{i,\\mathrm{O}}$ | $B_{i,\\mathrm{H}}$ |",
        "|:---:|:---:|:---:|:---:|",
    ]
    for i, b_o, b_h in viol:
        lines.append(f"| {i + 1} | {i} | {b_o:.6e} | {b_h:.6e} |")
    if not viol:
        lines.append("| — | — | *(no violations)* | |")
    lines.append("")
    lines.append("See also `docs/MASS_BALANCES.md` for the full 114-cell audit.")
    return "\n".join(lines) + "\n"
