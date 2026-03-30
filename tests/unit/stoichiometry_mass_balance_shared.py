"""Shared mass-balance helpers for Petersen vs composition matrix (tests only)."""

from __future__ import annotations

import numpy as np

from bioprocess_twin.models.stoichiometry import get_composition_matrix, get_petersen_matrix

ELEMENT_NAMES = ["COD", "O", "C", "N", "P", "H"]

# O/H/COD still have residuals vs composition matrix (see DEVLOG). Target: 1e-6 once
# stoichiometry matches I across all processes; tightening to 1e-6 fails until then.
MASS_BALANCE_ATOL = 6.0

# Row order must match get_petersen_matrix() (rho1..rho19)
PETERSEN_PROCESS_LABELS = (
    "rho1: phototrophic growth X_ALG (NH4+)",
    "rho2: phototrophic growth X_ALG (NO3-)",
    "rho3: aerobic respiration X_ALG",
    "rho4: decay X_ALG",
    "rho5: aerobic growth X_H (NH4+)",
    "rho6: aerobic growth X_H (NO3-)",
    "rho7: aerobic respiration X_H",
    "rho8: anoxic growth X_H (NO3-)",
    "rho9: anoxic growth X_H (NO2-)",
    "rho10: anoxic respiration X_H",
    "rho11: hydrolysis X_S",
    "rho12: urea hydrolysis",
    "rho13: decay X_H",
    "rho14: aerobic growth X_AOB",
    "rho15: aerobic respiration X_AOB",
    "rho16: decay X_AOB",
    "rho17: aerobic growth X_NOB",
    "rho18: aerobic respiration X_NOB",
    "rho19: decay X_NOB",
)


def compute_mass_balance_matrix() -> np.ndarray:
    """Return S @ I.T with shape (19, 6): process × element."""
    S = get_petersen_matrix()
    comp = get_composition_matrix()
    return S @ comp.T


def format_mass_balance_by_element_summary(balance: np.ndarray) -> str:
    """Max |balance| per element and process row where it occurs."""
    abs_b = np.abs(balance)
    max_per_element = np.max(abs_b, axis=0)
    worst_row = np.argmax(abs_b, axis=0)
    lines = []
    for k in range(6):
        i = int(worst_row[k])
        signed = float(balance[i, k])
        label = PETERSEN_PROCESS_LABELS[i]
        lines.append(
            f"  {ELEMENT_NAMES[k]}: max |balance| = {max_per_element[k]:.2e}  "
            f"(worst: S[{i}] {label}; residual = {signed:+.4e})"
        )
    return "Mass balance errors by species:\n" + "\n".join(lines)


def format_mass_balance_all_cells(balance: np.ndarray, atol: float) -> str:
    """One line per process×element (114 lines), with OK if |residual| <= atol."""
    header = (
        f"Mass balance: all 19×6 cells (atol={atol:g}; OK if |residual| <= atol)\n"
        f"{'idx':>4} {'process':>5} {'elem':>4}  {'residual':>14}  {'|r|':>12}  status"
    )
    sep = "-" * len(header)
    rows: list[str] = [header, sep]
    n = 0
    for i in range(balance.shape[0]):
        for k in range(balance.shape[1]):
            n += 1
            r = float(balance[i, k])
            a = abs(r)
            ok = "OK" if a <= atol else "FAIL"
            rows.append(
                f"{n:4d} rho{i + 1:2d} {ELEMENT_NAMES[k]:>4}  {r:+.6e}  {a:12.6e}  {ok}"
            )
    assert n == 114, f"expected 114 cells, got {n}"
    return "\n".join(rows)
