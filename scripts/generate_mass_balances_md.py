#!/usr/bin/env python3
"""Regenerate mass-balance Markdown from Petersen and composition matrices."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# Repo root on sys.path when run as `uv run python scripts/generate_mass_balances_md.py`
_ROOT = Path(__file__).resolve().parents[1]
_ARTIFACTS_DIR = _ROOT / "docs" / "mass_balances" / "artifacts"
_src = str(_ROOT / "src")
if _src not in sys.path:
    sys.path.insert(0, _src)

from bioprocess_twin.models.stoichiometry import (  # noqa: E402
    get_composition_matrix,
    get_petersen_matrix,
)
from bioprocess_twin.stoichiometry_validation import MASS_BALANCE_ATOL  # noqa: E402

STATE_LABELS = (
    "X_ALG",
    "X_AOB",
    "X_NOB",
    "X_H",
    "X_S",
    "X_I",
    "S_S",
    "S_I",
    "S_IC",
    "S_ND",
    "S_NH",
    "S_NO2",
    "S_NO3",
    "S_N2",
    "S_PO4",
    "S_O2",
    "S_H2O",
)

STATE_LABELS_EXTENDED = STATE_LABELS + ("S_H_PROTON",)

ELEMENT_NAMES = ("COD", "O", "C", "N", "P", "H")

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


def _fmt_float(x: float) -> str:
    if x == 0.0:
        return "0"
    ax = abs(x)
    if ax >= 1e4 or ax < 1e-4:
        return f"{x:.6e}"
    s = f"{x:.10f}".rstrip("0").rstrip(".")
    return s if s else "0"


def _anchor_id(rho: int, elem: str) -> str:
    return f"rho{rho}-{elem}"


def generate_markdown(
    S: np.ndarray,
    comp: np.ndarray,
    *,
    title: str,
    preamble_paragraphs: list[str],
    atol: float,
    regen_command: str,
    derivation_blurb: str,
) -> str:
    """Build full 114-cell Markdown for given ``S`` and ``comp``."""
    balance = S @ comp.T
    n_j = S.shape[1]
    if comp.shape[1] != n_j:
        raise ValueError(f"S has {n_j} columns but composition has {comp.shape[1]}")
    j_max = n_j - 1
    state_labels = STATE_LABELS_EXTENDED if n_j == len(STATE_LABELS_EXTENDED) else STATE_LABELS
    if len(state_labels) != n_j:
        raise ValueError(f"Need state_labels length {n_j}, got {len(state_labels)}")

    lines: list[str] = []
    lines.append(title)
    lines.append("")
    for p in preamble_paragraphs:
        lines.append(p)
        lines.append("")
    lines.append(rf"$$B_{{i,k}} = \sum_{{j=0}}^{{{j_max}}} S_{{i,j}}\, I_{{k,j}}\,.$$")
    lines.append("")
    lines.append(
        "Process index $i$ is **0-based** (`S[i, :]`, same as code). Casagli process numbering "
        "$\\rho_r$ uses $r = i + 1$."
    )
    lines.append("")
    lines.append(
        f"**Tolerance:** audit uses atol = `{atol:g}` (see `stoichiometry_validation.py`). "
        r"OK if $|B_{i,k}| \le$ atol."
    )
    lines.append("")
    lines.append(f"**Regeneration:** `{regen_command}`")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Table of contents")
    lines.append("")
    lines.append("- [Summary table (all 114 cells)](#summary-table-all-114-cells)")
    lines.append("- [Per-cell derivations](#per-cell-derivations)")
    lines.append("")
    for i in range(19):
        rho = i + 1
        lines.append(f"- **rho{rho}** ({PETERSEN_PROCESS_LABELS[i].split(':', 1)[1].strip()})")
        sub: list[str] = []
        for _k, elem in enumerate(ELEMENT_NAMES):
            aid = _anchor_id(rho, elem)
            sub.append(f"  - [{elem}](#{aid})")
        lines.append("\n".join(sub))
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Summary table (all 114 cells)")
    lines.append("")
    lines.append("| $\\rho$ | Element | Residual $B_{i,k}$ | $\\|B_{i,k}\\|$ | Status |")
    lines.append("|:---:|:---:|:---:|:---:|:---:|")
    for i in range(19):
        rho = i + 1
        for k, elem in enumerate(ELEMENT_NAMES):
            r = float(balance[i, k])
            a = abs(r)
            st = "OK" if a <= atol else "FAIL"
            lines.append(f"| {rho} | {elem} | {_fmt_float(r)} | {_fmt_float(a)} | **{st}** |")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Per-cell derivations")
    lines.append("")
    lines.append(derivation_blurb)
    lines.append("")

    cell_idx = 0
    for i in range(19):
        rho = i + 1
        short = PETERSEN_PROCESS_LABELS[i]
        lines.append(f"### rho{rho}")
        lines.append("")
        lines.append(f"*Full label:* {short}")
        lines.append("")
        for k, elem in enumerate(ELEMENT_NAMES):
            cell_idx += 1
            aid = _anchor_id(rho, elem)
            lines.append(f'<a id="{aid}"></a>')
            lines.append(f"#### rho{rho} — {elem} (cell {cell_idx}/114)")
            lines.append("")
            r = float(balance[i, k])
            a = abs(r)
            st = "**OK**" if a <= atol else "**FAIL**"
            lines.append(
                f"Balance for process row $i={i}$ ($\\rho_{{{rho}}}$) and composition row $k={k}$ (**{elem}**):"
            )
            lines.append("")
            lines.append(r"$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$")
            lines.append("")

            js = [j for j in range(S.shape[1]) if S[i, j] != 0.0]
            if not js:
                lines.append(r"No columns with $S_{i,j}\neq 0$; $B_{i,k}=0$ by convention.")
                lines.append("")
                lines.append(f"**Sum:** $B_{{{i},{k}}} = 0$  ")
                lines.append("")
                lines.append(f"**Status:** {st} (atol = {atol:g})")
                lines.append("")
                continue

            terms: list[str] = []
            for j in js:
                s_ij = float(S[i, j])
                i_kj = float(comp[k, j])
                prod = s_ij * i_kj
                terms.append(_fmt_float(prod))
            lines.append(f"Non-zero contributions ($S_{{i,j}}\\neq 0$), in column order $j = 0\\ldots {j_max}$:")
            lines.append("")
            lines.append(r"$$B_{i,k} = " + " + ".join(terms) + r"$$")
            lines.append("")

            lines.append("| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |")
            lines.append("|:---:|:---|:---:|:---:|:---:|")
            for j in js:
                s_ij = float(S[i, j])
                i_kj = float(comp[k, j])
                prod = s_ij * i_kj
                row = f"| {j} | `{state_labels[j]}` | {_fmt_float(s_ij)} | {_fmt_float(i_kj)} | {_fmt_float(prod)} |"
                lines.append(row)
            lines.append("")
            lines.append(f"**Sum:** $B_{{{i},{k}}} \\approx {_fmt_float(r)}$ (check: `numpy.sum` = {_fmt_float(r)})")
            lines.append("")
            lines.append(f"**Status:** {st} — $|B_{{{i},{k}}}| = {_fmt_float(a)}$ {'≤' if a <= atol else '>'} {atol:g}")
            lines.append("")

    assert cell_idx == 114, cell_idx
    return "\n".join(lines) + "\n"


def _generate_si() -> str:
    S = get_petersen_matrix()
    comp = get_composition_matrix()
    return generate_markdown(
        S,
        comp,
        title="# Mass balances: explicit $\\mathbf{S}\\mathbf{I}^\\top$ cells (ALBA stoichiometry)",
        preamble_paragraphs=[
            "**Artifact:** `docs/mass_balances/artifacts/audit-si-114cell.md` — SI-faithful "
            "`get_petersen_matrix()` (Casagli SI.3.3).",
            "This file documents the **elemental and COD-style mass-balance check** "
            "used in unit tests: for each biological process row $i$ of the Petersen "
            "matrix $\\mathbf{S}$ and each row $k$ of the composition matrix "
            "$\\mathbf{I}$ (COD, O, C, N, P, H). Coefficients follow **Casagli SI.3.3** "
            "via `get_petersen_matrix()`. Large $|B_{i,k}|$ for O/H indicates ASM-style "
            "implicit solvent (see ADR 007).",
            "**Source matrices:** `get_petersen_matrix()`, `get_composition_matrix()`. "
            "**Tolerance constant:** `MASS_BALANCE_ATOL` in "
            "`src/bioprocess_twin/stoichiometry_validation.py` "
            f"(currently `{MASS_BALANCE_ATOL:g}`; tests re-export via "
            "`tests/unit/stoichiometry_mass_balance_shared.py`).",
        ],
        atol=MASS_BALANCE_ATOL,
        regen_command="uv run python scripts/generate_mass_balances_md.py",
        derivation_blurb=(
            "For each cell: contributing columns are those with $S_{i,j} \\neq 0$. "
            "Products use numeric values from the same ``S`` and ``comp`` as this document."
        ),
    )


def _generate_closure() -> str:
    from bioprocess_twin.models.stoichiometry_closure import build_petersen_matrix_with_oh_closure

    S, _, detail = build_petersen_matrix_with_oh_closure()
    comp = get_composition_matrix()
    rows = ", ".join(str(i + 1) for i in detail.rows_adjusted)
    return generate_markdown(
        S,
        comp,
        title=(
            "# Mass balances: $\\mathbf{S}\\mathbf{I}^\\top$ with **stoichiometric** "
            "$S_{\\mathrm{H2O}}$ closure (residual O as water-O)"
        ),
        preamble_paragraphs=[
            "**Artifact:** `docs/mass_balances/artifacts/audit-oxygen-closure-114cell.md` — "
            "oxygen closure via stoichiometric $S_{\\mathrm{H2O}}$ only.",
            "This file uses **`build_petersen_matrix_with_oh_closure()`** "
            "(`stoichiometry_closure.py`): a **copy** of the SI Petersen matrix with "
            "column $S_{\\mathrm{H2O}}$ set to $\\alpha_i=-L_i/I_{\\mathrm{O},S_{\\mathrm{H2O}}}$ "
            "for failing rows ($L_i$ = elemental O from all non-water SI species). "
            "That is **stoichiometric water** under the convention that missing oxygen is "
            "water oxygen in this $\\mathbf{I}$ basis (see "
            "`docs/mass_balances/guides/stoichiometric_water_rationale.md`). "
            "**Hydrogen** may still deviate until a second closure species (e.g. $\\mathrm{H}^+$) "
            "is modeled (ADR 007).",
            f"**Adjusted process rows (1-based $\\rho$):** {rows or '*(none)*'}. "
            "All non-$S_{\\mathrm{H2O}}$ entries match `get_petersen_matrix()`.",
        ],
        atol=MASS_BALANCE_ATOL,
        regen_command="uv run python scripts/generate_mass_balances_md.py --closure-of-oxygen",
        derivation_blurb=(
            "Same layout as `audit-si-114cell.md`; coefficients come from the **oxygen-closure** "
            "Petersen matrix ``S`` (this document)."
        ),
    )


def _generate_closure_oxygen_and_protons() -> str:
    from bioprocess_twin.models.stoichiometry_closure import (
        build_petersen_matrix_with_oxygen_and_proton_closure,
        get_composition_matrix_proton_closure,
    )

    S, _, detail = build_petersen_matrix_with_oxygen_and_proton_closure()
    comp = get_composition_matrix_proton_closure()
    rows_o = ", ".join(str(i + 1) for i in detail.oxygen.rows_adjusted)
    rows_p = ", ".join(str(i + 1) for i in detail.rows_adjusted_proton)
    return generate_markdown(
        S,
        comp,
        title=(
            "# Mass balances: $\\mathbf{S}\\mathbf{I}^\\top$ with **O** (water) + **H** "
            "(proton) closure — extended **19×18** $\\mathbf{S}$"
        ),
        preamble_paragraphs=[
            "**Artifact:** `docs/mass_balances/artifacts/audit-oxygen-proton-closure-114cell.md` "
            "— audit layer only (**not** SI Casagli 17 states; **no** charge / SI.6 pH).",
            "This file uses **`build_petersen_matrix_with_oxygen_and_proton_closure()`**: "
            "oxygen closure on **`S_H2O`**, then **`S_H_PROTON`** with "
            "**`compute_stoichiometric_s_h_proton_total_for_row`** "
            "($\\beta_i=-R_i^{\\mathrm{H}}$: all H not in the 17 explicit columns booked as "
            "g H in free protons; not a separate tuning parameter). "
            "Non-zero only where $|\\beta_i| > \\texttt{atol}$. See "
            "`docs/mass_balances/guides/proton_closure_rationale.md`.",
            f"**O-adjusted process rows (1-based $\\rho$):** {rows_o or '*(none)*'}. "
            f"**Proton-adjusted rows:** {rows_p or '*(none)*'}. "
            "Columns **0–16** match `build_petersen_matrix_with_oh_closure()`.",
        ],
        atol=MASS_BALANCE_ATOL,
        regen_command=("uv run python scripts/generate_mass_balances_md.py --closure-of-oxygen-and-protons"),
        derivation_blurb=(
            "114 cells = 19 processes × 6 composition rows; sums run over **18** columns where "
            "applicable. Matrix ``S`` is **19×18**; ``comp`` is **6×18**."
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--closure-of-oxygen",
        action="store_true",
        help="Write docs/mass_balances/artifacts/audit-oxygen-closure-114cell.md.",
    )
    parser.add_argument(
        "--closure-of-oxygen-and-protons",
        action="store_true",
        help="Write docs/mass_balances/artifacts/audit-oxygen-proton-closure-114cell.md.",
    )
    parser.add_argument(
        "--closure",
        action="store_true",
        help="Deprecated alias for --closure-of-oxygen.",
    )
    args = parser.parse_args()

    n_closure = sum(
        [
            args.closure_of_oxygen,
            args.closure_of_oxygen_and_protons,
            args.closure,
        ]
    )
    if n_closure > 1:
        parser.error("Use at most one of --closure-of-oxygen, --closure-of-oxygen-and-protons, --closure.")

    _ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)
    if args.closure_of_oxygen or args.closure:
        text = _generate_closure()
        out = _ARTIFACTS_DIR / "audit-oxygen-closure-114cell.md"
    elif args.closure_of_oxygen_and_protons:
        text = _generate_closure_oxygen_and_protons()
        out = _ARTIFACTS_DIR / "audit-oxygen-proton-closure-114cell.md"
    else:
        text = _generate_si()
        out = _ARTIFACTS_DIR / "audit-si-114cell.md"

    out.write_text(text, encoding="utf-8")
    print(f"Wrote {out} ({len(text)} bytes)")


if __name__ == "__main__":
    main()
