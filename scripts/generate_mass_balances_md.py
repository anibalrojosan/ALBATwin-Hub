#!/usr/bin/env python3
"""Regenerate docs/MASS_BALANCES.md from stoichiometry matrices (19×17 × 6×17 → 114 cells)."""

from __future__ import annotations

import sys
from pathlib import Path

# Repo root on sys.path when run as `uv run python scripts/generate_mass_balances_md.py`
_ROOT = Path(__file__).resolve().parents[1]
_src = str(_ROOT / "src")
if _src not in sys.path:
    sys.path.insert(0, _src)

from bioprocess_twin.models.stoichiometry import (  # noqa: E402
    get_composition_matrix,
    get_petersen_matrix,
)

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

MASS_BALANCE_ATOL = 0.1


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


def generate() -> str:
    S = get_petersen_matrix()
    comp = get_composition_matrix()
    balance = S @ comp.T

    lines: list[str] = []
    lines.append(
        "# Mass balances: explicit $\\mathbf{S}\\mathbf{I}^\\top$ cells (ALBA stoichiometry)"
    )
    lines.append("")
    lines.append(
        "This file documents the **elemental and COD-style mass-balance check** "
        "used in unit tests: for each biological process row $i$ of the Petersen "
        "matrix $\\mathbf{S}$ and each row $k$ of the composition matrix "
        "$\\mathbf{I}$ (COD, O, C, N, P, H),"
    )
    lines.append("")
    lines.append(r"$$B_{i,k} = \sum_{j=0}^{16} S_{i,j}\, I_{k,j}\,.$$")
    lines.append("")
    lines.append(
        "Process index $i$ is **0-based** (`S[i, :]`, same as code). Casagli process numbering "
        "$\\rho_r$ uses $r = i + 1$. A value near zero means the row is consistent with that "
        "composition row; large $|B_{i,k}|$ indicates missing or implicit species in $\\mathbf{S}$ "
        "(see ADR 007 for O/H policy)."
    )
    lines.append("")
    lines.append(
        "**Tolerance** for OK/FAIL labels matches "
        "`tests/unit/stoichiometry_mass_balance_shared.py`: "
        f"`MASS_BALANCE_ATOL = {MASS_BALANCE_ATOL:g}` "
        r"(|residual| $\le$ atol $\Rightarrow$ OK)."
    )
    lines.append("")
    lines.append("**Regeneration:** `uv run python scripts/generate_mass_balances_md.py`")
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
            st = "OK" if a <= MASS_BALANCE_ATOL else "FAIL"
            lines.append(f"| {rho} | {elem} | {_fmt_float(r)} | {_fmt_float(a)} | **{st}** |")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Per-cell derivations")
    lines.append("")
    lines.append(
        "For each cell: contributing columns are those with $S_{i,j} \\neq 0$. "
        "Products use numeric values from `get_petersen_matrix()` and `get_composition_matrix()`."
    )
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
            st = "**OK**" if a <= MASS_BALANCE_ATOL else "**FAIL**"
            lines.append(
                f"Balance for process row $i={i}$ ($\\rho_{{{rho}}}$) and composition row "
                f"$k={k}$ (**{elem}**):"
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
                lines.append(f"**Status:** {st} (atol = {MASS_BALANCE_ATOL:g})")
                lines.append("")
                continue

            terms: list[str] = []
            for j in js:
                s_ij = float(S[i, j])
                i_kj = float(comp[k, j])
                prod = s_ij * i_kj
                terms.append(_fmt_float(prod))
            lines.append(
                "Non-zero contributions ($S_{i,j}\\neq 0$), in column order $j = 0\\ldots 16$:"
            )
            lines.append("")
            lines.append(r"$$B_{i,k} = " + " + ".join(terms) + r"$$")
            lines.append("")

            lines.append("| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |")
            lines.append("|:---:|:---|:---:|:---:|:---:|")
            for j in js:
                s_ij = float(S[i, j])
                i_kj = float(comp[k, j])
                prod = s_ij * i_kj
                row = (
                    f"| {j} | `{STATE_LABELS[j]}` | {_fmt_float(s_ij)} | "
                    f"{_fmt_float(i_kj)} | {_fmt_float(prod)} |"
                )
                lines.append(row)
            lines.append("")
            lines.append(
                f"**Sum:** $B_{{{i},{k}}} \\approx {_fmt_float(r)}$ "
                f"(check: `numpy.sum` = {_fmt_float(r)})"
            )
            lines.append("")
            lines.append(
                f"**Status:** {st} — $|B_{{{i},{k}}}| = {_fmt_float(a)}$ "
                f"{'≤' if a <= MASS_BALANCE_ATOL else '>'} {MASS_BALANCE_ATOL:g}"
            )
            lines.append("")

    assert cell_idx == 114, cell_idx
    return "\n".join(lines) + "\n"


def main() -> None:
    out = _ROOT / "docs" / "MASS_BALANCES.md"
    text = generate()
    out.write_text(text, encoding="utf-8")
    print(f"Wrote {out} ({len(text)} bytes)")


if __name__ == "__main__":
    main()
