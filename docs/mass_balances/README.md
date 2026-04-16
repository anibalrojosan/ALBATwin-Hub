# Mass balances and O/H closure documentation

This area documents **elemental audits** for the ALBA Petersen matrix \(\mathbf{S}\) and composition matrix \(\mathbf{I}\): for each biological process row and each composition row (COD, O, C, N, P, H), the check is \(\mathbf{B} = \mathbf{S}\mathbf{I}^\top\) with **114 cells** (19×6) for all closure policies below.

**Policy and architecture:** [ADR 007: Elemental mass-balance O/H closure](../adrs/007-elemental-mass-balance-oh-closure.md).

**Model matrices (SI vs extended):** [STOICHIOMETRY.md](../STOICHIOMETRY.md).

## Layout

| Kind | Path |
|------|------|
| Narrative & rationale | [`guides/`](guides/) |
| Case studies / worksheets | [`case-studies/`](case-studies/) |
| Generated 114-cell audits | [`artifacts/`](artifacts/) |

## Closure modes and artefacts

| Mode | \(\mathbf{S}\) shape | \(\mathbf{C}\) length | Generated audit | Regenerate |
|------|----------------------|------------------------|-----------------|------------|
| `si` | 19×17 | 17 | [`artifacts/audit-si-114cell.md`](artifacts/audit-si-114cell.md) | `uv run python scripts/generate_mass_balances_md.py` |
| `oxygen` | 19×17 (water column adjusted) | 17 | [`artifacts/audit-oxygen-closure-114cell.md`](artifacts/audit-oxygen-closure-114cell.md) | `uv run python scripts/generate_mass_balances_md.py --closure-of-oxygen` |
| `oxygen_and_protons` | 19×18 | 18 | [`artifacts/audit-oxygen-proton-closure-114cell.md`](artifacts/audit-oxygen-proton-closure-114cell.md) | `uv run python scripts/generate_mass_balances_md.py --closure-of-oxygen-and-protons` |

The biological rate vector \(\boldsymbol{\rho}\) always has **19** components (see [`MATH_MODEL.md`](../MATH_MODEL.md)); only \(\mathbf{S}\) and the state layout change with closure mode.

## Guides

- [OH_CLOSURE.md](guides/OH_CLOSURE.md) — overview of oxygen and proton closure layers.
- [stoichiometric_water_rationale.md](guides/stoichiometric_water_rationale.md) — closing elemental O via \(S_{\mathrm{H2O}}\).
- [proton_closure_rationale.md](guides/proton_closure_rationale.md) — \(S_{\mathrm{H\_PROTON}}\) column (audit layer).

## Case studies

- [rho14-aob-oxygen-balance.md](case-studies/rho14-aob-oxygen-balance.md) — oxygen residual for aerobic AOB growth (\(\rho_{14}\)).

## Root-level stubs

Index for all three audits (links + regeneration commands): [`docs/MASS_BALANCES.md`](../MASS_BALANCES.md).
