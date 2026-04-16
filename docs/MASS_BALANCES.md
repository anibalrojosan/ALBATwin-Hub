# Mass balance audits (114-cell checks)

This file is the **single entry point** for the three **generated** elemental audits. Each audit tabulates, for every biological process row and every composition row (COD, O, C, N, P, H), the balance residual with

$$
\mathbf{B} = \mathbf{S}\mathbf{I}^\top
$$

so there are always **114 cells** (19 × 6).

Full tables live under [`mass_balances/artifacts/`](mass_balances/artifacts/). **Do not edit those files by hand.**

**More context:** [mass_balances/README.md](mass_balances/README.md) (closure modes, guides, case studies, ADR 007).

---

## 1. SI-faithful audit (`si`)

**Artefact:** [`mass_balances/artifacts/audit-si-114cell.md`](mass_balances/artifacts/audit-si-114cell.md) — `get_petersen_matrix()`, Casagli SI.3.3.

Regenerate:

```bash
uv run python scripts/generate_mass_balances_md.py
```

---

## 2. Oxygen (water) closure (`oxygen`)

**Artefact:** [`mass_balances/artifacts/audit-oxygen-closure-114cell.md`](mass_balances/artifacts/audit-oxygen-closure-114cell.md) — stoichiometric **S_H2O** column; Petersen matrix **19×17**.

Regenerate:

```bash
uv run python scripts/generate_mass_balances_md.py --closure-of-oxygen
```

(`--closure` is a deprecated alias for `--closure-of-oxygen`.)

**Narrative:** [mass_balances/guides/OH_CLOSURE.md](mass_balances/guides/OH_CLOSURE.md), [mass_balances/guides/stoichiometric_water_rationale.md](mass_balances/guides/stoichiometric_water_rationale.md).

---

## 3. Oxygen + proton closure (`oxygen_and_protons`)

**Artefact:** [`mass_balances/artifacts/audit-oxygen-proton-closure-114cell.md`](mass_balances/artifacts/audit-oxygen-proton-closure-114cell.md) — extended Petersen matrix **19×18** (audit layer; no SI.6 pH coupling).

Regenerate:

```bash
uv run python scripts/generate_mass_balances_md.py --closure-of-oxygen-and-protons
```

**Narrative:** [mass_balances/guides/proton_closure_rationale.md](mass_balances/guides/proton_closure_rationale.md).
