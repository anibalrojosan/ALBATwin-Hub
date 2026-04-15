# Rationale: **`S_H_PROTON`** as stoichiometric hydrogen in free protons (audit layer)

This document defines the **18th Petersen column** used in **`build_petersen_matrix_with_oxygen_and_proton_closure()`** ([`stoichiometry_closure.py`](../../src/bioprocess_twin/models/stoichiometry_closure.py)). It does **not** implement charge balance, SI.6 speciation, or pH dynamics.

The coefficient on this column is **not** a free knob to “make the H test pass”: it is the **only** value consistent with (i) the oxygen-closed row **`S^{(O)}`**, (ii) the **17-column** composition matrix **`I`**, and (iii) the modelling convention that **all elemental hydrogen not already carried by those columns** is booked as **g H in free solvated protons** with **`I_H = 1`**.

## Layering

1. **Oxygen (water)** — `build_petersen_matrix_with_oh_closure`: column **`S_H2O`** is set stoichiometrically so **`B_{i,O} = 0`** (see [`stoichiometric_water_rationale.md`](stoichiometric_water_rationale.md)).

2. **Hydrogen (protons)** — over the **same** first 17 columns, define the **explicit** hydrogen load

$$
R_i^{\mathrm{H}} = \sum_{j=0}^{16} S_{i,j}^{(\mathrm{O})}\, I_{\mathrm{H},j}.
$$

Every term is fixed by the SI row (after O closure) and **`I`**: biomass, nutrients, inorganics, and **water-H** are all included in this sum.

## Stoichiometric assignment to the proton pool

Extend **`I`** to **6×18** with a column where only **`I_H = 1`** (g H per unit state in the g H·m⁻³ convention; other composition rows **0**). If **all** hydrogen that must still appear in the process balance is assigned to that **single** inventory — the acid–base carrier in this audit layer — then **`B_{i,\mathrm{H}} = 0`** requires

$$
R_i^{\mathrm{H}} + \beta_i \cdot 1 = 0
\quad\Rightarrow\quad
\beta_i = -R_i^{\mathrm{H}}.
$$

So **`β_i`** is the **total** Petersen entry on **`S_H_PROTON`**: net **g H** (per process normalization) attributed to free protons. Code: **`compute_stoichiometric_s_h_proton_total_for_row`**.

Algebraically **`β_i = -B_{i,\mathrm{H}}^{(\mathrm{O})}`** where **`B_{i,\mathrm{H}}^{(\mathrm{O})} = R_i^{\mathrm{H}}`** before the proton column exists; the **interpretation** is the **`R → β`** construction above, parallel to **`L_i → α_i`** for water-O.

Rows with **`|β_i| ≤ MASS_BALANCE_ATOL`** keep **`S_{i,S_{\mathrm{H}^+}} = 0`** (no numerical proton term when the explicit row already closes H within tolerance).

## State convention: **`S_H_PROTON`** (column index 17)

**Unit:** **g H·m⁻³** of hydrogen as **free solvated protons** (same hydrogen-mass basis as **`S_H2O`** in `MATH_MODEL.md`).

$$
(\text{g H from } \mathrm{H}^+\text{ per m}^3) \approx c_{\mathrm{H}^+}\, M_{\mathrm{H}},
$$

with **`M_{\mathrm{H}} \approx 1`** g·mol⁻¹ when concentration is in mol·m⁻³.

## Scope limits

- **Not** electroneutrality; **not** SI.6 Table SI.6.1 charge balance.
- **`StateVector`** supports **17** or **18** components via **`StateVectorVariant`** (`core/state.py`); the ODE RHS still must use matching **`S`** shape.

## Related

- Water-O closure: [`stoichiometric_water_rationale.md`](stoichiometric_water_rationale.md).
- Narrative: [`OH_CLOSURE.md`](OH_CLOSURE.md).
