# Rationale: auxiliary water column as stoichiometry, not a free tuning knob

State **`S_H2O`** in the Petersen matrix (see `MATH_MODEL.md`).

## Units

In ALBA / `MATH_MODEL.md`, **`S_H2O`** is a **hydrogen mass** inventory. In display form:

$$
\text{units of } S_{\mathrm{H2O}}:\ \mathrm{g\,H\cdot m^{-3}}.
$$

The composition row for that column gives:

$$
I_{\mathrm{H},S_{\mathrm{H2O}}}=1,
\qquad
I_{\mathrm{O},S_{\mathrm{H2O}}}=7.94
$$

(oxygen mass paired with that hydrogen in water: about **8 g O per g H**).

## What “closing oxygen” means

For process row *i*, the SI lists stoichiometric coefficients for every state **except** that some **atomic oxygen** appearing in real chemistry as **water** was left implicit (ASM-style).

Let

$$
L_i = \sum_{j \neq S_{\mathrm{H2O}}} S_{i,j}\, I_{\mathrm{O},j}
$$

be the **elemental oxygen** contributed by all **non-water** columns, using the same composition matrix **I** as the code. That sum is a **linear functional of the SI row**: it is fixed once the reaction is fixed.

If **all** oxygen that is still missing from the row must appear as oxygen bound in water, and water is tracked only through **`S_H2O`** with the composition above, then the **unique** total coefficient `α_i` on **`S_H2O`** that satisfies

$$
L_i + \alpha_i\, I_{\mathrm{O},S_{\mathrm{H2O}}} = 0
$$

is

$$
\alpha_i = -\frac{L_i}{I_{\mathrm{O},S_{\mathrm{H2O}}}}.
$$

So `α_i` is **not** chosen to “make the test pass” independently of chemistry: it is the **only** value consistent with (i) the SI stoichiometry for all explicit species and (ii) the modelling convention that **residual elemental O is water oxygen** in this column.

Algebraically, define the SI oxygen balance and the SI water coefficient:

$$
B_{i,\mathrm{O}}^{\mathrm{SI}}=\sum_j S_{i,j}^{\mathrm{SI}} I_{\mathrm{O},j},
\qquad
\text{SI coefficient on water: }\alpha_i^{\mathrm{SI}}.
$$

Then

$$
L_i = B_{i,\mathrm{O}}^{\mathrm{SI}} - \alpha_i^{\mathrm{SI}} I_{\mathrm{O},S_{\mathrm{H2O}}},
\qquad
\alpha_i = \alpha_i^{\mathrm{SI}} - \frac{B_{i,\mathrm{O}}^{\mathrm{SI}}}{I_{\mathrm{O},S_{\mathrm{H2O}}}}.
$$

The second term is what was previously called “Δ”. The **recommended reading** is: **total water-hydrogen production per process unit** satisfies

$$
\alpha_i = -\frac{L_i}{I_{\mathrm{O},S_{\mathrm{H2O}}}}.
$$

## Relation to ρ₁₄ (AOB)

The worksheet [`balance_of_oxygen_for_rho_14_Aerobic_growth_of_X_AOB_on_NH4+.md`](balance_of_oxygen_for_rho_14_Aerobic_growth_of_X_AOB_on_NH4+.md) computes, for biomass yield in g COD, approximately

$$
B_{14,\mathrm{O}}^{\mathrm{SI}}\approx -5.55\ \mathrm{g\,O}
$$

with **no** **`S_H2O`** term in the SI row. Interpreting that deficit as oxygen that must appear in water gives a **water-derived H mass** of order

$$
\frac{|B_{14,\mathrm{O}}^{\mathrm{SI}}|}{7.94} \approx 0.70\ \mathrm{g\,H}
$$

per g COD biomass — the same magnitude as `α_i` above.

## What this does **not** fix

One column cannot generally enforce **both** oxygen and hydrogen closure plus **charge** unless the true reaction also closes with water only. In display form, “O and H closed” means:

$$
B_{i,\mathrm{O}}=0,
\qquad
B_{i,\mathrm{H}}=0.
$$

Remaining **hydrogen** imbalance usually means an explicit **acid–base** species (e.g. **H⁺**) is still omitted (ADR 007).

## Code

`compute_stoichiometric_s_h2o_total_for_row` and `build_petersen_matrix_with_oh_closure` in
[`stoichiometry_closure.py`](../../src/bioprocess_twin/models/stoichiometry_closure.py).
