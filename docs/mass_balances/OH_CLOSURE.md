# Oxygen / hydrogen closure (experimental Petersen copy)

This note accompanies **[ADR 007](../adrs/007-elemental-mass-balance-oh-closure.md)** and the SI-faithful audit in [`MASS_BALANCES.md`](../MASS_BALANCES.md).

## Strategy

- **`get_petersen_matrix()`** remains the **literal Casagli SI.3.3** transcription.
- **`build_petersen_matrix_with_oh_closure()`** in
  [`src/bioprocess_twin/models/stoichiometry_closure.py`](../../src/bioprocess_twin/models/stoichiometry_closure.py)
  returns a **NumPy copy** of that matrix. Only column **`S_H2O`** may differ from the SI.
- For each process row $i$ with $|B_{i,\mathrm{O}}| > \texttt{atol}$ or
  $|B_{i,\mathrm{H}}| > \texttt{atol}$ (using `MASS_BALANCE_ATOL` in
  `stoichiometry_validation.py`), we set the **total** Petersen coefficient on
  $S_{\mathrm{H2O}}$ to the value implied by **elemental oxygen** once all
  non-water species are fixed by the SI row:

  $$L_i = \sum_{j \neq S_{\mathrm{H2O}}} S_{i,j}^{\mathrm{SI}}\, I_{\mathrm{O},j},
  \qquad
  \alpha_i = -\frac{L_i}{I_{\mathrm{O},S_{\mathrm{H2O}}}}.$$

  So $\alpha_i$ is the net **g H** (per process normalization) carried by water in
  this convention when **all “missing” oxygen** in the ASM-style row is assigned to
  water oxygen—not an independent tuning parameter. Algebraically,
  $\alpha_i = S_{i,S_{\mathrm{H2O}}}^{\mathrm{SI}} - B_{i,\mathrm{O}}^{\mathrm{SI}} / I_{\mathrm{O},S_{\mathrm{H2O}}}$
  (same as the former “$\Delta$ on top of SI” recipe).

  Default strategy name: **`stoichiometric_water`**. **`oxygen_exact`** is kept as an
  alias (identical numeric result). See
  [`stoichiometric_water_rationale.md`](stoichiometric_water_rationale.md).

- **`least_squares`** (optional) minimizes
  $(B_{i,\mathrm{O}} + \Delta I_{\mathrm{O}})^2 + (B_{i,\mathrm{H}} + \Delta I_{\mathrm{H}})^2$
  with a single $\Delta$; it can worsen one element while improving the sum of squares.

Hydrogen often **does not** reach zero with a single $S_{\mathrm{H2O}}$ degree of freedom;
full closure may require an explicit acid–base carrier (e.g. $\mathrm{H}^+$) and charge
bookkeeping (future work).

## Artefacts

- **[`MASS_BALANCES.md`](../MASS_BALANCES.md):** full 114-cell audit for **SI** $\mathbf{S}$.
- **[`MASS_BALANCES_CLOSURE_OF_OXYGEN.md`](../MASS_BALANCES_CLOSURE_OF_OXYGEN.md):** same 114-cell audit for the **oxygen (water) closure** $\mathbf{S}$.
- **`list_oh_mass_balance_violations()`:** programmatic inventory of SI O/H failures.

Regenerate Markdown:

```bash
uv run python scripts/generate_mass_balances_md.py
uv run python scripts/generate_mass_balances_md.py --closure
```

## Simulation integration (later)

`get_petersen_matrix_for_simulation(use_oh_closure=False)` returns the SI matrix by default.
Setting `use_oh_closure=True` selects the closure copy **without** changing the rest of the
ODE stack. Before enabling this in production runs, validate the effect on $S_{\mathrm{H2O}}$
dynamics and consistency with pH / gas–liquid submodels from the ALBA SI.
