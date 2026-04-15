# ADR 007: Elemental Mass-Balance Closure for Oxygen and Hydrogen (Water and Protons)

**Status:** Accepted  
**Date:** 2026-04-14  
**Related:** Sprint 2.5A (`phase1-02.5A`), `docs/mass_balances/`, `docs/STOICHIOMETRY.md`, Casagli et al. (2021) SI.3.1–SI.3.3

## Context

The ALBA stoichiometric backbone uses a **Petersen matrix** $\mathbf{S}$ and a **composition matrix** $\mathbf{I}$. A standard consistency check is that, for each biological process $i$ and each conserved quantity $k$ represented as a row of $\mathbf{I}$,

$$
\sum_j S_{i,j}\, I_{k,j} \approx 0
$$

(equivalently $\mathbf{S}\mathbf{I}^\top \approx \mathbf{0}$ in the layout used in this repository).

After aligning $\mathbf{S}$ and $\mathbf{I}$ with the Casagli et al. (2021) supplementary information, **COD**, **C**, **N**, and **P** are largely consistent with the intended ASM-style formulation, but **atomic oxygen (O)** and **atomic hydrogen (H)** exhibit **large residuals** on several **bacterial** processes (e.g. heterotrophic growth and respiration, AOB/NOB growth and respiration, urea hydrolysis). The worst-known case is aerobic growth of AOB on ammonium ($\rho_{14}$), with an oxygen residual of order **$-5.55$** (mass-based, in the units of the test), documented in `docs/mass_balances/balance_of_oxygen_for_rho_14_Aerobic_growth_of_X_AOB_on_NH4+.md` and internal synthesis.

The Casagli SI explicitly includes **$S_{\mathrm{H2O}}$** coefficients for **algal** processes in SI.3.3, while many **bacterial** rows follow **Activated Sludge Model (ASM)** convention: **water** is treated as an infinite solvent and **acid–base** (free protons, speciation) is handled outside the strict elemental bookkeeping of $\mathbf{S}$ for those rows. That hybrid is adequate for **electron (COD)** and **nutrient** tracking as in classical wastewater engineering, but it **does not** guarantee closure of **elemental O and H** under a strict $\mathbf{S}\mathbf{I}^\top$ audit unless **closure species** (e.g. $S_{\mathrm{H2O}}$, and where needed $\mathrm{H}^+$ or alkalinity) are added consistently.

The project must choose a **policy** for verification and for the **digital twin** roadmap: either extend the documented stoichiometry toward **RWQM-style** elemental rigor, or explicitly **scope out** strict O/H checks.

## Options Considered

### Option (a): Extend the model (RWQM-style elemental closure)

**Description:** Keep **SI-tabulated** coefficients for all columns that the SI defines. Add (or complete) entries on affected rows for **$S_{\mathrm{H2O}}$** and, where required for simultaneous O, H, and **charge** consistency, **acid–base carriers** (e.g. $\mathrm{H}^+$ or an alkalinity surrogate). Coefficients not given in the SI for bacterial processes are **derived** from elemental (and charge) balances and documented as a **project extension** traceable to first principles and peer references (ASM, RWQM1, elemental-balance methodology).

**Advantages (+)**

- Satisfies **strict elemental** mass balance for **O** and **H** in $\mathbf{S}\mathbf{I}^\top$, aligning the twin with **conservation-of-matter** expectations for mechanistic auditing.
- Supports future coupling to **pH**, **speciation**, or analyses that depend on consistent **atomic** bookkeeping, without known “missing oxygen” artifacts in $\mathbf{S}$.
- Allows tightening **test tolerances** (e.g. toward `1e-6`) for O/H once coefficients are stable.

**Disadvantages (-)**

- **Engineering effort:** every affected process row must be revisited; **charge** must be kept coherent with ions already in $\mathbf{I}$ ($\mathrm{NH_4^+}$, $\mathrm{NO_2^-}$, $\mathrm{PO_4^{3-}}$, etc.).
- **Extension vs publication:** new coefficients are **not** literal copies from Casagli tables for those rows; they must be clearly labeled as **twin-layer extensions** with derivation notes.
- **Dynamics:** if new species participate in $\mathrm{d}\mathbf{C}/\mathrm{d}t$, they must remain consistent with the **pH / equilibrium** and **gas–liquid** submodels described elsewhere in the SI (kinetics and transfer do not, by themselves, fix missing rows in $\mathbf{S}$).

### Option (b): Document ASM-style scope and narrow tests

**Description:** Declare that **strict elemental O/H closure** in $\mathbf{S}\mathbf{I}^\top$ is **out of scope** for the representation “as published” in the bacterial block. Document that **COD** (and selected elements) are the intended conservation checks. Relax or partition **unit tests** (e.g. exclude O/H for ASM-style rows or raise `atol` only for those cells).

**Advantages (+)**

- **Minimal change** and **maximum textual fidelity** to the SI tables for bacterial $\alpha$ entries.
- **Fast** path to a twin that reproduces **published** ALBA behavior for design variables (DO, nutrients, biomass) without re-deriving stoichiometry.

**Disadvantages (-)**

- The **digital twin** cannot claim **full elemental** consistency in $\mathbf{S}$ for O/H without qualification.
- Large **O/H residuals** remain a permanent feature of strict audits, which is weak for **long-term credibility** if the product sells **elemental** or **mass-balance-first** rigor.

## Decision

For **ALBATwinHub**, we adopt **Option (a)**.

We will **extend** the stoichiometric model so that **elemental O and H** (and charge, where required) **close** under the same $\mathbf{S}\mathbf{I}^\top$ convention used in tests, while **preserving SI-given** coefficients for all columns that the SI specifies. Any additional Petersen entries will be **derived, documented, and tested**, and treated as an explicit **extension layer** for a **digital twin** that prioritizes **conservation of matter** in the tabulated biological processes, not only COD-style electron accounting.

Option **(b)** is **rejected** as the long-term stance for this product, though it remains the correct historical reading of why the **published** SI alone does not close O/H on every bacterial row.

## Implementation Notes (Non-Exhaustive)

- Primary code: `src/bioprocess_twin/models/stoichiometry.py` (`get_petersen_matrix`, `get_composition_matrix` if new columns or rows are required).
- Tests: `tests/unit/stoichiometry_mass_balance_shared.py` (`MASS_BALANCE_ATOL`), `tests/unit/test_stoichiometry_mass_balance_cells.py`, `tests/unit/test_stoichiometry.py`.
- Documentation: update `docs/STOICHIOMETRY.md` and `docs/mass_balances/` as coefficients are fixed; keep cross-links to Casagli SI and this ADR.
- Explicit $\mathbf{S}\mathbf{I}^\top$ breakdown for all 114 process×element cells: [`docs/MASS_BALANCES.md`](../MASS_BALANCES.md) (regenerate with `uv run python scripts/generate_mass_balances_md.py`).
- Experimental **O-exact** $S_{\mathrm{H2O}}$ closure copy: [`stoichiometry_closure.py`](../../src/bioprocess_twin/models/stoichiometry_closure.py), narrative [`OH_CLOSURE.md`](../mass_balances/OH_CLOSURE.md), audit [`MASS_BALANCES_CLOSURE_OF_OXYGEN.md`](../MASS_BALANCES_CLOSURE_OF_OXYGEN.md) (`uv run python scripts/generate_mass_balances_md.py --closure`).

## Consequences

- (+) The twin’s stoichiometric core can meet **strict O/H (and charge-aware) closure** targets once implementation is complete.
- (+) Clear **separation of concerns**: “SI transcription” vs **twin extension** for closure species.
- (-) Higher **upfront and review** cost before O/H tests can be tightened to strict tolerances.
- (-) **Simulation and calibration** workflows may need revisiting if new terms materially affect **$S_{\mathrm{H2O}}$** or **pH-related** state dynamics.
