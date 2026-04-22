# Standard reaction enthalpies for aqueous acid–base equilibria (van’t Hoff)

This note summarizes the main technical points extracted from the literature on standard reaction enthalpies for aqueous acid–base equilibria in bioprocess models. It is **not** part of the Casagli et al. supplementary information; it documents **default** $\Delta H^\circ$ values used with the van’t Hoff form in [`SI.6`](../supporting_informations/SI.6%20Explicit%20chemical%20equilibria%2C%20their%20dissociation%20constants%20with%20temperature%20dependence.md) and implemented in [`src/bioprocess_twin/models/chemistry.py`](../../src/bioprocess_twin/models/chemistry.py).

---

## 1. Role in the model

Bioprocess simulators adjust $K_a$ (and $K_w$) with temperature using the integrated van’t Hoff equation (SI.6.4 shape):

$$
\ln\left(\frac{K_T}{K_{T_{\mathrm{ref}}}}\right) = \frac{\Delta H^\circ}{R} \left(\frac{1}{T_{\mathrm{ref}}} - \frac{1}{T}\right).
$$

Each equilibrium has its own $\Delta H^\circ$. Literature on bioprocess thermodynamics argues that treating $\Delta H^\circ$ as **constant** between roughly 273 K and 313 K is **adequate** for engineering bioreactors: heat-capacity corrections typically move $\log_{10} K$ by less than ~0.05 over ~40 K, smaller than typical kinetic and transport uncertainties.

---

## 2. Standard-state conventions (summary)

- **Conditions:** standard pressure **1 bar**, temperature **298.15 K** for the recommended $\Delta H^\circ$ (zero for nitric acid), **infinite dilution** (hypothetical ideal 1 mol·kg⁻¹ aqueous standard state extrapolated from infinite dilution — IUPAC / CODATA-style).
- **Scale:** evaluated databases usually report thermodynamics on a **molality** basis; ADM-style codes often track **molarity**. For dilute water, molality ≈ molarity; at higher $T$ or ionic strength, **density** and **activity coefficients** (e.g. Davies, extended Debye–Hückel) matter for full consistency.
- **Strong acids:** at infinite dilution, **nitric acid** is treated as fully dissociated; **$\Delta H^\circ$ for the formal dissociation equilibrium is 0 J·mol⁻¹** (enthalpy of formation of $\mathrm{HNO_3(aq)}$ aligned with $\mathrm{NO_3^-}$ + $\mathrm{H^+}$ in ATcT-style networks).

---

## 3. Recommended $\Delta H^\circ$ (J·mol⁻¹)

Reactions are written in the **acid dissociation** sense consistent with $K_a = [\text{products}]/[\text{acid}]$ (same convention as `chemistry.py`).

| Equilibrium (concept) | $\Delta H^\circ$ (J·mol⁻¹) | Notes |
|----------------------|------------------------------|--------|
| Lumped $\mathrm{CO_2(aq)/H_2CO_3^*}$ first dissociation → $\mathrm{HCO_3^-}$ + $\mathrm{H^+}$ | **9 155** | “Macro” convention used in biogeochemical / ADM-style codes (~9.1 kJ·mol⁻¹; Plummer & Busenberg ~9.109 kJ·mol⁻¹). |
| $\mathrm{HCO_3^-}$ → $\mathrm{CO_3^{2-}}$ + $\mathrm{H^+}$ | **14 700** | From CODATA-style formation enthalpies (Hess). |
| $\mathrm{NH_4^+}$ → $\mathrm{NH_3}$ + $\mathrm{H^+}$ | **52 201** | ATcT / Bates-style consensus; strongly endothermic → more free $\mathrm{NH_3}$ at higher $T$ at fixed pH. |
| $\mathrm{HNO_2}$ → $\mathrm{NO_2^-}$ + $\mathrm{H^+}$ | **11 400** | NEA-TDB / Grenthe & Fuger lineage; **larger uncertainty** (unstable species, kinetics). |
| $\mathrm{HNO_3}$ → $\mathrm{NO_3^-}$ + $\mathrm{H^+}$ | **0** | Strong-acid / infinite-dilution convention. |
| $\mathrm{H_3PO_4}$ → $\mathrm{H_2PO_4^-}$ + $\mathrm{H^+}$ | **−8 480** | Exothermic first step; CODATA-compatible / NEA-TDB–style network (historical NBS vs CODATA phosphate offsets are widely discussed in the literature). |
| $\mathrm{H_2PO_4^-}$ → $\mathrm{HPO_4^{2-}}$ + $\mathrm{H^+}$ | **3 600** | Mid-range of evaluated experimental band (~3.6–4.2 kJ·mol⁻¹). |
| $\mathrm{HPO_4^{2-}}$ → $\mathrm{PO_4^{3-}}$ + $\mathrm{H^+}$ | **14 600** | Derived consistently with cumulative protonation data in the same network. |
| Water autoprotolysis $\mathrm{H_2O}$ → $\mathrm{OH^-}$ + $\mathrm{H^+}$ ($K_w$) | **55 830** | CODATA key values; drives strong $T$ dependence of $K_w$. |

These numbers match `default_dissociation_enthalpy_j_per_mol()` in `chemistry.py`.

---

## 4. Literature / database lineage (high level)

The underlying literature traces values primarily to **CODATA** Key Values, **NIST**-related compilations, **ATcT** (Active Thermochemical Tables), and **NEA-TDB** (nuclear / environmental aqueous chemistry evaluations), with specific citations to carbonate electrochemistry (Harned-type cells), phosphate evaluations (e.g. Rard, Wolery), and nitrous acid evaluations.

---

## 5. Reference temperature for $K_{a,\mathrm{ref}}$ and $\Delta H^\circ$

Recommended $\Delta H^\circ$ in the table are anchored to **298.15 K** (infinite dilution, 1 bar). In [`src/bioprocess_twin/models/chemistry.py`](../../src/bioprocess_twin/models/chemistry.py), **`T_REF_K = 298.15`** and **`default_dissociation_constants_ref_molar()`** supplies $K_a$ and $K_w$ at that same temperature (from **`MATH_MODEL.md`** §1.2.7), so van’t Hoff uses a **matched** $(K_{a,\mathrm{ref}}, T_{\mathrm{ref}}, \Delta H^\circ)$ triple for the default bundle.

If you replace $K_{a,\mathrm{ref}}$ with values at another temperature (e.g. **293.15 K** from Table SI.6.1 alone) while keeping these $\Delta H^\circ$, re-check consistency or set **`T_REF_K`** and enthalpies to match your chosen baseline.

---

## 6. Related documentation

- Pedagogical context: [`HYDROCHEMISTRY.md`](../HYDROCHEMISTRY.md), Part E.  
- Tabulated $K_a$ and SI equation: [`SI.6`](../supporting_informations/SI.6%20Explicit%20chemical%20equilibria%2C%20their%20dissociation%20constants%20with%20temperature%20dependence.md).  
- $pK_a$ baseline in this repo: [`MATH_MODEL.md`](../MATH_MODEL.md) §1.2.7.
