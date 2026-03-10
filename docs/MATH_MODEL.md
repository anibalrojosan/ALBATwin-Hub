# Mathematical Model - ALBA (Algae-Bacteria)

## Table of Contents
1. [State Variables](#state-variables)
2. [Biological Kinetics](#biological-kinetics)
3. [Stoichiometry](#stoichiometry)
4. [pH and Chemical Equilibria](#ph-and-chemical-equilibria)
5. [Gas-Liquid Transfer](#gas-liquid-transfer)
6. [Temperature Dependence](#temperature-dependence)

---

## State Variables

The model tracks 17 state variables ($C_i$), typically expressed in $g \cdot m^{-3}$ (equivalent to $mg \cdot L^{-1}$).

| Symbol | Description | Unit |
| :--- | :--- | :--- |
| **Biomass** | | |
| $X_{ALG}$ | Algal biomass | $gCOD \cdot m^{-3}$ |
| $X_H$ | Heterotrophic bacteria | $gCOD \cdot m^{-3}$ |
| $X_{AOB}$ | Ammonia Oxidizing Bacteria | $gCOD \cdot m^{-3}$ |
| $X_{NOB}$ | Nitrite Oxidizing Bacteria | $gCOD \cdot m^{-3}$ |
| $X_I$ | Inert particulate organic matter | $gCOD \cdot m^{-3}$ |
| **Soluble Substrates** | | |
| $S_S$ | Readily biodegradable organic matter | $gCOD \cdot m^{-3}$ |
| $S_I$ | Inert soluble organic matter | $gCOD \cdot m^{-3}$ |
| $S_{NH}$ | Total Ammoniacal Nitrogen ($NH_4^+ + NH_3$) | $gN \cdot m^{-3}$ |
| $S_{NO2}$ | Nitrite Nitrogen | $gN \cdot m^{-3}$ |
| $S_{NO3}$ | Nitrate Nitrogen | $gN \cdot m^{-3}$ |
| $S_{PO4}$ | Total Inorganic Phosphorus | $gP \cdot m^{-3}$ |
| $S_{O2}$ | Dissolved Oxygen | $gO_2 \cdot m^{-3}$ |
| $S_{IC}$ | Total Inorganic Carbon ($CO_2 + HCO_3^- + CO_3^{2-}$) | $gC \cdot m^{-3}$ |
| $S_{ND}$ | Soluble Biodegradable Organic Nitrogen (Urea) | $gN \cdot m^{-3}$ |
| **Particulate Substrates** | | |
| $X_S$ | Slowly biodegradable organic matter | $gCOD \cdot m^{-3}$ |

---

## Biological Kinetics

The model uses **Liebig's Law of the Minimum** for nutrient limitation, rather than multiplicative Monod terms.

### General Growth Rate Equation
For a biomass $X_i$, the growth rate $\rho_{growth}$ is defined as:

$$ \rho_{growth} = \mu_{max} \cdot f(T) \cdot f(pH) \cdot f(I) \cdot \min \left( \frac{S_1}{K_1 + S_1}, \frac{S_2}{K_2 + S_2}, \dots \right) \cdot X_i $$

Where:
*   $\mu_{max}$: Maximum specific growth rate ($d^{-1}$)
*   $f(T)$: Temperature correction factor
*   $f(pH)$: pH correction factor
*   $f(I)$: Light correction factor (for algae only)
*   $S_j$: Limiting substrate concentration
*   $K_j$: Half-saturation constant

### 1. Algal Kinetics (Processes $\rho_1, \rho_2$)
**Light Dependence ($f_I$):**
Uses the Haldane model modified by Bernard & Remond (2012), integrated over depth $z$:

$$ f_I(I) = \frac{I}{I + \frac{\mu_{max}}{\alpha} \left( \frac{I}{I_{opt}} - 1 \right)^2} $$

The average growth rate in the reactor is obtained by integrating $f_I(I(z))$ where $I(z)$ follows Beer-Lambert law:
$$ I(z) = I_0 \cdot e^{-\epsilon \cdot TSS \cdot z} $$

**Oxygen Inhibition:**
Photosynthesis is inhibited by high $DO$ concentrations (Hill type):
$$ f_{O2, inh} = \frac{K_{O2, inh}^n}{S_{O2}^n + K_{O2, inh}^n} $$

### 2. Bacterial Kinetics (Processes $\rho_5 - \rho_{10}$)
Heterotrophs can grow aerobically ($\rho_5, \rho_6$) or anoxically ($\rho_8, \rho_9$).
*   **Aerobic:** Depends on $\frac{S_{O2}}{K_{O2} + S_{O2}}$
*   **Anoxic:** Depends on $\frac{K_{O2}}{K_{O2} + S_{O2}}$ (inhibition by oxygen) and presence of $NO_x$.

### 3. Nitrification (Processes $\rho_{14}, \rho_{17}$)
Two-step nitrification ($NH_4^+ \to NO_2^- \to NO_3^-$).
*   **AOB:** Growth on $S_{NH}$. Sensitive to $S_{O2}$ and $S_{IC}$.
*   **NOB:** Growth on $S_{NO2}$. Sensitive to $S_{O2}$ and $S_{IC}$.

---

## Stoichiometry

The model considers 19 biological processes. The mass balance is ensured via a stoichiometric matrix $\mathbf{S}$.

$$ \frac{d\mathbf{C}}{dt} = \mathbf{S}^T \cdot \boldsymbol{\rho} + \text{Transport} + \text{Transfer} $$

Key stoichiometric coefficients ($\alpha$) are derived from the elemental composition of biomass (e.g., $C_{100}H_{183}O_{48}N_{11}P$ for algae).

Example for Algal Growth on Ammonium ($\rho_1$):
*   Consumes: $S_{NH}$, $S_{PO4}$, $S_{IC}$
*   Produces: $X_{ALG}$, $S_{O2}$

---

## pH and Chemical Equilibria

The pH is calculated by solving the charge balance equation at each time step. The system considers the following equilibria:

1.  **Ammonia:** $NH_4^+ \leftrightarrow NH_3 + H^+$
2.  **Carbonate:** $CO_2 + H_2O \leftrightarrow HCO_3^- + H^+ \leftrightarrow CO_3^{2-} + 2H^+$
3.  **Phosphate:** $H_3PO_4 \leftrightarrow H_2PO_4^- \leftrightarrow HPO_4^{2-} \leftrightarrow PO_4^{3-}$
4.  **Nitrite/Nitrate:** $HNO_2 \leftrightarrow NO_2^-$, $HNO_3 \leftrightarrow NO_3^-$
5.  **Water:** $H_2O \leftrightarrow H^+ + OH^-$

**Charge Balance Equation:**
$$ [H^+] + [NH_4^+] + \dots - [OH^-] - [HCO_3^-] - 2[CO_3^{2-}] - \dots + Z_{net} = 0 $$

---

## Gas-Liquid Transfer

Gas transfer rates ($Q_j$) for $O_2$, $CO_2$, and $NH_3$ are modeled as:

$$ Q_j = k_L a_j \cdot (S_{j, sat} - S_j) $$

*   **Oxygen:** $k_L a_{O2}$ is a calibrated parameter (e.g., $34 d^{-1}$).
*   **Others:** Derived from diffusivity ratios:
    $$ k_L a_j = k_L a_{O2} \cdot \sqrt{\frac{D_j}{D_{O2}}} $$
*   **Saturation ($S_{sat}$):** Calculated using Henry's Law, dependent on temperature and partial pressure ($p_{gas}$).

---

## Temperature Dependence

### Growth Rates (CTMI Model)
The Cardinal Temperature Model with Inflection (Rosso et al., 1993) is used for all growth rates:

$$ f(T) = \frac{(T - T_{max})(T - T_{min})^2}{(T_{opt} - T_{min})[(T_{opt} - T_{min})(T - T_{opt}) - (T_{opt} - T_{max})(T_{opt} + T_{min} - 2T)]} $$

Defined by three parameters: $T_{min}$, $T_{opt}$, $T_{max}$.

### Decay and Hydrolysis (Arrhenius)
Simple exponential dependence:
$$ k(T) = k_{20} \cdot \theta^{(T - 20)} $$

---
*Last updated: March 10, 2026 by Anibal Rojo*
