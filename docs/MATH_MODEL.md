# ALBA Model Specification (Single Source of Truth)

This document serves as the technical specification for the ALBA (Algae-Bacteria) model implementation. It aggregates all parameters, equations, and logic required to build the Digital Twin.

**Source:** Casagli et al. (2021). *ALBA: A comprehensive growth model to optimize algae-bacteria wastewater treatment in raceway ponds*. Water Research, 190, 116734.

---

## 1. Model Constants & Parameters

### 1.1 Stoichiometric Parameters
Values defining the elemental composition of biomass and substrates.

| Symbol | Value | Unit | Description |
| :--- | :--- | :--- | :--- |
| **Biomass Composition** | | | |
| $i_{C,BM,ALG}$ | 0.327 | $gC \cdot gCOD_{BM}^{-1}$ | Carbon fraction in Algae |
| $i_{N,BM,ALG}$ | 0.042 | $gN \cdot gCOD_{BM}^{-1}$ | Nitrogen fraction in Algae |
| $i_{P,BM,ALG}$ | 0.008 | $gP \cdot gCOD_{BM}^{-1}$ | Phosphorus fraction in Algae |
| $i_{O,BM,ALG}$ | 0.209 | $gO \cdot gCOD_{BM}^{-1}$ | Oxygen fraction in Algae |
| $i_{H,BM,ALG}$ | 0.050 | $gH \cdot gCOD_{BM}^{-1}$ | Hydrogen fraction in Algae |
| $i_{C,BM}$ | 0.36 | $gC \cdot gCOD_{BM}^{-1}$ | Carbon fraction in Bacteria (H, AOB, NOB) |
| $i_{N,BM}$ | 0.084 | $gN \cdot gCOD_{BM}^{-1}$ | Nitrogen fraction in Bacteria |
| $i_{P,BM}$ | 0.016 | $gP \cdot gCOD_{BM}^{-1}$ | Phosphorus fraction in Bacteria |
| $i_{O,BM}$ | 0.184 | $gO \cdot gCOD_{BM}^{-1}$ | Oxygen fraction in Bacteria |
| $i_{H,BM}$ | 0.043 | $gH \cdot gCOD_{BM}^{-1}$ | Hydrogen fraction in Bacteria |
| **Yields & Fractions** | | | |
| $Y_H$ | 0.63 | $gCOD_{BM} \cdot gCOD_{S}^{-1}$ | Heterotrophic yield (Aerobic) |
| $Y_{H,NO2}$ | 0.3 | $gCOD_{BM} \cdot gCOD_{S}^{-1}$ | Heterotrophic yield (Anoxic NO2) |
| $Y_{H,NO3}$ | 0.5 | $gCOD_{BM} \cdot gCOD_{S}^{-1}$ | Heterotrophic yield (Anoxic NO3) |
| $Y_{AOB}$ | 0.2 | $gCOD_{BM} \cdot gN^{-1}$ | AOB yield |
| $Y_{NOB}$ | 0.05 | $gCOD_{BM} \cdot gN^{-1}$ | NOB yield |
| $f_{Xi,ALG}$ | 0.062 | $gCOD_{Xi} \cdot gCOD_{BM}^{-1}$ | Inert fraction from Algae decay |
| $f_{XI}$ | 0.1 | $gCOD_{XI} \cdot gCOD_{BM}^{-1}$ | Inert fraction from Bacteria decay |
| $f_{SI}$ | 0.1 | $gCOD_{SI} \cdot gCOD_{BM}^{-1}$ | Soluble inert fraction from hydrolysis |

### 1.2 Kinetic Parameters
Values governing the speed of biological reactions.

| Symbol | Value | Unit | Description |
| :--- | :--- | :--- | :--- |
| **Algae ($X_{ALG}$)** | | | |
| $\mu_{max,ALG}$ | 2.5 ± 0.055 | $d^{-1}$ | Max growth rate |
| $b_{max,r,ALG}$ | 0.1 | $d^{-1}$ | Respiration rate |
| $b_{max,d,ALG}$ | 0.03 | $d^{-1}$ | Decay rate |
| $K_{N,ALG}$ | 0.1 | $gN \cdot m^{-3}$ | Half-sat Ammonium |
| $K_{NO3,ALG}$ | 0.3 | $gN \cdot m^{-3}$ | Half-sat Nitrate |
| $K_{P,ALG}$ | 0.02 | $gP \cdot m^{-3}$ | Half-sat Phosphorus |
| $K_{C,ALG}$ | 0.004 | $gC \cdot m^{-3}$ | Half-sat Inorganic Carbon |
| $K_{O,ALG}$ | 0.2 | $gO_2 \cdot m^{-3}$ | Half-sat Oxygen (Respiration) |
| $EC_{50,O2}$ | 20 | $gO_2 \cdot m^{-3}$ | Oxygen inhibition constant ($k_{DO}$) |
| $n$ | 15 | - | Hill coefficient for O2 inhibition |
| $I_{opt}$ | 300 ± 3.8 | $\mu mol \cdot m^{-2} \cdot s^{-1}$ | Optimal Irradiance |
| $\alpha_{light}$ | 0.01 | $\mu mol^{-1} \cdot m^2 \cdot s$ | Initial slope of PI curve |
| **Heterotrophs ($X_H$)** | | | |
| $\mu_{max,H}$ | 6.0 | $d^{-1}$ | Max growth rate |
| $b_{max,r,H}$ | 0.3 | $d^{-1}$ | Respiration rate |
| $b_{max,d,H}$ | 0.9 | $d^{-1}$ | Decay rate |
| $K_{S,H}$ | 4.0 | $gCOD \cdot m^{-3}$ | Half-sat Organic Substrate |
| $K_{O,H}$ | 0.2 | $gO_2 \cdot m^{-3}$ | Half-sat Oxygen |
| $K_{N,H}$ | 0.05 | $gN \cdot m^{-3}$ | Half-sat Ammonium |
| $K_{P,H}$ | 0.01 | $gP \cdot m^{-3}$ | Half-sat Phosphorus |
| $\eta_{ANOX}$ | 0.6 | - | Anoxic reduction factor |
| **Nitrifiers ($X_{AOB}, X_{NOB}$)** | | | |
| $\mu_{max,AOB}$ | 0.72 ± 0.0005 | $d^{-1}$ | Max growth rate AOB |
| $\mu_{max,NOB}$ | 0.65 ± 0.023 | $d^{-1}$ | Max growth rate NOB |
| $K_{N,AOB}$ | 0.5 | $gN \cdot m^{-3}$ | Half-sat Ammonium (AOB) |
| $K_{NO2,NOB}$ | 0.5 | $gN \cdot m^{-3}$ | Half-sat Nitrite (NOB) |
| $K_{O,AOB}$ | 0.8 | $gO_2 \cdot m^{-3}$ | Half-sat Oxygen (AOB) |
| $K_{O,NOB}$ | 2.2 | $gO_2 \cdot m^{-3}$ | Half-sat Oxygen (NOB) |

---

## 2. State Vector Definition
The model state vector $\mathbf{y}$ consists of 17 variables.

| Index | Symbol | Unit | Description |
| :--- | :--- | :--- | :--- |
| 1 | $X_{ALG}$ | $gCOD \cdot m^{-3}$ | Algal biomass |
| 2 | $X_{AOB}$ | $gCOD \cdot m^{-3}$ | Ammonia Oxidizing Bacteria |
| 3 | $X_{NOB}$ | $gCOD \cdot m^{-3}$ | Nitrite Oxidizing Bacteria |
| 4 | $X_H$ | $gCOD \cdot m^{-3}$ | Heterotrophic Bacteria |
| 5 | $X_S$ | $gCOD \cdot m^{-3}$ | Slowly biodegradable organic matter |
| 6 | $X_I$ | $gCOD \cdot m^{-3}$ | Inert particulate organic matter |
| 7 | $S_S$ | $gCOD \cdot m^{-3}$ | Readily biodegradable organic matter |
| 8 | $S_I$ | $gCOD \cdot m^{-3}$ | Inert soluble organic matter |
| 9 | $S_{IC}$ | $gC \cdot m^{-3}$ | Total Inorganic Carbon ($CO_2 + HCO_3^- + CO_3^{2-}$) |
| 10 | $S_{ND}$ | $gN \cdot m^{-3}$ | Soluble Organic Nitrogen (Urea) |
| 11 | $S_{NH}$ | $gN \cdot m^{-3}$ | Total Ammoniacal Nitrogen ($NH_3 + NH_4^+$) |
| 12 | $S_{NO2}$ | $gN \cdot m^{-3}$ | Total Nitrite Nitrogen ($HNO_2 + NO_2^-$) |
| 13 | $S_{NO3}$ | $gN \cdot m^{-3}$ | Total Nitrate Nitrogen ($HNO_3 + NO_3^-$) |
| 14 | $S_{N2}$ | $gN \cdot m^{-3}$ | Nitrogen Gas |
| 15 | $S_{PO4}$ | $gP \cdot m^{-3}$ | Total Inorganic Phosphorus |
| 16 | $S_{O2}$ | $gO_2 \cdot m^{-3}$ | Dissolved Oxygen |
| 17 | $S_{H2O}$ | $gH \cdot m^{-3}$ | Water (Balance term) |

---

## 3. Algebraic Sub-models

### 3.1 Temperature Dependence ($f_T$)
*   **Growth:** Uses CTMI Model (Rosso et al., 1993).
    $$ f_T(T) = \frac{(T-T_{max})(T-T_{min})^2}{(T_{opt}-T_{min})[(T_{opt}-T_{min})(T-T_{opt}) - (T_{opt}-T_{max})(T_{opt}+T_{min}-2T)]} $$
*   **Decay:** Uses Arrhenius.
    $$ f_{T,decay}(T) = \theta^{(T-20)} $$

### 3.2 pH Dependence ($f_{pH}$)
Uses Cardinal pH Model (CPM).
$$ f_{pH}(pH) = \frac{(pH - pH_{min})(pH - pH_{max})}{(pH - pH_{min})(pH - pH_{max}) - (pH - pH_{opt})^2} $$

### 3.3 Light Dependence ($f_I$)
For Algae only. Uses Haldane model (Bernard & Rémond, 2012).
$$ f_I(I) = \frac{\mu_{max}}{1 + \frac{\mu_{max}}{\alpha} \left(\frac{I}{I_{opt}} - 1\right)^2} $$

### 3.4 Oxygen Inhibition ($f_{DO}$)
For Algae.
*   **Growth Inhibition:** $f_{DO,g} = \frac{EC_{50,O2}^n}{S_{O2}^n + EC_{50,O2}^n}$
*   **Decay Enhancement:** $f_{DO,d} = \frac{S_{O2}^n}{S_{O2}^n + EC_{50,O2}^n}$

---

## 4. Process Rate Equations ($\rho$)
The system defines 19 biological processes using **Liebig's Law of the Minimum** for nutrient limitation.

### Algae ($X_{ALG}$)
1.  **Growth on $NH_4$:**
    $$ \rho_1 = \mu_{max,ALG} \cdot f_I \cdot f_T \cdot f_{pH} \cdot f_{DO,g} \cdot \min\left(\frac{S_{IC}}{K_{C} + S_{IC}}, \frac{S_{NH}}{K_{N} + S_{NH}}, \frac{S_{PO4}}{K_{P} + S_{PO4}}\right) \cdot X_{ALG} $$
2.  **Growth on $NO_3$:**
    $$ \rho_2 = \mu_{max,ALG} \cdot f_I \cdot f_T \cdot f_{pH} \cdot f_{DO,g} \cdot \frac{K_{N}}{K_{N} + S_{NH}} \cdot \min\left(\frac{S_{IC}}{K_{C} + S_{IC}}, \frac{S_{NO3}}{K_{NO3} + S_{NO3}}, \frac{S_{PO4}}{K_{P} + S_{PO4}}\right) \cdot X_{ALG} $$
3.  **Respiration:**
    $$ \rho_3 = b_{max,r,ALG} \cdot f_T \cdot f_{pH} \cdot \frac{S_{O2}}{K_{O} + S_{O2}} \cdot X_{ALG} $$
4.  **Decay:**
    $$ \rho_4 = b_{max,d,ALG} \cdot \theta^{(T-20)} \cdot f_{pH} \cdot f_{DO,d} \cdot X_{ALG} $$

### Heterotrophs ($X_H$)
5.  **Aerobic Growth ($NH_4$):**
    $$ \rho_5 = \mu_{max,H} \cdot f_T \cdot f_{pH} \cdot \min\left(\frac{S_S}{K_S + S_S}, \frac{S_{O2}}{K_O + S_{O2}}, \frac{S_{NH}}{K_N + S_{NH}}, \frac{S_{PO4}}{K_P + S_{PO4}}\right) \cdot X_H $$
6.  **Aerobic Growth ($NO_3$):**
    $$ \rho_6 = \mu_{max,H} \cdot f_T \cdot f_{pH} \cdot \frac{K_N}{K_N + S_{NH}} \cdot \min\left(\frac{S_S}{K_S + S_S}, \frac{S_{O2}}{K_O + S_{O2}}, \frac{S_{NO3}}{K_{NO3} + S_{NO3}}, \frac{S_{PO4}}{K_P + S_{PO4}}\right) \cdot X_H $$
7.  **Aerobic Respiration:**
    $$ \rho_7 = b_{max,r,H} \cdot f_T \cdot f_{pH} \cdot \frac{S_{O2}}{K_O + S_{O2}} \cdot X_H $$
8.  **Anoxic Growth ($NO_2$):**
    $$ \rho_8 = \mu_{max,H} \cdot \eta_{ANOX} \cdot f_T \cdot f_{pH} \cdot \frac{K_O}{K_O + S_{O2}} \cdot \min\left(\frac{S_S}{K_S + S_S}, \frac{S_{NO2}}{K_{NO2} + S_{NO2}}, \frac{S_{PO4}}{K_P + S_{PO4}}\right) \cdot X_H $$
9.  **Anoxic Growth ($NO_3$):**
    $$ \rho_9 = \mu_{max,H} \cdot \eta_{ANOX} \cdot f_T \cdot f_{pH} \cdot \frac{K_O}{K_O + S_{O2}} \cdot \min\left(\frac{S_S}{K_S + S_S}, \frac{S_{NO3}}{K_{NO3} + S_{NO3}}, \frac{S_{PO4}}{K_P + S_{PO4}}\right) \cdot X_H $$
10. **Anoxic Respiration:**
    $$ \rho_{10} = b_{max,r,H} \cdot \eta_{ANOX} \cdot f_T \cdot f_{pH} \cdot \frac{K_O}{K_O + S_{O2}} \cdot \min\left(\frac{S_{NO2}}{K_{NO2} + S_{NO2}}, \frac{S_{NO3}}{K_{NO3} + S_{NO3}}\right) \cdot X_H $$

### Hydrolysis & Decay
11. **Hydrolysis ($X_S \to S_S$):**
    $$ \rho_{11} = \mu_{Hyd} \cdot \theta_{HYD}^{(T-20)} \cdot f_{pH} \cdot \frac{X_S/X_H}{K_{HYD} + (X_S/X_H)} \cdot X_H $$
12. **Ammonification (Urea $\to NH_4$):**
    $$ \rho_{12} = \mu_{a} \cdot \theta_{AMM}^{(T-20)} \cdot f_{pH} \cdot \frac{S_{ND}}{K_a + S_{ND}} \cdot X_H $$
13. **Decay ($X_H$):**
    $$ \rho_{13} = b_{max,d,H} \cdot \theta_H^{(T-20)} \cdot f_{pH} \cdot X_H $$

### Nitrifiers ($X_{AOB}, X_{NOB}$)
14. **AOB Growth:**
    $$ \rho_{14} = \mu_{max,AOB} \cdot f_T \cdot f_{pH} \cdot \min\left(\frac{S_{NH}}{K_N + S_{NH}}, \frac{S_{O2}}{K_O + S_{O2}}, \frac{S_{IC}}{K_C + S_{IC}}, \frac{S_{PO4}}{K_P + S_{PO4}}\right) \cdot X_{AOB} $$
15. **AOB Respiration:**
    $$ \rho_{15} = b_{max,r,AOB} \cdot f_T \cdot f_{pH} \cdot \frac{S_{O2}}{K_O + S_{O2}} \cdot X_{AOB} $$
16. **AOB Decay:**
    $$ \rho_{16} = b_{max,d,AOB} \cdot \theta_{AOB}^{(T-20)} \cdot f_{pH} \cdot X_{AOB} $$
17. **NOB Growth:**
    $$ \rho_{17} = \mu_{max,NOB} \cdot f_T \cdot f_{pH} \cdot \min\left(\frac{S_{NO2}}{K_{NO2} + S_{NO2}}, \frac{S_{O2}}{K_O + S_{O2}}, \frac{S_{IC}}{K_C + S_{IC}}, \frac{S_{PO4}}{K_P + S_{PO4}}\right) \cdot X_{NOB} $$
18. **NOB Respiration:**
    $$ \rho_{18} = b_{max,r,NOB} \cdot f_T \cdot f_{pH} \cdot \frac{S_{O2}}{K_O + S_{O2}} \cdot X_{NOB} $$
19. **NOB Decay:**
    $$ \rho_{19} = b_{max,d,NOB} \cdot \theta_{NOB}^{(T-20)} \cdot f_{pH} \cdot X_{NOB} $$

---

## 5. Stoichiometric Matrix Construction Rules

The stoichiometric coefficients ($\alpha_{i,j}$) for the matrix $\mathbf{S}$ ($19 \times 17$) are calculated dynamically based on mass and charge balances.

**Example Rule (Algal Growth on $NH_4$, $\rho_1$):**
*   $X_{ALG}$ (Production): $+1$
*   $S_{NH}$ (Consumption): $-i_{N,BM,ALG}$
*   $S_{PO4}$ (Consumption): $-i_{P,BM,ALG}$
*   $S_{IC}$ (Consumption): $-i_{C,BM,ALG}$
*   $S_{O2}$ (Production): $+ \frac{32}{12}i_{C,BM} + \frac{24}{14}i_{N,BM} + \frac{40}{31}i_{P,BM} + 8i_{H,BM} - i_{O,BM}$

**General Logic:**
For any growth process, the coefficient for $S_{O2}$ is derived from the COD balance:
$$ \alpha_{O2} = 1 - Y_{biomass} $$
(Adjusted for units of $gO_2$ vs $gCOD$).

> For the complete Petersen matrix of the ALBA model, see [STOICHIOMETRY.md](STOICHIOMETRY.md).