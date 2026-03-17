### ALBA model stoichiometric matrix

### 1. Petersen Matrix

This matrix shows the **dynamic interaction** between the processes and the variables. 
*   **Rows:** The 19 biological processes (e.g. "Phototrophic growth on NH₄⁺").
*   **Columns:** The 17 state variables (e.g. "Ammonium", "Oxygen").
*   **Cells ($\alpha_{i,j}$):** Indicate if a variable is consumed (negative value), produced (positive value) or not affected (empty/zero cell) by that specific process.

They are used to calculate the **differential equations (ODEs)**. The variation of a state variable in time ($dC/dt$) is the sum of the processes that affect it:
$$\frac{d\mathbf{C}}{dt} = \mathbf{S}^T \cdot \boldsymbol{\rho}$$
Where $\mathbf{S}$ is this matrix and $\boldsymbol{\rho}$ are the reaction rates.

| component j→ / process i↓ | X_ALG | X_AOB | X_NOB | X_H  | X_S  | X_I  | S_S   | S_I   | S_IC  | S_ND  | S_NH  | S_NO2 | S_NO3 | S_N2  | S_PO4 | S_O2  | S_H2O |
|--------------------------|-------|-------|-------|------|------|------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
|                          | gCOD/m³| gCOD/m³| gCOD/m³| gCOD/m³| gCOD/m³| gCOD/m³| gCOD/m³| gCOD/m³| gCm⁻³ | gNm⁻³ | gNm⁻³ | gNm⁻³ | gNm⁻³ | gNm⁻³ | gPm⁻³ | gO₂m⁻³| gHm⁻³ |
| **Algae** | | | | | | | | | | | | | | | | | | |
| 1. Phototrophic growth on NH₄⁺         | 1     |       |       |      |      |      |       |       | α₁,₉  |       | α₁,₁₁ |       |       |       | α₁,₁₅ | 1     | α₁,₁₇ |
| 2. Phototrophic growth on NO₃⁻         | 1     |       |       |      |      |      |       |       | α₂,₉  |       |       |       | α₂,₁₃ |       | α₂,₁₅ |       | α₂,₁₇ |
| 3. Aerobic respiration                 | -1    |       |       |      |      |      |       |       | α₃,₉  |       | α₃,₁₁ |       |       |       | α₃,₁₅ | -1    | α₃,₁₇ |
| 4. Decay                               | -1    |       |       |      | α₄,₅ | α₄,₆ |       |       | α₄,₉  |       | α₄,₁₁ |       |       |       | α₄,₁₅ |       |       |
**Heterotrophic bacteria** | | | | | | | | | | | | | | | | | | |
| 5. Aerobic growth on NH₄⁺              |       |       |       | 1    |      |      | α₅,₇  |       | α₅,₉  |       | α₅,₁₁ |       |       |       | α₅,₁₅ | α₅,₁₆|       |
| 6. Aerobic growth on NO₃⁻              |       |       |       | 1    |      |      | α₆,₇  |       | α₆,₉  |       |       |       | α₆,₁₃ |       | α₆,₁₅ | α₆,₁₆|       |
| 7. Aerobic respiration                 |       |       |       | -1   |      |      |       |       | α₇,₉  |       | α₇,₁₁ |       |       |       | α₇,₁₅ | -1    |       |
| 8. Anoxic growth on NO₃⁻               |       |       |       | 1    |      |      | α₈,₇  |       | α₈,₉  |       | α₈,₁₁ |       | α₈,₁₃ | α₈,₁₄ | α₈,₁₅ |       |       |
| 9. Anoxic growth on NO₂⁻               |       |       |       | 1    |      |      | α₉,₇  |       | α₉,₉  |       | α₉,₁₁ | α₉,₁₂ |       | α₉,₁₄ | α₉,₁₅ |       |       |
| 10. Anoxic respiration NO₂⁻ and NO₃⁻    |       |       |       | -1   |      |      |       |       | α₁₀,₉ |       | α₁₀,₁₁| α₁₀,₁₂| α₁₀,₁₃| α₁₀,₁₄| α₁₀,₁₅|       |       |
| 11. Hydrolysis of slowly biodegradable COD |       |       |       |      | -1   |      | α₁₁,₇ | α₁₁,₈ | α₁₁,₉ |       | α₁₁,₁₁|       |       |       | α₁₁,₁₅|       |       |
| 12. Hydrolysis of urea                  |       |       |       |      |      |      |       |       | α₁₂,₉ | -1    | 1     |       |       |       |       |       | α₁₂,₁₇|
| 13. Decay                               |       |       |       | -1   | α₁₃,₅| α₁₃,₆|       |       | α₁₃,₉ |       | α₁₃,₁₁|       |       |       | α₁₃,₁₅|       |       |
**Ammonium Oxydising Bacteria** | | | | | | | | | | | | | | | | | | |
| 14. Aerobic growth on NH₄⁺             |       | 1     |       |      |      |      |       |       | α₁₄,₉ |       | α₁₄,₁₁| α₁₄,₁₂|       |       | α₁₄,₁₅| α₁₄,₁₆|       |
| 15. Aerobic respiration                |       | -1    |       |      |      |      |       |       | α₁₅,₉ |       | α₁₅,₁₁|       |       |       | α₁₅,₁₅| α₁₅,₁₆|       |
| 16. Decay                              |       | -1    |       |      | α₁₆,₅| α₁₆,₆|       |       | α₁₆,₉ |       | α₁₆,₁₁|       |       |       | α₁₆,₁₅|       |       |
**Nitrite Oxydising Bacteria** | | | | | | | | | | | | | | | | | | |
| 17. Aerobic growth on NO₂⁻             |       |       | 1     |      |      |      |       |       | α₁₇,₉ |       | α₁₇,₁₁| α₁₇,₁₂| α₁₇,₁₃|       | α₁₇,₁₅| α₁₇,₁₆|       |
| 18. Aerobic respiration                |       |       | -1    |      |      |      |       |       | α₁₈,₉ |       | α₁₈,₁₁|       |       |       | α₁₈,₁₅| α₁₈,₁₆|       |
| 19. Decay                              |       |       | -1    |      | α₁₉,₅| α₁₉,₆|       |       | α₁₉,₉ |       | α₁₉,₁₁|       |       |       | α₁₉,₁₅|       |       |
**Equilibrium phase** | | | | | | | | | | | | | | | | | | |
| 20. Dissolution of O₂                  |       |       |       |      |      |      |       |       |      |      |      |       |       |      |       | 1    |       |
| 21. Dissolution of CO₂                 |       |       |       |      |      |      |       |       | 1    |      |      |       |       |      |       |      |       |
| 22. Dissolution of NH₃                 |       |       |       |      |      |      |       |       |      |      | 1    |       |       |      |       |      |       |

---

### 2. Composition matrix

This matrix shows the **fixed chemical composition** of each variable.
*   **Rows:** The basic conserved elements (Carbon, Nitrogen, Phosphorus, COD, etc.).
*   **Columns:** The 17 state variables.
*   **Cells ($i_{k,j}$):** Indicate the amount of a chemical element $k$ that contains a unit of the component $j$ (e.g. "How many grams of Nitrogen are in 1 gram of Algae biomass").

| k↓ | j→ | X_ALG | X_AOB | X_NOB | X_H  | X_S  | X_I  | S_S   | S_I   | S_IC | S_ND | S_NH | S_NO2 | S_NO3 | S_N2 | S_PO4 | S_O2 | S_H2O |
|----|----|-------|-------|-------|------|------|------|-------|-------|------|------|------|-------|-------|------|-------|------|-------|
| 1  | COD (4 → 19) | 1     | 1     | 1     | 1    | 1    | 1    | 1     | 1     | 1    | 0    | 0    | -3.43 | -4.57 | -1.71| 0     | -1   | 0     |
| 2  | O (1 → 3)    | $iO_{BM,ALG}$ | $iO_{BM}$ | $iO_{BM}$ | $iO_{BM}$ | $iO_{XS}$ | $iO_{XI}$ | $iO_{SS}$ | $iO_{SI}$ | 2.67 | 0.57 | 0    | 2.28  | 3.43  | 0     | 2.07  | 1    | 7.94  |
| 3  | C           | $iC_{BM,ALG}$ | $iC_{BM}$ | $iC_{BM}$ | $iC_{BM}$ | $iC_{XS}$ | $iC_{XI}$ | $iC_{SS}$ | $iC_{SI}$ | 1    | 0.43 | 0    | 0     | 0     | 0     | 0     | 0    | 0     |
| 4  | N           | $iN_{BM,ALG}$ | $iN_{BM}$ | $iN_{BM}$ | $iN_{BM}$ | $iN_{XS}$ | $iN_{XI}$ | $iN_{SS}$ | $iN_{SI}$ | 0    | 1    | 1    | 1     | 1     | 0     | 0     | 0    | 0     |
| 5  | P           | $iP_{BM,ALG}$ | $iP_{BM}$ | $iP_{BM}$ | $iP_{BM}$ | $iP_{XS}$ | $iP_{XI}$ | $iP_{SS}$ | $iP_{SI}$ | 0    | 0    | 0    | 0     | 0     | 0     | 1     | 0    | 0     |
| 6  | H (1 → 3)   | $iH_{BM,ALG}$ | $iH_{BM}$ | $iH_{BM}$ | $iH_{BM}$ | 0     | 0     | 0     | 0     | 0.14 | 0.22 | 0.07 | 0.07  | 0.07  | 0     | 0.10  | 0    | 1     |

Just like other models, like ASM1 and ASM3, ALBA assumes that all bacteria (AOB, NOB, H) have the same empirical formula of biomass. Instead of defining a different composition for each type of bacteria, a standard composition is used for "bacterial biomass" ($C_{100}H_{183}O_{48}N_{11}P$). Therefore:
*   $iC_{BM}$ (Carbon in biomass) is the same for the three.
*   $iN_{BM}$ (Nitrogen in biomass) is the same for the three.
*   $iO_{BM}$, $iP_{BM}$, $iH_{BM}$ also repeat.

The column of the **Algae** shows that their parameters are **different** (for example, $iC_{BM,ALG}$). This is because algae have a different bio-chemical composition than bacteria (generally richer in carbon and with different nitrogen requirements).

It has two critical uses:
1.  **Calculation of Coefficients:** It is used to find the numerical values of the Petersen matrix. For example, to know how much $CO_2$ a bacteria consumes when growing, the model uses the composition matrix to balance the Carbon atoms.
2.  **Validation (Mass Balance):** It is used in the **unit tests** to ensure that no reaction "creates" matter from nothing. If you multiply a row of the Petersen matrix by the composition matrix, the result **must be zero** for all elements (C, N, P, etc.).

---

**Note:** The values like $iO_{BM,ALG}$, $iC_{BM,ALG}$, etc., are specific parameters of the ALBA model and must be replaced by their numerical values according to the table SI.3.2 of the original document.