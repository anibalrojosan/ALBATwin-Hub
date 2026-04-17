### SI 1.2 Correlation between absorption measurements at 680 nm and algal biomass

Additional tests were performed to determine the correlation factor between the optical density measured at 680 nm (O.D.680 nm) (Helios Epsilon, Thermo Scientific) and the algal biomass ($X_{alg, meas}$, $gCOD \cdot m^{-3}$). Different batch with a volume of 500 mL, containing Tris-Acetate-Phosphate Medium (TAP) without organic carbon source and $NaHCO_3$ concentrations equal to 40 mM and 60 mM (Merck KGaA, Darmstadt, Germany) were inoculated with the biomass from the raceway mainly composed by *Chlorella sp.* The algal biomass ($gCOD \cdot m^{-3}$) was derived from dry weight measurements (TSS) using the conversion factor **1.57 gCOD · gALG⁻¹**.

![Figure SI.1.2: Correlation between optical density at 680 nm and algal biomass (gCOD m-3)](../figures/Figure%20SI%201.2%20Correlation%20between%20optical%20density%20at%20680%20nm%20and%20algal%20biomass%20gCOD%20m-3.png)

**Figure SI.1.2: Correlation between optical density at 680 nm and algal biomass ($gCOD \cdot m^{-3}$)**

The correlation equation obtained from the linear regression is:
$$y = 824.48x$$
Where:
*   $y$: Algal biomass ($X_{ALG}$) in $gCOD \cdot m^{-3}$.
*   $x$: Optical density at 680 nm ($OD_{680}$).

---

### Key data for implementation:

1.  **Conversion Factor Biomass-DQO:** The value of **1.57 gCOD / gTSS** is fundamental to convert laboratory measurements (solids) to the biomass units used in the ALBA model ($gCOD$).
2.  **Quick Estimation:** If in the future I need to integrate a turbidity or spectrophotometry sensor into the "Twin Hub", I can use the constant **824.48** to estimate the biomass in real time from the absorbance.
3.  **Species:** The model is calibrated specifically for *Chlorella sp.*
