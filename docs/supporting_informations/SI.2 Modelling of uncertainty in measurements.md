### SI 2. Modelling of uncertainty in measurements

The way standard deviation is estimated from the variation coefficient is presented in Table SI.2.1. For measurements lower than a threshold $\varpi$, the standard deviation is assumed to be constant.

**Table SI.2.1: measurement uncertainty modelling: standard deviation as a function of mean value $\omega$.**

| Measurement | Unit | Threshold $\varpi$ | Standard deviation ($\omega < \varpi$) | Standard deviation ($\varpi < \omega$) |
| :--- | :---: | :---: | :---: | :---: |
| DO | $mgO_2 \cdot L^{-1}$ | - | - | 5% $\omega$ |
| pH | - | - | - | 2% $\omega$ |
| sCOD | $mgCOD \cdot L^{-1}$ | 5 | 1 | 20% $\omega$ |
| $COD_{ALG}$ | $mgCOD \cdot L^{-1}$ | 5 | 1 | 20% $\omega$ |
| $P\text{-}PO_4^{3-}$ | $mgP \cdot L^{-1}$ | 5 | 1 | 20% $\omega$ |
| $N\text{-}NH_4$ | $mgN \cdot L^{-1}$ | 5 | 1 | 20% $\omega$ |
| $N\text{-}NO_3^-$ | $mgN \cdot L^{-1}$ | 5 | 1 | 20% $\omega$ |
| $N\text{-}NO_2$ | $mgN \cdot L^{-1}$ | 5 | 1 | 20% $\omega$ |

---

### Key implementation notes for the Digital Twin:

1.  **Sensor Error Modeling:** When I implement a Kalman Filter or any data assimilation method, these percentages (5% for DO, 2% for pH, 20% for nutrients) will define my **Measurement Noise Covariance Matrix ($R$)**.
2.  **Threshold Logic:** For nutrients and COD, if the concentration is very low (below 5 $mg/L$), the model should assume a fixed absolute error of 1 $mg/L$ instead of a percentage. This prevents the error from becoming zero when the substrate is depleted.
3.  **High Uncertainty:** Chemical measurements (sCOD, Nitrogen, Phosphorus) have a much higher uncertainty (20%) compared to physical probes like pH or DO.