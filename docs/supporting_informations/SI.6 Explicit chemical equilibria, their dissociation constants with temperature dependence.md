# SI.6 Explicit chemical equilibria, their dissociation constants with temperature dependence

## **Table SI.6.1: pH sub-model equation system**

| Description | Expression [mol m<sup>-3</sup>] | K<sub>A</sub> (293 K) [M] |
|-------------|--------------------------------|--------------------------|
| 1- Mass balance | $S_{\text{NH}_3} = \text{NH}_3 + \text{NH}_4^+$ | - |
| 2 - Dissociation $\text{NH}_4^+ \rightleftharpoons \text{NH}_3 + \text{H}^+$ | $\text{NH}_3 = \left(\frac{S_{\text{NH}_3}/14}{1 + \frac{K_{\text{a,NH}_4} \cdot 10^3}{\text{H}^+}}\right)$ | $K_{\text{a,NH}_4}: 5.62 \times 10^{-10}$ |
| 3- Mass balance | $\frac{S_{\text{NO}_2}}{14} = \text{NO}_2 + \text{HNO}_2$ | - |
| 4 - Dissociation $\text{HNO}_2 \rightleftharpoons \text{NO}_2^- + \text{H}^+$ | $\text{HNO}_2 = \left(\frac{S_{\text{NO}_2}/14}{1 + \frac{K_{\text{a,NO}_2} \cdot 10^3}{\text{H}^+}}\right)$ | $K_{\text{a,HNO}_2}: 4.47 \times 10^{-4}$ |
| 5- Mass balance | $\frac{S_{\text{NO}_3}}{14} = \text{NO}_3^- + \text{HNO}_3$ | - |
| 6 - Dissociation $\text{HNO}_3 \rightleftharpoons \text{NO}_3^- + \text{H}^+$ | $\text{HNO}_3 = \left(\frac{S_{\text{NO}_3}/14}{1 + \frac{K_{\text{a,NO}_3} \cdot 10^3}{\text{H}^+}}\right)$ | $K_{\text{a,HNO}_3}: 4.37 \times 10^{1}$ |
| 7- Mass balance | $\frac{S_{\text{IC}}}{12} = \text{CO}_2 + \text{HCO}_3^- + \text{CO}_3^{2-}$ | - |
| 8 - Dissociation $\text{H}_2\text{O} + \text{CO}_2 \rightleftharpoons \text{HCO}_3^- + \text{H}^+$ | $\text{CO}_2 = \frac{\frac{S_{\text{IC}}}{12}}{1 + \frac{K_{\text{a,CO}_2} \cdot 10^3}{\text{H}^+} + \frac{K_{\text{a,CO}_2} \cdot K_{\text{a,HCO}_3} \cdot 10^6}{\text{H}^2}}$ | $K_{\text{a,H}_2\text{CO}_3}: 4.27 \times 10^{-7}$ |
| 9 - Dissociation $\text{HCO}_3^- \rightleftharpoons \text{CO}_3^{2-} + \text{H}^+$ | $\text{HCO}_3^- = \frac{\frac{S_{\text{IC}}}{12}}{1 + \frac{\text{H}^+}{K_{\text{a,HCO}_3} \cdot 10^3} + \frac{K_{\text{a,HCO}_3} \cdot 10^3}{\text{H}^+}}$ | $K_{\text{a,HCO}_3}: 4.68 \times 10^{-11}$ |
| 10- Mass balance | $\frac{S_{\text{PO}_4}}{31} = \text{H}_3\text{PO}_4 + \text{H}_2\text{PO}_4^- + \text{HPO}_4^{2-} + \text{PO}_4^{3-}$ | - |
| 11 - Dissociation $\text{H}_3\text{PO}_4 \rightleftharpoons \text{H}_2\text{PO}_4^- + \text{H}^+$ | $\text{H}_3\text{PO}_4 = \frac{S_{\text{PO}_4}/31}{1 + \frac{K_{\text{a,H}_3\text{PO}_4} \cdot 10^3}{\text{H}^+} + \frac{K_{\text{a,H}_3\text{PO}_4} \cdot K_{\text{a,H}_2\text{PO}_4} \cdot 10^6}{\text{H}^2} + \frac{K_{\text{a,H}_3\text{PO}_4} \cdot K_{\text{a,H}_2\text{PO}_4} \cdot K_{\text{a,HPO}_4} \cdot 10^9}{\text{H}^3}}$ | $K_{\text{a,H}_3\text{PO}_4}: 7.24 \times 10^{-3}$ |
| 12 - Dissociation $\text{H}_2\text{PO}_4^- \rightleftharpoons \text{HPO}_4^{2-} + \text{H}^+$ | $\text{H}_2\text{PO}_4^- = \frac{S_{\text{PO}_4}/31}{1 + \frac{\text{H}^+}{K_{\text{a,H}_2\text{PO}_4} \cdot 10^3} + \frac{K_{\text{a,H}_2\text{PO}_4} \cdot 10^3}{\text{H}^+} + \frac{K_{\text{a,H}_2\text{PO}_4} \cdot K_{\text{a,HPO}_4} \cdot 10^6}{\text{H}^2}}$ | $K_{\text{a,H}_2\text{PO}_4}: 6.17 \times 10^{-8}$ |
| 13 - Dissociation $\text{HPO}_4^{2-} \rightleftharpoons \text{PO}_4^{3-} + \text{H}^+$ | $\text{HPO}_4^{2-} = \frac{S_{\text{PO}_4}/31}{1 + \frac{\text{H}^2}{K_{\text{a,H}_2\text{PO}_4} \cdot K_{\text{a,HPO}_4} \cdot 10^6} + \frac{\text{H}^+}{K_{\text{a,HPO}_4} \cdot 10^3} + \frac{K_{\text{a,HPO}_4} \cdot 10^3}{\text{H}^+}}$ | $K_{\text{a,HPO}_4}: 2.14 \times 10^{-13}$ |
| 14 - Dissociation $\text{H}_2\text{O} \rightleftharpoons \text{OH}^- + \text{H}^+$ | $\text{OH}^- = \frac{K_{\text{w}} \cdot 10^3}{\text{H}^+}$ | $K_{\text{w}}: 1.00 \times 10^{-14}$ |
| 15 - Charge balance | $\text{H}^+ + \text{NH}_4^+ + \Delta \text{CAT}_{\text{AN}} - \text{OH}^- - \text{NO}_2^- - \text{NO}_3^- - \text{HCO}_3^- - 2\text{CO}_3^{2-} - \text{H}_2\text{PO}_4^- - 2\text{HPO}_4^{2-} - 3\text{PO}_4^{3-} = 0$ | - |

---

As matter of illustration, the implemented equations for bicarbonate are shown below. Through the Henderson-Hasselbach formula, it is possible to compute the inorganic carbon ionic fractionation and then derive the amount of $\text{CO}_2$, $\text{HCO}_3^-$, $\text{CO}_3^{2-}$, according to the pH simulated. Conversion factors are required for every chemical equilibrium to transform the mass (or COD) concentrations into molar concentrations. Total inorganic carbon (S<sub>IC</sub> in equation SI6.2) is divided for carbon molecular weight to obtain the value in [molC m<sup>-3</sup>], while acidity constants ($K_{\text{a,CO}_2}$, $K_{\text{a,HCO}_3}$ in equation SI6.2 and SI6.3) are multiplied for $10^3$ [l m<sup>-3</sup>] and $10^6$ [l<sup>2</sup> m<sup>-6</sup>], since their value is typically reported in [M] in literature. The complete system of algebraic equations of pH sub-model can be found in Table SI.6.1.

$$
\frac{S_{\text{IC}}}{12} = \text{CO}_2 + \text{HCO}_3^- + \text{CO}_3^{2-} \tag{SI.6.1}
$$

$$
\text{CO}_2 = \frac{\frac{S_{\text{IC}}}{12}}{1 + \frac{K_{\text{a,CO}_2} \cdot 10^3}{\text{H}_{\text{ION}}^+} + \frac{K_{\text{a,CO}_2} \cdot K_{\text{a,HCO}_3} \cdot 10^6}{\text{H}_{\text{ION}}^2}} \tag{SI.6.2}
$$

$$
\text{HCO}_3^- = \frac{\frac{S_{\text{IC}}}{12}}{1 + \frac{\text{H}_{\text{ION}}^+}{K_{\text{a,HCO}_3} \cdot 10^3} + \frac{K_{\text{a,HCO}_3} \cdot 10^3}{\text{H}_{\text{ION}}^+}} \tag{SI.6.3}
$$

All the full set of equations are summarized in Table SI.6.1.

The temperature influence on the dissociation constants was taken into account by using the van't Hoff equation:

$$
\ln\left(\frac{K_{\text{a,T}}}{K_{\text{a,Tref}}}\right) = \frac{\Delta H^\circ}{R} \left(\frac{1}{T_{\text{ref}}} - \frac{1}{T + 273.15}\right) \tag{SI.6.4}
$$

In Equation SI.6.4, $T_{\text{ref}}$ is the standard temperature (298.15 K) for which the equilibrium coefficient value ($K_{\text{a,Tref}}$, [mol L<sup>-1</sup>]) is known, $T$ is the temperature at which we want to know the equilibrium coefficient value ($K_{\text{a,T}}$, [mol L<sup>-1</sup>]), $R$ is the gas law constant [J K<sup>-1</sup> mol<sup>-1</sup>] and $\Delta H^\circ$ is the heat of reaction at standard temperature and pressure [J]. 