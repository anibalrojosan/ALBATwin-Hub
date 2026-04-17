# SI.7 Gas-liquid mass transfer

The different mass transfer coefficients were expressed as a function of the oxygen one (Sperandio, 1997):

$$
\frac{k_{L}a_{j}}{k_{L}a_{\text{O}_2}} = \left(\frac{D_{S_j}}{D_{\text{O}_2}}\right)^{0.5} \tag{SI.7.1}
$$

where $D_{S_j}$ [m<sup>2</sup> s<sup>-1</sup>] represents the diffusivity coefficient for the gas $j$. Combining equations SI7.1, SI7.2 and 11 (see Section 3.2.2), the following expression for the kinetics was obtained:

$$
Q_j = k_{L}a_{\text{O}_2} \left(\frac{D_{S_j}}{D_{\text{O}_2}}\right)^{0.5} \left(H_{S_j}p_{S_j} - S_j\right) \tag{SI.7.2}
$$

In Equation SI.7.2 the mass transfer coefficient ($k_{L}a_{\text{O}_2}$), the Henry’s constant ($H_{S_j}$) and the diffusivity coefficient ($D_{S_j}$) are temperature dependent. Temperature dependence is expressed by the Arrhenius law. The temperature correction coefficient varies in the range 1.016-1.135, the value chosen in this study is 1.024 (Ginot et Hervé, 1994). Henry’s constant temperature dependence acts in an opposite direction. Lower temperatures correspond to higher gas solubility. The empirical functions proposed by Sander (2015) were implemented, as shown below (Equation SI7.3, SI7.4 and SI7.5):

$$
H_{\text{O}_2}(T) = 42.15 \cdot e^{1700 \left(\frac{1}{273.15 + T} - \frac{1}{298.15}\right)} \quad \left[\text{gO}_2\text{m}^{-3}\text{atm}^{-1}\right] \tag{SI.7.3}
$$

$$
H_{\text{CO}_2}(T) = \left[1511.13 \cdot e^{2400 \left(\frac{1}{273.15 + T} - \frac{1}{298.15}\right)}\right] \cdot \frac{12}{44} \quad \left[\text{gC-CO}_2\text{m}^{-3}\text{atm}^{-1}\right] \tag{SI.7.4}
$$

$$
H_{\text{NH}_3}(T) = \left[4.63 \cdot 10^5 \cdot e^{2100 \left(\frac{1}{273.15 + T} - \frac{1}{298.15}\right)}\right] \cdot \frac{14}{17} \quad \left[\text{gN-NH}_3\text{m}^{-3}\text{atm}^{-1}\right] \tag{SI.7.5}
$$

In Equation SI.7.2, the difference $(H_{S_j}p_{S_j} - S_j)$ for $\text{CO}_2$ and $\text{NH}_3$ must be written in order to consider only the form really subjected to stripping/dissolution (i.e. $\text{CO}_2$ and free ammonia). The complete expressions for gas-liquid mass transfer are reported in Table SI.7.1.

---

**Table SI.7.1. Gas-liquid mass transfer rates implemented in the ALBA model.**

| Process | Description | Unit |
|---------|-------------|------|
| $\rho_{20}$ oxygen stripping/dissolution | $\theta^{T-20} \cdot k_{L}a \cdot \left(H_{\text{O}_2}(T) \cdot p_{\text{O}_2} - S_{\text{O}_2}\right)$ | gO<sub>2</sub> m<sup>-3</sup> d<sup>-1</sup> |
| $\rho_{21}$ carbon dioxide stripping/dissolution | $\theta^{T-20} \cdot k_{L}a \cdot \left(\frac{D_{\text{CO}_2}}{D_{\text{O}_2}}\right)^{0.5} \cdot \left(H_{\text{CO}_2}(T) \cdot p_{\text{CO}_2} - \frac{S_{\text{IC}}}{1 + \frac{K_{\text{a,CO}_2} \cdot 10^3}{H_{\text{ION}}}}\right)$ | gC-CO<sub>2</sub> m<sup>-3</sup> d<sup>-1</sup> |
| $\rho_{22}$ ammonia stripping | $\theta^{T-20} \cdot k_{L}a \cdot \left(\frac{D_{\text{NH}_3}}{D_{\text{O}_2}}\right)^{0.5} \cdot \left(H_{\text{NH}_3}(T) \cdot p_{\text{NH}_3} - \frac{S_{\text{NH}}}{1 + \frac{K_{\text{a,NH}_3} \cdot 10^3}{H_{\text{ION}}}}\right)$ | gN-NH<sub>3</sub> m<sup>-3</sup> d<sup>-1</sup> |