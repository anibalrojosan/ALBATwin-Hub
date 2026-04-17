### SI 1.1 Estimation of the light extinction coefficient

The PAR was measured at different depths of the reactor. The probe was maintained at six points along depth (2, 5, 8, 12, 16 and 20 cm), registering the data for half an hour. The TSS concentration in the algal suspension was measured for each test. The Beer-Lambert equation was used to describe light extinction with the depth $z$ [m]:

$$\ln\left(\frac{I(z)}{I(0)}\right) = -\epsilon \cdot c \cdot z \tag{SI1.1}$$

where $I(z)$ is the PAR value [$\mu mol \cdot m^{-2} \cdot s^{-1}$] measured at depth $z$ [m]; $I(0)$ is the PAR value [$\mu mol \cdot m^{-2} \cdot s^{-1}$] measured at the pond surface ($z=0$); $c$ is the algal suspension TSS [$g \cdot m^{-3}$]; $\epsilon$ is the light extinction coefficient [$m^2 \cdot g^{-1}$]. The logarithm values of light data were evaluated for the different reactor depths (0 – 20 cm) during each test. The light extinction coefficient and its confidence interval were then estimated through linear regression (*fitlm* and *confCI* functions in MATLAB R2019b). Results are reported in table SI.1.1.

**Table SI.1.1: TSS measurements and estimated light extinction coefficients with corresponding 95% confidence intervals for each test**

| Test date | TSS | $\epsilon$ | Confidence interval on $\epsilon$ |
| :--- | :---: | :---: | :---: |
| | [$gTSS \cdot m^{-3}$] | [$m^2 \cdot gTSS^{-1}$] | [$m^2 \cdot gTSS^{-1}$] |
| 07/09/2018 | 183 | 0.113 | [0.094, 0.131] |
| 27/09/2018 | 300 | 0.095 | [0.072, 0.119] |
| 04/10/2018 | 292 | 0.083 | [0.062, 0.104] |
| 24/10/2018 | 212 | 0.112 | [0.082, 0.143] |
| **Avg. ± St.Dev** | **247 ± 58** | **0.101 ± 0.014** | - |

---

### Notes for the mathematical model:

* **Key Value**: Use the average value of $\epsilon = 0.101 \pm 0.014 \, m^2/gTSS$ as the default parameter for light attenuation in your simulation.
* **Relationship with TSS**: Note that attenuation depends directly on the concentration of total suspended solids ($c$); this means that as algae grow, light penetration in the reactor decreases (**self-shading**).

¿Necesitas que ajuste alguna otra variable o que integremos esto en un bloque de código específico?

