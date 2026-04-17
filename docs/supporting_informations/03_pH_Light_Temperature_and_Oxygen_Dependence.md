# 3. pH, Light, Temperature and Oxygen Dependence

### Light

Light is a crucial factor for algae growth, driving a large fraction of the energy and carbon fluxes in the system. Describing its effect on photosynthesis in a turbid system is challenging since it is affected by many different factors and it is species dependent (Martínez et al. 2018). As stated before, light penetration in the HRABP was estimated through the Lambert-Beer law (see Equation SI.1.1) and the light extinction coefficient ($\varepsilon$) was experimentally determined, as reported in supplementary information (SI.1).

The light dependence of algal growth was described by a Haldane-type function (Eq. (5)), choosing the parametrization proposed by Bernard and Rémond (2012):

$$
f_I = \frac{\mu_{MAX}}{1 + \frac{\mu_{MAX}}{\alpha} \left(\frac{I}{I_{OPT}} - 1\right)^2}
\tag{5}
$$

### Temperature

Temperature deeply affects biological process rates, and this influence must definitely be considered for outdoor systems. In this study, temperature fluctuates within large ranges along the campaign period, through daily oscillations and seasonal changes.

The model chosen for simulating the temperature dependence of growth and respiration rates, both for algae and bacteria, is the CTMI (Cardinal Temperature Model with Inflection) proposed by Rosso et al. (1993), shown in Eq. (6). This function has been shown to efficiently describe biomass growth, especially at high temperatures. It requires three parameters (the cardinal temperatures: T$_{MAX}$, T$_{OPT}$, T$_{MIN}$), which define the optimal working range for each functional group.

An Arrhenius function, requiring only one parameter ($\theta$ in Eq. (7)), was implemented for modelling the decay rate dependence on temperature for both algae and bacteria. With this function, the decay rate increases with temperature. Nominal and calibrated cardinal temperature values are shown later on in Table SI.8.1.

$$
f_T =
  \begin{cases}
   0 & T < T_{min} \\
   \frac{(T-T_{max}) \cdot (T-T_{min})^2}{(T_{opt}-T_{min}) \cdot (T_{opt}-T_{min}) \cdot [(T-T_{opt})-(T_{opt}-T_{max}) \cdot ((T_{opt}-T_{min})-(T+T_{min}-2T_{opt}))]} & T_{min} \leq T \leq T_{max} \\
   0 & T > T_{max}
  \end{cases}
\tag{6}
$$

$$
f_T = \frac{\mu_{Decay}(T)}{\mu_{Decay}(20^{\circ}C)} = \theta^{(T-20)}
\tag{7}
$$
---

### pH

The pH strongly influences system dynamics, since it directly affects the speciation of soluble compounds ($S_{IC}$, $S_{NH}$, $S_{NO2}$, $S_{NO3}$, $S_{PO4}$) and their availability.

The pH of the raceway was not controlled, so that the system exhibited large daily pH fluctuations (up to 10.5 during day and down to 7 during night). The pH dependence was modelled using the CPM (Cardinal pH Model, without inflection, Eq. (8)) function proposed by Rosso et al. (1995).

The CPM requires three parameters (the cardinal pH: pH$_{MAX}$, pH$_{OPT}$, pH$_{MIN}$), defining the growing range for each biomass.

$$
f_{pH} =
  \begin{cases}
   0 & \text{pH} < \text{pH}_{min} \\
   \frac{(\text{pH} - \text{pH}_{min})(\text{pH} - \text{pH}_{max})}{(\text{pH} - \text{pH}_{min})(\text{pH} - \text{pH}_{max}) - (\text{pH} - \text{pH}_{opt})^2} & \text{pH}_{min} \leq \text{pH} \leq \text{pH}_{max} \\
   0 & \text{pH} > \text{pH}_{max}
  \end{cases}
\tag{8}
$$

Nominal and calibrated cardinal pH values are reported in Table SI.8.1.

---

### Oxygen

High dissolved oxygen concentrations can negatively affect the photosynthetic activity of phototrophic microorganisms (Peng et al., 2013). The reduction of photosynthetic activity at high DO concentrations can be described with an inhibition Hill-type model (Eq. (9)) in the growth rate (Di Veroli et al., 2015):

$$
f_{DO,g} = \frac{k_{DO}^n}{S_{O2}^n + k_{DO}^n}
\tag{9}
$$

where $k_{DO}$ is the inhibition parameter of the model and $n$ is the dimensionless Hill coefficient [-]. Oxygen is the substrate of algae respiration and its limiting effect is classically represented with a Monod-type function (see Table 3).

The effect of high DO concentration on algal decay was represented with a Hill-type model (Eq. (10)), as reported in Table 3. It represents the increase in the decay rate above a certain DO concentration ($k_{DO}$):

$$
f_{DO,d} = \frac{S_{O2}^n}{S_{O2}^n + k_{DO}^n}
\tag{10}
$$

---

