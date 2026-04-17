# SI.8 Sensitivity analysis and calibration strategy

The model most sensitive parameters were determined using the available AQUASIM toolboxes for sensitivity analysis. The absolute-relative sensitivity function of model output $y_i$ to parameter $p_j$ is defined as below:

$$
\frac{\partial y_i^{o,r}}{\partial p_j^{o,r}} = \frac{p_j}{y_i} \frac{\partial y_i}{\partial p_j} \tag{SI.8.1}
$$

The sensitivity analysis was carried out accounting for the environmental conditions defining each season pattern (in terms of light, temperature and evaporation rate, see Fig. 1), therefore the parameters reported in Tab. SI.8.1 are the resulting most sensitive ones in every season investigated. The sensitivity functions were estimated running simulations under established periodic regime (see description in section 2.2).

These parameters were then calibrated, using the procedure described in Section 2.2 (see Eq. 1). The experimental data of dissolved oxygen and pH collected by online probes in the periods 02-21/10/2018 and 01-10/01/2019 were considered in the criterion defined by equation (1) to calibrate the model. The model was then run with the new set of parameters and validated on 414 days of monitoring campaign, covering therefore all the seasons (see section 2.2 for details and Fig. 2 for simulation results).

---

## Table SI.8.1. Most sensible parameters identified from the sensitivity analysis performed under periodic regime, the nominal and calibrated values with their standard deviation, including the most affected model variables.

| Parameter | Nominal value | Reference | Calibrated value ± std | Most affected variables |
|-----------|---------------|-----------|------------------------|--------------------------|
| Algae maximum specific growth rate [$\mu_{\text{max,g,ALG}}$] | 1.5 d<sup>-1</sup> | [Solimeno et al. 2019] | 2.5 ± 0.055 d<sup>-1</sup> | X<sub>ALG</sub>, S<sub>O2</sub>, pH |
| AOB maximum specific growth rate [$\mu_{\text{max,g,AOB}}$] | 0.9 d<sup>-1</sup> | [Arashiro et al. 2017] | 0.72 ± 0.005 d<sup>-1</sup> | X<sub>AOB</sub>, S<sub>O2</sub>, S<sub>NH</sub>, S<sub>NO2</sub>, pH |
| NOB maximum specific growth rate [$\mu_{\text{max,g,NOB}}$] | 0.67 d<sup>-1</sup> | [Arashiro et al. 2017] | 0.65 ± 0.023 d<sup>-1</sup> | X<sub>NOB</sub>, S<sub>O2</sub>, S<sub>NO3</sub>, S<sub>NO2</sub>, pH |
| Light optimal value for growth [I<sub>OPT</sub>] | 275 μmol m<sup>-2</sup> s<sup>-1</sup> | [Martinez et al. 2018] | 300 ± 3.814 μmol m<sup>-2</sup> s<sup>-1</sup> | X<sub>ALG</sub>, S<sub>O2</sub>, pH |
| Light extinction coefficient [ε] | 0.067 ± 0.001 m<sup>2</sup> gCOD<sup>-1</sup> | [measured] | - | X<sub>ALG</sub>, S<sub>O2</sub>, pH |
| Initial slope of PI curve [α] | 0.027 | [Martinez et al. 2018] | 0.01 ± 0.0003 μmol<sup>-1</sup> m<sup>2</sup> s<sup>1</sup> d<sup>-1</sup> | X<sub>ALG</sub>, S<sub>O2</sub>, pH, X<sub>AOB</sub>, X<sub>NOB</sub> |
| Mass transfer coefficient [k<sub>L</sub>a] | 25 d<sup>-1</sup> | [Decostere, 2016] | 34 ± 0.115 d<sup>-1</sup> | S<sub>O2</sub>, pH |
| Coefficient for temperature correction for hydrolysis [θ<sub>HYD</sub>] | 1.07 | [Reichert, 2001] | 1.04 ± 0.005 | S<sub>S</sub>, X<sub>H</sub>, S<sub>O2</sub> |
| Coefficient for temperature correction for ammonification [θ<sub>AMM</sub>] | 1.07 | [Reichert, 2001] | 1.12 ± 0.002 | S<sub>NH</sub>, X<sub>AOB</sub>, X<sub>NOB</sub>, pH, S<sub>O2</sub> |

---

| Parameter | Nominal Min | Nominal Opt | Nominal Max | Reference | Calibrated Min | Calibrated Opt | Calibrated Max | Most affected variables |
|-----------|-------------|-------------|-------------|-----------|----------------|----------------|----------------|-------------------------|
| Cardinal temperature values for X<sub>ALG</sub> [CTMI] | 1.1 | 32.5 | 39.3 |[Bernard & Rémond, 2012]| -10±1.524 | 20±0.148 | 42±0.513 | Biomass concentration |
| Cardinal temperature values for X<sub>AOB</sub> [CTMI] | 5 | 25-35 | 35 | [Jubany, 2007] | -8±0.741 | 24.5±0.232 | 40±0.817 | nutrient removal rates, S<sub>O2</sub>, pH |
| Cardinal temperature values for X<sub>NOB</sub> [CTMI] | 5 | 25-30 | 37 | [Jubany 2007] | -8±9.734 | 20±0.940 | 38.5±6.090 | - |
| Cardinal temperature values for X<sub>H</sub> [CTMI] | 5 | 40 | 47 | [Rosso et al. 1995] | -3±0.335 | 25±0.634 | 42±1.919 | - |
| Cardinal pH values for X<sub>ALG</sub> [CPM] | 2.24 | 7.34 | 10 | [Ippoliti et al. 2016] | 2±0.562 | 8.4±0.066 | 12±0.039 | - |
| Cardinal pH values for X<sub>AOB</sub> [CPM] | 5.8 | 7.8-8 | 9 | [Jubany 2007 ] | 5.8±0.355 | 8.1±0.078 | 12.4±0.115 | Biomass concentration |
| Cardinal pH values for X<sub>NOB</sub> [CPM] | 6.5 | 7.6-8 | 8.6 | [Jubany 2007 ] | 5±0.568 | 7.9±0.320 | 12.1±0.463 | nutrient removal rates, S<sub>O2</sub>, pH |
| Cardinal pH values for X<sub>H</sub> [CPM] | 4 | 7 | 9 | [Rosso et al. 1995] | 2±0.344 | 7±0.066 | 11.5±0.022 | - |