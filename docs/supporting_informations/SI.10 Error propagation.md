# SI.10 Error propagation

After performing the sensitivity analysis and estimating the parameter standard error, as described in SI.9, the error propagation $\sigma_{y_i}$ of the model predictions for $y_i$ was computed as:

$$
\sigma_{y_i}(t) = \sqrt{\sum_{j=1}^m \left(\frac{\partial y_i}{\partial p_j}(t)\right)^2 \sigma_{p_j}^2} \tag{SI.10.1}
$$

Where $p_j$ are the model parameters, $\sigma_{p_j}$ their standard deviations, $y_i(p_1, ..., p_m)$ is the model solution for each predicted state $y_i$ at a given time $t$ and $\sigma_{y_i}$ is the prediction standard deviation of the model result.

Then, the 95% confidence intervals on model predictions (TSS, COD<sub>S</sub>, X<sub>ALG</sub>, S<sub>NH</sub>, S<sub>NO2</sub>, S<sub>NO3</sub>, S<sub>O2</sub> and pH) were estimated by the interval $[y_i - 1.96 \sigma_{y_i}, y_i + 1.96 \sigma_{y_i}]$ shown in Fig. 2, Section 4.2.2.