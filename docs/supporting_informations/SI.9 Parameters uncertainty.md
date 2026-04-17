# SI.9 Parameters uncertainty

Once the model was calibrated and validated, a dynamic sensitivity analysis was run, accounting for all the period covered from the monitoring campaign (15/05/2018 - 01/08/2019) and therefore using the actual environmental conditions. The sensitivity functions were then computed in these real conditions for all the parameters reported in Table SI.8.1:

$$
\frac{\partial y_i^{a,r}}{\partial p_j} = \frac{p_j}{y_i} \frac{\partial y_i}{\partial p_j} \tag{SI.9.1}
$$

The parameter standard deviation was then derived from the Fisher Information Matrix $F$. The Fisher analysis is based on the local sensitivity functions $\frac{\partial y_i^{a,r}}{\partial p_j}$, and turned out to be efficient for biological dynamic systems (Ejiofor et al., 1994; Vatcheva et al., 2006).

The matrix $F$ was computed from the sensitivity matrix $\Delta_{Y_p}$ (Eq. 24) and covariance matrix of measured standard deviation $C$:

$$
\Delta = \begin{bmatrix}
\frac{\partial y_1}{\partial p_1} & \cdots & \frac{\partial y_1}{\partial p_m} \\
\vdots & \ddots & \vdots \\
\frac{\partial y_n}{\partial p_1} & \cdots & \frac{\partial y_n}{\partial p_m}
\end{bmatrix} \tag{SI.9.2}
$$

$$
F = \sum_{k=1}^K \Delta_k^T C^{-1} \Delta_k \tag{SI.9.3}
$$

The standard deviation $\delta_j$ associated to parameters $p_j$ is then computed as:

$$
\delta_j^2 = \left(F^{-1}\right)_{jj} \tag{SI.9.4}
$$