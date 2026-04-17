# Table SI.3.3: Stoichiometric coefficient expressions

**Stoichiometric coefficients**


| Symbol | Affected variable | Expression | Unit |
| ------ | :---------------: | :----------: | :-----: |
|$\rho_1$ - Growth of $X_{ALG}$ on $NH_4^+$ | | | |
| $\alpha_{1,1}$  | $X_{ALG}$         | $1$                                                                                                      | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{1,9}$  | $S_{IC}$          | $-iC_{BM}^{ALG}$                                                                                         | $gC / gCOD_{BM}$     |
| $\alpha_{1,11}$ | $S_{NH}$          | $-iN_{BM}^{ALG}$                                                                                         | $gN / gCOD_{BM}$     |
| $\alpha_{1,15}$ | $S_{PO4}$         | $-iP_{BM}^{ALG}$                                                                                         | $gP / gCOD_{BM}$     |
| $\alpha_{1,16}$ | $S_{O2}$          | $-iO_{BM}^{ALG} + (32/12)iC_{BM}^{ALG} - (24/14)iN_{BM}^{ALG} + (40/31)iP_{BM}^{ALG} + (8)iH_{BM}^{ALG}$ | $gO_2 / gCOD_{BM}$    |
| $\alpha_{1,17}$ | $S_{H2O}$         | $-0.0404$                                                                                                | $gH / gCOD_{BM}$     |
|$\rho_2$ - Growth of $X_{ALG}$ on $NO_3^-$ | | | |
| $\alpha_{2,1}$  | $X_{ALG}$         | $1$                                                                                                      | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{2,9}$  | $S_{IC}$          | $-iC_{BM}^{ALG}$                                                                                         | $gC / gCOD_{BM}$     |
| $\alpha_{2,13}$ | $S_{NO3}$         | $-iN_{BM}^{ALG}$                                                                                         | $gN / gCOD_{BM}$     |
| $\alpha_{2,15}$ | $S_{PO4}$         | $-iP_{BM}^{ALG}$                                                                                         | $gP / gCOD_{BM}$     |
| $\alpha_{2,16}$ | $S_{O2}$          | $-iO_{BM}^{ALG} + (32/12)iC_{BM}^{ALG} + (40/14)iN_{BM}^{ALG} + (40/31)iP_{BM}^{ALG} + (8)iH_{BM}^{ALG}$ | $gO_2 / gCOD_{BM}$     |
| $\alpha_{2,17}$ | $S_{H2O}$         | $-0.0464$                                                                                                | $gH / gCOD_{BM}$     |
|$\rho_3$ - Aerobic respiration of $X_{ALG}$ | | | |
| $\alpha_{3,1}$  | $X_{ALG}$         | $-1$                                                                                                     | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{3,9}$  | $S_{IC}$          | $iC_{BM}^{ALG}$                                                                                          | $gC / gCOD_{BM}$     |
| $\alpha_{3,11}$ | $S_{NH}$          | $iN_{BM}^{ALG}$                                                                                          | $gN / gCOD_{BM}$     |
| $\alpha_{3,15}$ | $S_{PO4}$         | $iP_{BM}^{ALG}$                                                                                          | $gP / gCOD_{BM}$     |
| $\alpha_{3,16}$ | $S_{O2}$          | $-iO_{BM}^{ALG} - (32/12)iC_{BM}^{ALG} + (24/14)iN_{BM}^{ALG} - (40/31)iP_{BM}^{ALG} - (8)iH_{BM}^{ALG}$ | $gO_2 / gCOD_{BM}$    |
| $\alpha_{3,17}$ | $S_{H2O}$         | $0.0404$                                                                                                 | $gH / gCOD_{BM}$     |
|$\rho_4$ - Decay of $X_{ALG}$ | | | |
| $\alpha_{4,1}$  | $X_{ALG}$         | $-1$                                                            | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{4,5}$  | $X_S$             | $(1-f_{XI,ALG})$                                                | $gCOD_{XS} / gCOD_{XS}$ |
| $\alpha_{4,6}$  | $X_I$             | $f_{XI,ALG}$                                                    | $gCOD_{XI} / gCOD_{XI}$ |
| $\alpha_{4,9}$  | $S_{IC}$          | $iC_{BM}^{ALG}-(1-f_{XI,ALG}) * iC_{XS} - f_{XI,ALG} * iC_{XI}$ | $gC / gCOD_{BM}$     |
| $\alpha_{4,11}$ | $S_{NH}$          | $iN_{BM}^{ALG}-(1-f_{XI,ALG}) * iN_{XS} - f_{XI,ALG} * iN_{XI}$ | $gN / gCOD_{BM}$     |
| $\alpha_{4,15}$ | $S_{PO4}$         | $iP_{BM}^{ALG}-(1-f_{XI,ALG}) * iP_{XS} - f_{XI,ALG} * iP_{XI}$ | $gP / gCOD_{BM}$     |
|$\rho_5$ - Aerobic growth of $X_H$ on $NH_4^+$ | | | |
| $\alpha_{5,4}$  | $X_H$             | $1$                     | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{5,7}$  | $S_S$             | $-1/Y_H$                | $gCOD_{SS} / gCOD_{SS}$ |
| $\alpha_{5,9}$  | $S_{IC}$          | $iC_{SS}/Y_H - iC_{BM}$ | $gC / gCOD_{BM}$     |
| $\alpha_{5,11}$ | $S_{NH}$          | $iN_{SS}/Y_H - iN_{BM}$ | $gN / gCOD_{BM}$     |
| $\alpha_{5,15}$ | $S_{PO4}$         | $iP_{SS}/Y_H - iP_{BM}$ | $gP / gCOD_{BM}$     |
| $\alpha_{5,16}$ | $S_{O2}$          | $-(1/Y_H - 1)$          | $gO_2 / gCOD_{BM}$    |
|$\rho_6$ - Aerobic growth of $X_H$ on $NO_3^-$ | | | |
| $\alpha_{6,4}$  | $X_H$             | $1$                                            | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{6,7}$  | $S_S$             | $-1/Y_H$                                       | $gCOD_{SS} / gCOD_{SS}$ |
| $\alpha_{6,9}$  | $S_{IC}$          | $iC_{SS}/Y_H - iC_{BM}$                        | $gC / gCOD_{BM}$     |
| $\alpha_{6,13}$ | $S_{NO3}$         | $iN_{SS}/Y_H - iN_{BM}$                        | $gN / gCOD_{BM}$     |
| $\alpha_{6,15}$ | $S_{PO4}$         | $iP_{SS}/Y_H - iP_{BM}$                        | $gP / gCOD_{BM}$     |
| $\alpha_{6,16}$ | $S_{O2}$          | $-(1/Y_H - 1) - 64/14*(iN_{SS}/Y_H - iN_{BM})$ | $gO_2 / gCOD_{BM}$    |
|$\rho_7$ - Aerobic respiration of $X_H$ | | | |
| $\alpha_{7,4}$  | $X_H$             | $-1$       | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{7,9}$  | $S_{IC}$          | $iC_{BM}$  | $gC / gCOD_{BM}$     |
| $\alpha_{7,11}$ | $S_{NH}$          | $iN_{BM}$  | $gN / gCOD_{BM}$     |
| $\alpha_{7,15}$ | $S_{PO4}$         | $iP_{BM}$  | $gP / gCOD_{BM}$     |
| $\alpha_{7,16}$ | $S_{O2}$          | $-1$       | $gO_2 / gCOD_{BM}$    |
|$\rho_8$ - Anoxic growth of $X_H$ on $NO_3^-$ | | | |
| $\alpha_{8,4}$  | $X_H$             | $1$                          | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{8,7}$  | $S_S$             | $-1/Y_{HNO3}$                | $gCOD_{SS} / gCOD_{SS}$ |
| $\alpha_{8,9}$  | $S_{IC}$          | $iC_{SS}/Y_{HNO3} - iC_{BM}$ | $gC / gCOD_{BM}$     |
| $\alpha_{8,11}$ | $S_{NH}$          | $iN_{SS}/Y_{HNO3} - iN_{BM}$ | $gN / gCOD_{BM}$     |
| $\alpha_{8,13}$ | $S_{NO3}$         | $-28/80*(1/Y_{HNO3} - 1)$    | $gN / gCOD_{BM}$     |
| $\alpha_{8,14}$ | $S_{N2}$          | $28/80*(1/Y_{HNO3} - 1)$     | $gN / gCOD_{BM}$     |
| $\alpha_{8,15}$ | $S_{PO4}$         | $iP_{SS}/Y_{HNO3} - iP_{BM}$ | $gP / gCOD_{BM}$     |
|$\rho_9$ - Anoxic growth of $X_H$ on $NO_2^-$ | | | |
| $\alpha_{9,4}$  | $X_H$             | $1$                          | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{9,7}$  | $S_S$             | $-1/Y_{HNO2}$                | $gCOD_{SS} / gCOD_{SS}$ |
| $\alpha_{9,9}$  | $S_{IC}$          | $iC_{SS}/Y_{HNO2} - iC_{BM}$ | $gC / gCOD_{BM}$     |
| $\alpha_{9,11}$ | $S_{NH}$          | $iN_{SS}/Y_{HNO2} - iN_{BM}$ | $gN / gCOD_{BM}$     |
| $\alpha_{9,12}$ | $S_{NO2}$         | $-28/48*(1/Y_{HNO2} - 1)$    | $gN / gCOD_{BM}$     |
| $\alpha_{9,14}$ | $S_{N2}$          | $28/48*(1/Y_{HNO2} - 1)$     | $gN / gCOD_{BM}$     |
| $\alpha_{9,15}$ | $S_{PO4}$         | $iP_{SS}/Y_{HNO2} - iP_{BM}$ | $gP / gCOD_{BM}$     |
|$\rho_{10}$ - Anoxic respiration of $X_H$ on $NO_2^-$ and $NO_3^-$ | | | |
| $\alpha_{10,4}$  | $X_H$             | $-1$       | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{10,9}$  | $S_{IC}$          | $iC_{BM}$  | $gC / gCOD_{BM}$     |
| $\alpha_{10,11}$ | $S_{NH}$          | $iN_{BM}$  | $gN / gCOD_{BM}$     |
| $\alpha_{10,12}$ | $S_{NO2}$         | $-14/64$   | $gN / gCOD_{BM}$     |
| $\alpha_{10,13}$ | $S_{NO3}$         | $-14/64$   | $gN / gCOD_{BM}$     |
| $\alpha_{10,14}$ | $S_{N2}$          | $28/64$    | $gN / gCOD_{BM}$     |
| $\alpha_{10,15}$ | $S_{PO4}$         | $iP_{BM}$  | $gP / gCOD_{BM}$     |
|$\rho_{11}$ - Hydrolysis of slowly biodegradable COD | | | |
| $\alpha_{11,5}$  | $X_S$             | $-1$                                               | $gCOD_{XS} / gCOD_{XS}$ |
| $\alpha_{11,7}$  | $S_S$             | $1-f_{SI}$                                        | $gCOD_{SS} / gCOD_{XS}$ |
| $\alpha_{11,8}$  | $S_I$             | $f_{SI}$                                          | $gCOD_{SI} / gCOD_{XS}$ |
| $\alpha_{11,9}$  | $S_{IC}$          | $iC_{XS}-(1-f_{SI}) * iC_{SS} - f_{SI} * iC_{SI}$ | $gC / gCOD_{XS}$     |
| $\alpha_{11,11}$ | $S_{NH}$          | $iN_{XS}-(1-f_{SI}) * iN_{SS} - f_{SI} * iN_{SI}$ | $gN / gCOD_{XS}$     |
| $\alpha_{11,15}$ | $S_{PO4}$         | $iP_{XS}-(1-f_{SI}) * iP_{SS} - f_{SI} * iP_{SI}$ | $gP / gCOD_{XS}$     |
|$\rho_{12}$ - Hydrolysis of urea | | | |
| $\alpha_{12,9}$  | $S_{IC}$          | $iC_{ND}$  | $gC / gN_{urea}$        |
| $\alpha_{12,10}$ | $S_{ND}$          | $-1$       | $gN_{urea} / gN_{urea}$       |
| $\alpha_{12,11}$ | $S_{NH}$          | $1$        | $gN_{ammonia} / gN_{urea}$       |
| $\alpha_{12,17}$ | $S_{H2O}$         | $iH_{ND}$  | $gH / gN_{urea}$        |
|$\rho_{13}$ - Decay of $X_H$ | | | |
| $\alpha_{13,4}$  | $X_H$             | $-1$                                              | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{13,5}$  | $X_S$             | $1-f_{XI}$                                        | $gCOD_{XS} / gCOD_{XS}$ |
| $\alpha_{13,6}$  | $X_I$             | $f_{XI}$                                          | $gCOD_{XI} / gCOD_{XI}$ |
| $\alpha_{13,9}$  | $S_{IC}$          | $iC_{BM}-(1-f_{XI}) * iC_{XS} - f_{XI} * iC_{XI}$ | $gC / gCOD_{BM}$     |
| $\alpha_{13,11}$ | $S_{NH}$          | $iN_{BM}-(1-f_{XI}) * iN_{XS} - f_{XI} * iN_{XI}$ | $gN / gCOD_{BM}$     |
| $\alpha_{13,15}$ | $S_{PO4}$         | $iP_{BM}-(1-f_{XI}) * iP_{XS} - f_{XI} * iP_{XI}$ | $gP / gCOD_{BM}$     |
|$\rho_{14}$ - Aerobic growth of $X_{AOB}$ on $NH_4^+$ | | | |
| $\alpha_{14,2}$  | $X_{AOB}$         | $1$                    | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{14,9}$  | $S_{IC}$          | $-iC_{BM}$             | $gC / gCOD_{BM}$     |
| $\alpha_{14,11}$ | $S_{NH}$          | $-1/Y_{AOB} - iN_{BM}$ | $gN / gCOD_{BM}$     |
| $\alpha_{14,12}$ | $S_{NO2}$         | $1/Y_{AOB}$            | $gN / gCOD_{BM}$     |
| $\alpha_{14,15}$ | $S_{PO4}$         | $-iP_{BM}$             | $gP / gCOD_{BM}$     |
| $\alpha_{14,16}$ | $S_{O2}$          | $1-48/14*1/Y_{AOB}$    | $gO_2 / gCOD_{BM}$    |
$\rho_{15}$ - Aerobic respiration of $X_{AOB}$ | | | |
| $\alpha_{15,2}$  | $X_{AOB}$         | $-1$       | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{15,9}$  | $S_{IC}$          | $iC_{BM}$  | $gC / gCOD_{BM}$     |
| $\alpha_{15,11}$ | $S_{NH}$          | $iN_{BM}$  | $gN / gCOD_{BM}$     |
| $\alpha_{15,15}$ | $S_{PO4}$         | $iP_{BM}$  | $gP / gCOD_{BM}$     |
| $\alpha_{15,16}$ | $S_{O2}$          | $-1$       | $gO_2 / gCOD_{BM}$    |
$\rho_{16}$ - Decay of $X_{AOB}$ | | | |
| $\alpha_{16,2}$  | $X_{AOB}$         | $-1$                                              | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{16,5}$  | $X_S$             | $1-f_{XI}$                                        | $gCOD_{XS} / gCOD_{XS}$ |
| $\alpha_{16,6}$  | $X_I$             | $f_{XI}$                                          | $gCOD_{XI} / gCOD_{XI}$ |
| $\alpha_{16,9}$  | $S_{IC}$          | $iC_{BM}-(1-f_{XI}) * iC_{XS} - f_{XI} * iC_{XI}$ | $gC / gCOD_{BM}$     |
| $\alpha_{16,11}$ | $S_{NH}$          | $iN_{BM}-(1-f_{XI}) * iN_{XS} - f_{XI} * iN_{XI}$ | $gN / gCOD_{BM}$     |
| $\alpha_{16,15}$ | $S_{PO4}$         | $iP_{BM}-(1-f_{XI}) * iP_{XS} - f_{XI} * iP_{XI}$ | $gP / gCOD_{BM}$     |
$\rho_{17}$ - Aerobic growth of $X_{NOB}$ on $NO_2^-$ | | | |
| $\alpha_{17,3}$  | $X_{NOB}$         | $1$                 | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{17,9}$  | $S_{IC}$          | $-iC_{BM}$          | $gC / gCOD_{BM}$     |
| $\alpha_{17,11}$ | $S_{NH}$          | $-iN_{BM}$          | $gN / gCOD_{BM}$     |
| $\alpha_{17,12}$ | $S_{NO2}$         | $-1/Y_{NOB}$        | $gN / gCOD_{BM}$     |
| $\alpha_{17,13}$ | $S_{NO3}$         | $1/Y_{NOB}$         | $gN / gCOD_{BM}$     |
| $\alpha_{17,15}$ | $S_{PO4}$         | $-iP_{BM}$          | $gP / gCOD_{BM}$     |
| $\alpha_{17,16}$ | $S_{O2}$          | $1-16/14*1/Y_{NOB}$ | $gO_2 / gCOD_{BM}$    |
$\rho_{18}$ - Aerobic respiration of $X_{NOB}$ | | | |
| $\alpha_{18,3}$  | $X_{NOB}$         | $-1$       | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{18,9}$  | $S_{IC}$          | $iC_{BM}$  | $gC / gCOD_{BM}$     |
| $\alpha_{18,11}$ | $S_{NH}$          | $iN_{BM}$  | $gN / gCOD_{BM}$     |
| $\alpha_{18,15}$ | $S_{PO4}$         | $iP_{BM}$  | $gP / gCOD_{BM}$     |
| $\alpha_{18,16}$ | $S_{O2}$          | $-1$       | $gO_2 / gCOD_{BM}$    |
$\rho_{19}$ - Decay of $X_{NOB}$ | | | |
| $\alpha_{19,3}$  | $X_{NOB}$         | $-1$                                              | $gCOD_{BM} / gCOD_{BM}$ |
| $\alpha_{19,5}$  | $X_S$             | $1-f_{XI}$                                        | $gCOD_{XS} / gCOD_{XS}$ |
| $\alpha_{19,6}$  | $X_I$             | $f_{XI}$                                          | $gCOD_{XI} / gCOD_{XI}$ |
| $\alpha_{19,9}$  | $S_{IC}$          | $iC_{BM}-(1-f_{XI}) * iC_{XS} - f_{XI} * iC_{XI}$ | $gC / gCOD_{BM}$     |
| $\alpha_{19,11}$ | $S_{NH}$          | $iN_{BM}-(1-f_{XI}) * iN_{XS} - f_{XI} * iN_{XI}$ | $gN / gCOD_{BM}$     |
| $\alpha_{19,15}$ | $S_{PO4}$         | $iP_{BM}-(1-f_{XI}) * iP_{XS} - f_{XI} * iP_{XI}$ | $gP / gCOD_{BM}$     |
$\rho_{20}$ - Dissolution of $O_2$ | | | |
| $\alpha_{20,15}$ | $S_{O2}$          | $1$        | $[-]$ |
$\rho_{21}$ - Dissolution of $CO_2$ | | | |
| $\alpha_{21,9}$ | $S_{IC}$          | $1$        | $[-]$ |
$\rho_{22}$ - Dissolution of $NH_3$ | | | |
| $\alpha_{22,11}$ | $S_{NH}$          | $1$        | $[-]$ |
