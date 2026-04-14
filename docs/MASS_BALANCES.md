# Mass balances: explicit $\mathbf{S}\mathbf{I}^\top$ cells (ALBA stoichiometry)

This file documents the **elemental and COD-style mass-balance check** used in unit tests: for each biological process row $i$ of the Petersen matrix $\mathbf{S}$ and each row $k$ of the composition matrix $\mathbf{I}$ (COD, O, C, N, P, H),

$$B_{i,k} = \sum_{j=0}^{16} S_{i,j}\, I_{k,j}\,.$$

Process index $i$ is **0-based** (`S[i, :]`, same as code). Casagli process numbering $\rho_r$ uses $r = i + 1$. A value near zero means the row is consistent with that composition row; large $|B_{i,k}|$ indicates missing or implicit species in $\mathbf{S}$ (see ADR 007 for O/H policy).

**Tolerance** for OK/FAIL labels matches `tests/unit/stoichiometry_mass_balance_shared.py`: `MASS_BALANCE_ATOL = 0.1` (|residual| $\le$ atol $\Rightarrow$ OK).

**Regeneration:** `uv run python scripts/generate_mass_balances_md.py`

---

## Table of contents

- [Summary table (all 114 cells)](#summary-table-all-114-cells)
- [Per-cell derivations](#per-cell-derivations)

- **rho1** (phototrophic growth X_ALG (NH4+))
  - [COD](#rho1-COD)
  - [O](#rho1-O)
  - [C](#rho1-C)
  - [N](#rho1-N)
  - [P](#rho1-P)
  - [H](#rho1-H)
- **rho2** (phototrophic growth X_ALG (NO3-))
  - [COD](#rho2-COD)
  - [O](#rho2-O)
  - [C](#rho2-C)
  - [N](#rho2-N)
  - [P](#rho2-P)
  - [H](#rho2-H)
- **rho3** (aerobic respiration X_ALG)
  - [COD](#rho3-COD)
  - [O](#rho3-O)
  - [C](#rho3-C)
  - [N](#rho3-N)
  - [P](#rho3-P)
  - [H](#rho3-H)
- **rho4** (decay X_ALG)
  - [COD](#rho4-COD)
  - [O](#rho4-O)
  - [C](#rho4-C)
  - [N](#rho4-N)
  - [P](#rho4-P)
  - [H](#rho4-H)
- **rho5** (aerobic growth X_H (NH4+))
  - [COD](#rho5-COD)
  - [O](#rho5-O)
  - [C](#rho5-C)
  - [N](#rho5-N)
  - [P](#rho5-P)
  - [H](#rho5-H)
- **rho6** (aerobic growth X_H (NO3-))
  - [COD](#rho6-COD)
  - [O](#rho6-O)
  - [C](#rho6-C)
  - [N](#rho6-N)
  - [P](#rho6-P)
  - [H](#rho6-H)
- **rho7** (aerobic respiration X_H)
  - [COD](#rho7-COD)
  - [O](#rho7-O)
  - [C](#rho7-C)
  - [N](#rho7-N)
  - [P](#rho7-P)
  - [H](#rho7-H)
- **rho8** (anoxic growth X_H (NO3-))
  - [COD](#rho8-COD)
  - [O](#rho8-O)
  - [C](#rho8-C)
  - [N](#rho8-N)
  - [P](#rho8-P)
  - [H](#rho8-H)
- **rho9** (anoxic growth X_H (NO2-))
  - [COD](#rho9-COD)
  - [O](#rho9-O)
  - [C](#rho9-C)
  - [N](#rho9-N)
  - [P](#rho9-P)
  - [H](#rho9-H)
- **rho10** (anoxic respiration X_H)
  - [COD](#rho10-COD)
  - [O](#rho10-O)
  - [C](#rho10-C)
  - [N](#rho10-N)
  - [P](#rho10-P)
  - [H](#rho10-H)
- **rho11** (hydrolysis X_S)
  - [COD](#rho11-COD)
  - [O](#rho11-O)
  - [C](#rho11-C)
  - [N](#rho11-N)
  - [P](#rho11-P)
  - [H](#rho11-H)
- **rho12** (urea hydrolysis)
  - [COD](#rho12-COD)
  - [O](#rho12-O)
  - [C](#rho12-C)
  - [N](#rho12-N)
  - [P](#rho12-P)
  - [H](#rho12-H)
- **rho13** (decay X_H)
  - [COD](#rho13-COD)
  - [O](#rho13-O)
  - [C](#rho13-C)
  - [N](#rho13-N)
  - [P](#rho13-P)
  - [H](#rho13-H)
- **rho14** (aerobic growth X_AOB)
  - [COD](#rho14-COD)
  - [O](#rho14-O)
  - [C](#rho14-C)
  - [N](#rho14-N)
  - [P](#rho14-P)
  - [H](#rho14-H)
- **rho15** (aerobic respiration X_AOB)
  - [COD](#rho15-COD)
  - [O](#rho15-O)
  - [C](#rho15-C)
  - [N](#rho15-N)
  - [P](#rho15-P)
  - [H](#rho15-H)
- **rho16** (decay X_AOB)
  - [COD](#rho16-COD)
  - [O](#rho16-O)
  - [C](#rho16-C)
  - [N](#rho16-N)
  - [P](#rho16-P)
  - [H](#rho16-H)
- **rho17** (aerobic growth X_NOB)
  - [COD](#rho17-COD)
  - [O](#rho17-O)
  - [C](#rho17-C)
  - [N](#rho17-N)
  - [P](#rho17-P)
  - [H](#rho17-H)
- **rho18** (aerobic respiration X_NOB)
  - [COD](#rho18-COD)
  - [O](#rho18-O)
  - [C](#rho18-C)
  - [N](#rho18-N)
  - [P](#rho18-P)
  - [H](#rho18-H)
- **rho19** (decay X_NOB)
  - [COD](#rho19-COD)
  - [O](#rho19-O)
  - [C](#rho19-C)
  - [N](#rho19-N)
  - [P](#rho19-P)
  - [H](#rho19-H)

---

## Summary table (all 114 cells)

| $\rho$ | Element | Residual $B_{i,k}$ | $\|B_{i,k}\|$ | Status |
|:---:|:---:|:---:|:---:|:---:|
| 1 | COD | -0.0013225806 | 0.0013225806 | **OK** |
| 1 | O | -0.0001034194 | 0.0001034194 | **OK** |
| 1 | C | 0 | 0 | **OK** |
| 1 | N | 0 | 0 | **OK** |
| 1 | P | 0 | 0 | **OK** |
| 1 | H | -0.00044 | 0.00044 | **OK** |
| 2 | COD | -0.0013825806 | 0.0013825806 | **OK** |
| 2 | O | 0.0001965806 | 0.0001965806 | **OK** |
| 2 | C | 0 | 0 | **OK** |
| 2 | N | 0 | 0 | **OK** |
| 2 | P | 0 | 0 | **OK** |
| 2 | H | -0.00014 | 0.00014 | **OK** |
| 3 | COD | 0.0013225806 | 0.0013225806 | **OK** |
| 3 | O | 0.0001034194 | 0.0001034194 | **OK** |
| 3 | C | 0 | 0 | **OK** |
| 3 | N | 0 | 0 | **OK** |
| 3 | P | 0 | 0 | **OK** |
| 3 | H | 0.00044 | 0.00044 | **OK** |
| 4 | COD | -5.551115e-17 | 5.551115e-17 | **OK** |
| 4 | O | -0.03072638 | 0.03072638 | **OK** |
| 4 | C | -7.806256e-18 | 7.806256e-18 | **OK** |
| 4 | N | -2.602085e-18 | 2.602085e-18 | **OK** |
| 4 | P | 0 | 0 | **OK** |
| 4 | H | -0.04832564 | 0.04832564 | **OK** |
| 5 | COD | 0 | 0 | **OK** |
| 5 | O | -0.2810977778 | 0.2810977778 | **FAIL** |
| 5 | C | 0 | 0 | **OK** |
| 5 | N | 0 | 0 | **OK** |
| 5 | P | 0 | 0 | **OK** |
| 5 | H | 0.028951746 | 0.028951746 | **OK** |
| 6 | COD | -8.598639e-05 | 8.598639e-05 | **OK** |
| 6 | O | -0.2123946485 | 0.2123946485 | **FAIL** |
| 6 | C | 0 | 0 | **OK** |
| 6 | N | 0 | 0 | **OK** |
| 6 | P | 0 | 0 | **OK** |
| 6 | H | 0.0379803175 | 0.0379803175 | **OK** |
| 7 | COD | 0 | 0 | **OK** |
| 7 | O | -0.18968 | 0.18968 | **FAIL** |
| 7 | C | 0 | 0 | **OK** |
| 7 | N | 0 | 0 | **OK** |
| 7 | P | 0 | 0 | **OK** |
| 7 | H | -0.02292 | 0.02292 | **OK** |
| 8 | COD | 0.001 | 0.001 | **OK** |
| 8 | O | -0.604 | 0.604 | **FAIL** |
| 8 | C | 0 | 0 | **OK** |
| 8 | N | 0 | 0 | **OK** |
| 8 | P | 0 | 0 | **OK** |
| 8 | H | 0.00602 | 0.00602 | **OK** |
| 9 | COD | 0.0077777778 | 0.0077777778 | **OK** |
| 9 | O | -1.5689533333 | 1.5689533333 | **FAIL** |
| 9 | C | 0 | 0 | **OK** |
| 9 | N | 0 | 0 | **OK** |
| 9 | P | -1.301043e-18 | 1.301043e-18 | **OK** |
| 9 | H | -0.0596911111 | 0.0596911111 | **OK** |
| 10 | COD | 0.001875 | 0.001875 | **OK** |
| 10 | O | -0.4387425 | 0.4387425 | **FAIL** |
| 10 | C | 0 | 0 | **OK** |
| 10 | N | 0 | 0 | **OK** |
| 10 | P | 0 | 0 | **OK** |
| 10 | H | -0.053545 | 0.053545 | **OK** |
| 11 | COD | 2.775558e-17 | 2.775558e-17 | **OK** |
| 11 | O | -0.011814 | 0.011814 | **OK** |
| 11 | C | 3.469447e-18 | 3.469447e-18 | **OK** |
| 11 | N | 0 | 0 | **OK** |
| 11 | P | -2.975051e-19 | 2.975051e-19 | **OK** |
| 11 | H | 0.00319 | 0.00319 | **OK** |
| 12 | COD | 0 | 0 | **OK** |
| 12 | O | 1.14711 | 1.14711 | **FAIL** |
| 12 | C | -0.001 | 0.001 | **OK** |
| 12 | N | 0 | 0 | **OK** |
| 12 | P | 0 | 0 | **OK** |
| 12 | H | 0.152 | 0.152 | **FAIL** |
| 13 | COD | 2.775558e-17 | 2.775558e-17 | **OK** |
| 13 | O | 0.094061 | 0.094061 | **OK** |
| 13 | C | 0 | 0 | **OK** |
| 13 | N | 0 | 0 | **OK** |
| 13 | P | 0 | 0 | **OK** |
| 13 | H | -0.031522 | 0.031522 | **OK** |
| 14 | COD | -0.0071428571 | 0.0071428571 | **OK** |
| 14 | O | -5.5531771429 | 5.5531771429 | **FAIL** |
| 14 | C | 0 | 0 | **OK** |
| 14 | N | 0 | 0 | **OK** |
| 14 | P | 0 | 0 | **OK** |
| 14 | H | -0.72708 | 0.72708 | **FAIL** |
| 15 | COD | 0 | 0 | **OK** |
| 15 | O | -0.18968 | 0.18968 | **FAIL** |
| 15 | C | 0 | 0 | **OK** |
| 15 | N | 0 | 0 | **OK** |
| 15 | P | 0 | 0 | **OK** |
| 15 | H | -0.02292 | 0.02292 | **OK** |
| 16 | COD | 2.775558e-17 | 2.775558e-17 | **OK** |
| 16 | O | 0.094061 | 0.094061 | **OK** |
| 16 | C | 0 | 0 | **OK** |
| 16 | N | 0 | 0 | **OK** |
| 16 | P | 0 | 0 | **OK** |
| 16 | H | -0.031522 | 0.031522 | **OK** |
| 17 | COD | 0.0571428571 | 0.0571428571 | **OK** |
| 17 | O | 0.3325371429 | 0.3325371429 | **FAIL** |
| 17 | C | 0 | 0 | **OK** |
| 17 | N | 0 | 0 | **OK** |
| 17 | P | 0 | 0 | **OK** |
| 17 | H | 0.02292 | 0.02292 | **OK** |
| 18 | COD | 0 | 0 | **OK** |
| 18 | O | -0.18968 | 0.18968 | **FAIL** |
| 18 | C | 0 | 0 | **OK** |
| 18 | N | 0 | 0 | **OK** |
| 18 | P | 0 | 0 | **OK** |
| 18 | H | -0.02292 | 0.02292 | **OK** |
| 19 | COD | 0 | 0 | **OK** |
| 19 | O | 0.094061 | 0.094061 | **OK** |
| 19 | C | 0 | 0 | **OK** |
| 19 | N | 0 | 0 | **OK** |
| 19 | P | 0 | 0 | **OK** |
| 19 | H | -0.031522 | 0.031522 | **OK** |

---

## Per-cell derivations

For each cell: contributing columns are those with $S_{i,j} \neq 0$. Products use numeric values from `get_petersen_matrix()` and `get_composition_matrix()`.

### rho1

*Full label:* rho1: phototrophic growth X_ALG (NH4+)

<a id="rho1-COD"></a>
#### rho1 — COD (cell 1/114)

Balance for process row $i=0$ ($\rho_{1}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 1 + 0 + 0 + 0 + -1.0013225806 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | 1 | 1 | 1 |
| 8 | `S_IC` | -0.327 | 0 | 0 |
| 10 | `S_NH` | -0.042 | 0 | 0 |
| 14 | `S_PO4` | -0.008 | 0 | 0 |
| 15 | `S_O2` | 1.0013225806 | -1 | -1.0013225806 |
| 16 | `S_H2O` | -0.0404 | 0 | 0 |

**Sum:** $B_{0,0} \approx -0.0013225806$ (check: `numpy.sum` = -0.0013225806)

**Status:** **OK** — $|B_{0,0}| = 0.0013225806$ ≤ 0.1

<a id="rho1-O"></a>
#### rho1 — O (cell 2/114)

Balance for process row $i=0$ ($\rho_{1}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.209 + -0.87309 + 0 + -0.01656 + 1.0013225806 + -0.320776$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | 1 | 0.209 | 0.209 |
| 8 | `S_IC` | -0.327 | 2.67 | -0.87309 |
| 10 | `S_NH` | -0.042 | 0 | 0 |
| 14 | `S_PO4` | -0.008 | 2.07 | -0.01656 |
| 15 | `S_O2` | 1.0013225806 | 1 | 1.0013225806 |
| 16 | `S_H2O` | -0.0404 | 7.94 | -0.320776 |

**Sum:** $B_{0,1} \approx -0.0001034194$ (check: `numpy.sum` = -0.0001034194)

**Status:** **OK** — $|B_{0,1}| = 0.0001034194$ ≤ 0.1

<a id="rho1-C"></a>
#### rho1 — C (cell 3/114)

Balance for process row $i=0$ ($\rho_{1}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.327 + -0.327 + 0 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | 1 | 0.327 | 0.327 |
| 8 | `S_IC` | -0.327 | 1 | -0.327 |
| 10 | `S_NH` | -0.042 | 0 | 0 |
| 14 | `S_PO4` | -0.008 | 0 | 0 |
| 15 | `S_O2` | 1.0013225806 | 0 | 0 |
| 16 | `S_H2O` | -0.0404 | 0 | 0 |

**Sum:** $B_{0,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{0,2}| = 0$ ≤ 0.1

<a id="rho1-N"></a>
#### rho1 — N (cell 4/114)

Balance for process row $i=0$ ($\rho_{1}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.042 + 0 + -0.042 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | 1 | 0.042 | 0.042 |
| 8 | `S_IC` | -0.327 | 0 | 0 |
| 10 | `S_NH` | -0.042 | 1 | -0.042 |
| 14 | `S_PO4` | -0.008 | 0 | 0 |
| 15 | `S_O2` | 1.0013225806 | 0 | 0 |
| 16 | `S_H2O` | -0.0404 | 0 | 0 |

**Sum:** $B_{0,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{0,3}| = 0$ ≤ 0.1

<a id="rho1-P"></a>
#### rho1 — P (cell 5/114)

Balance for process row $i=0$ ($\rho_{1}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.008 + 0 + 0 + -0.008 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | 1 | 0.008 | 0.008 |
| 8 | `S_IC` | -0.327 | 0 | 0 |
| 10 | `S_NH` | -0.042 | 0 | 0 |
| 14 | `S_PO4` | -0.008 | 1 | -0.008 |
| 15 | `S_O2` | 1.0013225806 | 0 | 0 |
| 16 | `S_H2O` | -0.0404 | 0 | 0 |

**Sum:** $B_{0,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{0,4}| = 0$ ≤ 0.1

<a id="rho1-H"></a>
#### rho1 — H (cell 6/114)

Balance for process row $i=0$ ($\rho_{1}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.05 + 0 + -0.00924 + -0.0008 + 0 + -0.0404$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | 1 | 0.05 | 0.05 |
| 8 | `S_IC` | -0.327 | 0 | 0 |
| 10 | `S_NH` | -0.042 | 0.22 | -0.00924 |
| 14 | `S_PO4` | -0.008 | 0.1 | -0.0008 |
| 15 | `S_O2` | 1.0013225806 | 0 | 0 |
| 16 | `S_H2O` | -0.0404 | 1 | -0.0404 |

**Sum:** $B_{0,5} \approx -0.00044$ (check: `numpy.sum` = -0.00044)

**Status:** **OK** — $|B_{0,5}| = 0.00044$ ≤ 0.1

### rho2

*Full label:* rho2: phototrophic growth X_ALG (NO3-)

<a id="rho2-COD"></a>
#### rho2 — COD (cell 7/114)

Balance for process row $i=1$ ($\rho_{2}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 1 + 0 + 0.19194 + 0 + -1.1933225806 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | 1 | 1 | 1 |
| 8 | `S_IC` | -0.327 | 0 | 0 |
| 12 | `S_NO3` | -0.042 | -4.57 | 0.19194 |
| 14 | `S_PO4` | -0.008 | 0 | 0 |
| 15 | `S_O2` | 1.1933225806 | -1 | -1.1933225806 |
| 16 | `S_H2O` | -0.0464 | 0 | 0 |

**Sum:** $B_{1,0} \approx -0.0013825806$ (check: `numpy.sum` = -0.0013825806)

**Status:** **OK** — $|B_{1,0}| = 0.0013825806$ ≤ 0.1

<a id="rho2-O"></a>
#### rho2 — O (cell 8/114)

Balance for process row $i=1$ ($\rho_{2}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.209 + -0.87309 + -0.14406 + -0.01656 + 1.1933225806 + -0.368416$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | 1 | 0.209 | 0.209 |
| 8 | `S_IC` | -0.327 | 2.67 | -0.87309 |
| 12 | `S_NO3` | -0.042 | 3.43 | -0.14406 |
| 14 | `S_PO4` | -0.008 | 2.07 | -0.01656 |
| 15 | `S_O2` | 1.1933225806 | 1 | 1.1933225806 |
| 16 | `S_H2O` | -0.0464 | 7.94 | -0.368416 |

**Sum:** $B_{1,1} \approx 0.0001965806$ (check: `numpy.sum` = 0.0001965806)

**Status:** **OK** — $|B_{1,1}| = 0.0001965806$ ≤ 0.1

<a id="rho2-C"></a>
#### rho2 — C (cell 9/114)

Balance for process row $i=1$ ($\rho_{2}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.327 + -0.327 + 0 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | 1 | 0.327 | 0.327 |
| 8 | `S_IC` | -0.327 | 1 | -0.327 |
| 12 | `S_NO3` | -0.042 | 0 | 0 |
| 14 | `S_PO4` | -0.008 | 0 | 0 |
| 15 | `S_O2` | 1.1933225806 | 0 | 0 |
| 16 | `S_H2O` | -0.0464 | 0 | 0 |

**Sum:** $B_{1,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{1,2}| = 0$ ≤ 0.1

<a id="rho2-N"></a>
#### rho2 — N (cell 10/114)

Balance for process row $i=1$ ($\rho_{2}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.042 + 0 + -0.042 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | 1 | 0.042 | 0.042 |
| 8 | `S_IC` | -0.327 | 0 | 0 |
| 12 | `S_NO3` | -0.042 | 1 | -0.042 |
| 14 | `S_PO4` | -0.008 | 0 | 0 |
| 15 | `S_O2` | 1.1933225806 | 0 | 0 |
| 16 | `S_H2O` | -0.0464 | 0 | 0 |

**Sum:** $B_{1,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{1,3}| = 0$ ≤ 0.1

<a id="rho2-P"></a>
#### rho2 — P (cell 11/114)

Balance for process row $i=1$ ($\rho_{2}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.008 + 0 + 0 + -0.008 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | 1 | 0.008 | 0.008 |
| 8 | `S_IC` | -0.327 | 0 | 0 |
| 12 | `S_NO3` | -0.042 | 0 | 0 |
| 14 | `S_PO4` | -0.008 | 1 | -0.008 |
| 15 | `S_O2` | 1.1933225806 | 0 | 0 |
| 16 | `S_H2O` | -0.0464 | 0 | 0 |

**Sum:** $B_{1,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{1,4}| = 0$ ≤ 0.1

<a id="rho2-H"></a>
#### rho2 — H (cell 12/114)

Balance for process row $i=1$ ($\rho_{2}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.05 + 0 + -0.00294 + -0.0008 + 0 + -0.0464$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | 1 | 0.05 | 0.05 |
| 8 | `S_IC` | -0.327 | 0 | 0 |
| 12 | `S_NO3` | -0.042 | 0.07 | -0.00294 |
| 14 | `S_PO4` | -0.008 | 0.1 | -0.0008 |
| 15 | `S_O2` | 1.1933225806 | 0 | 0 |
| 16 | `S_H2O` | -0.0464 | 1 | -0.0464 |

**Sum:** $B_{1,5} \approx -0.00014$ (check: `numpy.sum` = -0.00014)

**Status:** **OK** — $|B_{1,5}| = 0.00014$ ≤ 0.1

### rho3

*Full label:* rho3: aerobic respiration X_ALG

<a id="rho3-COD"></a>
#### rho3 — COD (cell 13/114)

Balance for process row $i=2$ ($\rho_{3}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -1 + 0 + 0 + 0 + 1.0013225806 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | -1 | 1 | -1 |
| 8 | `S_IC` | 0.327 | 0 | 0 |
| 10 | `S_NH` | 0.042 | 0 | 0 |
| 14 | `S_PO4` | 0.008 | 0 | 0 |
| 15 | `S_O2` | -1.0013225806 | -1 | 1.0013225806 |
| 16 | `S_H2O` | 0.0404 | 0 | 0 |

**Sum:** $B_{2,0} \approx 0.0013225806$ (check: `numpy.sum` = 0.0013225806)

**Status:** **OK** — $|B_{2,0}| = 0.0013225806$ ≤ 0.1

<a id="rho3-O"></a>
#### rho3 — O (cell 14/114)

Balance for process row $i=2$ ($\rho_{3}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.209 + 0.87309 + 0 + 0.01656 + -1.0013225806 + 0.320776$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | -1 | 0.209 | -0.209 |
| 8 | `S_IC` | 0.327 | 2.67 | 0.87309 |
| 10 | `S_NH` | 0.042 | 0 | 0 |
| 14 | `S_PO4` | 0.008 | 2.07 | 0.01656 |
| 15 | `S_O2` | -1.0013225806 | 1 | -1.0013225806 |
| 16 | `S_H2O` | 0.0404 | 7.94 | 0.320776 |

**Sum:** $B_{2,1} \approx 0.0001034194$ (check: `numpy.sum` = 0.0001034194)

**Status:** **OK** — $|B_{2,1}| = 0.0001034194$ ≤ 0.1

<a id="rho3-C"></a>
#### rho3 — C (cell 15/114)

Balance for process row $i=2$ ($\rho_{3}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.327 + 0.327 + 0 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | -1 | 0.327 | -0.327 |
| 8 | `S_IC` | 0.327 | 1 | 0.327 |
| 10 | `S_NH` | 0.042 | 0 | 0 |
| 14 | `S_PO4` | 0.008 | 0 | 0 |
| 15 | `S_O2` | -1.0013225806 | 0 | 0 |
| 16 | `S_H2O` | 0.0404 | 0 | 0 |

**Sum:** $B_{2,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{2,2}| = 0$ ≤ 0.1

<a id="rho3-N"></a>
#### rho3 — N (cell 16/114)

Balance for process row $i=2$ ($\rho_{3}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.042 + 0 + 0.042 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | -1 | 0.042 | -0.042 |
| 8 | `S_IC` | 0.327 | 0 | 0 |
| 10 | `S_NH` | 0.042 | 1 | 0.042 |
| 14 | `S_PO4` | 0.008 | 0 | 0 |
| 15 | `S_O2` | -1.0013225806 | 0 | 0 |
| 16 | `S_H2O` | 0.0404 | 0 | 0 |

**Sum:** $B_{2,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{2,3}| = 0$ ≤ 0.1

<a id="rho3-P"></a>
#### rho3 — P (cell 17/114)

Balance for process row $i=2$ ($\rho_{3}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.008 + 0 + 0 + 0.008 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | -1 | 0.008 | -0.008 |
| 8 | `S_IC` | 0.327 | 0 | 0 |
| 10 | `S_NH` | 0.042 | 0 | 0 |
| 14 | `S_PO4` | 0.008 | 1 | 0.008 |
| 15 | `S_O2` | -1.0013225806 | 0 | 0 |
| 16 | `S_H2O` | 0.0404 | 0 | 0 |

**Sum:** $B_{2,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{2,4}| = 0$ ≤ 0.1

<a id="rho3-H"></a>
#### rho3 — H (cell 18/114)

Balance for process row $i=2$ ($\rho_{3}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.05 + 0 + 0.00924 + 0.0008 + 0 + 0.0404$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | -1 | 0.05 | -0.05 |
| 8 | `S_IC` | 0.327 | 0 | 0 |
| 10 | `S_NH` | 0.042 | 0.22 | 0.00924 |
| 14 | `S_PO4` | 0.008 | 0.1 | 0.0008 |
| 15 | `S_O2` | -1.0013225806 | 0 | 0 |
| 16 | `S_H2O` | 0.0404 | 1 | 0.0404 |

**Sum:** $B_{2,5} \approx 0.00044$ (check: `numpy.sum` = 0.00044)

**Status:** **OK** — $|B_{2,5}| = 0.00044$ ≤ 0.1

### rho4

*Full label:* rho4: decay X_ALG

<a id="rho4-COD"></a>
#### rho4 — COD (cell 19/114)

Balance for process row $i=3$ ($\rho_{4}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -1 + 0.938 + 0.062 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | -1 | 1 | -1 |
| 4 | `X_S` | 0.938 | 1 | 0.938 |
| 5 | `X_I` | 0.062 | 1 | 0.062 |
| 8 | `S_IC` | 0.006396 | 0 | 0 |
| 10 | `S_NH` | 0.006388 | 0 | 0 |
| 14 | `S_PO4` | 0.00269 | 0 | 0 |

**Sum:** $B_{3,0} \approx -5.551115e-17$ (check: `numpy.sum` = -5.551115e-17)

**Status:** **OK** — $|B_{3,0}| = 5.551115e-17$ ≤ 0.1

<a id="rho4-O"></a>
#### rho4 — O (cell 20/114)

Balance for process row $i=3$ ($\rho_{4}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.209 + 0.146328 + 0.0093 + 0.01707732 + 0 + 0.0055683$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | -1 | 0.209 | -0.209 |
| 4 | `X_S` | 0.938 | 0.156 | 0.146328 |
| 5 | `X_I` | 0.062 | 0.15 | 0.0093 |
| 8 | `S_IC` | 0.006396 | 2.67 | 0.01707732 |
| 10 | `S_NH` | 0.006388 | 0 | 0 |
| 14 | `S_PO4` | 0.00269 | 2.07 | 0.0055683 |

**Sum:** $B_{3,1} \approx -0.03072638$ (check: `numpy.sum` = -0.03072638)

**Status:** **OK** — $|B_{3,1}| = 0.03072638$ ≤ 0.1

<a id="rho4-C"></a>
#### rho4 — C (cell 21/114)

Balance for process row $i=3$ ($\rho_{4}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.327 + 0.298284 + 0.02232 + 0.006396 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | -1 | 0.327 | -0.327 |
| 4 | `X_S` | 0.938 | 0.318 | 0.298284 |
| 5 | `X_I` | 0.062 | 0.36 | 0.02232 |
| 8 | `S_IC` | 0.006396 | 1 | 0.006396 |
| 10 | `S_NH` | 0.006388 | 0 | 0 |
| 14 | `S_PO4` | 0.00269 | 0 | 0 |

**Sum:** $B_{3,2} \approx -7.806256e-18$ (check: `numpy.sum` = -7.806256e-18)

**Status:** **OK** — $|B_{3,2}| = 7.806256e-18$ ≤ 0.1

<a id="rho4-N"></a>
#### rho4 — N (cell 22/114)

Balance for process row $i=3$ ($\rho_{4}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.042 + 0.031892 + 0.00372 + 0 + 0.006388 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | -1 | 0.042 | -0.042 |
| 4 | `X_S` | 0.938 | 0.034 | 0.031892 |
| 5 | `X_I` | 0.062 | 0.06 | 0.00372 |
| 8 | `S_IC` | 0.006396 | 0 | 0 |
| 10 | `S_NH` | 0.006388 | 1 | 0.006388 |
| 14 | `S_PO4` | 0.00269 | 0 | 0 |

**Sum:** $B_{3,3} \approx -2.602085e-18$ (check: `numpy.sum` = -2.602085e-18)

**Status:** **OK** — $|B_{3,3}| = 2.602085e-18$ ≤ 0.1

<a id="rho4-P"></a>
#### rho4 — P (cell 23/114)

Balance for process row $i=3$ ($\rho_{4}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.008 + 0.00469 + 0.00062 + 0 + 0 + 0.00269$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | -1 | 0.008 | -0.008 |
| 4 | `X_S` | 0.938 | 0.005 | 0.00469 |
| 5 | `X_I` | 0.062 | 0.01 | 0.00062 |
| 8 | `S_IC` | 0.006396 | 0 | 0 |
| 10 | `S_NH` | 0.006388 | 0 | 0 |
| 14 | `S_PO4` | 0.00269 | 1 | 0.00269 |

**Sum:** $B_{3,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{3,4}| = 0$ ≤ 0.1

<a id="rho4-H"></a>
#### rho4 — H (cell 24/114)

Balance for process row $i=3$ ($\rho_{4}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.05 + 0 + 0 + 0 + 0.00140536 + 0.000269$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 0 | `X_ALG` | -1 | 0.05 | -0.05 |
| 4 | `X_S` | 0.938 | 0 | 0 |
| 5 | `X_I` | 0.062 | 0 | 0 |
| 8 | `S_IC` | 0.006396 | 0 | 0 |
| 10 | `S_NH` | 0.006388 | 0.22 | 0.00140536 |
| 14 | `S_PO4` | 0.00269 | 0.1 | 0.000269 |

**Sum:** $B_{3,5} \approx -0.04832564$ (check: `numpy.sum` = -0.04832564)

**Status:** **OK** — $|B_{3,5}| = 0.04832564$ ≤ 0.1

### rho5

*Full label:* rho5: aerobic growth X_H (NH4+)

<a id="rho5-COD"></a>
#### rho5 — COD (cell 25/114)

Balance for process row $i=4$ ($\rho_{5}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 1 + -1.5873015873 + 0 + 0 + 0 + 0.5873015873$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 1 | 1 |
| 6 | `S_S` | -1.5873015873 | 1 | -1.5873015873 |
| 8 | `S_IC` | 0.1447619048 | 0 | 0 |
| 10 | `S_NH` | -0.0601904762 | 0 | 0 |
| 14 | `S_PO4` | -0.0080634921 | 0 | 0 |
| 15 | `S_O2` | -0.5873015873 | -1 | 0.5873015873 |

**Sum:** $B_{4,0} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{4,0}| = 0$ ≤ 0.1

<a id="rho5-O"></a>
#### rho5 — O (cell 26/114)

Balance for process row $i=4$ ($\rho_{5}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.184 + -0.2476190476 + 0.3865142857 + 0 + -0.0166914286 + -0.5873015873$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.184 | 0.184 |
| 6 | `S_S` | -1.5873015873 | 0.156 | -0.2476190476 |
| 8 | `S_IC` | 0.1447619048 | 2.67 | 0.3865142857 |
| 10 | `S_NH` | -0.0601904762 | 0 | 0 |
| 14 | `S_PO4` | -0.0080634921 | 2.07 | -0.0166914286 |
| 15 | `S_O2` | -0.5873015873 | 1 | -0.5873015873 |

**Sum:** $B_{4,1} \approx -0.2810977778$ (check: `numpy.sum` = -0.2810977778)

**Status:** **FAIL** — $|B_{4,1}| = 0.2810977778$ > 0.1

<a id="rho5-C"></a>
#### rho5 — C (cell 27/114)

Balance for process row $i=4$ ($\rho_{5}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.36 + -0.5047619048 + 0.1447619048 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.36 | 0.36 |
| 6 | `S_S` | -1.5873015873 | 0.318 | -0.5047619048 |
| 8 | `S_IC` | 0.1447619048 | 1 | 0.1447619048 |
| 10 | `S_NH` | -0.0601904762 | 0 | 0 |
| 14 | `S_PO4` | -0.0080634921 | 0 | 0 |
| 15 | `S_O2` | -0.5873015873 | 0 | 0 |

**Sum:** $B_{4,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{4,2}| = 0$ ≤ 0.1

<a id="rho5-N"></a>
#### rho5 — N (cell 28/114)

Balance for process row $i=4$ ($\rho_{5}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.084 + -0.0238095238 + 0 + -0.0601904762 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.084 | 0.084 |
| 6 | `S_S` | -1.5873015873 | 0.015 | -0.0238095238 |
| 8 | `S_IC` | 0.1447619048 | 0 | 0 |
| 10 | `S_NH` | -0.0601904762 | 1 | -0.0601904762 |
| 14 | `S_PO4` | -0.0080634921 | 0 | 0 |
| 15 | `S_O2` | -0.5873015873 | 0 | 0 |

**Sum:** $B_{4,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{4,3}| = 0$ ≤ 0.1

<a id="rho5-P"></a>
#### rho5 — P (cell 29/114)

Balance for process row $i=4$ ($\rho_{5}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.016 + -0.0079365079 + 0 + 0 + -0.0080634921 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.016 | 0.016 |
| 6 | `S_S` | -1.5873015873 | 0.005 | -0.0079365079 |
| 8 | `S_IC` | 0.1447619048 | 0 | 0 |
| 10 | `S_NH` | -0.0601904762 | 0 | 0 |
| 14 | `S_PO4` | -0.0080634921 | 1 | -0.0080634921 |
| 15 | `S_O2` | -0.5873015873 | 0 | 0 |

**Sum:** $B_{4,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{4,4}| = 0$ ≤ 0.1

<a id="rho5-H"></a>
#### rho5 — H (cell 30/114)

Balance for process row $i=4$ ($\rho_{5}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.043 + 0 + 0 + -0.0132419048 + -0.0008063492 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.043 | 0.043 |
| 6 | `S_S` | -1.5873015873 | 0 | 0 |
| 8 | `S_IC` | 0.1447619048 | 0 | 0 |
| 10 | `S_NH` | -0.0601904762 | 0.22 | -0.0132419048 |
| 14 | `S_PO4` | -0.0080634921 | 0.1 | -0.0008063492 |
| 15 | `S_O2` | -0.5873015873 | 0 | 0 |

**Sum:** $B_{4,5} \approx 0.028951746$ (check: `numpy.sum` = 0.028951746)

**Status:** **OK** — $|B_{4,5}| = 0.028951746$ ≤ 0.1

### rho6

*Full label:* rho6: aerobic growth X_H (NO3-)

<a id="rho6-COD"></a>
#### rho6 — COD (cell 31/114)

Balance for process row $i=5$ ($\rho_{6}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 1 + -1.5873015873 + 0 + 0.2750704762 + 0 + 0.3121451247$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 1 | 1 |
| 6 | `S_S` | -1.5873015873 | 1 | -1.5873015873 |
| 8 | `S_IC` | 0.1447619048 | 0 | 0 |
| 12 | `S_NO3` | -0.0601904762 | -4.57 | 0.2750704762 |
| 14 | `S_PO4` | -0.0080634921 | 0 | 0 |
| 15 | `S_O2` | -0.3121451247 | -1 | 0.3121451247 |

**Sum:** $B_{5,0} \approx -8.598639e-05$ (check: `numpy.sum` = -8.598639e-05)

**Status:** **OK** — $|B_{5,0}| = 8.598639e-05$ ≤ 0.1

<a id="rho6-O"></a>
#### rho6 — O (cell 32/114)

Balance for process row $i=5$ ($\rho_{6}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.184 + -0.2476190476 + 0.3865142857 + -0.2064533333 + -0.0166914286 + -0.3121451247$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.184 | 0.184 |
| 6 | `S_S` | -1.5873015873 | 0.156 | -0.2476190476 |
| 8 | `S_IC` | 0.1447619048 | 2.67 | 0.3865142857 |
| 12 | `S_NO3` | -0.0601904762 | 3.43 | -0.2064533333 |
| 14 | `S_PO4` | -0.0080634921 | 2.07 | -0.0166914286 |
| 15 | `S_O2` | -0.3121451247 | 1 | -0.3121451247 |

**Sum:** $B_{5,1} \approx -0.2123946485$ (check: `numpy.sum` = -0.2123946485)

**Status:** **FAIL** — $|B_{5,1}| = 0.2123946485$ > 0.1

<a id="rho6-C"></a>
#### rho6 — C (cell 33/114)

Balance for process row $i=5$ ($\rho_{6}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.36 + -0.5047619048 + 0.1447619048 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.36 | 0.36 |
| 6 | `S_S` | -1.5873015873 | 0.318 | -0.5047619048 |
| 8 | `S_IC` | 0.1447619048 | 1 | 0.1447619048 |
| 12 | `S_NO3` | -0.0601904762 | 0 | 0 |
| 14 | `S_PO4` | -0.0080634921 | 0 | 0 |
| 15 | `S_O2` | -0.3121451247 | 0 | 0 |

**Sum:** $B_{5,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{5,2}| = 0$ ≤ 0.1

<a id="rho6-N"></a>
#### rho6 — N (cell 34/114)

Balance for process row $i=5$ ($\rho_{6}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.084 + -0.0238095238 + 0 + -0.0601904762 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.084 | 0.084 |
| 6 | `S_S` | -1.5873015873 | 0.015 | -0.0238095238 |
| 8 | `S_IC` | 0.1447619048 | 0 | 0 |
| 12 | `S_NO3` | -0.0601904762 | 1 | -0.0601904762 |
| 14 | `S_PO4` | -0.0080634921 | 0 | 0 |
| 15 | `S_O2` | -0.3121451247 | 0 | 0 |

**Sum:** $B_{5,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{5,3}| = 0$ ≤ 0.1

<a id="rho6-P"></a>
#### rho6 — P (cell 35/114)

Balance for process row $i=5$ ($\rho_{6}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.016 + -0.0079365079 + 0 + 0 + -0.0080634921 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.016 | 0.016 |
| 6 | `S_S` | -1.5873015873 | 0.005 | -0.0079365079 |
| 8 | `S_IC` | 0.1447619048 | 0 | 0 |
| 12 | `S_NO3` | -0.0601904762 | 0 | 0 |
| 14 | `S_PO4` | -0.0080634921 | 1 | -0.0080634921 |
| 15 | `S_O2` | -0.3121451247 | 0 | 0 |

**Sum:** $B_{5,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{5,4}| = 0$ ≤ 0.1

<a id="rho6-H"></a>
#### rho6 — H (cell 36/114)

Balance for process row $i=5$ ($\rho_{6}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.043 + 0 + 0 + -0.0042133333 + -0.0008063492 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.043 | 0.043 |
| 6 | `S_S` | -1.5873015873 | 0 | 0 |
| 8 | `S_IC` | 0.1447619048 | 0 | 0 |
| 12 | `S_NO3` | -0.0601904762 | 0.07 | -0.0042133333 |
| 14 | `S_PO4` | -0.0080634921 | 0.1 | -0.0008063492 |
| 15 | `S_O2` | -0.3121451247 | 0 | 0 |

**Sum:** $B_{5,5} \approx 0.0379803175$ (check: `numpy.sum` = 0.0379803175)

**Status:** **OK** — $|B_{5,5}| = 0.0379803175$ ≤ 0.1

### rho7

*Full label:* rho7: aerobic respiration X_H

<a id="rho7-COD"></a>
#### rho7 — COD (cell 37/114)

Balance for process row $i=6$ ($\rho_{7}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -1 + 0 + 0 + 0 + 1$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 1 | -1 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 0 | 0 |
| 15 | `S_O2` | -1 | -1 | 1 |

**Sum:** $B_{6,0} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{6,0}| = 0$ ≤ 0.1

<a id="rho7-O"></a>
#### rho7 — O (cell 38/114)

Balance for process row $i=6$ ($\rho_{7}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.184 + 0.9612 + 0 + 0.03312 + -1$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 0.184 | -0.184 |
| 8 | `S_IC` | 0.36 | 2.67 | 0.9612 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 2.07 | 0.03312 |
| 15 | `S_O2` | -1 | 1 | -1 |

**Sum:** $B_{6,1} \approx -0.18968$ (check: `numpy.sum` = -0.18968)

**Status:** **FAIL** — $|B_{6,1}| = 0.18968$ > 0.1

<a id="rho7-C"></a>
#### rho7 — C (cell 39/114)

Balance for process row $i=6$ ($\rho_{7}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.36 + 0.36 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 0.36 | -0.36 |
| 8 | `S_IC` | 0.36 | 1 | 0.36 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 0 | 0 |
| 15 | `S_O2` | -1 | 0 | 0 |

**Sum:** $B_{6,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{6,2}| = 0$ ≤ 0.1

<a id="rho7-N"></a>
#### rho7 — N (cell 40/114)

Balance for process row $i=6$ ($\rho_{7}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.084 + 0 + 0.084 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 0.084 | -0.084 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 1 | 0.084 |
| 14 | `S_PO4` | 0.016 | 0 | 0 |
| 15 | `S_O2` | -1 | 0 | 0 |

**Sum:** $B_{6,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{6,3}| = 0$ ≤ 0.1

<a id="rho7-P"></a>
#### rho7 — P (cell 41/114)

Balance for process row $i=6$ ($\rho_{7}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.016 + 0 + 0 + 0.016 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 0.016 | -0.016 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 1 | 0.016 |
| 15 | `S_O2` | -1 | 0 | 0 |

**Sum:** $B_{6,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{6,4}| = 0$ ≤ 0.1

<a id="rho7-H"></a>
#### rho7 — H (cell 42/114)

Balance for process row $i=6$ ($\rho_{7}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.043 + 0 + 0.01848 + 0.0016 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 0.043 | -0.043 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 0.22 | 0.01848 |
| 14 | `S_PO4` | 0.016 | 0.1 | 0.0016 |
| 15 | `S_O2` | -1 | 0 | 0 |

**Sum:** $B_{6,5} \approx -0.02292$ (check: `numpy.sum` = -0.02292)

**Status:** **OK** — $|B_{6,5}| = 0.02292$ ≤ 0.1

### rho8

*Full label:* rho8: anoxic growth X_H (NO3-)

<a id="rho8-COD"></a>
#### rho8 — COD (cell 43/114)

Balance for process row $i=7$ ($\rho_{8}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 1 + -2 + 0 + 0 + 1.5995 + -0.5985 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 1 | 1 |
| 6 | `S_S` | -2 | 1 | -2 |
| 8 | `S_IC` | 0.276 | 0 | 0 |
| 10 | `S_NH` | -0.054 | 0 | 0 |
| 12 | `S_NO3` | -0.35 | -4.57 | 1.5995 |
| 13 | `S_N2` | 0.35 | -1.71 | -0.5985 |
| 14 | `S_PO4` | -0.006 | 0 | 0 |

**Sum:** $B_{7,0} \approx 0.001$ (check: `numpy.sum` = 0.001)

**Status:** **OK** — $|B_{7,0}| = 0.001$ ≤ 0.1

<a id="rho8-O"></a>
#### rho8 — O (cell 44/114)

Balance for process row $i=7$ ($\rho_{8}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.184 + -0.312 + 0.73692 + 0 + -1.2005 + 0 + -0.01242$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.184 | 0.184 |
| 6 | `S_S` | -2 | 0.156 | -0.312 |
| 8 | `S_IC` | 0.276 | 2.67 | 0.73692 |
| 10 | `S_NH` | -0.054 | 0 | 0 |
| 12 | `S_NO3` | -0.35 | 3.43 | -1.2005 |
| 13 | `S_N2` | 0.35 | 0 | 0 |
| 14 | `S_PO4` | -0.006 | 2.07 | -0.01242 |

**Sum:** $B_{7,1} \approx -0.604$ (check: `numpy.sum` = -0.604)

**Status:** **FAIL** — $|B_{7,1}| = 0.604$ > 0.1

<a id="rho8-C"></a>
#### rho8 — C (cell 45/114)

Balance for process row $i=7$ ($\rho_{8}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.36 + -0.636 + 0.276 + 0 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.36 | 0.36 |
| 6 | `S_S` | -2 | 0.318 | -0.636 |
| 8 | `S_IC` | 0.276 | 1 | 0.276 |
| 10 | `S_NH` | -0.054 | 0 | 0 |
| 12 | `S_NO3` | -0.35 | 0 | 0 |
| 13 | `S_N2` | 0.35 | 0 | 0 |
| 14 | `S_PO4` | -0.006 | 0 | 0 |

**Sum:** $B_{7,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{7,2}| = 0$ ≤ 0.1

<a id="rho8-N"></a>
#### rho8 — N (cell 46/114)

Balance for process row $i=7$ ($\rho_{8}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.084 + -0.03 + 0 + -0.054 + -0.35 + 0.35 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.084 | 0.084 |
| 6 | `S_S` | -2 | 0.015 | -0.03 |
| 8 | `S_IC` | 0.276 | 0 | 0 |
| 10 | `S_NH` | -0.054 | 1 | -0.054 |
| 12 | `S_NO3` | -0.35 | 1 | -0.35 |
| 13 | `S_N2` | 0.35 | 1 | 0.35 |
| 14 | `S_PO4` | -0.006 | 0 | 0 |

**Sum:** $B_{7,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{7,3}| = 0$ ≤ 0.1

<a id="rho8-P"></a>
#### rho8 — P (cell 47/114)

Balance for process row $i=7$ ($\rho_{8}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.016 + -0.01 + 0 + 0 + 0 + 0 + -0.006$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.016 | 0.016 |
| 6 | `S_S` | -2 | 0.005 | -0.01 |
| 8 | `S_IC` | 0.276 | 0 | 0 |
| 10 | `S_NH` | -0.054 | 0 | 0 |
| 12 | `S_NO3` | -0.35 | 0 | 0 |
| 13 | `S_N2` | 0.35 | 0 | 0 |
| 14 | `S_PO4` | -0.006 | 1 | -0.006 |

**Sum:** $B_{7,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{7,4}| = 0$ ≤ 0.1

<a id="rho8-H"></a>
#### rho8 — H (cell 48/114)

Balance for process row $i=7$ ($\rho_{8}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.043 + 0 + 0 + -0.01188 + -0.0245 + 0 + -0.0006$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.043 | 0.043 |
| 6 | `S_S` | -2 | 0 | 0 |
| 8 | `S_IC` | 0.276 | 0 | 0 |
| 10 | `S_NH` | -0.054 | 0.22 | -0.01188 |
| 12 | `S_NO3` | -0.35 | 0.07 | -0.0245 |
| 13 | `S_N2` | 0.35 | 0 | 0 |
| 14 | `S_PO4` | -0.006 | 0.1 | -0.0006 |

**Sum:** $B_{7,5} \approx 0.00602$ (check: `numpy.sum` = 0.00602)

**Status:** **OK** — $|B_{7,5}| = 0.00602$ ≤ 0.1

### rho9

*Full label:* rho9: anoxic growth X_H (NO2-)

<a id="rho9-COD"></a>
#### rho9 — COD (cell 49/114)

Balance for process row $i=8$ ($\rho_{9}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 1 + -3.3333333333 + 0 + 0 + 4.6686111111 + -2.3275 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 1 | 1 |
| 6 | `S_S` | -3.3333333333 | 1 | -3.3333333333 |
| 8 | `S_IC` | 0.7 | 0 | 0 |
| 10 | `S_NH` | -0.034 | 0 | 0 |
| 11 | `S_NO2` | -1.3611111111 | -3.43 | 4.6686111111 |
| 13 | `S_N2` | 1.3611111111 | -1.71 | -2.3275 |
| 14 | `S_PO4` | 0.0006666667 | 0 | 0 |

**Sum:** $B_{8,0} \approx 0.0077777778$ (check: `numpy.sum` = 0.0077777778)

**Status:** **OK** — $|B_{8,0}| = 0.0077777778$ ≤ 0.1

<a id="rho9-O"></a>
#### rho9 — O (cell 50/114)

Balance for process row $i=8$ ($\rho_{9}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.184 + -0.52 + 1.869 + 0 + -3.1033333333 + 0 + 0.00138$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.184 | 0.184 |
| 6 | `S_S` | -3.3333333333 | 0.156 | -0.52 |
| 8 | `S_IC` | 0.7 | 2.67 | 1.869 |
| 10 | `S_NH` | -0.034 | 0 | 0 |
| 11 | `S_NO2` | -1.3611111111 | 2.28 | -3.1033333333 |
| 13 | `S_N2` | 1.3611111111 | 0 | 0 |
| 14 | `S_PO4` | 0.0006666667 | 2.07 | 0.00138 |

**Sum:** $B_{8,1} \approx -1.5689533333$ (check: `numpy.sum` = -1.5689533333)

**Status:** **FAIL** — $|B_{8,1}| = 1.5689533333$ > 0.1

<a id="rho9-C"></a>
#### rho9 — C (cell 51/114)

Balance for process row $i=8$ ($\rho_{9}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.36 + -1.06 + 0.7 + 0 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.36 | 0.36 |
| 6 | `S_S` | -3.3333333333 | 0.318 | -1.06 |
| 8 | `S_IC` | 0.7 | 1 | 0.7 |
| 10 | `S_NH` | -0.034 | 0 | 0 |
| 11 | `S_NO2` | -1.3611111111 | 0 | 0 |
| 13 | `S_N2` | 1.3611111111 | 0 | 0 |
| 14 | `S_PO4` | 0.0006666667 | 0 | 0 |

**Sum:** $B_{8,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{8,2}| = 0$ ≤ 0.1

<a id="rho9-N"></a>
#### rho9 — N (cell 52/114)

Balance for process row $i=8$ ($\rho_{9}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.084 + -0.05 + 0 + -0.034 + -1.3611111111 + 1.3611111111 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.084 | 0.084 |
| 6 | `S_S` | -3.3333333333 | 0.015 | -0.05 |
| 8 | `S_IC` | 0.7 | 0 | 0 |
| 10 | `S_NH` | -0.034 | 1 | -0.034 |
| 11 | `S_NO2` | -1.3611111111 | 1 | -1.3611111111 |
| 13 | `S_N2` | 1.3611111111 | 1 | 1.3611111111 |
| 14 | `S_PO4` | 0.0006666667 | 0 | 0 |

**Sum:** $B_{8,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{8,3}| = 0$ ≤ 0.1

<a id="rho9-P"></a>
#### rho9 — P (cell 53/114)

Balance for process row $i=8$ ($\rho_{9}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.016 + -0.0166666667 + 0 + 0 + 0 + 0 + 0.0006666667$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.016 | 0.016 |
| 6 | `S_S` | -3.3333333333 | 0.005 | -0.0166666667 |
| 8 | `S_IC` | 0.7 | 0 | 0 |
| 10 | `S_NH` | -0.034 | 0 | 0 |
| 11 | `S_NO2` | -1.3611111111 | 0 | 0 |
| 13 | `S_N2` | 1.3611111111 | 0 | 0 |
| 14 | `S_PO4` | 0.0006666667 | 1 | 0.0006666667 |

**Sum:** $B_{8,4} \approx -1.301043e-18$ (check: `numpy.sum` = -1.301043e-18)

**Status:** **OK** — $|B_{8,4}| = 1.301043e-18$ ≤ 0.1

<a id="rho9-H"></a>
#### rho9 — H (cell 54/114)

Balance for process row $i=8$ ($\rho_{9}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.043 + 0 + 0 + -0.00748 + -0.0952777778 + 0 + 6.666667e-05$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | 1 | 0.043 | 0.043 |
| 6 | `S_S` | -3.3333333333 | 0 | 0 |
| 8 | `S_IC` | 0.7 | 0 | 0 |
| 10 | `S_NH` | -0.034 | 0.22 | -0.00748 |
| 11 | `S_NO2` | -1.3611111111 | 0.07 | -0.0952777778 |
| 13 | `S_N2` | 1.3611111111 | 0 | 0 |
| 14 | `S_PO4` | 0.0006666667 | 0.1 | 6.666667e-05 |

**Sum:** $B_{8,5} \approx -0.0596911111$ (check: `numpy.sum` = -0.0596911111)

**Status:** **OK** — $|B_{8,5}| = 0.0596911111$ ≤ 0.1

### rho10

*Full label:* rho10: anoxic respiration X_H

<a id="rho10-COD"></a>
#### rho10 — COD (cell 55/114)

Balance for process row $i=9$ ($\rho_{10}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -1 + 0 + 0 + 0.7503125 + 0.9996875 + -0.748125 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 1 | -1 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 11 | `S_NO2` | -0.21875 | -3.43 | 0.7503125 |
| 12 | `S_NO3` | -0.21875 | -4.57 | 0.9996875 |
| 13 | `S_N2` | 0.4375 | -1.71 | -0.748125 |
| 14 | `S_PO4` | 0.016 | 0 | 0 |

**Sum:** $B_{9,0} \approx 0.001875$ (check: `numpy.sum` = 0.001875)

**Status:** **OK** — $|B_{9,0}| = 0.001875$ ≤ 0.1

<a id="rho10-O"></a>
#### rho10 — O (cell 56/114)

Balance for process row $i=9$ ($\rho_{10}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.184 + 0.9612 + 0 + -0.49875 + -0.7503125 + 0 + 0.03312$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 0.184 | -0.184 |
| 8 | `S_IC` | 0.36 | 2.67 | 0.9612 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 11 | `S_NO2` | -0.21875 | 2.28 | -0.49875 |
| 12 | `S_NO3` | -0.21875 | 3.43 | -0.7503125 |
| 13 | `S_N2` | 0.4375 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 2.07 | 0.03312 |

**Sum:** $B_{9,1} \approx -0.4387425$ (check: `numpy.sum` = -0.4387425)

**Status:** **FAIL** — $|B_{9,1}| = 0.4387425$ > 0.1

<a id="rho10-C"></a>
#### rho10 — C (cell 57/114)

Balance for process row $i=9$ ($\rho_{10}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.36 + 0.36 + 0 + 0 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 0.36 | -0.36 |
| 8 | `S_IC` | 0.36 | 1 | 0.36 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 11 | `S_NO2` | -0.21875 | 0 | 0 |
| 12 | `S_NO3` | -0.21875 | 0 | 0 |
| 13 | `S_N2` | 0.4375 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 0 | 0 |

**Sum:** $B_{9,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{9,2}| = 0$ ≤ 0.1

<a id="rho10-N"></a>
#### rho10 — N (cell 58/114)

Balance for process row $i=9$ ($\rho_{10}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.084 + 0 + 0.084 + -0.21875 + -0.21875 + 0.4375 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 0.084 | -0.084 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 1 | 0.084 |
| 11 | `S_NO2` | -0.21875 | 1 | -0.21875 |
| 12 | `S_NO3` | -0.21875 | 1 | -0.21875 |
| 13 | `S_N2` | 0.4375 | 1 | 0.4375 |
| 14 | `S_PO4` | 0.016 | 0 | 0 |

**Sum:** $B_{9,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{9,3}| = 0$ ≤ 0.1

<a id="rho10-P"></a>
#### rho10 — P (cell 59/114)

Balance for process row $i=9$ ($\rho_{10}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.016 + 0 + 0 + 0 + 0 + 0 + 0.016$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 0.016 | -0.016 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 11 | `S_NO2` | -0.21875 | 0 | 0 |
| 12 | `S_NO3` | -0.21875 | 0 | 0 |
| 13 | `S_N2` | 0.4375 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 1 | 0.016 |

**Sum:** $B_{9,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{9,4}| = 0$ ≤ 0.1

<a id="rho10-H"></a>
#### rho10 — H (cell 60/114)

Balance for process row $i=9$ ($\rho_{10}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.043 + 0 + 0.01848 + -0.0153125 + -0.0153125 + 0 + 0.0016$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 0.043 | -0.043 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 0.22 | 0.01848 |
| 11 | `S_NO2` | -0.21875 | 0.07 | -0.0153125 |
| 12 | `S_NO3` | -0.21875 | 0.07 | -0.0153125 |
| 13 | `S_N2` | 0.4375 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 0.1 | 0.0016 |

**Sum:** $B_{9,5} \approx -0.053545$ (check: `numpy.sum` = -0.053545)

**Status:** **OK** — $|B_{9,5}| = 0.053545$ ≤ 0.1

### rho11

*Full label:* rho11: hydrolysis X_S

<a id="rho11-COD"></a>
#### rho11 — COD (cell 61/114)

Balance for process row $i=10$ ($\rho_{11}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -1 + 0.9 + 0.1 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 4 | `X_S` | -1 | 1 | -1 |
| 6 | `S_S` | 0.9 | 1 | 0.9 |
| 7 | `S_I` | 0.1 | 1 | 0.1 |
| 8 | `S_IC` | -0.0042 | 0 | 0 |
| 10 | `S_NH` | 0.0145 | 0 | 0 |
| 14 | `S_PO4` | -4.336809e-19 | 0 | 0 |

**Sum:** $B_{10,0} \approx 2.775558e-17$ (check: `numpy.sum` = 2.775558e-17)

**Status:** **OK** — $|B_{10,0}| = 2.775558e-17$ ≤ 0.1

<a id="rho11-O"></a>
#### rho11 — O (cell 62/114)

Balance for process row $i=10$ ($\rho_{11}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.156 + 0.1404 + 0.015 + -0.011214 + 0 + -8.977194e-19$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 4 | `X_S` | -1 | 0.156 | -0.156 |
| 6 | `S_S` | 0.9 | 0.156 | 0.1404 |
| 7 | `S_I` | 0.1 | 0.15 | 0.015 |
| 8 | `S_IC` | -0.0042 | 2.67 | -0.011214 |
| 10 | `S_NH` | 0.0145 | 0 | 0 |
| 14 | `S_PO4` | -4.336809e-19 | 2.07 | -8.977194e-19 |

**Sum:** $B_{10,1} \approx -0.011814$ (check: `numpy.sum` = -0.011814)

**Status:** **OK** — $|B_{10,1}| = 0.011814$ ≤ 0.1

<a id="rho11-C"></a>
#### rho11 — C (cell 63/114)

Balance for process row $i=10$ ($\rho_{11}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.318 + 0.2862 + 0.036 + -0.0042 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 4 | `X_S` | -1 | 0.318 | -0.318 |
| 6 | `S_S` | 0.9 | 0.318 | 0.2862 |
| 7 | `S_I` | 0.1 | 0.36 | 0.036 |
| 8 | `S_IC` | -0.0042 | 1 | -0.0042 |
| 10 | `S_NH` | 0.0145 | 0 | 0 |
| 14 | `S_PO4` | -4.336809e-19 | 0 | 0 |

**Sum:** $B_{10,2} \approx 3.469447e-18$ (check: `numpy.sum` = 3.469447e-18)

**Status:** **OK** — $|B_{10,2}| = 3.469447e-18$ ≤ 0.1

<a id="rho11-N"></a>
#### rho11 — N (cell 64/114)

Balance for process row $i=10$ ($\rho_{11}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.034 + 0.0135 + 0.006 + 0 + 0.0145 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 4 | `X_S` | -1 | 0.034 | -0.034 |
| 6 | `S_S` | 0.9 | 0.015 | 0.0135 |
| 7 | `S_I` | 0.1 | 0.06 | 0.006 |
| 8 | `S_IC` | -0.0042 | 0 | 0 |
| 10 | `S_NH` | 0.0145 | 1 | 0.0145 |
| 14 | `S_PO4` | -4.336809e-19 | 0 | 0 |

**Sum:** $B_{10,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{10,3}| = 0$ ≤ 0.1

<a id="rho11-P"></a>
#### rho11 — P (cell 65/114)

Balance for process row $i=10$ ($\rho_{11}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.005 + 0.0045 + 0.0005 + 0 + 0 + -4.336809e-19$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 4 | `X_S` | -1 | 0.005 | -0.005 |
| 6 | `S_S` | 0.9 | 0.005 | 0.0045 |
| 7 | `S_I` | 0.1 | 0.005 | 0.0005 |
| 8 | `S_IC` | -0.0042 | 0 | 0 |
| 10 | `S_NH` | 0.0145 | 0 | 0 |
| 14 | `S_PO4` | -4.336809e-19 | 1 | -4.336809e-19 |

**Sum:** $B_{10,4} \approx -2.975051e-19$ (check: `numpy.sum` = -2.975051e-19)

**Status:** **OK** — $|B_{10,4}| = 2.975051e-19$ ≤ 0.1

<a id="rho11-H"></a>
#### rho11 — H (cell 66/114)

Balance for process row $i=10$ ($\rho_{11}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0 + 0 + 0 + 0 + 0.00319 + -4.336809e-20$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 4 | `X_S` | -1 | 0 | 0 |
| 6 | `S_S` | 0.9 | 0 | 0 |
| 7 | `S_I` | 0.1 | 0 | 0 |
| 8 | `S_IC` | -0.0042 | 0 | 0 |
| 10 | `S_NH` | 0.0145 | 0.22 | 0.00319 |
| 14 | `S_PO4` | -4.336809e-19 | 0.1 | -4.336809e-20 |

**Sum:** $B_{10,5} \approx 0.00319$ (check: `numpy.sum` = 0.00319)

**Status:** **OK** — $|B_{10,5}| = 0.00319$ ≤ 0.1

### rho12

*Full label:* rho12: urea hydrolysis

<a id="rho12-COD"></a>
#### rho12 — COD (cell 67/114)

Balance for process row $i=11$ ($\rho_{12}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 8 | `S_IC` | 0.429 | 0 | 0 |
| 9 | `S_ND` | -1 | 0 | 0 |
| 10 | `S_NH` | 1 | 0 | 0 |
| 16 | `S_H2O` | 0.072 | 0 | 0 |

**Sum:** $B_{11,0} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{11,0}| = 0$ ≤ 0.1

<a id="rho12-O"></a>
#### rho12 — O (cell 68/114)

Balance for process row $i=11$ ($\rho_{12}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 1.14543 + -0.57 + 0 + 0.57168$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 8 | `S_IC` | 0.429 | 2.67 | 1.14543 |
| 9 | `S_ND` | -1 | 0.57 | -0.57 |
| 10 | `S_NH` | 1 | 0 | 0 |
| 16 | `S_H2O` | 0.072 | 7.94 | 0.57168 |

**Sum:** $B_{11,1} \approx 1.14711$ (check: `numpy.sum` = 1.14711)

**Status:** **FAIL** — $|B_{11,1}| = 1.14711$ > 0.1

<a id="rho12-C"></a>
#### rho12 — C (cell 69/114)

Balance for process row $i=11$ ($\rho_{12}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.429 + -0.43 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 8 | `S_IC` | 0.429 | 1 | 0.429 |
| 9 | `S_ND` | -1 | 0.43 | -0.43 |
| 10 | `S_NH` | 1 | 0 | 0 |
| 16 | `S_H2O` | 0.072 | 0 | 0 |

**Sum:** $B_{11,2} \approx -0.001$ (check: `numpy.sum` = -0.001)

**Status:** **OK** — $|B_{11,2}| = 0.001$ ≤ 0.1

<a id="rho12-N"></a>
#### rho12 — N (cell 70/114)

Balance for process row $i=11$ ($\rho_{12}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0 + -1 + 1 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 8 | `S_IC` | 0.429 | 0 | 0 |
| 9 | `S_ND` | -1 | 1 | -1 |
| 10 | `S_NH` | 1 | 1 | 1 |
| 16 | `S_H2O` | 0.072 | 0 | 0 |

**Sum:** $B_{11,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{11,3}| = 0$ ≤ 0.1

<a id="rho12-P"></a>
#### rho12 — P (cell 71/114)

Balance for process row $i=11$ ($\rho_{12}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 8 | `S_IC` | 0.429 | 0 | 0 |
| 9 | `S_ND` | -1 | 0 | 0 |
| 10 | `S_NH` | 1 | 0 | 0 |
| 16 | `S_H2O` | 0.072 | 0 | 0 |

**Sum:** $B_{11,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{11,4}| = 0$ ≤ 0.1

<a id="rho12-H"></a>
#### rho12 — H (cell 72/114)

Balance for process row $i=11$ ($\rho_{12}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0 + -0.14 + 0.22 + 0.072$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 8 | `S_IC` | 0.429 | 0 | 0 |
| 9 | `S_ND` | -1 | 0.14 | -0.14 |
| 10 | `S_NH` | 1 | 0.22 | 0.22 |
| 16 | `S_H2O` | 0.072 | 1 | 0.072 |

**Sum:** $B_{11,5} \approx 0.152$ (check: `numpy.sum` = 0.152)

**Status:** **FAIL** — $|B_{11,5}| = 0.152$ > 0.1

### rho13

*Full label:* rho13: decay X_H

<a id="rho13-COD"></a>
#### rho13 — COD (cell 73/114)

Balance for process row $i=12$ ($\rho_{13}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -1 + 0.9 + 0.1 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 1 | -1 |
| 4 | `X_S` | 0.9 | 1 | 0.9 |
| 5 | `X_I` | 0.1 | 1 | 0.1 |
| 8 | `S_IC` | 0.0378 | 0 | 0 |
| 10 | `S_NH` | 0.0474 | 0 | 0 |
| 14 | `S_PO4` | 0.0105 | 0 | 0 |

**Sum:** $B_{12,0} \approx 2.775558e-17$ (check: `numpy.sum` = 2.775558e-17)

**Status:** **OK** — $|B_{12,0}| = 2.775558e-17$ ≤ 0.1

<a id="rho13-O"></a>
#### rho13 — O (cell 74/114)

Balance for process row $i=12$ ($\rho_{13}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.184 + 0.1404 + 0.015 + 0.100926 + 0 + 0.021735$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 0.184 | -0.184 |
| 4 | `X_S` | 0.9 | 0.156 | 0.1404 |
| 5 | `X_I` | 0.1 | 0.15 | 0.015 |
| 8 | `S_IC` | 0.0378 | 2.67 | 0.100926 |
| 10 | `S_NH` | 0.0474 | 0 | 0 |
| 14 | `S_PO4` | 0.0105 | 2.07 | 0.021735 |

**Sum:** $B_{12,1} \approx 0.094061$ (check: `numpy.sum` = 0.094061)

**Status:** **OK** — $|B_{12,1}| = 0.094061$ ≤ 0.1

<a id="rho13-C"></a>
#### rho13 — C (cell 75/114)

Balance for process row $i=12$ ($\rho_{13}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.36 + 0.2862 + 0.036 + 0.0378 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 0.36 | -0.36 |
| 4 | `X_S` | 0.9 | 0.318 | 0.2862 |
| 5 | `X_I` | 0.1 | 0.36 | 0.036 |
| 8 | `S_IC` | 0.0378 | 1 | 0.0378 |
| 10 | `S_NH` | 0.0474 | 0 | 0 |
| 14 | `S_PO4` | 0.0105 | 0 | 0 |

**Sum:** $B_{12,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{12,2}| = 0$ ≤ 0.1

<a id="rho13-N"></a>
#### rho13 — N (cell 76/114)

Balance for process row $i=12$ ($\rho_{13}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.084 + 0.0306 + 0.006 + 0 + 0.0474 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 0.084 | -0.084 |
| 4 | `X_S` | 0.9 | 0.034 | 0.0306 |
| 5 | `X_I` | 0.1 | 0.06 | 0.006 |
| 8 | `S_IC` | 0.0378 | 0 | 0 |
| 10 | `S_NH` | 0.0474 | 1 | 0.0474 |
| 14 | `S_PO4` | 0.0105 | 0 | 0 |

**Sum:** $B_{12,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{12,3}| = 0$ ≤ 0.1

<a id="rho13-P"></a>
#### rho13 — P (cell 77/114)

Balance for process row $i=12$ ($\rho_{13}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.016 + 0.0045 + 0.001 + 0 + 0 + 0.0105$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 0.016 | -0.016 |
| 4 | `X_S` | 0.9 | 0.005 | 0.0045 |
| 5 | `X_I` | 0.1 | 0.01 | 0.001 |
| 8 | `S_IC` | 0.0378 | 0 | 0 |
| 10 | `S_NH` | 0.0474 | 0 | 0 |
| 14 | `S_PO4` | 0.0105 | 1 | 0.0105 |

**Sum:** $B_{12,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{12,4}| = 0$ ≤ 0.1

<a id="rho13-H"></a>
#### rho13 — H (cell 78/114)

Balance for process row $i=12$ ($\rho_{13}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.043 + 0 + 0 + 0 + 0.010428 + 0.00105$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 3 | `X_H` | -1 | 0.043 | -0.043 |
| 4 | `X_S` | 0.9 | 0 | 0 |
| 5 | `X_I` | 0.1 | 0 | 0 |
| 8 | `S_IC` | 0.0378 | 0 | 0 |
| 10 | `S_NH` | 0.0474 | 0.22 | 0.010428 |
| 14 | `S_PO4` | 0.0105 | 0.1 | 0.00105 |

**Sum:** $B_{12,5} \approx -0.031522$ (check: `numpy.sum` = -0.031522)

**Status:** **OK** — $|B_{12,5}| = 0.031522$ ≤ 0.1

### rho14

*Full label:* rho14: aerobic growth X_AOB

<a id="rho14-COD"></a>
#### rho14 — COD (cell 79/114)

Balance for process row $i=13$ ($\rho_{14}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 1 + 0 + 0 + -17.15 + 0 + 16.1428571429$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | 1 | 1 | 1 |
| 8 | `S_IC` | -0.36 | 0 | 0 |
| 10 | `S_NH` | -5.084 | 0 | 0 |
| 11 | `S_NO2` | 5 | -3.43 | -17.15 |
| 14 | `S_PO4` | -0.016 | 0 | 0 |
| 15 | `S_O2` | -16.1428571429 | -1 | 16.1428571429 |

**Sum:** $B_{13,0} \approx -0.0071428571$ (check: `numpy.sum` = -0.0071428571)

**Status:** **OK** — $|B_{13,0}| = 0.0071428571$ ≤ 0.1

<a id="rho14-O"></a>
#### rho14 — O (cell 80/114)

Balance for process row $i=13$ ($\rho_{14}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.184 + -0.9612 + 0 + 11.4 + -0.03312 + -16.1428571429$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | 1 | 0.184 | 0.184 |
| 8 | `S_IC` | -0.36 | 2.67 | -0.9612 |
| 10 | `S_NH` | -5.084 | 0 | 0 |
| 11 | `S_NO2` | 5 | 2.28 | 11.4 |
| 14 | `S_PO4` | -0.016 | 2.07 | -0.03312 |
| 15 | `S_O2` | -16.1428571429 | 1 | -16.1428571429 |

**Sum:** $B_{13,1} \approx -5.5531771429$ (check: `numpy.sum` = -5.5531771429)

**Status:** **FAIL** — $|B_{13,1}| = 5.5531771429$ > 0.1

<a id="rho14-C"></a>
#### rho14 — C (cell 81/114)

Balance for process row $i=13$ ($\rho_{14}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.36 + -0.36 + 0 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | 1 | 0.36 | 0.36 |
| 8 | `S_IC` | -0.36 | 1 | -0.36 |
| 10 | `S_NH` | -5.084 | 0 | 0 |
| 11 | `S_NO2` | 5 | 0 | 0 |
| 14 | `S_PO4` | -0.016 | 0 | 0 |
| 15 | `S_O2` | -16.1428571429 | 0 | 0 |

**Sum:** $B_{13,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{13,2}| = 0$ ≤ 0.1

<a id="rho14-N"></a>
#### rho14 — N (cell 82/114)

Balance for process row $i=13$ ($\rho_{14}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.084 + 0 + -5.084 + 5 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | 1 | 0.084 | 0.084 |
| 8 | `S_IC` | -0.36 | 0 | 0 |
| 10 | `S_NH` | -5.084 | 1 | -5.084 |
| 11 | `S_NO2` | 5 | 1 | 5 |
| 14 | `S_PO4` | -0.016 | 0 | 0 |
| 15 | `S_O2` | -16.1428571429 | 0 | 0 |

**Sum:** $B_{13,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{13,3}| = 0$ ≤ 0.1

<a id="rho14-P"></a>
#### rho14 — P (cell 83/114)

Balance for process row $i=13$ ($\rho_{14}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.016 + 0 + 0 + 0 + -0.016 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | 1 | 0.016 | 0.016 |
| 8 | `S_IC` | -0.36 | 0 | 0 |
| 10 | `S_NH` | -5.084 | 0 | 0 |
| 11 | `S_NO2` | 5 | 0 | 0 |
| 14 | `S_PO4` | -0.016 | 1 | -0.016 |
| 15 | `S_O2` | -16.1428571429 | 0 | 0 |

**Sum:** $B_{13,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{13,4}| = 0$ ≤ 0.1

<a id="rho14-H"></a>
#### rho14 — H (cell 84/114)

Balance for process row $i=13$ ($\rho_{14}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.043 + 0 + -1.11848 + 0.35 + -0.0016 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | 1 | 0.043 | 0.043 |
| 8 | `S_IC` | -0.36 | 0 | 0 |
| 10 | `S_NH` | -5.084 | 0.22 | -1.11848 |
| 11 | `S_NO2` | 5 | 0.07 | 0.35 |
| 14 | `S_PO4` | -0.016 | 0.1 | -0.0016 |
| 15 | `S_O2` | -16.1428571429 | 0 | 0 |

**Sum:** $B_{13,5} \approx -0.72708$ (check: `numpy.sum` = -0.72708)

**Status:** **FAIL** — $|B_{13,5}| = 0.72708$ > 0.1

### rho15

*Full label:* rho15: aerobic respiration X_AOB

<a id="rho15-COD"></a>
#### rho15 — COD (cell 85/114)

Balance for process row $i=14$ ($\rho_{15}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -1 + 0 + 0 + 0 + 1$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | -1 | 1 | -1 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 0 | 0 |
| 15 | `S_O2` | -1 | -1 | 1 |

**Sum:** $B_{14,0} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{14,0}| = 0$ ≤ 0.1

<a id="rho15-O"></a>
#### rho15 — O (cell 86/114)

Balance for process row $i=14$ ($\rho_{15}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.184 + 0.9612 + 0 + 0.03312 + -1$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | -1 | 0.184 | -0.184 |
| 8 | `S_IC` | 0.36 | 2.67 | 0.9612 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 2.07 | 0.03312 |
| 15 | `S_O2` | -1 | 1 | -1 |

**Sum:** $B_{14,1} \approx -0.18968$ (check: `numpy.sum` = -0.18968)

**Status:** **FAIL** — $|B_{14,1}| = 0.18968$ > 0.1

<a id="rho15-C"></a>
#### rho15 — C (cell 87/114)

Balance for process row $i=14$ ($\rho_{15}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.36 + 0.36 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | -1 | 0.36 | -0.36 |
| 8 | `S_IC` | 0.36 | 1 | 0.36 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 0 | 0 |
| 15 | `S_O2` | -1 | 0 | 0 |

**Sum:** $B_{14,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{14,2}| = 0$ ≤ 0.1

<a id="rho15-N"></a>
#### rho15 — N (cell 88/114)

Balance for process row $i=14$ ($\rho_{15}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.084 + 0 + 0.084 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | -1 | 0.084 | -0.084 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 1 | 0.084 |
| 14 | `S_PO4` | 0.016 | 0 | 0 |
| 15 | `S_O2` | -1 | 0 | 0 |

**Sum:** $B_{14,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{14,3}| = 0$ ≤ 0.1

<a id="rho15-P"></a>
#### rho15 — P (cell 89/114)

Balance for process row $i=14$ ($\rho_{15}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.016 + 0 + 0 + 0.016 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | -1 | 0.016 | -0.016 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 1 | 0.016 |
| 15 | `S_O2` | -1 | 0 | 0 |

**Sum:** $B_{14,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{14,4}| = 0$ ≤ 0.1

<a id="rho15-H"></a>
#### rho15 — H (cell 90/114)

Balance for process row $i=14$ ($\rho_{15}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.043 + 0 + 0.01848 + 0.0016 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | -1 | 0.043 | -0.043 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 0.22 | 0.01848 |
| 14 | `S_PO4` | 0.016 | 0.1 | 0.0016 |
| 15 | `S_O2` | -1 | 0 | 0 |

**Sum:** $B_{14,5} \approx -0.02292$ (check: `numpy.sum` = -0.02292)

**Status:** **OK** — $|B_{14,5}| = 0.02292$ ≤ 0.1

### rho16

*Full label:* rho16: decay X_AOB

<a id="rho16-COD"></a>
#### rho16 — COD (cell 91/114)

Balance for process row $i=15$ ($\rho_{16}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -1 + 0.9 + 0.1 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | -1 | 1 | -1 |
| 4 | `X_S` | 0.9 | 1 | 0.9 |
| 5 | `X_I` | 0.1 | 1 | 0.1 |
| 8 | `S_IC` | 0.0378 | 0 | 0 |
| 10 | `S_NH` | 0.0474 | 0 | 0 |
| 14 | `S_PO4` | 0.0105 | 0 | 0 |

**Sum:** $B_{15,0} \approx 2.775558e-17$ (check: `numpy.sum` = 2.775558e-17)

**Status:** **OK** — $|B_{15,0}| = 2.775558e-17$ ≤ 0.1

<a id="rho16-O"></a>
#### rho16 — O (cell 92/114)

Balance for process row $i=15$ ($\rho_{16}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.184 + 0.1404 + 0.015 + 0.100926 + 0 + 0.021735$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | -1 | 0.184 | -0.184 |
| 4 | `X_S` | 0.9 | 0.156 | 0.1404 |
| 5 | `X_I` | 0.1 | 0.15 | 0.015 |
| 8 | `S_IC` | 0.0378 | 2.67 | 0.100926 |
| 10 | `S_NH` | 0.0474 | 0 | 0 |
| 14 | `S_PO4` | 0.0105 | 2.07 | 0.021735 |

**Sum:** $B_{15,1} \approx 0.094061$ (check: `numpy.sum` = 0.094061)

**Status:** **OK** — $|B_{15,1}| = 0.094061$ ≤ 0.1

<a id="rho16-C"></a>
#### rho16 — C (cell 93/114)

Balance for process row $i=15$ ($\rho_{16}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.36 + 0.2862 + 0.036 + 0.0378 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | -1 | 0.36 | -0.36 |
| 4 | `X_S` | 0.9 | 0.318 | 0.2862 |
| 5 | `X_I` | 0.1 | 0.36 | 0.036 |
| 8 | `S_IC` | 0.0378 | 1 | 0.0378 |
| 10 | `S_NH` | 0.0474 | 0 | 0 |
| 14 | `S_PO4` | 0.0105 | 0 | 0 |

**Sum:** $B_{15,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{15,2}| = 0$ ≤ 0.1

<a id="rho16-N"></a>
#### rho16 — N (cell 94/114)

Balance for process row $i=15$ ($\rho_{16}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.084 + 0.0306 + 0.006 + 0 + 0.0474 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | -1 | 0.084 | -0.084 |
| 4 | `X_S` | 0.9 | 0.034 | 0.0306 |
| 5 | `X_I` | 0.1 | 0.06 | 0.006 |
| 8 | `S_IC` | 0.0378 | 0 | 0 |
| 10 | `S_NH` | 0.0474 | 1 | 0.0474 |
| 14 | `S_PO4` | 0.0105 | 0 | 0 |

**Sum:** $B_{15,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{15,3}| = 0$ ≤ 0.1

<a id="rho16-P"></a>
#### rho16 — P (cell 95/114)

Balance for process row $i=15$ ($\rho_{16}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.016 + 0.0045 + 0.001 + 0 + 0 + 0.0105$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | -1 | 0.016 | -0.016 |
| 4 | `X_S` | 0.9 | 0.005 | 0.0045 |
| 5 | `X_I` | 0.1 | 0.01 | 0.001 |
| 8 | `S_IC` | 0.0378 | 0 | 0 |
| 10 | `S_NH` | 0.0474 | 0 | 0 |
| 14 | `S_PO4` | 0.0105 | 1 | 0.0105 |

**Sum:** $B_{15,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{15,4}| = 0$ ≤ 0.1

<a id="rho16-H"></a>
#### rho16 — H (cell 96/114)

Balance for process row $i=15$ ($\rho_{16}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.043 + 0 + 0 + 0 + 0.010428 + 0.00105$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 1 | `X_AOB` | -1 | 0.043 | -0.043 |
| 4 | `X_S` | 0.9 | 0 | 0 |
| 5 | `X_I` | 0.1 | 0 | 0 |
| 8 | `S_IC` | 0.0378 | 0 | 0 |
| 10 | `S_NH` | 0.0474 | 0.22 | 0.010428 |
| 14 | `S_PO4` | 0.0105 | 0.1 | 0.00105 |

**Sum:** $B_{15,5} \approx -0.031522$ (check: `numpy.sum` = -0.031522)

**Status:** **OK** — $|B_{15,5}| = 0.031522$ ≤ 0.1

### rho17

*Full label:* rho17: aerobic growth X_NOB

<a id="rho17-COD"></a>
#### rho17 — COD (cell 97/114)

Balance for process row $i=16$ ($\rho_{17}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 1 + 0 + 0 + 68.6 + -91.4 + 0 + 21.8571428571$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | 1 | 1 | 1 |
| 8 | `S_IC` | -0.36 | 0 | 0 |
| 10 | `S_NH` | -0.084 | 0 | 0 |
| 11 | `S_NO2` | -20 | -3.43 | 68.6 |
| 12 | `S_NO3` | 20 | -4.57 | -91.4 |
| 14 | `S_PO4` | -0.016 | 0 | 0 |
| 15 | `S_O2` | -21.8571428571 | -1 | 21.8571428571 |

**Sum:** $B_{16,0} \approx 0.0571428571$ (check: `numpy.sum` = 0.0571428571)

**Status:** **OK** — $|B_{16,0}| = 0.0571428571$ ≤ 0.1

<a id="rho17-O"></a>
#### rho17 — O (cell 98/114)

Balance for process row $i=16$ ($\rho_{17}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.184 + -0.9612 + 0 + -45.6 + 68.6 + -0.03312 + -21.8571428571$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | 1 | 0.184 | 0.184 |
| 8 | `S_IC` | -0.36 | 2.67 | -0.9612 |
| 10 | `S_NH` | -0.084 | 0 | 0 |
| 11 | `S_NO2` | -20 | 2.28 | -45.6 |
| 12 | `S_NO3` | 20 | 3.43 | 68.6 |
| 14 | `S_PO4` | -0.016 | 2.07 | -0.03312 |
| 15 | `S_O2` | -21.8571428571 | 1 | -21.8571428571 |

**Sum:** $B_{16,1} \approx 0.3325371429$ (check: `numpy.sum` = 0.3325371429)

**Status:** **FAIL** — $|B_{16,1}| = 0.3325371429$ > 0.1

<a id="rho17-C"></a>
#### rho17 — C (cell 99/114)

Balance for process row $i=16$ ($\rho_{17}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.36 + -0.36 + 0 + 0 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | 1 | 0.36 | 0.36 |
| 8 | `S_IC` | -0.36 | 1 | -0.36 |
| 10 | `S_NH` | -0.084 | 0 | 0 |
| 11 | `S_NO2` | -20 | 0 | 0 |
| 12 | `S_NO3` | 20 | 0 | 0 |
| 14 | `S_PO4` | -0.016 | 0 | 0 |
| 15 | `S_O2` | -21.8571428571 | 0 | 0 |

**Sum:** $B_{16,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{16,2}| = 0$ ≤ 0.1

<a id="rho17-N"></a>
#### rho17 — N (cell 100/114)

Balance for process row $i=16$ ($\rho_{17}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.084 + 0 + -0.084 + -20 + 20 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | 1 | 0.084 | 0.084 |
| 8 | `S_IC` | -0.36 | 0 | 0 |
| 10 | `S_NH` | -0.084 | 1 | -0.084 |
| 11 | `S_NO2` | -20 | 1 | -20 |
| 12 | `S_NO3` | 20 | 1 | 20 |
| 14 | `S_PO4` | -0.016 | 0 | 0 |
| 15 | `S_O2` | -21.8571428571 | 0 | 0 |

**Sum:** $B_{16,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{16,3}| = 0$ ≤ 0.1

<a id="rho17-P"></a>
#### rho17 — P (cell 101/114)

Balance for process row $i=16$ ($\rho_{17}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.016 + 0 + 0 + 0 + 0 + -0.016 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | 1 | 0.016 | 0.016 |
| 8 | `S_IC` | -0.36 | 0 | 0 |
| 10 | `S_NH` | -0.084 | 0 | 0 |
| 11 | `S_NO2` | -20 | 0 | 0 |
| 12 | `S_NO3` | 20 | 0 | 0 |
| 14 | `S_PO4` | -0.016 | 1 | -0.016 |
| 15 | `S_O2` | -21.8571428571 | 0 | 0 |

**Sum:** $B_{16,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{16,4}| = 0$ ≤ 0.1

<a id="rho17-H"></a>
#### rho17 — H (cell 102/114)

Balance for process row $i=16$ ($\rho_{17}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = 0.043 + 0 + -0.01848 + -1.4 + 1.4 + -0.0016 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | 1 | 0.043 | 0.043 |
| 8 | `S_IC` | -0.36 | 0 | 0 |
| 10 | `S_NH` | -0.084 | 0.22 | -0.01848 |
| 11 | `S_NO2` | -20 | 0.07 | -1.4 |
| 12 | `S_NO3` | 20 | 0.07 | 1.4 |
| 14 | `S_PO4` | -0.016 | 0.1 | -0.0016 |
| 15 | `S_O2` | -21.8571428571 | 0 | 0 |

**Sum:** $B_{16,5} \approx 0.02292$ (check: `numpy.sum` = 0.02292)

**Status:** **OK** — $|B_{16,5}| = 0.02292$ ≤ 0.1

### rho18

*Full label:* rho18: aerobic respiration X_NOB

<a id="rho18-COD"></a>
#### rho18 — COD (cell 103/114)

Balance for process row $i=17$ ($\rho_{18}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -1 + 0 + 0 + 0 + 1$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | -1 | 1 | -1 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 0 | 0 |
| 15 | `S_O2` | -1 | -1 | 1 |

**Sum:** $B_{17,0} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{17,0}| = 0$ ≤ 0.1

<a id="rho18-O"></a>
#### rho18 — O (cell 104/114)

Balance for process row $i=17$ ($\rho_{18}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.184 + 0.9612 + 0 + 0.03312 + -1$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | -1 | 0.184 | -0.184 |
| 8 | `S_IC` | 0.36 | 2.67 | 0.9612 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 2.07 | 0.03312 |
| 15 | `S_O2` | -1 | 1 | -1 |

**Sum:** $B_{17,1} \approx -0.18968$ (check: `numpy.sum` = -0.18968)

**Status:** **FAIL** — $|B_{17,1}| = 0.18968$ > 0.1

<a id="rho18-C"></a>
#### rho18 — C (cell 105/114)

Balance for process row $i=17$ ($\rho_{18}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.36 + 0.36 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | -1 | 0.36 | -0.36 |
| 8 | `S_IC` | 0.36 | 1 | 0.36 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 0 | 0 |
| 15 | `S_O2` | -1 | 0 | 0 |

**Sum:** $B_{17,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{17,2}| = 0$ ≤ 0.1

<a id="rho18-N"></a>
#### rho18 — N (cell 106/114)

Balance for process row $i=17$ ($\rho_{18}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.084 + 0 + 0.084 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | -1 | 0.084 | -0.084 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 1 | 0.084 |
| 14 | `S_PO4` | 0.016 | 0 | 0 |
| 15 | `S_O2` | -1 | 0 | 0 |

**Sum:** $B_{17,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{17,3}| = 0$ ≤ 0.1

<a id="rho18-P"></a>
#### rho18 — P (cell 107/114)

Balance for process row $i=17$ ($\rho_{18}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.016 + 0 + 0 + 0.016 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | -1 | 0.016 | -0.016 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 0 | 0 |
| 14 | `S_PO4` | 0.016 | 1 | 0.016 |
| 15 | `S_O2` | -1 | 0 | 0 |

**Sum:** $B_{17,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{17,4}| = 0$ ≤ 0.1

<a id="rho18-H"></a>
#### rho18 — H (cell 108/114)

Balance for process row $i=17$ ($\rho_{18}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.043 + 0 + 0.01848 + 0.0016 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | -1 | 0.043 | -0.043 |
| 8 | `S_IC` | 0.36 | 0 | 0 |
| 10 | `S_NH` | 0.084 | 0.22 | 0.01848 |
| 14 | `S_PO4` | 0.016 | 0.1 | 0.0016 |
| 15 | `S_O2` | -1 | 0 | 0 |

**Sum:** $B_{17,5} \approx -0.02292$ (check: `numpy.sum` = -0.02292)

**Status:** **OK** — $|B_{17,5}| = 0.02292$ ≤ 0.1

### rho19

*Full label:* rho19: decay X_NOB

<a id="rho19-COD"></a>
#### rho19 — COD (cell 109/114)

Balance for process row $i=18$ ($\rho_{19}$) and composition row $k=0$ (**COD**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -1 + 0.9 + 0.1 + 0 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | -1 | 1 | -1 |
| 4 | `X_S` | 0.9 | 1 | 0.9 |
| 5 | `X_I` | 0.1 | 1 | 0.1 |
| 8 | `S_IC` | 0.0378 | 0 | 0 |
| 10 | `S_NH` | 0.0474 | 0 | 0 |
| 14 | `S_PO4` | 0.0105 | 0 | 0 |

**Sum:** $B_{18,0} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{18,0}| = 0$ ≤ 0.1

<a id="rho19-O"></a>
#### rho19 — O (cell 110/114)

Balance for process row $i=18$ ($\rho_{19}$) and composition row $k=1$ (**O**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.184 + 0.1404 + 0.015 + 0.100926 + 0 + 0.021735$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | -1 | 0.184 | -0.184 |
| 4 | `X_S` | 0.9 | 0.156 | 0.1404 |
| 5 | `X_I` | 0.1 | 0.15 | 0.015 |
| 8 | `S_IC` | 0.0378 | 2.67 | 0.100926 |
| 10 | `S_NH` | 0.0474 | 0 | 0 |
| 14 | `S_PO4` | 0.0105 | 2.07 | 0.021735 |

**Sum:** $B_{18,1} \approx 0.094061$ (check: `numpy.sum` = 0.094061)

**Status:** **OK** — $|B_{18,1}| = 0.094061$ ≤ 0.1

<a id="rho19-C"></a>
#### rho19 — C (cell 111/114)

Balance for process row $i=18$ ($\rho_{19}$) and composition row $k=2$ (**C**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.36 + 0.2862 + 0.036 + 0.0378 + 0 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | -1 | 0.36 | -0.36 |
| 4 | `X_S` | 0.9 | 0.318 | 0.2862 |
| 5 | `X_I` | 0.1 | 0.36 | 0.036 |
| 8 | `S_IC` | 0.0378 | 1 | 0.0378 |
| 10 | `S_NH` | 0.0474 | 0 | 0 |
| 14 | `S_PO4` | 0.0105 | 0 | 0 |

**Sum:** $B_{18,2} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{18,2}| = 0$ ≤ 0.1

<a id="rho19-N"></a>
#### rho19 — N (cell 112/114)

Balance for process row $i=18$ ($\rho_{19}$) and composition row $k=3$ (**N**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.084 + 0.0306 + 0.006 + 0 + 0.0474 + 0$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | -1 | 0.084 | -0.084 |
| 4 | `X_S` | 0.9 | 0.034 | 0.0306 |
| 5 | `X_I` | 0.1 | 0.06 | 0.006 |
| 8 | `S_IC` | 0.0378 | 0 | 0 |
| 10 | `S_NH` | 0.0474 | 1 | 0.0474 |
| 14 | `S_PO4` | 0.0105 | 0 | 0 |

**Sum:** $B_{18,3} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{18,3}| = 0$ ≤ 0.1

<a id="rho19-P"></a>
#### rho19 — P (cell 113/114)

Balance for process row $i=18$ ($\rho_{19}$) and composition row $k=4$ (**P**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.016 + 0.0045 + 0.001 + 0 + 0 + 0.0105$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | -1 | 0.016 | -0.016 |
| 4 | `X_S` | 0.9 | 0.005 | 0.0045 |
| 5 | `X_I` | 0.1 | 0.01 | 0.001 |
| 8 | `S_IC` | 0.0378 | 0 | 0 |
| 10 | `S_NH` | 0.0474 | 0 | 0 |
| 14 | `S_PO4` | 0.0105 | 1 | 0.0105 |

**Sum:** $B_{18,4} \approx 0$ (check: `numpy.sum` = 0)

**Status:** **OK** — $|B_{18,4}| = 0$ ≤ 0.1

<a id="rho19-H"></a>
#### rho19 — H (cell 114/114)

Balance for process row $i=18$ ($\rho_{19}$) and composition row $k=5$ (**H**):

$$B_{i,k} = \sum_{j} S_{i,j}\, I_{k,j}.$$

Non-zero contributions ($S_{i,j}\neq 0$), in column order $j = 0\ldots 16$:

$$B_{i,k} = -0.043 + 0 + 0 + 0 + 0.010428 + 0.00105$$

| $j$ | Variable | $S_{i,j}$ | $I_{k,j}$ | $S_{i,j} I_{k,j}$ |
|:---:|:---|:---:|:---:|:---:|
| 2 | `X_NOB` | -1 | 0.043 | -0.043 |
| 4 | `X_S` | 0.9 | 0 | 0 |
| 5 | `X_I` | 0.1 | 0 | 0 |
| 8 | `S_IC` | 0.0378 | 0 | 0 |
| 10 | `S_NH` | 0.0474 | 0.22 | 0.010428 |
| 14 | `S_PO4` | 0.0105 | 0.1 | 0.00105 |

**Sum:** $B_{18,5} \approx -0.031522$ (check: `numpy.sum` = -0.031522)

**Status:** **OK** — $|B_{18,5}| = 0.031522$ ≤ 0.1

