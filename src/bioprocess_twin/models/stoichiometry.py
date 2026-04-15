"""
ALBA model stoichiometry: Petersen matrix and composition matrix.

Implements the 19x17 Petersen matrix (biological processes) and 6x17 composition
matrix for mass balance validation. Based on Casagli et al. (2021) and ASM conventions.

**Supplementary material (SI) mapping**

- **Table SI.3.1** (`docs/papers/supplementary_material_ALBA/...SI.3.1...md`): process
  layout and composition matrix structure. Numeric composition cells here must match
  `get_composition_matrix()` (the repo table is kept aligned with this module).
- **Table SI.3.2**: all stoichiometric *parameters* (yields, i_C, i_N, f_XI, …) used in
  both the composition matrix and Petersen coefficients. Kinetic parameters from
  other SI tables do **not** appear in **S**; they only scale **ρ** in **dC/dt = Sᵀρ**.
- **Table SI.3.3**: closed-form expressions for Petersen entries **α**; implemented in
  `get_petersen_matrix()`.
"""

import numpy as np

# ---------------------------------------------------------------------------
# State variable column indices (must match StateVector.to_array order for j=0..16;
# extended layout appends S_H_PROTON at j=17 — see core.state.StateVectorVariant)
# ---------------------------------------------------------------------------
X_ALG, X_AOB, X_NOB, X_H, X_S, X_I = 0, 1, 2, 3, 4, 5
S_S, S_I, S_IC, S_ND, S_NH, S_NO2, S_NO3, S_N2, S_PO4, S_O2, S_H2O = (
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
)

N_STATE = 17
N_PROCESSES = 19
N_ELEMENTS = 6  # COD, O, C, N, P, H

# ---------------------------------------------------------------------------
# Stoichiometric constants (MATH_MODEL.md Section 1.1, SI.3.2)
# ---------------------------------------------------------------------------
# Algae biomass composition
I_C_ALG = 0.327
I_N_ALG = 0.042
I_P_ALG = 0.008
I_O_ALG = 0.209
I_H_ALG = 0.050
F_XI_ALG = 0.062  # Inert fraction from algae decay

# Bacterial biomass composition (AOB, NOB, H)
I_C_BM = 0.36
I_N_BM = 0.084
I_P_BM = 0.016
I_O_BM = 0.184
I_H_BM = 0.043
F_XI = 0.1  # Inert fraction from bacteria decay
F_SI = 0.1  # Soluble inert from hydrolysis

# Particulate organics (X_S, X_I)
I_C_XS, I_N_XS, I_P_XS, I_O_XS = 0.318, 0.034, 0.005, 0.156
I_C_XI, I_N_XI, I_P_XI, I_O_XI = 0.36, 0.06, 0.01, 0.15

# Soluble organics (S_S, S_I)
I_C_SS, I_N_SS, I_P_SS, I_O_SS = 0.318, 0.015, 0.005, 0.156
I_C_SI, I_N_SI, I_P_SI, I_O_SI = 0.36, 0.06, 0.005, 0.15

# Urea (S_ND)
I_C_ND = 0.429
I_H_ND = 0.072

# Yields
Y_H = 0.63
Y_H_NO2 = 0.3
Y_H_NO3 = 0.5
Y_AOB = 0.2
Y_NOB = 0.05

# Table SI.3.3: fixed S_H2O coefficients for algal processes (g H / g COD_BM)
ALPHA_ALG_H2O_RHO1 = -0.0404
ALPHA_ALG_H2O_RHO2 = -0.0464
ALPHA_ALG_H2O_RHO3 = 0.0404

# Composition matrix fixed values (SI.3.1; cross-check against Casagli et al. 2021 SI)
# S_IC: O=2.67, C=1, H=0 (COD row for S_IC is 0; see docs/stoichiometry-external-comparison.md)
# S_ND: O=0.57, C=0.43 (I_C_S_ND), H=0.14
# S_NH: H=0.22
# S_NO2, S_NO3: O=2.28, 3.43; H=0.07
# S_PO4: O=2.07, H=0.10
# S_O2: O=1
# S_H2O: O=7.94, H=1
I_O_S_IC = 2.67
I_H_S_IC = 0
I_O_S_ND = 0.57
I_C_S_ND = 0.43
I_H_S_ND = 0.14
I_H_S_NH = 0.22
I_O_S_NO2 = 2.28
I_O_S_NO3 = 3.43
I_H_S_NO2 = 0.07
I_H_S_NO3 = 0.07
I_O_S_PO4 = 2.07
I_H_S_PO4 = 0.10
I_O_S_O2 = 1.0
I_O_S_H2O = 7.94
I_H_S_H2O = 1.0


def get_composition_matrix() -> np.ndarray:
    """
    Return the 6x17 composition matrix.

    Rows: COD (0), O (1), C (2), N (3), P (4), H (5).
    Columns: 17 state variables in StateVector order.
    """
    comp = np.zeros((N_ELEMENTS, N_STATE))

    # Row 0: COD
    comp[0, :] = [
        1,
        1,
        1,
        1,
        1,
        1,  # X_ALG, X_AOB, X_NOB, X_H, X_S, X_I
        1,
        1,
        0,  # S_S, S_I, S_IC
        0,
        0,  # S_ND, S_NH (no COD)
        -3.43,
        -4.57,
        -1.71,  # S_NO2, S_NO3, S_N2 (electron acceptors)
        0,
        -1,
        0,  # S_PO4, S_O2, S_H2O
    ]

    # Row 1: O
    comp[1, :] = [
        I_O_ALG,
        I_O_BM,
        I_O_BM,
        I_O_BM,
        I_O_XS,
        I_O_XI,
        I_O_SS,
        I_O_SI,
        I_O_S_IC,
        I_O_S_ND,
        0,
        I_O_S_NO2,
        I_O_S_NO3,
        0,
        I_O_S_PO4,
        I_O_S_O2,
        I_O_S_H2O,
    ]

    # Row 2: C
    comp[2, :] = [
        I_C_ALG,
        I_C_BM,
        I_C_BM,
        I_C_BM,
        I_C_XS,
        I_C_XI,
        I_C_SS,
        I_C_SI,
        1,
        I_C_S_ND,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ]

    # Row 3: N
    comp[3, :] = [
        I_N_ALG,
        I_N_BM,
        I_N_BM,
        I_N_BM,
        I_N_XS,
        I_N_XI,
        I_N_SS,
        I_N_SI,
        0,
        1,
        1,
        1,
        1,
        1,
        0,
        0,
        0,
    ]

    # Row 4: P
    comp[4, :] = [
        I_P_ALG,
        I_P_BM,
        I_P_BM,
        I_P_BM,
        I_P_XS,
        I_P_XI,
        I_P_SS,
        I_P_SI,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
    ]

    # Row 5: H
    comp[5, :] = [
        I_H_ALG,
        I_H_BM,
        I_H_BM,
        I_H_BM,
        0,
        0,
        0,
        0,
        I_H_S_IC,
        I_H_S_ND,
        I_H_S_NH,
        I_H_S_NO2,
        I_H_S_NO3,
        0,
        I_H_S_PO4,
        0,
        I_H_S_H2O,
    ]

    return comp


def _alpha_alg_o2_growth_nh4() -> float:
    """SI.3.3 rho1: S_O2 coefficient for phototrophic growth on NH4+."""
    return (
        -I_O_ALG
        + (32.0 / 12.0) * I_C_ALG
        - (24.0 / 14.0) * I_N_ALG
        + (40.0 / 31.0) * I_P_ALG
        + 8.0 * I_H_ALG
    )


def _alpha_alg_o2_growth_no3() -> float:
    """SI.3.3 rho2: S_O2 coefficient for phototrophic growth on NO3-."""
    return (
        -I_O_ALG
        + (32.0 / 12.0) * I_C_ALG
        + (40.0 / 14.0) * I_N_ALG
        + (40.0 / 31.0) * I_P_ALG
        + 8.0 * I_H_ALG
    )


def _alpha_alg_o2_respiration() -> float:
    """SI.3.3 rho3: S_O2 coefficient for algal aerobic respiration."""
    return (
        + I_O_ALG
        - (32.0 / 12.0) * I_C_ALG
        + (24.0 / 14.0) * I_N_ALG
        - (40.0 / 31.0) * I_P_ALG
        - 8.0 * I_H_ALG
    )


def get_petersen_matrix() -> np.ndarray:
    """
    Return the 19x17 Petersen matrix (biological processes only).

    Stoichiometric coefficients follow Table SI.3.3 in the ALBA supplementary
    material (docs/papers/supplementary_material_ALBA/). Processes rho20-22
    (equilibrium) are excluded (hydrochemistry sprint).
    """
    S = np.zeros((N_PROCESSES, N_STATE))

    # --- rho1: Phototrophic growth on NH4+ ---
    S[0, X_ALG] = 1
    S[0, S_IC] = -I_C_ALG
    S[0, S_NH] = -I_N_ALG
    S[0, S_PO4] = -I_P_ALG
    S[0, S_O2] = _alpha_alg_o2_growth_nh4()
    S[0, S_H2O] = ALPHA_ALG_H2O_RHO1

    # --- rho2: Phototrophic growth on NO3- ---
    S[1, X_ALG] = 1
    S[1, S_IC] = -I_C_ALG
    S[1, S_NO3] = -I_N_ALG
    S[1, S_PO4] = -I_P_ALG
    S[1, S_O2] = _alpha_alg_o2_growth_no3()
    S[1, S_H2O] = ALPHA_ALG_H2O_RHO2

    # --- rho3: Aerobic respiration of X_ALG ---
    S[2, X_ALG] = -1
    S[2, S_IC] = I_C_ALG
    S[2, S_NH] = I_N_ALG
    S[2, S_PO4] = I_P_ALG
    S[2, S_O2] = _alpha_alg_o2_respiration()
    S[2, S_H2O] = ALPHA_ALG_H2O_RHO3

    # --- rho4: Decay of X_ALG ---
    S[3, X_ALG] = -1
    S[3, X_S] = 1 - F_XI_ALG
    S[3, X_I] = F_XI_ALG
    S[3, S_IC] = I_C_ALG - (1 - F_XI_ALG) * I_C_XS - F_XI_ALG * I_C_XI
    S[3, S_NH] = I_N_ALG - (1 - F_XI_ALG) * I_N_XS - F_XI_ALG * I_N_XI
    S[3, S_PO4] = I_P_ALG - (1 - F_XI_ALG) * I_P_XS - F_XI_ALG * I_P_XI

    # --- rho5: Aerobic growth of X_H on NH4+ ---
    S[4, X_H] = 1
    S[4, S_S] = -1 / Y_H
    S[4, S_IC] = I_C_SS / Y_H - I_C_BM
    S[4, S_NH] = I_N_SS / Y_H - I_N_BM
    S[4, S_PO4] = I_P_SS / Y_H - I_P_BM
    S[4, S_O2] = -(1 / Y_H - 1)

    # --- rho6: Aerobic growth of X_H on NO3- ---
    n_no3_r6 = I_N_SS / Y_H - I_N_BM
    S[5, X_H] = 1
    S[5, S_S] = -1 / Y_H
    S[5, S_IC] = I_C_SS / Y_H - I_C_BM
    S[5, S_NO3] = n_no3_r6
    S[5, S_PO4] = I_P_SS / Y_H - I_P_BM
    S[5, S_O2] = -(1 / Y_H - 1) - (64.0 / 14.0) * n_no3_r6

    # --- rho7: Aerobic respiration of X_H ---
    S[6, X_H] = -1
    S[6, S_IC] = I_C_BM
    S[6, S_NH] = I_N_BM
    S[6, S_PO4] = I_P_BM
    S[6, S_O2] = -1

    # --- rho8: Anoxic growth of X_H on NO3- (SI: 28/80 factor) ---
    yld = Y_H_NO3
    oex = 1 / yld - 1
    S[7, X_H] = 1
    S[7, S_S] = -1 / yld
    S[7, S_IC] = I_C_SS / yld - I_C_BM
    S[7, S_NH] = I_N_SS / yld - I_N_BM
    S[7, S_NO3] = -(28.0 / 80.0) * oex
    S[7, S_N2] = (28.0 / 80.0) * oex
    S[7, S_PO4] = I_P_SS / yld - I_P_BM

    # --- rho9: Anoxic growth of X_H on NO2- (SI: 28/48 factor) ---
    yld2 = Y_H_NO2
    oex2 = 1 / yld2 - 1
    S[8, X_H] = 1
    S[8, S_S] = -1 / yld2
    S[8, S_IC] = I_C_SS / yld2 - I_C_BM
    S[8, S_NH] = I_N_SS / yld2 - I_N_BM
    S[8, S_NO2] = -(28.0 / 48.0) * oex2
    S[8, S_N2] = (28.0 / 48.0) * oex2
    S[8, S_PO4] = I_P_SS / yld2 - I_P_BM

    # --- rho10: Anoxic respiration of X_H on NO2- and NO3- ---
    S[9, X_H] = -1
    S[9, S_IC] = I_C_BM
    S[9, S_NH] = I_N_BM
    S[9, S_NO2] = -14.0 / 64.0
    S[9, S_NO3] = -14.0 / 64.0
    S[9, S_N2] = 28.0 / 64.0
    S[9, S_PO4] = I_P_BM

    # --- rho11: Hydrolysis of slowly biodegradable COD ---
    S[10, X_S] = -1
    S[10, S_S] = 1 - F_SI
    S[10, S_I] = F_SI
    S[10, S_IC] = I_C_XS - (1 - F_SI) * I_C_SS - F_SI * I_C_SI
    S[10, S_NH] = I_N_XS - (1 - F_SI) * I_N_SS - F_SI * I_N_SI
    S[10, S_PO4] = I_P_XS - (1 - F_SI) * I_P_SS - F_SI * I_P_SI

    # --- rho12: Hydrolysis of urea ---
    S[11, S_IC] = I_C_ND
    S[11, S_ND] = -1
    S[11, S_NH] = 1
    S[11, S_H2O] = I_H_ND

    # --- rho13: Decay of X_H ---
    S[12, X_H] = -1
    S[12, X_S] = 1 - F_XI
    S[12, X_I] = F_XI
    S[12, S_IC] = I_C_BM - (1 - F_XI) * I_C_XS - F_XI * I_C_XI
    S[12, S_NH] = I_N_BM - (1 - F_XI) * I_N_XS - F_XI * I_N_XI
    S[12, S_PO4] = I_P_BM - (1 - F_XI) * I_P_XS - F_XI * I_P_XI

    # --- rho14: Aerobic growth of X_AOB on NH4+ ---
    S[13, X_AOB] = 1
    S[13, S_IC] = -I_C_BM
    S[13, S_NH] = -I_N_BM - (1 / Y_AOB) 
    S[13, S_NO2] = 1 / Y_AOB
    S[13, S_PO4] = -I_P_BM
    S[13, S_O2] = 1 - (48.0 / 14.0) / Y_AOB

    # --- rho15: Aerobic respiration of X_AOB ---
    S[14, X_AOB] = -1
    S[14, S_IC] = I_C_BM
    S[14, S_NH] = I_N_BM
    S[14, S_PO4] = I_P_BM
    S[14, S_O2] = -1

    # --- rho16: Decay of X_AOB ---
    S[15, X_AOB] = -1
    S[15, X_S] = 1 - F_XI
    S[15, X_I] = F_XI
    S[15, S_IC] = I_C_BM - (1 - F_XI) * I_C_XS - F_XI * I_C_XI
    S[15, S_NH] = I_N_BM - (1 - F_XI) * I_N_XS - F_XI * I_N_XI
    S[15, S_PO4] = I_P_BM - (1 - F_XI) * I_P_XS - F_XI * I_P_XI

    # --- rho17: Aerobic growth of X_NOB on NO2- ---
    S[16, X_NOB] = 1
    S[16, S_IC] = -I_C_BM
    S[16, S_NH] = -I_N_BM
    S[16, S_NO2] = -1 / Y_NOB
    S[16, S_NO3] = 1 / Y_NOB
    S[16, S_PO4] = -I_P_BM
    S[16, S_O2] = 1 - (16.0 / 14.0) / Y_NOB

    # --- rho18: Aerobic respiration of X_NOB ---
    S[17, X_NOB] = -1
    S[17, S_IC] = I_C_BM
    S[17, S_NH] = I_N_BM
    S[17, S_PO4] = I_P_BM
    S[17, S_O2] = -1

    # --- rho19: Decay of X_NOB ---
    S[18, X_NOB] = -1
    S[18, X_S] = 1 - F_XI
    S[18, X_I] = F_XI
    S[18, S_IC] = I_C_BM - (1 - F_XI) * I_C_XS - F_XI * I_C_XI
    S[18, S_NH] = I_N_BM - (1 - F_XI) * I_N_XS - F_XI * I_N_XI
    S[18, S_PO4] = I_P_BM - (1 - F_XI) * I_P_XS - F_XI * I_P_XI

    return S
