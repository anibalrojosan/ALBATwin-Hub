"""
ALBA model stoichiometry: Petersen matrix and composition matrix.

Implements the 19x17 Petersen matrix (biological processes) and 6x17 composition
matrix for mass balance validation. Based on Casagli et al. (2021) and ASM conventions.
"""

import numpy as np

# ---------------------------------------------------------------------------
# State variable column indices (must match StateVector.to_array order)
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

# Composition matrix fixed values (SI.3.1)
# S_IC: O=2.67, C=1, H=0.14
# S_ND: O=0.57, C=0.43 (i_C_ND), H=0.22 (SI.3.1 table)
# S_NH: H=0.07
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


def _solve_alpha_h2o_from_h_balance(row: np.ndarray, comp: np.ndarray) -> float:
    """Solve for S_H2O coefficient from H balance: row @ comp[5,:] = 0."""
    contrib = np.dot(row, comp[5, :])
    return -contrib / comp[5, S_H2O]


def _solve_alpha_from_balance(
    row: np.ndarray, comp: np.ndarray, unknown_idx: int, element_row: int
) -> float:
    """Solve for single unknown from one balance equation."""
    contrib = np.dot(row, comp[element_row, :])
    return -contrib / comp[element_row, unknown_idx]


def get_petersen_matrix() -> np.ndarray:
    """
    Return the 19x17 Petersen matrix (biological processes only).

    Processes 20-22 (equilibrium) are excluded (hydrochemistry sprint).
    """
    comp = get_composition_matrix()
    S = np.zeros((N_PROCESSES, N_STATE))

    # --- Process 1: Phototrophic growth on NH4 ---
    S[0, X_ALG] = 1
    S[0, S_IC] = -I_C_ALG
    S[0, S_NH] = -I_N_ALG
    S[0, S_PO4] = -I_P_ALG
    S[0, S_O2] = 1
    S[0, S_H2O] = _solve_alpha_h2o_from_h_balance(S[0], comp)

    # --- Process 2: Phototrophic growth on NO3 ---
    S[1, X_ALG] = 1
    S[1, S_IC] = -I_C_ALG
    S[1, S_NO3] = -I_N_ALG
    S[1, S_PO4] = -I_P_ALG
    S[1, S_H2O] = _solve_alpha_h2o_from_h_balance(S[1], comp)

    # --- Process 3: Algae aerobic respiration ---
    S[2, X_ALG] = -1
    S[2, S_O2] = -1
    S[2, S_IC] = _solve_alpha_from_balance(S[2], comp, S_IC, 2)
    S[2, S_NH] = _solve_alpha_from_balance(S[2], comp, S_NH, 3)
    S[2, S_PO4] = _solve_alpha_from_balance(S[2], comp, S_PO4, 4)
    S[2, S_H2O] = _solve_alpha_h2o_from_h_balance(S[2], comp)

    # --- Process 4: Algae decay ---
    # -1 X_ALG -> (1-f_Xi_ALG) X_S, f_Xi_ALG X_I; N, P, C released
    S[3, X_ALG] = -1
    S[3, X_S] = 1 - F_XI_ALG
    S[3, X_I] = F_XI_ALG
    S[3, S_IC] = _solve_alpha_from_balance(S[3], comp, S_IC, 2)
    S[3, S_NH] = _solve_alpha_from_balance(S[3], comp, S_NH, 3)
    S[3, S_PO4] = _solve_alpha_from_balance(S[3], comp, S_PO4, 4)

    # --- Process 5: Heterotrophic aerobic growth on NH4 ---
    S[4, X_H] = 1
    S[4, S_S] = -1 / Y_H
    S[4, S_IC] = -I_C_BM + I_C_SS / Y_H
    S[4, S_NH] = -I_N_BM + I_N_SS / Y_H
    S[4, S_PO4] = -I_P_BM + I_P_SS / Y_H
    S[4, S_O2] = _solve_alpha_from_balance(S[4], comp, S_O2, 0)  # COD balance

    # --- Process 6: Heterotrophic aerobic growth on NO3 ---
    S[5, X_H] = 1
    S[5, S_S] = -1 / Y_H
    S[5, S_IC] = -I_C_BM + I_C_SS / Y_H
    S[5, S_NO3] = -I_N_BM + I_N_SS / Y_H
    S[5, S_PO4] = -I_P_BM + I_P_SS / Y_H
    S[5, S_O2] = _solve_alpha_from_balance(S[5], comp, S_O2, 0)  # COD balance

    # --- Process 7: Heterotrophic aerobic respiration ---
    S[6, X_H] = -1
    S[6, S_IC] = _solve_alpha_from_balance(S[6], comp, S_IC, 2)
    S[6, S_NH] = _solve_alpha_from_balance(S[6], comp, S_NH, 3)
    S[6, S_PO4] = _solve_alpha_from_balance(S[6], comp, S_PO4, 4)
    S[6, S_O2] = _solve_alpha_from_balance(S[6], comp, S_O2, 0)  # COD balance

    # --- Process 8: Anoxic growth on NO3 ---
    # 2.86 gCOD/gN for NO3 reduction (ASM convention)
    S[7, X_H] = 1
    S[7, S_S] = -1 / Y_H_NO3
    S[7, S_IC] = -I_C_BM + I_C_SS / Y_H_NO3
    S[7, S_NH] = -I_N_BM + I_N_SS / Y_H_NO3
    S[7, S_NO3] = -(1 - Y_H_NO3) / Y_H_NO3 / 2.86
    S[7, S_N2] = (1 - Y_H_NO3) / Y_H_NO3 / 2.86
    S[7, S_PO4] = -I_P_BM + I_P_SS / Y_H_NO3

    # --- Process 9: Anoxic growth on NO2 ---
    # 1.71 gCOD/gN for NO2 reduction (ASM convention)
    S[8, X_H] = 1
    S[8, S_S] = -1 / Y_H_NO2
    S[8, S_IC] = -I_C_BM + I_C_SS / Y_H_NO2
    S[8, S_NH] = -I_N_BM + I_N_SS / Y_H_NO2
    S[8, S_NO2] = -(1 - Y_H_NO2) / Y_H_NO2 / 1.71
    S[8, S_N2] = (1 - Y_H_NO2) / Y_H_NO2 / 1.71
    S[8, S_PO4] = -I_P_BM + I_P_SS / Y_H_NO2

    # --- Process 10: Anoxic respiration (NO2 and NO3) ---
    # Endogenous respiration: N to S_NH, use NO3 as electron acceptor
    S[9, X_H] = -1
    S[9, S_NH] = I_N_BM
    S[9, S_NO3] = -1 / 2.86
    S[9, S_N2] = 0.5 / 2.86
    S[9, S_IC] = _solve_alpha_from_balance(S[9], comp, S_IC, 2)
    S[9, S_PO4] = _solve_alpha_from_balance(S[9], comp, S_PO4, 4)

    # --- Process 11: Hydrolysis ---
    S[10, X_S] = -1
    S[10, S_S] = 1 - F_SI
    S[10, S_I] = F_SI
    S[10, S_IC] = _solve_alpha_from_balance(S[10], comp, S_IC, 2)
    S[10, S_NH] = _solve_alpha_from_balance(S[10], comp, S_NH, 3)
    S[10, S_PO4] = _solve_alpha_from_balance(S[10], comp, S_PO4, 4)

    # --- Process 12: Urea hydrolysis ---
    S[11, S_ND] = -1
    S[11, S_NH] = 1
    S[11, S_IC] = I_C_ND  # Carbon from urea
    S[11, S_H2O] = _solve_alpha_h2o_from_h_balance(S[11], comp)

    # --- Process 13: Heterotrophic decay ---
    S[12, X_H] = -1
    S[12, X_S] = 1 - F_XI
    S[12, X_I] = F_XI
    S[12, S_IC] = _solve_alpha_from_balance(S[12], comp, S_IC, 2)
    S[12, S_NH] = _solve_alpha_from_balance(S[12], comp, S_NH, 3)
    S[12, S_PO4] = _solve_alpha_from_balance(S[12], comp, S_PO4, 4)

    # --- Process 14: AOB growth ---
    S[13, X_AOB] = 1
    S[13, S_IC] = -I_C_BM
    S[13, S_NH] = -1 / Y_AOB
    S[13, S_NO2] = 1 / Y_AOB
    S[13, S_PO4] = -I_P_BM
    S[13, S_O2] = _solve_alpha_from_balance(S[13], comp, S_O2, 0)  # COD balance

    # --- Process 15: AOB respiration ---
    S[14, X_AOB] = -1
    S[14, S_IC] = _solve_alpha_from_balance(S[14], comp, S_IC, 2)
    S[14, S_NH] = _solve_alpha_from_balance(S[14], comp, S_NH, 3)
    S[14, S_PO4] = _solve_alpha_from_balance(S[14], comp, S_PO4, 4)
    S[14, S_O2] = _solve_alpha_from_balance(S[14], comp, S_O2, 0)

    # --- Process 16: AOB decay ---
    S[15, X_AOB] = -1
    S[15, X_S] = 1 - F_XI
    S[15, X_I] = F_XI
    S[15, S_IC] = _solve_alpha_from_balance(S[15], comp, S_IC, 2)
    S[15, S_NH] = _solve_alpha_from_balance(S[15], comp, S_NH, 3)
    S[15, S_PO4] = _solve_alpha_from_balance(S[15], comp, S_PO4, 4)

    # --- Process 17: NOB growth ---
    S[16, X_NOB] = 1
    S[16, S_IC] = -I_C_BM
    S[16, S_NH] = 0
    S[16, S_NO2] = -1 / Y_NOB
    S[16, S_NO3] = 1 / Y_NOB
    S[16, S_PO4] = -I_P_BM
    S[16, S_O2] = _solve_alpha_from_balance(S[16], comp, S_O2, 0)  # COD balance

    # --- Process 18: NOB respiration ---
    S[17, X_NOB] = -1
    S[17, S_IC] = _solve_alpha_from_balance(S[17], comp, S_IC, 2)
    S[17, S_NH] = _solve_alpha_from_balance(S[17], comp, S_NH, 3)
    S[17, S_PO4] = _solve_alpha_from_balance(S[17], comp, S_PO4, 4)
    S[17, S_O2] = _solve_alpha_from_balance(S[17], comp, S_O2, 0)

    # --- Process 19: NOB decay ---
    S[18, X_NOB] = -1
    S[18, X_S] = 1 - F_XI
    S[18, X_I] = F_XI
    S[18, S_IC] = _solve_alpha_from_balance(S[18], comp, S_IC, 2)
    S[18, S_NH] = _solve_alpha_from_balance(S[18], comp, S_NH, 3)
    S[18, S_PO4] = _solve_alpha_from_balance(S[18], comp, S_PO4, 4)

    return S
