# Balance of oxygen for ρ₁₄: Aerobic growth of X_AOB on NH₄⁺ (AOB)


## 1. Row **O** of **`I`** (`comp[1, :]`)

Only the columns with a coefficient ≠ 0 in ρ₁₄ are relevant. For oxygen:

| Variable | `I_O` in code | Constant |
|----------|------------------|-----------|
| `X_AOB`  | `I_O_BM`         | 0.184 |
| `S_IC`   | `I_O_S_IC`       | 2.67 |
| `S_NH`   | **0**            | (no O in this row) |
| `S_NO2`  | `I_O_S_NO2`      | 2.28 |
| `S_PO4`  | `I_O_S_PO4`      | 2.07 |
| `S_O2`   | `I_O_S_O2`       | 1.0 |

Definition in code:

```142:161:src/bioprocess_twin/models/stoichiometry.py
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
```

(`S_NH` is the column with `0` between `S_ND` and `S_NO2`.)

## 2. Row **ρ₁₄** of **`S`** (index 13)

```402:408:src/bioprocess_twin/models/stoichiometry.py
    # --- rho14: Aerobic growth of X_AOB on NH4+ ---
    S[13, X_AOB] = 1
    S[13, S_IC] = -I_C_BM
    S[13, S_NH] = -I_N_BM - (1 / Y_AOB) 
    S[13, S_NO2] = 1 / Y_AOB
    S[13, S_PO4] = -I_P_BM
    S[13, S_O2] = 1 - (48.0 / 14.0) / Y_AOB
```

With `Y_AOB = 0.2` → `1/Y_AOB = 5`.

## 3. Expression of the balance (what the test evaluates)

Para el proceso 14 y el elemento O:

$$
B_{14,\mathrm{O}} = \sum_{j} S_{13,j}\, I_{\mathrm{O},j}.
$$

Solo contribuyen columnas con $S_{13,j}\neq 0$ y $I_{\mathrm{O},j}$ relevante:

$$
\begin{aligned}
B_{14,\mathrm{O}} &=
(1)\,I_{\mathrm{O,BM}} \\
&\quad + (-I_{C,\mathrm{BM}})\,I_{\mathrm{O},S_{\mathrm{IC}}} \\
&\quad + \underbrace{(\ldots)\cdot 0}_{S_{\mathrm{NH}}} \\
&\quad + \frac{1}{Y_{\mathrm{AOB}}}\,I_{\mathrm{O},S_{\mathrm{NO2}}} \\
&\quad + (-I_{P,\mathrm{BM}})\,I_{\mathrm{O},S_{\mathrm{PO4}}} \\
&\quad + \Bigl(1 - \frac{48/14}{Y_{\mathrm{AOB}}}\Bigr)\,I_{\mathrm{O},S_{\mathrm{O2}}}.
\end{aligned}
$$

Substituting numbers from the code (`I_C_BM=0.36`, `I_P_BM=0.016`, etc.):

| Term | Calculation | Approximate value |
|--------|---------|----------------|
| Biomasa | $1 \times 0.184$ | $+0.184$ |
| $S_{\mathrm{IC}}$ | $-0.36 \times 2.67$ | $-0.9612$ |
| $S_{\mathrm{NO2}}$ | $5 \times 2.28$ | $+11.4$ |
| $S_{\mathrm{PO4}}$ | $-0.016 \times 2.07$ | $-0.03312$ |
| $S_{\mathrm{O2}}$ | $\bigl(1 - \frac{48/14}{0.2}\bigr) \times 1 = 1 - \frac{120}{7}$ | $\approx -16.142857$ |

Sum:

$$
0.184 - 0.9612 + 11.4 - 0.03312 - 16.142857 \approx -5.553.
$$

That matches **`−5.553177e+00`** from the test.

---

> Further investigation is needed to understand this large residual.