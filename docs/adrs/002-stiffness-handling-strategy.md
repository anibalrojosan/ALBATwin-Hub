# ADR 002: Stiffness Handling Strategy

**Context:**
Bioprocesses are inherently "stiff" systems due to disparate time scales:
*   Chemical equilibria (pH): Milliseconds ($10^{-3}$ s).
*   Gas transfer ($k_La$): Minutes ($10^2$ s).
*   Biomass growth ($\mu$): Days ($10^5$ s).

Explicit solvers (like Runge-Kutta 45) fail or become extremely slow as they must take tiny steps to maintain stability.

**Decision:**
Use **Implicit Solvers** specifically designed for stiff systems.
*   **Primary:** `LSODA` (Livermore Solver for Ordinary Differential Equations) - automatically switches between non-stiff (Adams) and stiff (BDF) methods.
*   **Secondary:** `Radau` (Implicit Runge-Kutta) or `BDF` (Backward Differentiation Formula).

**Consequences:**
*   (+) Stability: Can take large time steps when the system is slow, without blowing up.
*   (+) Speed: Orders of magnitude faster for long simulations (e.g., 365 days).
*   (-) Complexity: Requires computing the Jacobian matrix (can be approximated internally by the solver).
