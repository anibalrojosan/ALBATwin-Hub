# ADR 003: DAE Resolution Strategy (Nested Solver)

**Context:**
The model is a Differential-Algebraic Equation (DAE) system:
*   Differential: $dX/dt = \dots$ (Biology)
*   Algebraic: $ChargeBalance(H^+) = 0$ (pH)

Solving DAEs directly in Python (`scikits.odes`) is complex and dependency-heavy.

**Decision:**
Use a **Nested Algebraic Solver (Operator Splitting)** approach.
*   At each time step of the ODE solver, pause and solve the algebraic pH system using **Newton-Raphson** (`scipy.optimize.root`).
*   The pH calculation is treated as a "sub-routine" within the derivative function.

**Consequences:**
*   (+) Robustness: Decouples the fast chemical equilibrium from the slow biological dynamics.
*   (+) Compatibility: Allows using standard ODE solvers (SciPy) instead of specialized DAE solvers.
*   (-) Cost: Adds computational overhead to every function evaluation.
