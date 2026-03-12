# ADR 001: Simulation Engine Stack

**Context:**
The ALBA model involves complex biological kinetics and requires high-performance numerical integration. While MATLAB is traditional in academia, the industry standard for production software is Python. However, pure Python is too slow for stiff ODE systems.

**Decision:**
Use **Python + SciPy + Numba (JIT)** as the core simulation engine.
*   **SciPy (`scipy.integrate.solve_ivp`):** Provides robust, battle-tested ODE solvers.
*   **Numba (`@jit`):** Compiles the derivative function (`dydt`) to optimized machine code (LLVM) at runtime, achieving C++-like performance.
*   **JAX (Fallback):** Reserved as a "Plan B" only if gradient-based optimization or GPU acceleration becomes strictly necessary.

**Consequences:**
*   (+) Performance: 10-100x speedup over pure Python.
*   (+) Maintainability: Code remains readable Python, avoiding complex C++ wrappers.
*   (-) Constraint: Numba-compatible functions cannot use some Python features (e.g., dynamic dictionaries inside loops).
