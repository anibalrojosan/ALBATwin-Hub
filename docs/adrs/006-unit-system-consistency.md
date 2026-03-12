# ADR 006: Unit System Consistency

**Context:**
Bioprocess engineering mixes units dangerously ($mg/L$, $g/m^3$, $M$, $d^{-1}$, $s^{-1}$). Unit conversion errors are a common cause of model failures.

**Decision:**
Standardize strictly on **SI Base Units + Days** for all internal calculations.
*   **Mass:** Grams ($g$) or Grams COD ($gCOD$).
*   **Volume:** Cubic Meters ($m^3$).
*   **Time:** Days ($d$).
*   **Temperature:** Celsius ($^\circ C$).
*   **Exception:** Light in $\mu mol \cdot m^{-2} \cdot s^{-1}$ (standard in photosynthesis models).

**Implementation:**
*   All inputs must be converted to this base system at the boundary (Input Layer).
*   No "magic number" conversions inside the derivative functions.
*   Output data is saved in base units; conversion for visualization happens in the Dashboard layer.

**Consequences:**
*   (+) Clarity: Developers don't need to guess if a value is in $L$ or $m^3$.
*   (+) Safety: Eliminates a whole class of bugs.
