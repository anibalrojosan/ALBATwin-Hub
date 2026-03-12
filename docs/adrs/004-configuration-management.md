# ADR 004: Configuration Management

**Context:**
The model relies on ~50 parameters ($\mu_{max}$, $K_S$, yields). Hardcoding these in Python makes the system rigid and error-prone.

**Decision:**
Use **YAML files with Pydantic validation**.
*   **YAML:** Human-readable format for configuration (`constants.yaml`, `reactor.yaml`).
*   **Pydantic:** Strict schema validation at load time. Ensures types (float vs int), ranges ($pH > 0$), and required fields.

**Consequences:**
*   (+) Safety: Prevents silent errors (e.g., negative growth rates).
*   (+) Flexibility: Easy to swap parameter sets for different algae strains without changing code.
