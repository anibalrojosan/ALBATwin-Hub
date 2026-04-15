"""Tolerances for Petersen vs composition matrix mass-balance audits (see docs/MASS_BALANCES.md)."""

# Target: tighten toward 1e-6 once O/H closure (ADR 007) is implemented in stoichiometry.
MASS_BALANCE_ATOL = 0.01
