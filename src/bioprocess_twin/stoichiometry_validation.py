"""Tolerances for Petersen vs composition matrix mass-balance audits (see docs/mass_balances/README.md)."""

# Target: tighten toward 1e-6 once O/H closure (ADR 007) is implemented in stoichiometry.
MASS_BALANCE_ATOL = 0.01
