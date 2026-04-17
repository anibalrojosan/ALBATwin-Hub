# Testing

How to install dev dependencies and run the test suite for **bioprocess-twin-hub**.

## Prerequisites

From the repository root, with [uv](https://docs.astral.sh/uv/) installed:

```bash
uv sync --group dev
uv pip install -e .
```

Python **3.12+** is required (see `.python-version`).

## Default test run (CI-equivalent gate)

[`pyproject.toml`](../../pyproject.toml) configures pytest with:

- Coverage on `src/bioprocess_twin` and a terminal missing-lines report.
- **Exclusion** of tests marked `strict_si_mass_balance` (literal SI elemental balance on `S @ Iᵀ` at `MASS_BALANCE_ATOL`; many cells fail until extended O/H closure is fully enforced).

```bash
uv run pytest
```

This is what the **`quality`** job runs in [`.github/workflows/ci.yml`](../../.github/workflows/ci.yml).

## Strict SI mass balance audit

To run only the marked tests and see which process×element cells violate the tolerance (useful for closure policy work):

```bash
uv run pytest -m strict_si_mass_balance -o addopts=''
```

`-o addopts=''` overrides the default `addopts` for this invocation so the `-m "not strict_si_mass_balance"` filter from `pyproject.toml` does not combine with `-m strict_si_mass_balance` (which would select nothing).

GitHub Actions also runs this as a separate **`mass-balance-si-audit`** job with `continue-on-error: true` so the branch stays mergeable while the audit remains visible.

## Run a subset of tests

Examples:

```bash
uv run pytest tests/unit/test_kinetics.py -v
uv run pytest tests/unit/test_kinetic_parameters.py tests/unit/test_kinetic_modifiers.py -q
```

To run **all** tests including those normally deselected by the default marker filter, either use the audit command above or temporarily adjust `addopts` / use `-o addopts=''` and pass an explicit file or expression.

## Coverage focused on one module

Override the default `addopts` so only the kinetics module is measured:

```bash
uv run pytest tests/unit/test_kinetics.py \
  -o addopts='--cov=bioprocess_twin.models.kinetics --cov-report=term-missing'
```

## Linting

Not part of pytest, but required before merge in CI:

```bash
uv run ruff check .
uv run ruff format --check .
```

## Related documentation

- [`DEVLOG.md`](DEVLOG.md) — sprint notes and decisions.
- [`docs/MATH_MODEL.md`](../MATH_MODEL.md) — process rates $\rho_1\ldots\rho_{19}$ and kinetic SSOT.
- [`docs/STOICHIOMETRY.md`](../STOICHIOMETRY.md) — Petersen matrix and closure modes.
