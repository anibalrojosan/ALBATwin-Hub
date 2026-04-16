# Mass balances (oxygen / water closure): moved

The **oxygen-closure** 114-cell audit now lives at:

**[mass_balances/artifacts/audit-oxygen-closure-114cell.md](mass_balances/artifacts/audit-oxygen-closure-114cell.md)**

Regenerate:

```bash
uv run python scripts/generate_mass_balances_md.py --closure-of-oxygen
```

(`--closure` remains an alias for `--closure-of-oxygen`.)

See **[mass_balances/README.md](mass_balances/README.md)** for context and narrative docs.
