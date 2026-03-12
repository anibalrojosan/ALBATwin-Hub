# ADR 005: Data Persistence Format

**Context:**
Simulations generate large time-series datasets (e.g., 1 year at 1-minute intervals = ~500k rows $\times$ 17 columns). CSV files are slow to read/write and lose type information.

**Decision:**
Use **Apache Parquet** with Snappy compression.

**Consequences:**
*   (+) Performance: 10-50x faster I/O than CSV.
*   (+) Efficiency: Significantly smaller file size (columnar compression).
*   (+) Typing: Preserves data types (float32/64), critical for ML training data.
