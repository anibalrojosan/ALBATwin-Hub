"""
Typical daily environmental series by season, aligned with Casagli et al. (2021) Fig. 1.

Construction: visual control points plus monotone piecewise cubic
interpolation (`scipy.interpolate.PchipInterpolator`, Fritsch–Carlson). 

Irradiance values at t ≥ 20 h are forced to 0 (night). Control knots at 20 h use 0 so the daytime curve 
meets night without relying on the sparse 20–24 h segment alone. 

Non-negative variables use `np.clip(..., 0, None)` after evaluation for tiny negative numerical noise.

Hourly tables for the simulator are obtained by evaluating the interpolants at `t_h = 0 … 23`.

Units:
  - irradiance: µmol m⁻² s⁻¹ (figure axis often prints µE m⁻² s⁻¹; same numeric scale for PAR),
  - temperature: °C,
  - evaporation: m³ h⁻¹.

Usage
-----
    uv run python -m bioprocess_twin.forcing.typical_daily_forcing_per_season

writes ``data/forcing/typical_daily_patterns_of_temperature_irradiance_and_evaporation_rates.csv``

Import from notebooks or code::

    from bioprocess_twin.forcing import typical_daily_forcing_per_season as forcing

See ``notebooks/fig1_forcing_visual_check.ipynb`` for plots vs control points.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

Season = Literal["spring", "summer", "autumn", "winter"]

_REPO_ROOT = Path(__file__).resolve().parents[3]

DEFAULT_EXPORT_CSV = "typical_daily_patterns_of_temperature_irradiance_and_evaporation_rates.csv"

# Irradiance is forced to zero for t ≥ this hour (night).
IRRADIANCE_NIGHT_FROM_H = 20.0

# --- Control-point knots [h] (visual extraction from Fig. 1) ---
T_CTRL_IRR = np.array([0.0, 4.0, 6.0, 8.0, 10.0, 12.0, 13.0, 14.0, 16.0, 18.0, 20.0, 24.0])
T_CTRL_TEMP = np.array([0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0])
T_CTRL_EVAP = np.array([0.0, 4.0, 8.0, 12.0, 14.0, 16.0, 20.0, 24.0])

# Irradiance [µmol m⁻² s⁻¹] at T_CTRL_IRR
IRRADIANCE_CONTROL: dict[Season, np.ndarray] = {
    "spring": np.array([0, 0, 0, 100, 750, 1300, 1380, 1250, 750, 100, 0, 0], dtype=np.float64),
    "summer": np.array([0, 0, 50, 250, 950, 1550, 1650, 1500, 950, 250, 0, 0], dtype=np.float64),
    "autumn": np.array([0, 0, 0, 50, 350, 750, 800, 600, 300, 0, 0, 0], dtype=np.float64),
    "winter": np.array([0, 0, 0, 0, 300, 750, 800, 650, 250, 0, 0, 0], dtype=np.float64),
}

# Temperature [°C] at T_CTRL_TEMP
TEMPERATURE_CONTROL: dict[Season, np.ndarray] = {
    "spring": np.array([20.8, 19.5, 18.8, 20.5, 24.2, 22.5, 20.8], dtype=np.float64),
    "summer": np.array([23.5, 22.0, 21.2, 23.0, 27.2, 25.5, 23.5], dtype=np.float64),
    "autumn": np.array([14.6, 14.0, 13.5, 14.2, 16.3, 15.2, 14.6], dtype=np.float64),
    "winter": np.array([8.0, 7.0, 6.5, 7.5, 10.1, 9.0, 8.0], dtype=np.float64),
}

# Evaporation [m³ h⁻¹] at T_CTRL_EVAP
EVAPORATION_CONTROL: dict[Season, np.ndarray] = {
    "spring": np.array([0.005, 0.004, 0.009, 0.026, 0.029, 0.025, 0.008, 0.005], dtype=np.float64),
    "summer": np.array([0.008, 0.006, 0.015, 0.042, 0.051, 0.040, 0.015, 0.008], dtype=np.float64),
    "autumn": np.array([0.004, 0.003, 0.004, 0.013, 0.017, 0.012, 0.005, 0.004], dtype=np.float64),
    "winter": np.array([0.004, 0.003, 0.003, 0.014, 0.019, 0.013, 0.004, 0.004], dtype=np.float64),
}

# Matplotlib styles approximating the paper (colour + linestyle); keys = Season.
FIG1_LINE_STYLES: dict[Season, dict[str, object]] = {
    "spring": {"color": "#cc0000", "linestyle": (0, (5, 5)), "linewidth": 2.0},
    "summer": {"color": "#ff9999", "linestyle": "--", "linewidth": 2.0},
    "autumn": {"color": "#0000cc", "linestyle": ":", "linewidth": 2.5},
    "winter": {"color": "#9999ff", "linestyle": "-", "linewidth": 2.0},
}


def _interp_clip(
    t_ctrl: np.ndarray,
    y_ctrl: np.ndarray,
    t_eval: np.ndarray,
    *,
    clip_nonnegative: bool,
) -> np.ndarray:
    """Piecewise cubic Hermite (PCHIP) preserving monotonicity between consecutive knots."""
    tx = np.asarray(t_ctrl, dtype=np.float64)
    yy = np.asarray(y_ctrl, dtype=np.float64)
    pchip = PchipInterpolator(tx, yy, extrapolate=False)
    y = np.asarray(pchip(t_eval), dtype=np.float64)
    if clip_nonnegative:
        y = np.clip(y, 0.0, None)
    return y


def evaluate_forcing(
    season: Season,
    t_hours: np.ndarray,
    *,
    clip_irradiance: bool = True,
    clip_evaporation: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Evaluate irradiance, temperature, evaporation at arbitrary times in [0, 24].

    Uses monotone PCHIP interpolants through the Fig. 1 control points for each variable.
    Irradiance is set to **0** for **t ≥ IRRADIANCE_NIGHT_FROM_H** (default 20 h).
    """
    t_eval = np.asarray(t_hours, dtype=np.float64)
    irr = _interp_clip(T_CTRL_IRR, IRRADIANCE_CONTROL[season], t_eval, clip_nonnegative=clip_irradiance)
    irr = np.where(t_eval >= IRRADIANCE_NIGHT_FROM_H, 0.0, irr)
    temp = _interp_clip(T_CTRL_TEMP, TEMPERATURE_CONTROL[season], t_eval, clip_nonnegative=False)
    evap = _interp_clip(T_CTRL_EVAP, EVAPORATION_CONTROL[season], t_eval, clip_nonnegative=clip_evaporation)
    return irr, temp, evap


def dense_time_grid(n_points: int = 200) -> np.ndarray:
    """Uniformly spaced times in [0, 24] for smooth plotting."""
    return np.linspace(0.0, 24.0, int(n_points))


def hourly_array_t_h() -> np.ndarray:
    return np.arange(0, 24, dtype=np.float64)


def typical_daily_forcing_table() -> pd.DataFrame:
    """Build 96 rows: 4 seasons × 24 hours (t_h = 0 … 23)."""
    hours = hourly_array_t_h()
    rows: list[dict[str, float | str]] = []
    for season in ("spring", "summer", "autumn", "winter"):
        irr, temp, evap = evaluate_forcing(season, hours)
        for i, t_h in enumerate(hours):
            rows.append(
                {
                    "season": season,
                    "t_h": float(t_h),
                    "irradiance_umol_m2_s": float(irr[i]),
                    "temperature_C": float(temp[i]),
                    "evaporation_m3_h": float(evap[i]),
                }
            )
    return pd.DataFrame(rows)


def typical_daily_forcing_numpy(season: Season) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (t_h, irradiance_umol_m2_s, temperature_C, evaporation_m3_h), each shape (24,)."""
    df = typical_daily_forcing_table()
    sub = df.loc[df["season"] == season].sort_values(by="t_h")
    return (
        sub["t_h"].to_numpy(),
        sub["irradiance_umol_m2_s"].to_numpy(),
        sub["temperature_C"].to_numpy(),
        sub["evaporation_m3_h"].to_numpy(),
    )


def write_csv(path: Path | None = None) -> Path:
    """Write the full 4×24 table to CSV under ``data/forcing/`` by default."""
    out = path if path is not None else _REPO_ROOT / "data" / "forcing" / DEFAULT_EXPORT_CSV
    out.parent.mkdir(parents=True, exist_ok=True)
    df = typical_daily_forcing_table().round(
        {"t_h": 0, "irradiance_umol_m2_s": 2, "temperature_C": 4, "evaporation_m3_h": 6}
    )
    df.to_csv(out, index=False)
    return out


def main() -> None:
    out = write_csv()
    print(f"Wrote {out} ({len(typical_daily_forcing_table())} rows)")


if __name__ == "__main__":
    main()
