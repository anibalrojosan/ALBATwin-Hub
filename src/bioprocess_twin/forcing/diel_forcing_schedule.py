"""
Diel forcing schedules built on Casagli et al. (2021) Fig. 1 series plus optional drivers.

EnvConditions (kinetics) only carries temperature, pH, and irradiance. Additional
quantities (RH, wind, rain, evaporation override, inflow Q) are carried on
ForcingSample / DielForcingArrays for upcoming hydraulics (e.g. phase1-04c) and
are not consumed by calculate_rates or calculate_gas_transfer yet.

Units:
  - rain_mm_h: mm per hour (rain depth rate).
  - evaporation_m3_h: m³ h⁻¹ (from Fig. 1 digitization when not overridden).
  - inflow_m3_h: volumetric inflow Q [m³ h⁻¹].
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import TypeAlias

import numpy as np

from bioprocess_twin.models.kinetics import EnvConditions

from .typical_daily_forcing_per_season import (
    Season,
    evaluate_forcing,
)

OptionalDriver: TypeAlias = float | Callable[[float], float] | None


def _wrap_hours_scalar(t_hours: float) -> float:
    """Map clock time to $[0, 24)$ for periodic daily evaluation."""
    return float(np.mod(t_hours, 24.0))


def _wrap_hours_array(t_hours: np.ndarray) -> np.ndarray:
    return np.mod(np.asarray(t_hours, dtype=np.float64), 24.0)


def _eval_optional(opt: OptionalDriver, t_hours: float) -> float | None:
    if opt is None:
        return None
    if callable(opt):
        return float(opt(_wrap_hours_scalar(t_hours)))
    return float(opt)


def _eval_optional_array(
    opt: OptionalDriver,
    t_wrapped: np.ndarray,
) -> np.ndarray | None:
    if opt is None:
        return None
    if callable(opt):
        return np.vectorize(lambda tx: float(opt(float(tx))), otypes=[np.float64])(t_wrapped)
    return np.full_like(t_wrapped, float(opt), dtype=np.float64)


@dataclass(frozen=True, slots=True)
class ForcingSample:
    """Environmental snapshot at one clock time within a repeating daily cycle."""

    t_hours: float
    temperature_C: float
    irradiance_umol_m2_s: float
    evaporation_m3_h: float
    relative_humidity_percent: float | None = None
    wind_speed_m_s: float | None = None
    rain_mm_h: float | None = None
    inflow_m3_h: float | None = None


@dataclass(frozen=True, slots=True)
class DielForcingArrays:
    """Vectorized diel series on t_hours (same shape for all array fields)."""

    t_hours: np.ndarray
    temperature_C: np.ndarray
    irradiance_umol_m2_s: np.ndarray
    evaporation_m3_h: np.ndarray
    relative_humidity_percent: np.ndarray | None = None
    wind_speed_m_s: np.ndarray | None = None
    rain_mm_h: np.ndarray | None = None
    inflow_m3_h: np.ndarray | None = None


@dataclass(frozen=True, slots=True)
class DielForcingSchedule:
    """
    Fig. 1–based T(t), I0(t), and default evaporation; optional extra drivers.

    Optional keyword arguments accept None (not provided), a constant float,
    or a callable f(t_hours) -> float for a time-varying schedule on [0,24).
    evaporation_m3_h when set overrides the Fig. 1 evaporation curve; when
    None, evaporation comes from evaluate_forcing.
    """

    season: Season
    relative_humidity_percent: OptionalDriver = None
    wind_speed_m_s: OptionalDriver = None
    rain_mm_h: OptionalDriver = None
    evaporation_m3_h: OptionalDriver = None
    inflow_m3_h: OptionalDriver = None

    def at(self, t_hours: float) -> ForcingSample:
        tw = _wrap_hours_scalar(t_hours)
        irr, temp, evap_base = evaluate_forcing(self.season, np.array([tw]))
        irr_v, temp_v, evap_base_v = float(irr[0]), float(temp[0]), float(evap_base[0])
        evap_use = _eval_optional(self.evaporation_m3_h, t_hours)
        evap_final = evap_base_v if evap_use is None else evap_use
        return ForcingSample(
            t_hours=tw,
            temperature_C=temp_v,
            irradiance_umol_m2_s=irr_v,
            evaporation_m3_h=evap_final,
            relative_humidity_percent=_eval_optional(self.relative_humidity_percent, t_hours),
            wind_speed_m_s=_eval_optional(self.wind_speed_m_s, t_hours),
            rain_mm_h=_eval_optional(self.rain_mm_h, t_hours),
            inflow_m3_h=_eval_optional(self.inflow_m3_h, t_hours),
        )

    def at_many(self, t_hours: np.ndarray) -> DielForcingArrays:
        """Evaluate schedules on a dense grid (e.g. for a future time integrator)."""
        t = np.asarray(t_hours, dtype=np.float64)
        tw = _wrap_hours_array(t)
        irr, temp, evap_base = evaluate_forcing(self.season, tw)
        evap_override = _eval_optional_array(self.evaporation_m3_h, tw)
        if evap_override is None:
            evap = np.asarray(evap_base, dtype=np.float64)
        else:
            evap = evap_override
        return DielForcingArrays(
            t_hours=t,
            temperature_C=np.asarray(temp, dtype=np.float64),
            irradiance_umol_m2_s=np.asarray(irr, dtype=np.float64),
            evaporation_m3_h=evap,
            relative_humidity_percent=_eval_optional_array(self.relative_humidity_percent, tw),
            wind_speed_m_s=_eval_optional_array(self.wind_speed_m_s, tw),
            rain_mm_h=_eval_optional_array(self.rain_mm_h, tw),
            inflow_m3_h=_eval_optional_array(self.inflow_m3_h, tw),
        )


def to_env_conditions(sample: ForcingSample, *, ph: float) -> EnvConditions:
    """Map kinetic-relevant fields to EnvConditions; ph is required (not from Fig. 1)."""
    return EnvConditions(
        temperature_C=float(sample.temperature_C),
        pH=float(ph),
        irradiance_umol_m2_s=float(sample.irradiance_umol_m2_s),
    )
