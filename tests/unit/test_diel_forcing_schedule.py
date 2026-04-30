"""Tests for diel forcing schedules and ``EnvConditions`` mapping."""

from __future__ import annotations

import numpy as np
import pytest

from bioprocess_twin.forcing import (
    DielForcingSchedule,
    dense_time_grid,
    evaluate_forcing,
    to_env_conditions,
)
from bioprocess_twin.models.kinetics import EnvConditions


def test_summer_irradiance_zero_at_night_and_nonneg_dense_grid() -> None:
    sched = DielForcingSchedule(season="summer")
    for t in (0.0, 12.0, 20.0):
        s = sched.at(t)
        assert np.isfinite(s.temperature_C)
        if t >= 20.0:
            assert s.irradiance_umol_m2_s == pytest.approx(0.0)
        assert s.irradiance_umol_m2_s >= 0.0

    t_dense = dense_time_grid(200)
    arr = sched.at_many(t_dense)
    assert arr.temperature_C.shape == t_dense.shape
    assert np.all(arr.irradiance_umol_m2_s >= 0.0)


def test_at_many_matches_evaluate_forcing_without_optionals() -> None:
    t = np.linspace(0.0, 24.0, 50, endpoint=False)
    sched = DielForcingSchedule(season="autumn")
    arr = sched.at_many(t)
    irr_e, temp_e, evap_e = evaluate_forcing("autumn", t)
    np.testing.assert_allclose(arr.irradiance_umol_m2_s, irr_e)
    np.testing.assert_allclose(arr.temperature_C, temp_e)
    np.testing.assert_allclose(arr.evaporation_m3_h, evap_e)


def test_optional_constants_broadcast() -> None:
    sched = DielForcingSchedule(
        season="spring",
        relative_humidity_percent=65.0,
        inflow_m3_h=1.0,
    )
    a = sched.at(3.0)
    b = sched.at(15.0)
    assert a.relative_humidity_percent == pytest.approx(65.0)
    assert b.relative_humidity_percent == pytest.approx(65.0)
    assert a.inflow_m3_h == pytest.approx(1.0)
    assert b.inflow_m3_h == pytest.approx(1.0)


def test_optional_callable_wind() -> None:
    sched = DielForcingSchedule(season="winter", wind_speed_m_s=lambda t: float(t))
    s = sched.at(6.0)
    assert s.wind_speed_m_s == pytest.approx(6.0)


def test_evaporation_override_and_optional_arrays_callable() -> None:
    sched_flat = DielForcingSchedule(season="spring", evaporation_m3_h=0.05)
    s = sched_flat.at(8.0)
    assert s.evaporation_m3_h == pytest.approx(0.05)

    t = np.linspace(0.0, 24.0, 13, endpoint=False)
    sched_rain = DielForcingSchedule(season="summer", rain_mm_h=lambda tt: 0.01 * tt)
    arr = sched_rain.at_many(t)
    assert arr.rain_mm_h is not None
    np.testing.assert_allclose(arr.rain_mm_h, 0.01 * t)


def test_to_env_conditions_mapping() -> None:
    sched = DielForcingSchedule(season="summer")
    sample = sched.at(10.0)
    env = to_env_conditions(sample, ph=7.2)
    assert isinstance(env, EnvConditions)
    assert env.pH == pytest.approx(7.2)
    assert env.temperature_C == pytest.approx(sample.temperature_C)
    assert env.irradiance_umol_m2_s == pytest.approx(sample.irradiance_umol_m2_s)
