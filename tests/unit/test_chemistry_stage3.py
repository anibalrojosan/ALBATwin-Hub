"""Stage 3 tests: SI.6.1 row 15 residual and pH solving."""

from __future__ import annotations

import math

import pytest

from bioprocess_twin.models.chemistry import (
    ChargeBalanceInputs,
    PHSolveError,
    PHSolverOptions,
    SpeciationTotals,
    charge_residual,
    default_dissociation_constants_ref_molar,
    default_dissociation_enthalpy_j_per_mol,
    ph_from_h_plus_mol_per_m3,
    scale_dissociation_constants_at_t,
    solve_pH,
    speciate_aqueous,
)


def test_charge_residual_matches_row15_manual_assembly() -> None:
    """charge_residual reproduces SI.6 row 15 sign convention and valence weights."""
    totals = SpeciationTotals(
        c_tot_nh3=0.015,
        c_tot_no2=0.002,
        c_tot_no3=0.006,
        c_tot_ic=0.035,
        c_tot_po4=0.0012,
    )
    inputs = ChargeBalanceInputs(
        totals=totals,
        t_celsius=25.0,
        delta_cat_an_mol_per_m3=0.0075,
    )
    h = 3.0e-5
    k_t = scale_dissociation_constants_at_t(
        default_dissociation_constants_ref_molar(),
        default_dissociation_enthalpy_j_per_mol(),
        inputs.t_celsius,
    )
    sp = speciate_aqueous(h, totals, k_t)
    expected = (
        h
        + sp.nh4
        + inputs.delta_cat_an_mol_per_m3
        - sp.oh
        - sp.no2
        - sp.no3
        - sp.hco3
        - 2.0 * sp.co3
        - sp.h2po4
        - 2.0 * sp.hpo4
        - 3.0 * sp.po4
    )
    assert math.isclose(charge_residual(h, inputs), expected, rel_tol=0.0, abs_tol=1e-15)


def test_solve_ph_external_reference_pure_water() -> None:
    """External-style reference: SI row-14 convention yields pH ~8.5 for pure water at 25 °C."""
    inputs = ChargeBalanceInputs(
        totals=SpeciationTotals(0.0, 0.0, 0.0, 0.0, 0.0),
        t_celsius=25.0,
        delta_cat_an_mol_per_m3=0.0,
    )
    result = solve_pH(inputs, initial_ph=6.0)
    assert result.converged
    assert result.iterations <= 10
    assert math.isclose(result.ph, 8.5, rel_tol=0.0, abs_tol=5.0e-3)
    assert abs(result.residual) <= 1.0e-10


def test_solve_ph_converges_from_multiple_initial_guesses() -> None:
    """Different initial pH guesses converge to the same root for a fixed chemistry case."""
    inputs = ChargeBalanceInputs(
        totals=SpeciationTotals(
            c_tot_nh3=0.010,
            c_tot_no2=0.0015,
            c_tot_no3=0.0040,
            c_tot_ic=0.040,
            c_tot_po4=0.0010,
        ),
        t_celsius=25.0,
        delta_cat_an_mol_per_m3=0.010,
    )
    result_a = solve_pH(inputs, initial_ph=5.0)
    result_b = solve_pH(inputs, initial_ph=6.5)
    assert result_a.converged and result_b.converged
    assert result_a.iterations <= 10
    assert result_b.iterations <= 10
    assert math.isclose(result_a.ph, result_b.ph, rel_tol=0.0, abs_tol=2.0e-5)
    assert math.isclose(result_a.h_plus_mol_per_m3, result_b.h_plus_mol_per_m3, rel_tol=0.0, abs_tol=1.0e-8)


def test_solve_ph_internal_regression_case() -> None:
    """Pinned synthetic case to detect drifts in charge residual assembly or solver behavior."""
    inputs = ChargeBalanceInputs(
        totals=SpeciationTotals(
            c_tot_nh3=0.008,
            c_tot_no2=0.0005,
            c_tot_no3=0.0025,
            c_tot_ic=0.025,
            c_tot_po4=0.0008,
        ),
        t_celsius=25.0,
        delta_cat_an_mol_per_m3=0.006,
    )
    result = solve_pH(inputs, initial_ph=7.2)
    assert result.converged
    assert result.method_used == "newton_logh"
    assert math.isclose(result.ph, 5.784138, rel_tol=0.0, abs_tol=1.0e-4)


def test_solve_ph_raises_clear_error_when_no_bracket() -> None:
    """Extreme cation offset without balancing anions fails with a clear bracket error."""
    inputs = ChargeBalanceInputs(
        totals=SpeciationTotals(0.0, 0.0, 0.0, 0.0, 0.0),
        t_celsius=25.0,
        delta_cat_an_mol_per_m3=1.0e6,
    )
    with pytest.raises(PHSolveError, match="No sign change in pH bounds"):
        solve_pH(inputs, initial_ph=7.0)


def test_solve_ph_uses_bounded_fallback_when_newton_is_disabled() -> None:
    """Fallback bracket converges when Newton iterations are intentionally skipped."""
    inputs = ChargeBalanceInputs(
        totals=SpeciationTotals(0.0, 0.0, 0.0, 0.0, 0.0),
        t_celsius=25.0,
        delta_cat_an_mol_per_m3=0.0,
    )
    options = PHSolverOptions(max_iterations=0, residual_tol=1.0e-12, x_step_tol=1.0e-12)
    result = solve_pH(inputs, initial_ph=8.5, options=options)
    assert result.converged
    assert result.method_used == "fallback_bracket"
    assert math.isclose(result.ph, 8.5, rel_tol=0.0, abs_tol=5.0e-3)
    assert math.isclose(result.ph, ph_from_h_plus_mol_per_m3(result.h_plus_mol_per_m3), rel_tol=0.0, abs_tol=1e-12)
