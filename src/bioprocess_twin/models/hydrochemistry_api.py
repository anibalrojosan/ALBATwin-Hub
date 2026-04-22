"""
Stage 5: simulator-facing hydrochemistry API (SI.6 pH + SI.7 gas transfer).

These functions are the thin entry points intended for the future ODE RHS
(phase1-04 Orchestration). They read StateVector and EnvConditions only.

Contract: SI.6 charge balance and solve_pH use env.temperature_C for van't Hoff
and speciation totals derived from state. env.irradiance_umol_m2_s does not
enter SI.6. EnvConditions.pH is not used here; it remains the cardinal pH input
for calculate_rates until Sprint 4 wires solved pH into kinetics each step.
"""

from __future__ import annotations

from dataclasses import dataclass

from bioprocess_twin.core.state import StateVector
from bioprocess_twin.models.chemistry import (
    AlbaDissociationConstantsRef,
    AlbaDissociationEnthalpy,
    ChargeBalanceInputs,
    PHSolveResult,
    PHSolverOptions,
    SpeciationTotals,
    mol_c_per_m3_carbon_from_s_ic,
    mol_n_per_m3_nitrogen_from_s_nh,
    mol_n_per_m3_nitrogen_from_s_no2,
    mol_n_per_m3_nitrogen_from_s_no3,
    mol_p_per_m3_phosphorus_from_s_po4,
    solve_pH,
)
from bioprocess_twin.models.gas_transfer import GasTransferConditions, GasTransferRates, calculate_gas_transfer
from bioprocess_twin.models.kinetics import EnvConditions


def speciation_totals_from_state(state: StateVector) -> SpeciationTotals:
    """
    Map ALBA soluble totals (g m^-3) to SI.6.1 molar totals (mol m^-3) for charge balance.

    Same mapping as speciate_from_alba_totals uses internally (S_IC/12, S_NH/14, ...).
    """
    return SpeciationTotals(
        c_tot_nh3=mol_n_per_m3_nitrogen_from_s_nh(float(state.S_NH)),
        c_tot_no2=mol_n_per_m3_nitrogen_from_s_no2(float(state.S_NO2)),
        c_tot_no3=mol_n_per_m3_nitrogen_from_s_no3(float(state.S_NO3)),
        c_tot_ic=mol_c_per_m3_carbon_from_s_ic(float(state.S_IC)),
        c_tot_po4=mol_p_per_m3_phosphorus_from_s_po4(float(state.S_PO4)),
    )


def solve_pH_from_state(
    state: StateVector,
    env_conditions: EnvConditions,
    *,
    delta_cat_an_mol_per_m3: float = 0.0,
    initial_ph: float = 7.0,
    options: PHSolverOptions | None = None,
    k_ref: AlbaDissociationConstantsRef | None = None,
    dh: AlbaDissociationEnthalpy | None = None,
) -> PHSolveResult:
    """
    Solve SI.6.1 row 15 for pH from StateVector and operating temperature.

    Uses env_conditions.temperature_C only (not env_conditions.pH).
    """
    totals = speciation_totals_from_state(state)
    inputs = ChargeBalanceInputs(
        totals=totals,
        t_celsius=env_conditions.temperature_C,
        delta_cat_an_mol_per_m3=delta_cat_an_mol_per_m3,
        k_ref=k_ref,
        dh=dh,
    )
    return solve_pH(inputs, initial_ph=initial_ph, options=options)


@dataclass(frozen=True, slots=True)
class HydrochemistryStepResult:
    """One evaluation: algebraic pH plus SI.7 gas transfer rates at that [H+]."""

    ph_result: PHSolveResult
    gas_rates: GasTransferRates


def hydrochemistry_step(
    state: StateVector,
    env_conditions: EnvConditions,
    *,
    theta_kla: float,
    delta_cat_an_mol_per_m3: float = 0.0,
    initial_ph: float = 7.0,
    options: PHSolverOptions | None = None,
    k_ref: AlbaDissociationConstantsRef | None = None,
    dh: AlbaDissociationEnthalpy | None = None,
    gas_conditions: GasTransferConditions | None = None,
) -> HydrochemistryStepResult:
    """
    Run solve_pH_from_state then calculate_gas_transfer with the solved [H+].

    Gas Henry and K_a(T) use env_conditions.temperature_C. Pass theta_kla from
    KineticParameters.theta_kla when wiring the simulator.
    """
    ph_result = solve_pH_from_state(
        state,
        env_conditions,
        delta_cat_an_mol_per_m3=delta_cat_an_mol_per_m3,
        initial_ph=initial_ph,
        options=options,
        k_ref=k_ref,
        dh=dh,
    )
    gas_rates = calculate_gas_transfer(
        state,
        env_conditions,
        ph_result.h_plus_mol_per_m3,
        theta_kla=theta_kla,
        gas_conditions=gas_conditions,
    )
    return HydrochemistryStepResult(ph_result=ph_result, gas_rates=gas_rates)
