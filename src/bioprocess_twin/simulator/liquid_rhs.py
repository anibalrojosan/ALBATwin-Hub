"""Stage 6 liquid-phase RHS assembler for ALBA SI 17-state layout."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from bioprocess_twin.core.state import StateVector
from bioprocess_twin.models.chemistry import (
    AlbaDissociationConstantsRef,
    AlbaDissociationEnthalpy,
    AqueousSpeciesMolar,
    PHSolveResult,
    PHSolverOptions,
    default_dissociation_constants_ref_molar,
    default_dissociation_enthalpy_j_per_mol,
    scale_dissociation_constants_at_t,
    speciate_from_alba_totals,
)
from bioprocess_twin.models.gas_transfer import GasTransferConditions, GasTransferRates
from bioprocess_twin.models.hydrochemistry_api import hydrochemistry_step
from bioprocess_twin.models.kinetic_parameters import KineticParameters, default_alba
from bioprocess_twin.models.kinetics import EnvConditions, calculate_rates
from bioprocess_twin.models.stoichiometry import (
    N_PROCESSES,
    N_PROCESSES_WITH_GAS_TRANSFER,
    N_STATE,
    get_petersen_matrix_with_gas_transfer,
)


@dataclass(frozen=True, slots=True)
class BiomassConcentrations:
    """Biomass pools from StateVector in gCOD m^-3."""

    x_alg: float
    x_aob: float
    x_nob: float
    x_h: float
    x_s: float
    x_i: float

    @property
    def total_particulate_cod(self) -> float:
        """Total particulate COD-like mass proxy from biomass and particulate pools."""
        return self.x_alg + self.x_aob + self.x_nob + self.x_h + self.x_s + self.x_i


@dataclass(frozen=True, slots=True)
class EnvironmentalSnapshot:
    """Environment and pH values used for one RHS evaluation."""

    temperature_c: float
    irradiance_umol_m2_s: float
    ph_input: float
    ph_solved: float


@dataclass(frozen=True, slots=True)
class LiquidRhsDiagnostics:
    """Diagnostics payload for plotting and verification of one RHS call."""

    state: StateVector
    biomass: BiomassConcentrations
    environment: EnvironmentalSnapshot
    h_plus_mol_per_m3: float
    aqueous_species_molar: AqueousSpeciesMolar


@dataclass(frozen=True, slots=True)
class AlbaLiquidRhsResult:
    """RHS outputs for one liquid-phase evaluation at fixed inputs."""

    dcdt_g_m3_d: np.ndarray
    rho_bio_g_m3_d: np.ndarray
    rho_gas_g_m3_d: np.ndarray
    rho_full_g_m3_d: np.ndarray
    ph_result: PHSolveResult
    gas_rates: GasTransferRates
    diagnostics: LiquidRhsDiagnostics


def state_vector_from_y(state: StateVector | np.ndarray) -> StateVector:
    """Coerce SI-layout length-17 ``ndarray`` (or ``StateVector``) to ``StateVector``."""
    if isinstance(state, StateVector):
        return state
    arr = np.asarray(state, dtype=np.float64).ravel()
    if arr.size != N_STATE:
        raise ValueError(f"state must have length {N_STATE} for Stage 6 SI layout, got {arr.size}")
    return StateVector.from_array(arr)


def _build_biomass_snapshot(state: StateVector) -> BiomassConcentrations:
    """Extract biomass and particulate pools from StateVector."""
    return BiomassConcentrations(
        x_alg=float(state.X_ALG),
        x_aob=float(state.X_AOB),
        x_nob=float(state.X_NOB),
        x_h=float(state.X_H),
        x_s=float(state.X_S),
        x_i=float(state.X_I),
    )


def evaluate_liquid_rhs(
    state: StateVector | np.ndarray,
    env_conditions: EnvConditions,
    *,
    kinetic_parameters: KineticParameters | None = None,
    theta_kla: float | None = None,
    gas_conditions: GasTransferConditions | None = None,
    delta_cat_an_mol_per_m3: float = 0.0,
    initial_ph: float = 7.0,
    options: PHSolverOptions | None = None,
    k_ref: AlbaDissociationConstantsRef | None = None,
    dh: AlbaDissociationEnthalpy | None = None,
) -> AlbaLiquidRhsResult:
    """
    Evaluate one Stage 6 RHS with biological and gas-transfer processes.

    Workflow:
    - solve pH from state and temperature
    - align kinetics pH to solved pH
    - compute rho_bio (19) and rho_gas (3)
    - assemble rho_full (22) and dC/dt (17) via S_full^T @ rho_full
    """
    st = state_vector_from_y(state)
    params = kinetic_parameters if kinetic_parameters is not None else default_alba()
    theta_value = float(theta_kla) if theta_kla is not None else float(params.theta_kla)

    step = hydrochemistry_step(
        st,
        env_conditions,
        theta_kla=theta_value,
        delta_cat_an_mol_per_m3=delta_cat_an_mol_per_m3,
        initial_ph=initial_ph,
        options=options,
        k_ref=k_ref,
        dh=dh,
        gas_conditions=gas_conditions,
    )

    env_aligned = EnvConditions(
        temperature_C=env_conditions.temperature_C,
        pH=step.ph_result.ph,
        irradiance_umol_m2_s=env_conditions.irradiance_umol_m2_s,
    )
    rho_bio = calculate_rates(st, env_aligned, params)

    rho_gas = np.array(
        [
            step.gas_rates.rho_o2_g_m3_d,
            step.gas_rates.rho_co2_gC_m3_d,
            step.gas_rates.rho_nh3_gN_m3_d,
        ],
        dtype=np.float64,
    )
    rho_full = np.concatenate((rho_bio, rho_gas))
    if rho_full.size != N_PROCESSES_WITH_GAS_TRANSFER:
        raise RuntimeError(f"rho_full must have length {N_PROCESSES_WITH_GAS_TRANSFER}, got {rho_full.size}")
    if rho_bio.size != N_PROCESSES:
        raise RuntimeError(f"rho_bio must have length {N_PROCESSES}, got {rho_bio.size}")

    s_full = get_petersen_matrix_with_gas_transfer()
    dcdt = s_full.T @ rho_full

    k_ref_eval = k_ref if k_ref is not None else default_dissociation_constants_ref_molar()
    dh_eval = dh if dh is not None else default_dissociation_enthalpy_j_per_mol()
    k_at_t = scale_dissociation_constants_at_t(k_ref_eval, dh_eval, env_conditions.temperature_C)
    aqueous_species = speciate_from_alba_totals(
        h_plus_mol_per_m3=step.ph_result.h_plus_mol_per_m3,
        s_ic_g_per_m3=float(st.S_IC),
        s_nh_g_per_m3=float(st.S_NH),
        s_no2_g_per_m3=float(st.S_NO2),
        s_no3_g_per_m3=float(st.S_NO3),
        s_po4_g_per_m3=float(st.S_PO4),
        k=k_at_t,
    )

    diagnostics = LiquidRhsDiagnostics(
        state=st,
        biomass=_build_biomass_snapshot(st),
        environment=EnvironmentalSnapshot(
            temperature_c=env_conditions.temperature_C,
            irradiance_umol_m2_s=env_conditions.irradiance_umol_m2_s,
            ph_input=env_conditions.pH,
            ph_solved=step.ph_result.ph,
        ),
        h_plus_mol_per_m3=step.ph_result.h_plus_mol_per_m3,
        aqueous_species_molar=aqueous_species,
    )
    return AlbaLiquidRhsResult(
        dcdt_g_m3_d=dcdt,
        rho_bio_g_m3_d=rho_bio,
        rho_gas_g_m3_d=rho_gas,
        rho_full_g_m3_d=rho_full,
        ph_result=step.ph_result,
        gas_rates=step.gas_rates,
        diagnostics=diagnostics,
    )
