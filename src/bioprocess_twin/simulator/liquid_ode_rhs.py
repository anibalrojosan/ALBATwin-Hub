"""ODE-sized liquid RHS: Stage 6 only, with diel forcing (phase1-04a) and no transport (04c)."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from bioprocess_twin.forcing.diel_forcing_schedule import DielForcingSchedule, to_env_conditions
from bioprocess_twin.models.chemistry import AlbaDissociationConstantsRef, AlbaDissociationEnthalpy, PHSolverOptions
from bioprocess_twin.models.gas_transfer import GasTransferConditions
from bioprocess_twin.models.kinetic_parameters import KineticParameters

from .liquid_rhs import evaluate_liquid_rhs, state_vector_from_y


@dataclass(frozen=True, slots=True)
class LiquidOdeRhsProblem:
    """
    Immutable bundle for evaluate_liquid_ode_rhs.

    Time: t_hours passed to the RHS is hours of day (clock time). Values outside
    [0,24) are reduced modulo 24 inside DielForcingSchedule.at (same convention as forcing).
    """

    schedule: DielForcingSchedule
    initial_ph: float = 7.0
    placeholder_ph_for_env: float = 7.0
    kinetic_parameters: KineticParameters | None = None
    theta_kla: float | None = None
    gas_conditions: GasTransferConditions | None = None
    delta_cat_an_mol_per_m3: float = 0.0
    options: PHSolverOptions | None = None
    k_ref: AlbaDissociationConstantsRef | None = None
    dh: AlbaDissociationEnthalpy | None = None


def evaluate_liquid_ode_rhs(t_hours: float, y: np.ndarray, *, problem: LiquidOdeRhsProblem) -> np.ndarray:
    """
    SciPy-style vector field: dy/dt in g m⁻³ d⁻¹, shape (N_STATE,).

    Builds EnvConditions from problem.schedule.at(t_hours) and
    to_env_conditions(..., ph=problem.placeholder_ph_for_env). Solved pH inside
    evaluate_liquid_rhs still comes from SI.6 charge balance (initial_ph is only the
    solver guess).
    """
    st = state_vector_from_y(y)
    sample = problem.schedule.at(t_hours)
    env = to_env_conditions(sample, ph=problem.placeholder_ph_for_env)
    out = evaluate_liquid_rhs(
        st,
        env,
        kinetic_parameters=problem.kinetic_parameters,
        theta_kla=problem.theta_kla,
        gas_conditions=problem.gas_conditions,
        delta_cat_an_mol_per_m3=problem.delta_cat_an_mol_per_m3,
        initial_ph=problem.initial_ph,
        options=problem.options,
        k_ref=problem.k_ref,
        dh=problem.dh,
    )
    return np.asarray(out.dcdt_g_m3_d, dtype=np.float64).ravel()


def make_liquid_rhs(problem: LiquidOdeRhsProblem) -> Callable[[float, np.ndarray], np.ndarray]:
    """Return rhs(t, y) closed over problem for use with solve_ivp."""

    def rhs(t_hours: float, y: np.ndarray) -> np.ndarray:
        return evaluate_liquid_ode_rhs(t_hours, y, problem=problem)

    return rhs
