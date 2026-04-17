"""
ALBA biological process rates ρ (stub).

Process order matches ``get_petersen_matrix()`` row index ``i`` = ρ_{i+1}.
Full rate expressions: ``docs/MATH_MODEL.md`` §4; environment scalars §3.
"""

import numpy as np
from pydantic import BaseModel, ConfigDict, Field

from bioprocess_twin.core.state import StateVector, StateVectorVariant

from .stoichiometry import N_PROCESSES, N_STATE


class EnvConditions(BaseModel):
    """Environmental inputs for kinetic modifiers (§3) and future rate laws."""

    model_config = ConfigDict(frozen=True)

    temperature_C: float = Field(..., description="Temperature [°C]; may be negative (CTMI T_min for algae).")
    pH: float = Field(..., description="pH for cardinal pH factor (CPM).")
    irradiance_umol_m2_s: float = Field(
        ...,
        ge=0.0,
        description="Photosynthetically active irradiance [μmol m⁻² s⁻¹] for algal f_I (§3.3).",
    )


def calculate_rates(
    state: np.ndarray | StateVector,
    env_conditions: EnvConditions,
) -> np.ndarray:
    """
    Return the vector of **19** process rates ρ₁…ρ₁₉ [d⁻¹] × concentrations as in §4.

    **Stub:** returns zeros; implements contract only.

    State layout: first **17** entries are the Casagli SI variables (same order as
    ``StateVector.to_array`` indices 0..16). A length-**18** ndarray uses only
    ``state[:17]`` (proton closure compartment is ignored for ρ; see ``MATH_MODEL.md`` §2.1).

    ``StateVector`` is converted via ``to_array(variant=SI)`` (17 components).

    Args:
        state: Concentrations ``(17,)`` or ``(18,)``, or a ``StateVector`` instance.
        env_conditions: Temperature, pH, irradiance (unused in stub).

    Returns:
        ``ndarray`` of shape ``(19,)``, ``dtype`` float64.

    Raises:
        ValueError: If ``state`` array size is not 17 or 18.
    """
    _ = env_conditions  # used in full implementation
    if isinstance(state, StateVector):
        conc = state.to_array(variant=StateVectorVariant.SI)
    else:
        flat = np.asarray(state, dtype=np.float64).ravel()
        n = flat.size
        if n == N_STATE:
            conc = flat
        elif n == N_STATE + 1:
            conc = flat[:N_STATE]
        else:
            raise ValueError(f"state as ndarray must have length {N_STATE} or {N_STATE + 1} (proton layout), got {n}")

    _ = conc  # SI concentrations; used when rate laws are implemented

    return np.zeros(N_PROCESSES, dtype=np.float64)
