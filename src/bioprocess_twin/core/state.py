from enum import StrEnum
from typing import Literal, Self

import numpy as np
from pydantic import BaseModel, ConfigDict, Field


class StateVectorVariant(StrEnum):
    """
    Which simulation stoichiometry layout this state is paired with.

    Values align with ``PetersenClosureMode`` in ``stoichiometry_closure.py``:
    SI and oxygen closure use **17** ODE components; oxygen+proton closure uses **18**
    (append ``S_H_PROTON`` after ``S_H2O``).
    """

    SI = "si"
    OXYGEN_CLOSURE = "oxygen"
    OXYGEN_AND_PROTON_CLOSURE = "oxygen_and_protons"


def state_array_len(variant: StateVectorVariant) -> Literal[17, 18]:
    """Length of ``StateVector.to_array(..., variant=variant)``."""
    if variant == StateVectorVariant.OXYGEN_AND_PROTON_CLOSURE:
        return 18
    return 17


class StateVector(BaseModel):
    """
    State vector for the ALBA bioprocess model.

    Holds **18** named fields; the first **17** match Casagli SI order. ``S_H_PROTON``
    is the optional 18th compartment (g H as free-proton inventory) used only when
    ``variant=StateVectorVariant.OXYGEN_AND_PROTON_CLOSURE`` in ``to_array`` /
    ``from_array``.

    For SI and oxygen-closure runs, use ``to_array(variant=...)`` with a **17** layout;
    ``S_H_PROTON`` is omitted from the array and ``from_array`` with 17 values sets it
    to **0**.
    """

    model_config = ConfigDict(frozen=True)

    X_ALG: float = Field(..., ge=0.0, description="Algal biomass [gCOD/m3]")
    X_AOB: float = Field(..., ge=0.0, description="Ammonia Oxidizing Bacteria [gCOD/m3]")
    X_NOB: float = Field(..., ge=0.0, description="Nitrite Oxidizing Bacteria [gCOD/m3]")
    X_H: float = Field(..., ge=0.0, description="Heterotrophic Bacteria [gCOD/m3]")
    X_S: float = Field(..., ge=0.0, description="Slowly biodegradable organic matter [gCOD/m3]")
    X_I: float = Field(..., ge=0.0, description="Inert particulate organic matter [gCOD/m3]")

    S_S: float = Field(..., ge=0.0, description="Readily biodegradable organic matter [gCOD/m3]")
    S_I: float = Field(..., ge=0.0, description="Inert soluble organic matter [gCOD/m3]")
    S_IC: float = Field(..., ge=0.0, description="Total Inorganic Carbon [gC/m3]")
    S_ND: float = Field(..., ge=0.0, description="Soluble Organic Nitrogen (Urea) [gN/m3]")
    S_NH: float = Field(..., ge=0.0, description="Total Ammoniacal Nitrogen [gN/m3]")
    S_NO2: float = Field(..., ge=0.0, description="Total Nitrite Nitrogen [gN/m3]")
    S_NO3: float = Field(..., ge=0.0, description="Total Nitrate Nitrogen [gN/m3]")
    S_N2: float = Field(..., ge=0.0, description="Nitrogen Gas [gN/m3]")
    S_PO4: float = Field(..., ge=0.0, description="Total Inorganic Phosphorus [gP/m3]")
    S_O2: float = Field(..., ge=0.0, description="Dissolved Oxygen [gO2/m3]")
    S_H2O: float = Field(..., ge=0.0, description="Water (Balance term) [gH/m3]")
    S_H_PROTON: float = Field(
        default=0.0,
        ge=0.0,
        description="Free proton pool (balance term) [gH/m3]; 18th ODE component when using proton closure",
    )

    def to_array(self, *, variant: StateVectorVariant = StateVectorVariant.SI) -> np.ndarray:
        """
        Convert to a 1D array of length **17** or **18** depending on ``variant``.

        Order: indices ``0..16`` match ``stoichiometry.N_STATE`` / Petersen columns;
        index **17** is ``S_H_PROTON`` when ``variant`` is ``OXYGEN_AND_PROTON_CLOSURE``.
        """
        base = [
            self.X_ALG,
            self.X_AOB,
            self.X_NOB,
            self.X_H,
            self.X_S,
            self.X_I,
            self.S_S,
            self.S_I,
            self.S_IC,
            self.S_ND,
            self.S_NH,
            self.S_NO2,
            self.S_NO3,
            self.S_N2,
            self.S_PO4,
            self.S_O2,
            self.S_H2O,
        ]
        if state_array_len(variant) == 18:
            return np.array(base + [self.S_H_PROTON], dtype=float)
        return np.array(base, dtype=float)

    @classmethod
    def from_array(
        cls,
        arr: np.ndarray,
        *,
        variant: StateVectorVariant | None = None,
    ) -> Self:
        """
        Build from a 1D array. If ``variant`` is ``None``, infer **17** vs **18** from
        ``arr.size`` (``17`` or ``18`` only).

        17-component loads set ``S_H_PROTON`` to **0**.
        """
        flat = np.asarray(arr, dtype=float).ravel()
        if variant is None:
            if flat.size == 18:
                variant = StateVectorVariant.OXYGEN_AND_PROTON_CLOSURE
            elif flat.size == 17:
                variant = StateVectorVariant.SI
            else:
                raise ValueError(
                    f"Expected array length 17 or 18 when variant is omitted, got {flat.size}"
                )
        n = state_array_len(variant)
        if flat.size != n:
            raise ValueError(f"Expected array length {n} for variant {variant!r}, got {flat.size}")

        kw = {
            "X_ALG": float(flat[0]),
            "X_AOB": float(flat[1]),
            "X_NOB": float(flat[2]),
            "X_H": float(flat[3]),
            "X_S": float(flat[4]),
            "X_I": float(flat[5]),
            "S_S": float(flat[6]),
            "S_I": float(flat[7]),
            "S_IC": float(flat[8]),
            "S_ND": float(flat[9]),
            "S_NH": float(flat[10]),
            "S_NO2": float(flat[11]),
            "S_NO3": float(flat[12]),
            "S_N2": float(flat[13]),
            "S_PO4": float(flat[14]),
            "S_O2": float(flat[15]),
            "S_H2O": float(flat[16]),
            "S_H_PROTON": float(flat[17]) if n == 18 else 0.0,
        }
        return cls(**kw)
