from typing import Self

import numpy as np
from pydantic import BaseModel, ConfigDict, Field


class StateVector(BaseModel):
    """
    State vector for the ALBA bioprocess model.
    Encapsulates the 17 state variables, ensuring non-negative concentrations
    and providing easy conversion to/from NumPy arrays for the ODE solver.
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

    def to_array(self) -> np.ndarray:
        """
        Convert the state vector to a 1D NumPy array for numerical integration.
        The order of variables is strictly maintained.
        """
        return np.array(
            [
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
            ],
            dtype=float,
        )

    @classmethod
    def from_array(cls, arr: np.ndarray) -> Self:
        """
        Create a StateVector instance from a 1D NumPy array.
        This is typically used when receiving the state from the ODE solver.

        Args:
            arr: 1D NumPy array of length 17

        Returns:
            StateVector instance
        """
        if len(arr) != 17:
            raise ValueError(f"Expected array of length 17, got {len(arr)}")

        return cls(
            X_ALG=arr[0],
            X_AOB=arr[1],
            X_NOB=arr[2],
            X_H=arr[3],
            X_S=arr[4],
            X_I=arr[5],
            S_S=arr[6],
            S_I=arr[7],
            S_IC=arr[8],
            S_ND=arr[9],
            S_NH=arr[10],
            S_NO2=arr[11],
            S_NO3=arr[12],
            S_N2=arr[13],
            S_PO4=arr[14],
            S_O2=arr[15],
            S_H2O=arr[16],
        )
