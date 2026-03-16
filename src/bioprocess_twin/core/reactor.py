from pathlib import Path
from typing import Self

import yaml
from pydantic import BaseModel, ConfigDict, Field


class GeometryConfig(BaseModel):
    """Physical dimensions of the reactor."""
    model_config = ConfigDict(frozen=True)

    surface_area: float = Field(..., gt=0, description="Total surface area in m2")
    length: float = Field(..., gt=0, description="Total length in m")
    nominal_depth: float = Field(..., gt=0, description="Water depth in m")
    total_volume: float = Field(..., gt=0, description="Total volume in m3")


class OperationalConfig(BaseModel):
    """Operational parameters."""
    model_config = ConfigDict(frozen=True)

    linear_velocity: float = Field(..., ge=0, description="Velocity in m/s")
    pump_flow_rate: float = Field(..., ge=0, description="Recirculation flow in m3/d")
    kla_20: float = Field(..., ge=0, description="Mass transfer coefficient at 20C in 1/d")


class LocationConfig(BaseModel):
    """Geographic metadata."""
    model_config = ConfigDict(frozen=True)

    name: str
    latitude: float
    longitude: float


class ReactorConfig(BaseModel):
    """
    Main configuration for the reactor.
    Aggregates geometry, operation, and location.
    """
    model_config = ConfigDict(frozen=True)

    name: str
    geometry: GeometryConfig
    operational: OperationalConfig
    location: LocationConfig | None = None

    @classmethod
    def from_yaml(cls, path: str | Path) -> Self:
        """Load configuration from a YAML file."""
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Config file not found: {path}")

        with open(path, "r") as f:
            data = yaml.safe_load(f)
        
        # Expects a 'reactor' root key
        if "reactor" not in data:
            raise ValueError("YAML must contain a 'reactor' root key")
            
        return cls(**data["reactor"])

    def calculate_hrt(self, influent_flow: float) -> float:
        """
        Calculate Hydraulic Retention Time (HRT).
        HRT = Volume / Q_in
        
        Args:
            influent_flow: Influent flow rate in m3/d
            
        Returns:
            HRT in days
        """
        if influent_flow <= 0:
            raise ValueError("Influent flow must be positive")
        return self.geometry.total_volume / influent_flow

    @property
    def theoretical_volume(self) -> float:
        """Calculate volume based on surface area and depth (A * h)."""
        return self.geometry.surface_area * self.geometry.nominal_depth

    @property
    def volume_discrepancy(self) -> float:
        """
        Calculate the difference between reported total volume and theoretical volume.
        Returns: Absolute difference in m3.
        """
        return abs(self.geometry.total_volume - self.theoretical_volume)
