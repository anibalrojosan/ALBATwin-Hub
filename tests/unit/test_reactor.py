from pathlib import Path

import pytest
from pydantic import ValidationError

from bioprocess_twin.core.reactor import ReactorConfig


def test_reactor_config_loading():
    """
    Test that ReactorConfig correctly loads data from the actual YAML file.
    This ensures the link between the config file and the Python object is solid.
    """
    config_path = Path("config/reactor_setup.yaml")

    # Load the config
    config = ReactorConfig.from_yaml(config_path)

    # Verify exact values from the Narbonne HRABP paper
    assert config.name == "Narbonne HRABP Raceway"
    assert config.geometry.surface_area == 56.0
    assert config.geometry.nominal_depth == 0.3
    assert config.geometry.total_volume == 17.0
    assert config.operational.kla_20 == 34.0


def test_reactor_geometry_validation():
    """
    Test that Pydantic validation works.
    It should raise a ValidationError if we try to create a reactor with invalid data.
    """
    # Invalid data (negative surface area)
    invalid_data = {
        "name": "Invalid Reactor",
        "geometry": {"surface_area": -10.0, "nominal_depth": 0.3, "total_volume": 17.0},
        "operational": {"kla_20": 34.0, "linear_velocity": 0.2, "pump_flow_rate": 182.0},
    }

    with pytest.raises(ValidationError):
        ReactorConfig(**invalid_data)


def test_reactor_config_missing_fields():
    """Test that the ReactorConfig raises a ValidationError if the YAML is incomplete."""
    incomplete_data = {
        "name": "Incomplete Reactor",
        "geometry": {"surface_area": 56.0},  # Missing depth and volume
    }
    with pytest.raises(ValidationError):
        ReactorConfig(**incomplete_data)


def test_hrt_calculation():
    """Test the calculation of the Hydraulic Retention Time."""
    config_path = Path("config/reactor_setup.yaml")
    config = ReactorConfig.from_yaml(config_path)

    influent_flow = 4.25  # m3/d
    expected_hrt = 17.0 / 4.25
    assert config.calculate_hrt(influent_flow) == expected_hrt


def test_from_yaml_file_not_found():
    """Test that the ReactorConfig raises a FileNotFoundError if the file does not exist."""
    with pytest.raises(FileNotFoundError):
        ReactorConfig.from_yaml("non_existent_file.yaml")


def test_config_is_immutable():
    """Test that the ReactorConfig is immutable once loaded (frozen=True)."""
    config_path = Path("config/reactor_setup.yaml")
    config = ReactorConfig.from_yaml(config_path)

    # Attempting to modify a value should raise a ValidationError or AttributeError
    with pytest.raises((ValidationError, AttributeError)):
        config.geometry.surface_area = 100.0


def test_volume_calculation_logic():
    """
    Test derived properties for volume consistency.
    The class should expose theoretical volume and discrepancy check.
    """
    config_path = Path("config/reactor_setup.yaml")
    config = ReactorConfig.from_yaml(config_path)

    # Theoretical volume: 56m2 * 0.3m = 16.8 m3
    assert abs(config.theoretical_volume - 16.8) < 1e-6

    # Reported volume is 17.0 m3, so discrepancy should be 0.2 m3
    assert abs(config.volume_discrepancy - 0.2) < 1e-6
