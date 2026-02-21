"""Heat exchanger component and geometry config. Properties from pcd/ipd/thermo.

System structure:
- system["components"]: list of all components (defines property package).
- stream["composition"]: mole fractions by CAS; all system components must have a fraction.
- stream flow_mol, temperature, pressure: edit these in the system preset below; they are
  used directly (no override from INPUT_PARAMETERS). INPUT_PARAMETERS overrides geometry only.
"""

# =============================================================================
# INPUT PARAMETERS — geometry overrides only. Stream data (flow, T, P, composition)
# is taken from the system presets below (EASY_SYSTEM, HYDROCARBON_SYSTEM, etc.).
# Edit the system's hot_stream/cold_stream directly to change flowrates and conditions.
# =============================================================================

INPUT_PARAMETERS = {
    "geometry": {
        "0D": {
            "area": 50.0,
            "overall_U": 600.0,
            "hot_side_deltaP": -8000,
            "cold_side_deltaP": -12000,
            "hot_outlet_temperature": None,
        },
        "1D": {
            "area": 50.0,
            "length": 5.0,
            "heat_transfer_coefficient": 600.0,
            "finite_elements": 10,
            "collocation_points": 5,
        },
        "ShellAndTube1D": {
            "length": 5.0,
            "shell_diameter": 0.5,
            "tube_inner_diameter": 0.020,
            "tube_outer_diameter": 0.025,
            "number_of_tubes": 100,
            "hot_side_htc": 1000.0,
            "cold_side_htc": 1200.0,
            "finite_elements": 10,
            "collocation_points": 5,
        },
    },
}

# =============================================================================
# END INPUT PARAMETERS
# =============================================================================


def _deep_merge(base, override):
    """Merge override into base; return new dict. None in override skips that key."""
    if override is None:
        return dict(base)
    result = dict(base)
    for k, v in override.items():
        if v is None:
            continue
        if k in result and isinstance(result[k], dict) and isinstance(v, dict):
            result[k] = _deep_merge(result[k], v)
        else:
            result[k] = v
    return result


EASY_SYSTEM = {
    "name": "Water-Toluene-Nitrogen Heat Exchange",
    "components": [
        {"identifier": "7732-18-5", "name": "water"},
        {"identifier": "108-88-3", "name": "toluene"},
        {"identifier": "7727-37-9", "name": "nitrogen"},
    ],
    "hot_stream": {
        "flow_mol": 100.0,
        "temperature": 373.15,
        "pressure": 501325,
        "composition": {"7732-18-5": 0.9, "108-88-3": 0.001, "7727-37-9": 0.099},
    },
    "cold_stream": {
        "flow_mol": 100.0,
        "temperature": 323.15,
        "pressure": 202650,
        "composition": {"7732-18-5": 0.001, "108-88-3": 0.999},
    },
}

WATER_TOLUENE_LIQUID = {
    "name": "Water-Toluene Liquid (Liquid Phase Only)",
    "components": [
        {"identifier": "7732-18-5", "name": "water"},
        {"identifier": "108-88-3", "name": "toluene"},
    ],
    "hot_stream": {
        "flow_mol": 100.0,
        "temperature": 573.15,
        "pressure": 101325,
        "composition": {"7732-18-5": 0.999, "108-88-3": 0.001},
    },
    "cold_stream": {
        "flow_mol": 100.0,
        "temperature": 323.15,
        "pressure": 202650,
        "composition": {"7732-18-5": 0.001, "108-88-3": 0.999},
    },
}

HYDROCARBON_SYSTEM = {
    "name": "Hydrocarbon Heat Exchange",
    "components": [
        {"identifier": "74-82-8", "name": "methane"},
        {"identifier": "74-84-0", "name": "ethane"},
        {"identifier": "74-98-6", "name": "propane"},
        {"identifier": "106-97-8", "name": "butane"},
    ],
    "hot_stream": {
        "flow_mol": 50.0,
        "temperature": 700.0,
        "pressure": 500000,
        "composition": {
            "74-82-8": 0.01,
            "74-84-0": 0.49,
            "74-98-6": 0.49,
            "106-97-8": 0.01,
        },
    },
    "cold_stream": {
        "flow_mol": 200.0,
        "temperature": 280.0,
        "pressure": 600000,
        "composition": {
            "74-82-8": 0.70,
            "74-84-0": 0.29,
            "74-98-6": 0.005,
            "106-97-8": 0.005,
        },
    },
}

CO2_CAPTURE_SYSTEM = {
    "name": "CO2-MEA System",
    "components": [
        {"identifier": "124-38-9", "name": "co2"},
        {"identifier": "7732-18-5", "name": "water"},  # Water
        {"identifier": "141-43-5", "name": "mea"},  # Monoethanolamine
    ],
    "hot_stream": {
        "flow_mol": 60.55,
        "temperature": 392.23,
        "pressure": 202650,
        "composition": {"124-38-9": 0.0158, "7732-18-5": 0.8747, "141-43-5": 0.1095},
    },
    "cold_stream": {
        "flow_mol": 63.02,
        "temperature": 326.36,
        "pressure": 202650,
        "composition": {"124-38-9": 0.0414, "7732-18-5": 0.8509, "141-43-5": 0.1077},
    },
}

POWER_PLANT_SYSTEM = {
    "name": "Steam-Flue Gas (Power Plant)",
    "note": "Uses IDAES built-in iapws95 and FlueGasParameterBlock",
    "cold_stream": {
        "property_package": "iapws95",
        "flow_mol": 24678.26,
        "temperature": 773.15,
        "pressure": 2.5449e7,
    },
    "hot_stream": {
        "property_package": "FlueGasParameterBlock",
        "total_flow_mol": 28.3876e3 * 0.18,
        "temperature": 1102.335,  # K
        "pressure": 100145,  # Pa
        "composition": {
            "H2O": 0.0869,
            "CO2": 0.1449,
            "N2": 0.7434,
            "O2": 0.0247,
            "NO": 0.0006,
            "SO2": 0.002,
        },
    },
}

DEFAULT_SYSTEM = EASY_SYSTEM
#geometry defaults are advanced settings, make sure you know what you are doing before changing them.
GEOMETRY_DEFAULTS = {
    "0D": {
        "area": 50.0,
        "overall_U": 600.0,
        "hot_side_deltaP": -8000,
        "cold_side_deltaP": -12000,
        "hot_outlet_temperature": None,
    },
    "1D": {
        "area": 25.0,
        "length": 5.0,
        "heat_transfer_coefficient": 600.0,
        "finite_elements": 10,
        "collocation_points": 5,
    },
    "ShellAndTube1D": {
        "length": 5.0,
        "shell_diameter": 0.5,
        "tube_inner_diameter": 0.020,
        "tube_outer_diameter": 0.025,
        "number_of_tubes": 100,
        "tube_pitch": 0.035,
        "tube_layout": "triangular",
        "hot_side_htc": 1000.0,
        "cold_side_htc": 1200.0,
        "finite_elements": 10,
        "collocation_points": 5,
    },
    "Plate": {
        "area": 114.3,
        "plate_length": 1.270,
        "plate_width": 0.508,
        "plate_thickness": 0.0006,
        "plate_pact_length": 0.381,
        "port_diameter": 0.2045,
        "plate_thermal_conductivity": 16.2,
        "passes": 4,
        "channels_per_pass": 12,
        "number_of_divider_plates": 2,
    },
    "Boiler": {
        "tube_di": 0.0587,
        "tube_thickness": 0.00419,
        "pitch_x": 0.0762,
        "pitch_y": 0.1524,
        "tube_length": 16.16,
        "tube_nrow": 40,
        "tube_ncol": 108,
        "nrow_inlet": 4,
        "delta_elevation": 50,
        "tube_r_fouling": 0.003131,
        "emissivity_wall": 0.7,
        "fcorrection_htc": 1.5,
        "fcorrection_dp_tube": 1.0,
        "fcorrection_dp_shell": 1.0,
    },
    "Boiler2D": {
        "tube_inner_diameter": 0.0587,
        "tube_thickness": 0.00419,
        "header_inner_diameter": 0.30,
        "header_wall_thickness": 0.02,
        "tube_ncol": 108,
        "tube_nseg": 10,
        "tube_inlet_nrow": 4,
        "pitch_x": 0.0762,
        "pitch_y": 0.1524,
        "tube_length_seg": 1.616,
        "delta_elevation": 50,
        "finite_elements": 10,
        "collocation_points": 3,
        "radial_elements": 8,
        "header_radial_elements": 6,
        "tube_r_fouling": 0.003131,
        "shell_r_fouling": 0.0001,
        "emissivity_wall": 0.7,
        "fcorrection_htc_shell": 1.5,
        "fcorrection_htc_tube": 1.5,
        "fcorrection_dp_tube": 1.0,
        "fcorrection_dp_shell": 1.0,
    },
}


def _filter_composition(comp_dict, component_ids):
    """Keep only composition keys that exist in component_ids; renormalize to sum=1."""
    if comp_dict is None:
        return comp_dict
    ids_set = {c["identifier"] for c in component_ids}
    filtered = {k: v for k, v in comp_dict.items() if k in ids_set}
    total = sum(filtered.values())
    if total > 1e-12:
        return {k: v / total for k, v in filtered.items()}
    return filtered


def get_system(system_name: str = None):
    """Return system for database-driven models (lmtd, ntu, lc, 1d, shell_tube).
    power_plant is not supported; use idaes_hex_boiler or idaes_hex_boiler_2d instead."""
    systems = {
        "water_toluene_nitrogen": EASY_SYSTEM,
        "water_toluene": WATER_TOLUENE_LIQUID,
        "liquid": WATER_TOLUENE_LIQUID,
        "hydrocarbon": HYDROCARBON_SYSTEM,
        "co2_capture": CO2_CAPTURE_SYSTEM,
    }
    key = (system_name or "").strip().lower()
    if key == "power_plant":
        raise ValueError(
            "power_plant uses iapws95+FlueGas, not generic properties. "
            "Use idaes_hex_boiler.py or idaes_hex_boiler_2d.py instead."
        )
    base = systems.get(key, DEFAULT_SYSTEM)
    if "components" not in base:
        return base
    hot = dict(base["hot_stream"])
    cold = dict(base["cold_stream"])
    comp_ids = base["components"]
    if "composition" in hot:
        hot["composition"] = _filter_composition(hot["composition"], comp_ids)
    if "composition" in cold:
        cold["composition"] = _filter_composition(cold["composition"], comp_ids)
    return {
        "name": base["name"],
        "components": base["components"],
        "hot_stream": hot,
        "cold_stream": cold,
    }


def get_geometry(model_type: str):
    base = GEOMETRY_DEFAULTS.get(model_type, GEOMETRY_DEFAULTS["0D"])
    geom_overrides = INPUT_PARAMETERS.get("geometry", {}).get(model_type)
    return _deep_merge(base, geom_overrides) if geom_overrides else dict(base)


def validate_system(system: dict) -> bool:
    required_keys = ["components", "hot_stream", "cold_stream"]
    for key in required_keys:
        if key not in system:
            raise ValueError(f"System configuration missing required key: {key}")
    
    if not system["components"]:
        raise ValueError("System must have at least one component")
    
    for comp in system["components"]:
        if "identifier" not in comp:
            raise ValueError(f"Component missing 'identifier' field: {comp}")
    
    return True
