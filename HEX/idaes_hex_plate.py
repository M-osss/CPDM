import pyomo.environ as pyo
from pyomo.environ import units as pyunits

# =============================================================================
# INPUT PARAMETERS — edit this section only (MEA/CO2 capture PHE)
# =============================================================================

# --- Hot side inlet ---
HOT_FLOW_MOL       = 60.54879   # mol/s
HOT_TEMPERATURE_IN = 392.23     # K
HOT_PRESSURE_IN    = 202650     # Pa
HOT_MOLE_FRACS     = {"CO2": 0.0158, "H2O": 0.8747, "MEA": 0.1095}

# --- Cold side inlet ---
COLD_FLOW_MOL       = 63.01910  # mol/s
COLD_TEMPERATURE_IN = 326.36    # K
COLD_PRESSURE_IN    = 202650    # Pa
COLD_MOLE_FRACS     = {"CO2": 0.0414, "H2O": 0.8509, "MEA": 0.1077}

# --- Geometry (uses IDAES defaults if not overridden) ---
INIT_DUTY_GUESS = 245000  # W (for initialization)

# =============================================================================
# END INPUT PARAMETERS
# =============================================================================

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.column_models.plate_heat_exchanger import PlateHeatExchanger
from idaes.models_extra.column_models.properties.MEA_solvent import (
    configuration as aqueous_mea,
)


def main():
    """
    Plate heat exchanger (PHE) model in IDAES.

    Property provision note:
    - This script uses IDAES-provided MEA solvent property configuration
      (`idaes.models_extra.column_models.properties.MEA_solvent`), which includes
      viscosity and thermal conductivity required by PHE correlations.
    - No ChemSep/IPD extraction is used here because the PHE model expects
      transport properties in the IDAES modular property framework.
    """
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.hotside_properties = GenericParameterBlock(**aqueous_mea)
    m.fs.coldside_properties = GenericParameterBlock(**aqueous_mea)

    m.fs.phe = PlateHeatExchanger(
        passes=4,
        channels_per_pass=12,
        number_of_divider_plates=2,
        hot_side={"property_package": m.fs.hotside_properties},
        cold_side={"property_package": m.fs.coldside_properties},
    )

    t0 = m.fs.time.first()

    m.fs.phe.hot_side_inlet.flow_mol[t0].fix(HOT_FLOW_MOL)
    m.fs.phe.hot_side_inlet.temperature[t0].fix(HOT_TEMPERATURE_IN)
    m.fs.phe.hot_side_inlet.pressure[t0].fix(HOT_PRESSURE_IN)
    for comp, xf in HOT_MOLE_FRACS.items():
        m.fs.phe.hot_side_inlet.mole_frac_comp[t0, comp].fix(xf)

    m.fs.phe.cold_side_inlet.flow_mol[t0].fix(COLD_FLOW_MOL)
    m.fs.phe.cold_side_inlet.temperature[t0].fix(COLD_TEMPERATURE_IN)
    m.fs.phe.cold_side_inlet.pressure[t0].fix(COLD_PRESSURE_IN)
    for comp, xf in COLD_MOLE_FRACS.items():
        m.fs.phe.cold_side_inlet.mole_frac_comp[t0, comp].fix(xf)

    # Fix unit geometry (IDAES defaults are reasonable; fixing makes DOF=0)
    m.fs.phe.plate_length.fix()
    m.fs.phe.plate_width.fix()
    m.fs.phe.plate_thickness.fix()
    m.fs.phe.plate_pact_length.fix()
    m.fs.phe.port_diameter.fix()
    m.fs.phe.plate_therm_cond.fix()
    m.fs.phe.area.fix()

    m.fs.phe.initialize(outlvl=idaeslog.INFO, duty=(INIT_DUTY_GUESS, pyunits.W))
    solver = get_solver()
    res = solver.solve(m, tee=False)

    print("\n--- IDAES PlateHeatExchanger (PHE) ---")
    print(f"Solver termination: {res.solver.termination_condition}")
    print(f"Heat duty: {pyo.value(m.fs.phe.heat_duty[t0]) / 1000:.3f} kW")
    print(f"Effectiveness: {pyo.value(m.fs.phe.effectiveness[t0]):.6f}")
    print(f"NTU: {pyo.value(m.fs.phe.NTU[t0]):.3f}")
    print(f"U (overall): {pyo.value(m.fs.phe.heat_transfer_coefficient[t0]):.2f} W/m^2/K")
    print(f"Re_hot: {pyo.value(m.fs.phe.Re_hot[t0]):.3f} | Re_cold: {pyo.value(m.fs.phe.Re_cold[t0]):.3f}")
    print(f"Pr_hot: {pyo.value(m.fs.phe.Pr_hot[t0]):.3f} | Pr_cold: {pyo.value(m.fs.phe.Pr_cold[t0]):.3f}")
    print(
        f"Hot out P: {pyo.value(m.fs.phe.hot_side_outlet.pressure[t0]):.1f} Pa | "
        f"Cold out P: {pyo.value(m.fs.phe.cold_side_outlet.pressure[t0]):.1f} Pa"
    )
    print(
        f"Hot out T: {pyo.value(m.fs.phe.hot_side_outlet.temperature[t0]):.3f} K | "
        f"Cold out T: {pyo.value(m.fs.phe.cold_side_outlet.temperature[t0]):.3f} K"
    )


if __name__ == "__main__":
    main()

