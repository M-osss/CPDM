import pyomo.environ as pyo
from pyomo.environ import units as pyunits

# =============================================================================
# INPUT PARAMETERS — edit this section only (steam/flue gas boiler)
# =============================================================================

# --- Cold side (tube: steam/water) inlet ---
COLD_FLOW_MOL       = 24678.26   # mol/s
COLD_TEMPERATURE_IN = 773.15     # K
COLD_PRESSURE_IN    = 2.5449e7   # Pa

# --- Hot side (shell: flue gas) inlet ---
FG_TOTAL_FLOW_MOL   = 28.3876e3 * 0.18  # mol/s
FG_TEMPERATURE_IN   = 1102.335   # K
FG_PRESSURE_IN      = 100145     # Pa
FG_COMPOSITION      = {"H2O": 0.0869, "CO2": 0.1449, "N2": 0.7434, "O2": 0.0247, "NO": 0.0006, "SO2": 0.002}

# --- Geometry (inch-to-meter: ITM = 0.0254) ---
TUBE_DI             = (2.5 - 2 * 0.165) * 0.0254   # m
TUBE_THICKNESS      = 0.165 * 0.0254                # m
PITCH_X             = 3 * 0.0254                    # m
PITCH_Y             = 54.78 / 108 * 12 * 0.0254    # m
TUBE_LENGTH         = 53.13 * 12 * 0.0254          # m
TUBE_NROW           = 40
TUBE_NCOL           = 108
NROW_INLET          = 4
DELTA_ELEVATION     = 50
EMISSIVITY_WALL     = 0.7
FCORRECTION_HTC     = 1.5
TUBE_R_FOULING      = 0.003131

# =============================================================================
# END INPUT PARAMETERS
# =============================================================================

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver

from idaes.models.properties import iapws95
from idaes.models_extra.power_generation.properties import FlueGasParameterBlock

from idaes.models_extra.power_generation.unit_models.boiler_heat_exchanger import (
    BoilerHeatExchanger,
    TubeArrangement,
    HeatExchangerFlowPattern,
    delta_temperature_underwood_callback,
)


def main():
    """
    BoilerHeatExchanger (0D) model in IDAES (power generation library).

    Property provision note:
    - Steam/water properties: IDAES `iapws95` (built-in, no external database needed).
    - Flue gas properties: IDAES `FlueGasParameterBlock` (built-in).
    """
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.prop_steam = iapws95.Iapws95ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()

    m.fs.bhx = BoilerHeatExchanger(
        delta_temperature_callback=delta_temperature_underwood_callback,
        cold_side={"property_package": m.fs.prop_steam, "has_pressure_change": True},
        hot_side={"property_package": m.fs.prop_fluegas, "has_pressure_change": True},
        has_holdup=False,
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
        tube_arrangement=TubeArrangement.inLine,
        cold_side_water_phase="Liq",
        has_radiation=True,
    )

    t0 = m.fs.time.first()

    h = pyo.value(iapws95.htpx(COLD_TEMPERATURE_IN * pyunits.K, COLD_PRESSURE_IN * pyunits.Pa))
    m.fs.bhx.cold_side_inlet.flow_mol[t0].fix(COLD_FLOW_MOL)
    m.fs.bhx.cold_side_inlet.enth_mol[t0].fix(h)
    m.fs.bhx.cold_side_inlet.pressure[t0].fix(COLD_PRESSURE_IN)

    for comp, frac in FG_COMPOSITION.items():
        m.fs.bhx.hot_side_inlet.flow_mol_comp[t0, comp].fix(FG_TOTAL_FLOW_MOL * frac)
    m.fs.bhx.hot_side_inlet.temperature[t0].fix(FG_TEMPERATURE_IN)
    m.fs.bhx.hot_side_inlet.pressure[t0].fix(FG_PRESSURE_IN)

    m.fs.bhx.tube_di.fix(TUBE_DI)
    m.fs.bhx.tube_thickness.fix(TUBE_THICKNESS)
    m.fs.bhx.pitch_x.fix(PITCH_X)
    m.fs.bhx.pitch_y.fix(PITCH_Y)
    m.fs.bhx.tube_length.fix(TUBE_LENGTH)
    m.fs.bhx.tube_nrow.fix(TUBE_NROW)
    m.fs.bhx.tube_ncol.fix(TUBE_NCOL)
    m.fs.bhx.nrow_inlet.fix(NROW_INLET)
    m.fs.bhx.delta_elevation.fix(DELTA_ELEVATION)
    m.fs.bhx.tube_r_fouling = TUBE_R_FOULING
    m.fs.bhx.emissivity_wall.fix(EMISSIVITY_WALL)
    m.fs.bhx.fcorrection_htc.fix(FCORRECTION_HTC)
    m.fs.bhx.fcorrection_dp_tube.fix(1.0)
    m.fs.bhx.fcorrection_dp_shell.fix(1.0)

    # Initialize and solve
    m.fs.bhx.initialize()
    solver = get_solver()
    res = solver.solve(m, tee=False)

    print("\n--- IDAES BoilerHeatExchanger (0D, power-generation) ---")
    print(f"Solver termination: {res.solver.termination_condition}")
    print(f"Cold out T: {pyo.value(m.fs.bhx.cold_side.properties_out[t0].temperature):.2f} K")
    print(f"Hot out T: {pyo.value(m.fs.bhx.hot_side.properties_out[t0].temperature):.2f} K")
    print(f"Heat duty: {pyo.value(m.fs.bhx.heat_duty[t0]) / 1e6:.3f} MW")
    print(f"Area: {pyo.value(m.fs.bhx.area):.2f} m^2")
    print(f"U: {pyo.value(m.fs.bhx.overall_heat_transfer_coefficient[t0]):.2f} W/m^2/K")


if __name__ == "__main__":
    main()

