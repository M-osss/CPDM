import pyomo.environ as pyo
from pyomo.environ import units as pyunits

# =============================================================================
# INPUT PARAMETERS — edit this section only (2D crossflow boiler)
# =============================================================================

ITM = 0.0254

# --- Tube-side (steam/water) inlet ---
TUBE_FLOW_MOL       = 24678.26   # mol/s
TUBE_TEMPERATURE_IN = 773.15     # K
TUBE_PRESSURE_IN    = 2.5449e7   # Pa

# --- Shell-side (flue gas) inlet ---
SHELL_TOTAL_FLOW   = 28.3876e3 * 0.18  # mol/s
SHELL_TEMPERATURE   = 1102.335   # K
SHELL_PRESSURE     = 100145     # Pa
SHELL_COMPOSITION  = {"H2O": 0.0869, "CO2": 0.1449, "N2": 0.7434, "O2": 0.0247, "NO": 0.0006, "SO2": 0.002}

# --- Geometry ---
TUBE_DI            = (2.5 - 2 * 0.165) * ITM
TUBE_THICKNESS      = 0.165 * ITM
HEADER_DI          = 0.30
HEADER_THICKNESS   = 0.02
TUBE_NCOL          = 108
TUBE_NSEG          = 10
TUBE_INLET_NROW    = 4
PITCH_X            = 3 * ITM
PITCH_Y            = 54.78 / 108 * 12 * ITM
TUBE_LENGTH_SEG    = (53.13 * 12 * ITM) / 10.0
DELTA_ELEVATION    = 50
FINITE_ELEMENTS    = 10
COLLOCATION_POINTS = 3
RADIAL_ELEMENTS    = 8
HEADER_RADIAL_EL   = 6
EMISSIVITY_WALL    = 0.7
TUBE_R_FOULING     = 0.003131
SHELL_R_FOULING   = 0.0001
FCORRECTION_HTC_SHELL = 1.5
FCORRECTION_HTC_TUBE  = 1.5

# =============================================================================
# END INPUT PARAMETERS
# =============================================================================

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver

import idaes.logger as idaeslog

from idaes.models.properties import iapws95
from idaes.models_extra.power_generation.properties import FlueGasParameterBlock
from idaes.models_extra.power_generation.unit_models.boiler_heat_exchanger_2D import (
    HeatExchangerCrossFlow2D_Header,
)


def main(system_override: dict = None):
    """
    2D crossflow boiler heat exchanger with tube-wall radial discretization.

    Model: `HeatExchangerCrossFlow2D_Header` (IDAES power_generation library).

    Args:
        system_override: Optional system dict with custom steam/flue-gas conditions

    Property provision note:
    - Tube-side (water/steam): IDAES `iapws95` (built-in).
    - Shell-side (flue gas): IDAES `FlueGasParameterBlock` (built-in).
    - Tube-wall material properties are internal Params (therm_cond_wall, cp_wall, dens_wall, etc.).
    - For custom components beyond steam/flue-gas, users must provide
      compatible property packages.
    """
    print("\n" + "="*80)
    print("IDAES HeatExchangerCrossFlow2D_Header (Power Plant System, 2D)")
    print("Property Packages: iapws95 (steam) + FlueGas (built-in)")
    print("="*80 + "\n")
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.prop_steam = iapws95.Iapws95ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()

    m.fs.hx2d = HeatExchangerCrossFlow2D_Header(
        shell_side={"property_package": m.fs.prop_fluegas, "has_pressure_change": True},
        tube_side={"property_package": m.fs.prop_steam, "has_pressure_change": True},
        flow_type="counter_current",
        tube_arrangement="in-line",
        tube_side_water_phase="Liq",
        has_radiation=True,
        finite_elements=FINITE_ELEMENTS,
        collocation_points=COLLOCATION_POINTS,
        radial_elements=RADIAL_ELEMENTS,
        tube_inner_diameter=TUBE_DI,
        tube_thickness=TUBE_THICKNESS,
        has_header=True,
        header_inner_diameter=HEADER_DI,
        header_wall_thickness=HEADER_THICKNESS,
        header_radial_elements=HEADER_RADIAL_EL,
    )

    t0 = m.fs.time.first()

    h = pyo.value(iapws95.htpx(TUBE_TEMPERATURE_IN * pyunits.K, TUBE_PRESSURE_IN * pyunits.Pa))
    m.fs.hx2d.tube_inlet.flow_mol[t0].fix(TUBE_FLOW_MOL)
    m.fs.hx2d.tube_inlet.enth_mol[t0].fix(h)
    m.fs.hx2d.tube_inlet.pressure[t0].fix(TUBE_PRESSURE_IN)

    for comp, frac in SHELL_COMPOSITION.items():
        m.fs.hx2d.shell_inlet.flow_mol_comp[t0, comp].fix(SHELL_TOTAL_FLOW * frac)
    m.fs.hx2d.shell_inlet.temperature[t0].fix(SHELL_TEMPERATURE)
    m.fs.hx2d.shell_inlet.pressure[t0].fix(SHELL_PRESSURE)

    m.fs.hx2d.tube_ncol.fix(TUBE_NCOL)
    m.fs.hx2d.tube_nseg.fix(TUBE_NSEG)
    m.fs.hx2d.tube_inlet_nrow.fix(TUBE_INLET_NROW)
    m.fs.hx2d.pitch_x.fix(PITCH_X)
    m.fs.hx2d.pitch_y.fix(PITCH_Y)
    m.fs.hx2d.tube_length_seg.fix(TUBE_LENGTH_SEG)
    m.fs.hx2d.delta_elevation.fix(DELTA_ELEVATION)

    m.fs.hx2d.fcorrection_htc_shell.fix(FCORRECTION_HTC_SHELL)
    m.fs.hx2d.fcorrection_htc_tube.fix(FCORRECTION_HTC_TUBE)
    m.fs.hx2d.fcorrection_dp_tube.fix(1.0)
    m.fs.hx2d.fcorrection_dp_shell.fix(1.0)
    m.fs.hx2d.tube_r_fouling.set_value(TUBE_R_FOULING)
    m.fs.hx2d.shell_r_fouling.set_value(SHELL_R_FOULING)
    m.fs.hx2d.emissivity_wall.fix(EMISSIVITY_WALL)

    # Initialize and solve
    m.fs.hx2d.initialize(outlvl=idaeslog.INFO)  # type: ignore[name-defined]
    solver = get_solver()
    res = solver.solve(m, tee=False)

    z_shell_out = m.fs.hx2d.shell.length_domain.last()
    # Tube length domain direction depends on flow_type (counter_current uses backward),
    # so the outlet is at the *first* point when tube is backward.
    z_tube_out = (
        m.fs.hx2d.tube.length_domain.first()
        if m.fs.hx2d.config.flow_type == "counter_current"
        else m.fs.hx2d.tube.length_domain.last()
    )

    print("\n--- IDAES HeatExchangerCrossFlow2D_Header (boiler 2D) ---")
    print(f"Solver termination: {res.solver.termination_condition}")
    print(f"Tube out T: {pyo.value(m.fs.hx2d.tube.properties[t0, z_tube_out].temperature):.2f} K")
    print(f"Shell out T: {pyo.value(m.fs.hx2d.shell.properties[t0, z_shell_out].temperature):.2f} K")
    print(f"Total heat (from shell): {pyo.value(m.fs.hx2d.total_heat[t0]) / 1e6:.3f} MW")


if __name__ == "__main__":
    main()

