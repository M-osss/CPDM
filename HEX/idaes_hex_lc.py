import pyomo.environ as pyo
from pyomo.environ import units as pyunits

# =============================================================================
# INPUT PARAMETERS — edit to override defaults from hex_component_properties.
# Stream/geometry from system presets in hex_component_properties.
# LC-specific parameters below (or set geometry_overrides in main()).
# =============================================================================

LC_AREA              = 70.0      # m^2
LC_UA_HOT_SIDE       = 2.0e4     # W/K
LC_UA_COLD_SIDE      = 2.0e4     # W/K
LC_HEAT_CAPACITY_WALL = 5.0e5    # J/K
LC_THERMAL_FOULING_HOT = 0.0
LC_THERMAL_FOULING_COLD = 0.0

# =============================================================================
# END INPUT PARAMETERS
# =============================================================================

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

from idaes.models.unit_models import HeatExchangerLumpedCapacitance
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

import idaes_property_builder as ipb
import hex_component_properties as hcp


def main(system_name: str = "liquid", geometry_overrides: dict = None):
    system = hcp.get_system(system_name)
    hcp.validate_system(system)
    
    geom = hcp.get_geometry("0D")
    if geometry_overrides:
        geom.update(geometry_overrides)
    geom.setdefault("area", LC_AREA)
    
    print(f"\n{'='*80}")
    print(f"IDAES HeatExchangerLC: {system['name']}")
    print(f"Components: {', '.join(c['name'] for c in system['components'])}")
    print(f"{'='*80}\n")
    
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 10], time_units=pyunits.s)

    components_list = [c["identifier"] for c in system["components"]]
    props_cfg, built = ipb.build_generic_parameter_block_config(
        components_list, phases=("Liq",), eos="PR", prefer_rpp_ig_cp=True
    )
    ipb.validate_config_dict(props_cfg)
    m.fs.props = GenericParameterBlock(**props_cfg)

    m.fs.hx = HeatExchangerLumpedCapacitance(
        hot_side={"property_package": m.fs.props, "has_pressure_change": True},
        cold_side={"property_package": m.fs.props, "has_pressure_change": True},
        dynamic=False,
        dynamic_heat_balance=False,
    )

    comp_id_map = {}
    for cas, built_comp in zip(components_list, built):
        comp_id_map[cas] = built_comp.comp_id

    hot = system["hot_stream"]
    cold = system["cold_stream"]
    
    for t in m.fs.time:
        m.fs.hx.hot_side.properties_in[t].flow_mol.fix(hot["flow_mol"])
        m.fs.hx.hot_side.properties_in[t].pressure.fix(hot["pressure"])
        m.fs.hx.hot_side.properties_in[t].temperature.fix(hot["temperature"])
        for cas, mole_frac in hot["composition"].items():
            m.fs.hx.hot_side.properties_in[t].mole_frac_comp[comp_id_map[cas]].fix(mole_frac)

        m.fs.hx.cold_side.properties_in[t].flow_mol.fix(cold["flow_mol"])
        m.fs.hx.cold_side.properties_in[t].pressure.fix(cold["pressure"])
        m.fs.hx.cold_side.properties_in[t].temperature.fix(cold["temperature"])
        for cas, mole_frac in cold["composition"].items():
            m.fs.hx.cold_side.properties_in[t].mole_frac_comp[comp_id_map[cas]].fix(mole_frac)

        m.fs.hx.hot_side.deltaP[t].fix(geom["hot_side_deltaP"])
        m.fs.hx.cold_side.deltaP[t].fix(geom["cold_side_deltaP"])

    m.fs.hx.area.fix(geom["area"])
    m.fs.hx.heat_capacity_wall.set_value(LC_HEAT_CAPACITY_WALL)
    m.fs.hx.thermal_resistance_wall.set_value(0.0)
    m.fs.hx.thermal_fouling_hot_side.set_value(LC_THERMAL_FOULING_HOT)
    m.fs.hx.thermal_fouling_cold_side.set_value(LC_THERMAL_FOULING_COLD)

    for t in m.fs.time:
        m.fs.hx.ua_hot_side[t].fix(LC_UA_HOT_SIDE)
        m.fs.hx.ua_cold_side[t].fix(LC_UA_COLD_SIDE)

    m.fs.hx.initialize(outlvl=idaeslog.INFO)
    solver = get_solver()
    res = solver.solve(m, tee=False)

    print("\n--- IDAES HeatExchangerLumpedCapacitance ---")
    print(f"Solver termination: {res.solver.termination_condition}")
    t0 = m.fs.time.first()
    print(f"Area: {pyo.value(m.fs.hx.area):.3f} m^2")
    print(f"UA hot/cold: {pyo.value(m.fs.hx.ua_hot_side[t0]):.2f} / {pyo.value(m.fs.hx.ua_cold_side[t0]):.2f} W/K")
    print(f"Wall T: {pyo.value(m.fs.hx.temperature_wall[t0]):.2f} K")
    print(f"Heat duty: {pyo.value(m.fs.hx.heat_duty[t0]) / 1000:.3f} kW")
    print(
        f"Hot out T: {pyo.value(m.fs.hx.hot_side.properties_out[t0].temperature):.2f} K | "
        f"Cold out T: {pyo.value(m.fs.hx.cold_side.properties_out[t0].temperature):.2f} K"
    )


if __name__ == "__main__":
    import sys
    
    system_name = sys.argv[1] if len(sys.argv) > 1 else "liquid"
    main(system_name=system_name)

