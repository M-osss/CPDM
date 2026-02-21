"""IDAES NTU heat exchanger. Edit system presets in hex_component_properties for flow/T/P/composition; INPUT_PARAMETERS for geometry only."""
#WARNING: This model only considers liquid phase, use others if you want to simulate vapor phase or VLE.
import pyomo.environ as pyo
from pyomo.environ import units as pyunits

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

from idaes.models.unit_models.heat_exchanger_ntu import HeatExchangerNTU
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

import idaes_property_builder as ipb
import hex_component_properties as hcp


def main(system_name: str = "water_toluene_nitrogen", geometry_overrides: dict = None):
    system = hcp.get_system(system_name)
    hcp.validate_system(system)
    
    geom = hcp.get_geometry("0D")
    if geometry_overrides:
        geom.update(geometry_overrides)
    
    print(f"\n{'='*80}")
    print(f"IDAES HeatExchangerNTU: {system['name']}")
    print(f"Components: {', '.join(c['name'] for c in system['components'])}")
    print(f"{'='*80}\n")
    
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    components_list = [c["identifier"] for c in system["components"]]
    props_cfg, built = ipb.build_generic_parameter_block_config(
        components_list, phases=("Liq",), eos="PR", prefer_rpp_ig_cp=True
    )
    ipb.validate_config_dict(props_cfg)
    m.fs.props = GenericParameterBlock(**props_cfg)

    m.fs.hx = HeatExchangerNTU(
        hot_side={"property_package": m.fs.props, "has_pressure_change": True},
        cold_side={"property_package": m.fs.props, "has_pressure_change": True},
    )

    t0 = m.fs.time.first()

    comp_id_map = {}
    for cas, built_comp in zip(components_list, built):
        comp_id_map[cas] = built_comp.comp_id
    
    hot = system["hot_stream"]
    m.fs.hx.hot_side.properties_in[t0].flow_mol.fix(hot["flow_mol"])
    m.fs.hx.hot_side.properties_in[t0].temperature.fix(hot["temperature"])
    m.fs.hx.hot_side.properties_in[t0].pressure.fix(hot["pressure"])
    for cas, mole_frac in hot["composition"].items():
        comp_id = comp_id_map[cas]
        m.fs.hx.hot_side.properties_in[t0].mole_frac_comp[comp_id].fix(mole_frac)

    cold = system["cold_stream"]
    m.fs.hx.cold_side.properties_in[t0].flow_mol.fix(cold["flow_mol"])
    m.fs.hx.cold_side.properties_in[t0].temperature.fix(cold["temperature"])
    m.fs.hx.cold_side.properties_in[t0].pressure.fix(cold["pressure"])
    for cas, mole_frac in cold["composition"].items():
        comp_id = comp_id_map[cas]
        m.fs.hx.cold_side.properties_in[t0].mole_frac_comp[comp_id].fix(mole_frac)

    m.fs.hx.hot_side.deltaP[t0].fix(geom["hot_side_deltaP"])
    m.fs.hx.cold_side.deltaP[t0].fix(geom["cold_side_deltaP"])

    m.fs.hx.area.fix(geom["area"])
    m.fs.hx.heat_transfer_coefficient[t0].fix(geom["overall_U"])
    m.fs.hx.effectiveness[t0].fix(0.80)

    duty_guess = 0.8 * min(hot["flow_mol"], cold["flow_mol"]) * 100 * (hot["temperature"] - cold["temperature"])
    duty_guess = max(duty_guess, 1e5)
    m.fs.hx.initialize(outlvl=idaeslog.INFO, duty=(duty_guess, pyunits.W))
    solver = get_solver()
    res = solver.solve(m, tee=False)

    print("\n--- IDAES HeatExchangerNTU ---")
    print(f"Solver termination: {res.solver.termination_condition}")
    print(f"Area: {pyo.value(m.fs.hx.area):.3f} m^2")
    print(f"U: {pyo.value(m.fs.hx.heat_transfer_coefficient[t0]):.3f} W/m^2/K")
    print(f"Effectiveness: {pyo.value(m.fs.hx.effectiveness[t0]):.4f}")
    print(f"NTU: {pyo.value(m.fs.hx.NTU[t0]):.4f}")
    print(f"Cratio: {pyo.value(m.fs.hx.Cratio[t0]):.4f}")
    print(f"Heat duty: {pyo.value(m.fs.hx.hot_side.heat[t0]) / 1000:.3f} kW")
    print(
        f"Hot in/out T: {pyo.value(m.fs.hx.hot_side.properties_in[t0].temperature):.2f} -> "
        f"{pyo.value(m.fs.hx.hot_side.properties_out[t0].temperature):.2f} K"
    )
    print(
        f"Cold in/out T: {pyo.value(m.fs.hx.cold_side.properties_in[t0].temperature):.2f} -> "
        f"{pyo.value(m.fs.hx.cold_side.properties_out[t0].temperature):.2f} K"
    )
    print(
        f"Hot dP: {-pyo.value(m.fs.hx.hot_side.deltaP[t0]) / 1000:.2f} kPa | "
        f"Cold dP: {-pyo.value(m.fs.hx.cold_side.deltaP[t0]) / 1000:.2f} kPa"
    )


if __name__ == "__main__":
    import sys
    
    system_name = sys.argv[1] if len(sys.argv) > 1 else "water_toluene_nitrogen"
    main(system_name=system_name)

