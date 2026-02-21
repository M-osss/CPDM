"""IDAES LMTD heat exchanger. Edit system presets (e.g. HYDROCARBON_SYSTEM) in hex_component_properties for flow/T/P/composition; INPUT_PARAMETERS for geometry only."""
import pyomo.environ as pyo

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

from idaes.models.unit_models import HeatExchanger
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

import idaes_property_builder as ipb
import hex_component_properties as hcp

#specify the desired system here and make sure it matches with what you want to simulate in the hex_component_properties.py file
#matching starts on line 300 of hex_component_properties.py file. 
def main(system_name: str = "hydrocarbon", geometry_overrides: dict = None):
    system = hcp.get_system(system_name)
    hcp.validate_system(system)
    
    geom = hcp.get_geometry("0D")
    if geometry_overrides:
        geom.update(geometry_overrides)
    
    print(f"\n{'='*80}")
    print(f"IDAES HeatExchanger (LMTD): {system['name']}")
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

    m.fs.hx = HeatExchanger(
        hot_side={"property_package": m.fs.props, "has_pressure_change": True},
        cold_side={"property_package": m.fs.props, "has_pressure_change": True},
    )

    t0 = m.fs.time.first()

    comp_id_map = {}
    for cas, built_comp in zip(components_list, built):
        comp_id_map[cas] = built_comp.comp_id

    hot = system["hot_stream"]
    m.fs.hx.hot_side.properties_in[t0].flow_mol.fix(hot["flow_mol"])
    m.fs.hx.hot_side.properties_in[t0].pressure.fix(hot["pressure"])
    m.fs.hx.hot_side.properties_in[t0].temperature.fix(hot["temperature"])
    for cas, mole_frac in hot["composition"].items():
        comp_id = comp_id_map[cas]
        m.fs.hx.hot_side.properties_in[t0].mole_frac_comp[comp_id].fix(mole_frac)

    cold = system["cold_stream"]
    m.fs.hx.cold_side.properties_in[t0].flow_mol.fix(cold["flow_mol"])
    m.fs.hx.cold_side.properties_in[t0].pressure.fix(cold["pressure"])
    m.fs.hx.cold_side.properties_in[t0].temperature.fix(cold["temperature"])
    for cas, mole_frac in cold["composition"].items():
        comp_id = comp_id_map[cas]
        m.fs.hx.cold_side.properties_in[t0].mole_frac_comp[comp_id].fix(mole_frac)

    m.fs.hx.hot_side.deltaP[t0].fix(geom["hot_side_deltaP"])
    m.fs.hx.cold_side.deltaP[t0].fix(geom["cold_side_deltaP"])

    hot_out_T = geom.get("hot_outlet_temperature")
    if hot_out_T is None:
        hot_out_T = hot["temperature"] - 30.0
    m.fs.hx.hot_side.properties_out[t0].temperature.fix(hot_out_T)
    m.fs.hx.overall_heat_transfer_coefficient[t0].fix(geom["overall_U"])

    m.fs.hx.initialize(outlvl=idaeslog.INFO)
    solver = get_solver()
    res = solver.solve(m, tee=False)

    comp_id_to_name = {comp_id_map[c["identifier"]]: c["name"] for c in system["components"]}

    def _fmt_comp(props_block):
        return ", ".join(
            f"{comp_id_to_name.get(c, c)}={pyo.value(props_block.mole_frac_comp[c]):.4f}"
            for c in m.fs.props.component_list
        )

    print("\n--- IDAES HeatExchanger (LMTD, database-driven properties) ---")
    print(f"Solver termination: {res.solver.termination_condition}")
    print(f"Heat duty: {pyo.value(m.fs.hx.heat_duty[t0]) / 1000:.3f} kW")
    print(f"Area: {pyo.value(m.fs.hx.area):.3f} m^2")
    print(f"LMTD (delta_temperature): {pyo.value(m.fs.hx.delta_temperature[t0]):.3f} K")

    print("\n--- HOT STREAM ---")
    T_hi = pyo.value(m.fs.hx.hot_side.properties_in[t0].temperature)
    T_ho = pyo.value(m.fs.hx.hot_side.properties_out[t0].temperature)
    P_hi = pyo.value(m.fs.hx.hot_side.properties_in[t0].pressure)
    P_ho = pyo.value(m.fs.hx.hot_side.properties_out[t0].pressure)
    print(f"Inlet:  T={T_hi:.2f} K  P={P_hi/1000:.1f} kPa  flow={pyo.value(m.fs.hx.hot_side.properties_in[t0].flow_mol):.2f} mol/s")
    print(f"        composition: {_fmt_comp(m.fs.hx.hot_side.properties_in[t0])}")
    print(f"Outlet: T={T_ho:.2f} K  P={P_ho/1000:.1f} kPa  dP={-pyo.value(m.fs.hx.hot_side.deltaP[t0])/1000:.1f} kPa")
    print(f"        composition: {_fmt_comp(m.fs.hx.hot_side.properties_out[t0])}")

    print("\n--- COLD STREAM ---")
    T_ci = pyo.value(m.fs.hx.cold_side.properties_in[t0].temperature)
    T_co = pyo.value(m.fs.hx.cold_side.properties_out[t0].temperature)
    P_ci = pyo.value(m.fs.hx.cold_side.properties_in[t0].pressure)
    P_co = pyo.value(m.fs.hx.cold_side.properties_out[t0].pressure)
    print(f"Inlet:  T={T_ci:.2f} K  P={P_ci/1000:.1f} kPa  flow={pyo.value(m.fs.hx.cold_side.properties_in[t0].flow_mol):.2f} mol/s")
    print(f"        composition: {_fmt_comp(m.fs.hx.cold_side.properties_in[t0])}")
    print(f"Outlet: T={T_co:.2f} K  P={P_co/1000:.1f} kPa  dP={-pyo.value(m.fs.hx.cold_side.deltaP[t0])/1000:.1f} kPa")
    print(f"        composition: {_fmt_comp(m.fs.hx.cold_side.properties_out[t0])}")


if __name__ == "__main__":
    import sys
    
    system_name = sys.argv[1] if len(sys.argv) > 1 else "hydrocarbon"
    main(system_name=system_name)

