"""IDAES shell-and-tube 1D. Edit system presets in hex_component_properties for flow/T/P/composition; INPUT_PARAMETERS for geometry only."""
import pyomo.environ as pyo

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

from idaes.models.unit_models import ShellAndTube1D, HeatExchangerFlowPattern
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

import idaes_property_builder as ipb
import property_database_interface as pdb
import hex_component_properties as hcp


def main(system_name: str = "water_toluene", geometry_overrides: dict = None):
    system = hcp.get_system(system_name)
    hcp.validate_system(system)
    
    geom = hcp.get_geometry("ShellAndTube1D")
    if geometry_overrides:
        geom.update(geometry_overrides)
    
    print(f"\n{'='*80}")
    print(f"IDAES ShellAndTube1D: {system['name']}")
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

    m.fs.hx = ShellAndTube1D(
        hot_side={"property_package": m.fs.props, "has_pressure_change": False},
        cold_side={"property_package": m.fs.props, "has_pressure_change": False},
        shell_is_hot=True,
        flow_type=HeatExchangerFlowPattern.countercurrent,
        finite_elements=geom["finite_elements"],
        collocation_points=geom["collocation_points"],
    )

    t0 = m.fs.time.first()

    comp_id_map = {}
    for cas, built_comp in zip(components_list, built):
        comp_id_map[cas] = built_comp.comp_id

    m.fs.hx.length.fix(geom["length"])
    m.fs.hx.shell_diameter.fix(geom["shell_diameter"])
    m.fs.hx.tube_outer_diameter.fix(geom["tube_outer_diameter"])
    m.fs.hx.tube_inner_diameter.fix(geom["tube_inner_diameter"])
    m.fs.hx.number_of_tubes.fix(geom["number_of_tubes"])

    m.fs.hx.hot_side_heat_transfer_coefficient.fix(geom["hot_side_htc"])
    m.fs.hx.cold_side_heat_transfer_coefficient.fix(geom["cold_side_htc"])

    hot = system["hot_stream"]
    m.fs.hx.hot_side_inlet.flow_mol[t0].fix(hot["flow_mol"])
    m.fs.hx.hot_side_inlet.temperature[t0].fix(hot["temperature"])
    m.fs.hx.hot_side_inlet.pressure[t0].fix(hot["pressure"])
    for cas, mole_frac in hot["composition"].items():
        m.fs.hx.hot_side_inlet.mole_frac_comp[t0, comp_id_map[cas]].fix(mole_frac)

    cold = system["cold_stream"]
    m.fs.hx.cold_side_inlet.flow_mol[t0].fix(cold["flow_mol"])
    m.fs.hx.cold_side_inlet.temperature[t0].fix(cold["temperature"])
    m.fs.hx.cold_side_inlet.pressure[t0].fix(cold["pressure"])
    for cas, mole_frac in cold["composition"].items():
        m.fs.hx.cold_side_inlet.mole_frac_comp[t0, comp_id_map[cas]].fix(mole_frac)

    m.fs.hx.initialize(outlvl=idaeslog.INFO)
    solver = get_solver()
    res = solver.solve(m, tee=False)

    x0 = m.fs.hx.hot_side.length_domain.first()
    x1 = m.fs.hx.hot_side.length_domain.last()
    xm = sorted(list(m.fs.hx.hot_side.length_domain))[len(list(m.fs.hx.hot_side.length_domain)) // 2]

    print("\n--- IDAES ShellAndTube1D ---")
    print(f"Solver termination: {res.solver.termination_condition}")
    print(f"Length: {pyo.value(m.fs.hx.length):.3f} m")
    print(f"Shell D: {pyo.value(m.fs.hx.shell_diameter):.3f} m")
    print(
        f"Tubes: N={pyo.value(m.fs.hx.number_of_tubes):.0f} | "
        f"Do={pyo.value(m.fs.hx.tube_outer_diameter):.4f} m | "
        f"Di={pyo.value(m.fs.hx.tube_inner_diameter):.4f} m"
    )
    print(f"Hot out T: {pyo.value(m.fs.hx.hot_side_outlet.temperature[t0]):.2f} K")
    print(f"Cold out T: {pyo.value(m.fs.hx.cold_side_outlet.temperature[t0]):.2f} K")
    print(
        "Wall T profile (x=0, mid, 1): "
        f"{pyo.value(m.fs.hx.temperature_wall[t0, x0]):.2f}, "
        f"{pyo.value(m.fs.hx.temperature_wall[t0, xm]):.2f}, "
        f"{pyo.value(m.fs.hx.temperature_wall[t0, x1]):.2f} K"
    )

    kij = pdb.get_pr_binary_interaction("7732-18-5", "108-88-3")
    print(f"PR k_ij (ChemSep ipd/pr.ipd): {kij if kij is not None else 0.0}")


if __name__ == "__main__":
    import sys
    
    system_name = sys.argv[1] if len(sys.argv) > 1 else "water_toluene"
    main(system_name=system_name)

