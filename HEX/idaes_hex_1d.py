"""IDAES 1D heat exchanger. Edit system presets in hex_component_properties for flow/T/P/composition; INPUT_PARAMETERS for geometry only."""
import pyomo.environ as pyo
from pyomo.environ import units as pyunits

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

from idaes.models.unit_models import HeatExchanger1D, HeatExchangerFlowPattern
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

import idaes_property_builder as ipb
import property_database_interface as pdb
import hex_component_properties as hcp


def _check_two_phase_warning(m, t0, x0, xm, x1):
    has_vap_phase = False
    try:
        if hasattr(m.fs.props, "phase_list"):
            has_vap_phase = "Vap" in [str(p) for p in m.fs.props.phase_list]
    except Exception:
        pass
    two_phase_detected = False
    warning_locations = []
    
    if has_vap_phase:
        check_points = [
            ("hot_inlet", m.fs.hx.hot_side.properties[t0, x0]),
            ("hot_mid", m.fs.hx.hot_side.properties[t0, xm]),
            ("hot_outlet", m.fs.hx.hot_side.properties[t0, x1]),
            ("cold_inlet", m.fs.hx.cold_side.properties[t0, x1]),
            ("cold_mid", m.fs.hx.cold_side.properties[t0, xm]),
            ("cold_outlet", m.fs.hx.cold_side.properties[t0, x0]),
        ]
        
        for loc_name, props_block in check_points:
            try:
                if hasattr(props_block, "vapor_fraction"):
                    vf = pyo.value(props_block.vapor_fraction)
                    if vf > 1e-6 and vf < (1.0 - 1e-6):
                        two_phase_detected = True
                        warning_locations.append(f"{loc_name} (VF={vf:.4f})")
                elif hasattr(props_block, "phase_frac"):
                    if "Vap" in props_block.phase_frac:
                        vf = pyo.value(props_block.phase_frac["Vap"])
                        if vf > 1e-6 and vf < (1.0 - 1e-6):
                            two_phase_detected = True
                            warning_locations.append(f"{loc_name} (Vap frac={vf:.4f})")
            except (AttributeError, KeyError, ValueError, TypeError):
                continue
    else:
        try:
            tc_values = []
            if hasattr(m.fs.props, "component_list"):
                for comp in m.fs.props.component_list:
                    try:
                        tc = pyo.value(m.fs.props.temperature_crit[comp])
                        tc_values.append(tc)
                    except (AttributeError, KeyError):
                        continue
            
            if tc_values:
                avg_tc = sum(tc_values) / len(tc_values)
                check_points = [
                    ("hot_inlet", m.fs.hx.hot_side.properties[t0, x0]),
                    ("hot_outlet", m.fs.hx.hot_side.properties[t0, x1]),
                    ("cold_inlet", m.fs.hx.cold_side.properties[t0, x1]),
                    ("cold_outlet", m.fs.hx.cold_side.properties[t0, x0]),
                ]
                
                for loc_name, props_block in check_points:
                    try:
                        T = pyo.value(props_block.temperature)
                        P = pyo.value(props_block.pressure)
                        if avg_tc > 0:
                            if (T > 0.7 * avg_tc and P < 1e6) or T > 500:
                                warning_locations.append(
                                    f"{loc_name} (T={T:.1f}K, P={P/1e5:.2f}bar)"
                                )
                    except (AttributeError, ValueError, TypeError):
                        continue
        except Exception:
            pass

    check_points = [
        ("hot_inlet", m.fs.hx.hot_side.properties[t0, x0]),
        ("hot_mid", m.fs.hx.hot_side.properties[t0, xm]),
        ("hot_outlet", m.fs.hx.hot_side.properties[t0, x1]),
        ("cold_inlet", m.fs.hx.cold_side.properties[t0, x1]),
        ("cold_mid", m.fs.hx.cold_side.properties[t0, xm]),
        ("cold_outlet", m.fs.hx.cold_side.properties[t0, x0]),
    ]
    
    for loc_name, props_block in check_points:
        try:
            # Access temperature and pressure from the properties block
            # props_block is already indexed by [t0, x], so access directly
            T = pyo.value(props_block.temperature)
            P = pyo.value(props_block.pressure)
            # Warn if temperature is very high (>500K) and pressure is moderate
            # This is a conservative check for potential vaporization
            if T > 500 and P < 2e6:  # >500K and <20 bar
                if not any(loc_name in loc for loc in warning_locations):
                    warning_locations.append(
                        f"{loc_name} (T={T:.1f}K, P={P/1e5:.2f}bar)"
                    )
        except (AttributeError, ValueError, TypeError, KeyError):
            # Skip if property not accessible
            continue
    
    if two_phase_detected or warning_locations:
        print("\n" + "!"*80)
        print("WARNING: Two-phase conditions detected in heat exchanger!")
        print("!"*80)
        print("Heat exchangers are typically designed for single-phase operation.")
        print("Two-phase flow can cause:")
        print("  - Incorrect heat transfer calculations")
        print("  - Flow maldistribution")
        print("  - Pressure drop issues")
        print("  - Mechanical problems (vibration, erosion)")
        if warning_locations:
            print(f"\nDetected at locations: {', '.join(warning_locations)}")
        if not has_vap_phase:
            print("\nNote: Using liquid-only property package.")
            print("For accurate two-phase detection, use VLE (vapor-liquid equilibrium) property package.")
        print("!"*80 + "\n")


def main(system_name: str = "water_toluene", geometry_overrides: dict = None):
    system = hcp.get_system(system_name)
    hcp.validate_system(system)
    
    geom = hcp.get_geometry("1D")
    if geometry_overrides:
        geom.update(geometry_overrides)
    
    print(f"\n{'='*80}")
    print(f"IDAES HeatExchanger1D: {system['name']}")
    print(f"Components: {', '.join(c['name'] for c in system['components'])}")
    print(f"{'='*80}\n")
    
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    components_list = [c["identifier"] for c in system["components"]]
    props_cfg, built = ipb.build_generic_parameter_block_config(
        components_list,
        phases=("Liq",),
        eos="PR",
        prefer_rpp_ig_cp=True,
    )
    ipb.validate_config_dict(props_cfg)
    m.fs.props = GenericParameterBlock(**props_cfg)

    m.fs.hx = HeatExchanger1D(
        hot_side={"property_package": m.fs.props, "has_pressure_change": False},
        cold_side={"property_package": m.fs.props, "has_pressure_change": False},
        flow_type=HeatExchangerFlowPattern.countercurrent,
        finite_elements=geom["finite_elements"],
        collocation_points=geom["collocation_points"],
    )

    t0 = m.fs.time.first()

    comp_id_map = {}
    for cas, built_comp in zip(components_list, built):
        comp_id_map[cas] = built_comp.comp_id

    m.fs.hx.area.fix(geom["area"])
    m.fs.hx.length.fix(geom["length"])
    m.fs.hx.heat_transfer_coefficient.fix(geom["heat_transfer_coefficient"])

    hot = system["hot_stream"]
    m.fs.hx.hot_side_inlet.flow_mol[t0].fix(hot["flow_mol"])
    m.fs.hx.hot_side_inlet.temperature[t0].fix(hot["temperature"])
    m.fs.hx.hot_side_inlet.pressure[t0].fix(hot["pressure"])
    hot_mole_fracs = {}
    hot_sum = sum(hot["composition"].values())
    missing_count = 0
    for cas in components_list:
        comp_id = comp_id_map[cas]
        if cas in hot["composition"]:
            hot_mole_fracs[comp_id] = hot["composition"][cas]
        else:
            hot_mole_fracs[comp_id] = 1e-20
            missing_count += 1
    
    total = sum(hot_mole_fracs.values())
    if missing_count > 0 and abs(total - 1.0) > 1e-6:
        adjustment = (1.0 - hot_sum) / missing_count
        for cas in components_list:
            comp_id = comp_id_map[cas]
            if cas not in hot["composition"]:
                hot_mole_fracs[comp_id] = max(1e-20, adjustment)
    
    for comp_id, mole_frac in hot_mole_fracs.items():
        m.fs.hx.hot_side_inlet.mole_frac_comp[t0, comp_id].fix(mole_frac)

    cold = system["cold_stream"]
    m.fs.hx.cold_side_inlet.flow_mol[t0].fix(cold["flow_mol"])
    m.fs.hx.cold_side_inlet.temperature[t0].fix(cold["temperature"])
    m.fs.hx.cold_side_inlet.pressure[t0].fix(cold["pressure"])
    cold_mole_fracs = {}
    cold_sum = sum(cold["composition"].values())
    missing_count = 0
    for cas in components_list:
        comp_id = comp_id_map[cas]
        if cas in cold["composition"]:
            cold_mole_fracs[comp_id] = cold["composition"][cas]
        else:
            cold_mole_fracs[comp_id] = 1e-20
            missing_count += 1
    
    # Normalize only if we added missing components and sum != 1.0
    total = sum(cold_mole_fracs.values())
    if missing_count > 0 and abs(total - 1.0) > 1e-6:
        # Adjust the small values proportionally to maintain sum = 1.0
        adjustment = (1.0 - cold_sum) / missing_count
        for cas in components_list:
            comp_id = comp_id_map[cas]
            if cas not in cold["composition"]:
                cold_mole_fracs[comp_id] = max(1e-20, adjustment)
    
    for comp_id, mole_frac in cold_mole_fracs.items():
        m.fs.hx.cold_side_inlet.mole_frac_comp[t0, comp_id].fix(mole_frac)

    # Initialize and solve
    # Provide state arguments to help initialization
    hot_state_args = {
        "flow_mol": hot["flow_mol"],
        "temperature": hot["temperature"],
        "pressure": hot["pressure"],
    }
    cold_state_args = {
        "flow_mol": cold["flow_mol"],
        "temperature": cold["temperature"],
        "pressure": cold["pressure"],
    }
    
    # Estimate initial heat duty: Q = U*A*LMTD (rough estimate)
    # LMTD ≈ (Th_in - Tc_out) - (Th_out - Tc_in) / ln((Th_in - Tc_out)/(Th_out - Tc_in))
    # Simplified: assume Th_out ≈ Th_in - 50K, Tc_out ≈ Tc_in + 50K
    Th_in = hot["temperature"]
    Tc_in = cold["temperature"]
    Th_out_guess = Th_in - 50.0
    Tc_out_guess = Tc_in + 50.0
    dT1 = Th_in - Tc_out_guess
    dT2 = Th_out_guess - Tc_in
    if dT1 > 0 and dT2 > 0 and abs(dT1 - dT2) > 1e-6:
        lmtd = (dT1 - dT2) / pyo.log(dT1 / dT2)
    else:
        lmtd = (dT1 + dT2) / 2.0
    duty_guess = geom["heat_transfer_coefficient"] * geom["area"] * lmtd  # W
    
    try:
        m.fs.hx.initialize(
            outlvl=idaeslog.INFO,
            hot_side_state_args=hot_state_args,
            cold_side_state_args=cold_state_args,
            duty=(duty_guess, pyunits.W),
        )
    except Exception as e:
        print(f"\nWarning: Initialization had issues: {type(e).__name__}")
        print("This may occur with gas-liquid mixtures using liquid-only property package.")
        print("Attempting to solve anyway...")
    solver = get_solver()
    res = solver.solve(m, tee=False)

    x0 = m.fs.hx.hot_side.length_domain.first()
    x1 = m.fs.hx.hot_side.length_domain.last()
    xm = sorted(list(m.fs.hx.hot_side.length_domain))[len(list(m.fs.hx.hot_side.length_domain)) // 2]

    # Check for two-phase conditions (run regardless of solver status)
    try:
        _check_two_phase_warning(m, t0, x0, xm, x1)
    except Exception as e:
        # Log but don't fail - the check is advisory
        print(f"\nNote: Two-phase check encountered an issue (non-critical): {type(e).__name__}")

    print("\n--- IDAES HeatExchanger1D (database-driven properties) ---")
    print(f"Solver termination: {res.solver.termination_condition}")
    print(f"Area: {pyo.value(m.fs.hx.area):.3f} m^2 | Length: {pyo.value(m.fs.hx.length):.3f} m")
    print(f"h (fixed): {pyo.value(m.fs.hx.heat_transfer_coefficient[t0, x0]):.2f} W/m^2/K")
    print(f"Hot out T: {pyo.value(m.fs.hx.hot_side_outlet.temperature[t0]):.2f} K")
    print(f"Cold out T: {pyo.value(m.fs.hx.cold_side_outlet.temperature[t0]):.2f} K")
    print(
        "Hot-side T profile (x=0, mid, 1): "
        f"{pyo.value(m.fs.hx.hot_side.properties[t0, x0].temperature):.2f}, "
        f"{pyo.value(m.fs.hx.hot_side.properties[t0, xm].temperature):.2f}, "
        f"{pyo.value(m.fs.hx.hot_side.properties[t0, x1].temperature):.2f} K"
    )
    print(
        "Cold-side T profile (x=0, mid, 1): "
        f"{pyo.value(m.fs.hx.cold_side.properties[t0, x0].temperature):.2f}, "
        f"{pyo.value(m.fs.hx.cold_side.properties[t0, xm].temperature):.2f}, "
        f"{pyo.value(m.fs.hx.cold_side.properties[t0, x1].temperature):.2f} K"
    )

    # Transport correlations availability (ChemSep XML) – for correlation-based h/dP work.
    print("\nTransport correlation presence (ChemSep XML):")
    for spec in ["7732-18-5", "108-88-3"]:
        av = pdb.check_property_availability(spec)
        print(
            f"- {spec}: "
            f"liq_visc={av.get('transport_present', {}).get('liq_visc_present')}, "
            f"vap_visc={av.get('transport_present', {}).get('vap_visc_present')}, "
            f"liq_k={av.get('transport_present', {}).get('liq_k_present')}, "
            f"vap_k={av.get('transport_present', {}).get('vap_k_present')}"
        )


if __name__ == "__main__":
    import sys
    
    system_name = sys.argv[1] if len(sys.argv) > 1 else "water_toluene_nitrogen"
    main(system_name=system_name)

