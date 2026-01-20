import pyomo.environ as pyo
from pyomo.environ import (Constraint,
                           log,
                           units as pyunits)

from idaes.core import FlowsheetBlock, LiquidPhase
from idaes.models.unit_models import HeatExchanger
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock, Component)
from idaes.models.properties.modular_properties.pure import Perrys, NIST
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

def main():
    
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Using Peng-Robinson equation of state for realistic liquid behavior
    # and temperature-dependent property correlations
    props_config = {
        "components": {
            "h2o": {
                "type": Component,
                "dens_mol_liq_comp": Perrys.dens_mol_liq_comp,
                "enth_mol_liq_comp": Perrys.enth_mol_liq_comp,
                "enth_mol_ig_comp": NIST.enth_mol_ig_comp,
                "cp_mol_liq_comp": Perrys.cp_mol_liq_comp,
                "cp_mol_ig_comp": NIST.cp_mol_ig_comp,
                "entr_mol_liq_comp": Perrys.entr_mol_liq_comp,
                "entr_mol_ig_comp": NIST.entr_mol_ig_comp,
                "parameter_data": {
                    "mw": (18.015e-3, pyunits.kg/pyunits.mol),  # molecular weight
                    "pressure_crit": (220.64e5, pyunits.Pa),
                    "temperature_crit": (647.14, pyunits.K),
                    "omega": 0.344,  # acentric factor
                    # Perry's parameters for liquid density
                    "dens_mol_liq_comp_coeff": {
                        'eqn_type': 1,
                        '1': (5.459, pyunits.kmol/pyunits.m**3),
                        '2': (0.30542, pyunits.dimensionless),
                        '3': (647.13, pyunits.K),
                        '4': (0.081, pyunits.dimensionless),
                    },
                    # Perry's parameters for liquid heat capacity
                    "cp_mol_liq_comp_coeff": {
                        '1': (2.7637e5, pyunits.J/pyunits.mol/pyunits.K),
                        '2': (-2.0901e3, pyunits.J/pyunits.mol/pyunits.K**2),
                        '3': (8.1250, pyunits.J/pyunits.mol/pyunits.K**3),
                        '4': (-1.4116e-2, pyunits.J/pyunits.mol/pyunits.K**4),
                        '5': (9.3701e-6, pyunits.J/pyunits.mol/pyunits.K**5),
                    },
                    # NIST Shomate parameters for ideal gas heat capacity
                    "cp_mol_ig_comp_coeff": {
                        'A': 30.092,
                        'B': 6.832514,
                        'C': 6.793435,
                        'D': -2.53448,
                        'E': 0.082139,
                        'F': -250.881,
                        'G': 223.3967,
                        'H': -241.8264,
                    },
                    # Enthalpy of formation
                    "enth_mol_form_liq_comp_ref": (-285.83e3, pyunits.J/pyunits.mol),
                    "enth_mol_form_vap_comp_ref": (-241.82e3, pyunits.J/pyunits.mol),
                    # Standard entropy
                    "entr_mol_form_liq_comp_ref": (69.95, pyunits.J/pyunits.mol/pyunits.K),
                    "entr_mol_form_vap_comp_ref": (188.84, pyunits.J/pyunits.mol/pyunits.K),
                },
            },
            "toluene": {
                "type": Component,
                "dens_mol_liq_comp": Perrys.dens_mol_liq_comp,
                "enth_mol_liq_comp": Perrys.enth_mol_liq_comp,
                "enth_mol_ig_comp": NIST.enth_mol_ig_comp,
                "cp_mol_liq_comp": Perrys.cp_mol_liq_comp,
                "cp_mol_ig_comp": NIST.cp_mol_ig_comp,
                "entr_mol_liq_comp": Perrys.entr_mol_liq_comp,
                "entr_mol_ig_comp": NIST.entr_mol_ig_comp,
                "parameter_data": {
                    "mw": (92.138e-3, pyunits.kg/pyunits.mol),  # molecular weight
                    "pressure_crit": (41.06e5, pyunits.Pa),
                    "temperature_crit": (591.8, pyunits.K),
                    "omega": 0.262,  # acentric factor
                    # Perry's parameters for liquid density
                    "dens_mol_liq_comp_coeff": {
                        'eqn_type': 1,
                        '1': (0.8488, pyunits.kmol/pyunits.m**3),
                        '2': (0.26655, pyunits.dimensionless),
                        '3': (591.8, pyunits.K),
                        '4': (0.2878, pyunits.dimensionless),
                    },
                    # Perry's parameters for liquid heat capacity
                    "cp_mol_liq_comp_coeff": {
                        '1': (1.4014e5, pyunits.J/pyunits.mol/pyunits.K),
                        '2': (-1.5230e2, pyunits.J/pyunits.mol/pyunits.K**2),
                        '3': (6.9500e-1, pyunits.J/pyunits.mol/pyunits.K**3),
                        '4': (0.0, pyunits.J/pyunits.mol/pyunits.K**4),
                        '5': (0.0, pyunits.J/pyunits.mol/pyunits.K**5),
                    },
                    # NIST Shomate parameters for ideal gas heat capacity
                    "cp_mol_ig_comp_coeff": {
                        'A': -24.35,
                        'B': 512.5,
                        'C': -276.5,
                        'D': 49.11,
                        'E': 0.0,
                        'F': 0.0,
                        'G': 0.0,
                        'H': 50.0,
                    },
                    # Enthalpy of formation
                    "enth_mol_form_liq_comp_ref": (12.0e3, pyunits.J/pyunits.mol),
                    "enth_mol_form_vap_comp_ref": (50.0e3, pyunits.J/pyunits.mol),
                    # Standard entropy
                    "entr_mol_form_liq_comp_ref": (220.96, pyunits.J/pyunits.mol/pyunits.K),
                    "entr_mol_form_vap_comp_ref": (320.77, pyunits.J/pyunits.mol/pyunits.K),
                },
            }
        },
        "phases": {
            "Liq": {
                "type": LiquidPhase,
                "equation_of_state": Cubic,
                "equation_of_state_options": {"type": CubicType.PR},  # Peng-Robinson EoS
            }
        },
        "base_units": {
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        "state_definition": FTPx,
        "state_bounds": {
            "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
            "temperature": (273.15, 300, 600, pyunits.K),
            "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
        },
        "pressure_ref": (101325, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
        "parameter_data": {
            "PR_kappa": {
                ("h2o", "h2o"): 0.000,
                ("h2o", "toluene"): 0.000,
                ("toluene", "h2o"): 0.000,
                ("toluene", "toluene"): 0.000,
            }
        },
    }
    m.fs.props = GenericParameterBlock(**props_config)

    # Create heat exchanger with realistic configuration
    m.fs.hex = HeatExchanger(
        hot_side={"property_package": m.fs.props, 
                  "has_pressure_change": True},  # Account for pressure drop
        cold_side={"property_package": m.fs.props,
                   "has_pressure_change": True},  # Account for pressure drop
    )

    # Get time index
    t0 = m.fs.time.first()
    
    # Fix hot side inlet conditions
    m.fs.hex.hot_side.properties_in[t0].flow_mol.fix(100)  # mol/s
    m.fs.hex.hot_side.properties_in[t0].pressure.fix(101325)  # Pa
    m.fs.hex.hot_side.properties_in[t0].temperature.fix(473.15) # K (200°C)
    m.fs.hex.hot_side.properties_in[t0].mole_frac_comp['h2o'].fix(0.999)
    m.fs.hex.hot_side.properties_in[t0].mole_frac_comp['toluene'].fix(0.001)

    # Fix cold side inlet conditions
    m.fs.hex.cold_side.properties_in[t0].flow_mol.fix(100)  # mol/s
    m.fs.hex.cold_side.properties_in[t0].pressure.fix(202650)  # Pa
    m.fs.hex.cold_side.properties_in[t0].temperature.fix(323.15)  # K (50°C)
    m.fs.hex.cold_side.properties_in[t0].mole_frac_comp['toluene'].fix(0.999)
    m.fs.hex.cold_side.properties_in[t0].mole_frac_comp['h2o'].fix(0.001)
    
    # Fix hot side outlet temperature
    m.fs.hex.hot_side.properties_out[t0].temperature.fix(373.15) # K (100°C)
    
    # Fix pressure changes (realistic pressure drops)
    m.fs.hex.hot_side.deltaP[t0].fix(-5000)  # 5 kPa pressure drop
    m.fs.hex.cold_side.deltaP[t0].fix(-10000)  # 10 kPa pressure drop
    
    # Fix overall heat transfer coefficient (realistic value)
    m.fs.hex.overall_heat_transfer_coefficient[t0].fix(800)  # W/m2/K
    
    # Initialize the heat exchanger
    m.fs.hex.initialize(outlvl=idaeslog.INFO)

    # Solve the model
    solver = get_solver()
    results = solver.solve(m, tee=True)

    # Print results
    print("\n--- IDAES HEAT EXCHANGER RESULTS (REALISTIC MODEL) ---")
    print(f"Heat Duty: {pyo.value(m.fs.hex.heat_duty[t0]) / 1000:.2f} kW")
    print(f"Heat Transfer Area: {pyo.value(m.fs.hex.area):.2f} m^2")
    print("\n--- HOT STREAM (Water) ---")
    print(f"Inlet Temperature: {pyo.value(m.fs.hex.hot_side.properties_in[t0].temperature):.2f} K ({pyo.value(m.fs.hex.hot_side.properties_in[t0].temperature) - 273.15:.1f}°C)")
    print(f"Outlet Temperature: {pyo.value(m.fs.hex.hot_side.properties_out[t0].temperature):.2f} K ({pyo.value(m.fs.hex.hot_side.properties_out[t0].temperature) - 273.15:.1f}°C)")
    print(f"Inlet Pressure: {pyo.value(m.fs.hex.hot_side.properties_in[t0].pressure)/1000:.1f} kPa")
    print(f"Outlet Pressure: {pyo.value(m.fs.hex.hot_side.properties_out[t0].pressure)/1000:.1f} kPa")
    print(f"Pressure Drop: {-pyo.value(m.fs.hex.hot_side.deltaP[t0])/1000:.1f} kPa")
    
    # Print hot stream properties
    print(f"Inlet Density: {pyo.value(m.fs.hex.hot_side.properties_in[t0].dens_mol):.1f} mol/m³")
    print(f"Inlet Cp: {pyo.value(m.fs.hex.hot_side.properties_in[t0].cp_mol):.1f} J/mol/K")
    
    print("\n--- COLD STREAM (Toluene) ---")
    print(f"Inlet Temperature: {pyo.value(m.fs.hex.cold_side.properties_in[t0].temperature):.2f} K ({pyo.value(m.fs.hex.cold_side.properties_in[t0].temperature) - 273.15:.1f}°C)")
    print(f"Outlet Temperature: {pyo.value(m.fs.hex.cold_side.properties_out[t0].temperature):.2f} K ({pyo.value(m.fs.hex.cold_side.properties_out[t0].temperature) - 273.15:.1f}°C)")
    print(f"Inlet Pressure: {pyo.value(m.fs.hex.cold_side.properties_in[t0].pressure)/1000:.1f} kPa")
    print(f"Outlet Pressure: {pyo.value(m.fs.hex.cold_side.properties_out[t0].pressure)/1000:.1f} kPa")
    print(f"Pressure Drop: {-pyo.value(m.fs.hex.cold_side.deltaP[t0])/1000:.1f} kPa")
    
    # Print cold stream properties
    print(f"Inlet Density: {pyo.value(m.fs.hex.cold_side.properties_in[t0].dens_mol):.1f} mol/m³")
    print(f"Inlet Cp: {pyo.value(m.fs.hex.cold_side.properties_in[t0].cp_mol):.1f} J/mol/K")
    
    print("\n--- HEAT EXCHANGER ---")
    print(f"LMTD: {pyo.value(m.fs.hex.delta_temperature[t0]):.2f} K")
    print(f"Overall HTC: {pyo.value(m.fs.hex.overall_heat_transfer_coefficient[t0]):.1f} W/m^2/K")
    
    # Calculate and display effectiveness
    C_hot = m.fs.hex.hot_side.properties_in[t0].flow_mol * \
            m.fs.hex.hot_side.properties_in[t0].cp_mol
    C_cold = m.fs.hex.cold_side.properties_in[t0].flow_mol * \
             m.fs.hex.cold_side.properties_in[t0].cp_mol
    C_min = min(pyo.value(C_hot), pyo.value(C_cold))
    Q_max = C_min * (pyo.value(m.fs.hex.hot_side.properties_in[t0].temperature) - 
                     pyo.value(m.fs.hex.cold_side.properties_in[t0].temperature))
    effectiveness = pyo.value(m.fs.hex.heat_duty[t0]) / Q_max
    print(f"Effectiveness: {effectiveness:.3f}")
    
    print("------------------------------------")
    
    # Check solver status
    if results.solver.termination_condition == pyo.TerminationCondition.optimal:
        print("\nSolver converged to optimal solution!")
    else:
        print(f"\nSolver terminated with condition: {results.solver.termination_condition}")


if __name__ == '__main__':
    main() 
