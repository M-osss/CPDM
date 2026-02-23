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


# edit these
HOT_FLOW_MOL       = 100
HOT_TEMPERATURE_IN  = 473.15
HOT_PRESSURE_IN     = 101325
HOT_MOLE_FRACS      = {"h2o": 0.999, "toluene": 0.001}

COLD_FLOW_MOL       = 100
COLD_TEMPERATURE_IN = 323.15
COLD_PRESSURE_IN    = 101325
COLD_MOLE_FRACS     = {"toluene": 0.999, "h2o": 0.001}

HOT_TEMPERATURE_OUT = 373.15
HOT_DELTA_P         = -5000
COLD_DELTA_P        = -10000
OVERALL_HTC         = 800


def main():

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

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
                    "mw": (18.015e-3, pyunits.kg/pyunits.mol),
                    "pressure_crit": (220.64e5, pyunits.Pa),
                    "temperature_crit": (647.14, pyunits.K),
                    "omega": 0.344,
                    "dens_mol_liq_comp_coeff": {
                        'eqn_type': 1,
                        '1': (5.459, pyunits.kmol/pyunits.m**3),
                        '2': (0.30542, pyunits.dimensionless),
                        '3': (647.13, pyunits.K),
                        '4': (0.081, pyunits.dimensionless),
                    },
                    "cp_mol_liq_comp_coeff": {
                        '1': (2.7637e5, pyunits.J/pyunits.mol/pyunits.K),
                        '2': (-2.0901e3, pyunits.J/pyunits.mol/pyunits.K**2),
                        '3': (8.1250, pyunits.J/pyunits.mol/pyunits.K**3),
                        '4': (-1.4116e-2, pyunits.J/pyunits.mol/pyunits.K**4),
                        '5': (9.3701e-6, pyunits.J/pyunits.mol/pyunits.K**5),
                    },
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
                    "enth_mol_form_liq_comp_ref": (-285.83e3, pyunits.J/pyunits.mol),
                    "enth_mol_form_vap_comp_ref": (-241.82e3, pyunits.J/pyunits.mol),
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
                    "mw": (92.138e-3, pyunits.kg/pyunits.mol),
                    "pressure_crit": (41.06e5, pyunits.Pa),
                    "temperature_crit": (591.8, pyunits.K),
                    "omega": 0.262,
                    "dens_mol_liq_comp_coeff": {
                        'eqn_type': 1,
                        '1': (0.8488, pyunits.kmol/pyunits.m**3),
                        '2': (0.26655, pyunits.dimensionless),
                        '3': (591.8, pyunits.K),
                        '4': (0.2878, pyunits.dimensionless),
                    },
                    "cp_mol_liq_comp_coeff": {
                        '1': (1.4014e5, pyunits.J/pyunits.mol/pyunits.K),
                        '2': (-1.5230e2, pyunits.J/pyunits.mol/pyunits.K**2),
                        '3': (6.9500e-1, pyunits.J/pyunits.mol/pyunits.K**3),
                        '4': (0.0, pyunits.J/pyunits.mol/pyunits.K**4),
                        '5': (0.0, pyunits.J/pyunits.mol/pyunits.K**5),
                    },
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
                    "enth_mol_form_liq_comp_ref": (12.0e3, pyunits.J/pyunits.mol),
                    "enth_mol_form_vap_comp_ref": (50.0e3, pyunits.J/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (220.96, pyunits.J/pyunits.mol/pyunits.K),
                    "entr_mol_form_vap_comp_ref": (320.77, pyunits.J/pyunits.mol/pyunits.K),
                },
            }
        },
        "phases": {
            "Liq": {
                "type": LiquidPhase,
                "equation_of_state": Cubic,
                "equation_of_state_options": {"type": CubicType.PR},
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

    m.fs.hex = HeatExchanger(
        hot_side={"property_package": m.fs.props, "has_pressure_change": True},
        cold_side={"property_package": m.fs.props, "has_pressure_change": True},
    )

    t0 = m.fs.time.first()

    m.fs.hex.hot_side.properties_in[t0].flow_mol.fix(HOT_FLOW_MOL)
    m.fs.hex.hot_side.properties_in[t0].pressure.fix(HOT_PRESSURE_IN)
    m.fs.hex.hot_side.properties_in[t0].temperature.fix(HOT_TEMPERATURE_IN)
    for comp, xf in HOT_MOLE_FRACS.items():
        m.fs.hex.hot_side.properties_in[t0].mole_frac_comp[comp].fix(xf)

    m.fs.hex.cold_side.properties_in[t0].flow_mol.fix(COLD_FLOW_MOL)
    m.fs.hex.cold_side.properties_in[t0].pressure.fix(COLD_PRESSURE_IN)
    m.fs.hex.cold_side.properties_in[t0].temperature.fix(COLD_TEMPERATURE_IN)
    for comp, xf in COLD_MOLE_FRACS.items():
        m.fs.hex.cold_side.properties_in[t0].mole_frac_comp[comp].fix(xf)

    m.fs.hex.hot_side.properties_out[t0].temperature.fix(HOT_TEMPERATURE_OUT)
    m.fs.hex.hot_side.deltaP[t0].fix(HOT_DELTA_P)
    m.fs.hex.cold_side.deltaP[t0].fix(COLD_DELTA_P)
    m.fs.hex.overall_heat_transfer_coefficient[t0].fix(OVERALL_HTC)

    m.fs.hex.initialize(outlvl=idaeslog.INFO)

    solver = get_solver()
    results = solver.solve(m, tee=True)

    print("\n--- IDAES HEAT EXCHANGER RESULTS ---")
    print(f"Heat Duty: {pyo.value(m.fs.hex.heat_duty[t0]) / 1000:.2f} kW")
    print(f"Heat Transfer Area: {pyo.value(m.fs.hex.area):.2f} m^2")

    print("\n--- HOT STREAM ---")
    T_hi = pyo.value(m.fs.hex.hot_side.properties_in[t0].temperature)
    T_ho = pyo.value(m.fs.hex.hot_side.properties_out[t0].temperature)
    P_hi = pyo.value(m.fs.hex.hot_side.properties_in[t0].pressure)
    P_ho = pyo.value(m.fs.hex.hot_side.properties_out[t0].pressure)
    print(f"Inlet Temperature:  {T_hi:.2f} K  ({T_hi - 273.15:.1f} C)")
    print(f"Outlet Temperature: {T_ho:.2f} K  ({T_ho - 273.15:.1f} C)")
    print(f"Inlet Pressure:     {P_hi/1000:.1f} kPa")
    print(f"Outlet Pressure:    {P_ho/1000:.1f} kPa")
    print(f"Pressure Drop:      {-pyo.value(m.fs.hex.hot_side.deltaP[t0])/1000:.1f} kPa")
    print(f"Inlet Density:      {pyo.value(m.fs.hex.hot_side.properties_in[t0].dens_mol):.1f} mol/m3")
    print(f"Inlet Cp:           {pyo.value(m.fs.hex.hot_side.properties_in[t0].cp_mol):.1f} J/mol/K")

    print("\n--- COLD STREAM ---")
    T_ci = pyo.value(m.fs.hex.cold_side.properties_in[t0].temperature)
    T_co = pyo.value(m.fs.hex.cold_side.properties_out[t0].temperature)
    P_ci = pyo.value(m.fs.hex.cold_side.properties_in[t0].pressure)
    P_co = pyo.value(m.fs.hex.cold_side.properties_out[t0].pressure)
    print(f"Inlet Temperature:  {T_ci:.2f} K  ({T_ci - 273.15:.1f} C)")
    print(f"Outlet Temperature: {T_co:.2f} K  ({T_co - 273.15:.1f} C)")
    print(f"Inlet Pressure:     {P_ci/1000:.1f} kPa")
    print(f"Outlet Pressure:    {P_co/1000:.1f} kPa")
    print(f"Pressure Drop:      {-pyo.value(m.fs.hex.cold_side.deltaP[t0])/1000:.1f} kPa")
    print(f"Inlet Density:      {pyo.value(m.fs.hex.cold_side.properties_in[t0].dens_mol):.1f} mol/m3")
    print(f"Inlet Cp:           {pyo.value(m.fs.hex.cold_side.properties_in[t0].cp_mol):.1f} J/mol/K")

    print("\n--- HEAT EXCHANGER ---")
    print(f"LMTD:           {pyo.value(m.fs.hex.delta_temperature[t0]):.2f} K")
    print(f"Overall HTC:    {pyo.value(m.fs.hex.overall_heat_transfer_coefficient[t0]):.1f} W/m2/K")

    C_hot = pyo.value(
        m.fs.hex.hot_side.properties_in[t0].flow_mol
        * m.fs.hex.hot_side.properties_in[t0].cp_mol)
    C_cold = pyo.value(
        m.fs.hex.cold_side.properties_in[t0].flow_mol
        * m.fs.hex.cold_side.properties_in[t0].cp_mol)
    C_min = min(C_hot, C_cold)
    Q_max = C_min * (T_hi - T_ci)
    effectiveness = pyo.value(m.fs.hex.heat_duty[t0]) / Q_max
    print(f"Effectiveness:  {effectiveness:.3f}")
    print("------------------------------------")

    if results.solver.termination_condition == pyo.TerminationCondition.optimal:
        print("\nSolver converged to optimal solution.")
    else:
        print(f"\nSolver terminated: {results.solver.termination_condition}")


if __name__ == '__main__':
    main()
