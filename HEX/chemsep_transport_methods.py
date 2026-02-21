from __future__ import annotations

import pyomo.environ as pyo
from pyomo.environ import Var, Param, units as pyunits, exp, log

from idaes.core.util.misc import set_param_from_config


class ChemSepLiquidViscosity:
    """
    Pure-component liquid viscosity correlation from ChemSep XML.

    Implemented eqno:
    - 101: ln(mu) = A + B/T + C*ln(T) + D*T^E

    Notes:
    - This is treated as a *pure-component* correlation and is used via
      the generic property package hook `visc_d_phase_comp`.
    - Units are handled by non-dimensionalizing T inside the exponent.
    """

    class visc_d_phase_comp:
        @staticmethod
        def build_parameters(cobj):
            cobj.visc_d_phase_comp_eqno = Param(
                mutable=True, initialize=101, doc="ChemSep eqno for liquid viscosity"
            )
            set_param_from_config(cobj, param="visc_d_phase_comp_coeff", index="eqno")

            cobj.visc_d_phase_comp_coeff_A = Var(
                initialize=0.0, units=pyunits.dimensionless
            )
            set_param_from_config(cobj, param="visc_d_phase_comp_coeff", index="A")

            cobj.visc_d_phase_comp_coeff_B = Var(
                initialize=0.0, units=pyunits.K
            )
            set_param_from_config(cobj, param="visc_d_phase_comp_coeff", index="B")

            cobj.visc_d_phase_comp_coeff_C = Var(
                initialize=0.0, units=pyunits.dimensionless
            )
            set_param_from_config(cobj, param="visc_d_phase_comp_coeff", index="C")

            cobj.visc_d_phase_comp_coeff_D = Var(
                initialize=0.0, units=pyunits.dimensionless
            )
            set_param_from_config(cobj, param="visc_d_phase_comp_coeff", index="D")

            cobj.visc_d_phase_comp_coeff_E = Var(
                initialize=2.0, units=pyunits.dimensionless
            )
            set_param_from_config(cobj, param="visc_d_phase_comp_coeff", index="E")

        @staticmethod
        def return_expression(b, cobj, p, T):  # pylint: disable=unused-argument
            # Only liquid phase is supported here.
            if str(p) != "Liq":
                return pyo.Expression.Skip

            # Ensure parameters exist
            if not hasattr(cobj, "visc_d_phase_comp_eqno"):
                ChemSepLiquidViscosity.visc_d_phase_comp.build_parameters(cobj)

            if int(pyo.value(cobj.visc_d_phase_comp_eqno)) != 101:
                raise ValueError(
                    f"Unsupported ChemSep liquid viscosity eqno={pyo.value(cobj.visc_d_phase_comp_eqno)}"
                )

            Tk = pyunits.convert(T, to_units=pyunits.K)
            t = Tk / pyunits.K  # dimensionless

            expo = (
                cobj.visc_d_phase_comp_coeff_A
                + cobj.visc_d_phase_comp_coeff_B / Tk
                + cobj.visc_d_phase_comp_coeff_C * log(t)
                + cobj.visc_d_phase_comp_coeff_D * t ** cobj.visc_d_phase_comp_coeff_E
            )

            mu = exp(expo) * pyunits.Pa * pyunits.s
            units = b.params.get_metadata().derived_units
            return pyunits.convert(mu, units.VISCOSITY)


class ChemSepLiquidThermalConductivity:
    """
    Pure-component liquid thermal conductivity correlation from ChemSep XML.

    Implemented eqno:
    - 16: k = A + B*T + C*T^2 + D*T^3 + E*T^4
    """

    class therm_cond_phase_comp:
        @staticmethod
        def build_parameters(cobj):
            cobj.therm_cond_phase_comp_eqno = Param(
                mutable=True, initialize=16, doc="ChemSep eqno for liquid thermal conductivity"
            )
            set_param_from_config(cobj, param="therm_cond_phase_comp_coeff", index="eqno")

            cobj.therm_cond_phase_comp_coeff_A = Var(
                initialize=0.1, units=pyunits.W / pyunits.m / pyunits.K
            )
            set_param_from_config(cobj, param="therm_cond_phase_comp_coeff", index="A")

            cobj.therm_cond_phase_comp_coeff_B = Var(
                initialize=0.0, units=pyunits.W / pyunits.m / pyunits.K**2
            )
            set_param_from_config(cobj, param="therm_cond_phase_comp_coeff", index="B")

            cobj.therm_cond_phase_comp_coeff_C = Var(
                initialize=0.0, units=pyunits.W / pyunits.m / pyunits.K**3
            )
            set_param_from_config(cobj, param="therm_cond_phase_comp_coeff", index="C")

            cobj.therm_cond_phase_comp_coeff_D = Var(
                initialize=0.0, units=pyunits.W / pyunits.m / pyunits.K**4
            )
            set_param_from_config(cobj, param="therm_cond_phase_comp_coeff", index="D")

            cobj.therm_cond_phase_comp_coeff_E = Var(
                initialize=0.0, units=pyunits.W / pyunits.m / pyunits.K**5
            )
            set_param_from_config(cobj, param="therm_cond_phase_comp_coeff", index="E")

        @staticmethod
        def return_expression(b, cobj, p, T):  # pylint: disable=unused-argument
            if str(p) != "Liq":
                return pyo.Expression.Skip

            if not hasattr(cobj, "therm_cond_phase_comp_eqno"):
                ChemSepLiquidThermalConductivity.therm_cond_phase_comp.build_parameters(cobj)

            if int(pyo.value(cobj.therm_cond_phase_comp_eqno)) != 16:
                raise ValueError(
                    f"Unsupported ChemSep liquid thermal conductivity eqno={pyo.value(cobj.therm_cond_phase_comp_eqno)}"
                )

            Tk = pyunits.convert(T, to_units=pyunits.K)
            k = (
                cobj.therm_cond_phase_comp_coeff_A
                + cobj.therm_cond_phase_comp_coeff_B * Tk
                + cobj.therm_cond_phase_comp_coeff_C * Tk**2
                + cobj.therm_cond_phase_comp_coeff_D * Tk**3
                + cobj.therm_cond_phase_comp_coeff_E * Tk**4
            )

            units = b.params.get_metadata().derived_units
            return pyunits.convert(k, units.THERMAL_CONDUCTIVITY)

