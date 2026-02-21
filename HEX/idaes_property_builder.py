from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import pyomo.environ as pyo
from pyomo.environ import units as pyunits

from idaes.core import LiquidPhase, VaporPhase
from idaes.models.properties.modular_properties.base.generic_property import Component
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.pure import Perrys, RPP4
from idaes.models.properties.modular_properties.pure.ConstantProperties import Constant

# Local (same-directory) import by design; DWSIM_HEX is not a Python package.
import property_database_interface as pdb


@dataclass(frozen=True)
class BuiltComponent:
    comp_id: str
    display_name: str
    cas: Optional[str]
    component_config: Dict[str, Any]


def _safe_id(s: str) -> str:
    s = (s or "").strip().lower()
    s = s.replace("/", "_").replace(" ", "_").replace("-", "_")
    s = "".join(ch for ch in s if (ch.isalnum() or ch == "_"))
    if not s:
        s = "comp"
    if s[0].isdigit():
        s = f"c_{s}"
    return s


def _value_or_raise(v: Any, msg: str) -> Any:
    if v is None:
        raise ValueError(msg)
    return v


def build_component_from_databases(
    cas_or_name: str,
    *,
    prefer_rpp_ig_cp: bool = True,
    allow_constant_ig_cp_fallback: bool = True,
    allow_constant_liq_cp_fallback: bool = True,
    allow_constant_liq_density_fallback: bool = True,
) -> BuiltComponent:
    """
    Build an IDAES component config entry from local databases.

    Critical properties are REQUIRED for cubic EOS use. If missing, this raises.
    """
    crit = pdb.get_critical_properties(cas_or_name)
    tc = crit.get("Tc_K")
    pc = crit.get("Pc_Pa")
    omega = crit.get("omega")
    mw_kg_per_kmol = crit.get("MW_kg_per_kmol")

    _value_or_raise(tc, f"Missing Tc for {cas_or_name}")
    _value_or_raise(pc, f"Missing Pc for {cas_or_name}")
    _value_or_raise(mw_kg_per_kmol, f"Missing MW for {cas_or_name}")
    if omega is None:
        # Cubic EOS can still run with omega, but treat as required for realism.
        raise ValueError(f"Missing omega for {cas_or_name}")

    display_name = crit.get("name") or cas_or_name
    cas = crit.get("cas")
    comp_id = _safe_id(cas_or_name if cas_or_name else (cas or display_name))

    # --- Liquid density (Perry's eqn_type 1 preferred)
    dens = pdb.get_perrys_density_coefficients(cas_or_name)
    if dens:
        dens_method = Perrys.dens_mol_liq_comp
        dens_param = {
            "eqn_type": 1,
            "1": (dens["1"], pyunits.kmol / pyunits.m**3),
            "2": (dens["2"], pyunits.dimensionless),
            "3": (dens["3"], pyunits.K),
            "4": (dens["4"], pyunits.dimensionless),
        }
    elif allow_constant_liq_density_fallback:
        # Fallback: constant liquid molar density is a blunt approximation.
        # Use 55.5 kmol/m3 as a generic "waterlike" default unless user overrides.
        dens_method = Constant.dens_mol_liq_comp
        dens_param = (55.5, pyunits.kmol / pyunits.m**3)
    else:
        dens_method = None
        dens_param = None

    # --- Liquid Cp (Perry's 5-term polynomial preferred)
    cp_liq = pdb.get_perrys_cp_liq_coefficients(cas_or_name)
    if cp_liq:
        cp_liq_method = Perrys.cp_mol_liq_comp
        cp_liq_param = {
            "1": (cp_liq["1"], pyunits.J / pyunits.kmol / pyunits.K),
            "2": (cp_liq["2"], pyunits.J / pyunits.kmol / pyunits.K**2),
            "3": (cp_liq["3"], pyunits.J / pyunits.kmol / pyunits.K**3),
            "4": (cp_liq["4"], pyunits.J / pyunits.kmol / pyunits.K**4),
            "5": (cp_liq["5"], pyunits.J / pyunits.kmol / pyunits.K**5),
        }
    elif allow_constant_liq_cp_fallback:
        cp_liq_method = Constant.cp_mol_liq_comp
        cp_liq_param = (75e3, pyunits.J / pyunits.kmol / pyunits.K)  # crude default
    else:
        cp_liq_method = None
        cp_liq_param = None

    # --- Ideal gas Cp (RPP4 from ChemSep RPPHeatCapacityCp preferred)
    ig_cp_rpp = pdb.get_rpp_ig_cp_coefficients(cas_or_name) if prefer_rpp_ig_cp else None
    if ig_cp_rpp:
        cp_ig_method = RPP4.cp_mol_ig_comp
        enth_ig_method = RPP4.enth_mol_ig_comp
        entr_ig_method = RPP4.entr_mol_ig_comp

        # ChemSep coefficients are in J/kmol/K^n; RPP4 expects J/mol/K^n.
        cp_ig_param = {
            "A": (ig_cp_rpp["A"] / 1000.0, pyunits.J / pyunits.mol / pyunits.K),
            "B": (ig_cp_rpp["B"] / 1000.0, pyunits.J / pyunits.mol / pyunits.K**2),
            "C": (ig_cp_rpp["C"] / 1000.0, pyunits.J / pyunits.mol / pyunits.K**3),
            "D": (ig_cp_rpp["D"] / 1000.0, pyunits.J / pyunits.mol / pyunits.K**4),
        }
        # Entropy/formation refs: use 0 by default; offsets cancel in HX duties.
        enth_form = 0.0
        entr_form = 0.0
    elif allow_constant_ig_cp_fallback:
        cp_ig_method = Constant.cp_mol_ig_comp
        enth_ig_method = Constant.enth_mol_ig_comp
        entr_ig_method = Constant.entr_mol_ig_comp
        cp_ig_param = (35.0, pyunits.J / pyunits.mol / pyunits.K)
        enth_form = 0.0
        entr_form = 0.0
    else:
        cp_ig_method = None
        enth_ig_method = None
        entr_ig_method = None
        cp_ig_param = None
        enth_form = None
        entr_form = None

    parameter_data: Dict[str, Any] = {
        "mw": (mw_kg_per_kmol / 1000.0, pyunits.kg / pyunits.mol),
        "pressure_crit": (pc, pyunits.Pa),
        "temperature_crit": (tc, pyunits.K),
        "omega": float(omega),
    }

    if dens_method is Perrys.dens_mol_liq_comp and isinstance(dens_param, dict):
        parameter_data["dens_mol_liq_comp_coeff"] = dens_param
    elif dens_method is Constant.dens_mol_liq_comp:
        parameter_data["dens_mol_liq_comp_coeff"] = dens_param

    if cp_liq_method is Perrys.cp_mol_liq_comp and isinstance(cp_liq_param, dict):
        parameter_data["cp_mol_liq_comp_coeff"] = cp_liq_param
        # Required by Perrys.entr_mol_liq_comp (offset is arbitrary for HX work)
        parameter_data["entr_mol_form_liq_comp_ref"] = (0.0, pyunits.J / pyunits.mol / pyunits.K)
        # Only used if include_enthalpy_of_formation=True (we keep False by default)
        parameter_data["enth_mol_form_liq_comp_ref"] = (0.0, pyunits.J / pyunits.mol)
    elif cp_liq_method is Constant.cp_mol_liq_comp:
        parameter_data["cp_mol_liq_comp_coeff"] = cp_liq_param
        parameter_data["entr_mol_form_liq_comp_ref"] = (0.0, pyunits.J / pyunits.mol / pyunits.K)
        parameter_data["enth_mol_form_liq_comp_ref"] = (0.0, pyunits.J / pyunits.mol)

    if cp_ig_method is RPP4.cp_mol_ig_comp and isinstance(cp_ig_param, dict):
        parameter_data["cp_mol_ig_comp_coeff"] = cp_ig_param
        parameter_data["enth_mol_form_vap_comp_ref"] = (enth_form, pyunits.J / pyunits.mol)
        parameter_data["entr_mol_form_vap_comp_ref"] = (entr_form, pyunits.J / pyunits.mol / pyunits.K)
    elif cp_ig_method is Constant.cp_mol_ig_comp:
        parameter_data["cp_mol_ig_comp_coeff"] = cp_ig_param
        parameter_data["enth_mol_form_ig_comp_ref"] = (0.0, pyunits.J / pyunits.mol)
        parameter_data["entr_mol_form_ig_comp_ref"] = (0.0, pyunits.J / pyunits.mol / pyunits.K)

    component_config: Dict[str, Any] = {
        "type": Component,
        "parameter_data": parameter_data,
    }

    # Attach methods (only when available)
    if dens_method is not None:
        component_config["dens_mol_liq_comp"] = dens_method
    if cp_liq_method is not None:
        component_config["cp_mol_liq_comp"] = cp_liq_method
        component_config["enth_mol_liq_comp"] = Perrys.enth_mol_liq_comp if cp_liq_method is Perrys.cp_mol_liq_comp else Constant.enth_mol_liq_comp
        component_config["entr_mol_liq_comp"] = Perrys.entr_mol_liq_comp if cp_liq_method is Perrys.cp_mol_liq_comp else Constant.entr_mol_liq_comp

    if cp_ig_method is not None:
        component_config["cp_mol_ig_comp"] = cp_ig_method
    if enth_ig_method is not None:
        component_config["enth_mol_ig_comp"] = enth_ig_method
    if entr_ig_method is not None:
        component_config["entr_mol_ig_comp"] = entr_ig_method

    return BuiltComponent(comp_id=comp_id, display_name=display_name, cas=cas, component_config=component_config)


def _build_pr_kappa_matrix(components: Sequence[BuiltComponent]) -> Dict[Tuple[str, str], float]:
    """
    Build `PR_kappa` for CubicType.PR using `ipd/pr.ipd` as primary source.
    Default is 0.0 for missing pairs.
    """
    cas_by_id: Dict[str, Optional[str]] = {c.comp_id: c.cas for c in components}
    kij: Dict[Tuple[str, str], float] = {}
    for i in components:
        for j in components:
            cas_i = cas_by_id[i.comp_id]
            cas_j = cas_by_id[j.comp_id]
            if i.comp_id == j.comp_id:
                kij[(i.comp_id, j.comp_id)] = 0.0
                continue
            if cas_i and cas_j:
                k = pdb.get_pr_binary_interaction(cas_i, cas_j)
                if k is None:
                    k = 0.0
            else:
                k = 0.0
            kij[(i.comp_id, j.comp_id)] = float(k)
    return kij


def build_generic_parameter_block_config(
    component_specs: Sequence[str],
    *,
    phases: Sequence[str] = ("Liq",),
    eos: str = "PR",
    prefer_rpp_ig_cp: bool = True,
    include_enthalpy_of_formation: bool = False,
) -> Tuple[Dict[str, Any], List[BuiltComponent]]:
    """
    Build a config dict suitable for:
        GenericParameterBlock(**config_dict)

    Returns:
        (config_dict, built_components)
    """
    built: List[BuiltComponent] = [build_component_from_databases(s, prefer_rpp_ig_cp=prefer_rpp_ig_cp) for s in component_specs]

    phases_cfg: Dict[str, Any] = {}
    for ph in phases:
        if ph == "Liq":
            phases_cfg["Liq"] = {
                "type": LiquidPhase,
                "equation_of_state": Cubic,
                "equation_of_state_options": {"type": CubicType.PR if eos.upper() == "PR" else CubicType.SRK},
            }
        elif ph == "Vap":
            phases_cfg["Vap"] = {
                "type": VaporPhase,
                "equation_of_state": Cubic,
                "equation_of_state_options": {"type": CubicType.PR if eos.upper() == "PR" else CubicType.SRK},
            }
        else:
            raise ValueError(f"Unsupported phase tag: {ph}")

    components_cfg = {c.comp_id: c.component_config for c in built}

    parameter_data: Dict[str, Any] = {}
    if eos.upper() == "PR":
        parameter_data["PR_kappa"] = _build_pr_kappa_matrix(built)

    config_dict: Dict[str, Any] = {
        "components": components_cfg,
        "phases": phases_cfg,
        "state_definition": FTPx,
        "state_bounds": {
            "flow_mol": (0, 10, 1000, pyunits.mol / pyunits.s),
            "temperature": (250, 300, 2000, pyunits.K),
            "pressure": (1e4, 1e5, 5e7, pyunits.Pa),
        },
        "pressure_ref": (101325, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
        "base_units": {
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        "parameter_data": parameter_data,
        "include_enthalpy_of_formation": include_enthalpy_of_formation,
    }
    return config_dict, built


def validate_config_dict(config_dict: Dict[str, Any]) -> None:
    # Minimal structural validation (avoid deep IDAES internals here).
    for k in ("components", "phases", "state_definition", "base_units"):
        if k not in config_dict:
            raise ValueError(f"Invalid property config: missing key '{k}'")
    if not config_dict["components"]:
        raise ValueError("Invalid property config: empty components")

