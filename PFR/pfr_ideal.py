"""Stage-1: Ideal 1-D Plug-Flow Reactor model for ethane pyrolysis.

State vector y = [Y_0 .. Y_{n-1}, T, P]
  dY_i/dz = (1/v_z) * sum_j nu_ij * Omega_j
  dT/dz   = (1/(rho Cp v_z)) * [(-dH_R) sum_j Omega_j + 4 q''/D]
  dP/dz   = -(4f/D) * rho v_z^2 / 2  (Darcy-Weisbach, Blasius f)
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from scipy.integrate import solve_ivp

try:
    from .kinetics_parser import (
        parse_csv_mechanism,
        Reaction,
        resolve_mechanism_path,
        MECHANISM_REDUCED,
        MECHANISM_ALL,
        get_skipped_rows,
    )
    from . import props
    from .props import cp_molar, sensible_enthalpy, gas_viscosity, resolve_species
    from .equilibrium import (
        EquilibriumCalculator,
        build_product_orders,
        net_reaction_rates,
        delta_G_reaction,
        delta_H_reaction_thermo,
    )
    from .eos_ppr78 import (
        evaluate_eos,
        get_Z,
        get_density,
        get_concentration,
        compressibility_factor,
        fugacity_coefficients,
        residual_enthalpy,
        residual_Cp,
        CRITICAL_PROPERTIES,
    )
except ImportError:
    from kinetics_parser import (
        parse_csv_mechanism,
        Reaction,
        resolve_mechanism_path,
        MECHANISM_REDUCED,
        MECHANISM_ALL,
        get_skipped_rows,
    )
    import props
    from props import cp_molar, sensible_enthalpy, gas_viscosity, resolve_species
    from equilibrium import (
        EquilibriumCalculator,
        build_product_orders,
        net_reaction_rates,
        delta_G_reaction,
        delta_H_reaction_thermo,
    )
    from eos_ppr78 import (
        evaluate_eos,
        get_Z,
        get_density,
        get_concentration,
        compressibility_factor,
        fugacity_coefficients,
        residual_enthalpy,
        residual_Cp,
        CRITICAL_PROPERTIES,
    )

ROOT = Path(__file__).resolve().parent

R_GAS = 8.314462618  # J/mol/K
D_TUBE = 0.005       # m
L_TUBE = 15.0        # m
T_SET = 1000.0       # K
P_IN = 2.0 * 101_325.0  # Pa
V_Z0 = 20.0         # m/s
TARGET_CONVERSION = 0.65

DEFAULT_N2_DILUTION = 0.1
FEED_FORMULA = {"C2H6": 1.0 - DEFAULT_N2_DILUTION, "N2": DEFAULT_N2_DILUTION}

MECHANISM_DEFAULT = MECHANISM_REDUCED
_MECHANISM_PATH = resolve_mechanism_path(MECHANISM_DEFAULT)


def set_mechanism_path(path: Path | str | None) -> None:
    global _MECHANISM_PATH, _KINETICS, _MW_ARRAY, _EQUILIBRIUM_CALC
    new_path = resolve_mechanism_path(path)
    if new_path != _MECHANISM_PATH:
        _MECHANISM_PATH = new_path
        _KINETICS = None
        _MW_ARRAY = None
        _EQUILIBRIUM_CALC = None


def get_mechanism_path() -> Path:
    return _MECHANISM_PATH


def _build_kinetics(mechanism_path: Path) -> tuple:
    """Parse mechanism and build numpy arrays for rate evaluation."""
    rxns = parse_csv_mechanism(mechanism_path)
    
    valid_rxns: List[Reaction] = []
    missing_rxns: List[str] = []
    
    for r in rxns:
        missing_species = []
        for sp in r["stoich"]:
            if sp != "M":
                try:
                    resolve_species(sp)
                except KeyError:
                    missing_species.append(sp)
        
        if missing_species:
            missing_rxns.append(f"{r['label']} (missing: {', '.join(missing_species)})")
        else:
            valid_rxns.append(r)
    
    global _SKIPPED_REACTIONS_MISSING
    _SKIPPED_REACTIONS_MISSING = missing_rxns
    
    species_formula: List[str] = []
    for r in valid_rxns:
        for sp in r["stoich"]:
            if sp not in species_formula:
                species_formula.append(sp)
    
    INERT_SPECIES = ["N2"]
    for inert in INERT_SPECIES:
        if inert not in species_formula:
            species_formula.append(inert)
    
    species_keys = [resolve_species(sp) for sp in species_formula]
    spec_idx = {sp: i for i, sp in enumerate(species_formula)}
    
    n = len(species_formula)
    m = len(valid_rxns)
    
    nu = np.zeros((n, m))
    A = np.zeros(m)
    n_exp = np.zeros(m)
    Ea = np.zeros(m)
    reactant_orders: List[Dict[int, float]] = []
    product_orders: List[Dict[int, float]] = []
    has_third_body = np.zeros(m, dtype=bool)
    
    for j, r in enumerate(valid_rxns):
        A[j] = r["A_si"]
        n_exp[j] = r["n"]
        Ea[j] = r["Ea_si"]
        has_third_body[j] = r.get("has_third_body", False)
        
        for sp, coeff in r["stoich"].items():
            nu[species_formula.index(sp), j] = coeff
        
        rord = {}
        for sp, order in r.get("reactants", {}).items():
            rord[species_formula.index(sp)] = order
        reactant_orders.append(rord)
        
        pord = build_product_orders(r["stoich"], spec_idx)
        product_orders.append(pord)
    
    return (species_formula, species_keys, nu, A, n_exp, Ea, 
            reactant_orders, has_third_body, product_orders, valid_rxns)


_KINETICS = None
_SKIPPED_REACTIONS_MISSING: List[str] = []
_EQUILIBRIUM_CALC = None


def _get_kinetics(mechanism_path: Path | str | None = None):
    global _KINETICS
    if mechanism_path is not None:
        set_mechanism_path(mechanism_path)
    if _KINETICS is None:
        _KINETICS = _build_kinetics(_MECHANISM_PATH)
    return _KINETICS


def get_skipped_reactions() -> Dict[str, List[str]]:
    return {
        "empty_kinetics": get_skipped_rows(),
        "missing_species": list(_SKIPPED_REACTIONS_MISSING),
    }


def report_skipped_reactions(context: str | None = None) -> None:
    skipped = get_skipped_reactions()
    empty_rows = skipped["empty_kinetics"]
    missing = skipped["missing_species"]
    total = len(empty_rows) + len(missing)
    prefix = f"[{context}] " if context else ""
    if total == 0:
        print(f"{prefix}Skipped reactions: none")
        return
    print(f"{prefix}Skipped reactions:")
    if empty_rows:
        print(f"  Empty kinetics ({len(empty_rows)}):")
        for item in empty_rows:
            print(f"    - {item}")
    if missing:
        print(f"  Missing species ({len(missing)}):")
        for item in missing:
            print(f"    - {item}")

def get_species_list() -> List[str]:
    return _get_kinetics()[0]

def get_species_index() -> Dict[str, int]:
    sp = get_species_list()
    return {s: i for i, s in enumerate(sp)}


def _get_MW_array() -> np.ndarray:
    """Molecular weights [kg/mol] for all species."""
    _, species_keys, *_ = _get_kinetics()
    db = props._get_db()
    mw = []
    for key in species_keys:
        mw.append(db[key]["identifiers"]["MW_kg_per_mol"])
    return np.array(mw)

_MW_ARRAY = None

def get_MW_array() -> np.ndarray:
    global _MW_ARRAY
    if _MW_ARRAY is None:
        _MW_ARRAY = _get_MW_array()
    return _MW_ARRAY


def mixture_viscosity(Y: np.ndarray, T: float) -> float:
    """Molar-fraction weighted viscosity [Pa s]."""
    _, species_keys, *_ = _get_kinetics()
    mu = np.array([gas_viscosity(T, key) for key in species_keys])
    return float(np.dot(Y, mu))


def mixture_cp(Y: np.ndarray, T: float, P: float = None, use_eos: bool = True) -> float:
    """Molar Cp [J/mol/K] = ideal (NASA7) + residual (PPR78)."""
    species_formula, species_keys, *_ = _get_kinetics()
    cps = np.array([cp_molar(T, key) for key in species_keys])
    Cp_ideal = float(np.dot(Y, cps))
    
    if use_eos and P is not None:
        Cp_res = residual_Cp(species_formula, Y, T, P)
        return Cp_ideal + Cp_res
    
    return Cp_ideal


def mixture_density(Y: np.ndarray, T: float, P: float, use_eos: bool = True) -> float:
    """Mass density [kg/m3] via PPR78 EOS (or ideal gas if use_eos=False)."""
    species_formula, species_keys, *_ = _get_kinetics()
    MW = get_MW_array()
    MW_mix = float(np.dot(Y, MW))
    
    if use_eos:
        Z = compressibility_factor(species_formula, Y, T, P)
    else:
        Z = 1.0
    
    return P * MW_mix / (Z * R_GAS * T)


def mixture_molar_mass(Y: np.ndarray) -> float:
    return float(np.dot(Y, get_MW_array()))


def mixture_concentration(Y: np.ndarray, T: float, P: float, use_eos: bool = True) -> float:
    """Total molar concentration [mol/m3] via PPR78 EOS."""
    if use_eos:
        species_formula, *_ = _get_kinetics()
        Z = compressibility_factor(species_formula, Y, T, P)
    else:
        Z = 1.0
    
    return P / (Z * R_GAS * T)


def reaction_rate_constants(T: float) -> np.ndarray:
    """Forward Arrhenius rate constants k_f(T). Reference T = 298 K."""
    kinetics = _get_kinetics()
    A = kinetics[3]
    n_exp = kinetics[4]
    Ea = kinetics[5]
    return A * (T / 298.0) ** n_exp * np.exp(-Ea / (R_GAS * T))


def _get_equilibrium_calculator() -> EquilibriumCalculator:
    global _EQUILIBRIUM_CALC
    if _EQUILIBRIUM_CALC is None:
        kinetics = _get_kinetics()
        species_formula = kinetics[0]
        species_keys = kinetics[1]
        reactions_raw = kinetics[9]
        _EQUILIBRIUM_CALC = EquilibriumCalculator(species_formula, species_keys, reactions_raw)
    return _EQUILIBRIUM_CALC


def omega_rates(C: np.ndarray, T: float, P: float, 
                include_reverse: bool = True) -> np.ndarray:
    """Net reaction rates Omega_j [mol/m3/s].
    
    Forward: k_f * prod(C_reactants^order)
    Reverse: k_r * prod(C_products^order), k_r = k_f / Kc
    Third-body: multiply by [M] = P/(Z*R*T)
    """
    kinetics = _get_kinetics()
    reactant_orders = kinetics[6]
    has_third_body = kinetics[7]
    product_orders = kinetics[8]
    
    k_forward = reaction_rate_constants(T)
    m = len(k_forward)
    
    if not include_reverse:
        omega = np.zeros(m)
        species_formula, *_ = _get_kinetics()
        C_sum = max(C.sum(), 1e-30)
        Y_approx = C / C_sum
        Z = compressibility_factor(species_formula, Y_approx, T, P)
        C_total = P / (Z * R_GAS * T)
        for j in range(m):
            rate = k_forward[j]
            for i_sp, order in reactant_orders[j].items():
                rate *= max(C[i_sp], 1e-30) ** order
            if has_third_body[j]:
                rate *= C_total
            omega[j] = rate
        return omega
    
    eq_calc = _get_equilibrium_calculator()
    C_sum = max(C.sum(), 1e-30)
    Y_approx = C / C_sum
    k_reverse = eq_calc.compute_reverse_rate_constants(
        T, k_forward, y=Y_approx, P=P, use_fugacity=True
    )
    
    _, _, omega_net = net_reaction_rates(
        C, T, P, k_forward, k_reverse,
        reactant_orders, product_orders, has_third_body
    )
    
    return omega_net


def omega_rates_detailed(C: np.ndarray, T: float, P: float) -> tuple:
    """Return (omega_forward, omega_reverse, omega_net) for analysis."""
    kinetics = _get_kinetics()
    reactant_orders = kinetics[6]
    has_third_body = kinetics[7]
    product_orders = kinetics[8]
    
    k_forward = reaction_rate_constants(T)
    
    C_sum = max(C.sum(), 1e-30)
    Y_approx = C / C_sum
    
    eq_calc = _get_equilibrium_calculator()
    k_reverse = eq_calc.compute_reverse_rate_constants(
        T, k_forward, y=Y_approx, P=P, use_fugacity=True
    )
    
    return net_reaction_rates(
        C, T, P, k_forward, k_reverse,
        reactant_orders, product_orders, has_third_body
    )


def delta_H_reaction(T: float) -> np.ndarray:
    """dH_j(T) for each reaction [J/mol] from NASA7 enthalpies."""
    kinetics = _get_kinetics()
    species_formula = kinetics[0]
    species_keys = kinetics[1]
    reactions_raw = kinetics[9]
    
    formula_to_key = {f: k for f, k in zip(species_formula, species_keys)}
    
    dH = np.zeros(len(reactions_raw))
    for j, rxn in enumerate(reactions_raw):
        stoich = rxn.get('stoich', {})
        try:
            dH[j] = delta_H_reaction_thermo(T, stoich, formula_to_key)
        except (KeyError, ValueError):
            h = np.array([sensible_enthalpy(T, key) for key in species_keys])
            nu = kinetics[2]
            dH[j] = np.dot(nu[:, j], h)
    
    return dH


def estimate_heat_flux(T: float = T_SET, P: float = P_IN,
                       Y0: np.ndarray | None = None,
                       target_conversion: float = TARGET_CONVERSION,
                       mechanism_path: Path | str | None = None) -> float:
    """Estimate constant q'' [W/m2] for near-isothermal operation."""
    species_formula, species_keys, nu, *_ = _get_kinetics(mechanism_path)
    spec_idx = get_species_index()
    n_spec = len(species_formula)
    
    if Y0 is None:
        Y0 = np.zeros(n_spec)
        for sp, frac in FEED_FORMULA.items():
            Y0[spec_idx[sp]] = frac
        Y0 /= Y0.sum()
    
    dH_main = 137e3  # J/mol, C2H6 -> C2H4 + H2
    
    C_total = mixture_concentration(Y0, T, P, use_eos=True)
    A_cross = math.pi * (D_TUBE / 2) ** 2
    F_in = C_total * V_Z0 * A_cross
    
    Q_total = F_in * target_conversion * dH_main
    A_wall = math.pi * D_TUBE * L_TUBE
    
    return Q_total / A_wall


class PFRState:
    def __init__(self, T_in: float = T_SET, P_in: float = P_IN,
                 v_z0: float = V_Z0, D: float = D_TUBE,
                 q_flux: float | None = None,
                 mechanism_path: Path | str | None = None,
                 rho_in: float | None = None):
        self.T_in = T_in
        self.P_in = P_in
        self.v_z0 = v_z0
        self.D = D
        self.mechanism_path = mechanism_path
        self.q_flux = q_flux if q_flux is not None else estimate_heat_flux(
            T_in, P_in, mechanism_path=mechanism_path
        )
        self.rho_in = rho_in
        self.mu_in: float | None = None


def make_diluted_feed(hydrocarbon: str = "C2H6", n2_fraction: float = 0.5) -> Dict[str, float]:
    if not 0.0 <= n2_fraction <= 1.0:
        raise ValueError(f"n2_fraction must be between 0 and 1, got {n2_fraction}")
    
    feed = {hydrocarbon: 1.0 - n2_fraction}
    if n2_fraction > 0:
        feed["N2"] = n2_fraction
    return feed


def ode_pfr(z: float, y: np.ndarray, state: PFRState) -> np.ndarray:
    """RHS of PFR ODE. State vector: y = [Y_0..Y_{n-1}, T, P]."""
    species_formula, _, nu, *_ = _get_kinetics(state.mechanism_path)
    n_spec = len(species_formula)
    
    Y = y[:n_spec].copy()
    Y = np.clip(Y, 1e-15, 1.0)
    Y /= Y.sum()
    T = max(y[n_spec], 300.0)
    P = max(y[n_spec + 1], 1e3)
    
    mu_mix = mixture_viscosity(Y, T)
    if state.mu_in is None:
        state.mu_in = mu_mix
    
    rho = mixture_density(Y, T, P, use_eos=True)
    Cp_mix = mixture_cp(Y, T, P, use_eos=True)
    MW_mix = mixture_molar_mass(Y)
    
    if state.rho_in is None:
        state.rho_in = rho
    v_z = state.v_z0 * (state.rho_in / rho)
    v_z = max(v_z, 0.01)
    
    C_total = mixture_concentration(Y, T, P, use_eos=True)
    C = Y * C_total
    
    omega = omega_rates(C, T, P)
    
    dYdz = (nu @ omega) / (v_z * C_total)
    
    dH = delta_H_reaction(T)
    Q_rxn = -np.dot(dH, omega)
    Cp_mass = Cp_mix / MW_mix
    dTdz = (1.0 / (rho * Cp_mass * v_z)) * (Q_rxn + 4.0 * state.q_flux / state.D)
    
    Re = rho * v_z * state.D / mu_mix
    Re = max(Re, 100.0)
    f = 0.0791 * Re ** (-0.25)  # Blasius
    dPdz = -4.0 * f / state.D * 0.5 * rho * v_z ** 2
    
    return np.concatenate([dYdz, [dTdz, dPdz]])


def run_pfr(L: float = L_TUBE, T_in: float = T_SET, P_in: float = P_IN,
            v_z0: float = V_Z0, D: float = D_TUBE,
            feed: Dict[str, float] | None = None,
            q_flux: float | None = None,
            max_step: float = 0.01,
            mechanism_path: Path | str | None = None) -> Any:
    """Integrate PFR from z=0 to z=L. Returns scipy OdeSolution."""
    if mechanism_path is not None:
        set_mechanism_path(mechanism_path)
    species_formula = get_species_list()
    spec_idx = get_species_index()
    report_skipped_reactions("run_pfr")
    n_spec = len(species_formula)
    
    Y0 = np.zeros(n_spec)
    feed = feed or FEED_FORMULA
    for sp, frac in feed.items():
        if sp in spec_idx:
            Y0[spec_idx[sp]] = frac
        else:
            raise ValueError(f"Unknown species in feed: {sp}")
    Y0 /= Y0.sum()
    
    y0 = np.concatenate([Y0, [T_in, P_in]])
    
    rho_in = mixture_density(Y0, T_in, P_in, use_eos=True)
    state = PFRState(
        T_in=T_in,
        P_in=P_in,
        v_z0=v_z0,
        D=D,
        q_flux=q_flux,
        mechanism_path=mechanism_path,
        rho_in=rho_in,
    )
    
    sol = solve_ivp(
        lambda z, y: ode_pfr(z, y, state),
        (0.0, L),
        y0,
        method="BDF",
        atol=1e-10,
        rtol=1e-8,
        max_step=max_step,
    )
    
    sol.species = species_formula
    sol.n_spec = n_spec
    sol.state = state
    
    return sol


def print_results(sol) -> None:
    spec_idx = get_species_index()
    n_spec = sol.n_spec
    
    print(f"Integration success: {sol.success}")
    print(f"Message: {sol.message}")
    print(f"Number of z-steps: {len(sol.t)}")
    print(f"z range: {sol.t[0]:.3f} to {sol.t[-1]:.3f} m")
    print()
    
    Y_out = sol.y[:n_spec, -1].copy()
    Y_out = np.clip(Y_out, 0, None)
    Y_out /= Y_out.sum()
    
    Y_in = sol.y[:n_spec, 0].copy()
    Y_in = np.clip(Y_in, 0, None)
    Y_in /= Y_in.sum()
    
    print("Inlet -> Outlet:")
    print(f"  T: {sol.y[n_spec, 0]:.1f} -> {sol.y[n_spec, -1]:.1f} K")
    print(f"  P: {sol.y[n_spec+1, 0]/1e5:.4f} -> {sol.y[n_spec+1, -1]/1e5:.4f} bar")
    dP = (sol.y[n_spec+1, 0] - sol.y[n_spec+1, -1]) / 1e5
    print(f"  dP: {dP*1000:.2f} mbar")
    print()
    
    print("Mole fractions (inlet -> outlet):")
    for sp in ["C2H6", "C2H4", "C2H2", "C2H5", "C2H3", "C1H4", "C1H3", "H2", "H", "C3H6", "C3H7", "N2"]:
        if sp in spec_idx:
            i = spec_idx[sp]
            if Y_out[i] > 1e-6 or Y_in[i] > 1e-6:
                print(f"  {sp:6s}: {Y_in[i]:.6f} -> {Y_out[i]:.6f}")
    
    print(f"\n  Sum Y: {Y_out.sum():.6f}")
    
    if "C2H6" in spec_idx:
        i_eth = spec_idx["C2H6"]
        conv = 1.0 - Y_out[i_eth] / Y_in[i_eth]
        print(f"\nEthane conversion: {conv*100:.1f}%")
        
        if "C2H4" in spec_idx and conv > 0.01:
            i_ene = spec_idx["C2H4"]
            ethylene_formed = Y_out[i_ene] - Y_in[i_ene]
            ethane_reacted = Y_in[i_eth] - Y_out[i_eth]
            selectivity = ethylene_formed / ethane_reacted if ethane_reacted > 1e-10 else 0
            print(f"Ethylene selectivity: {selectivity*100:.1f}%")
        
        if "C2H2" in spec_idx and conv > 0.01:
            i_ace = spec_idx["C2H2"]
            acetylene_formed = Y_out[i_ace] - Y_in[i_ace]
            ethane_reacted = Y_in[i_eth] - Y_out[i_eth]
            selectivity_ace = acetylene_formed / ethane_reacted if ethane_reacted > 1e-10 else 0
            print(f"Acetylene selectivity: {selectivity_ace*100:.1f}%")
    
    if hasattr(sol, 'state'):
        print(f"Wall heat flux: {sol.state.q_flux/1e3:.2f} kW/m2")


if __name__ == "__main__":
    print("Stage-1 Ideal PFR - Ethane Pyrolysis with N2 Dilution")
    print("=" * 60)
    print(f"Tube: D={D_TUBE*1000:.0f} mm, L={L_TUBE:.1f} m")
    print(f"Inlet: T={T_SET:.0f} K, P={P_IN/1e5:.1f} bar, v_z={V_Z0:.1f} m/s")
    print(f"N2 dilution: {DEFAULT_N2_DILUTION*100:.0f}% (mole basis)")
    print()
    
    sol = run_pfr()
    print_results(sol)
