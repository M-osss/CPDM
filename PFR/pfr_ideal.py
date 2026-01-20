"""Stage-1: Ideal 1-D Plug-Flow Reactor model for ethane pyrolysis.

State vector y = [Y_0 .. Y_{n-1}, T, P]
  Y_i   : molar fraction of species i (sum=1)
  T     : temperature [K]
  P     : pressure [Pa]

ODEs (no inflow terms):
  dY_i/dz = (1/v_z) * Σ_j ν_ij * Ω_j
  dT/dz   = (1/(ρ C_p v_z)) * [ (-ΔH_R) Σ_j Ω_j + 4 q''/D ]
  dP/dz   = − (4 f / D) * ρ v_z² / 2    (Darcy–Weisbach, f via Blasius)

Uses kinetics_parser & props helpers.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from scipy.integrate import solve_ivp

# Handle imports for both module and standalone execution
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

# -----------------------------------------------------------------------------
# Physical constants & reactor geometry
# -----------------------------------------------------------------------------
R_GAS = 8.314462618  # J/mol/K
D_TUBE = 0.005        # m tube diameter
L_TUBE = 15.0         # m tube length

# Operating conditions (Stage-1 defaults)
T_SET = 1000.0       # K target temperature
P_IN = 2.0 * 101_325.0  # Pa (2 bar)
V_Z0 = 20.0           # m/s feed superficial velocity
TARGET_CONVERSION = 0.65  # target ethane conversion for heat flux estimation

# Feed composition: ethane with N2 dilution (default 50% N2 by mole)
# N2 is inert - does not participate in reactions but affects:
# - Heat capacity (thermal ballast)
# - Partial pressures (shifts equilibrium)
# - Density and velocity
DEFAULT_N2_DILUTION = 0.1  # mole fraction of N2 in feed
FEED_FORMULA = {"C2H6": 1.0 - DEFAULT_N2_DILUTION, "N2": DEFAULT_N2_DILUTION}

# -----------------------------------------------------------------------------
# Build kinetics data structures at module load
# -----------------------------------------------------------------------------

MECHANISM_DEFAULT = MECHANISM_REDUCED

# Active mechanism path (mutable)
_MECHANISM_PATH = resolve_mechanism_path(MECHANISM_DEFAULT)


def set_mechanism_path(path: Path | str | None) -> None:
    """Set active mechanism CSV path for kinetics evaluation."""
    global _MECHANISM_PATH, _KINETICS, _MW_ARRAY, _EQUILIBRIUM_CALC
    new_path = resolve_mechanism_path(path)
    if new_path != _MECHANISM_PATH:
        _MECHANISM_PATH = new_path
        _KINETICS = None
        _MW_ARRAY = None
        _EQUILIBRIUM_CALC = None


def get_mechanism_path() -> Path:
    """Return active mechanism CSV path."""
    return _MECHANISM_PATH


def _build_kinetics(mechanism_path: Path) -> tuple:
    """Parse mechanism and build numpy arrays for fast rate evaluation.
    
    Returns:
        species_formula: list of formula tokens (e.g. 'C2H6')
        species_keys: list of resolved db keys (e.g. 'C2H6::ethane::74-84-0')
        nu: stoichiometry matrix (n_species, n_reactions)
        A: pre-exponential factors array
        n_exp: temperature exponents array
        Ea: activation energies array (J/mol)
        reactant_orders: list of dicts {species_idx: order} for each reaction
        has_third_body: boolean array
        product_orders: list of dicts {species_idx: order} for products (reverse rxn)
        reactions_raw: list of raw reaction dicts for equilibrium calculations
    """
    rxns = parse_csv_mechanism(mechanism_path)
    
    # Drop reactions that reference species missing from the thermo database
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
    
    # Collect all species from mechanism
    species_formula: List[str] = []
    for r in valid_rxns:
        for sp in r["stoich"]:
            if sp not in species_formula:
                species_formula.append(sp)
    
    # Add inert species (N2 for dilution) - these have zero stoichiometry
    # but are tracked for mass/energy balance and property calculations
    INERT_SPECIES = ["N2"]
    for inert in INERT_SPECIES:
        if inert not in species_formula:
            species_formula.append(inert)
    
    # Resolve to database keys
    species_keys = [resolve_species(sp) for sp in species_formula]
    
    # Build species index map
    spec_idx = {sp: i for i, sp in enumerate(species_formula)}
    
    n = len(species_formula)
    m = len(valid_rxns)
    
    nu = np.zeros((n, m))
    A = np.zeros(m)
    n_exp = np.zeros(m)
    Ea = np.zeros(m)
    reactant_orders: List[Dict[int, float]] = []
    product_orders: List[Dict[int, float]] = []  # For reverse reactions
    has_third_body = np.zeros(m, dtype=bool)
    
    for j, r in enumerate(valid_rxns):
        A[j] = r["A_si"]
        n_exp[j] = r["n"]
        Ea[j] = r["Ea_si"]
        has_third_body[j] = r.get("has_third_body", False)
        
        for sp, coeff in r["stoich"].items():
            nu[species_formula.index(sp), j] = coeff
        
        # Build reactant order dict for this reaction (negative stoich = reactants)
        rord = {}
        for sp, order in r.get("reactants", {}).items():
            rord[species_formula.index(sp)] = order
        reactant_orders.append(rord)
        
        # Build product order dict for reverse reaction (positive stoich = products)
        pord = build_product_orders(r["stoich"], spec_idx)
        product_orders.append(pord)
    
    # Note: N2 (and other inerts) have nu[i_N2, :] = 0 (no reaction)
    return (species_formula, species_keys, nu, A, n_exp, Ea, 
            reactant_orders, has_third_body, product_orders, valid_rxns)

# Module-level cached kinetics
_KINETICS = None
_SKIPPED_REACTIONS_MISSING: List[str] = []
_EQUILIBRIUM_CALC = None  # EquilibriumCalculator instance

def _get_kinetics(mechanism_path: Path | str | None = None):
    global _KINETICS
    if mechanism_path is not None:
        set_mechanism_path(mechanism_path)
    if _KINETICS is None:
        _KINETICS = _build_kinetics(_MECHANISM_PATH)
    return _KINETICS


def get_skipped_reactions() -> Dict[str, List[str]]:
    """Return skipped reactions grouped by reason."""
    return {
        "empty_kinetics": get_skipped_rows(),
        "missing_species": list(_SKIPPED_REACTIONS_MISSING),
    }


def report_skipped_reactions(context: str | None = None) -> None:
    """Print skipped reactions for the current mechanism."""
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
    """Return list of species formula tokens."""
    return _get_kinetics()[0]

def get_species_index() -> Dict[str, int]:
    """Return mapping formula -> index."""
    sp = get_species_list()
    return {s: i for i, s in enumerate(sp)}

# -----------------------------------------------------------------------------
# Molecular weights lookup
# -----------------------------------------------------------------------------

def _get_MW_array() -> np.ndarray:
    """Return molecular weights in kg/mol for all species."""
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

# -----------------------------------------------------------------------------
# Mixture properties using resolved species keys
# -----------------------------------------------------------------------------

def mixture_viscosity(Y: np.ndarray, T: float) -> float:
    """Molar-fraction weighted viscosity [Pa·s]."""
    _, species_keys, *_ = _get_kinetics()
    mu = np.array([gas_viscosity(T, key) for key in species_keys])
    return float(np.dot(Y, mu))


def mixture_cp(Y: np.ndarray, T: float, P: float = None, use_eos: bool = True) -> float:
    """Molar-fraction weighted heat capacity [J/mol/K].
    
    Cp_total = Cp_ideal + Cp_residual
    
    where Cp_ideal is from NASA7 polynomials and Cp_residual is from PPR78 EOS.
    
    Args:
        Y: molar fractions
        T: temperature [K]
        P: pressure [Pa] (required if use_eos=True)
        use_eos: if True, include residual Cp from EOS
    """
    species_formula, species_keys, *_ = _get_kinetics()
    cps = np.array([cp_molar(T, key) for key in species_keys])
    Cp_ideal = float(np.dot(Y, cps))
    
    if use_eos and P is not None:
        Cp_res = residual_Cp(species_formula, Y, T, P)
        return Cp_ideal + Cp_res
    
    return Cp_ideal


def mixture_density(Y: np.ndarray, T: float, P: float, use_eos: bool = True) -> float:
    """Mixture mass density [kg/m³] using PPR78 EOS or ideal gas law.
    
    rho = P * MW_mix / (Z * R * T)
    where Z is the compressibility factor from PPR78 EOS.
    
    Args:
        Y: molar fractions
        T: temperature [K]
        P: pressure [Pa]
        use_eos: if True, use PPR78 EOS; if False, use ideal gas (Z=1)
    """
    species_formula, species_keys, *_ = _get_kinetics()
    MW = get_MW_array()
    MW_mix = float(np.dot(Y, MW))  # kg/mol
    
    if use_eos:
        Z = compressibility_factor(species_formula, Y, T, P)
    else:
        Z = 1.0
    
    return P * MW_mix / (Z * R_GAS * T)


def mixture_molar_mass(Y: np.ndarray) -> float:
    """Average molar mass [kg/mol]."""
    return float(np.dot(Y, get_MW_array()))


def mixture_concentration(Y: np.ndarray, T: float, P: float, use_eos: bool = True) -> float:
    """Total molar concentration [mol/m³] using PPR78 EOS or ideal gas law.
    
    C_total = P / (Z * R * T)
    where Z is the compressibility factor from PPR78 EOS.
    
    Args:
        Y: molar fractions
        T: temperature [K]
        P: pressure [Pa]
        use_eos: if True, use PPR78 EOS; if False, use ideal gas (Z=1)
    """
    if use_eos:
        species_formula, *_ = _get_kinetics()
        Z = compressibility_factor(species_formula, Y, T, P)
    else:
        Z = 1.0
    
    return P / (Z * R_GAS * T)

# -----------------------------------------------------------------------------
# Reaction rate calculations
# -----------------------------------------------------------------------------

def reaction_rate_constants(T: float) -> np.ndarray:
    """Forward Arrhenius rate constants k_f(T).
    
    k = A * (T/298)^n * exp(-Ea/(R*T))
    
    Units depend on reaction order (mechanism uses molecule-based units,
    converted at parse time).
    Reference temperature: 298 K (standard form).
    """
    kinetics = _get_kinetics()
    A = kinetics[3]
    n_exp = kinetics[4]
    Ea = kinetics[5]
    return A * (T / 298.0) ** n_exp * np.exp(-Ea / (R_GAS * T))


def _get_equilibrium_calculator() -> EquilibriumCalculator:
    """Get or create equilibrium calculator for current mechanism."""
    global _EQUILIBRIUM_CALC
    if _EQUILIBRIUM_CALC is None:
        kinetics = _get_kinetics()
        species_formula = kinetics[0]
        species_keys = kinetics[1]
        reactions_raw = kinetics[9]  # Raw reaction dicts
        _EQUILIBRIUM_CALC = EquilibriumCalculator(species_formula, species_keys, reactions_raw)
    return _EQUILIBRIUM_CALC


def omega_rates(C: np.ndarray, T: float, P: float, 
                include_reverse: bool = True) -> np.ndarray:
    """Compute net reaction rates Ω_j [mol/m³/s].
    
    For each reaction:
        Ω_net = Ω_forward - Ω_reverse
    
    where:
        Ω_forward = k_f * Π[C_reactants]^order
        Ω_reverse = k_r * Π[C_products]^order
        k_r = k_f / Kc (detailed balance)
    
    For third-body reactions: multiply by [M] = P/(R*T).
    
    Args:
        C: species concentrations [mol/m³]
        T: temperature [K]
        P: pressure [Pa]
        include_reverse: if True, compute net rates; if False, forward only
    
    Returns:
        omega: array of net reaction rates (n_reactions,)
    """
    kinetics = _get_kinetics()
    reactant_orders = kinetics[6]
    has_third_body = kinetics[7]
    product_orders = kinetics[8]
    
    k_forward = reaction_rate_constants(T)
    m = len(k_forward)
    
    if not include_reverse:
        # Original forward-only calculation with EOS correction
        omega = np.zeros(m)
        species_formula, *_ = _get_kinetics()
        # Use EOS for total concentration [M] in third-body reactions
        # Y is approximated from C for this purpose
        C_sum = max(C.sum(), 1e-30)
        Y_approx = C / C_sum
        Z = compressibility_factor(species_formula, Y_approx, T, P)
        C_total = P / (Z * R_GAS * T)  # Non-ideal concentration
        for j in range(m):
            rate = k_forward[j]
            for i_sp, order in reactant_orders[j].items():
                rate *= max(C[i_sp], 1e-30) ** order
            if has_third_body[j]:
                rate *= C_total
            omega[j] = rate
        return omega
    
    # Compute reverse rate constants from equilibrium with fugacity corrections
    eq_calc = _get_equilibrium_calculator()
    # Get molar fractions from concentrations for fugacity correction
    C_sum = max(C.sum(), 1e-30)
    Y_approx = C / C_sum
    k_reverse = eq_calc.compute_reverse_rate_constants(
        T, k_forward, y=Y_approx, P=P, use_fugacity=True
    )
    
    # Compute net rates using the equilibrium module
    _, _, omega_net = net_reaction_rates(
        C, T, P, k_forward, k_reverse,
        reactant_orders, product_orders, has_third_body
    )
    
    return omega_net


def omega_rates_detailed(C: np.ndarray, T: float, P: float) -> tuple:
    """Compute forward, reverse, and net reaction rates.
    
    Useful for analysis and debugging.
    
    Args:
        C: species concentrations [mol/m³]
        T: temperature [K]
        P: pressure [Pa]
    
    Returns:
        omega_forward: forward rates [mol/m³/s]
        omega_reverse: reverse rates [mol/m³/s]
        omega_net: net rates [mol/m³/s]
    """
    kinetics = _get_kinetics()
    reactant_orders = kinetics[6]
    has_third_body = kinetics[7]
    product_orders = kinetics[8]
    
    k_forward = reaction_rate_constants(T)
    
    # Get molar fractions from concentrations for fugacity correction
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
    """Compute ΔH_j(T) for each reaction [J/mol].
    
    Uses thermodynamically consistent calculation:
    ΔH_j = Σ_i ν_ij * H_i(T)  (products - reactants)
    
    This uses the NASA7 polynomial enthalpy which includes
    the standard enthalpy of formation.
    """
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
            # Fallback to matrix calculation
            h = np.array([sensible_enthalpy(T, key) for key in species_keys])
            nu = kinetics[2]
            dH[j] = np.dot(nu[:, j], h)
    
    return dH

# -----------------------------------------------------------------------------
# Heat flux estimation (Stage-1: constant q'' to maintain near-isothermal)
# -----------------------------------------------------------------------------

def estimate_heat_flux(T: float = T_SET, P: float = P_IN,
                       Y0: np.ndarray | None = None,
                       target_conversion: float = TARGET_CONVERSION,
                       mechanism_path: Path | str | None = None) -> float:
    """Estimate constant heat flux q'' [W/m²] to maintain near-isothermal operation.
    
    Uses an iterative approach: run short integration segments to estimate
    average heat demand along the reactor, then compute required q''.
    
    For Stage-1, we use a simplified estimate based on expected conversion
    and reaction enthalpy.
    """
    species_formula, species_keys, nu, *_ = _get_kinetics(mechanism_path)
    spec_idx = get_species_index()
    n_spec = len(species_formula)
    
    if Y0 is None:
        Y0 = np.zeros(n_spec)
        for sp, frac in FEED_FORMULA.items():
            Y0[spec_idx[sp]] = frac
        Y0 /= Y0.sum()
    
    # Estimate based on main reaction C2H6 -> C2H4 + H2
    # ΔH_rxn ≈ +137 kJ/mol (endothermic)
    dH_main = 137e3  # J/mol
    
    # Molar flow rate estimate [mol/s] using EOS
    C_total = mixture_concentration(Y0, T, P, use_eos=True)  # mol/m³ (non-ideal)
    A_cross = math.pi * (D_TUBE / 2) ** 2  # m²
    F_in = C_total * V_Z0 * A_cross  # mol/s
    
    # Heat required for target conversion [W]
    Q_total = F_in * target_conversion * dH_main
    
    # Wall area [m²]
    A_wall = math.pi * D_TUBE * L_TUBE
    
    # Heat flux [W/m²]
    q_flux = Q_total / A_wall
    
    return q_flux

# -----------------------------------------------------------------------------
# ODE system for ideal PFR
# -----------------------------------------------------------------------------

class PFRState:
    """Container for PFR operating parameters."""
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
        self.mu_in: float | None = None  # set after first call


def make_diluted_feed(hydrocarbon: str = "C2H6", n2_fraction: float = 0.5) -> Dict[str, float]:
    """Create a feed composition with N2 dilution.
    
    Args:
        hydrocarbon: primary hydrocarbon species (default "C2H6" for ethane)
        n2_fraction: mole fraction of N2 diluent (0 to 1)
    
    Returns:
        dict of {species: mole_fraction}
    
    Examples:
        make_diluted_feed("C2H6", 0.5)  -> {"C2H6": 0.5, "N2": 0.5}
        make_diluted_feed("C2H6", 0.0)  -> {"C2H6": 1.0}  (pure feed)
        make_diluted_feed("C2H6", 0.8)  -> {"C2H6": 0.2, "N2": 0.8}  (highly diluted)
    """
    if not 0.0 <= n2_fraction <= 1.0:
        raise ValueError(f"n2_fraction must be between 0 and 1, got {n2_fraction}")
    
    feed = {hydrocarbon: 1.0 - n2_fraction}
    if n2_fraction > 0:
        feed["N2"] = n2_fraction
    return feed


def ode_pfr(z: float, y: np.ndarray, state: PFRState) -> np.ndarray:
    """Right-hand side of PFR ODE system.
    
    State vector: y = [Y_0, ..., Y_{n-1}, T, P]
    
    Returns: dy/dz
    """
    species_formula, _, nu, *_ = _get_kinetics(state.mechanism_path)
    n_spec = len(species_formula)
    
    # Unpack state
    Y = y[:n_spec].copy()
    Y = np.clip(Y, 1e-15, 1.0)
    Y /= Y.sum()  # renormalize
    T = max(y[n_spec], 300.0)
    P = max(y[n_spec + 1], 1e3)
    
    # Mixture properties
    mu_mix = mixture_viscosity(Y, T)
    if state.mu_in is None:
        state.mu_in = mu_mix
    
    rho = mixture_density(Y, T, P, use_eos=True)  # PPR78 EOS density
    Cp_mix = mixture_cp(Y, T, P, use_eos=True)  # Ideal + residual Cp
    MW_mix = mixture_molar_mass(Y)  # kg/mol
    
    # Velocity from continuity using constant mass flux:
    # v_z = v_z0 * (rho_in / rho)
    if state.rho_in is None:
        state.rho_in = rho
    v_z = state.v_z0 * (state.rho_in / rho)
    v_z = max(v_z, 0.01)  # prevent zero velocity
    
    # Concentrations [mol/m³] using PPR78 EOS
    C_total = mixture_concentration(Y, T, P, use_eos=True)
    C = Y * C_total
    
    # Reaction rates
    omega = omega_rates(C, T, P)
    
    # Species ODEs: dY_i/dz = (1/v_z) * (1/C_total) * Σ_j ν_ij * Ω_j
    # (divide by C_total because Y_i = C_i / C_total)
    dYdz = (nu @ omega) / (v_z * C_total)
    
    # Energy ODE: dT/dz = (1/(ρ*Cp_mass*v_z)) * [Q_rxn + 4*q''/D]
    # where Cp_mass = Cp_mix / MW_mix [J/kg/K]
    # Q_rxn = -Σ_j ΔH_j * Ω_j [W/m³]
    dH = delta_H_reaction(T)
    Q_rxn = -np.dot(dH, omega)  # W/m³
    Cp_mass = Cp_mix / MW_mix  # J/kg/K
    dTdz = (1.0 / (rho * Cp_mass * v_z)) * (Q_rxn + 4.0 * state.q_flux / state.D)
    
    # Pressure ODE: Darcy-Weisbach
    # dP/dz = -4*f/D * ρ*v_z²/2
    Re = rho * v_z * state.D / mu_mix
    Re = max(Re, 100.0)  # avoid numerical issues
    f = 0.0791 * Re ** (-0.25)  # Blasius for turbulent
    dPdz = -4.0 * f / state.D * 0.5 * rho * v_z ** 2
    
    return np.concatenate([dYdz, [dTdz, dPdz]])

# -----------------------------------------------------------------------------
# Integration wrapper
# -----------------------------------------------------------------------------

def run_pfr(L: float = L_TUBE, T_in: float = T_SET, P_in: float = P_IN,
            v_z0: float = V_Z0, D: float = D_TUBE,
            feed: Dict[str, float] | None = None,
            q_flux: float | None = None,
            max_step: float = 0.01,
            mechanism_path: Path | str | None = None) -> Any:
    """Integrate PFR from z=0 to z=L.
    
    Args:
        L: reactor length [m]
        T_in: inlet temperature [K]
        P_in: inlet pressure [Pa]
        v_z0: inlet superficial velocity [m/s]
        D: tube diameter [m]
        feed: feed composition dict {formula: mole_frac}
        q_flux: wall heat flux [W/m²], None for auto-estimate
        max_step: max integration step [m]
    
    Returns:
        scipy OdeSolution object
    """
    if mechanism_path is not None:
        set_mechanism_path(mechanism_path)
    species_formula = get_species_list()
    spec_idx = get_species_index()
    report_skipped_reactions("run_pfr")
    n_spec = len(species_formula)
    
    # Build initial mole fraction vector
    Y0 = np.zeros(n_spec)
    feed = feed or FEED_FORMULA
    for sp, frac in feed.items():
        if sp in spec_idx:
            Y0[spec_idx[sp]] = frac
        else:
            raise ValueError(f"Unknown species in feed: {sp}")
    Y0 /= Y0.sum()
    
    # Build initial state vector
    y0 = np.concatenate([Y0, [T_in, P_in]])
    
    # Create state container using PPR78 EOS for inlet density
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
    
    # Integrate
    sol = solve_ivp(
        lambda z, y: ode_pfr(z, y, state),
        (0.0, L),
        y0,
        method="BDF",
        atol=1e-10,
        rtol=1e-8,
        max_step=max_step,
    )
    
    # Attach metadata
    sol.species = species_formula
    sol.n_spec = n_spec
    sol.state = state
    
    return sol


def print_results(sol) -> None:
    """Print summary of PFR solution."""
    spec_idx = get_species_index()
    n_spec = sol.n_spec
    
    print(f"Integration success: {sol.success}")
    print(f"Message: {sol.message}")
    print(f"Number of z-steps: {len(sol.t)}")
    print(f"z range: {sol.t[0]:.3f} to {sol.t[-1]:.3f} m")
    print()
    
    # Renormalize outlet mole fractions for display
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
    
    print("Mole fractions (inlet -> outlet, normalized):")
    for sp in ["C2H6", "C2H4", "C2H2", "C2H5", "C2H3", "C1H4", "C1H3", "H2", "H", "C3H6", "C3H7", "N2"]:
        if sp in spec_idx:
            i = spec_idx[sp]
            if Y_out[i] > 1e-6 or Y_in[i] > 1e-6:
                print(f"  {sp:6s}: {Y_in[i]:.6f} -> {Y_out[i]:.6f}")
    
    print(f"\n  Sum Y: {Y_out.sum():.6f}")
    
    # Conversion based on normalized fractions
    if "C2H6" in spec_idx:
        i_eth = spec_idx["C2H6"]
        conv = 1.0 - Y_out[i_eth] / Y_in[i_eth]
        print(f"\nEthane conversion: {conv*100:.1f}%")
        
        # Selectivity to ethylene (C2H4)
        if "C2H4" in spec_idx and conv > 0.01:
            i_ene = spec_idx["C2H4"]
            ethylene_formed = Y_out[i_ene] - Y_in[i_ene]
            ethane_reacted = Y_in[i_eth] - Y_out[i_eth]
            selectivity = ethylene_formed / ethane_reacted if ethane_reacted > 1e-10 else 0
            print(f"Ethylene selectivity: {selectivity*100:.1f}%")
        
        # Selectivity to acetylene (C2H2)  
        if "C2H2" in spec_idx and conv > 0.01:
            i_ace = spec_idx["C2H2"]
            acetylene_formed = Y_out[i_ace] - Y_in[i_ace]
            ethane_reacted = Y_in[i_eth] - Y_out[i_eth]
            selectivity_ace = acetylene_formed / ethane_reacted if ethane_reacted > 1e-10 else 0
            print(f"Acetylene selectivity: {selectivity_ace*100:.1f}%")
    
    # Heat flux used
    if hasattr(sol, 'state'):
        print(f"Wall heat flux: {sol.state.q_flux/1e3:.2f} kW/m²")


# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    print("Stage-1 Ideal PFR - Ethane Pyrolysis with N2 Dilution")
    print("=" * 60)
    print(f"Tube: D={D_TUBE*1000:.0f} mm, L={L_TUBE:.1f} m")
    print(f"Inlet: T={T_SET:.0f} K, P={P_IN/1e5:.1f} bar, v_z={V_Z0:.1f} m/s")
    print(f"N2 dilution: {DEFAULT_N2_DILUTION*100:.0f}% (mole basis)")
    print()
    
    # Run with default N2 dilution
    sol = run_pfr()
    print_results(sol)
