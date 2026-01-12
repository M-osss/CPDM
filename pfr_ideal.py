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
if __name__ == "__main__":
    from kinetics_parser import parse_csv_mechanism, Reaction
    import props
    from props import cp_molar, sensible_enthalpy, gas_viscosity, resolve_species
else:
    from .kinetics_parser import parse_csv_mechanism, Reaction
    from . import props
    from .props import cp_molar, sensible_enthalpy, gas_viscosity, resolve_species

ROOT = Path(__file__).resolve().parent

# -----------------------------------------------------------------------------
# Physical constants & reactor geometry
# -----------------------------------------------------------------------------
R_GAS = 8.314462618  # J/mol/K
D_TUBE = 0.1        # m tube diameter
L_TUBE = 10.0         # m tube length

# Operating conditions (Stage-1 defaults)
T_SET = 1000.0       # K target temperature
P_IN = 2.0 * 101_325.0  # Pa (2 bar)
V_Z0 = 1.0           # m/s feed superficial velocity

# Feed composition: pure ethane
FEED_FORMULA = {"C2H6": 1.0}

# -----------------------------------------------------------------------------
# Build kinetics data structures at module load
# -----------------------------------------------------------------------------

def _build_kinetics():
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
    """
    rxns = parse_csv_mechanism()
    
    # Collect all species from mechanism
    species_formula: List[str] = []
    for r in rxns:
        for sp in r["stoich"]:
            if sp not in species_formula:
                species_formula.append(sp)
    
    # Resolve to database keys
    species_keys = [resolve_species(sp) for sp in species_formula]
    
    n = len(species_formula)
    m = len(rxns)
    
    nu = np.zeros((n, m))
    A = np.zeros(m)
    n_exp = np.zeros(m)
    Ea = np.zeros(m)
    reactant_orders: List[Dict[int, float]] = []
    has_third_body = np.zeros(m, dtype=bool)
    
    for j, r in enumerate(rxns):
        A[j] = r["A_si"]
        n_exp[j] = r["n"]
        Ea[j] = r["Ea_si"]
        has_third_body[j] = r.get("has_third_body", False)
        
        for sp, coeff in r["stoich"].items():
            nu[species_formula.index(sp), j] = coeff
        
        # Build reactant order dict for this reaction
        rord = {}
        for sp, order in r.get("reactants", {}).items():
            rord[species_formula.index(sp)] = order
        reactant_orders.append(rord)
    
    return species_formula, species_keys, nu, A, n_exp, Ea, reactant_orders, has_third_body

# Module-level cached kinetics
_KINETICS = None

def _get_kinetics():
    global _KINETICS
    if _KINETICS is None:
        _KINETICS = _build_kinetics()
    return _KINETICS

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


def mixture_cp(Y: np.ndarray, T: float) -> float:
    """Molar-fraction weighted heat capacity [J/mol/K]."""
    _, species_keys, *_ = _get_kinetics()
    cps = np.array([cp_molar(T, key) for key in species_keys])
    return float(np.dot(Y, cps))


def mixture_density(Y: np.ndarray, T: float, P: float) -> float:
    """Mixture mass density [kg/m³] from ideal gas law."""
    MW = get_MW_array()
    MW_mix = float(np.dot(Y, MW))  # kg/mol
    return P * MW_mix / (R_GAS * T)


def mixture_molar_mass(Y: np.ndarray) -> float:
    """Average molar mass [kg/mol]."""
    return float(np.dot(Y, get_MW_array()))

# -----------------------------------------------------------------------------
# Reaction rate calculations
# -----------------------------------------------------------------------------

def reaction_rate_constants(T: float) -> np.ndarray:
    """Arrhenius rate constants k_j(T).
    
    k = A * (T/298.15)^n * exp(-Ea/(R*T))
    
    Units depend on reaction order (mechanism uses molecule-based units,
    converted at parse time).
    """
    _, _, _, A, n_exp, Ea, *_ = _get_kinetics()
    return A * (T / 298.15) ** n_exp * np.exp(-Ea / (R_GAS * T))


def omega_rates(C: np.ndarray, T: float, P: float) -> np.ndarray:
    """Compute reaction rates Ω_j [mol/m³/s].
    
    For elementary reactions: Ω_j = k_j * Π_i C_i^order_i
    For third-body reactions: multiply by total concentration [M] = P/(R*T).
    
    Args:
        C: species concentrations [mol/m³]
        T: temperature [K]
        P: pressure [Pa]
    
    Returns:
        omega: array of reaction rates (n_reactions,)
    """
    _, _, _, _, _, _, reactant_orders, has_third_body = _get_kinetics()
    k = reaction_rate_constants(T)
    m = len(k)
    omega = np.zeros(m)
    
    C_total = P / (R_GAS * T)  # total molar concentration [mol/m³]
    
    for j in range(m):
        rate = k[j]
        for i_sp, order in reactant_orders[j].items():
            rate *= max(C[i_sp], 1e-30) ** order
        if has_third_body[j]:
            rate *= C_total
        omega[j] = rate
    
    return omega


def delta_H_reaction(T: float) -> np.ndarray:
    """Compute ΔH_j(T) for each reaction [J/mol].
    
    ΔH_j = Σ_i ν_ij * h_i(T)  (products - reactants)
    """
    _, species_keys, nu, *_ = _get_kinetics()
    h = np.array([sensible_enthalpy(T, key) for key in species_keys])
    return nu.T @ h  # shape (m,)

# -----------------------------------------------------------------------------
# Heat flux estimation (Stage-1: constant q'' to maintain near-isothermal)
# -----------------------------------------------------------------------------

def estimate_heat_flux(T: float = T_SET, P: float = P_IN, 
                       Y0: np.ndarray | None = None,
                       target_conversion: float = 0.6) -> float:
    """Estimate constant heat flux q'' [W/m²] to maintain near-isothermal operation.
    
    Uses an iterative approach: run short integration segments to estimate
    average heat demand along the reactor, then compute required q''.
    
    For Stage-1, we use a simplified estimate based on expected conversion
    and reaction enthalpy.
    """
    species_formula, species_keys, nu, *_ = _get_kinetics()
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
    
    # Molar flow rate estimate [mol/s]
    C_total = P / (R_GAS * T)  # mol/m³
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
                 q_flux: float | None = None):
        self.T_in = T_in
        self.P_in = P_in
        self.v_z0 = v_z0
        self.D = D
        self.q_flux = q_flux if q_flux is not None else estimate_heat_flux(T_in, P_in)
        self.mu_in: float | None = None  # set after first call


def ode_pfr(z: float, y: np.ndarray, state: PFRState) -> np.ndarray:
    """Right-hand side of PFR ODE system.
    
    State vector: y = [Y_0, ..., Y_{n-1}, T, P]
    
    Returns: dy/dz
    """
    species_formula, _, nu, *_ = _get_kinetics()
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
    
    rho = mixture_density(Y, T, P)
    Cp_mix = mixture_cp(Y, T)  # J/mol/K
    MW_mix = mixture_molar_mass(Y)  # kg/mol
    
    # Velocity from continuity: v_z = v_z0 * (P0/P) * (T/T0) * (μ0/μ)
    # The viscosity ratio accounts for composition changes affecting flow
    v_z = state.v_z0 * (state.P_in / P) * (T / state.T_in)
    v_z = max(v_z, 0.01)  # prevent zero velocity
    
    # Concentrations [mol/m³]
    C_total = P / (R_GAS * T)
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
            max_step: float = 0.01) -> Any:
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
    species_formula = get_species_list()
    spec_idx = get_species_index()
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
    
    # Create state container
    state = PFRState(T_in=T_in, P_in=P_in, v_z0=v_z0, D=D, q_flux=q_flux)
    
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
    for sp in ["C2H6", "C2H4", "C2H2", "C2H5", "C2H3", "C1H4", "C1H3", "H2", "H", "C3H6", "C3H7"]:
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
    
    # Heat flux used
    if hasattr(sol, 'state'):
        print(f"Wall heat flux: {sol.state.q_flux/1e3:.2f} kW/m²")


# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    print("Stage-1 Ideal PFR - Ethane Pyrolysis")
    print("=" * 50)
    print(f"Tube: D={D_TUBE*1000:.0f} mm, L={L_TUBE:.1f} m")
    print(f"Inlet: T={T_SET:.0f} K, P={P_IN/1e5:.1f} bar, v_z={V_Z0:.1f} m/s")
    print()
    
    sol = run_pfr()
    print_results(sol)
