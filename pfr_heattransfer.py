"""Stage-2: 2-D Heat Transfer Augmentation for Ethane Pyrolysis PFR.

Extends Stage-1 ideal PFR with:
  • Tube wall radial conduction: ρ_w C_{p,w} ∂T_w/∂t = (1/r)∂/∂r(k_w r ∂T_w/∂r)
  • Gas-phase radiative transfer via WSGG (Weighted Sum of Gray Gases)
  • Coupled boundary conditions:
    - Inner: q''_inner = h_i(T_fluid - T_w,inner) + q_rad(fluid→wall)
    - Outer: q''_outer = h_o(T_w,outer - T_furnace)

Operator-split coupling: advance species/P with current wall T, then update wall T.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Tuple

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_banded

# Handle imports for both module and standalone execution
try:
    from .pfr_ideal import (
        run_pfr, get_species_list, get_species_index, get_MW_array,
        mixture_viscosity, mixture_cp, mixture_density, mixture_molar_mass,
        omega_rates, delta_H_reaction, _get_kinetics, R_GAS,
        FEED_FORMULA, make_diluted_feed
    )
    from . import props
except ImportError:
    from pfr_ideal import (
        run_pfr, get_species_list, get_species_index, get_MW_array,
        mixture_viscosity, mixture_cp, mixture_density, mixture_molar_mass,
        omega_rates, delta_H_reaction, _get_kinetics, R_GAS,
        FEED_FORMULA, make_diluted_feed
    )
    import props

ROOT = Path(__file__).resolve().parent

# =============================================================================
# Physical Constants
# =============================================================================
STEFAN_BOLTZMANN = 5.670374419e-8  # W/(m²·K⁴)
PI = math.pi

# =============================================================================
# Wall Material Properties (Inconel 800H - typical for ethylene furnaces)
# =============================================================================
@dataclass
class WallMaterial:
    """Tube wall material properties."""
    name: str = "Inconel 800H"
    rho: float = 7940.0          # kg/m³ density
    k: float = 25.0              # W/(m·K) thermal conductivity at ~1000K
    cp: float = 500.0            # J/(kg·K) specific heat
    emissivity: float = 0.85     # surface emissivity (oxidized)
    
    def thermal_diffusivity(self) -> float:
        """α = k/(ρ·Cp) [m²/s]"""
        return self.k / (self.rho * self.cp)


# =============================================================================
# WSGG (Weighted Sum of Gray Gases) Model
# =============================================================================
@dataclass
class WSGGModel:
    """Weighted Sum of Gray Gases for radiative properties.
    
    For ethane pyrolysis, participating species are primarily:
    - H2O (if steam dilution)
    - CO2 (minor from coke oxidation, typically negligible in pyrolysis)
    - CH4 (weak absorber)
    
    Using Smith et al. (1982) coefficients for H2O-CO2 mixtures.
    For pure hydrocarbon pyrolysis without steam, gas radiation is weak.
    """
    # 3-gray-gas + 1-clear model coefficients (simplified)
    # kappa_i [1/(m·atm)], a_i(T) weights
    # From Smith et al. (1982) for H2O/CO2 = 2
    kappa: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.4, 6.5, 130.0]))
    
    # Weight polynomial coefficients: a_i = b_{i,0} + b_{i,1}*T + b_{i,2}*T²
    # Simplified for T in [1000, 2000] K
    b_matrix: np.ndarray = field(default_factory=lambda: np.array([
        [0.364, -1.01e-4, 0.0],      # clear gas
        [0.252,  1.21e-4, -2.1e-8],  # gray 1
        [0.272,  7.53e-5, -1.5e-8],  # gray 2
        [0.112, -8.72e-5,  1.8e-8],  # gray 3
    ]))
    
    def weights(self, T: float) -> np.ndarray:
        """Compute WSGG weights a_i(T)."""
        T_arr = np.array([1.0, T, T**2])
        return self.b_matrix @ T_arr
    
    def absorption_coefficient(self, T: float, P_partial_atm: float, 
                               path_length: float) -> float:
        """Compute effective absorption coefficient [1/m].
        
        For ethane pyrolysis without steam, P_partial (H2O+CO2) ≈ 0.
        Returns a minimal value for gas-phase radiation.
        
        Args:
            T: gas temperature [K]
            P_partial_atm: partial pressure of participating species [atm]
            path_length: mean beam length [m]
        """
        if P_partial_atm < 1e-6:
            # No participating species - optically thin
            return 0.01  # minimal value [1/m]
        
        a = self.weights(T)
        # Planck mean approximation
        kappa_eff = 0.0
        for i in range(len(self.kappa)):
            kappa_eff += a[i] * self.kappa[i]
        
        return kappa_eff * P_partial_atm
    
    def emissivity_gas(self, T: float, P_partial_atm: float, 
                       path_length: float) -> float:
        """Gas emissivity using WSGG model.
        
        ε_g = Σ_i a_i(T) * [1 - exp(-κ_i * P * L)]
        """
        if P_partial_atm < 1e-6:
            return 0.01  # optically thin limit
        
        a = self.weights(T)
        eps = 0.0
        for i in range(len(self.kappa)):
            tau = self.kappa[i] * P_partial_atm * path_length
            eps += a[i] * (1.0 - math.exp(-tau))
        
        return max(0.01, min(eps, 0.99))


# =============================================================================
# Heat Transfer Correlations
# =============================================================================

def nusselt_dittus_boelter(Re: float, Pr: float, heating: bool = True) -> float:
    """Dittus-Boelter correlation for turbulent pipe flow.
    
    Nu = 0.023 * Re^0.8 * Pr^n
    n = 0.4 for heating, 0.3 for cooling
    
    Valid: Re > 10000, 0.6 < Pr < 160, L/D > 10
    """
    n = 0.4 if heating else 0.3
    return 0.023 * Re**0.8 * Pr**n


def nusselt_gnielinski(Re: float, Pr: float) -> float:
    """Gnielinski correlation for pipe flow.
    
    Nu = (f/8)(Re - 1000)Pr / [1 + 12.7(f/8)^0.5 (Pr^(2/3) - 1)]
    
    Valid: 3000 < Re < 5e6, 0.5 < Pr < 2000
    """
    if Re < 3000:
        Re = 3000  # clamp to validity range
    
    f = (0.79 * math.log(Re) - 1.64)**(-2)  # Petukhov friction
    num = (f / 8) * (Re - 1000) * Pr
    denom = 1.0 + 12.7 * (f / 8)**0.5 * (Pr**(2/3) - 1.0)
    return num / denom


def prandtl_number(mu: float, cp_mass: float, k: float) -> float:
    """Prandtl number: Pr = μ·Cp/k"""
    return mu * cp_mass / k


def gas_thermal_conductivity(T: float, Y: np.ndarray) -> float:
    """Estimate mixture thermal conductivity [W/(m·K)].
    
    Uses Eucken correlation: k = (Cp + 5R/4) * μ / M
    Simplified for light hydrocarbons.
    """
    # Simple correlation for hydrocarbon gases
    # k ≈ 0.03 + 6e-5 * T [W/m/K] at 1 bar
    return 0.02 + 7e-5 * T


# =============================================================================
# Tube Wall Model (Radial Conduction)
# =============================================================================

@dataclass
class TubeWall:
    """Radial wall conduction model with finite-volume discretization.
    
    Uses Nr radial nodes through wall thickness.
    """
    r_inner: float           # inner radius [m]
    r_outer: float           # outer radius [m]
    Nr: int = 5              # number of radial nodes
    material: WallMaterial = field(default_factory=WallMaterial)
    
    # State
    T_wall: np.ndarray = field(init=False)
    r_nodes: np.ndarray = field(init=False)
    
    def __post_init__(self):
        """Initialize radial grid and temperature."""
        self.r_nodes = np.linspace(self.r_inner, self.r_outer, self.Nr)
        self.T_wall = np.ones(self.Nr) * 1100.0  # initial guess [K]
    
    @property
    def thickness(self) -> float:
        return self.r_outer - self.r_inner
    
    @property
    def dr(self) -> float:
        return (self.r_outer - self.r_inner) / (self.Nr - 1)
    
    def steady_state_conduction(self, q_inner: float, T_outer: float) -> np.ndarray:
        """Solve steady-state radial conduction.
        
        1/r d/dr(r k dT/dr) = 0 with:
          - q_inner = -k dT/dr|_{r=r_i} (inner heat flux)
          - T(r_outer) = T_outer (outer boundary)
        
        Analytical solution for constant k:
        T(r) = T_outer + (q_inner * r_inner / k) * ln(r_outer / r)
        """
        k = self.material.k
        r_i = self.r_inner
        r_o = self.r_outer
        
        T = np.zeros(self.Nr)
        for i, r in enumerate(self.r_nodes):
            T[i] = T_outer + (q_inner * r_i / k) * math.log(r_o / r)
        
        self.T_wall = T
        return T
    
    def solve_transient(self, q_inner: float, h_outer: float, T_furnace: float,
                        dt: float) -> np.ndarray:
        """Advance wall temperature by dt using implicit scheme.
        
        ρ·Cp·∂T/∂t = (1/r)∂/∂r(k·r·∂T/∂r)
        
        BCs:
          - Inner: -k ∂T/∂r = q_inner
          - Outer: -k ∂T/∂r = h_outer (T - T_furnace)
        
        Uses Crank-Nicolson (θ = 0.5).
        """
        Nr = self.Nr
        dr = self.dr
        alpha = self.material.thermal_diffusivity()
        k = self.material.k
        rho_cp = self.material.rho * self.material.cp
        
        Fo = alpha * dt / dr**2  # Fourier number
        theta = 0.5  # Crank-Nicolson
        
        # Build tridiagonal system: A·T^{n+1} = B·T^n + source
        # Using banded storage for scipy.linalg.solve_banded
        
        # Coefficients for interior nodes using cylindrical coords
        r = self.r_nodes
        
        # Main diagonal, sub-diagonal, super-diagonal
        diag = np.ones(Nr)
        sub = np.zeros(Nr)
        sup = np.zeros(Nr)
        rhs = self.T_wall.copy()
        
        for i in range(1, Nr - 1):
            r_m = (r[i-1] + r[i]) / 2  # r_{i-1/2}
            r_p = (r[i] + r[i+1]) / 2  # r_{i+1/2}
            
            c_m = theta * Fo * r_m / (r[i] * dr)
            c_p = theta * Fo * r_p / (r[i] * dr)
            
            sub[i] = -c_m
            sup[i] = -c_p
            diag[i] = 1.0 + c_m + c_p
            
            # Explicit part
            c_m_exp = (1 - theta) * Fo * r_m / (r[i] * dr)
            c_p_exp = (1 - theta) * Fo * r_p / (r[i] * dr)
            
            rhs[i] = (1.0 - c_m_exp - c_p_exp) * self.T_wall[i] \
                   + c_m_exp * self.T_wall[i-1] \
                   + c_p_exp * self.T_wall[i+1]
        
        # Inner BC: -k dT/dr = q_inner → ghost point method
        # T[0] - T[1] ≈ -q_inner * dr / k
        diag[0] = 1.0
        sup[0] = -1.0
        rhs[0] = -q_inner * dr / k
        
        # Outer BC: -k dT/dr = h_outer (T - T_furnace)
        # k(T[-2] - T[-1])/dr = h_outer (T[-1] - T_furnace)
        Bi = h_outer * dr / k  # Biot number
        diag[-1] = 1.0 + Bi
        sub[-1] = -1.0
        rhs[-1] = Bi * T_furnace
        
        # Solve banded system
        ab = np.zeros((3, Nr))
        ab[0, 1:] = sup[:-1]   # super-diagonal
        ab[1, :] = diag        # main diagonal
        ab[2, :-1] = sub[1:]   # sub-diagonal
        
        self.T_wall = solve_banded((1, 1), ab, rhs)
        return self.T_wall
    
    @property
    def T_inner(self) -> float:
        return self.T_wall[0]
    
    @property
    def T_outer(self) -> float:
        return self.T_wall[-1]


# =============================================================================
# Stage-2 PFR State
# =============================================================================

@dataclass
class Stage2State:
    """Operating state for Stage-2 PFR with heat transfer."""
    # Geometry
    D_inner: float = 0.05           # m inner diameter
    wall_thickness: float = 0.008  # m (8 mm typical for cracker tubes)
    L: float = 10.0                # m length
    
    # Operating conditions
    T_in: float = 1000.0           # K inlet gas temperature
    P_in: float = 20.0 * 101325.0   # Pa inlet pressure
    v_z0: float = 5              # m/s inlet velocity
    T_furnace: float = 1500.0      # K furnace temperature (radiant section)
    
    # External heat transfer
    h_outer: float = 30.0          # W/(m²·K) outer convective coefficient
    
    # Discretization
    Nz: int = 100                  # axial nodes
    Nr_wall: int = 5               # radial wall nodes
    
    # Models
    wall: TubeWall = field(init=False)
    wsgg: WSGGModel = field(default_factory=WSGGModel)
    wall_material: WallMaterial = field(default_factory=WallMaterial)
    
    # Solution storage
    z_nodes: np.ndarray = field(init=False)
    T_gas: np.ndarray = field(init=False)      # (Nz,) gas temperature
    T_wall_inner: np.ndarray = field(init=False)  # (Nz,) inner wall T
    T_wall_outer: np.ndarray = field(init=False)  # (Nz,) outer wall T
    q_inner: np.ndarray = field(init=False)    # (Nz,) inner wall heat flux
    
    def __post_init__(self):
        r_inner = self.D_inner / 2
        r_outer = r_inner + self.wall_thickness
        self.wall = TubeWall(r_inner, r_outer, self.Nr_wall, self.wall_material)
        
        self.z_nodes = np.linspace(0, self.L, self.Nz)
        self.T_gas = np.ones(self.Nz) * self.T_in
        self.T_wall_inner = np.ones(self.Nz) * (self.T_in + 50)
        self.T_wall_outer = np.ones(self.Nz) * (self.T_in + 100)
        self.q_inner = np.zeros(self.Nz)
    
    @property
    def dz(self) -> float:
        return self.L / (self.Nz - 1)
    
    @property
    def D_outer(self) -> float:
        return self.D_inner + 2 * self.wall_thickness


# =============================================================================
# Heat Flux Calculations
# =============================================================================

def compute_inner_heat_flux(T_gas: float, T_wall_inner: float,
                           Y: np.ndarray, P: float, D: float,
                           v_z: float, rho: float, mu: float,
                           wsgg: WSGGModel) -> Tuple[float, float, float]:
    """Compute heat flux from wall to gas at inner surface.
    
    q''_inner = h_i(T_wall - T_gas) + q_rad(wall→gas) - q_rad(gas→wall)
    
    Returns:
        q_total: total heat flux [W/m²]
        q_conv: convective component [W/m²]
        q_rad: net radiative component [W/m²]
    """
    # Gas properties
    MW = get_MW_array()
    MW_mix = float(np.dot(Y, MW))
    species_keys = _get_kinetics()[1]
    
    # Cp estimation (simplified - use mixture rule)
    Cp_mix = 0.0
    for i, key in enumerate(species_keys):
        Cp_mix += Y[i] * props.cp_molar(T_gas, key)
    Cp_mass = Cp_mix / MW_mix  # J/(kg·K)
    
    # Thermal conductivity
    k_gas = gas_thermal_conductivity(T_gas, Y)
    
    # Prandtl number
    Pr = prandtl_number(mu, Cp_mass, k_gas)
    Pr = max(0.6, min(Pr, 100.0))  # clamp to validity range
    
    # Reynolds number
    Re = rho * v_z * D / mu
    Re = max(3000, Re)  # ensure turbulent
    
    # Nusselt number (Gnielinski)
    Nu = nusselt_gnielinski(Re, Pr)
    
    # Convective heat transfer coefficient
    h_i = Nu * k_gas / D
    
    # Convective heat flux (wall → gas is positive)
    q_conv = h_i * (T_wall_inner - T_gas)
    
    # Radiative heat flux
    # Mean beam length for cylinder: L_m ≈ 0.9 * D
    L_beam = 0.9 * D
    
    # Partial pressure of participating species (H2O, CO2)
    # In ethane pyrolysis without steam, these are negligible
    # Approximate: P_H2O ≈ 0, P_CO2 ≈ 0
    P_atm = P / 101325.0
    P_participating = 0.0  # negligible for dry pyrolysis
    
    # Gas emissivity
    eps_g = wsgg.emissivity_gas(T_gas, P_participating * P_atm, L_beam)
    
    # Wall emissivity (oxidized Inconel)
    eps_w = 0.85
    
    # Net radiative exchange (gray enclosure approximation)
    # q_rad = σ * (eps_w * T_wall^4 - eps_g * T_gas^4) * correction_factor
    # For nearly transparent gas, q_rad ≈ σ * eps_w * (T_wall^4 - T_gas^4)
    
    q_rad = STEFAN_BOLTZMANN * eps_w * (T_wall_inner**4 - T_gas**4)
    
    # Total heat flux to gas
    q_total = q_conv + q_rad
    
    return q_total, q_conv, q_rad


def compute_outer_heat_flux(T_wall_outer: float, T_furnace: float,
                           h_outer: float, eps_wall: float = 0.85) -> Tuple[float, float, float]:
    """Compute heat flux from furnace to outer wall.
    
    q''_outer = h_o(T_furnace - T_wall) + ε_w * σ * (T_furnace^4 - T_wall^4)
    
    Returns:
        q_total: total flux [W/m²]
        q_conv: convective component [W/m²]
        q_rad: radiative component [W/m²]
    """
    # Convective
    q_conv = h_outer * (T_furnace - T_wall_outer)
    
    # Radiative from furnace (assume furnace as blackbody enclosure)
    q_rad = eps_wall * STEFAN_BOLTZMANN * (T_furnace**4 - T_wall_outer**4)
    
    return q_conv + q_rad, q_conv, q_rad


# =============================================================================
# Coupled ODE System for Species + Energy + Pressure
# =============================================================================

def stage2_ode(z: float, y: np.ndarray, state: Stage2State,
               q_flux_func: Callable[[float], float]) -> np.ndarray:
    """Right-hand side of Stage-2 PFR ODE system.
    
    Same as Stage-1 but q''(z) comes from wall model.
    
    State vector: y = [Y_0, ..., Y_{n-1}, T, P]
    """
    species_formula, _, nu, *_ = _get_kinetics()
    n_spec = len(species_formula)
    
    # Unpack state
    Y = y[:n_spec].copy()
    Y = np.clip(Y, 1e-15, 1.0)
    Y /= Y.sum()
    T = max(y[n_spec], 300.0)
    P = max(y[n_spec + 1], 1e3)
    
    # Mixture properties
    mu_mix = mixture_viscosity(Y, T)
    rho = mixture_density(Y, T, P)
    Cp_mix = mixture_cp(Y, T)  # J/mol/K
    MW_mix = mixture_molar_mass(Y)  # kg/mol
    
    # Velocity from continuity
    v_z = state.v_z0 * (state.P_in / P) * (T / state.T_in)
    v_z = max(v_z, 0.01)
    
    # Concentrations [mol/m³]
    C_total = P / (R_GAS * T)
    C = Y * C_total
    
    # Reaction rates
    omega = omega_rates(C, T, P)
    
    # Species ODEs
    dYdz = (nu @ omega) / (v_z * C_total)
    
    # Energy ODE with spatially-varying heat flux
    dH = delta_H_reaction(T)
    Q_rxn = -np.dot(dH, omega)  # W/m³
    Cp_mass = Cp_mix / MW_mix  # J/kg/K
    
    q_flux = q_flux_func(z)  # W/m²
    dTdz = (1.0 / (rho * Cp_mass * v_z)) * (Q_rxn + 4.0 * q_flux / state.D_inner)
    
    # Pressure ODE
    Re = rho * v_z * state.D_inner / mu_mix
    Re = max(Re, 100.0)
    f = 0.0791 * Re ** (-0.25)
    dPdz = -4.0 * f / state.D_inner * 0.5 * rho * v_z ** 2
    
    return np.concatenate([dYdz, [dTdz, dPdz]])


# =============================================================================
# Operator-Split Solver
# =============================================================================

def solve_stage2(state: Stage2State,
                feed: Dict[str, float] | None = None,
                max_outer_iter: int = 50,
                tol_T: float = 1.0,
                relax: float = 0.3,
                verbose: bool = True) -> Dict[str, Any]:
    """Solve Stage-2 PFR with coupled heat transfer.
    
    Operator-split iteration:
    1. Initialize wall temperatures from furnace BC
    2. Compute inner heat flux from convection + radiation
    3. Integrate species/T/P ODEs with current q''(z)
    4. Update wall temperatures from new gas T profile (with under-relaxation)
    5. Repeat until convergence
    
    Args:
        state: Stage2State with geometry and operating conditions
        feed: feed composition {formula: mole_frac}
        max_outer_iter: max coupling iterations
        tol_T: convergence tolerance on gas temperature [K]
        relax: under-relaxation factor (0 < relax <= 1)
        verbose: print iteration info
    
    Returns:
        dict with solution arrays and convergence info
    """
    species_formula = get_species_list()
    spec_idx = get_species_index()
    n_spec = len(species_formula)
    
    # Build initial mole fraction vector
    Y0 = np.zeros(n_spec)
    feed = feed or FEED_FORMULA  # Default: diluted feed from pfr_ideal
    for sp, frac in feed.items():
        if sp in spec_idx:
            Y0[spec_idx[sp]] = frac
        else:
            raise ValueError(f"Unknown species in feed: {sp}")
    Y0 /= Y0.sum()
    
    # Initialize
    Nz = state.Nz
    z = state.z_nodes
    dz = state.dz
    
    # Initial gas T profile (from Stage-1 or guess)
    T_gas = np.ones(Nz) * state.T_in
    P_gas = np.ones(Nz) * state.P_in
    Y_gas = np.zeros((n_spec, Nz))
    for i in range(Nz):
        Y_gas[:, i] = Y0
    
    # Initialize wall from steady-state with estimated flux
    # Use furnace temperature and estimate heat balance
    for iz in range(Nz):
        # Start with wall close to gas temperature
        # This is more stable than starting from furnace temperature
        T_wo = state.T_in + 50  # slightly above gas T
        T_wi = state.T_in + 30
        state.T_wall_outer[iz] = T_wo
        state.T_wall_inner[iz] = T_wi
    
    # Initial heat flux estimate (conservative)
    q_inner_profile = np.ones(Nz) * 20000.0  # 20 kW/m2 initial guess
    state.q_inner = q_inner_profile.copy()
    
    # Outer iteration loop
    converged = False
    for outer_iter in range(max_outer_iter):
        T_gas_old = T_gas.copy()
        q_old = state.q_inner.copy()
        
        # Step 1: Compute heat flux profile from current wall state
        q_inner_new = np.zeros(Nz)
        
        for iz in range(Nz):
            Y = Y_gas[:, iz]
            T = max(T_gas[iz], 300.0)  # clamp to physical values
            P = P_gas[iz]
            T_wi = state.T_wall_inner[iz]
            
            # Ensure wall is hotter than gas (heating)
            T_wi = max(T_wi, T + 1.0)
            
            # Gas properties
            mu = mixture_viscosity(Y, T)
            rho = mixture_density(Y, T, P)
            v_z = state.v_z0 * (state.P_in / P) * (T / state.T_in)
            v_z = max(v_z, 0.1)  # prevent zero velocity
            
            # Heat flux to gas
            q_total, _, _ = compute_inner_heat_flux(
                T, T_wi, Y, P, state.D_inner, v_z, rho, mu, state.wsgg
            )
            
            # Clamp heat flux to reasonable range
            q_total = max(1000.0, min(q_total, 200000.0))  # 1-200 kW/m2
            q_inner_new[iz] = q_total
        
        # Apply under-relaxation to heat flux
        state.q_inner = relax * q_inner_new + (1 - relax) * q_old
        q_inner_profile = state.q_inner.copy()
        
        # Interpolation function for ODE solver
        def q_flux_interp(z_val: float) -> float:
            if z_val <= 0:
                return q_inner_profile[0]
            if z_val >= state.L:
                return q_inner_profile[-1]
            idx = z_val / dz
            i0 = int(idx)
            i1 = min(i0 + 1, Nz - 1)
            frac = idx - i0
            return (1 - frac) * q_inner_profile[i0] + frac * q_inner_profile[i1]
        
        # Step 2: Integrate species/T/P ODEs
        y0 = np.concatenate([Y0, [state.T_in, state.P_in]])
        
        sol = solve_ivp(
            lambda zz, yy: stage2_ode(zz, yy, state, q_flux_interp),
            (0.0, state.L),
            y0,
            method="BDF",
            t_eval=z,
            atol=1e-10,
            rtol=1e-8,
            max_step=dz,
        )
        
        if not sol.success:
            if verbose:
                print(f"ODE integration failed: {sol.message}")
            break
        
        # Extract solution with under-relaxation on temperature
        T_gas_new = sol.y[n_spec, :]
        T_gas = relax * T_gas_new + (1 - relax) * T_gas_old
        T_gas = np.clip(T_gas, 300.0, 2500.0)  # physical limits
        
        P_gas = sol.y[n_spec + 1, :]
        
        for iz in range(Nz):
            Y_gas[:, iz] = sol.y[:n_spec, iz]
            Y_gas[:, iz] = np.clip(Y_gas[:, iz], 1e-15, 1.0)
            Y_gas[:, iz] /= Y_gas[:, iz].sum()
        
        # Step 3: Update wall temperatures from new gas profile
        for iz in range(Nz):
            Y = Y_gas[:, iz]
            T = T_gas[iz]
            P = P_gas[iz]
            
            # Solve for wall temperatures using energy balance
            T_wo_old = state.T_wall_outer[iz]
            T_wi_old = state.T_wall_inner[iz]
            
            # Target: find wall T such that q_out = q_in at steady state
            # q_out = h_o*(T_furn - T_wo) + eps*sigma*(T_furn^4 - T_wo^4)
            # q_in = h_i*(T_wi - T_gas) + radiation
            # q_wall = k*(T_wo - T_wi)/thickness
            
            # Iterate for wall temperature
            T_wo = T_wo_old
            T_wi = T_wi_old
            
            for wall_iter in range(30):
                # Outer heat flux from furnace
                q_out, _, _ = compute_outer_heat_flux(T_wo, state.T_furnace,
                                                       state.h_outer,
                                                       state.wall_material.emissivity)
                
                # Inner surface temperature from wall conduction
                # q_wall = k*(T_wo - T_wi)/t  (thin wall approximation)
                T_wi_calc = T_wo - q_out * state.wall_thickness / state.wall_material.k
                
                # Ensure T_wi > T_gas
                T_wi_calc = max(T_wi_calc, T + 1.0)
                
                # Inner heat flux to gas
                mu = mixture_viscosity(Y, T)
                rho = mixture_density(Y, T, P)
                v_z = state.v_z0 * (state.P_in / P) * (T / state.T_in)
                v_z = max(v_z, 0.1)
                
                q_in, _, _ = compute_inner_heat_flux(
                    T, T_wi_calc, Y, P, state.D_inner, v_z, rho, mu, state.wsgg
                )
                
                # Residual: q_out should equal q_in for steady-state wall
                res = q_out - q_in
                
                if abs(res) < 50:  # 50 W/m2 tolerance
                    T_wi = T_wi_calc
                    break
                
                # Adjust outer wall temperature
                # If q_out > q_in, wall is receiving more heat than it transfers
                # -> wall temperature should increase
                dT = 0.05 * res / (state.h_outer + 20)
                T_wo += dT
                
                # Clamp to physical range
                T_wo = max(T + 10.0, min(T_wo, state.T_furnace - 10.0))
                T_wi = T_wi_calc
            
            # Apply under-relaxation to wall temperatures
            state.T_wall_outer[iz] = relax * T_wo + (1 - relax) * T_wo_old
            state.T_wall_inner[iz] = relax * T_wi + (1 - relax) * T_wi_old
        
        # Check convergence
        dT_max = np.max(np.abs(T_gas - T_gas_old))
        
        if verbose:
            i_eth = spec_idx.get("C2H6", 0)
            conv = 1.0 - Y_gas[i_eth, -1] / Y_gas[i_eth, 0]
            print(f"Iter {outer_iter+1}: dT_max={dT_max:.2f} K, "
                  f"T_out={T_gas[-1]:.1f} K, conv={conv*100:.1f}%")
        
        if dT_max < tol_T:
            converged = True
            if verbose:
                print(f"Converged after {outer_iter+1} iterations")
            break
    
    # Build result dictionary
    result = {
        "z": z,
        "Y": Y_gas,
        "T_gas": T_gas,
        "P": P_gas,
        "T_wall_inner": state.T_wall_inner.copy(),
        "T_wall_outer": state.T_wall_outer.copy(),
        "q_inner": state.q_inner.copy(),
        "species": species_formula,
        "converged": converged,
        "n_iter": outer_iter + 1,
    }
    
    return result


# =============================================================================
# Results Display
# =============================================================================

def print_stage2_results(result: Dict[str, Any]) -> None:
    """Print summary of Stage-2 solution."""
    z = result["z"]
    Y = result["Y"]
    T_gas = result["T_gas"]
    P = result["P"]
    T_wi = result["T_wall_inner"]
    T_wo = result["T_wall_outer"]
    q = result["q_inner"]
    species = result["species"]
    
    spec_idx = {s: i for i, s in enumerate(species)}
    
    print(f"\nStage-2 PFR Results")
    print("=" * 60)
    print(f"Converged: {result['converged']} ({result['n_iter']} iterations)")
    print(f"z range: {z[0]:.3f} to {z[-1]:.3f} m")
    print()
    
    print("Axial profiles (inlet -> outlet):")
    print(f"  T_gas:       {T_gas[0]:.1f} -> {T_gas[-1]:.1f} K")
    print(f"  T_wall_in:   {T_wi[0]:.1f} -> {T_wi[-1]:.1f} K")
    print(f"  T_wall_out:  {T_wo[0]:.1f} -> {T_wo[-1]:.1f} K")
    print(f"  P:           {P[0]/1e5:.4f} -> {P[-1]/1e5:.4f} bar")
    print(f"  q_inner:     {q[0]/1e3:.2f} -> {q[-1]/1e3:.2f} kW/m2")
    print(f"  q_avg:       {np.mean(q)/1e3:.2f} kW/m²")
    print()
    
    print("Mole fractions (inlet -> outlet):")
    Y_in = Y[:, 0] / Y[:, 0].sum()
    Y_out = Y[:, -1] / Y[:, -1].sum()
    
    for sp in ["C2H6", "C2H4", "C2H2", "C1H4", "H2", "N2"]:
        if sp in spec_idx:
            i = spec_idx[sp]
            if Y_out[i] > 1e-6 or Y_in[i] > 1e-6:
                print(f"  {sp:6s}: {Y_in[i]:.6f} -> {Y_out[i]:.6f}")
    
    # Conversion
    if "C2H6" in spec_idx:
        i_eth = spec_idx["C2H6"]
        conv = 1.0 - Y_out[i_eth] / Y_in[i_eth]
        print(f"\nEthane conversion: {conv*100:.1f}%")
    
    # Selectivity
    if "C2H4" in spec_idx and "C2H6" in spec_idx:
        i_eth = spec_idx["C2H6"]
        i_ene = spec_idx["C2H4"]
        ethane_reacted = Y_in[i_eth] - Y_out[i_eth]
        ethylene_formed = Y_out[i_ene] - Y_in[i_ene]
        if ethane_reacted > 1e-10:
            selectivity = ethylene_formed / ethane_reacted
            print(f"C2H4 selectivity: {selectivity*100:.1f}%")


# =============================================================================
# Main Entry Point
# =============================================================================

if __name__ == "__main__":
    print("Stage-2 PFR with Heat Transfer - Ethane Pyrolysis")
    print("=" * 60)
    
    # Create state with typical industrial conditions
    # Uses Stage2State defaults (can be overridden)
    state = Stage2State(
        # D_inner=0.1,           # 100 mm tube (default)
        # wall_thickness=0.008,  # 8 mm wall (default)
        # L=10.0,                # 10 m length (default)
        # T_in=1000.0,           # 1000 K inlet (default)
        # P_in=2.0 * 101325,     # 2 bar (default)
        # v_z0=1.0,              # 1 m/s (default)
        # T_furnace=1300.0,      # 1300 K furnace (default)
        # h_outer=60.0,          # 60 W/(m²·K) (default)
        Nz=50,                 # 50 axial nodes
        Nr_wall=5,             # 5 radial wall nodes
    )
    
    print(f"Tube: D_inner={state.D_inner*1000:.0f} mm, "
          f"wall={state.wall_thickness*1000:.0f} mm, L={state.L:.1f} m")
    print(f"Inlet: T={state.T_in:.0f} K, P={state.P_in/1e5:.1f} bar")
    print(f"Furnace: T={state.T_furnace:.0f} K, h_outer={state.h_outer:.0f} W/(m²·K)")
    print()
    
    result = solve_stage2(state, feed={"C2H6": 1.0}, verbose=True)
    print_stage2_results(result)

