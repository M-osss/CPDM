"""Stage-2: wall conduction, WSGG radiation, operator-split coupling."""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Tuple

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_banded

try:
    from .pfr_ideal import (
        run_pfr, get_species_list, get_species_index, get_MW_array,
        mixture_viscosity, mixture_cp, mixture_density, mixture_molar_mass,
        mixture_concentration,
        omega_rates, delta_H_reaction, _get_kinetics, R_GAS,
        set_mechanism_path, report_skipped_reactions,
    )
    from . import props
except ImportError:
    from pfr_ideal import (
        run_pfr, get_species_list, get_species_index, get_MW_array,
        mixture_viscosity, mixture_cp, mixture_density, mixture_molar_mass,
        mixture_concentration,
        omega_rates, delta_H_reaction, _get_kinetics, R_GAS,
        set_mechanism_path, report_skipped_reactions,
    )
    import props

ROOT = Path(__file__).resolve().parent

DEFAULT_N2_DILUTION = 0.1
FEED_FORMULA = {"C2H6": 1.0 - DEFAULT_N2_DILUTION, "N2": DEFAULT_N2_DILUTION}

STEFAN_BOLTZMANN = 5.670374419e-8
PI = math.pi


@dataclass
class WallMaterial:
    name: str = "Inconel 800H"
    rho: float = 7940.0
    k: float = 25.0
    cp: float = 500.0
    emissivity: float = 0.85
    
    def thermal_diffusivity(self) -> float:
        return self.k / (self.rho * self.cp)


@dataclass
class WSGGModel:
    kappa: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.4, 6.5, 130.0]))
    
    b_matrix: np.ndarray = field(default_factory=lambda: np.array([
        [0.364, -1.01e-4, 0.0],
        [0.252,  1.21e-4, -2.1e-8],
        [0.272,  7.53e-5, -1.5e-8],
        [0.112, -8.72e-5,  1.8e-8],
    ]))
    
    def weights(self, T: float) -> np.ndarray:
        T_arr = np.array([1.0, T, T**2])
        return self.b_matrix @ T_arr
    
    def absorption_coefficient(self, T: float, P_partial_atm: float, 
                               path_length: float) -> float:
        if P_partial_atm < 1e-6:
            return 0.01
        
        a = self.weights(T)
        kappa_eff = 0.0
        for i in range(len(self.kappa)):
            kappa_eff += a[i] * self.kappa[i]
        
        return kappa_eff * P_partial_atm
    
    def emissivity_gas(self, T: float, P_partial_atm: float, 
                       path_length: float) -> float:
        if P_partial_atm < 1e-6:
            return 0.01
        
        a = self.weights(T)
        eps = 0.0
        for i in range(len(self.kappa)):
            tau = self.kappa[i] * P_partial_atm * path_length
            eps += a[i] * (1.0 - math.exp(-tau))
        
        return max(0.01, min(eps, 0.99))


def nusselt_dittus_boelter(Re: float, Pr: float, heating: bool = True) -> float:
    n = 0.4 if heating else 0.3
    return 0.023 * Re**0.8 * Pr**n


def nusselt_gnielinski(Re: float, Pr: float) -> float:
    if Re < 3000:
        Re = 3000
    
    f = (0.79 * math.log(Re) - 1.64)**(-2)
    num = (f / 8) * (Re - 1000) * Pr
    denom = 1.0 + 12.7 * (f / 8)**0.5 * (Pr**(2/3) - 1.0)
    return num / denom


def prandtl_number(mu: float, cp_mass: float, k: float) -> float:
    return mu * cp_mass / k


def gas_thermal_conductivity(T: float, Y: np.ndarray) -> float:
    return 0.02 + 7e-5 * T


@dataclass
class TubeWall:
    r_inner: float
    r_outer: float
    Nr: int = 5
    material: WallMaterial = field(default_factory=WallMaterial)
    
    T_wall: np.ndarray = field(init=False)
    r_nodes: np.ndarray = field(init=False)
    
    def __post_init__(self):
        self.r_nodes = np.linspace(self.r_inner, self.r_outer, self.Nr)
        self.T_wall = np.ones(self.Nr) * 1100.0
    
    @property
    def thickness(self) -> float:
        return self.r_outer - self.r_inner
    
    @property
    def dr(self) -> float:
        return (self.r_outer - self.r_inner) / (self.Nr - 1)
    
    def steady_state_conduction(self, q_inner: float, T_outer: float) -> np.ndarray:
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
        Nr = self.Nr
        dr = self.dr
        alpha = self.material.thermal_diffusivity()
        k = self.material.k
        rho_cp = self.material.rho * self.material.cp
        
        Fo = alpha * dt / dr**2
        theta = 0.5
        
        r = self.r_nodes
        
        diag = np.ones(Nr)
        sub = np.zeros(Nr)
        sup = np.zeros(Nr)
        rhs = self.T_wall.copy()
        
        for i in range(1, Nr - 1):
            r_m = (r[i-1] + r[i]) / 2
            r_p = (r[i] + r[i+1]) / 2
            
            c_m = theta * Fo * r_m / (r[i] * dr)
            c_p = theta * Fo * r_p / (r[i] * dr)
            
            sub[i] = -c_m
            sup[i] = -c_p
            diag[i] = 1.0 + c_m + c_p
            
            c_m_exp = (1 - theta) * Fo * r_m / (r[i] * dr)
            c_p_exp = (1 - theta) * Fo * r_p / (r[i] * dr)
            
            rhs[i] = (1.0 - c_m_exp - c_p_exp) * self.T_wall[i] \
                   + c_m_exp * self.T_wall[i-1] \
                   + c_p_exp * self.T_wall[i+1]
        
        diag[0] = 1.0
        sup[0] = -1.0
        rhs[0] = -q_inner * dr / k
        
        Bi = h_outer * dr / k
        diag[-1] = 1.0 + Bi
        sub[-1] = -1.0
        rhs[-1] = Bi * T_furnace
        
        ab = np.zeros((3, Nr))
        ab[0, 1:] = sup[:-1]
        ab[1, :] = diag
        ab[2, :-1] = sub[1:]
        
        self.T_wall = solve_banded((1, 1), ab, rhs)
        return self.T_wall
    
    @property
    def T_inner(self) -> float:
        return self.T_wall[0]
    
    @property
    def T_outer(self) -> float:
        return self.T_wall[-1]


@dataclass
class Stage2State:
    D_inner: float = 0.025
    wall_thickness: float = 0.008
    L: float = 12.0
    T_in: float = 1000.0
    P_in: float = 10.0 * 101325.0
    v_z0: float = 20
    T_furnace: float = 1450.0
    mechanism_path: Path | str | None = "reduced"
    h_outer: float = 30.0
    Nz: int = 100
    Nr_wall: int = 5
    
    wall: TubeWall = field(init=False)
    wsgg: WSGGModel = field(default_factory=WSGGModel)
    wall_material: WallMaterial = field(default_factory=WallMaterial)
    
    z_nodes: np.ndarray = field(init=False)
    T_gas: np.ndarray = field(init=False)
    T_wall_inner: np.ndarray = field(init=False)
    T_wall_outer: np.ndarray = field(init=False)
    q_inner: np.ndarray = field(init=False)
    
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


def compute_inner_heat_flux(T_gas: float, T_wall_inner: float,
                           Y: np.ndarray, P: float, D: float,
                           v_z: float, rho: float, mu: float,
                           wsgg: WSGGModel) -> Tuple[float, float, float]:
    MW = get_MW_array()
    MW_mix = float(np.dot(Y, MW))
    species_keys = _get_kinetics()[1]
    
    Cp_mix = 0.0
    for i, key in enumerate(species_keys):
        Cp_mix += Y[i] * props.cp_molar(T_gas, key)
    Cp_mass = Cp_mix / MW_mix
    
    k_gas = gas_thermal_conductivity(T_gas, Y)
    
    Pr = prandtl_number(mu, Cp_mass, k_gas)
    Pr = max(0.6, min(Pr, 100.0))
    
    Re = rho * v_z * D / mu
    Re = max(3000, Re)
    
    Nu = nusselt_gnielinski(Re, Pr)
    h_i = Nu * k_gas / D
    q_conv = h_i * (T_wall_inner - T_gas)
    
    L_beam = 0.9 * D
    P_atm = P / 101325.0
    P_participating = 0.0
    eps_g = wsgg.emissivity_gas(T_gas, P_participating * P_atm, L_beam)
    eps_w = 0.85
    
    q_rad = STEFAN_BOLTZMANN * eps_w * (T_wall_inner**4 - T_gas**4)
    q_total = q_conv + q_rad
    
    return q_total, q_conv, q_rad


def compute_outer_heat_flux(T_wall_outer: float, T_furnace: float,
                           h_outer: float, eps_wall: float = 0.85) -> Tuple[float, float, float]:
    q_conv = h_outer * (T_furnace - T_wall_outer)
    q_rad = eps_wall * STEFAN_BOLTZMANN * (T_furnace**4 - T_wall_outer**4)
    return q_conv + q_rad, q_conv, q_rad


def stage2_ode(z: float, y: np.ndarray, state: Stage2State,
               q_flux_func: Callable[[float], float]) -> np.ndarray:
    set_mechanism_path(state.mechanism_path)
    species_formula, _, nu, *_ = _get_kinetics()
    n_spec = len(species_formula)
    
    Y = y[:n_spec].copy()
    Y = np.clip(Y, 1e-15, 1.0)
    Y /= Y.sum()
    T = max(y[n_spec], 300.0)
    P = max(y[n_spec + 1], 1e3)
    
    mu_mix = mixture_viscosity(Y, T)
    rho = mixture_density(Y, T, P, use_eos=True)
    Cp_mix = mixture_cp(Y, T, P, use_eos=True)
    MW_mix = mixture_molar_mass(Y)
    
    v_z = state.v_z0 * (state.P_in / P) * (T / state.T_in)
    v_z = max(v_z, 0.01)
    
    C_total = mixture_concentration(Y, T, P, use_eos=True)
    C = Y * C_total
    
    omega = omega_rates(C, T, P)
    dYdz = (nu @ omega) / (v_z * C_total)
    
    dH = delta_H_reaction(T)
    Q_rxn = -np.dot(dH, omega)
    Cp_mass = Cp_mix / MW_mix
    
    q_flux = q_flux_func(z)
    dTdz = (1.0 / (rho * Cp_mass * v_z)) * (Q_rxn + 4.0 * q_flux / state.D_inner)
    
    Re = rho * v_z * state.D_inner / mu_mix
    Re = max(Re, 100.0)
    f = 0.0791 * Re ** (-0.25)
    dPdz = -4.0 * f / state.D_inner * 0.5 * rho * v_z ** 2
    
    return np.concatenate([dYdz, [dTdz, dPdz]])


def solve_stage2(state: Stage2State,
                feed: Dict[str, float] | None = None,
                max_outer_iter: int = 50,
                tol_T: float = 1.0,
                relax: float = 0.3,
                verbose: bool = True) -> Dict[str, Any]:
    set_mechanism_path(state.mechanism_path)
    species_formula = get_species_list()
    spec_idx = get_species_index()
    report_skipped_reactions("solve_stage2")
    n_spec = len(species_formula)
    
    Y0 = np.zeros(n_spec)
    feed = feed or FEED_FORMULA
    for sp, frac in feed.items():
        if sp in spec_idx:
            Y0[spec_idx[sp]] = frac
        else:
            raise ValueError(f"Unknown species in feed: {sp}")
    Y0 /= Y0.sum()
    
    Nz = state.Nz
    z = state.z_nodes
    dz = state.dz
    
    T_gas = np.ones(Nz) * state.T_in
    P_gas = np.ones(Nz) * state.P_in
    Y_gas = np.zeros((n_spec, Nz))
    for i in range(Nz):
        Y_gas[:, i] = Y0
    
    for iz in range(Nz):
        state.T_wall_outer[iz] = state.T_in + 50
        state.T_wall_inner[iz] = state.T_in + 30
    
    q_inner_profile = np.ones(Nz) * 20000.0
    state.q_inner = q_inner_profile.copy()
    
    converged = False
    for outer_iter in range(max_outer_iter):
        T_gas_old = T_gas.copy()
        q_old = state.q_inner.copy()
        
        q_inner_new = np.zeros(Nz)
        
        for iz in range(Nz):
            Y = Y_gas[:, iz]
            T = max(T_gas[iz], 300.0)
            P = P_gas[iz]
            T_wi = float(np.clip(state.T_wall_inner[iz], 250.0, 3000.0))
            
            mu = mixture_viscosity(Y, T)
            rho = mixture_density(Y, T, P, use_eos=True)
            v_z = state.v_z0 * (state.P_in / P) * (T / state.T_in)
            v_z = max(v_z, 0.1)
            
            q_total, _, _ = compute_inner_heat_flux(
                T, T_wi, Y, P, state.D_inner, v_z, rho, mu, state.wsgg
            )
            
            q_total = float(np.clip(q_total, -200000.0, 200000.0))
            q_inner_new[iz] = q_total
        
        state.q_inner = relax * q_inner_new + (1 - relax) * q_old
        q_inner_profile = state.q_inner.copy()
        
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
        
        T_gas_new = sol.y[n_spec, :]
        T_gas = relax * T_gas_new + (1 - relax) * T_gas_old
        T_gas = np.clip(T_gas, 300.0, 2500.0)
        
        P_gas = sol.y[n_spec + 1, :]
        
        for iz in range(Nz):
            Y_gas[:, iz] = sol.y[:n_spec, iz]
            Y_gas[:, iz] = np.clip(Y_gas[:, iz], 1e-15, 1.0)
            Y_gas[:, iz] /= Y_gas[:, iz].sum()
        
        for iz in range(Nz):
            Y = Y_gas[:, iz]
            T = T_gas[iz]
            P = P_gas[iz]
            
            T_wo_old = state.T_wall_outer[iz]
            T_wi_old = state.T_wall_inner[iz]
            
            T_wo = T_wo_old
            T_wi = T_wi_old
            
            for wall_iter in range(30):
                q_out, _, _ = compute_outer_heat_flux(T_wo, state.T_furnace,
                                                       state.h_outer,
                                                       state.wall_material.emissivity)
                
                T_wi_calc = T_wo - q_out * state.wall_thickness / state.wall_material.k
                
                mu = mixture_viscosity(Y, T)
                rho = mixture_density(Y, T, P, use_eos=True)
                v_z = state.v_z0 * (state.P_in / P) * (T / state.T_in)
                v_z = max(v_z, 0.1)
                
                q_in, _, _ = compute_inner_heat_flux(
                    T, T_wi_calc, Y, P, state.D_inner, v_z, rho, mu, state.wsgg
                )
                
                res = q_out - q_in
                
                if abs(res) < 50:
                    T_wi = T_wi_calc
                    break
                
                dT = 0.05 * res / (state.h_outer + 20)
                T_wo += dT
                T_wo = float(np.clip(T_wo, 250.0, 3000.0))
                T_wi = T_wi_calc
            
            state.T_wall_outer[iz] = relax * T_wo + (1 - relax) * T_wo_old
            state.T_wall_inner[iz] = relax * T_wi + (1 - relax) * T_wi_old
        
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


def print_stage2_results(result: Dict[str, Any]) -> None:
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
    print(f"  q_avg:       {np.mean(q)/1e3:.2f} kW/m2")
    print()
    
    print("Mole fractions (inlet -> outlet):")
    Y_in = Y[:, 0] / Y[:, 0].sum()
    Y_out = Y[:, -1] / Y[:, -1].sum()
    
    for sp in ["C2H6", "C2H4", "C2H2", "C1H4", "H2", "N2"]:
        if sp in spec_idx:
            i = spec_idx[sp]
            if Y_out[i] > 1e-6 or Y_in[i] > 1e-6:
                print(f"  {sp:6s}: {Y_in[i]:.6f} -> {Y_out[i]:.6f}")
    
    if "C2H6" in spec_idx:
        i_eth = spec_idx["C2H6"]
        conv = 1.0 - Y_out[i_eth] / Y_in[i_eth]
        print(f"\nEthane conversion: {conv*100:.1f}%")
    
    if "C2H4" in spec_idx and "C2H6" in spec_idx:
        i_eth = spec_idx["C2H6"]
        i_ene = spec_idx["C2H4"]
        ethane_reacted = Y_in[i_eth] - Y_out[i_eth]
        ethylene_formed = Y_out[i_ene] - Y_in[i_ene]
        if ethane_reacted > 1e-10:
            selectivity = ethylene_formed / ethane_reacted
            print(f"C2H4 selectivity: {selectivity*100:.1f}%")


if __name__ == "__main__":
    print("Stage-2 PFR with Heat Transfer - Ethane Pyrolysis")
    print("=" * 60)
    
    feed_composition = FEED_FORMULA
    
    state = Stage2State(
        Nz=50,
        Nr_wall=5,
    )
    
    print(f"Tube: D_inner={state.D_inner*1000:.0f} mm, "
          f"wall={state.wall_thickness*1000:.0f} mm, L={state.L:.1f} m")
    print(f"Inlet: T={state.T_in:.0f} K, P={state.P_in/1e5:.1f} bar")
    print(f"Furnace: T={state.T_furnace:.0f} K, h_outer={state.h_outer:.0f} W/(m2 K)")
    print(f"Feed: {feed_composition} (N2 dilution: {DEFAULT_N2_DILUTION*100:.0f}%)")
    print()
    
    result = solve_stage2(state, feed=feed_composition, verbose=True)
    print_stage2_results(result)
