"""Stage-3: axial dispersion, Danckwerts BCs, Taylor-Aris."""
from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy.integrate import solve_bvp
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

try:
    from .pfr_ideal import (
        get_species_list, get_species_index, get_MW_array,
        mixture_viscosity, mixture_cp, mixture_density, mixture_molar_mass,
        mixture_concentration,
        omega_rates, delta_H_reaction, _get_kinetics, R_GAS,
        set_mechanism_path, report_skipped_reactions,
        D_TUBE, L_TUBE, T_SET, P_IN, V_Z0,
        FEED_FORMULA, make_diluted_feed
    )
    from .pfr_heattransfer import (
        Stage2State, solve_stage2, gas_thermal_conductivity
    )
    from . import props
except ImportError:
    from pfr_ideal import (
        get_species_list, get_species_index, get_MW_array,
        mixture_viscosity, mixture_cp, mixture_density, mixture_molar_mass,
        mixture_concentration,
        omega_rates, delta_H_reaction, _get_kinetics, R_GAS,
        set_mechanism_path, report_skipped_reactions,
        D_TUBE, L_TUBE, T_SET, P_IN, V_Z0,
        FEED_FORMULA, make_diluted_feed
    )
    from pfr_heattransfer import (
        Stage2State, solve_stage2, gas_thermal_conductivity
    )
    import props

ROOT = Path(__file__).resolve().parent

PI = math.pi


def taylor_aris_dispersion(D_mol: float, v_z: float, D_tube: float) -> float:
    return D_mol + (v_z**2 * D_tube**2) / (192 * D_mol)


def turbulent_dispersion(v_z: float, D_tube: float, Re: float) -> float:
    if Re < 2300:
        D_mol = 1e-5
        return taylor_aris_dispersion(D_mol, v_z, D_tube)
    
    inv_Pe = 0.3 / Re + 0.5 / (1.0 + 3.8 / Re)
    return inv_Pe * v_z * D_tube


def estimate_molecular_diffusivity(T: float, P: float, MW1: float, MW2: float) -> float:
    V_A = 40.0
    V_B = 6.12
    
    P_atm = P / 101325.0
    M_A = MW1 * 1000
    M_B = MW2 * 1000

    D_AB = (1.013e-2 * T**1.75 * math.sqrt(1/M_A + 1/M_B)) / \
           (P_atm * (V_A**(1/3) + V_B**(1/3))**2)
    
    return D_AB * 1e-4


@dataclass
class Stage3State:
    D: float = 0.025
    L: float = 10.0
    T_in: float = 1100.0
    P_in: float = 5.0 * 101325
    v_z0: float = 5
    mechanism_path: Path | str | None = "reduced"
    Nz: int = 100
    use_taylor_aris: bool = False
    D_ax_constant: float | None = None
    alpha_ax_constant: float | None = None
    q_flux_profile: np.ndarray | None = None
    q_flux_constant: float | None = None
    Pe_min: float = 10.0
    Pe_max: float = 1000.0
    
    z_nodes: np.ndarray = field(init=False)
    dz: float = field(init=False)
    _q_flux_auto: float = field(init=False)
    
    def __post_init__(self):
        self.z_nodes = np.linspace(0, self.L, self.Nz)
        self.dz = self.L / (self.Nz - 1)
        
        if self.q_flux_constant is None:
            try:
                from .pfr_ideal import estimate_heat_flux
            except ImportError:
                from pfr_ideal import estimate_heat_flux
            self._q_flux_auto = estimate_heat_flux(self.T_in, self.P_in)
        else:
            self._q_flux_auto = self.q_flux_constant
    
    def get_q_flux(self, z: float) -> float:
        if self.q_flux_profile is not None:
            idx = z / self.dz
            i0 = int(idx)
            i1 = min(i0 + 1, self.Nz - 1)
            frac = idx - i0
            return (1 - frac) * self.q_flux_profile[i0] + frac * self.q_flux_profile[i1]
        return self._q_flux_auto


def compute_D_ax_profile(state: Stage3State, T: np.ndarray, P: np.ndarray,
                         Y: np.ndarray) -> np.ndarray:
    Nz = state.Nz
    D_ax = np.zeros(Nz)
    
    if state.D_ax_constant is not None:
        return np.ones(Nz) * state.D_ax_constant
    
    for iz in range(Nz):
        Y_loc = Y[:, iz]
        T_loc = T[iz]
        P_loc = P[iz]
        
        v_z = state.v_z0 * (state.P_in / P_loc) * (T_loc / state.T_in)
        mu = mixture_viscosity(Y_loc, T_loc)
        rho = mixture_density(Y_loc, T_loc, P_loc, use_eos=True)
        Re = rho * v_z * state.D / mu
        
        if state.use_taylor_aris or Re < 2300:
            MW = get_MW_array()
            MW_mix = float(np.dot(Y_loc, MW))
            D_mol = estimate_molecular_diffusivity(T_loc, P_loc, MW_mix, 0.002)
            D_ax[iz] = taylor_aris_dispersion(D_mol, v_z, state.D)
        else:
            D_ax[iz] = turbulent_dispersion(v_z, state.D, Re)
    
    return D_ax


def compute_alpha_ax_profile(state: Stage3State, T: np.ndarray, P: np.ndarray,
                            Y: np.ndarray, D_ax: np.ndarray) -> np.ndarray:
    if state.alpha_ax_constant is not None:
        return np.ones(state.Nz) * state.alpha_ax_constant
    
    Pr_t = 0.85
    return D_ax / Pr_t


def solve_dispersion_fd(state: Stage3State,
                       feed: Dict[str, float] | None = None,
                       max_iter: int = 800,
                       tol: float = 1e-5,
                       relax: float = 0.3,
                       verbose: bool = True) -> Dict[str, Any]:
    try:
        from .pfr_ideal import run_pfr
    except ImportError:
        from pfr_ideal import run_pfr
    
    set_mechanism_path(state.mechanism_path)
    species_formula = get_species_list()
    spec_idx = get_species_index()
    report_skipped_reactions("solve_dispersion_fd")
    n_spec = len(species_formula)
    
    feed = feed or FEED_FORMULA
    
    Nz = state.Nz
    z = state.z_nodes
    dz = state.dz
    
    q_flux_used = state._q_flux_auto
    
    if verbose:
        print("Initializing from Stage-1 plug-flow solution...")
        print(f"Using heat flux: {q_flux_used/1e3:.2f} kW/m2")
    
    sol = run_pfr(
        L=state.L,
        T_in=state.T_in,
        P_in=state.P_in,
        v_z0=state.v_z0,
        D=state.D,
        feed=feed,
        q_flux=q_flux_used,
        max_step=dz / 2,
        mechanism_path=state.mechanism_path,
    )
    
    if not sol.success:
        raise RuntimeError(f"Stage-1 solver failed: {sol.message}")
    
    from scipy.interpolate import interp1d
    
    T = np.zeros(Nz)
    P = np.zeros(Nz)
    Y = np.zeros((n_spec, Nz))
    
    for iz, z_val in enumerate(z):
        idx = np.searchsorted(sol.t, z_val)
        if idx == 0:
            T[iz] = sol.y[n_spec, 0]
            P[iz] = sol.y[n_spec + 1, 0]
            Y[:, iz] = sol.y[:n_spec, 0]
        elif idx >= len(sol.t):
            T[iz] = sol.y[n_spec, -1]
            P[iz] = sol.y[n_spec + 1, -1]
            Y[:, iz] = sol.y[:n_spec, -1]
        else:
            z0, z1 = sol.t[idx-1], sol.t[idx]
            frac = (z_val - z0) / (z1 - z0)
            T[iz] = (1 - frac) * sol.y[n_spec, idx-1] + frac * sol.y[n_spec, idx]
            P[iz] = (1 - frac) * sol.y[n_spec+1, idx-1] + frac * sol.y[n_spec+1, idx]
            Y[:, iz] = (1 - frac) * sol.y[:n_spec, idx-1] + frac * sol.y[:n_spec, idx]
    
    for iz in range(Nz):
        Y[:, iz] = np.clip(Y[:, iz], 1e-15, 1.0)
        Y[:, iz] /= Y[:, iz].sum()
    
    if verbose:
        i_eth = spec_idx.get("C2H6", 0)
        conv_pf = 1.0 - Y[i_eth, -1] / Y[i_eth, 0]
        print(f"Plug-flow solution: T_out={T[-1]:.1f} K, conv={conv_pf*100:.1f}%")
    
    D_ax = compute_D_ax_profile(state, T, P, Y)
    alpha_ax = compute_alpha_ax_profile(state, T, P, Y, D_ax)
    
    v_z_avg = state.v_z0 * np.mean(T / state.T_in)
    Pe_est = v_z_avg * state.L / np.mean(D_ax)
    
    if verbose:
        print(f"Estimated Peclet number: {Pe_est:.1f}")
    
    if Pe_est < 10:
        print("WARNING: Low Peclet number - dispersion is significant.")
        print("         Results may be inaccurate; consider full BVP solver.")
    
    converged = False
    T_min = max(300.0, state.T_in - 200)
    T_max = min(2500.0, state.T_in + 500)
    
    for it in range(max_iter):
        Y_old = Y.copy()
        T_old = T.copy()
        
        D_ax = compute_D_ax_profile(state, T, P, Y)
        alpha_ax = compute_alpha_ax_profile(state, T, P, Y, D_ax)
        
        for i_sp in range(n_spec):
            d2Ydz2 = np.zeros(Nz)
            for iz in range(1, Nz - 1):
                d2Ydz2[iz] = (Y[i_sp, iz+1] - 2*Y[i_sp, iz] + Y[i_sp, iz-1]) / dz**2
            
            v_z_arr = state.v_z0 * (state.P_in / P) * (T / state.T_in)
            disp_corr = (D_ax / v_z_arr) * d2Ydz2 * dz
            
            for iz in range(1, Nz - 1):
                Y[i_sp, iz] += relax * disp_corr[iz]
        
        for iz in range(Nz):
            Y[:, iz] = np.clip(Y[:, iz], 1e-15, 1.0)
            Y[:, iz] /= Y[:, iz].sum()
        
        d2Tdz2 = np.zeros(Nz)
        for iz in range(1, Nz - 1):
            d2Tdz2[iz] = (T[iz+1] - 2*T[iz] + T[iz-1]) / dz**2
        
        v_z_arr = state.v_z0 * (state.P_in / P) * (T / state.T_in)
        disp_corr_T = (alpha_ax / v_z_arr) * d2Tdz2 * dz
        
        for iz in range(1, Nz - 1):
            T[iz] += relax * disp_corr_T[iz]
        
        T = np.clip(T, T_min, T_max)
        
        dY_max = np.max(np.abs(Y - Y_old))
        dT_max = np.max(np.abs(T - T_old))
        
        if verbose and (it % 20 == 0 or it < 5):
            i_eth = spec_idx.get("C2H6", 0)
            conv = 1.0 - Y[i_eth, -1] / max(Y[i_eth, 0], 1e-10)
            print(f"Dispersion iter {it+1}: dY={dY_max:.2e}, dT={dT_max:.2f} K, "
                  f"conv={conv*100:.1f}%")
        
        if dY_max < tol and dT_max < 0.5:
            converged = True
            if verbose:
                print(f"Dispersion corrections converged after {it+1} iterations")
            break
    
    if not converged and verbose:
        print(f"Dispersion corrections did not fully converge after {max_iter} iterations")
    
    v_z_avg = state.v_z0 * np.mean(T / state.T_in)
    Pe_mass = v_z_avg * state.L / np.mean(D_ax)
    Pe_heat = v_z_avg * state.L / np.mean(alpha_ax)
    
    return {
        "z": z,
        "Y": Y,
        "T": T,
        "P": P,
        "D_ax": D_ax,
        "alpha_ax": alpha_ax,
        "Pe_mass": Pe_mass,
        "Pe_heat": Pe_heat,
        "species": species_formula,
        "converged": converged,
        "n_iter": it + 1,
    }


def _solve_species_fd(state: Stage3State, Y: np.ndarray, T: np.ndarray,
                     P: np.ndarray, i_sp: int, Y_in: float,
                     D_ax: np.ndarray) -> np.ndarray:
    Nz = state.Nz
    dz = state.dz
    
    a = np.zeros(Nz)
    b = np.zeros(Nz)
    c = np.zeros(Nz)
    d = np.zeros(Nz)
    
    set_mechanism_path(state.mechanism_path)
    species_formula, _, nu, *_ = _get_kinetics()
    
    for iz in range(1, Nz - 1):
        Y_loc = Y[:, iz]
        T_loc = T[iz]
        P_loc = P[iz]
        
        C_total = mixture_concentration(Y_loc, T_loc, P_loc, use_eos=True)
        C = Y_loc * C_total
        
        v_z = state.v_z0 * (state.P_in / P_loc) * (T_loc / state.T_in)
        D_loc = max(D_ax[iz], 1e-6)
        
        omega = omega_rates(C, T_loc, P_loc)
        R_i = np.dot(nu[i_sp, :], omega)
        
        Pe_local = v_z * dz / D_loc
        
        a[iz] = 1.0 / dz**2 + v_z / (D_loc * dz)
        b[iz] = -2.0 / dz**2 - v_z / (D_loc * dz)
        c[iz] = 1.0 / dz**2
        d[iz] = -R_i / (C_total * D_loc)
    
    iz = 0
    v_z = state.v_z0
    D_loc = max(D_ax[iz], 1e-6)
    
    b[0] = 1.0 + v_z * dz / D_loc
    c[0] = -1.0
    d[0] = v_z * dz / D_loc * Y_in
    
    a[-1] = -1.0
    b[-1] = 1.0
    d[-1] = 0.0
    
    Y_new = _solve_tridiag(a, b, c, d)
    return Y_new


def _solve_energy_fd(state: Stage3State, Y: np.ndarray, T: np.ndarray,
                    P: np.ndarray, alpha_ax: np.ndarray) -> np.ndarray:
    Nz = state.Nz
    dz = state.dz
    
    a = np.zeros(Nz)
    b = np.zeros(Nz)
    c = np.zeros(Nz)
    d = np.zeros(Nz)
    
    for iz in range(1, Nz - 1):
        Y_loc = Y[:, iz]
        T_loc = T[iz]
        P_loc = P[iz]
        
        C_total = mixture_concentration(Y_loc, T_loc, P_loc, use_eos=True)
        C = Y_loc * C_total
        
        v_z = state.v_z0 * (state.P_in / P_loc) * (T_loc / state.T_in)
        alpha_loc = max(alpha_ax[iz], 1e-6)
        
        rho = mixture_density(Y_loc, T_loc, P_loc, use_eos=True)
        Cp_mix = mixture_cp(Y_loc, T_loc, P_loc, use_eos=True)
        MW_mix = mixture_molar_mass(Y_loc)
        Cp_mass = Cp_mix / MW_mix
        
        omega = omega_rates(C, T_loc, P_loc)
        dH = delta_H_reaction(T_loc)
        Q_rxn = -np.dot(dH, omega)
        
        q_flux = state.get_q_flux(state.z_nodes[iz])
        
        S = (Q_rxn + 4.0 * q_flux / state.D) / (rho * Cp_mass * alpha_loc)
        
        a[iz] = 1.0 / dz**2 + v_z / (alpha_loc * dz)
        b[iz] = -2.0 / dz**2 - v_z / (alpha_loc * dz)
        c[iz] = 1.0 / dz**2
        d[iz] = -S
    
    b[0] = 1.0
    c[0] = 0.0
    d[0] = state.T_in
    
    a[-1] = -1.0
    b[-1] = 1.0
    d[-1] = 0.0
    
    T_new = _solve_tridiag(a, b, c, d)
    return T_new


def _solve_tridiag(a: np.ndarray, b: np.ndarray, c: np.ndarray,
                  d: np.ndarray) -> np.ndarray:
    n = len(d)
    
    c_prime = np.zeros(n)
    d_prime = np.zeros(n)
    
    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]
    
    for i in range(1, n):
        denom = b[i] - a[i] * c_prime[i-1]
        if abs(denom) < 1e-15:
            denom = 1e-15
        c_prime[i] = c[i] / denom
        d_prime[i] = (d[i] - a[i] * d_prime[i-1]) / denom
    
    x = np.zeros(n)
    x[-1] = d_prime[-1]
    for i in range(n - 2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i+1]
    
    return x


def generate_openfoam_case(state: Stage3State,
                          output_dir: str | Path,
                          stage2_result: Dict[str, Any] | None = None) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    (output_dir / "0").mkdir(exist_ok=True)
    (output_dir / "constant").mkdir(exist_ok=True)
    (output_dir / "system").mkdir(exist_ok=True)
    
    R = state.D / 2
    L = state.L
    n_radial = 20
    n_axial = int(state.Nz)
    wedge_angle = 5.0
    
    _write_blockmesh(output_dir / "system" / "blockMeshDict",
                    R, L, n_radial, n_axial, wedge_angle)
    _write_controldict(output_dir / "system" / "controlDict")
    _write_fvschemes(output_dir / "system" / "fvSchemes")
    _write_fvsolution(output_dir / "system" / "fvSolution")
    _write_thermophysical(output_dir / "constant" / "thermophysicalProperties")
    _write_transport(output_dir / "constant" / "transportProperties")
    _write_turbulence(output_dir / "constant" / "turbulenceProperties")
    
    species_list = get_species_list()
    _write_reaction_sources(output_dir / "constant" / "reactionSources.json",
                           species_list, state, stage2_result)
    _write_boundary_conditions(output_dir / "0", state, stage2_result)
    
    print(f"OpenFOAM case generated in: {output_dir}")


def _write_blockmesh(filepath: Path, R: float, L: float,
                    n_radial: int, n_axial: int, wedge_angle: float) -> None:
    theta = math.radians(wedge_angle / 2)
    
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}

convertToMeters 1.0;

R {R};
L {L};
theta {theta};

vertices
(
    (0 0 0)
    ({R} 0 0)
    (#calc "$R*cos($theta)" #calc "$R*sin($theta)" 0)
    (#calc "$R*cos($theta)" #calc "-$R*sin($theta)" 0)
    
    (0 0 {L})
    ({R} 0 {L})
    (#calc "$R*cos($theta)" #calc "$R*sin($theta)" {L})
    (#calc "$R*cos($theta)" #calc "-$R*sin($theta)" {L})
);

blocks
(
    hex (0 1 2 0 4 5 6 4) ({n_radial} 1 {n_axial}) simpleGrading (1 1 1)
);

edges
(
    arc 1 2 (#calc "$R*cos($theta/2)" #calc "$R*sin($theta/2)" 0)
    arc 5 6 (#calc "$R*cos($theta/2)" #calc "$R*sin($theta/2)" {L})
);

boundary
(
    inlet
    {{
        type patch;
        faces
        (
            (0 1 2 0)
        );
    }}
    outlet
    {{
        type patch;
        faces
        (
            (4 5 6 4)
        );
    }}
    wall
    {{
        type wall;
        faces
        (
            (1 5 6 2)
        );
    }}
    wedgeFront
    {{
        type wedge;
        faces
        (
            (0 1 5 4)
        );
    }}
    wedgeBack
    {{
        type wedge;
        faces
        (
            (0 2 6 4)
        );
    }}
    axis
    {{
        type empty;
        faces
        (
            (0 4 4 0)
        );
    }}
);

mergePatchPairs
(
);
"""
    filepath.write_text(content)


def _write_controldict(filepath: Path) -> None:
    content = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}

application     reactingFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         1000;

deltaT          0.001;
writeControl    adjustableRunTime;
writeInterval   100;
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
writeCompression off;

timeFormat      general;
timePrecision   6;

runTimeModifiable true;

adjustTimeStep  yes;
maxCo           0.5;

functions
{
    fieldAverage
    {
        type            fieldAverage;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
        fields
        (
            U
            {
                mean        on;
                prime2Mean  on;
                base        time;
            }
            T
            {
                mean        on;
                prime2Mean  on;
                base        time;
            }
        );
    }
}
"""
    filepath.write_text(content)


def _write_fvschemes(filepath: Path) -> None:
    content = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}

ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,Yi_h)   bounded Gauss linearUpwind grad(Yi_h);
    div(phi,h)      bounded Gauss linearUpwind grad(h);
    div(phi,K)      bounded Gauss linearUpwind grad(K);
    div(phi,k)      bounded Gauss linearUpwind grad(k);
    div(phi,omega)  bounded Gauss linearUpwind grad(omega);
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}
"""
    filepath.write_text(content)


def _write_fvsolution(filepath: Path) -> None:
    content = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}

solvers
{
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-6;
        relTol          0.1;
    }

    pFinal
    {
        $p;
        relTol          0;
    }

    "(U|h|k|omega|Yi)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-6;
        relTol          0.1;
    }

    "(U|h|k|omega|Yi)Final"
    {
        $U;
        relTol          0;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
    residualControl
    {
        p               1e-4;
        U               1e-4;
        h               1e-4;
        "(k|omega)"     1e-4;
        "Yi.*"          1e-4;
    }
}

relaxationFactors
{
    fields
    {
        p               0.3;
    }
    equations
    {
        U               0.7;
        h               0.7;
        k               0.7;
        omega           0.7;
        "Yi.*"          0.7;
    }
}
"""
    filepath.write_text(content)


def _write_thermophysical(filepath: Path) -> None:
    content = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      thermophysicalProperties;
}

thermoType
{
    type            heRhoThermo;
    mixture         reactingMixture;
    transport       sutherland;
    thermo          janaf;
    energy          sensibleEnthalpy;
    equationOfState perfectGas;
    specie          specie;
}

chemistryReader foamChemistryReader;
foamChemistryFile "$FOAM_CASE/constant/reactions";
foamChemistryThermoFile "$FOAM_CASE/constant/thermo.compressibleGas";

inertSpecie N2;
"""
    filepath.write_text(content)


def _write_transport(filepath: Path) -> None:
    content = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}
"""
    filepath.write_text(content)


def _write_turbulence(filepath: Path) -> None:
    content = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}

simulationType  RAS;

RAS
{
    RASModel        kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}
"""
    filepath.write_text(content)


def _write_reaction_sources(filepath: Path, species_list: List[str],
                           state: Stage3State,
                           stage2_result: Dict[str, Any] | None) -> None:
    Nz = state.Nz
    z = state.z_nodes.tolist()
    
    sources = {
        "z": z,
        "species": species_list,
        "description": "Reaction source terms [mol/m3/s] from Python PFR model",
        "rates": {}
    }
    
    if stage2_result is not None:
        Y = stage2_result["Y"]
        T = stage2_result["T_gas"]
        P = stage2_result["P"]
        
        set_mechanism_path(state.mechanism_path)
        species_formula, _, nu, *_ = _get_kinetics()
        
        for iz in range(Nz):
            Y_loc = Y[:, iz]
            T_loc = T[iz]
            P_loc = P[iz]
            
            C_total = mixture_concentration(Y_loc, T_loc, P_loc, use_eos=True)
            C = Y_loc * C_total
            
            omega = omega_rates(C, T_loc, P_loc)
            R = nu @ omega
            
            sources["rates"][f"z_{iz}"] = {
                "z": z[iz],
                "T": T_loc,
                "P": P_loc,
                "rates": {sp: R[i] for i, sp in enumerate(species_list)}
            }
    
    with open(filepath, "w") as f:
        json.dump(sources, f, indent=2)


def _write_boundary_conditions(bc_dir: Path, state: Stage3State,
                              stage2_result: Dict[str, Any] | None) -> None:
    
    content_U = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}}

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 {state.v_z0});

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform (0 0 {state.v_z0});
    }}
    outlet
    {{
        type            zeroGradient;
    }}
    wall
    {{
        type            noSlip;
    }}
    wedgeFront
    {{
        type            wedge;
    }}
    wedgeBack
    {{
        type            wedge;
    }}
    axis
    {{
        type            empty;
    }}
}}
"""
    (bc_dir / "U").write_text(content_U)
    
    content_p = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}}

dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform {state.P_in};

boundaryField
{{
    inlet
    {{
        type            zeroGradient;
    }}
    outlet
    {{
        type            fixedValue;
        value           uniform {state.P_in};
    }}
    wall
    {{
        type            zeroGradient;
    }}
    wedgeFront
    {{
        type            wedge;
    }}
    wedgeBack
    {{
        type            wedge;
    }}
    axis
    {{
        type            empty;
    }}
}}
"""
    (bc_dir / "p").write_text(content_p)
    
    content_T = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      T;
}}

dimensions      [0 0 0 1 0 0 0];

internalField   uniform {state.T_in};

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform {state.T_in};
    }}
    outlet
    {{
        type            zeroGradient;
    }}
    wall
    {{
        type            fixedGradient;
        gradient        uniform {state._q_flux_auto / 25.0};  // q''/k
    }}
    wedgeFront
    {{
        type            wedge;
    }}
    wedgeBack
    {{
        type            wedge;
    }}
    axis
    {{
        type            empty;
    }}
}}
"""
    (bc_dir / "T").write_text(content_T)


def import_openfoam_dispersion(case_dir: str | Path,
                              time_dir: str = "1000") -> Dict[str, np.ndarray]:
    case_dir = Path(case_dir)
    time_path = case_dir / time_dir
    
    result = {
        "z": None,
        "D_ax": None,
        "U": None,
        "T": None,
    }
    
    print(f"Would read CFD data from: {time_path}")
    print("Note: OpenFOAM parsing not implemented in this version.")
    print("Use foamToVTK and VTK/ParaView for post-processing.")
    
    return result


def print_stage3_results(result: Dict[str, Any]) -> None:
    z = result["z"]
    Y = result["Y"]
    T = result["T"]
    P = result["P"]
    D_ax = result["D_ax"]
    species = result["species"]
    
    spec_idx = {s: i for i, s in enumerate(species)}
    
    print(f"\nStage-3 PFR with Axial Dispersion")
    print("=" * 60)
    print(f"Converged: {result['converged']} ({result['n_iter']} iterations)")
    print(f"z range: {z[0]:.3f} to {z[-1]:.3f} m")
    print()
    
    print("Peclet numbers:")
    print(f"  Pe_mass: {result['Pe_mass']:.1f}")
    print(f"  Pe_heat: {result['Pe_heat']:.1f}")
    print()
    
    print("Axial dispersion coefficient:")
    print(f"  D_ax:    {D_ax[0]*1e4:.4f} -> {D_ax[-1]*1e4:.4f} cm2/s")
    print(f"  D_ax_avg: {np.mean(D_ax)*1e4:.4f} cm2/s")
    print()
    
    print("Axial profiles (inlet -> outlet):")
    print(f"  T:       {T[0]:.1f} -> {T[-1]:.1f} K")
    print(f"  P:       {P[0]/1e5:.4f} -> {P[-1]/1e5:.4f} bar")
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


def compare_with_plug_flow(dispersion_result: Dict[str, Any],
                          plug_flow_result: Any) -> Dict[str, float]:
    spec_idx = get_species_index()
    i_eth = spec_idx.get("C2H6", 0)
    i_ene = spec_idx.get("C2H4", 0)
    
    Y_disp = dispersion_result["Y"]
    conv_disp = 1.0 - Y_disp[i_eth, -1] / Y_disp[i_eth, 0]
    sel_disp = (Y_disp[i_ene, -1] - Y_disp[i_ene, 0]) / (Y_disp[i_eth, 0] - Y_disp[i_eth, -1])
    
    if hasattr(plug_flow_result, 'y'):
        n_spec = len(get_species_list())
        Y_pf_in = plug_flow_result.y[:n_spec, 0]
        Y_pf_out = plug_flow_result.y[:n_spec, -1]
    else:
        Y_pf = plug_flow_result["Y"]
        Y_pf_in = Y_pf[:, 0]
        Y_pf_out = Y_pf[:, -1]
    
    conv_pf = 1.0 - Y_pf_out[i_eth] / Y_pf_in[i_eth]
    sel_pf = (Y_pf_out[i_ene] - Y_pf_in[i_ene]) / (Y_pf_in[i_eth] - Y_pf_out[i_eth])
    
    comparison = {
        "conv_plug_flow": conv_pf,
        "conv_dispersion": conv_disp,
        "conv_diff_pct": (conv_disp - conv_pf) / conv_pf * 100,
        "sel_plug_flow": sel_pf,
        "sel_dispersion": sel_disp,
        "Pe_mass": dispersion_result["Pe_mass"],
        "Pe_heat": dispersion_result["Pe_heat"],
    }
    
    print("\nComparison: Plug Flow vs Dispersion Model")
    print("=" * 50)
    print(f"Peclet number (mass): {comparison['Pe_mass']:.1f}")
    print(f"Peclet number (heat): {comparison['Pe_heat']:.1f}")
    print()
    print(f"Ethane conversion:")
    print(f"  Plug flow:   {conv_pf*100:.2f}%")
    print(f"  Dispersion:  {conv_disp*100:.2f}%")
    print(f"  Difference:  {comparison['conv_diff_pct']:+.2f}%")
    print()
    print(f"C2H4 selectivity:")
    print(f"  Plug flow:   {sel_pf*100:.2f}%")
    print(f"  Dispersion:  {sel_disp*100:.2f}%")
    
    return comparison


if __name__ == "__main__":
    print("Stage-3 PFR with Axial Dispersion - Ethane Pyrolysis")
    print("=" * 60)
    
    state = Stage3State(
        Nz=100,
        use_taylor_aris=False,
    )
    
    print(f"Tube: D={state.D*1000:.0f} mm, L={state.L:.1f} m")
    print(f"Inlet: T={state.T_in:.0f} K, P={state.P_in/1e5:.1f} bar")
    print(f"Heat flux (auto): q''={state._q_flux_auto/1e3:.2f} kW/m2")
    print()
    
    result = solve_dispersion_fd(state, feed=FEED_FORMULA, verbose=True)
    print_stage3_results(result)
    
    print("\n" + "=" * 60)
    print("Generating OpenFOAM case...")
    generate_openfoam_case(state, ROOT / "openfoam_case")
