# DWSIM_PFR: Complete Guide - Plug-Flow Reactor Models for Ethane Pyrolysis

## Table of Contents

1. [Overview](#overview)
2. [Installation & Setup](#installation--setup)
3. [Stage-1: Ideal 1-D PFR](#stage-1-ideal-1-d-pfr)
4. [Stage-2: 2-D Heat Transfer Model](#stage-2-2-d-heat-transfer-model)
5. [Stage-3: Axial Dispersion & CFD Validation](#stage-3-axial-dispersion--cfd-validation)
6. [Kinetics & Thermophysical Properties](#kinetics--thermophysical-properties)
7. [Model Comparison & Selection Guide](#model-comparison--selection-guide)
8. [Troubleshooting](#troubleshooting)
9. [References](#references)

---

## Overview

**DWSIM_PFR** is a comprehensive multi-stage Plug-Flow Reactor (PFR) model for thermal cracking of ethane to ethylene. The package implements three progressively sophisticated models:

- **Stage-1**: Ideal 1-D PFR with constant heat flux (fast, preliminary design)
- **Stage-2**: 2-D heat transfer with wall conduction and radiation (industrial furnaces)
- **Stage-3**: Axial dispersion with CFD validation tools (non-ideal mixing)

### Key Features

- **18-reaction mechanism** from reduced Glarborg-Sendt kinetics
- **NASA-7 thermochemistry** from RMG database
- **Coupled species/energy/momentum** equations
- **Reversible kinetics** with detailed balance ($k_r = k_f / K_c$)
- **Thermodynamically consistent** reaction enthalpies and Gibbs energy
- **Industrial furnace modeling** with wall conduction and radiation
- **OpenFOAM integration** for CFD validation
- **Robust numerical solvers** (BDF for stiff kinetics)

### Typical Applications

- Ethylene plant design and optimization
- Furnace tube sizing and heat flux estimation
- Conversion and selectivity prediction
- Process sensitivity analysis
- CFD model validation

---

## Installation & Setup

### Requirements

- Python 3.8+
- numpy >= 1.20.0
- scipy >= 1.7.0

### Installation

```bash
pip install numpy scipy
```

Or using conda:

```bash
conda install numpy scipy
```

### Quick Start

```python
from DWSIM_PFR import run_pfr, print_results

# Run Stage-1 with default conditions
sol = run_pfr()
print_results(sol)
```

See `SETUP.md` for detailed installation instructions.

---

## Stage-1: Ideal 1-D PFR

### Model Description

Stage-1 implements an ideal plug-flow reactor solving coupled species, energy, and pressure equations. It features:
- **Reversible kinetics**: Net rates $\Omega_{net} = \Omega_f - \Omega_r$ with $k_r$ via detailed balance
- **Thermodynamic consistency**: Enthalpy and Gibbs energy via NASA-7 polynomials
- **PPR78 Equation of State**: Non-ideal gas behavior with compressibility factor and fugacity corrections
- Perfect plug flow (no radial gradients)
- Constant heat flux along the reactor
- No axial dispersion

### Governing Equations

**Species balance** (molar fractions):
```
dY_i/dz = (1/v_z) * (1/C_total) * Σ_j(ν_ij * Ω_net,j)
```

Where:
- Y_i: mole fraction of species i
- v_z: axial velocity [m/s]
- C_total: total molar concentration [mol/m³]
- ν_ij: stoichiometric coefficient
- Ω_net,j = Ω_f,j - Ω_r,j: net reaction rate [mol/m³/s]

**Energy balance**:
```
dT/dz = (1/(ρ * Cp_mass * v_z)) * [Q_rxn + 4*q''/D]
```

Where:
- Q_rxn = -Σ_j(ΔH_j * Ω_net,j): volumetric reaction heat [W/m³]
- ΔH_j: thermodynamically consistent reaction enthalpy [J/mol]
- q'': wall heat flux [W/m²]
- D: tube diameter [m]
- ρ: density [kg/m³]
- Cp_mass: mass heat capacity [J/(kg·K)]

**Pressure drop** (Darcy-Weisbach with Blasius friction):
```
dP/dz = -4*f/D * ρ*v_z²/2
f = 0.0791 * Re^(-0.25)
```

**Velocity** from continuity:
```
v_z = v_z0 * (P0/P) * (T/T0)
```

### Usage

#### Basic Usage

```python
from DWSIM_PFR import run_pfr, print_results

# Default conditions: 1300 K, 2 bar, 10 m tube, pure ethane
sol = run_pfr()
print_results(sol)
```

#### Custom Conditions

```python
sol = run_pfr(
    L=10.0,           # reactor length [m]
    T_in=1000.0,      # inlet temperature [K]
    P_in=2e5,         # inlet pressure [Pa]
    v_z0=1.0,         # superficial velocity [m/s]
    D=0.1,            # tube diameter [m]
    q_flux=50000.0,   # wall heat flux [W/m²]
    feed={"C2H6": 1.0},  # feed composition
    max_step=0.01,    # max integration step [m]
)
print_results(sol)
```

#### Accessing Results

```python
# Solution object attributes
z = sol.t                    # axial positions [m]
Y = sol.y[:n_spec, :]       # mole fractions (n_species, n_points)
T = sol.y[n_spec, :]        # temperature [K]
P = sol.y[n_spec+1, :]      # pressure [Pa]

# Species names
species = sol.species        # list of species formula tokens

# Conversion calculation
from DWSIM_PFR import get_species_index
spec_idx = get_species_index()
i_eth = spec_idx["C2H6"]
conv = 1.0 - Y[i_eth, -1] / Y[i_eth, 0]
print(f"Ethane conversion: {conv*100:.1f}%")
```

### Heat Flux Estimation

The model can auto-estimate heat flux for near-isothermal operation:

```python
from DWSIM_PFR.pfr_ideal import estimate_heat_flux, T_SET, P_IN

q_flux = estimate_heat_flux(
    T=T_SET,              # target temperature [K]
    P=P_IN,               # pressure [Pa]
    target_conversion=0.9  # target conversion
)
print(f"Estimated heat flux: {q_flux/1e3:.2f} kW/m²")
```

### Strengths

1. **Fast computation**: Typically < 1 second for 10 m tube
2. **Robust solver**: BDF method handles stiff kinetics
3. **Reversible kinetics**: Rigorous treatment of chemical equilibrium via detailed balance
4. **Thermodynamic consistency**: NASA-7 polynomials ensure accurate $K_c$ and $\Delta H_{rxn}$
5. **Proper unit conversion**: Pre-exponential factors correctly converted
6. **Third-body reactions**: Handled with total concentration [M] = P/(Z·R·T) using PPR78 EOS
7. **Property database**: Uses resolved species keys for NASA-7 lookups
8. **Non-ideal gas**: PPR78 EOS with fugacity-corrected equilibrium constants

### Limitations

1. **Linear mixing rules**: Viscosity and thermal conductivity use simple molar-fraction weighting instead of Wilke/Wassiljewa rules
2. **Constant q''**: Heat flux uniform along reactor in Stage-1
3. **No radial gradients**: Ideal plug flow assumption (significant in large D tubes)
4. **No dispersion**: Perfect plug flow (no back-mixing) in Stage-1/Stage-2
5. **Property dependency**: Accuracy is strictly bound by `species_properties.json` data quality
6. **Radical critical properties**: Estimated from parent molecules (uncertainty in EOS for radicals)

### Typical Results

For 10 m tube, 100 mm diameter, 1000 K inlet, 2 bar:
- Outlet temperature: ~1097 K (with auto heat flux ~5 kW/m²)
- Ethane conversion: ~53% (depends on heat flux and residence time)
- Pressure drop: < 0.1 mbar (negligible)

---

## Stage-2: 2-D Heat Transfer Model

### Model Description

Stage-2 extends Stage-1 with **coupled 2-D heat transfer** through the tube wall and **gas-phase radiation**. This enables realistic simulation of industrial ethane cracking furnaces where heat is supplied from an external firebox.

### Key Differences from Stage-1

| Feature | Stage-1 | Stage-2 |
|---------|---------|---------|
| Heat flux | Constant q''(z) | Variable q''(z) from wall-gas T difference |
| Wall model | None | Radial conduction through tube wall |
| Radiation | None | WSGG model for gas + wall radiation |
| Boundary conditions | Simple | Coupled inner (gas-wall) + outer (furnace-wall) |
| Computational cost | Fast (< 1 s) | Moderate (~10-30 s) |

### Physical Model

#### 1. Tube Wall Conduction

The tube wall is modeled as a **cylindrical shell** with radial conduction:

```
ρ_w * Cp_w * ∂T_w/∂t = (1/r) * ∂/∂r(k_w * r * ∂T_w/∂r)
```

For steady-state operation:
```
(1/r) * d/dr(k_w * r * dT_w/dr) = 0
```

**Solution** (constant k_w):
```
T_w(r) = T_outer + (q''_inner * r_inner / k_w) * ln(r_outer / r)
```

#### 2. Inner Boundary: Gas-Wall Interface

Heat flux from wall to gas has **convective and radiative** components:

```
q''_inner = h_i * (T_wall - T_gas) + ε_w * σ * (T_wall^4 - T_gas^4)
```

**Convective coefficient h_i:**
- Computed from **Gnielinski correlation** for turbulent pipe flow:
  ```
  Nu = (f/8)(Re - 1000)Pr / [1 + 12.7(f/8)^0.5(Pr^(2/3) - 1)]
  h_i = Nu * k_gas / D
  ```
- Valid for: 3000 < Re < 5×10⁶, 0.5 < Pr < 2000

**Radiative component:**
- Wall emissivity: ε_w = 0.85 (oxidized Inconel)
- Stefan-Boltzmann: σ = 5.67×10⁻⁸ W/(m²·K⁴)
- Gas radiation is weak for dry pyrolysis (no H₂O/CO₂)

#### 3. Outer Boundary: Furnace-Wall Interface

Heat flux from furnace to outer wall:

```
q''_outer = h_o * (T_furnace - T_wall) + ε_w * σ * (T_furnace^4 - T_wall^4)
```

**Typical values:**
- h_o = 60 W/(m²·K) (external convection in firebox)
- T_furnace = 1300 K (radiant section temperature)
- Radiative term dominates at high temperatures (T⁴ scaling)

#### 4. Gas-Phase Radiation (WSGG Model)

**Weighted Sum of Gray Gases (WSGG)** model for participating species:

```
κ_eff = Σ_i a_i(T) * κ_i * P_participating
ε_g = Σ_i a_i(T) * [1 - exp(-κ_i * P * L_beam)]
```

Where:
- κ_i: absorption coefficients for gray gases
- a_i(T): temperature-dependent weights
- P_participating: partial pressure of H₂O + CO₂
- L_beam ≈ 0.9×D (mean beam length for cylinder)

**For ethane pyrolysis:**
- No steam dilution → P_participating ≈ 0
- Gas is optically thin (minimal self-absorption)
- Model returns minimal κ ≈ 0.01 m⁻¹

### Solution Algorithm: Operator-Split Iteration

The coupled system is solved using **operator-splitting** with under-relaxation:

1. **Initialize**: Set initial gas T profile and wall temperatures
2. **Compute heat flux profile**: From current wall-gas temperature differences
3. **Integrate species/T/P ODEs**: Use Stage-1 ODE system with variable q''(z)
4. **Update wall temperatures**: From heat balance (q_out = q_wall = q_in)
5. **Check convergence**: Repeat with under-relaxation until ΔT < 1 K

**Under-relaxation:**
```
T_gas = α * T_gas_new + (1 - α) * T_gas_old
q'' = α * q''_new + (1 - α) * q''_old
```
Where α = 0.3 (default) prevents oscillation.

### Usage

#### Basic Usage

```python
from DWSIM_PFR import Stage2State, solve_stage2, print_stage2_results

# Configure reactor
state = Stage2State(
    D_inner=0.1,           # 100 mm inner diameter
    wall_thickness=0.008,  # 8 mm wall
    L=10.0,                # 10 m length
    T_in=1000.0,           # inlet gas temperature [K]
    P_in=2e5,              # inlet pressure [Pa]
    v_z0=1.0,              # inlet velocity [m/s]
    T_furnace=1300.0,      # furnace temperature [K]
    h_outer=60.0,          # outer convection [W/(m²·K)]
    Nz=50,                 # axial nodes
    Nr_wall=5,             # radial wall nodes
)

# Solve
result = solve_stage2(
    state,
    feed={"C2H6": 1.0},
    max_outer_iter=50,
    tol_T=1.0,             # convergence tolerance [K]
    relax=0.3,             # under-relaxation factor
    verbose=True
)

# Display results
print_stage2_results(result)
```

#### Accessing Results

```python
# Solution arrays
z = result["z"]                    # axial positions [m]
T_gas = result["T_gas"]            # gas temperature [K]
T_wall_inner = result["T_wall_inner"]  # inner wall T [K]
T_wall_outer = result["T_wall_outer"]  # outer wall T [K]
q_inner = result["q_inner"]        # heat flux [W/m²]
Y = result["Y"]                    # mole fractions (n_species, Nz)
P = result["P"]                    # pressure [Pa]

# Convergence info
converged = result["converged"]
n_iter = result["n_iter"]
```

#### Custom Wall Material

```python
from DWSIM_PFR import WallMaterial, Stage2State

# Define custom material
material = WallMaterial(
    name="Custom Alloy",
    rho=8000.0,       # kg/m³
    k=30.0,           # W/(m·K)
    cp=450.0,         # J/(kg·K)
    emissivity=0.80,  # surface emissivity
)

# Use in Stage2State
state = Stage2State(
    D_inner=0.1,
    wall_thickness=0.008,
    L=10.0,
    wall_material=material,
    # ... other parameters
)
```

### Material Properties

#### Wall Material: Inconel 800H (Default)

Typical for ethylene cracking furnaces:
- **Density**: ρ = 7940 kg/m³
- **Thermal conductivity**: k = 25 W/(m·K) at ~1000 K
- **Specific heat**: Cp = 500 J/(kg·K)
- **Emissivity**: ε = 0.85 (oxidized surface)
- **Thermal diffusivity**: α = k/(ρ·Cp) ≈ 6.3×10⁻⁶ m²/s

#### Gas Properties

**Thermal conductivity** (Eucken correlation):
```
k_gas ≈ 0.02 + 7×10⁻⁵ * T  [W/(m·K)]
```

**Viscosity**: From fitted ln-polynomial in `species_properties.json`

**Heat capacity**: NASA-7 polynomials from RMG database

### Typical Results

For a 10 m tube, 100 mm diameter, 8 mm wall:

**Operating conditions:**
- Inlet: T = 1000 K, P = 2 bar, v_z = 1 m/s
- Furnace: T = 1300 K, h_outer = 60 W/(m²·K)

**Results:**
- Outlet gas temperature: **1240 K** (vs 1097 K in Stage-1)
- Ethane conversion: **90%** (vs 53% in Stage-1)
- Heat flux profile: **46 → 11 kW/m²** (decreasing axially)
- Wall temperatures: **1178 → 1269 K** (inner), **1193 → 1273 K** (outer)

**Why higher conversion?**
1. **Higher heat input**: Average q'' ≈ 17.4 kW/m² vs 5.0 kW/m² in Stage-1
2. **Higher temperatures**: 1240 K vs 1097 K → faster reaction rates (Arrhenius)
3. **Variable flux**: High flux at inlet where ethane concentration is highest

### Numerical Considerations

**Convergence:**
- Typical iterations: 15-25 for ΔT < 1 K
- Under-relaxation: α = 0.3 prevents oscillation
- Tolerance: 1 K on gas temperature

**Stability:**
- Heat flux clamping: Limited to ±200 kW/m² (magnitude clamp; negative q'' allowed for gas→wall cooling)
- Temperature clamping: T_gas ∈ [300, 2500] K
- No forced heat-flow direction: the sign of q'' follows from (T_wall − T_gas) and (T_furnace − T_wall)

**Discretization:**
- Axial nodes: Nz = 50-100 (default 50)
- Radial wall nodes: Nr = 5 (sufficient for thin wall)
- ODE solver: BDF method with adaptive step size

### Limitations

1. **Steady-state wall**: No transient wall response (quasi-steady assumption)
2. **1-D gas flow**: Still plug flow, no radial gradients in gas phase
3. **Simplified WSGG**: Generic H₂O/CO₂ coefficients, not tuned for hydrocarbons
4. **No coking model**: No carbon deposition on inner wall
5. **No steam dilution**: Pure ethane feed assumed (industrial uses steam)
6. **Constant wall properties**: k, Cp, ρ assumed constant (weak temperature dependence)

### Comparison with Industrial Data

**Typical industrial ethane cracker:**
- Tube diameter: 80-120 mm
- Wall thickness: 6-10 mm
- Furnace temperature: 1200-1400 K
- Inlet temperature: 800-1000 K
- Outlet temperature: 1100-1300 K
- Conversion: 60-90% (depending on residence time)

**Stage-2 predictions:**
- Outlet temperature: 1240 K ✓ (within range)
- Conversion: 90% ✓ (high end, typical for long tubes)
- Heat flux: 11-46 kW/m² ✓ (typical range: 10-50 kW/m²)

---

## Stage-3: Axial Dispersion & CFD Validation

### Model Description

Stage-3 extends Stage-2 with **axial dispersion** effects and provides tools for **CFD validation** using OpenFOAM. This accounts for non-ideal mixing caused by velocity gradients and turbulent diffusion.

### Key Features

1. **Axial dispersion model**: Adds second-derivative terms to species and energy equations
2. **Danckwerts boundary conditions**: Proper inlet/outlet BCs for dispersion
3. **Dispersion correlations**: Taylor-Aris (laminar) and Levenspiel (turbulent)
4. **OpenFOAM case generation**: Complete CFD setup for validation

### Governing Equations with Dispersion

**Species equation:**
```
dY_i/dz = (1/v_z) * Σ_j(ν_ij * Ω_j) + (D_ax/v_z) * d²Y_i/dz²
```

**Energy equation:**
```
dT/dz = (1/(ρ*Cp*v_z)) * [Q_rxn + 4q''/D] + (alpha_ax/v_z) * d²T/dz²
```

Where:
- D_ax: axial mass dispersion coefficient [m²/s]
- alpha_ax: axial thermal dispersivity [m²/s]

### Danckwerts Boundary Conditions

**Inlet (z=0):**
```
v_z * (Y_in - Y(0)) = -D_ax * dY/dz |_{z=0}
```

This ensures continuity of flux at the inlet.

**Outlet (z=L):**
```
dY/dz |_{z=L} = 0
```

Zero gradient (fully developed flow).

### Dispersion Coefficient Correlations

#### Taylor-Aris (Laminar Flow)

For Re < 2300:
```
D_ax = D_mol + (v_z² * D²) / (192 * D_mol)
```

Where D_mol is molecular diffusivity.

#### Turbulent (Levenspiel Correlation)

For Re > 2300:
```
1/Pe_ax = 0.3/Re + 0.5/(1 + 3.8/Re)
D_ax = (1/Pe_ax) * v_z * D
```

Where Pe_ax = v_z * D / D_ax is the axial Peclet number.

### Physical Background

#### Peclet Number

The Peclet number characterizes the ratio of convection to dispersion:
```
Pe = v_z * L / D_ax
```

- **Pe > 100**: Plug flow dominates, dispersion is a small correction
- **Pe ~ 10-100**: Moderate dispersion effects
- **Pe < 10**: Strong back-mixing, approaches CSTR behavior

#### Effect on Conversion

Axial dispersion typically **reduces conversion** compared to ideal plug flow because:
1. Back-mixing dilutes reactant concentrations
2. Products are transported upstream, shifting equilibrium
3. Temperature gradients are smoothed

For high Pe (> 50), the effect is usually < 2% on conversion.

### Solution Algorithm

#### Two-Step Approach

1. **Initialize with Stage-1 plug flow**: Uses implicit BDF solver (numerically stable)
2. **Apply dispersion corrections**: Iteratively add second-derivative corrections with under-relaxation

This approach is efficient for high Peclet flows where dispersion is a perturbation.

#### Numerical Details

- **Discretization**: Central differences for d²/dz², upwind for d/dz
- **Under-relaxation**: Default α = 0.3 for stability
- **Grid**: Uniform spacing, typically 100 nodes for 10 m tube
- **Convergence**: dY_max < 1e-5, dT_max < 0.5 K

### Usage

#### Basic Usage

```python
from DWSIM_PFR import Stage3State, solve_dispersion_fd, print_stage3_results

# Configure reactor (heat flux auto-estimated)
state = Stage3State(
    D=0.1,              # 100 mm tube diameter
    L=10.0,             # 10 m length
    T_in=1000.0,        # inlet temperature [K]
    P_in=2e5,           # inlet pressure [Pa]
    v_z0=1.0,           # inlet velocity [m/s]
    Nz=100,             # 100 axial nodes
    use_taylor_aris=False,  # use turbulent correlation
    # q_flux_constant=None uses auto-estimate
)

# Solve
result = solve_dispersion_fd(
    state,
    feed={"C2H6": 1.0},
    max_iter=200,
    tol=1e-5,
    relax=0.3,
    verbose=True
)

# Display results
print_stage3_results(result)
```

#### Using Stage-2 Heat Flux Profile

```python
from DWSIM_PFR import Stage2State, solve_stage2, Stage3State, solve_dispersion_fd

# First solve Stage-2
state2 = Stage2State(D_inner=0.1, L=10.0, T_in=1000.0, T_furnace=1300.0)
result2 = solve_stage2(state2)

# Use heat flux profile in Stage-3
state3 = Stage3State(D=0.1, L=10.0, T_in=1000.0)
state3.q_flux_profile = result2["q_inner"]

result3 = solve_dispersion_fd(state3)
```

#### Comparing with Plug Flow

```python
from DWSIM_PFR import (
    run_pfr, solve_dispersion_fd, Stage3State, compare_with_plug_flow
)

# Solve both models
pf_result = run_pfr(L=10.0, T_in=1000.0)
state = Stage3State(D=0.1, L=10.0, T_in=1000.0)
disp_result = solve_dispersion_fd(state)

# Compare
comparison = compare_with_plug_flow(disp_result, pf_result)
```

### OpenFOAM Case Generation

#### Generate Case

```python
from DWSIM_PFR import Stage3State, generate_openfoam_case, solve_stage2

# Generate case with Stage-2 reaction sources
state = Stage3State(D=0.1, L=10.0, T_in=1000.0)
result2 = solve_stage2(Stage2State(...))

generate_openfoam_case(state, "./openfoam_case", stage2_result=result2)
```

#### Case Structure

```
openfoam_case/
├── 0/
│   ├── U           # Velocity BC
│   ├── p           # Pressure BC
│   └── T           # Temperature BC
├── constant/
│   ├── thermophysicalProperties
│   ├── transportProperties
│   ├── turbulenceProperties
│   └── reactionSources.json  # Tabulated reaction rates
└── system/
    ├── blockMeshDict
    ├── controlDict
    ├── fvSchemes
    └── fvSolution
```

#### Running OpenFOAM

```bash
cd openfoam_case
blockMesh
reactingFoam
```

#### Post-Processing

```bash
foamToVTK
# Open in ParaView for visualization
```

### Configuration Options

#### Stage3State Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `D` | 0.1 | Tube diameter [m] |
| `L` | 10.0 | Tube length [m] |
| `T_in` | 1000.0 | Inlet temperature [K] |
| `P_in` | 2e5 | Inlet pressure [Pa] |
| `v_z0` | 1.0 | Inlet velocity [m/s] |
| `Nz` | 100 | Number of axial nodes |
| `use_taylor_aris` | False | Use laminar (True) or turbulent (False) correlation |
| `D_ax_constant` | None | Optional constant D_ax [m²/s] |
| `q_flux_constant` | None | Heat flux [W/m²], None = auto-estimate |
| `q_flux_profile` | None | Optional q''(z) array |

### Typical Results

For 10 m tube, 100 mm diameter, 1000 K inlet (auto heat flux ~7.5 kW/m²):

**Temperature profile:**
- Inlet: 1000 K
- Outlet: ~1137 K (temperature rises due to heat input exceeding reaction endotherm initially)

**Dispersion parameters:**
- D_ax: 250-500 cm²/s (varies with T and composition)
- Pe_mass: 300-400 (high, near plug flow)
- Pe_heat: 250-350

**Effect on conversion:**
- Plug flow conversion: ~65%
- With dispersion: ~65% (< 0.5% difference for Pe > 100)

**Observation:**
For typical industrial conditions (Re > 10000, Pe > 100), axial dispersion has minimal impact on overall conversion. The main effect is smoothing of concentration and temperature gradients.

**Important Notes on Heat Flux:**
- The `q_flux_constant` parameter defaults to `None`, which triggers auto-estimation
- Auto-estimate uses the reaction enthalpy and target conversion to compute heat input
- For higher heat flux (industrial furnaces), use Stage-2 with coupled wall model
- Manually setting very high heat flux (> 50 kW/m²) without proper reaction heat consumption leads to unphysical temperature rise

### Limitations

1. **1-D approximation**: No radial gradients in dispersion model
2. **Steady-state only**: No transient effects
3. **Simplified correlations**: D_ax from empirical fits, not first-principles
4. **No coking effects**: No fouling or carbon deposition
5. **OpenFOAM setup**: Basic case generation, may need manual tuning

### Validation Workflow

1. **Generate OpenFOAM case** from Python Stage-3
2. **Run CFD simulation** with full 3D or axisymmetric mesh
3. **Extract effective D_ax(z)** from CFD results:
   - Compute flux: J = -D_eff * dC/dz + v*C
   - Solve for D_eff at each z location
4. **Import D_ax(z) into Python** using `import_openfoam_dispersion()`
5. **Re-solve Stage-3** with CFD-derived D_ax
6. **Compare** conversion and profiles between Python and CFD

---

## Kinetics & Thermophysical Properties

### Reaction Mechanism

The model uses **18 reactions** from the reduced Glarborg-Sendt mechanism for ethane pyrolysis:

- Primary cracking: C2H6 → C2H4 + H2
- Secondary reactions: C2H4 → C2H2 + H2, C2H6 → CH4 + CH4, etc.
- Radical chain reactions: H + C2H6 → H2 + C2H5, etc.

**Reaction format:**
```
k = A * (T/298.15)^n * exp(-Ea/(R*T))
```

**Unit conversion:**
- Pre-exponential factors converted from molecule/cm³ basis to mol/m³
- Activation energies in J/mol
- Third-body reactions (M) handled with [M] = P/(Z·R·T) using PPR78 EOS

### Thermophysical Properties

#### Heat Capacity and Enthalpy

**NASA-7 polynomials** from RMG database:
```
Cp(T) = R * (a1 + a2*T + a3*T² + a4*T³ + a5*T⁴)
H(T) = R*T * (a1 + a2*T/2 + a3*T²/3 + a4*T³/4 + a5*T⁴/5 + a6/T)
```

Stored in `species_properties.json` with temperature ranges (low/high).

#### Viscosity

**Fitted ln-polynomial correlations:**
```
ln(μ) = c0*ln³(T) + c1*ln²(T) + c2*ln(T) + c3
```

Fitted to tabulated data at reference pressure.

#### Thermal Conductivity

**Eucken correlation** (simplified):
```
k_gas ≈ 0.02 + 7×10⁻⁵ * T  [W/(m·K)]
```

#### Mixture Properties

**Simple mixing rules:**
- Heat capacity: Cp_mix = Σ_i(Y_i * Cp_i) + Cp_residual (EOS)
- Viscosity: μ_mix = Σ_i(Y_i * μ_i)
- Density: ρ = P * MW_mix / (Z * R * T) (PPR78 EOS)

#### PPR78 Equation of State

The models use the **Predictive Peng-Robinson 1978 (PPR78)** equation of state for non-ideal gas behavior:

```
Z³ - (1-B)Z² + (A - 3B² - 2B)Z - (AB - B² - B³) = 0
```

Where:
- Z = compressibility factor (Z=1 for ideal gas)
- A, B = reduced EOS parameters from critical properties

**Features:**
- Temperature-dependent binary interaction parameters via group contribution
- Fugacity coefficients for equilibrium calculations
- Residual enthalpy, entropy, and heat capacity

**Critical properties database** (`eos_ppr78.py`):
- Stable species: from ChemSep/literature (Tc, Pc, ω)
- Radicals: estimated from parent molecules

**Fugacity corrections** for equilibrium:
- Kc_corrected = Kc / K_φ
- K_φ = Π(φ_products)^ν / Π(φ_reactants)^|ν|

**Note:** At typical pyrolysis conditions (T > 800K, P < 5 bar), non-ideality corrections are small (<1%), but the framework supports higher pressures where deviations become significant.

### Data Files

| File | Description |
|------|-------------|
| `reduced_by_GS_mechanism_reactions.csv` | 18-reaction mechanism with Arrhenius parameters |
| `species_properties.json` | NASA-7 coefficients, transport properties, molecular weights |
| `ethane_mechanism_complete.json` | Full RMG mechanism (reference) |
| `ethane_pyrolysis_rmg.json` | RMG kinetics export (reference) |

---

## Model Comparison & Selection Guide

### When to Use Each Stage

#### Stage-1: Ideal PFR

**Use when:**
- Preliminary design and sizing
- Fast parametric studies
- Validation of kinetics
- Educational purposes

**Advantages:**
- Fast computation (< 1 s)
- Simple to use
- Robust numerical solver

**Limitations:**
- Constant heat flux assumption
- No wall model
- No dispersion effects

#### Stage-2: Heat Transfer Model

**Use when:**
- Industrial furnace design
- Heat flux optimization
- Wall temperature prediction
- Realistic operating conditions

**Advantages:**
- Variable heat flux from wall-gas coupling
- Wall conduction model
- Radiative transfer included
- Matches industrial operation

**Limitations:**
- Slower computation (~10-30 s)
- More parameters to configure
- Still assumes plug flow

#### Stage-3: Axial Dispersion

**Use when:**
- Non-ideal mixing effects important
- CFD validation needed
- Low Peclet number flows
- Research/validation studies

**Advantages:**
- Accounts for back-mixing
- OpenFOAM integration
- Can use CFD-derived D_ax

**Limitations:**
- Slowest computation (~30-60 s)
- Most complex setup
- Dispersion effects small for high Pe

### Performance Comparison

| Metric | Stage-1 | Stage-2 | Stage-3 |
|--------|---------|---------|---------|
| Computation time | < 1 s | 10-30 s | 30-60 s |
| Memory usage | Low | Moderate | Moderate |
| Accuracy (high Pe) | Good | Excellent | Excellent |
| Accuracy (low Pe) | Good | Good | Excellent |
| Industrial realism | Low | High | High |

### Typical Results Comparison

For 10 m tube, 100 mm diameter, 1000 K inlet:

| Result | Stage-1 | Stage-2 | Stage-3 |
|--------|---------|---------|---------|
| T_out [K] | 1097 | 1240 | 1137 |
| Conversion [%] | 53 | 90 | 65 |
| Reversibility | Net rates | Net rates | Net rates |
| Heat flux [kW/m²] | 5 (const) | 11-46 (var) | 7.5 (auto) |
| Wall T [K] | N/A | 1178-1269 | N/A |
| Pe number | N/A | N/A | 316 |

**Note:** Results depend strongly on heat flux and operating conditions.

---

## Equilibrium & Reversibility Model

The model now incorporates a rigorous thermodynamic equilibrium module (`equilibrium.py`) to handle reversible reactions.

### Detailed Balance

For every reaction $j$, the net rate is computed as:
$$\Omega_{net,j} = \Omega_{f,j} - \Omega_{r,j}$$

The reverse rate constant $k_{r,j}$ is derived from the forward rate constant $k_{f,j}$ and the concentration-based equilibrium constant $K_{c,j}$:
$$k_{r,j} = \frac{k_{f,j}}{K_{c,j}}$$

### Thermodynamic Calculations

1. **Gibbs Free Energy**: $\Delta G_{rxn}^\circ(T) = \sum \nu_i G_i^\circ(T)$
2. **Equilibrium Constant ($K_p$)**: $K_p(T) = \exp\left(-\frac{\Delta G_{rxn}^\circ}{RT}\right)$
3. **Conversion to $K_c$**: $K_c(T) = K_p(T) \cdot \left(\frac{P_{std}}{RT}\right)^{-\Delta n}$

Where $\Delta n$ is the net change in moles (excluding third bodies).

### Harsh Analysis of the Reversible Model

**Strengths:**
- **Closer to Physics**: Eliminates the "infinite conversion" artifact of irreversible models, especially relevant for secondary reactions like $C_2H_2 + H_2 \rightleftharpoons C_2H_4$.
- **Thermodynamic Closure**: Energy balance ($Q_{rxn}$) and kinetics ($\Omega_{net}$) are now coupled through the same NASA-7 property basis, preventing violations of the first and second laws of thermodynamics.

**Weaknesses:**
- **Increased Stiffness**: The presence of both very fast forward and reverse rates near equilibrium can lead to extreme numerical stiffness, challenging the BDF integrator.
- **Data Sensitivity**: A 1% error in Gibbs energy coefficients can lead to orders of magnitude error in $K_c$ due to the exponential dependence, making the model highly sensitive to the RMG database accuracy.
- **Inert Species Handling**: N2 dilution affects equilibrium through partial pressure reduction and fugacity corrections via the PPR78 EOS.

---

## Troubleshooting

### Common Issues

#### 1. Temperature Drops to 300 K

**Symptom:** Outlet temperature clamped at minimum (300 K)

**Cause:** Heat flux too low or reaction endotherm too large

**Solution:**
- Increase heat flux
- Use auto-estimation: `q_flux=None` in Stage3State
- Check energy balance: Q_rxn + 4*q''/D should be positive

#### 2. Temperature Rises Unphysically (> 2500 K)

**Symptom:** Outlet temperature exceeds 2500 K

**Cause:** Heat flux too high without sufficient reaction heat consumption

**Solution:**
- Reduce heat flux
- Use auto-estimation for near-isothermal operation
- Check that reactions are active (not all ethane consumed)

#### 3. Stage-2 Not Converging

**Symptom:** Iterations continue without convergence

**Cause:** Under-relaxation too high or initial guess poor

**Solution:**
- Reduce relaxation factor: `relax=0.2`
- Increase tolerance: `tol_T=2.0`
- Check furnace temperature is reasonable (1200-1400 K)

#### 4. Import Errors

**Symptom:** `ImportError: attempted relative import with no known parent package`

**Cause:** Running script directly instead of as module

**Solution:**
- Run from package: `python -m DWSIM_PFR.pfr_ideal`
- Or ensure proper PYTHONPATH
- Scripts handle both module and standalone execution

#### 5. OpenFOAM Case Generation Fails

**Symptom:** Error writing OpenFOAM files

**Cause:** Directory permissions or path issues

**Solution:**
- Check write permissions
- Use absolute paths
- Ensure output directory exists or can be created

### Numerical Issues

#### Stiff Kinetics

If solver fails with "stiff system" error:
- Reduce `max_step` parameter
- Use BDF method (default)
- Check for very fast reactions (high A, low Ea)

#### Convergence Problems

If iterations don't converge:
- Reduce under-relaxation factor
- Increase tolerance
- Check for unphysical parameter values
- Verify energy balance is reasonable

### Performance Optimization

#### Speed Up Stage-2

- Reduce `Nz` (axial nodes): 50 instead of 100
- Reduce `Nr_wall`: 3 instead of 5
- Increase `tol_T`: 2.0 instead of 1.0
- Use fewer iterations: `max_outer_iter=30`

#### Speed Up Stage-3

- Reduce `Nz`: 50 instead of 100
- Increase `tol`: 1e-4 instead of 1e-5
- Reduce `max_iter`: 100 instead of 200

---

## References

### Stage-1 References

1. **Plug-flow reactor theory**: O. Levenspiel, "Chemical Reaction Engineering," 3rd ed., Wiley (1999)
2. **BDF solver**: L.F. Shampine and M.W. Reichelt, "The MATLAB ODE Suite," SIAM J. Sci. Comput., 18, 1-22 (1997)

### Stage-2 References

1. **Gnielinski correlation**: V. Gnielinski, "New equations for heat and mass transfer in turbulent pipe and channel flow," Int. Chem. Eng., 16, 359-368 (1976)
2. **WSGG model**: T.F. Smith et al., "Evaluation of coefficients for the weighted sum of gray gases model," J. Heat Transfer, 104, 602-608 (1982)
3. **Wall conduction**: Standard cylindrical coordinate heat conduction (see any heat transfer textbook)
4. **Industrial furnaces**: S. Matar and L.F. Hatch, "Chemistry of Petrochemical Processes," 2nd ed., Gulf Publishing (2001)

### Stage-3 References

1. **Taylor dispersion**: G.I. Taylor, "Dispersion of soluble matter in solvent flowing slowly through a tube," Proc. R. Soc. London A, 219, 186-203 (1953)
2. **Aris extension**: R. Aris, "On the dispersion of a solute in a fluid flowing through a tube," Proc. R. Soc. London A, 235, 67-77 (1956)
3. **Turbulent dispersion**: O. Levenspiel, "Chemical Reaction Engineering," 3rd ed., Wiley (1999), Ch. 13
4. **Danckwerts BCs**: P.V. Danckwerts, "Continuous flow systems: Distribution of residence times," Chem. Eng. Sci., 2, 1-13 (1953)

### Kinetics References

1. **Glarborg-Sendt mechanism**: P. Glarborg and K. Sendt, "Reaction mechanism for ethane pyrolysis," Combust. Flame, 123, 452-468 (2000)
2. **NASA polynomials**: B.J. McBride et al., "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species," NASA/TP-2002-211556 (2002)
3. **RMG database**: W.H. Green et al., "RMG - Reaction Mechanism Generator," https://rmg.mit.edu/

---

## File Structure

```
DWSIM_PFR/
├── __init__.py                    # Package initialization
├── pfr_ideal.py                  # Stage-1: Ideal PFR
├── pfr_heattransfer.py           # Stage-2: Heat transfer
├── pfr_dispersion.py             # Stage-3: Axial dispersion
├── equilibrium.py                # Equilibrium calculations
├── eos_ppr78.py                  # PPR78 Equation of State
├── kinetics_parser.py             # CSV mechanism parser
├── props.py                      # Thermophysical properties
├── validate_eos_integration.py   # EOS validation script
├── species_properties.json        # NASA-7 + transport data
├── reduced_by_GS_mechanism_reactions.csv  # Reaction mechanism
├── ethane_mechanism_complete.json # Full RMG mechanism (ref)
├── ethane_pyrolysis_rmg.json     # RMG kinetics export (ref)
├── README.md                     # Basic overview
├── README_STAGE2.md              # Stage-2 technical docs
├── README_STAGE3.md              # Stage-3 technical docs
├── README_COMPLETE.md            # This file
├── SETUP.md                      # Installation guide
└── openfoam_case/                # Generated OpenFOAM case
    ├── 0/                        # Boundary conditions
    ├── constant/                 # Properties and sources
    └── system/                   # Solver settings
```

---

## Quick Reference

### Stage-1 Quick Start

```python
from DWSIM_PFR import run_pfr, print_results
sol = run_pfr(L=10.0, T_in=1000.0, P_in=2e5)
print_results(sol)
```

### Stage-2 Quick Start

```python
from DWSIM_PFR import Stage2State, solve_stage2, print_stage2_results
state = Stage2State(D_inner=0.1, L=10.0, T_in=1000.0, T_furnace=1300.0)
result = solve_stage2(state)
print_stage2_results(result)
```

### Stage-3 Quick Start

```python
from DWSIM_PFR import Stage3State, solve_dispersion_fd, print_stage3_results
state = Stage3State(D=0.1, L=10.0, T_in=1000.0)
result = solve_dispersion_fd(state)
print_stage3_results(result)
```

---

## License & Citation

If you use this code in your research, please cite:

```
DWSIM_PFR: Multi-Stage Plug-Flow Reactor Models for Ethane Pyrolysis
[Your citation information here]
```

---

**End of Complete Guide**

