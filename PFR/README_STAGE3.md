# Stage-3 PFR: Axial Dispersion & CFD Validation

## Overview

Stage-3 extends Stage-2 with **axial dispersion** effects and provides tools for **CFD validation** using OpenFOAM. This accounts for non-ideal mixing caused by velocity gradients and turbulent diffusion.

---

## Key Features

### 1. Axial Dispersion Model

**Species equation with dispersion:**
```
dY_i/dz = (1/v_z) * sum_j(nu_ij * Omega_j) + (D_ax/v_z) * d²Y_i/dz²
```

**Energy equation with dispersion:**
```
dT/dz = (1/(rho*Cp*v_z)) * [Q_rxn + 4q''/D] + (alpha_ax/v_z) * d²T/dz²
```

### 2. Danckwerts Boundary Conditions

**Inlet (z=0):**
```
v_z * (Y_in - Y(0)) = -D_ax * dY/dz |_{z=0}
```

**Outlet (z=L):**
```
dY/dz |_{z=L} = 0
```

### 3. Dispersion Coefficient Correlations

**Taylor-Aris (laminar flow):**
```
D_ax = D_mol + (v_z² * D²) / (192 * D_mol)
```

**Turbulent (Levenspiel correlation):**
```
1/Pe_ax = 0.3/Re + 0.5/(1 + 3.8/Re)
D_ax = (1/Pe_ax) * v_z * D
```

### 4. OpenFOAM Case Generation

Generates complete axisymmetric CFD case for validation:
- Wedge mesh (5° sector)
- k-omega SST turbulence model
- Species transport with tabulated reaction sources
- Heat flux boundary conditions

---

## Physical Background

### Peclet Number

The Peclet number characterizes the ratio of convection to dispersion:
```
Pe = v_z * L / D_ax
```

- **Pe > 100**: Plug flow dominates, dispersion is a small correction
- **Pe ~ 10-100**: Moderate dispersion effects
- **Pe < 10**: Strong back-mixing, approaches CSTR behavior

### Effect on Conversion

Axial dispersion typically **reduces conversion** compared to ideal plug flow because:
1. Back-mixing dilutes reactant concentrations
2. Products are transported upstream, shifting equilibrium
3. Temperature gradients are smoothed

For high Pe (> 50), the effect is usually < 2% on conversion.

---

## Usage

### Basic Usage

```python
from DWSIM_PFR import Stage3State, solve_dispersion_fd, print_stage3_results

# Configure reactor (heat flux auto-estimated for near-isothermal operation)
state = Stage3State(
    D=0.1,              # 100 mm tube diameter
    L=10.0,             # 10 m length
    T_in=1000.0,        # inlet temperature [K]
    P_in=2e5,           # inlet pressure [Pa]
    v_z0=1.0,           # inlet velocity [m/s]
    Nz=100,             # 100 axial nodes
    use_taylor_aris=False,  # use turbulent correlation
    # q_flux_constant=None uses auto-estimate from energy balance
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

### Using Stage-2 Heat Flux Profile

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

### Generating OpenFOAM Case

```python
from DWSIM_PFR import Stage3State, generate_openfoam_case, solve_stage2

# Generate case with Stage-2 reaction sources
state = Stage3State(D=0.1, L=10.0, T_in=1000.0)
result2 = solve_stage2(Stage2State(...))

generate_openfoam_case(state, "./openfoam_case", stage2_result=result2)
```

### Comparing with Plug Flow

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

---

## Solution Algorithm

### Marching + Correction Method

1. **Initialize with plug flow**: March forward solving species/energy/pressure without dispersion

2. **Compute dispersion coefficients**: Calculate D_ax(z) and alpha_ax(z) from local conditions

3. **Apply dispersion corrections**: Iteratively add second-derivative corrections with under-relaxation

4. **Converge**: Repeat until temperature and composition changes are below tolerance

### Numerical Details

- **Discretization**: Central differences for d²/dz², upwind for d/dz
- **Under-relaxation**: Default α = 0.3 for stability
- **Grid**: Uniform spacing, typically 100 nodes for 10 m tube
- **Convergence**: dY_max < 1e-5, dT_max < 1 K

---

## OpenFOAM Case Structure

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

### Running OpenFOAM

```bash
cd openfoam_case
blockMesh
reactingFoam
```

### Post-Processing

```bash
foamToVTK
# Open in ParaView for visualization
```

---

## Configuration Options

### Stage3State Parameters

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
| `q_flux_constant` | 50000 | Heat flux [W/m²] if no profile |
| `q_flux_profile` | None | Optional q''(z) array |

---

## Typical Results

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

---

## Limitations

1. **1-D approximation**: No radial gradients in dispersion model
2. **Steady-state only**: No transient effects
3. **Simplified correlations**: D_ax from empirical fits, not first-principles
4. **No coking effects**: No fouling or carbon deposition
5. **OpenFOAM setup**: Basic case generation, may need manual tuning

---

## Validation Workflow

1. **Generate OpenFOAM case** from Python Stage-3

2. **Run CFD simulation** with full 3D or axisymmetric mesh

3. **Extract effective D_ax(z)** from CFD results:
   - Compute flux: J = -D_eff * dC/dz + v*C
   - Solve for D_eff at each z location

4. **Import D_ax(z) into Python** using `import_openfoam_dispersion()`

5. **Re-solve Stage-3** with CFD-derived D_ax

6. **Compare** conversion and profiles between Python and CFD

---

## References

1. **Taylor dispersion**: G.I. Taylor, "Dispersion of soluble matter in solvent flowing slowly through a tube," Proc. R. Soc. London A, 219, 186-203 (1953)

2. **Aris extension**: R. Aris, "On the dispersion of a solute in a fluid flowing through a tube," Proc. R. Soc. London A, 235, 67-77 (1956)

3. **Turbulent dispersion**: O. Levenspiel, "Chemical Reaction Engineering," 3rd ed., Wiley (1999), Ch. 13

4. **Danckwerts BCs**: P.V. Danckwerts, "Continuous flow systems: Distribution of residence times," Chem. Eng. Sci., 2, 1-13 (1953)

