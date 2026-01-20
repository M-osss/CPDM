# Stage-2 PFR: Heat Transfer Model - Technical Documentation

## Overview

Stage-2 extends the ideal PFR (Stage-1) with **coupled 2-D heat transfer** through the tube wall and **gas-phase radiation**. This enables realistic simulation of industrial ethane cracking furnaces where heat is supplied from an external firebox.

---

## Key Differences from Stage-1

### Stage-1 (Ideal PFR)
- **Constant heat flux** q''(z) = constant
- **No wall model**: Heat flux is a user input or estimated from simple energy balance
- **Isothermal assumption**: Designed to maintain near-constant temperature
- **Simplified**: Suitable for preliminary design and validation

### Stage-2 (Heat Transfer Model)
- **Variable heat flux** q''(z) computed from wall-gas temperature difference
- **Wall conduction**: Radial temperature gradient through tube wall
- **Coupled boundaries**: Inner (gas-wall) and outer (furnace-wall) heat transfer
- **Radiative transfer**: Gas and wall radiation included
- **Realistic**: Matches industrial furnace operation

---

## Physical Model

### 1. Tube Wall Conduction

The tube wall is modeled as a **cylindrical shell** with radial conduction:

```
ρ_w * Cp_w * ∂T_w/∂t = (1/r) * ∂/∂r(k_w * r * ∂T_w/∂r)
```

For steady-state operation (assumed here):
```
(1/r) * d/dr(k_w * r * dT_w/dr) = 0
```

**Boundary Conditions:**
- **Inner surface** (r = r_inner): Heat flux from wall to gas
- **Outer surface** (r = r_outer): Heat flux from furnace to wall

**Solution:**
For constant thermal conductivity k_w:
```
T_w(r) = T_outer + (q''_inner * r_inner / k_w) * ln(r_outer / r)
```

### 2. Inner Boundary: Gas-Wall Interface

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
- Gas radiation is weak for dry pyrolysis (no H₂O/CO₂), so dominated by wall emission

### 3. Outer Boundary: Furnace-Wall Interface

Heat flux from furnace to outer wall:

```
q''_outer = h_o * (T_furnace - T_wall) + ε_w * σ * (T_furnace^4 - T_wall^4)
```

**Typical values:**
- h_o = 60 W/(m²·K) (external convection in firebox)
- T_furnace = 1300 K (radiant section temperature)
- Radiative term dominates at high temperatures (T⁴ scaling)

### 4. Gas-Phase Radiation (WSGG Model)

**Weighted Sum of Gray Gases (WSGG)** model for participating species:

```
κ_eff = Σ_i a_i(T) * κ_i * P_participating
```

Where:
- κ_i: absorption coefficients for gray gases
- a_i(T): temperature-dependent weights
- P_participating: partial pressure of H₂O + CO₂

**For ethane pyrolysis:**
- No steam dilution → P_participating ≈ 0
- Gas is optically thin (minimal self-absorption)
- Model returns minimal κ ≈ 0.01 m⁻¹

**Gas emissivity:**
```
ε_g = Σ_i a_i(T) * [1 - exp(-κ_i * P * L_beam)]
```

Where L_beam ≈ 0.9×D (mean beam length for cylinder).

### 5. Energy Balance Coupling

The gas-phase energy equation uses the computed heat flux:

```
dT/dz = (1/(ρ * Cp_mass * v_z)) * [Q_rxn + 4*q''(z)/D]
```

Where q''(z) is **spatially varying** and computed from:
1. Wall temperature profile (from conduction)
2. Gas temperature (from previous iteration)
3. Convective and radiative heat transfer

---

## Solution Algorithm: Operator-Split Iteration

The coupled system is solved using **operator-splitting** with under-relaxation:

### Step 1: Initialize
- Set initial gas temperature profile: T_gas(z) = T_in
- Set initial wall temperatures: T_wall ≈ T_in + 50 K
- Set initial heat flux estimate: q''(z) = 20 kW/m²

### Step 2: Compute Heat Flux Profile
For each axial position z:
1. Get current gas temperature T_gas(z) and composition Y(z)
2. Get current inner wall temperature T_wall_inner(z)
3. Compute gas properties: μ, ρ, Cp, k, Pr, Re
4. Compute Nusselt number (Gnielinski)
5. Compute convective heat flux: q_conv = h_i × (T_wall - T_gas)
6. Compute radiative heat flux: q_rad = ε_w × σ × (T_wall⁴ - T_gas⁴)
7. Total: q''(z) = q_conv + q_rad

### Step 3: Integrate Species/T/P ODEs
- Use Stage-1 ODE system with **variable** q''(z) from Step 2
- Integrate from z=0 to z=L using `solve_ivp` (BDF method)
- Extract new profiles: T_gas(z), Y(z), P(z)

### Step 4: Update Wall Temperatures
For each axial position z:
1. Compute outer heat flux from furnace:
   ```
   q_out = h_o*(T_furnace - T_wo) + ε_w*σ*(T_furnace⁴ - T_wo⁴)
   ```
2. Compute inner heat flux to gas (from Step 2)
3. Solve wall conduction to match: q_out = q_wall = q_in
4. Update T_wall_outer and T_wall_inner

### Step 5: Check Convergence
- Compute maximum temperature change: ΔT_max = max|T_gas_new - T_gas_old|
- If ΔT_max < tolerance (1 K): **converged**
- Otherwise: apply under-relaxation and repeat from Step 2

**Under-relaxation:**
```
T_gas = α * T_gas_new + (1 - α) * T_gas_old
q'' = α * q''_new + (1 - α) * q''_old
```
Where α = 0.3 (default) prevents oscillation.

---

## Material Properties

### Wall Material: Inconel 800H

Typical for ethylene cracking furnaces:
- **Density**: ρ = 7940 kg/m³
- **Thermal conductivity**: k = 25 W/(m·K) at ~1000 K
- **Specific heat**: Cp = 500 J/(kg·K)
- **Emissivity**: ε = 0.85 (oxidized surface)
- **Thermal diffusivity**: α = k/(ρ·Cp) ≈ 6.3×10⁻⁶ m²/s

### Gas Properties

**Thermal conductivity** (Eucken correlation):
```
k_gas ≈ 0.02 + 7×10⁻⁵ * T  [W/(m·K)]
```

**Viscosity**: From fitted ln-polynomial in `species_properties.json`

**Heat capacity**: NASA-7 polynomials from RMG database

---

## Typical Results

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

---

## Numerical Considerations

### Convergence
- **Typical iterations**: 15-25 for ΔT < 1 K
- **Under-relaxation**: α = 0.3 prevents oscillation
- **Tolerance**: 1 K on gas temperature

### Stability
- **Heat flux clamping**: Limited to 1-200 kW/m² to prevent unphysical values
- **Temperature clamping**: T_gas ∈ [300, 2500] K
- **Wall-gas constraint**: T_wall > T_gas (heating only)

### Discretization
- **Axial nodes**: Nz = 50-100 (default 50)
- **Radial wall nodes**: Nr = 5 (sufficient for thin wall)
- **ODE solver**: BDF method with adaptive step size

---

## Limitations and Assumptions

1. **Steady-state wall**: No transient response (quasi-steady assumption)
2. **1-D gas flow**: Plug flow, no radial gradients in gas phase
3. **Simplified WSGG**: Generic H₂O/CO₂ coefficients, not tuned for hydrocarbons
4. **No coking**: No carbon deposition on inner wall
5. **No steam dilution**: Pure ethane feed (industrial uses steam)
6. **Constant wall properties**: k, Cp, ρ assumed constant (weak temperature dependence)
7. **Thin wall approximation**: Wall conduction uses simplified 1-D model

---

## Usage Example

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

# Solve with custom relaxation
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

# Access solution arrays
z = result["z"]                    # axial positions [m]
T_gas = result["T_gas"]            # gas temperature [K]
T_wall_inner = result["T_wall_inner"]  # inner wall T [K]
T_wall_outer = result["T_wall_outer"]  # outer wall T [K]
q_inner = result["q_inner"]        # heat flux [W/m²]
Y = result["Y"]                    # mole fractions (n_species, Nz)
```

---

## Comparison with Industrial Data

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

## References

1. **Gnielinski correlation**: V. Gnielinski, "New equations for heat and mass transfer in turbulent pipe and channel flow," Int. Chem. Eng., 16, 359-368 (1976)

2. **WSGG model**: T.F. Smith et al., "Evaluation of coefficients for the weighted sum of gray gases model," J. Heat Transfer, 104, 602-608 (1982)

3. **Wall conduction**: Standard cylindrical coordinate heat conduction (see any heat transfer textbook)

4. **Industrial furnaces**: S. Matar and L.F. Hatch, "Chemistry of Petrochemical Processes," 2nd ed., Gulf Publishing (2001)

