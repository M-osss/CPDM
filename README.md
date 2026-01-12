# PFR: Stage-1 Ideal Plug-Flow Reactor for Ethane Pyrolysis (Stage 2 to be improved upon that, see Stage 2)

## Model Overview

1-D ideal PFR solving coupled species, energy, and pressure ODEs for thermal cracking of ethane to ethylene. Implements the Stage-1 requirements from the implementation plan.

## Governing Equations

**Species balance** (molar fractions):
```
dY_i/dz = (1/v_z) * (1/C_total) * sum_j(nu_ij * Omega_j)
```

**Energy balance**:
```
dT/dz = (1/(rho * Cp_mass * v_z)) * [Q_rxn + 4*q''/D]
```
where Q_rxn = -sum_j(dH_j * Omega_j) is volumetric reaction heat.

**Pressure drop** (Darcy-Weisbach with Blasius friction):
```
dP/dz = -4*f/D * rho*v_z^2/2
f = 0.0791 * Re^(-0.25)
```

**Velocity** from continuity:
```
v_z = v_z0 * (P0/P) * (T/T0)
```

## Kinetics

- 18 reactions from reduced Glarborg-Sendt mechanism
- Pre-exponential factors converted from molecule/cm^3 to mol/m^3 basis
- Third-body reactions (M) handled with total concentration [M] = P/(R*T)
- Arrhenius form: k = A * (T/298.15)^n * exp(-Ea/(R*T))

## Thermophysical Properties

- NASA-7 polynomials for Cp(T) and H(T) from RMG database
- Gas viscosity from fitted ln-polynomial correlations
- Ideal gas equation of state

## Usage

```python
from DWSIM_PFR import run_pfr, print_results

# Default: 1200 K, 2 bar, 5 m tube, pure ethane feed
sol = run_pfr()
print_results(sol)

# Custom conditions
sol = run_pfr(
    L=10.0,           # reactor length [m]
    T_in=1150.0,      # inlet temperature [K]
    P_in=1.5e5,       # inlet pressure [Pa]
    v_z0=2.0,         # superficial velocity [m/s]
    D=0.04,           # tube diameter [m]
    q_flux=5000.0,    # wall heat flux [W/m^2]
)
```

## Files

| File | Description |
|------|-------------|
| `pfr_ideal.py` | Main Stage-1 PFR model |
| `kinetics_parser.py` | CSV mechanism parser with unit conversion |
| `props.py` | Thermophysical property functions |
| `species_properties.json` | NASA-7 + transport data |
| `reduced_by_GS_mechanism_reactions.csv` | 18-reaction mechanism |

---

## Critical Analysis

### Strengths

1. **Proper unit conversion**: Pre-exponential factors correctly converted from molecule/cm^3 basis to mol/m^3 for all reaction orders (uni-, bi-, termolecular).

2. **Third-body handling**: Reactions with M are identified and multiplied by total concentration.

3. **Coupled ODEs**: Species, energy, and momentum equations solved simultaneously with BDF stiff solver.

4. **Heat flux estimation**: Automatic estimate based on target conversion and reaction enthalpy.

5. **Property database integration**: Uses resolved species keys to look up NASA-7 coefficients and transport properties.

### Weaknesses

1. **No reverse reactions**: All reactions treated as irreversible. At high T and partial equilibrium, reverse rates can be significant (especially for H2 + C2H4 <-> C2H6).

2. **Linear mixing rules**: Viscosity and Cp use simple molar-fraction weighting instead of Wilke or Wassiljewa rules.

3. **No radial gradients**: Ideal plug flow assumes flat velocity and concentration profiles. Industrial cracking furnaces have significant radial T gradients.

4. **Constant q''**: Heat flux is uniform along reactor. Real furnaces have axial q'' profiles based on firebox geometry.

5. **Pressure drop underestimated**: Blasius correlation valid only for smooth tubes. Real coils have roughness and bends.

6. **No coking model**: No carbon deposition kinetics or tube fouling effects.

7. **Temperature exponent form**: Uses (T/298.15)^n instead of standard T^n; may introduce small errors for high n values.

8. **Mole-fraction normalization**: Renormalizes Y after each step; can mask mass conservation errors.

### Recommendations for Stage-2

1. Add equilibrium constant K_eq(T) from NASA polynomials to enable reversible reactions.
2. Implement Wilke viscosity mixing rule.
3. Add radial energy balance with wall heat transfer coefficient.
4. Allow q''(z) profile input.
5. Consider adding simplified coking rate for long-term operation.


