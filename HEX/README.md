# DWSIM Heat Exchanger Models

Heat exchanger unit models for process simulation. All 10 models pass the integrated test suite (`python test_all_hex_models.py`).

## Models

| Model | File | Description |
|-------|------|-------------|
| IDAES HeatExchanger | `idaes_heat_exchanger.py` | 0D LMTD counter-current; all operating conditions in a single INPUT PARAMETERS block at module top; Peng-Robinson EoS; water-toluene default |
| MixtureHeatExchanger | `mixture_heat_exchanger.py` | Thermo library, phase-aware; INPUT PARAMETERS at module top |
| IDAES LMTD | `idaes_hex_lmtd.py` | Database-driven; edit system presets in `hex_component_properties` |
| IDAES NTU | `idaes_hex_ntu.py` | NTU method; edit system presets in `hex_component_properties` |
| IDAES LC | `idaes_hex_lc.py` | Lumped capacitance; edit system presets + LC params |
| IDAES 1D | `idaes_hex_1d.py` | Spatial discretization; edit system presets in `hex_component_properties` |
| ShellAndTube1D | `idaes_hex_shell_tube_1d.py` | Shell-and-tube; edit system presets in `hex_component_properties` |
| PlateHeatExchanger | `idaes_hex_plate.py` | MEA/CO2 PHE; INPUT PARAMETERS at module top |
| BoilerHeatExchanger | `idaes_hex_boiler.py` | Steam/flue gas; INPUT PARAMETERS at module top |
| Boiler2D | `idaes_hex_boiler_2d.py` | 2D crossflow boiler; INPUT PARAMETERS at module top |

## Strengths

- **idaes_heat_exchanger**: All independent parameters (flowrates, inlet T/P, mole fractions, pressure drops, outlet T spec, U) consolidated at module top; no hunting through code to change operating conditions. Self-contained, no external databases.
- **mixture_heat_exchanger**: INPUT PARAMETERS block at top (components, T/P, flow, mode); phase change handling; thermo library.
- **idaes_hex_lmtd/ntu/lc/1d/shell_tube**: Edit system presets (HYDROCARBON_SYSTEM, etc.) for flow/T/P/composition; `INPUT_PARAMETERS` for geometry only. Currently for liquid only, for vle change a package in props_cfg, built - see the IDAES documentation. Or just use mixture hex model or boiler models. 
- **idaes_hex_1d/shell_tube_1d**: Spatial resolution; transport correlations from ChemSep.
- **idaes_hex_plate**: INPUT PARAMETERS at top; MEA/CO2 capture; built-in transport properties.
- **idaes_hex_boiler/boiler_2d**: INPUT PARAMETERS at top; steam/flue gas; iapws95 + FlueGas.

## Weaknesses

- **idaes_heat_exchanger**: Component property data (Perry's, NIST coefficients) still hardcoded inside `main()`; adding a new component requires editing the property config dict, not just the INPUT PARAMETERS block. Liquid-only phase; no VLE. No validation that mole fractions sum to 1.
- **mixture_heat_exchanger**: Uses `np.trapz` fallback for NumPy <1.22; limited for azeotropes and highly non-ideal mixtures. Infeasible cases (e.g. cold flow too low to absorb hot duty) raise `ValueError` with diagnosis.
- **idaes_hex_lmtd**: Hot outlet T must be < hot inlet; defaults to `hot_inlet - 30 K` if not in geometry.
- **idaes_hex_ntu**: EASY_SYSTEM (water_toluene_nitrogen) uses 373 K hot inlet; liquid-only package fails at vapor conditions (e.g. 573 K).
- **idaes_hex_lc**: Dynamic wall mode (`dynamic_heat_balance=True`) can fail to initialize; run in steady-state mode for robustness. Default system is liquid (power_plant not supported).
- **idaes_hex_1d**: Two-phase conditions with liquid-only property package produce warnings; use VLE package for phase change.
- **idaes_hex_plate**: Deprecated in IDAES; may fail on Pyomo 6.7+.
- **Boiler models**: Require `idaes.models_extra.power_generation`.

## Usage

**Run all tests:**
```bash
cd DWSIM_HEX
python test_all_hex_models.py
```

**Run a single model:**
```bash
python idaes_hex_lmtd.py liquid          # liquid system
python idaes_hex_lmtd.py water_toluene_nitrogen
python idaes_hex_ntu.py hydrocarbon      # hydrocarbon mix
```

**Where to edit parameters:**
- **Self-contained** (idaes_heat_exchanger, mixture_heat_exchanger, idaes_hex_plate, idaes_hex_boiler, idaes_hex_boiler_2d): Edit the INPUT PARAMETERS block at the top of each file.
- **Database-driven** (idaes_hex_lmtd, ntu, lc, 1d, shell_tube_1d): Edit the system preset in `hex_component_properties.py` for flow, T, P, composition. Valid systems: liquid, water_toluene, water_toluene_nitrogen, hydrocarbon, co2_capture. power_plant is not supported; use idaes_hex_boiler or idaes_hex_boiler_2d.
