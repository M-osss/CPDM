# Setup Guide: DWSIM_PFR

This guide explains how to set up and run the DWSIM_PFR package from scratch.

## Required Packages

The package requires the following Python packages:

### External Dependencies
- **numpy** - Numerical computing
- **scipy** - Scientific computing (specifically `scipy.integrate.solve_ivp`)

### Standard Library (included with Python)
- `json` - JSON file handling
- `math` - Mathematical functions
- `pathlib` - Path operations
- `typing` - Type hints
- `csv` - CSV file parsing

## Installation

### Option 1: Using pip

```bash
pip install numpy scipy
```

### Option 2: Using conda

```bash
conda install numpy scipy
```

### Option 3: Using requirements.txt

Create a `requirements.txt` file in the DWSIM_PFR folder:

```
numpy>=1.20.0
scipy>=1.7.0
```

Then install:

```bash
pip install -r requirements.txt
```

## Required Data Files

Ensure the following data files are present in the `DWSIM_PFR` folder:

- `reduced_by_GS_mechanism_reactions.csv` - Reaction mechanism (18 reactions)
- `species_properties.json` - Species thermophysical properties database
- `ethane_mechanism_complete.json` - Reference RMG mechanism (optional, for metadata)
- `ethane_pyrolysis_rmg.json` - Reference RMG kinetics (optional, for reference)

## Running the Program

### Method 1: Run as Standalone Script

Navigate to the `DWSIM_PFR` folder and run:

```bash
python pfr_ideal.py
```

This will execute the default simulation with:
- Tube diameter: 0.1 m
- Tube length: 10.0 m
- Inlet temperature: 1000 K
- Inlet pressure: 2.0 bar
- Superficial velocity: 1.0 m/s

### Method 2: Import as Python Module

If the `DWSIM_PFR` folder is in your Python path, you can import and use it:

```python
from DWSIM_PFR import run_pfr, print_results

# Run with default parameters
sol = run_pfr()
print_results(sol)

# Run with custom parameters
sol = run_pfr(
    L=10.0,           # reactor length [m]
    T_in=1000.0,      # inlet temperature [K]
    P_in=1.5e5,       # inlet pressure [Pa]
    v_z0=2.0,         # superficial velocity [m/s]
    D=0.1,            # tube diameter [m]
    q_flux=5000.0,    # wall heat flux [W/m^2]
)
print_results(sol)
```

### Method 3: Install as Package (Optional)

To use `DWSIM_PFR` from any location, install it in development mode:

```bash
cd DWSIM_PFR
pip install -e .
```

(Note: This requires a `setup.py` or `pyproject.toml` file. If not present, use Method 2 with proper Python path configuration.)

## Verification

After installation, verify the setup by running:

```bash
python -c "import numpy; import scipy; print('All packages installed successfully')"
```

Then test the PFR model:

```bash
cd DWSIM_PFR
python pfr_ideal.py
```

You should see output showing:
- Integration success status
- Temperature and pressure profiles
- Species mole fractions at inlet and outlet
- Ethane conversion percentage

## Troubleshooting

### ImportError: No module named 'numpy' or 'scipy'
- Install missing packages: `pip install numpy scipy`

### FileNotFoundError for CSV or JSON files
- Ensure all data files are in the same folder as the Python scripts
- Check that file names match exactly (case-sensitive on Linux/Mac)

### Relative import errors when running as script
- The code handles both module and standalone execution automatically
- If issues persist, ensure you're running from the `DWSIM_PFR` directory

## Python Version

The code uses Python 3.7+ features:
- `from __future__ import annotations` (PEP 563)
- Type hints with `typing` module
- `pathlib.Path` for file operations

Recommended: Python 3.8 or higher.

