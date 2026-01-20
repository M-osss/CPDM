"""Property helper utility for Stage-1 PFR.

Provides:
    • load_species_db() → dict of species data
    • sensible_enthalpy(T, name)  [J/mol] using NASA7
    • cp_molar(T, name)           [J/mol/K]
    • gas_viscosity(T, name)      [Pa·s] (Lucas fit if present else None)
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, Any

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "species_properties.json"

R = 8.314462618  # J/mol/K

# -----------------------------------------------------------------------------
# JSON loading
# -----------------------------------------------------------------------------

def _nasa_cp(coeffs: list[float], T: float) -> float:
    a1, a2, a3, a4, a5, *_ = coeffs
    return R * (a1 + a2 * T + a3 * T ** 2 + a4 * T ** 3 + a5 * T ** 4)

def _nasa_h(coeffs: list[float], T: float) -> float:
    a1, a2, a3, a4, a5, a6, *_ = coeffs
    return R * T * (
        a1 + a2 * T / 2 + a3 * T ** 2 / 3 + a4 * T ** 3 / 4 + a5 * T ** 4 / 5 + a6 / T
    )

def load_species_db(path: Path | None = None) -> Dict[str, Any]:
    data = json.loads((path or DB_PATH).read_text())
    return data["species"]

_DB = None


def _get_db() -> Dict[str, Any]:
    global _DB
    if _DB is None:
        _DB = load_species_db(DB_PATH)
    return _DB

_FORMULA_MAP: Dict[str, str] | None = None

def _get_formula_map() -> Dict[str, str]:
    """Return mapping from formula token (e.g. 'C2H6') to full db key."""
    global _FORMULA_MAP
    if _FORMULA_MAP is None:
        data = json.loads(DB_PATH.read_text())
        f2e = data.get("formula_to_entries", {})
        _FORMULA_MAP = {k: v[0] for k, v in f2e.items() if v}
    return _FORMULA_MAP


def resolve_species(token: str) -> str:
    """Convert mechanism formula token to species_properties.json key.
    
    E.g. 'C2H6' -> 'C2H6::ethane::74-84-0'
    Returns token unchanged if already a full key or if 'M' (third body).
    """
    if "::" in token:
        return token
    if token == "M":
        return "M"  # third-body placeholder
    fm = _get_formula_map()
    if token in fm:
        return fm[token]
    raise KeyError(f"Species '{token}' not found in formula_to_entries")

# -----------------------------------------------------------------------------
# Thermo helpers
# -----------------------------------------------------------------------------

def cp_molar(T: float, species: str) -> float:
    sp = _get_db()[species]
    nasa = sp["thermo"]["nasa7"]
    if nasa is None:
        raise ValueError(f"No NASA data for {species}")
    T_mid = nasa["T_mid_K"]
    coeffs = nasa["nasa7_high"] if T >= T_mid else nasa["nasa7_low"]
    return _nasa_cp(coeffs, T)


def sensible_enthalpy(T: float, species: str) -> float:
    sp = _get_db()[species]
    nasa = sp["thermo"]["nasa7"]
    if nasa is None:
        raise ValueError(f"No NASA data for {species}")
    T_mid = nasa["T_mid_K"]
    coeffs = nasa["nasa7_high"] if T >= T_mid else nasa["nasa7_low"]
    return _nasa_h(coeffs, T)

# -----------------------------------------------------------------------------
# Viscosity – use fitted ln poly if available
# -----------------------------------------------------------------------------

def gas_viscosity(T: float, species: str) -> float:
    sp = _get_db()[species]
    gas = sp.get("transport", {}).get("gas")
    if not gas:
        raise ValueError("Missing transport block")
    fit = gas["fits"].get("mu_fit")
    if not fit:
        raise ValueError("No viscosity fit for species")
    c = fit["coeffs_high_to_low"]
    # ln(mu) = c0 ln^3 T + c1 ln^2 T + c2 ln T + c3  (if degree 3)
    lnT = math.log(T)
    deg = len(c) - 1
    ln_mu = 0.0
    for i, coeff in enumerate(c):
        ln_mu += coeff * lnT ** (deg - i)
    return math.exp(ln_mu)

if __name__ == "__main__":
    db = load_species_db()
    print("Loaded", len(db), "species")
    print("Cp_ethane(1000K)", cp_molar(1000, "C2H6::ethane::74-84-0"))
