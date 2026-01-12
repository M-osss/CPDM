"""Kinetics parser for reduced_by_GS_mechanism_reactions.csv

Parses the CSV mechanism file into a list of reaction dictionaries
suitable for the Stage-1 PFR model.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Dict, List

__all__ = ["Reaction", "parse_csv_mechanism"]

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "reduced_by_GS_mechanism_reactions.csv"

class Reaction(dict):
    """Small helper so mypy / IDE knows keys."""

    @property
    def A(self) -> float:  # 1/s or (m^3/mol/s) depending order
        return self["A_si"]

    @property
    def n(self) -> float:
        return self["n"]

    @property
    def Ea(self) -> float:  # J/mol
        return self["Ea_si"]

    @property
    def stoich(self) -> Dict[str, float]:
        return self["stoich"]


N_A = 6.02214076e23  # Avogadro's number

def _arrhenius_SI(A_raw: str, Ea_raw: str, energy_unit: str, 
                  order: int, has_third_body: bool) -> tuple[float, float]:
    """Return (A, Ea) in SI units: A [m³/mol]^(order-1)/s, Ea [J/mol].
    
    CSV uses molecule-based cm³ units. Conversion:
    - For unimolecular (order=1): A [1/s] unchanged
    - For bimolecular (order=2): A [cm³/molecule/s] → A [m³/mol/s]
      multiply by N_A * 1e-6
    - For third-body M reactions: effective order increases by 1
    """
    # Pre-exponential may be given like 1.37*10^(-12)
    A = float(eval(A_raw.replace("^", "**")))
    Ea = float(eval(Ea_raw.replace("^", "**")))
    
    # Convert Ea to J/mol
    unit = energy_unit.strip().lower()
    if unit == "kj":
        Ea *= 1e3
    elif unit == "j":
        pass
    else:
        Ea *= 1e3  # assume kJ default
    
    # Convert A from cm³/molecule to m³/mol basis
    # effective order = order + 1 if third body present
    eff_order = order + (1 if has_third_body else 0)
    
    if eff_order == 1:
        # Unimolecular: A in 1/s, no conversion needed
        pass
    elif eff_order == 2:
        # Bimolecular: cm³/molecule/s → m³/mol/s
        # 1 cm³ = 1e-6 m³, multiply by N_A to go from per-molecule to per-mol
        A = A * N_A * 1e-6
    elif eff_order == 3:
        # Termolecular: (cm³)²/(molecule)²/s → (m³)²/(mol)²/s
        A = A * (N_A ** 2) * (1e-6 ** 2)
    else:
        raise ValueError(f"Unsupported reaction order: {eff_order}")
    
    return A, Ea


def _parse_stoich(label: str) -> tuple[Dict[str, float], bool]:
    """Parse stoichiometry from reaction label.
    
    Returns (stoich_dict, has_third_body).
    Third body 'M' is excluded from stoichiometry but flagged.
    """
    left, right = label.split("=")
    has_M = False
    def _side(part: str, sign: int):
        nonlocal has_M
        out: Dict[str, float] = {}
        for token in part.split("+"):
            sp = token.strip()
            if not sp:
                continue
            if sp == "M":
                has_M = True
                continue  # exclude M from stoich
            out[sp] = out.get(sp, 0.0) + sign
        return out
    stoich = _side(left, -1)
    for sp, coeff in _side(right, 1).items():
        stoich[sp] = stoich.get(sp, 0.0) + coeff
    return stoich, has_M


def parse_csv_mechanism(path: Path | None = None) -> List[Reaction]:
    file = path or CSV_PATH
    reactions: List[Reaction] = []
    with open(file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            label = row["reaction"].strip()
            A_raw = row["A"].strip()
            n = float(row[" n"].strip())
            Ea_raw = row[" Ea"].strip()
            energy_unit = row[" units for Energy"].strip()
            stoich, has_M = _parse_stoich(label)
            # Extract reactant species (negative stoich coefficients)
            reactants = {sp: -coef for sp, coef in stoich.items() if coef < 0}
            # Reaction order = sum of reactant stoichiometric coefficients
            order = int(sum(reactants.values()))
            
            A_si, Ea_si = _arrhenius_SI(A_raw, Ea_raw, energy_unit, order, has_M)
            reactions.append(
                Reaction(
                    label=label,
                    A_si=A_si,
                    n=n,
                    Ea_si=Ea_si,
                    stoich=stoich,
                    reactants=reactants,
                    has_third_body=has_M,
                    order=order,
                )
            )
    return reactions

if __name__ == "__main__":
    rxs = parse_csv_mechanism()
    print(f"Parsed {len(rxs)} reactions, example:\n", rxs[0])
