from __future__ import annotations

import csv
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Optional, Tuple


CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")


def _workspace_root() -> Path:
    # DWSIM_HEX/property_database_interface.py -> workspace root
    return Path(__file__).resolve().parents[1]


def _norm(s: str) -> str:
    s = (s or "").strip().lower()
    s = re.sub(r"[\s_]+", " ", s)
    s = re.sub(r"[^a-z0-9 \-()+,.]", "", s)
    return s


def _as_float(val: Optional[str]) -> Optional[float]:
    if val is None:
        return None
    try:
        return float(val)
    except Exception:
        return None


@dataclass(frozen=True)
class ChemSepCorrelation:
    """Raw ChemSep correlation block as stored in chemsep XML."""

    eqno: Optional[int]
    coeffs: Dict[str, float]  # keys typically A..E (and sometimes others)
    t_min_k: Optional[float] = None
    t_max_k: Optional[float] = None
    units: Optional[str] = None


@dataclass(frozen=True)
class ChemSepCompound:
    name: str
    cas: Optional[str]
    mw_kg_per_kmol: Optional[float]
    tc_k: Optional[float]
    pc_pa: Optional[float]
    omega: Optional[float]

    # Correlations (subset we care about for IDAES HX work)
    liquid_density: Optional[ChemSepCorrelation] = None
    liquid_cp: Optional[ChemSepCorrelation] = None
    ig_cp: Optional[ChemSepCorrelation] = None
    rpp_ig_cp: Optional[ChemSepCorrelation] = None
    antoine_vp: Optional[ChemSepCorrelation] = None
    liq_visc: Optional[ChemSepCorrelation] = None
    vap_visc: Optional[ChemSepCorrelation] = None
    liq_k: Optional[ChemSepCorrelation] = None
    vap_k: Optional[ChemSepCorrelation] = None

    # Formation/reference data (ChemSep uses IG refs)
    ig_hf_j_per_kmol: Optional[float] = None
    ig_s_abs_j_per_kmol_k: Optional[float] = None


class ChemSepXMLDatabase:
    """
    Read-only, cached view of ChemSep `pcd/chemsep1.xml`.

    This is the authoritative local source for:
    - Critical properties (Tc, Pc), MW, omega
    - Liquid density correlation (LiquidDensity)
    - Liquid Cp (LiquidHeatCapacityCp)
    - Ideal-gas Cp (RPPHeatCapacityCp or IdealGasHeatCapacityCp)
    - Transport correlations (viscosity, thermal conductivity)

    NOTE:
    - This database provides **correlation coefficients**, not IDAES-ready formatting.
    - Formatting to IDAES config_dict is done in `idaes_property_builder.py`.
    """

    def __init__(self, xml_path: Optional[Path] = None):
        root = _workspace_root()
        self.xml_path = xml_path or (root / "pcd" / "chemsep1.xml")
        if not self.xml_path.exists():
            raise FileNotFoundError(f"ChemSep XML not found: {self.xml_path}")

        self._by_cas: Dict[str, ChemSepCompound] = {}
        self._by_name: Dict[str, str] = {}  # normalized name -> CAS (when available)
        self._loaded = False

    def ensure_loaded(self) -> None:
        if self._loaded:
            return
        self._load()
        self._loaded = True

    def _parse_corr(self, elem: ET.Element) -> ChemSepCorrelation:
        # correlation element has: <eqno>, <A>, <B>, ... and Tmin/Tmax with units="K"
        eqno = None
        coeffs: Dict[str, float] = {}
        tmin = None
        tmax = None
        units = elem.attrib.get("units")

        for child in list(elem):
            tag = child.tag
            if tag == "eqno":
                eqno = int(child.attrib.get("value"))
                continue
            if tag in {"A", "B", "C", "D", "E", "F", "G", "H", "a0", "a1", "a2", "a3", "a4"}:
                v = _as_float(child.attrib.get("value"))
                if v is not None:
                    coeffs[tag] = v
                continue
            if tag == "Tmin":
                tmin = _as_float(child.attrib.get("value"))
                continue
            if tag == "Tmax":
                tmax = _as_float(child.attrib.get("value"))
                continue

        return ChemSepCorrelation(eqno=eqno, coeffs=coeffs, t_min_k=tmin, t_max_k=tmax, units=units)

    def _load(self) -> None:
        # Stream parse to avoid loading entire XML into memory.
        # One compound element is processed and cleared at a time.
        ctx = ET.iterparse(self.xml_path, events=("end",))
        for event, elem in ctx:
            if elem.tag != "compound":
                continue

            name = None
            cas = None
            mw = None
            tc = None
            pc = None
            omega = None
            ig_hf = None
            ig_sabs = None

            liquid_density = None
            liquid_cp = None
            ig_cp = None
            rpp_ig_cp = None
            antoine_vp = None
            liq_visc = None
            vap_visc = None
            liq_k = None
            vap_k = None

            for child in list(elem):
                tag = child.tag
                if tag == "CompoundID":
                    name = child.attrib.get("value")
                elif tag == "CAS":
                    cas = child.attrib.get("value")
                elif tag == "MolecularWeight":
                    mw = _as_float(child.attrib.get("value"))
                elif tag == "CriticalTemperature":
                    tc = _as_float(child.attrib.get("value"))
                elif tag == "CriticalPressure":
                    pc = _as_float(child.attrib.get("value"))
                elif tag == "AcentricityFactor":
                    omega = _as_float(child.attrib.get("value"))
                elif tag == "HeatOfFormation":
                    ig_hf = _as_float(child.attrib.get("value"))
                elif tag == "AbsEntropy":
                    ig_sabs = _as_float(child.attrib.get("value"))
                elif tag == "LiquidDensity":
                    liquid_density = self._parse_corr(child)
                elif tag == "LiquidHeatCapacityCp":
                    liquid_cp = self._parse_corr(child)
                elif tag == "IdealGasHeatCapacityCp":
                    ig_cp = self._parse_corr(child)
                elif tag == "RPPHeatCapacityCp":
                    rpp_ig_cp = self._parse_corr(child)
                elif tag == "AntoineVaporPressure":
                    antoine_vp = self._parse_corr(child)
                elif tag == "LiquidViscosity":
                    liq_visc = self._parse_corr(child)
                elif tag == "VaporViscosity":
                    vap_visc = self._parse_corr(child)
                elif tag == "LiquidThermalConductivity":
                    liq_k = self._parse_corr(child)
                elif tag == "VaporThermalConductivity":
                    vap_k = self._parse_corr(child)

            comp = ChemSepCompound(
                name=name or "",
                cas=cas,
                mw_kg_per_kmol=mw,
                tc_k=tc,
                pc_pa=pc,
                omega=omega,
                liquid_density=liquid_density,
                liquid_cp=liquid_cp,
                ig_cp=ig_cp,
                rpp_ig_cp=rpp_ig_cp,
                antoine_vp=antoine_vp,
                liq_visc=liq_visc,
                vap_visc=vap_visc,
                liq_k=liq_k,
                vap_k=vap_k,
                ig_hf_j_per_kmol=ig_hf,
                ig_s_abs_j_per_kmol_k=ig_sabs,
            )

            if comp.cas and CAS_RE.match(comp.cas):
                self._by_cas[comp.cas] = comp
            if comp.name:
                n = _norm(comp.name)
                if comp.cas and CAS_RE.match(comp.cas):
                    self._by_name[n] = comp.cas
                else:
                    # still index name for find-by-name (CAS may be missing for pseudo-components)
                    self._by_name[n] = ""

            # free memory
            elem.clear()

    def get(self, cas_or_name: str) -> Optional[ChemSepCompound]:
        self.ensure_loaded()
        key = (cas_or_name or "").strip()
        if not key:
            return None

        if CAS_RE.match(key):
            return self._by_cas.get(key)

        n = _norm(key)
        cas = self._by_name.get(n)
        if cas:
            return self._by_cas.get(cas)

        # fallback: scan for exact normalized match with slight variants
        # (kept cheap; no fuzzy matching to avoid false positives)
        return None


class AntoineCSVDatabase:
    """Optional fallback vapor-pressure coefficients from `Antoine_Coefficients.csv` (name-based)."""

    def __init__(self, csv_path: Optional[Path] = None):
        root = _workspace_root()
        self.csv_path = csv_path or (root / "Antoine_Coefficients.csv")
        self._by_name: Dict[str, Dict[str, Any]] = {}
        self._loaded = False

    def ensure_loaded(self) -> None:
        if self._loaded:
            return
        if not self.csv_path.exists():
            self._loaded = True
            return
        with open(self.csv_path, "r", encoding="utf-8", errors="ignore", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = row.get("Compound Name") or row.get("Compound") or row.get("CompoundName")
                if not name:
                    continue
                A = _as_float(row.get("A"))
                B = _as_float(row.get("B"))
                C = _as_float(row.get("C"))
                tmin = _as_float(row.get("TMIN"))
                tmax = _as_float(row.get("TMAX"))
                if A is None or B is None or C is None:
                    continue
                self._by_name[_norm(name)] = {"A": A, "B": B, "C": C, "TMIN": tmin, "TMAX": tmax}
        self._loaded = True

    def get(self, name: str) -> Optional[Dict[str, Any]]:
        self.ensure_loaded()
        return self._by_name.get(_norm(name))


class IPDDatabase:
    """
    Minimal reader for ChemSep IPD files (e.g., PR, SRK interaction parameters).

    For HX property packages we only need a single k12 per pair.
    Many files contain multiple entries for same pair at different conditions;
    we select the **first** value encountered (conservative, deterministic).
    """

    def __init__(self, ipd_dir: Optional[Path] = None):
        root = _workspace_root()
        self.ipd_dir = ipd_dir or (root / "ipd")

    @staticmethod
    def _pair_key(cas1: str, cas2: str) -> Tuple[str, str]:
        return (cas1, cas2) if cas1 <= cas2 else (cas2, cas1)

    @lru_cache(maxsize=8)
    def _load_pr(self) -> Dict[Tuple[str, str], float]:
        fpath = self.ipd_dir / "pr.ipd"
        if not fpath.exists():
            return {}
        data: Dict[Tuple[str, str], float] = {}
        in_data = False
        with open(fpath, "r", encoding="latin-1", errors="ignore") as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                if s.startswith("[IPD]"):
                    in_data = True
                    continue
                if not in_data:
                    continue
                if s.startswith("#") or s.startswith("Comment="):
                    continue
                parts = s.split()
                if len(parts) < 3:
                    continue
                cas1, cas2 = parts[0], parts[1]
                k12 = _as_float(parts[2])
                if not (CAS_RE.match(cas1) and CAS_RE.match(cas2)) or k12 is None:
                    continue
                key = self._pair_key(cas1, cas2)
                if key not in data:
                    data[key] = k12
        return data

    def get_pr_kij(self, cas1: str, cas2: str) -> Optional[float]:
        data = self._load_pr()
        return data.get(self._pair_key(cas1, cas2))


@lru_cache(maxsize=1)
def _chemsep_db() -> ChemSepXMLDatabase:
    return ChemSepXMLDatabase()


@lru_cache(maxsize=1)
def _antoine_csv_db() -> AntoineCSVDatabase:
    return AntoineCSVDatabase()


@lru_cache(maxsize=1)
def _ipd_db() -> IPDDatabase:
    return IPDDatabase()


def resolve_chemsep_compound(cas_or_name: str) -> Optional[ChemSepCompound]:
    return _chemsep_db().get(cas_or_name)


def get_critical_properties(cas_or_name: str) -> Dict[str, Any]:
    """
    Return Tc, Pc, omega, MW with provenance.
    Prefers ChemSep XML; falls back to `thermo` library if installed.
    """
    comp = resolve_chemsep_compound(cas_or_name)
    if comp and comp.tc_k and comp.pc_pa and comp.mw_kg_per_kmol is not None:
        return {
            "source": "chemsep_xml",
            "name": comp.name,
            "cas": comp.cas,
            "Tc_K": comp.tc_k,
            "Pc_Pa": comp.pc_pa,
            "omega": comp.omega,
            "MW_kg_per_kmol": comp.mw_kg_per_kmol,
        }

    # Fallback: thermo library
    try:
        from thermo import Chemical  # type: ignore

        chem = Chemical(cas_or_name)
        return {
            "source": "thermo",
            "name": getattr(chem, "name", cas_or_name),
            "cas": getattr(chem, "CAS", None),
            "Tc_K": getattr(chem, "Tc", None),
            "Pc_Pa": getattr(chem, "Pc", None),
            "omega": getattr(chem, "omega", None),
            "MW_kg_per_kmol": (getattr(chem, "MW", None) * 1000.0) if getattr(chem, "MW", None) else None,
        }
    except Exception:
        return {"source": "none", "name": cas_or_name, "cas": None, "Tc_K": None, "Pc_Pa": None, "omega": None, "MW_kg_per_kmol": None}


def get_perrys_density_coefficients(cas_or_name: str) -> Optional[Dict[str, Any]]:
    """
    Return liquid molar density coefficients in a form directly mappable to
    `Perrys.dens_mol_liq_comp` with `eqn_type == 1`.

    ChemSep XML uses LiquidDensity eqno=105 with A,B,C,D:
      rho = A / B ** ( 1 + (1 - T/C) ** D )
    which matches IDAES Perrys dens eqn_type=1.
    """
    comp = resolve_chemsep_compound(cas_or_name)
    if not comp or not comp.liquid_density:
        return None
    corr = comp.liquid_density
    if corr.eqno not in (105, 1):  # 105 observed in chemsep XML; 1 reserved as alias
        return None
    A = corr.coeffs.get("A")
    B = corr.coeffs.get("B")
    C = corr.coeffs.get("C")
    D = corr.coeffs.get("D")
    if A is None or B is None or C is None or D is None:
        return None
    return {
        "source": "chemsep_xml",
        "eqn_type": 1,
        "1": A,
        "2": B,
        "3": C,
        "4": D,
        "Tmin_K": corr.t_min_k,
        "Tmax_K": corr.t_max_k,
        "units": corr.units,
    }


def get_perrys_cp_liq_coefficients(cas_or_name: str) -> Optional[Dict[str, Any]]:
    """
    Return 5-term polynomial coefficients for `Perrys.cp_mol_liq_comp`.

    IDAES expects:
      Cp = c1 + c2*T + c3*T^2 + c4*T^3 + c5*T^4
    with coefficients stored as indices "1".."5".

    ChemSep XML LiquidHeatCapacityCp provides A..E for a 5-term polynomial.
    We map: c1=A, c2=B, c3=C, c4=D, c5=E.
    """
    comp = resolve_chemsep_compound(cas_or_name)
    if not comp or not comp.liquid_cp:
        return None
    corr = comp.liquid_cp
    A = corr.coeffs.get("A")
    B = corr.coeffs.get("B")
    C = corr.coeffs.get("C")
    D = corr.coeffs.get("D")
    E = corr.coeffs.get("E")
    if A is None or B is None or C is None or D is None or E is None:
        return None
    return {
        "source": "chemsep_xml",
        "1": A,
        "2": B,
        "3": C,
        "4": D,
        "5": E,
        "Tmin_K": corr.t_min_k,
        "Tmax_K": corr.t_max_k,
        "units": corr.units,
        "eqno": corr.eqno,
    }


def get_rpp_ig_cp_coefficients(cas_or_name: str) -> Optional[Dict[str, Any]]:
    """
    Return A..D coefficients for IDAES `RPP4.cp_mol_ig_comp` (ideal-gas Cp polynomial).

    ChemSep XML RPPHeatCapacityCp provides A..D.
    """
    comp = resolve_chemsep_compound(cas_or_name)
    if not comp or not comp.rpp_ig_cp:
        return None
    corr = comp.rpp_ig_cp
    A = corr.coeffs.get("A")
    B = corr.coeffs.get("B")
    C = corr.coeffs.get("C")
    D = corr.coeffs.get("D")
    if A is None or B is None or C is None or D is None:
        return None
    return {
        "source": "chemsep_xml",
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "Tmin_K": corr.t_min_k,
        "Tmax_K": corr.t_max_k,
        "units": corr.units,
        "eqno": corr.eqno,
    }


def get_antoine_coefficients(cas_or_name: str) -> Optional[Dict[str, Any]]:
    """
    Vapor pressure coefficients.
    Prefer ChemSep AntoineVaporPressure; fallback to Antoine_Coefficients.csv (name-based).
    """
    comp = resolve_chemsep_compound(cas_or_name)
    if comp and comp.antoine_vp:
        corr = comp.antoine_vp
        A = corr.coeffs.get("A")
        B = corr.coeffs.get("B")
        C = corr.coeffs.get("C")
        if A is not None and B is not None and C is not None:
            return {
                "source": "chemsep_xml",
                "A": A,
                "B": B,
                "C": C,
                "Tmin_K": corr.t_min_k,
                "Tmax_K": corr.t_max_k,
                "units": corr.units,
                "eqno": corr.eqno,
            }

    if comp and comp.name:
        row = _antoine_csv_db().get(comp.name)
        if row:
            return {
                "source": "antoine_csv",
                "A": row["A"],
                "B": row["B"],
                "C": row["C"],
                "Tmin_C": row.get("TMIN"),
                "Tmax_C": row.get("TMAX"),
            }
    return None


def get_transport_correlations(cas_or_name: str) -> Optional[Dict[str, Any]]:
    """
    Return raw viscosity and thermal conductivity correlations from ChemSep XML.

    This function intentionally returns **raw correlation blocks** (eqno + A..E)
    because the functional form depends on `eqno` and is not yet standardized
    in our IDAES property packages. For the detailed HX scripts we typically
    fix/estimate heat transfer coefficients directly (and thus do not require
    transport correlations to be embedded into the property package).
    """
    comp = resolve_chemsep_compound(cas_or_name)
    if not comp:
        return None

    def _corr_to_dict(c: Optional[ChemSepCorrelation]) -> Optional[Dict[str, Any]]:
        if c is None:
            return None
        return {
            "eqno": c.eqno,
            "coeffs": dict(c.coeffs),
            "Tmin_K": c.t_min_k,
            "Tmax_K": c.t_max_k,
            "units": c.units,
        }

    return {
        "liq_visc": _corr_to_dict(comp.liq_visc),
        "vap_visc": _corr_to_dict(comp.vap_visc),
        "liq_k": _corr_to_dict(comp.liq_k),
        "vap_k": _corr_to_dict(comp.vap_k),
    }


def get_pr_binary_interaction(cas1: str, cas2: str) -> Optional[float]:
    """Return PR k_ij (k12) from `ipd/pr.ipd` if present."""
    if not (CAS_RE.match(cas1) and CAS_RE.match(cas2)):
        return None
    return _ipd_db().get_pr_kij(cas1, cas2)


def check_property_availability(cas_or_name: str) -> Dict[str, Any]:
    comp = resolve_chemsep_compound(cas_or_name)
    if not comp:
        crit = get_critical_properties(cas_or_name)
        return {
            "found_in_chemsep_xml": False,
            "chemsep_name": None,
            "chemsep_cas": None,
            "critical": crit,
            "liquid_density": None,
            "liquid_cp": None,
            "ig_cp_rpp": None,
            "antoine": None,
            "transport": None,
        }

    return {
        "found_in_chemsep_xml": True,
        "chemsep_name": comp.name,
        "chemsep_cas": comp.cas,
        "critical": get_critical_properties(comp.cas or comp.name),
        "liquid_density": get_perrys_density_coefficients(comp.cas or comp.name),
        "liquid_cp": get_perrys_cp_liq_coefficients(comp.cas or comp.name),
        "ig_cp_rpp": get_rpp_ig_cp_coefficients(comp.cas or comp.name),
        "antoine": get_antoine_coefficients(comp.cas or comp.name),
        "transport_present": {
            "liq_visc_present": comp.liq_visc is not None,
            "vap_visc_present": comp.vap_visc is not None,
            "liq_k_present": comp.liq_k is not None,
            "vap_k_present": comp.vap_k is not None,
        },
        "transport_raw": get_transport_correlations(comp.cas or comp.name),
    }

