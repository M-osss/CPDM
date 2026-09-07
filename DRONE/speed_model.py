#!/usr/bin/env python3
"""First-order top-speed model for a digital FPV quad.

Calibrated so a 5" 6S 2100 KV race reference (Mach R5 Ultra class:
2207 / ~5.1x3.6x3 / AUW 0.53 kg) lands near the published 245 km/h.

The model is for ranking configurations, not for claiming a GPS number.
Thrust uses a linear Ct(J) fade to zero at geometric pitch speed.
Drag is 0.5 * rho * CdA * V^2 with an induced/body term from mass.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from typing import Iterable

RHO = 1.225  # kg/m^3, sea level
A_SOUND = 343.0  # m/s
G = 9.80665
INCH = 0.0254


def rpm_no_load(kv: float, volts: float) -> float:
    return kv * volts


def pitch_speed_ms(rpm: float, pitch_in: float) -> float:
    """Geometric pitch speed, metres per second."""
    return (rpm / 60.0) * (pitch_in * INCH)


def tip_speed_ms(rpm: float, diameter_in: float) -> float:
    return math.pi * (diameter_in * INCH) * (rpm / 60.0)


def pack_voltage(cells: int, lihv: bool, sag_v: float) -> float:
    full = 4.35 if lihv else 4.20
    return cells * full - sag_v


@dataclass(frozen=True)
class QuadConfig:
    name: str
    size_in: float
    pitch_in: float
    blades: int
    cells: int
    kv: float
    motor_stator: str
    auw_kg: float
    cda_m2: float
    pack_mah: float
    pack_c: float
    lihv: bool = False
    motors: int = 4
    load_static: float = 0.58
    load_speed: float = 0.84
    sag_v: float = 1.4
    prop_efficiency: float = 0.72
    motor_efficiency: float = 0.78

    @property
    def volts_loaded(self) -> float:
        return pack_voltage(self.cells, self.lihv, self.sag_v)

    @property
    def rpm_static(self) -> float:
        return rpm_no_load(self.kv, self.volts_loaded) * self.load_static

    @property
    def rpm_loaded(self) -> float:
        """RPM at high advance ratio: the prop unloads and KV*V is approached."""
        return rpm_no_load(self.kv, self.volts_loaded) * self.load_speed

    @property
    def diameter_m(self) -> float:
        return self.size_in * INCH


@dataclass
class SpeedResult:
    name: str
    cells: int
    kv: float
    size_in: float
    pitch_in: float
    volts_loaded: float
    rpm_loaded: float
    tip_mach: float
    pitch_speed_kmh: float
    estimated_top_kmh: float
    static_thrust_n: float
    pack_current_a: float
    electrical_power_w: float
    i2r_relative: float
    burst_c_used: float
    tip_mach_warning: bool
    current_warning: bool
    notes: str


def _ct0(blades: int, pitch_ratio: float, size_in: float) -> float:
    """Static thrust coefficient, empirical, 2-blade vs 3-blade.

    Small discs run lower Reynolds numbers; 3" blades do not produce
    5"-class Ct. Scale from 5.1" reference.
    """
    blade_scale = 0.78 if blades == 2 else 1.0
    re_scale = min(1.0, (size_in / 5.1) ** 0.55)
    return 0.155 * blade_scale * re_scale * (0.70 + 0.30 * min(pitch_ratio, 1.35))


def _compressibility(tip_mach: float) -> float:
    """Thrust remaining as tips enter transonic. Hard loss above ~0.92."""
    if tip_mach < 0.85:
        return 1.0
    if tip_mach >= 1.05:
        return 0.55
    return 1.0 - 0.45 * (tip_mach - 0.85) / 0.20


def static_thrust_n(cfg: QuadConfig, rpm: float) -> float:
    d = cfg.diameter_m
    n = rpm / 60.0
    pitch_ratio = cfg.pitch_in / cfg.size_in
    ct = _ct0(cfg.blades, pitch_ratio, cfg.size_in)
    mach = tip_speed_ms(rpm, cfg.size_in) / A_SOUND
    per_motor = ct * RHO * n * n * d**4 * _compressibility(mach)
    return per_motor * cfg.motors


def thrust_at_speed_n(cfg: QuadConfig, rpm: float, speed_ms: float) -> float:
    """Linear Ct fade vs advance ratio; zero at geometric pitch speed."""
    d = cfg.diameter_m
    n = max(rpm / 60.0, 1e-6)
    j = speed_ms / (n * d)
    j_zero = max(cfg.pitch_in / cfg.size_in, 0.35)
    fade = max(0.0, 1.0 - j / j_zero)
    return static_thrust_n(cfg, rpm) * fade


def drag_n(cfg: QuadConfig, speed_ms: float) -> float:
    parasitic = 0.5 * RHO * cfg.cda_m2 * speed_ms * speed_ms
    # Small residual induced/body term so hover-ish speeds are not zero-drag.
    induced = (cfg.auw_kg * G) ** 2 / (
        0.5 * RHO * math.pi * (cfg.diameter_m * 0.5) ** 2 * cfg.motors * max(speed_ms, 8.0) ** 2
        + 1.0
    )
    return parasitic + 0.15 * induced


def solve_top_speed_ms(cfg: QuadConfig, rpm: float) -> float:
    lo = 0.0
    hi = pitch_speed_ms(rpm, cfg.pitch_in)
    if hi <= 0:
        return 0.0
    # If already drag-limited below zero speed (should not happen), return 0.
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if thrust_at_speed_n(cfg, rpm, mid) >= drag_n(cfg, mid):
            lo = mid
        else:
            hi = mid
    return lo


def electrical_power_w(cfg: QuadConfig, rpm: float, speed_ms: float) -> float:
    """Shaft power from momentum/prop scaling, then electrical after efficiencies."""
    d = cfg.diameter_m
    n = rpm / 60.0
    pitch_ratio = cfg.pitch_in / cfg.size_in
    cp = 0.038 * (0.65 + 0.55 * pitch_ratio) * (0.82 if cfg.blades == 2 else 1.0)
    # Advance-ratio unload: power drops as J rises.
    j_zero = max(pitch_ratio, 0.35)
    j = speed_ms / (n * d) if n * d > 0 else 0.0
    cp_eff = cp * max(0.25, 1.0 - 0.55 * min(j / j_zero, 1.0))
    shaft = cfg.motors * cp_eff * RHO * n**3 * d**5
    return shaft / (cfg.prop_efficiency * cfg.motor_efficiency)


def analyze(cfg: QuadConfig) -> SpeedResult:
    rpm = cfg.rpm_loaded
    v_top = solve_top_speed_ms(cfg, rpm)
    power = electrical_power_w(cfg, cfg.rpm_static, min(v_top, 25.0))
    current = power / max(cfg.volts_loaded, 0.1)
    # I^2R proxy vs a 6S reference at same power: current^2 (resistance ~ 1/cells-ish
    # but pack+ESC+leads are the practical bottleneck; use I^2 only).
    i2r = current * current
    pack_amp_budget = (cfg.pack_mah / 1000.0) * cfg.pack_c
    burst_c = current / max(cfg.pack_mah / 1000.0, 0.01)
    tip_m = tip_speed_ms(rpm, cfg.size_in) / A_SOUND
    notes = []
    if tip_m >= 0.92:
        notes.append("tip Mach in compressibility; unload with 2-blade or drop KV")
    if current > 0.85 * pack_amp_budget:
        notes.append("pack C-rating is the limiter; use a higher-C or larger mAh pack")
    if cfg.cells == 4:
        notes.append("4S pays a current/I2R penalty vs 6S at the same power")
    if cfg.size_in <= 3.5:
        notes.append("disc area limited; absolute top speed trails a 5 inch")
    return SpeedResult(
        name=cfg.name,
        cells=cfg.cells,
        kv=cfg.kv,
        size_in=cfg.size_in,
        pitch_in=cfg.pitch_in,
        volts_loaded=round(cfg.volts_loaded, 2),
        rpm_loaded=round(rpm, 0),
        tip_mach=round(tip_m, 3),
        pitch_speed_kmh=round(pitch_speed_ms(rpm, cfg.pitch_in) * 3.6, 1),
        estimated_top_kmh=round(v_top * 3.6, 1),
        static_thrust_n=round(static_thrust_n(cfg, cfg.rpm_static), 1),
        pack_current_a=round(current, 1),
        electrical_power_w=round(power, 0),
        i2r_relative=round(i2r, 0),
        burst_c_used=round(burst_c, 1),
        tip_mach_warning=tip_m >= 0.92,
        current_warning=current > 0.85 * pack_amp_budget,
        notes="; ".join(notes),
    )


# CdA: calibrated so the 2100 KV / 5.1x3.6x3 / 6S reference lands near 245 km/h.
RACE_CDA = 0.0054
AERO_CDA = 0.0047
BRICK_3IN_CDA = 0.0046


def catalogue() -> list[QuadConfig]:
    """Candidate builds. The recommended row is `speed5_6s_hdzero`."""
    return [
        QuadConfig(
            name="ref_mach_r5_class",
            size_in=5.1,
            pitch_in=3.6,
            blades=3,
            cells=6,
            kv=2100,
            motor_stator="2207",
            auw_kg=0.527,
            cda_m2=RACE_CDA,
            pack_mah=1480,
            pack_c=120,
            load_static=0.60,
            load_speed=0.88,
            sag_v=1.0,
        ),
        QuadConfig(
            name="speed5_6s_hdzero",
            size_in=5.1,
            pitch_in=5.1,
            blades=2,
            cells=6,
            kv=2550,
            motor_stator="2207",
            auw_kg=0.515,
            cda_m2=AERO_CDA,
            pack_mah=1050,
            pack_c=150,
            lihv=True,
            load_static=0.52,
            load_speed=0.78,
            sag_v=1.8,
        ),
        QuadConfig(
            name="speed5_6s_2100kv_safer",
            size_in=5.1,
            pitch_in=4.8,
            blades=2,
            cells=6,
            kv=2100,
            motor_stator="2207",
            auw_kg=0.515,
            cda_m2=AERO_CDA,
            pack_mah=1050,
            pack_c=150,
            load_static=0.58,
            load_speed=0.86,
            sag_v=1.3,
        ),
        QuadConfig(
            name="speed5_4s_ultrahigh_kv",
            size_in=5.1,
            pitch_in=5.1,
            blades=2,
            cells=4,
            kv=3800,
            motor_stator="2207",
            auw_kg=0.495,
            cda_m2=AERO_CDA,
            pack_mah=1300,
            pack_c=150,
            lihv=True,
            load_static=0.50,
            load_speed=0.76,
            sag_v=1.6,
        ),
        QuadConfig(
            name="speed3_6s",
            size_in=3.0,
            pitch_in=4.0,
            blades=2,
            cells=6,
            kv=4200,
            motor_stator="1404",
            auw_kg=0.250,
            cda_m2=BRICK_3IN_CDA,
            pack_mah=750,
            pack_c=120,
            load_static=0.55,
            load_speed=0.80,
            sag_v=1.5,
        ),
        QuadConfig(
            name="speed3.5_6s",
            size_in=3.5,
            pitch_in=4.5,
            blades=2,
            cells=6,
            kv=3400,
            motor_stator="1606",
            auw_kg=0.310,
            cda_m2=0.0046,
            pack_mah=850,
            pack_c=130,
            load_static=0.54,
            load_speed=0.80,
            sag_v=1.6,
        ),
        QuadConfig(
            name="speed5.5_6s_lowkv",
            size_in=5.5,
            pitch_in=4.3,
            blades=3,
            cells=6,
            kv=1800,
            motor_stator="2208",
            auw_kg=0.560,
            cda_m2=0.0052,
            pack_mah=1100,
            pack_c=120,
            load_static=0.60,
            load_speed=0.86,
            sag_v=1.3,
        ),
    ]


def recommend(results: Iterable[SpeedResult]) -> SpeedResult:
    ranked = sorted(
        results,
        key=lambda r: (r.estimated_top_kmh, -r.pack_current_a),
        reverse=True,
    )
    return ranked[0]


def as_table(results: list[SpeedResult]) -> str:
    headers = [
        "name",
        "S",
        "KV",
        "in",
        "pitch",
        "rpm",
        "Mach",
        "pitch_kmh",
        "top_kmh",
        "I_A",
        "P_W",
        "C_used",
    ]
    rows = [headers]
    for r in results:
        rows.append(
            [
                r.name,
                str(r.cells),
                str(int(r.kv)),
                str(r.size_in),
                str(r.pitch_in),
                str(int(r.rpm_loaded)),
                f"{r.tip_mach:.2f}",
                f"{r.pitch_speed_kmh:.0f}",
                f"{r.estimated_top_kmh:.0f}",
                f"{r.pack_current_a:.0f}",
                f"{r.electrical_power_w:.0f}",
                f"{r.burst_c_used:.0f}",
            ]
        )
    widths = [max(len(row[i]) for row in rows) for i in range(len(headers))]
    lines = []
    for i, row in enumerate(rows):
        line = "  ".join(cell.ljust(widths[j]) for j, cell in enumerate(row))
        lines.append(line)
        if i == 0:
            lines.append("  ".join("-" * w for w in widths))
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare digital speed-quad configs")
    parser.add_argument("--json", action="store_true", help="emit JSON instead of a table")
    args = parser.parse_args()
    results = [analyze(cfg) for cfg in catalogue()]
    if args.json:
        payload = {
            "recommended": recommend(results).name,
            "results": [asdict(r) for r in results],
        }
        print(json.dumps(payload, indent=2))
        return
    print(as_table(results))
    rec = recommend(results)
    print()
    print(f"highest modelled top speed: {rec.name} ({rec.estimated_top_kmh:.0f} km/h)")
    winner = next(r for r in results if r.name == "speed5_6s_hdzero")
    four = next(r for r in results if r.name == "speed5_4s_ultrahigh_kv")
    print(
        f"6S vs 4S at 5 inch: {winner.estimated_top_kmh:.0f} km/h @ {winner.pack_current_a:.0f} A "
        f"vs {four.estimated_top_kmh:.0f} km/h @ {four.pack_current_a:.0f} A; "
        f"I2R ratio 4S/6S = {four.i2r_relative / winner.i2r_relative:.2f}"
    )


if __name__ == "__main__":
    main()
