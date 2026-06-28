"""
Stair-cascade open-channel flow simulation
===========================================

Physical problem
----------------
Water (rho ~ 1000 kg/m3) is recirculated from a 400 kg tank by a pump that
delivers 200 kg/h and lifts it 1 m. From the top it flows down a stainless-steel
staircase of N = 33 steps. Each step: riser height h = 3 cm, tread length
l = 20 cm. Channel width W is 5, 10 or 15 cm (separate runs).

Goal
----
1. Classify the flow regime (nappe vs skimming).
2. Estimate the time for water to travel from the top step to the bottom step.
3. Estimate the liquid hold-up on the staircase and the required roof / wall
   clearance ("freeboard").
4. Support sizing / pitfall analysis of the pump and the recirculation loop.

Hydraulic model
---------------
For these flow rates yc/h << 1, so the flow is in the *nappe* regime: water
cascades from step to step as a thin film + free-falling nappe rather than
skimming over the step edges. Each horizontal tread is therefore treated as a
short horizontal open channel that ends in a *free overfall* (brink). The
control section is critical depth at the brink. The water-surface profile on
the tread is obtained by integrating the gradually-varied-flow (GVF) equation
for a horizontal bed (bed slope S0 = 0):

    dy/dx = -Sf / (1 - Fr^2)         (x measured in the flow direction)

    Sf  = n^2 * (q/y)^2 / R^(4/3)     (Manning friction slope)
    Fr^2 = q^2 / (g * y^3)            (rectangular channel)
    R   = W*y / (W + 2y)             (hydraulic radius, rectangular)

with q = Q/W the discharge per unit width. The profile is integrated *upstream*
from the brink (critical depth) over the wetted part of the tread.

The nappe leaves each brink horizontally at the brink velocity and falls h to
the next tread (projectile motion), which fixes the landing point and the
wetted tread length, and the per-step free-fall time.

The transit time of a water parcel across one tread equals the residence time
on that tread:  t_tread = integral(dx / V) = integral(y dx) / q.
Total transit time = sum over steps of (t_tread + t_fall).

All correlations are documented inline with their assumptions/limitations.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, field, asdict

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))


# --------------------------------------------------------------------------- #
# Fixed inputs (SI units unless noted)
# --------------------------------------------------------------------------- #
G = 9.81                      # gravitational acceleration [m/s^2]

# Water properties at ~20 C
RHO = 998.2                   # density [kg/m3]
MU = 1.002e-3                 # dynamic viscosity [Pa.s]
NU = MU / RHO                 # kinematic viscosity [m2/s]
SIGMA = 0.0728                # surface tension water/air [N/m]

# System
M_TANK = 400.0                # tank inventory [kg]
MDOT = 200.0 / 3600.0         # mass flow [kg/s]  (200 kg/h)
LIFT = 1.0                    # pump lift [m]

# Staircase geometry
N_STEPS = 33
H_STEP = 0.03                 # riser height [m]
L_STEP = 0.20                 # tread length [m]
WIDTHS = [0.05, 0.10, 0.15]   # channel widths to test [m]

# Stainless steel: hydraulically smooth metal.
# Manning n for clean smooth metal/steel ~ 0.011-0.013. Use 0.012 nominal.
MANNING_N = 0.012
N_RANGE = (0.011, 0.013)      # for sensitivity

Q_VOL = MDOT / RHO            # volumetric flow [m3/s]


# --------------------------------------------------------------------------- #
# Elementary hydraulics
# --------------------------------------------------------------------------- #
def critical_depth(q: float) -> float:
    """Critical depth of a wide rectangular channel, yc = (q^2/g)^(1/3)."""
    return (q * q / G) ** (1.0 / 3.0)


def hydraulic_radius(y: float, W: float) -> float:
    """Rectangular channel hydraulic radius R = W*y/(W + 2y)."""
    return W * y / (W + 2.0 * y)


def froude_sq(y: float, q: float) -> float:
    """Froude number squared for rectangular channel, Fr^2 = q^2/(g y^3)."""
    return q * q / (G * y ** 3)


def friction_slope(y: float, q: float, W: float, n: float) -> float:
    """Manning friction slope Sf = n^2 (q/y)^2 / R^(4/3)."""
    R = hydraulic_radius(y, W)
    V = q / y
    return n * n * V * V / R ** (4.0 / 3.0)


def reynolds_film(W: float) -> float:
    """Falling-film Reynolds number Re = 4*Gamma/mu, Gamma = mdot/W."""
    gamma = MDOT / W
    return 4.0 * gamma / MU


def wetting_rate(W: float) -> float:
    """Mass flow per unit width Gamma [kg/(m.s)]."""
    return MDOT / W


# --------------------------------------------------------------------------- #
# Gradually-varied-flow profile on a single horizontal tread
# --------------------------------------------------------------------------- #
def tread_profile(q: float, W: float, n: float, wetted_len: float,
                  npts: int = 400):
    """
    Integrate the horizontal-bed GVF equation upstream from the brink.

    Returns arrays (s, y) where s is distance measured *upstream* from the
    brink (s = 0 at brink, s = wetted_len at the landing point) and y is the
    flow depth. Critical depth is imposed as the brink control.
    """
    yc = critical_depth(q)
    # Start just on the subcritical side of critical to avoid the 1/(1-Fr^2)
    # singularity exactly at the brink.
    y0 = yc * (1.0 + 1.0e-3)

    def rhs(s, y):
        yy = max(y[0], 1e-9)
        Sf = friction_slope(yy, q, W, n)
        denom = 1.0 - froude_sq(yy, q)
        # dy/ds = -dy/dx ; on horizontal bed dy/dx = -Sf/(1-Fr^2)
        return [Sf / denom]

    s_eval = np.linspace(0.0, wetted_len, npts)
    sol = solve_ivp(rhs, (0.0, wetted_len), [y0], t_eval=s_eval,
                    method="RK45", rtol=1e-8, atol=1e-12, max_step=wetted_len / 50)
    s = sol.t
    y = sol.y[0]
    return s, y


def step_hydraulics(W: float, n: float = MANNING_N):
    """
    Full per-step hydraulic summary for a given channel width.
    """
    q = Q_VOL / W                       # unit discharge [m2/s]
    yc = critical_depth(q)
    Vc = q / yc                         # critical velocity = brink velocity (control)

    # Nappe projectile: leaves brink horizontally at Vb=Vc, falls h.
    t_fall = math.sqrt(2.0 * H_STEP / G)
    x_land = Vc * t_fall                 # horizontal landing distance from riser
    wetted_len = max(L_STEP - x_land, 1e-3)

    # Impact velocity (for splash height estimate)
    v_vert = math.sqrt(2.0 * G * H_STEP)
    v_impact = math.hypot(Vc, v_vert)
    splash_h = v_impact ** 2 / (2.0 * G)

    # GVF profile on the wetted tread
    s, y = tread_profile(q, W, n, wetted_len)

    # Hold-up per width = integral(y ds); residence/transit time = that / q
    holdup_per_width = np.trapezoid(y, s)        # [m2]
    t_tread = holdup_per_width / q               # [s]
    volume_step = holdup_per_width * W           # [m3]
    mass_step = volume_step * RHO                # [kg]
    y_max = float(np.max(y))                     # deepest point on tread [m]
    y_mean = holdup_per_width / wetted_len

    overshoot = x_land > L_STEP                  # nappe overshoots tread?

    return {
        "q": q,
        "yc": yc,
        "Vc": Vc,
        "t_fall": t_fall,
        "x_land": x_land,
        "wetted_len": wetted_len,
        "overshoot": overshoot,
        "v_impact": v_impact,
        "splash_h": splash_h,
        "holdup_per_width": holdup_per_width,
        "t_tread": t_tread,
        "volume_step": volume_step,
        "mass_step": mass_step,
        "y_max": y_max,
        "y_mean": y_mean,
        "profile_s": s,
        "profile_y": y,
    }


# --------------------------------------------------------------------------- #
# Flow-regime classification (nappe vs skimming)
# --------------------------------------------------------------------------- #
def regime(W: float):
    """
    Classify the cascade flow regime using Chanson's stepped-spillway criteria.

    Skimming flow onset:  (yc/h)_onset = 1.057 - 0.465 (h/l)   [Chanson 1994]
    Nappe flow if yc/h < onset (with a transition band just below).
    """
    q = Q_VOL / W
    yc = critical_depth(q)
    ratio = yc / H_STEP
    onset = 1.057 - 0.465 * (H_STEP / L_STEP)
    if ratio >= onset:
        label = "skimming"
    elif ratio >= 0.8 * onset:
        label = "transition"
    else:
        label = "nappe"
    return {"yc_over_h": ratio, "skimming_onset": onset, "regime": label}


# --------------------------------------------------------------------------- #
# Full system roll-up for one width
# --------------------------------------------------------------------------- #
@dataclass
class WidthResult:
    width_cm: float
    q: float
    Re_film: float
    Gamma: float
    regime: str
    yc_over_h: float
    yc_mm: float
    Vc: float
    t_tread: float
    t_fall: float
    transit_time_s: float
    holdup_mass_total_kg: float
    holdup_frac_of_tank: float
    y_max_mm: float
    splash_h_cm: float
    wetted_len_cm: float
    nappe_overshoot: bool
    freeboard_min_cm: float
    profile_s: list = field(default_factory=list, repr=False)
    profile_y: list = field(default_factory=list, repr=False)


def analyse_width(W: float, n: float = MANNING_N) -> WidthResult:
    h = step_hydraulics(W, n)
    reg = regime(W)

    transit = N_STEPS * (h["t_tread"] + h["t_fall"])
    holdup_mass = N_STEPS * h["mass_step"]

    # Freeboard / roof clearance: deepest water + splash rise + 50% margin
    freeboard = (h["y_max"] + h["splash_h"]) * 1.5

    return WidthResult(
        width_cm=W * 100,
        q=h["q"],
        Re_film=reynolds_film(W),
        Gamma=wetting_rate(W),
        regime=reg["regime"],
        yc_over_h=reg["yc_over_h"],
        yc_mm=h["yc"] * 1000,
        Vc=h["Vc"],
        t_tread=h["t_tread"],
        t_fall=h["t_fall"],
        transit_time_s=transit,
        holdup_mass_total_kg=holdup_mass,
        holdup_frac_of_tank=holdup_mass / M_TANK,
        y_max_mm=h["y_max"] * 1000,
        splash_h_cm=h["splash_h"] * 100,
        wetted_len_cm=h["wetted_len"] * 100,
        nappe_overshoot=h["overshoot"],
        freeboard_min_cm=freeboard * 100,
        profile_s=h["profile_s"].tolist(),
        profile_y=h["profile_y"].tolist(),
    )


# --------------------------------------------------------------------------- #
# Pump / recirculation analysis
# --------------------------------------------------------------------------- #
def pump_and_recirculation(passes: int = 5):
    Q = Q_VOL
    P_hyd = RHO * G * Q * LIFT                      # hydraulic power [W]
    eff = 0.35                                      # small-pump efficiency guess
    P_shaft = P_hyd / eff
    turnover_h = M_TANK / (MDOT * 3600.0)           # hours per pass
    total_h = passes * turnover_h
    # Temperature rise if all shaft power dumped to the 400 kg (worst case)
    cp = 4182.0
    dT_per_pass = P_shaft * (turnover_h * 3600.0) / (M_TANK * cp)
    return {
        "Q_m3_per_h": Q * 3600.0,
        "Q_L_per_min": Q * 1000.0 * 60.0,
        "hydraulic_power_W": P_hyd,
        "shaft_power_W_est": P_shaft,
        "turnover_h_per_pass": turnover_h,
        "passes": passes,
        "total_runtime_h": total_h,
        "dT_per_pass_C_worstcase": dT_per_pass,
    }


def skimming_onset_flow(W: float):
    """
    Mass flow (kg/h) at which yc/h reaches the skimming-flow onset for this
    geometry, i.e. the upper bound of the nappe regime. Above this the cascade
    would skim instead of cascade step-by-step.
    """
    onset = 1.057 - 0.465 * (H_STEP / L_STEP)
    yc = onset * H_STEP
    q = math.sqrt(G * yc ** 3)          # invert yc=(q^2/g)^(1/3)
    Q = q * W
    return Q * RHO * 3600.0             # kg/h


def min_wetting_rate():
    """
    Minimum wetting rate Gamma_min [kg/(m.s)] needed to sustain a continuous
    (non-rivulet) film. Closed-form values (Hartley & Murgatroyd 1964) depend
    strongly on contact angle and are uncertain, so the widely quoted
    experimental band for water on clean metal is returned as the anchor.
    Gamma below this band -> film breaks into rivulets / dry stripes.
    """
    return (0.07, 0.30)   # kg/(m.s), water on metal, typical experimental band


def evaporation_estimate(W: float):
    """
    Crude evaporation loss from exposed film area.
    Wetted tread area = N * l * W (ignoring risers/nappe; conservative-low).
    Use an indoor still-air mass-transfer coefficient ~ 1e-3 kg/m2/s/(kg/kg)
    is too high; instead use an engineering pan-style rate of ~0.05 kg/m2/h at
    20 C still air as a low anchor, and ~0.3 kg/m2/h with mild airflow.
    """
    area = N_STEPS * L_STEP * W
    rate_still = 0.05 * area      # kg/h
    rate_airflow = 0.3 * area     # kg/h
    return {"film_area_m2": area,
            "evap_kg_per_h_still": rate_still,
            "evap_kg_per_h_airflow": rate_airflow}


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    results = [analyse_width(W) for W in WIDTHS]

    # Sensitivity of transit time to Manning n at W = 0.05 m
    sens = {}
    for nval in N_RANGE:
        r = analyse_width(0.05, nval)
        sens[nval] = r.transit_time_s

    pump = pump_and_recirculation(passes=5)
    evap = {f"{int(W*100)}cm": evaporation_estimate(W) for W in WIDTHS}

    out = {
        "inputs": {
            "M_tank_kg": M_TANK,
            "mdot_kg_per_h": 200.0,
            "lift_m": LIFT,
            "N_steps": N_STEPS,
            "riser_h_cm": H_STEP * 100,
            "tread_l_cm": L_STEP * 100,
            "widths_cm": [W * 100 for W in WIDTHS],
            "manning_n": MANNING_N,
            "rho": RHO, "mu": MU, "sigma": SIGMA,
        },
        "per_width": [
            {k: v for k, v in asdict(r).items()
             if k not in ("profile_s", "profile_y")}
            for r in results
        ],
        "manning_sensitivity_W5cm_s": sens,
        "pump": pump,
        "evaporation": evap,
    }

    # Console report
    print("=" * 74)
    print("STAIR-CASCADE FLUID SIMULATION  -  water, stainless-steel staircase")
    print("=" * 74)
    print(f"Tank {M_TANK:.0f} kg | flow 200 kg/h ({pump['Q_L_per_min']:.2f} L/min, "
          f"{pump['Q_m3_per_h']:.3f} m3/h) | {N_STEPS} steps "
          f"({H_STEP*100:.0f} cm x {L_STEP*100:.0f} cm)\n")

    hdr = (f"{'W [cm]':>7}{'regime':>11}{'yc/h':>8}{'yc[mm]':>8}"
           f"{'Re_film':>9}{'transit[s]':>11}{'holdup[kg]':>11}"
           f"{'ymax[mm]':>9}{'freeboard[cm]':>14}")
    print(hdr)
    print("-" * len(hdr))
    for r in results:
        print(f"{r.width_cm:>7.0f}{r.regime:>11}{r.yc_over_h:>8.3f}"
              f"{r.yc_mm:>8.2f}{r.Re_film:>9.0f}{r.transit_time_s:>11.1f}"
              f"{r.holdup_mass_total_kg:>11.3f}{r.y_max_mm:>9.2f}"
              f"{r.freeboard_min_cm:>14.1f}")

    print("\nNappe per step: t_fall = %.3f s, "
          "landing %.1f-%.1f cm from riser (wetted tread %.1f-%.1f cm)"
          % (results[0].t_fall,
             min(L_STEP*100 - r.wetted_len_cm for r in results),
             max(L_STEP*100 - r.wetted_len_cm for r in results),
             min(r.wetted_len_cm for r in results),
             max(r.wetted_len_cm for r in results)))

    print("\nManning-n sensitivity of transit time (W=5 cm):")
    for nval, t in sens.items():
        print(f"   n={nval}: {t:.1f} s")

    print("\nPump / recirculation:")
    print(f"   Hydraulic power  : {pump['hydraulic_power_W']:.2f} W")
    print(f"   Shaft power (est): {pump['shaft_power_W_est']:.2f} W "
          f"(eta~0.35)")
    print(f"   Tank turnover    : {pump['turnover_h_per_pass']:.2f} h per pass")
    print(f"   {pump['passes']} passes      : {pump['total_runtime_h']:.1f} h")
    print(f"   Worst-case dT/pass from pump heat: "
          f"{pump['dT_per_pass_C_worstcase']:.4f} C")

    print("\nEvaporation (exposed tread film area):")
    for k, v in evap.items():
        print(f"   {k:>5}: area {v['film_area_m2']:.3f} m2 -> "
              f"{v['evap_kg_per_h_still']:.3f}-{v['evap_kg_per_h_airflow']:.3f} kg/h")

    print("\nRegime headroom and wetting:")
    g_lo, g_hi = min_wetting_rate()
    for r in results:
        W = r.width_cm / 100.0
        on = skimming_onset_flow(W)
        flag = "OK" if r.Gamma > g_hi else ("marginal" if r.Gamma > g_lo
                                            else "RIVULETS")
        print(f"   W={r.width_cm:>2.0f} cm: Gamma={r.Gamma:.3f} kg/m/s "
              f"({flag} vs wetting band {g_lo}-{g_hi}); "
              f"skimming onset at {on:,.0f} kg/h "
              f"(x{on/200.0:.0f} the design flow)")
    out["regime_headroom"] = {
        f"{int(r.width_cm)}cm": {
            "Gamma_kg_per_m_s": r.Gamma,
            "skimming_onset_kg_per_h": skimming_onset_flow(r.width_cm / 100.0),
        } for r in results
    }
    out["min_wetting_band_kg_per_m_s"] = list(min_wetting_rate())

    # Save JSON (without big profile arrays in the per_width block; keep separate)
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(out, f, indent=2)
    # Save profiles separately for plotting
    np.savez(os.path.join(HERE, "profiles.npz"),
             **{f"s_{int(r.width_cm)}": np.array(r.profile_s) for r in results},
             **{f"y_{int(r.width_cm)}": np.array(r.profile_y) for r in results})

    print("\nWrote results.json and profiles.npz to", HERE)
    return out, results


if __name__ == "__main__":
    main()
