"""
Plots for the stair-cascade simulation.

Generates:
  fig_profiles.png   - water-surface (depth) profile on a tread, per width
  fig_summary.png    - transit time and hold-up vs channel width
  fig_cascade.png    - scaled side view of steps + nappe trajectories (W=5 cm)
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stair_cascade_sim as sim


def _out(name):
    return os.path.join(HERE, name)


def fig_profiles():
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for W in sim.WIDTHS:
        h = sim.step_hydraulics(W)
        # plot depth vs distance measured from the riser (flow direction)
        x_from_riser = h["wetted_len"] - h["profile_s"]   # 0 at landing, L at brink
        x_plot = x_from_riser + h["x_land"]
        ax.plot(x_plot * 100, h["profile_y"] * 1000,
                label=f"W = {W*100:.0f} cm  (q={h['q']*1e3:.2f} L/s/m)")
        ax.axhline(h["yc"] * 1000, ls=":", lw=0.8, color="grey")
    ax.set_xlabel("distance along tread from riser [cm]")
    ax.set_ylabel("flow depth y [mm]")
    ax.set_title("Steady water-surface profile on one tread\n"
                 "(free-overfall control: y -> critical depth at the brink)")
    ax.axvline(sim.L_STEP * 100, color="k", lw=1.0)
    ax.text(sim.L_STEP * 100 - 0.2, 0.3, "brink", rotation=90,
            va="bottom", ha="right", fontsize=8)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(_out("fig_profiles.png"), dpi=140)
    plt.close(fig)


def fig_summary():
    Ws = sim.WIDTHS
    res = [sim.analyse_width(W) for W in Ws]
    x = [W * 100 for W in Ws]
    transit = [r.transit_time_s for r in res]
    ymax = [r.y_max_mm for r in res]

    fig, ax1 = plt.subplots(figsize=(6.5, 4.2))
    c1, c2 = "tab:blue", "tab:red"
    ax1.plot(x, transit, "o-", color=c1)
    for xi, ti, hi in zip(x, transit, [r.holdup_mass_total_kg for r in res]):
        ax1.annotate(f"{ti:.0f} s\n{hi:.2f} kg", (xi, ti),
                     textcoords="offset points", xytext=(6, -18), fontsize=8)
    ax1.set_xlabel("channel width W [cm]")
    ax1.set_ylabel("top-to-bottom transit time [s]", color=c1)
    ax1.tick_params(axis="y", labelcolor=c1)
    ax1.set_xticks(x)
    ax1.grid(alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(x, ymax, "s--", color=c2)
    ax2.set_ylabel("max film depth on a tread [mm]", color=c2)
    ax2.tick_params(axis="y", labelcolor=c2)

    ax1.set_title("Transit time vs channel width (hold-up annotated)\n"
                  "narrower channel -> higher unit flow -> faster, thicker film")
    fig.tight_layout()
    fig.savefig(_out("fig_summary.png"), dpi=140)
    plt.close(fig)


def fig_cascade(W=0.05, nshow=6):
    """Scaled side view of the first nshow steps with nappe trajectories."""
    h = sim.step_hydraulics(W)
    Vb = h["Vc"]
    fig, ax = plt.subplots(figsize=(7.5, 4.5))

    # Draw steps (side profile). Top step at top-left.
    for i in range(nshow):
        x0 = i * sim.L_STEP
        z_tread = -(i * sim.H_STEP)
        # tread (horizontal)
        ax.plot([x0, x0 + sim.L_STEP], [z_tread, z_tread], "k-", lw=2)
        # riser (vertical) down to next tread
        ax.plot([x0 + sim.L_STEP, x0 + sim.L_STEP],
                [z_tread, z_tread - sim.H_STEP], "k-", lw=2)

        # water film on tread (depth profile), thin layer above tread
        x_from_riser = h["wetted_len"] - h["profile_s"]
        xs = x0 + h["x_land"] + x_from_riser
        ys = z_tread + h["profile_y"]
        ax.fill_between(xs, z_tread, ys, color="tab:blue", alpha=0.5)

        # nappe trajectory from this brink to next tread
        brink_x = x0 + sim.L_STEP
        brink_z = z_tread + h["yc"]
        tt = np.linspace(0, h["t_fall"], 30)
        nx = brink_x + Vb * tt
        nz = brink_z - 0.5 * sim.G * tt ** 2
        ax.plot(nx, nz, color="tab:cyan", lw=1.5)

    ax.set_aspect("equal")
    ax.set_xlabel("horizontal distance [m]")
    ax.set_ylabel("elevation [m]")
    ax.set_title(f"Cascade side view, W = {W*100:.0f} cm "
                 f"(first {nshow} of {sim.N_STEPS} steps)\n"
                 "blue = film on tread, cyan = free-falling nappe")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(_out("fig_cascade.png"), dpi=140)
    plt.close(fig)


def fig_scenarios():
    """Transit time vs width for the two flow rates, with regime labels."""
    scen = [(200.0, "200 kg/h (design)", "tab:green"),
            (3000.0, "3000 kg/h = 50 kg/min", "tab:purple")]
    fig, ax = plt.subplots(figsize=(7, 4.4))
    x = [W * 100 for W in sim.WIDTHS]
    for mdot_h, lab, col in scen:
        mdot = mdot_h / 3600.0
        rs = [sim.analyse_scenario_width(W, mdot) for W in sim.WIDTHS]
        t = [r["transit_time_s"] for r in rs]
        ax.plot(x, t, "o-", color=col, label=lab)
        for xi, ti, r in zip(x, t, rs):
            tag = "skim" if "skim" in r["model_used"] else "nappe"
            ax.annotate(f"{ti:.0f}s\n{tag}\n{r['holdup_total_kg']:.0f}kg",
                        (xi, ti), textcoords="offset points",
                        xytext=(6, 4), fontsize=7.5, color=col)
    ax.set_xlabel("channel width W [cm]")
    ax.set_ylabel("top-to-bottom transit time [s]")
    ax.set_xticks(x)
    ax.set_title("Transit time vs width at two flow rates\n"
                 "(15x higher flow -> faster transit, but ~5x more hold-up; "
                 "5 cm goes skimming)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(_out("fig_scenarios.png"), dpi=140)
    plt.close(fig)


def fig_roof_lamp(W=0.10, mdot_h=3000.0):
    """Cross-section across the channel: water envelope, lamp, roof clearance."""
    mdot = mdot_h / 3600.0
    r = sim.analyse_scenario_width(W, mdot)
    env = r["water_envelope_cm"]
    rc = sim.roof_clearance(env / 100.0, lamp_dia_m=0.04)
    lamp_d = 4.0
    gap_below = 3.0
    lamp_cy = env + gap_below + lamp_d / 2.0
    roof_y = rc["clearance_with_lamp_cm"]
    Wc = W * 100

    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    # channel walls + floor
    ax.plot([0, 0], [0, roof_y], "k-", lw=2)
    ax.plot([Wc, Wc], [0, roof_y], "k-", lw=2)
    ax.plot([0, Wc], [0, 0], "k-", lw=3)            # tread floor
    # water (bulked) + spray hatch
    dbulk = r.get("d_bulk_mm", r["y_max_mm"]) / 10.0
    ax.add_patch(plt.Rectangle((0, 0), Wc, dbulk, color="tab:blue", alpha=0.6))
    ax.add_patch(plt.Rectangle((0, dbulk), Wc, env - dbulk,
                               color="tab:cyan", alpha=0.3, hatch="..."))
    # lamp
    lamp = plt.Circle((Wc / 2.0, lamp_cy), lamp_d / 2.0,
                      color="gold", ec="darkgoldenrod", lw=1.5)
    ax.add_patch(lamp)
    # roof
    ax.plot([-1, Wc + 1], [roof_y, roof_y], color="dimgray", lw=3)
    ax.text(Wc / 2.0, roof_y + 0.3, "roof", ha="center", fontsize=9)

    # dimension annotations
    ax.annotate("", (Wc + 1.5, 0), (Wc + 1.5, env),
                arrowprops=dict(arrowstyle="<->", color="tab:blue"))
    ax.text(Wc + 2.0, env / 2, f"water+spray\n{env:.1f} cm",
            color="tab:blue", fontsize=8, va="center")
    ax.annotate("", (Wc + 1.5, env), (Wc + 1.5, lamp_cy - lamp_d / 2),
                arrowprops=dict(arrowstyle="<->", color="green"))
    ax.text(Wc + 2.0, env + gap_below / 2, f"gap {gap_below:.0f} cm",
            color="green", fontsize=8, va="center")
    ax.annotate("", (-1.5, 0), (-1.5, roof_y),
                arrowprops=dict(arrowstyle="<->", color="red"))
    ax.text(-2.0, roof_y / 2, f"clear height\n{roof_y:.1f} cm",
            color="red", fontsize=9, va="center", ha="right", rotation=90)
    ax.text(Wc / 2, lamp_cy, "lamp\n4 cm", ha="center", va="center", fontsize=7)

    ax.set_xlim(-6, Wc + 8)
    ax.set_ylim(-1, roof_y + 2)
    ax.set_aspect("equal")
    ax.set_xlabel("across channel [cm]")
    ax.set_ylabel("height above tread [cm]")
    ax.set_title(f"Roof clearance with overhead lamp\n"
                 f"W={Wc:.0f} cm, {mdot_h:.0f} kg/h, {r['regime'][:20]}")
    fig.tight_layout()
    fig.savefig(_out("fig_roof_lamp.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    fig_profiles()
    fig_summary()
    fig_cascade()
    fig_scenarios()
    fig_roof_lamp()
    print("Wrote fig_profiles, fig_summary, fig_cascade, fig_scenarios, "
          "fig_roof_lamp .png")
