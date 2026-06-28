# StairCascade — water cascade over a stepped channel

Open-channel hydraulic model of water recirculated from a 400 kg tank, lifted
1 m by a pump, and cascaded down a 33-step stainless-steel staircase. Computes
flow regime, top-to-bottom transit time, liquid hold-up, required roof/wall
clearance, and pump / recirculation sizing for channel widths of 5, 10 and
15 cm.

## Files

- `stair_cascade_sim.py` — model and analysis (run this).
- `make_plots.py` — figures (`fig_profiles.png`, `fig_summary.png`, `fig_cascade.png`).
- `results.json` — machine-readable results.
- `requirements.txt` — `numpy`, `scipy`, `matplotlib`.

```
pip install -r StairCascade/requirements.txt
python StairCascade/stair_cascade_sim.py
python StairCascade/make_plots.py
```

## Inputs

| quantity | value |
|---|---|
| tank inventory | 400 kg |
| mass flow | 200 kg/h = 0.0556 kg/s = 3.34 L/min = 0.200 m³/h |
| pump lift | 1 m |
| steps | 33 |
| riser height `h` | 3 cm |
| tread length `l` | 20 cm |
| widths `W` | 5, 10, 15 cm |
| material | stainless steel, Manning `n` ≈ 0.012 |
| fluid | water at 20 °C (ρ = 998 kg/m³, μ = 1.00 mPa·s, σ = 0.073 N/m) |

Total vertical drop 0.99 m; total horizontal run 6.6 m; pseudo-slope 8.5°.

## Physics model

### 1. Flow regime — this is a *nappe* (cascade) flow, not skimming

The discharge per unit width is `q = Q/W`; critical depth `yc = (q²/g)^(1/3)`.
The nappe/skimming boundary for stepped channels (Chanson 1994) is
`(yc/h)_onset = 1.057 − 0.465·(h/l) = 0.99` here.

| W [cm] | q [L/s per m] | yc [mm] | yc/h | regime |
|---|---|---|---|---|
| 5  | 1.11 | 5.02 | 0.167 | nappe |
| 10 | 0.56 | 3.16 | 0.105 | nappe |
| 15 | 0.37 | 2.41 | 0.080 | nappe |

`yc/h` is far below 0.99 in every case, so water cascades step-by-step: a thin
film crosses each horizontal tread, springs off the edge (brink) as a free
nappe, and falls 3 cm onto the next tread. This is the controlling picture for
the whole analysis.

### 2. Per-tread water surface — gradually varied flow with a free overfall

Each tread is horizontal (bed slope `S0 = 0`) and ends in a brink. The control
section is critical depth at the brink. The depth profile is obtained by
integrating the gradually-varied-flow equation upstream from the brink:

```
dy/dx = −Sf / (1 − Fr²)
Sf  = n²·(q/y)² / R^(4/3)      (Manning friction slope)
Fr² = q² / (g·y³)             (rectangular channel)
R   = W·y / (W + 2y)
```

Depth draws down to `yc` at the brink and is deeper (sub-critical) upstream
(see `fig_profiles.png`). The deepest film on any tread is 3.9–7.2 mm, well
below the 30 mm riser, so water never overtops a step.

### 3. Nappe trajectory between steps

The nappe leaves the brink horizontally at the brink velocity `Vb ≈ Vc = q/yc`
and falls `h = 3 cm`:

```
t_fall = sqrt(2h/g) = 0.078 s
x_land = Vb·t_fall   = 1.2–1.7 cm from the riser
```

The jet lands 1–2 cm onto the next tread (no overshoot of the 20 cm tread), so
the cascade is well-formed and each tread is wetted over ≈ 18.3–18.8 cm.

### 4. Transit time = residence time

The time for a water parcel to cross one tread equals the local
residence time, `t_tread = ∫ dx/V = (∫ y dx)/q`. The total top-to-bottom
transit time is

```
T = Σ_steps ( t_tread + t_fall )
```

This is identical to (system hold-up)/(throughput), which is the physically
correct "time water spends on the staircase" at steady operation.

## Results

| W [cm] | transit time T | hold-up on stairs | max film depth | film Re |
|---|---|---|---|---|
| 5  | **38 s** | 1.96 kg (0.49 % of tank) | 7.2 mm | 4400 |
| 10 | **50 s** | 2.65 kg (0.66 % of tank) | 4.8 mm | 2200 |
| 15 | **60 s** | 3.18 kg (0.79 % of tank) | 3.9 mm | 1500 |

Trend: a **narrower** channel has a **higher** unit discharge → higher velocity
→ **shorter** transit and **less** hold-up (but a thicker film). Widening the
channel slows the water and increases hold-up.

The transit-time answer is robust: varying Manning `n` over the full
stainless-steel range (0.011–0.013) changes the 5 cm transit time by only
±0.6 s (37.3–38.5 s). The leading-edge (first-droplet) arrival on an initially
dry staircase is the same order of magnitude (tens of seconds); steady
operation is reached within ~1–2 transit times (≈1–2 min).

## Roof / wall clearance ("accommodate all the water")

Two distinct concerns:

1. **Hold-up cannot flood the enclosure.** Total liquid on the 33 steps is only
   2–3 kg (≤0.8 % of the 400 kg). The tank always holds ≥397 kg. The staircase
   physically cannot accumulate a meaningful fraction of the inventory at
   200 kg/h, so there is no flooding risk and no need to size the roof for bulk
   water. The deepest film (7 mm) sits far below the 30 mm risers.

2. **Splash / nappe clearance.** Each impacting nappe throws spray to roughly
   `v_impact²/(2g) ≈ 3.1–3.3 cm` above the tread. Adding the film depth and a
   50 % margin gives a per-tread envelope of ≈ 5–6 cm.

   **Recommended interior clear height ≥ 10 cm above each tread** (measured
   perpendicular to the treads), with the roof running parallel to the 8.5°
   stair slope. The extra margin over the 6 cm physics envelope:
   - keeps the roof clear of splash and wavy-film crests,
   - lets the air cavity under each nappe stay ventilated (an unventilated
     nappe clings to the riser and changes the flow),
   - allows cleaning/inspection.

   **Inlet caution:** if water is released from the full 1 m lift *directly
   above* the top step it strikes at √(2g·1) ≈ 4.4 m/s and splashes ~1 m high.
   Deliver it tangentially at top-step level (or via a short stilling
   box/weir). Only the head end then needs extra height, not the whole roof.

## Pump selection

Duty point: **0.200 m³/h against 1 m head** → hydraulic power
`ρgQH ≈ 0.55 W`, shaft power ≈ 1.5 W at ~35 % efficiency. This is a *tiny,
low-head, low-flow* duty.

- A **centrifugal pump is a poor match**. Catalogue centrifugals have best
  efficiency at far higher flow; at 0.2 m³/h one runs near shut-off, where it
  suffers internal recirculation, vibration and self-heating, and the delivered
  flow is sensitive to small head changes — bad for a process that needs a
  *steady, metered* 200 kg/h. Forcing it to work needs a continuous-recirculation
  (min-flow) bypass, which wastes energy and adds complexity.
- **Recommended: a positive-displacement pump** — a peristaltic or small gear /
  diaphragm pump. At 3.3 L/min these deliver a near-constant, head-insensitive
  flow, self-prime, dose accurately, and (peristaltic) keep the wetted path to
  tubing only. This directly satisfies "a known amount must pass each pass."
- If a centrifugal is mandatory, use a small magnetic-drive unit with a
  throttling/recirculation valve and accept the efficiency penalty.

Pump heat into the water is negligible: ≈0.007 °C per pass worst-case.

## Recirculation: "pump all the water through several times"

- One pass (turn the whole tank over once) = 400/200 = **2.0 h**. `N` passes =
  `2N` hours of continuous running (e.g. 5 passes = 10 h).
- The staircase transit (38–60 s) is **0.5–0.8 % of the 2 h turnover**. The
  staircase is *not* the rate-limiting element; the tank turnover (pump flow
  ÷ inventory) sets the cycle time. To get more passes in a given time, raise
  the pump flow, not change the stairs.
- **Regime headroom is large.** The flow would have to rise to ~2,900 kg/h
  (5 cm), ~5,700 kg/h (10 cm) or ~8,600 kg/h (15 cm) — **14× to 43× the design
  flow** — before the cascade transitions to skimming. So the pump rate can be
  increased substantially for faster turnover while staying firmly in the
  nappe regime.

## Pitfalls / further analysis (on "how much water to pump")

1. **Full-width wetting at low flow.** Mass flow per unit width Γ = ṁ/W is
   1.11 (5 cm), 0.56 (10 cm) and 0.37 kg/(m·s) (15 cm). The minimum wetting
   rate for water on metal is ≈0.07–0.30 kg/(m·s). The **15 cm channel sits just
   above this band**, so at 200 kg/h it is the most likely to break into
   rivulets / leave dry stripes — giving uneven residence time and a reduced
   effective contact area. **5 cm and 10 cm wet robustly.** If 15 cm width is
   required for capacity, either raise the flow or subdivide it into narrower
   parallel lanes.
2. **Films are turbulent / wavy** (film Re 1500–4400), not smooth laminar
   sheets. Expect surface waves, some splashing and strong air entrainment in
   the nappes. Good if the device is meant to aerate; it means real transit
   time scatters ±~20 % around the model value, and the Manning closure is an
   engineering approximation (validated here to be insensitive, ±0.6 s).
3. **Evaporation over many passes.** Exposed tread film area is 0.33–0.99 m².
   Estimated loss 0.02–0.30 kg/h (still air → mild airflow). Over a 10 h /
   5-pass run that is ~0.2–3 kg (≤0.75 %); it also cools the water and
   concentrates any dissolved solids. Provide make-up water / level control for
   long runs.
4. **Centrifugal minimum-flow** instability at this duty (see pump section).
5. **Inlet free-fall** splash if the 1 m lift discharges above the top step
   (see roof section).
6. **Nappe ventilation:** seal/vent design must keep air under each nappe;
   otherwise the jet clings to the riser (Coandă) and the regime/timing shift.
7. **Long-run water quality:** recirculating warm, aerated water invites
   biofilm; stainless steel resists corrosion but the loop should be drainable
   and cleanable.

## Key assumptions / limitations

- Steady, 1-D, depth-averaged open-channel flow per tread; brink = critical
  depth control; Manning friction (n = 0.012).
- Nappe modelled as a horizontal-launch projectile (brink velocity ≈ critical
  velocity); true brink velocity is ~1.4× higher, shortening the wetted tread
  by ≈1 cm — negligible effect on transit time (hold-up dominated by the deeper
  upstream film).
- Surface tension, three-dimensional rivulet break-up, and air entrainment are
  discussed but not resolved by the 1-D model; they set the ±~20 % uncertainty
  band and the 15 cm wetting caveat. A full free-surface CFD (VOF) run would be
  the next fidelity step if sub-10 % accuracy is needed.
