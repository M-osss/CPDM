# StairCascade — water cascade over a stepped channel

Open-channel hydraulic model of water recirculated from a 400 kg tank, lifted
1 m by a pump, and cascaded down a 33-step stainless-steel staircase. Computes
flow regime, top-to-bottom transit time, liquid hold-up, required roof/wall
clearance, and pump / recirculation sizing for channel widths of 5, 10 and
15 cm.

## Files

- `stair_cascade_sim.py` — model and analysis (run this; prints the design-flow
  block, the 50 kg/min scenario, and the wrap-around report).
- `make_plots.py` — figures (`fig_profiles.png`, `fig_summary.png`,
  `fig_cascade.png`, `fig_scenarios.png`, `fig_roof_lamp.png`).
- `results.json` — machine-readable results (both flow rates).
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

## Follow-up: 50 kg/min flow, wrap-around, hold-up, lamp + roof

Run `python StairCascade/stair_cascade_sim.py` — the high-flow scenario and the
wrap-around report print after the design-flow block. Figures
`fig_scenarios.png` and `fig_roof_lamp.png` cover this section.

### 1. Running at 50 kg/min (= 3000 kg/h, 15x the design flow)

The model auto-selects the correct regime. **The regime changes** at this
flow:

| W [cm] | yc/h | regime | model | velocity | transit time | hold-up | water+spray envelope |
|---|---|---|---|---|---|---|---|
| 5  | 1.02 | **skimming** | skimming sheet | 0.81 m/s | **8 s** | 11.8 kg | 5.4 cm |
| 10 | 0.64 | high nappe | per-tread cascade | 0.43 m/s | **17 s** | 12.0 kg | 6.3 cm |
| 15 | 0.49 | nappe | per-tread cascade | 0.38 m/s | **20 s** | 14.2 kg | 5.5 cm |

What changed versus 200 kg/h:

- **5 cm crosses into skimming flow** (yc/h = 1.02 > 0.99 onset). Water no
  longer cascades step-by-step; it shoots over the step edges as a fast,
  aerated supercritical sheet with recirculating vortices filling the step
  cavities. Modelled as uniform flow on the inclined pseudo-bottom (slope
  8.5°, skimming friction factor f ≈ 0.2), not as a per-tread overfall.
- **Transit time drops sharply** (38→8 s at 5 cm, 60→20 s at 15 cm): 15x more
  flow means much higher velocity.
- **Hold-up rises ~5x** to 12–14 kg (still only ≤3.5 % of the 400 kg tank — no
  flooding).
- **Tank turnover = 8 min/pass** (was 2 h). 10 passes = 1.3 h. This flow is
  what makes "pump everything through many times" practical.
- **Pump:** duty is now 3.0 m³/h @ 1 m (≈8 W hydraulic, ~23 W shaft). At this
  flow a **small centrifugal circulator is a reasonable choice** — the earlier
  objection (centrifugal starved at 0.2 m³/h) goes away. A positive-displacement
  pump still gives steadier metering if exact per-pass volume matters.
- **Caveat — model edge:** the 10 cm "nappe" film reaches 23 mm on a 30 mm
  riser, i.e. it is close to the point where consecutive pools merge into
  transition/skimming. Treat the 10 cm number as ±25 % and verify with CFD if
  it is the chosen design. The 5 cm skimming sheet is heavily aerated
  ("white water"); its bulked depth (~34 mm) and spray drive the roof sizing.

### 2. Wrapping the staircase around a 1 m³ cube

Geometry (auto-computed): **5 steps per face** (1.0 m flight), **7 flights**,
**6 corner turns**, **1.75 loops** around the cube. Total drop 0.99 m fits the
1 m cube height; the 6.6 m straight run folds into a 1 m × 1 m footprint.

- **Straight-run hydraulics are unchanged.** Transit time, film depth and
  hold-up depend on the slope and the flow per unit width, not on the plan
  layout. Folding the run around a cube does not change the per-flight
  numbers above.
- **The 6 corners do change things.** At each 90° turn the moving water climbs
  the outer wall (superelevation `Δz = V²·W/(g·r_bend)`):

  | W [cm] | V [m/s] | Δz at r=5 cm | Δz at r=10 cm |
  |---|---|---|---|
  | 5  | 0.81 | 6.7 cm | 3.4 cm |
  | 10 | 0.43 | 3.8 cm | 1.9 cm |
  | 15 | 0.38 | 4.4 cm | 2.2 cm |

  So at sharp corners the water surface rises **3–7 cm** on the outside, with
  splashing and secondary flow. **Put a small stilling landing/pool at each of
  the 6 corners** to kill momentum and redistribute the sheet, raise the wall
  and roof locally by the superelevation, and expect minor extra head loss
  (transit time up a few percent at most). Corners are also natural **lamp
  gaps** → UV/illumination dead zones (see below).

### 3. What the hold-up is, and where it sits

**Hold-up = the water resident on the staircase at any instant** (water "in
transit"), as opposed to the ~400 kg in the tank and the small amount in the
pipes. Numerically it equals `transit_time × mass_flow`.

- **At 200 kg/h (nappe):** 2–3 kg total, located as the **3.9–7.2 mm film
  flowing across each of the 33 treads** (drawing down to critical depth at
  every brink), plus a negligible amount in the falling nappes. ~60–95 g per
  step.
- **At 50 kg/min:** 12–14 kg total.
  - Nappe widths (10, 15 cm): a thicker **1.8–2.3 cm film on each tread**.
  - Skimming width (5 cm): it **splits in two** — the fast **sheet skimming
    over the step edges (~6.9 kg)** plus **recirculating vortices trapped in
    the triangular step cavities under the pseudo-bottom (~4.9 kg)**. The
    cavity water is a new component that only exists in skimming flow.
- Even at the high flow the hold-up is ≤3.5 % of the tank, so the staircase
  cannot accumulate a meaningful fraction of the inventory and there is no
  flooding risk; the tank always holds ≥386 kg.

### 4. Overhead lamp and required roof headroom

A 3–4 cm-diameter lamp running the length of each straight flight (≈1.0 m, 5
steps → ~7 lamps total) sits in the headroom. (A 3–4 cm quartz lamp dosing
recirculated water that passes "many times" is consistent with **UV
disinfection** — dose accumulates over passes.) The lamp must stay **dry and
un-fouled**: water on a hot quartz UV lamp causes thermal-shock cracking and
scale fouling kills UV transmission.

Required interior clear height (perpendicular to the treads):

```
clear height = water+spray envelope
             + 3 cm gap (keep lamp dry/clean)
             + 4 cm lamp diameter
             + 3 cm gap (mounting + lamp cooling)
```

| condition | envelope | clear height WITHOUT lamp | clear height WITH 4 cm lamp |
|---|---|---|---|
| 200 kg/h, W 5–15 cm | 4–6 cm | ~6–9 cm | **13–15 cm** |
| 50 kg/min, W 5 cm (skimming) | 5.4 cm | 8 cm | **15.4 cm** |
| 50 kg/min, W 10 cm (worst) | 6.3 cm | 9 cm | **16.3 cm** |
| 50 kg/min, W 15 cm | 5.5 cm | 8 cm | **15.5 cm** |

**Yes — you need more headroom.** Two compounding reasons: (a) 50 kg/min makes
the water thicker and splashier than the 4–7 mm design-flow film, and (b) the
lamp plus its dry-gap and mounting eat ~10 cm of height on their own.

Recommendation:

- **≈16–18 cm interior clear height along the flights** (up from the ~10 cm
  recommended without a lamp), measured perpendicular to the treads, roof
  parallel to the 8.5° slope.
- **≈20–22 cm at the 6 corners** to cover the 3–7 cm superelevation/splash.
- Centre the lamp above the channel, keep the ≥3 cm air gap above the splash
  crest, and fit a transparent splash guard if the envelope is uncertain. The
  6 corners have no lamp → accept UV/illumination dead bands there (the corner
  stilling pools also mix the flow, which helps even out dose between passes).

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
- Skimming model (high-flow 5 cm case): uniform normal-depth flow on the
  pseudo-bottom with an equivalent Darcy friction factor f = 0.2 (literature
  range 0.17–1.0; transit time scales with √f) and a depth-averaged air
  concentration of 0.4 for bulking; the step-cavity hold-up is taken as the
  full triangular cavity volume (upper bound — real cavities are part-aerated).
  The velocity is taken as the normal-depth value over the whole length
  (ignores the short acceleration reach near the crest), so the 8 s transit is
  a slight under-estimate.
- Corner superelevation uses the standard `V²W/(g·r)` bend formula; the actual
  rise depends on the corner detail (mitre, radius, landing) and is a sizing
  guide, not a precise value.
