# Engineering Report — Water Cascade Over a Stepped Stainless-Steel Staircase

**System:** 400 kg water tank → pump (1 m lift) → 33-step cascade → return to tank
**Scope:** flow-regime analysis, transit/cycle time, hold-up, roof clearance
(incl. overhead lamp), pump selection, recirculation, and operating pitfalls.
**Flow cases analysed:** 200 kg/h (design) and 50 kg/min = 3000 kg/h (high flow).
**Channel widths:** 5, 10, 15 cm.

All numbers in this report are produced by `stair_cascade_sim.py` and are
reproducible (`python StairCascade/stair_cascade_sim.py`). Figures referenced
are in the same folder.

---

## 1. Specifications

### 1.1 Process / fluid

| Quantity | Symbol | Value |
|---|---|---|
| Tank inventory | M | 400 kg |
| Design mass flow | ṁ | 200 kg/h = 0.0556 kg/s = 3.34 L/min = 0.200 m³/h |
| High-flow case | ṁ | 50 kg/min = 3000 kg/h = 0.833 kg/s = 50.1 L/min = 3.01 m³/h |
| Pump static lift | H | 1.0 m |
| Fluid | — | water, 20 °C |
| Density | ρ | 998.2 kg/m³ |
| Dynamic viscosity | μ | 1.002 × 10⁻³ Pa·s |
| Surface tension | σ | 0.0728 N/m |

### 1.2 Staircase geometry

| Quantity | Symbol | Value |
|---|---|---|
| Number of steps | N | 33 |
| Riser height | h | 3 cm |
| Tread length | l | 20 cm |
| Channel widths | W | 5 / 10 / 15 cm |
| Material | — | stainless steel (smooth, Manning n ≈ 0.012) |
| Total vertical drop | N·h | 0.99 m |
| Total horizontal run | N·l | 6.60 m |
| Pseudo-bottom slope | θ = atan(h/l) | 8.53° |
| Inclined (chord) length | — | 6.67 m |

Net elevation check: pump lifts +1.0 m, cascade drops 0.99 m → the bottom step
discharges ≈0.01 m above the start, i.e. essentially back at tank level.

---

## 2. Methodology

The discharge per unit width is `q = Q/W` with `Q = ṁ/ρ`. Critical depth
`yc = (q²/g)^(1/3)`. The flow regime is set by `yc/h` relative to the
skimming-flow onset (Chanson 1994):

```
(yc/h)_onset = 1.057 − 0.465·(h/l) = 1.057 − 0.465·0.15 = 0.987
```

- `yc/h < 0.987` → **nappe (cascade)** flow: thin film crosses each horizontal
  tread, springs off the brink as a free nappe, falls h to the next tread.
- `yc/h ≥ 0.987` → **skimming** flow: a fast aerated sheet skims the step edges
  over the pseudo-bottom while recirculating vortices fill the step cavities.

**Nappe model (per tread).** Horizontal-bed gradually-varied flow with a free
overfall (critical-depth control at the brink):

```
dy/dx = −Sf/(1 − Fr²),   Sf = n²(q/y)²/R^(4/3),   Fr² = q²/(g y³),   R = Wy/(W+2y)
```

integrated upstream from the brink over the wetted tread. Transit time across a
tread = residence time = ∫(dx/V) = (∫y dx)/q. The nappe is launched horizontally
at the brink velocity and falls h (projectile).

**Skimming model.** Uniform normal-depth flow on the inclined pseudo-bottom from
Darcy–Weisbach with a skimming friction factor f ≈ 0.2:

```
Sf = f·V²/(8 g R) = sin θ ,  V = q/d   →  solve for normal depth d
```

Velocity V = q/d; transit time = (inclined length)/V; depth bulked by a
depth-averaged air concentration of 0.4 for freeboard.

**Cycle time** (1st stair → bottom → tank) = staircase transit + free-fall of
the bottom nappe into the tank (`t = √(2·0.10/g) ≈ 0.14 s`, taking a nominal
0.10 m exit drop).

---

## 3. Sample calculations

### 3.1 Worked example A — W = 5 cm, 200 kg/h (nappe)

```
Q   = ṁ/ρ              = 0.0556 / 998.2          = 5.566e-5 m³/s
q   = Q/W              = 5.566e-5 / 0.05          = 1.113e-3 m²/s
yc  = (q²/g)^(1/3)     = (1.113e-3²/9.81)^(1/3)   = 5.02e-3 m   (5.02 mm)
yc/h= 5.02/30 = 0.167  < 0.987  → NAPPE flow                    ✓
Vc  = q/yc            = 1.113e-3 / 5.02e-3        = 0.222 m/s
t_fall = √(2h/g)      = √(2·0.03/9.81)            = 0.0782 s
x_land = Vc·t_fall    = 0.222·0.0782              = 0.0174 m  (lands 1.7 cm in)
wetted tread = l − x_land = 20 − 1.7              = 18.3 cm

GVF profile (numeric): depth 7.20 mm at landing → 5.02 mm (critical) at brink.
Hold-up per tread (mass) = 1.962 kg / 33 steps    = 0.0595 kg
t_tread = mass_step/ṁ = 0.0595/0.0556             = 1.070 s    (= residence)
Transit T = N·(t_tread + t_fall) = 33·(1.070 + 0.0782) = 37.9 s
Cycle = T + exit fall = 37.9 + 0.14               = 38.0 s

Splash: v_impact = √(Vc² + 2gh) = √(0.222² + 0.589) = 0.799 m/s
        splash height = v_impact²/(2g) = 0.638/19.62 = 0.0325 m (3.25 cm)
Freeboard (min) = (y_max + splash)·1.5 = (7.20 + 32.5 mm)·1.5 = 5.96 cm
Total hold-up = 1.96 kg (0.49 % of tank)
```

### 3.2 Worked example B — W = 5 cm, 50 kg/min (skimming)

```
Q   = 0.833/998.2 = 8.348e-4 m³/s ;  q = Q/W = 0.01670 m²/s
yc  = (q²/g)^(1/3) = 0.0305 m (30.5 mm) ;  yc/h = 1.017 > 0.987 → SKIMMING ✓
θ   = atan(3/20) = 8.53° ; sinθ = 0.1483 ; chord = 0.2022 m ; L_incl = 6.674 m

Normal depth (Darcy–Weisbach, f = 0.2):  f·V²/(8 g R) = sinθ , V = q/d
  → d (normal) = 20.6 mm ,  V = q/d = 0.01670/0.0206 = 0.811 m/s
  check: R = Wd/(W+2d) = 0.01129 m ; f·V²/(8gR) = 0.1483 = sinθ ✓
  Fr = V/√(g d cosθ) = 0.811/√(0.1998) = 1.81  (supercritical) ✓
  bulked depth = d/(1−0.4) = 34.3 mm

Transit T = L_incl/V = 6.674/0.811 = 8.23 s ;  Cycle = 8.23 + 0.14 = 8.4 s

Hold-up:
  mainstream sheet = d·L_incl·W·ρ = 0.0206·6.674·0.05·998.2 = 6.86 kg
  cavity vortices  = ½·h·l·W·N·ρ  = 0.5·0.03·0.20·0.05·33·998.2 = 4.94 kg
  TOTAL            = 11.80 kg  (3.0 % of tank)
Water+spray envelope = bulked + spray = 34.3 + 20 mm = 5.43 cm
```

### 3.3 Pump sample

```
Hydraulic power  P_hyd = ρ g Q H
  200 kg/h : 998.2·9.81·5.566e-5·1.0 = 0.55 W ;  shaft ≈ 0.55/0.35 = 1.6 W
  50 kg/min: 998.2·9.81·8.348e-4·1.0 = 8.18 W ;  shaft ≈ 8.18/0.35 = 23 W
Tank turnover = M/ṁ :  400/200 = 2.0 h ;  400/3000 = 0.133 h = 8.0 min
```

### 3.4 Corner superelevation sample (wrap-around, §8)

```
Δz = V²·W/(g·r_bend)  ;  W = 5 cm, V = 0.81 m/s, r_bend = 5 cm
   = 0.811²·0.05/(9.81·0.05) = 6.7 cm  (outer-wall rise at a 90° corner)
```

### 3.5 Calculation map for every reported column

This section explains how each number in the results tables is produced. It is
the same calculation chain used in `stair_cascade_sim.py`.

#### A. Unit conversions

All calculations are SI internally:

```
ṁ [kg/s] = ṁ [kg/h] / 3600
Q [m³/s] = ṁ / ρ
W [m]    = W [cm] / 100
h [m]    = 0.03
l [m]    = 0.20
```

Example:

```
200 kg/h  → ṁ = 200/3600 = 0.0556 kg/s
Q         = 0.0556/998.2 = 5.566×10⁻⁵ m³/s = 3.34 L/min

50 kg/min → ṁ = 50/60 = 0.833 kg/s
Q         = 0.833/998.2 = 8.348×10⁻⁴ m³/s = 50.1 L/min
```

#### B. Flow per unit width `q`

The same pump flow spread over a wider channel gives a smaller unit-width flow:

```
q = Q/W
```

At 200 kg/h:

```
W = 5 cm  → q = 5.566e-5/0.05 = 1.113e-3 m²/s
W = 10 cm → q = 5.566e-5/0.10 = 5.566e-4 m²/s
W = 15 cm → q = 5.566e-5/0.15 = 3.710e-4 m²/s
```

This is why the narrow channel runs faster: the same water is squeezed into a
smaller width, so the depth and velocity change.

#### C. Critical depth and regime

For a rectangular channel, critical depth is

```
yc = (q²/g)^(1/3)
```

Then compare `yc/h` with the stepped-channel skimming onset:

```
(yc/h)_onset = 1.057 − 0.465(h/l) = 0.987
```

Interpretation:

- `yc/h < 0.987` → nappe/cascade.
- `yc/h ≥ 0.987` → skimming.

At 50 kg/min, W = 5 cm:

```
q  = 8.348e-4/0.05 = 0.01670 m²/s
yc = (0.01670²/9.81)^(1/3) = 0.0305 m
yc/h = 0.0305/0.030 = 1.02 > 0.987 → skimming
```

At 50 kg/min, W = 10 cm:

```
q  = 8.348e-4/0.10 = 0.00835 m²/s
yc = 0.0192 m
yc/h = 0.64 < 0.987 → still nappe/high-nappe
```

#### D. Nappe-flow velocity, landing distance and wetted tread

In nappe flow the brink acts as the control section, so the brink depth is
critical and the brink velocity is

```
Vc = q/yc
```

The falling nappe drops one riser (`h = 0.03 m`):

```
t_fall = √(2h/g) = 0.0782 s
x_land = Vc·t_fall
wetted tread = l − x_land
```

At 200 kg/h, W = 5 cm:

```
Vc = 0.001113/0.00502 = 0.222 m/s
x_land = 0.222·0.0782 = 0.0174 m = 1.7 cm
wetted tread = 20 − 1.7 = 18.3 cm
```

This confirms the jet lands on the next tread and does not overshoot it.

#### E. Nappe-flow depth profile and hold-up

The water depth on each tread is not assumed constant. It is computed by
integrating the gradually-varied-flow equation:

```
dy/dx = −Sf/(1 − Fr²)
Sf  = n²(q/y)²/R^(4/3)
Fr² = q²/(g y³)
R   = Wy/(W+2y)
```

The boundary condition is `y = yc` at the brink. The model integrates upstream
from the brink to the landing point. That gives a depth profile `y(x)`.

The water volume on one tread is the cross-sectional water area integrated over
the wetted length:

```
volume_step = W ∫ y dx
mass_step   = ρ·volume_step
```

The total staircase hold-up in nappe flow is

```
hold-up = N·mass_step
```

The residence time on a tread is mass on that tread divided by mass flow:

```
t_tread = mass_step/ṁ = (∫y dx)/q
```

Total staircase descent time:

```
T_descent = N·(t_tread + t_fall)
```

The film-residence component is equivalent to

```
T_film = tread-film hold-up / ṁ
```

Example check, 200 kg/h, W = 5 cm:

```
hold-up = 1.962 kg
ṁ       = 0.0556 kg/s
film hold-up/ṁ = 1.962/0.0556 = 35.3 s
```

The difference between 35.3 s and the reported 37.9 s is the explicit free-fall
time in the nappes (`33·0.0782 = 2.58 s`). The film hold-up is on the treads;
the falling-jet hold-up is small but still contributes to parcel travel time.

#### F. Skimming-flow normal depth, velocity and hold-up

For skimming flow the water no longer stops and re-forms on every tread. It
flows over a pseudo-bottom at angle `θ = atan(h/l) = 8.53°`. The model solves
normal depth `d` from

```
Sf = sinθ = f V²/(8 g R)
V = q/d
R = Wd/(W+2d)
```

This equation is solved numerically because `R` contains `d`.

After `d` is known:

```
V = q/d
T_descent = L_incl/V
L_incl = N·√(h²+l²)
```

The total skimming hold-up has two parts:

```
sheet hold-up  = ρ·d·L_incl·W
cavity hold-up = ρ·(0.5·h·l·W)·N
total hold-up  = sheet + cavity
```

The cavity term is the water trapped in the triangular recirculation zones
inside each step. This term is absent in nappe flow.

#### G. Cycle time from first stair back to tank

The report uses:

```
cycle time = staircase descent time + bottom exit fall time
exit fall  = √(2·0.10/g) = 0.143 s
```

The exit fall is small. It is included so the stated cycle explicitly ends at
the bottom tank, not at the last stair edge.

The pipe/pump-up travel time is not included in the headline cycle because it
depends on the actual pipe length and diameter. For the recommended pipe sizes,
it is typically ~1–2 s and does not control the result. Tank turnover controls
the time between repeated passes of the same parcel.

#### H. Roof/freeboard and lamp clearance

For no lamp, the report uses:

```
minimum clear height = 1.5·(water depth + splash allowance)
```

For a 4 cm lamp above the water:

```
clear height = water+spray envelope
             + 3 cm dry gap below lamp
             + 4 cm lamp diameter
             + 3 cm mounting/cooling gap above lamp
```

Example, high-flow W = 10 cm:

```
water+spray envelope = 6.3 cm
clear height = 6.3 + 3 + 4 + 3 = 16.3 cm
```

Corners get extra headroom because the water surface rises on the outer wall:

```
corner rise Δz = V²W/(g r_bend)
```

The recommendation of 20–22 cm at corners is the flight clearance plus the
3–7 cm bend rise, rounded upward for splash.

---

## 4. Results

### 4.1 Design flow — 200 kg/h (all widths nappe)

| W [cm] | regime | yc [mm] | yc/h | Vc [m/s] | **transit** | **cycle** | hold-up [kg] | max film [mm] | film Re |
|---|---|---|---|---|---|---|---|---|---|
| 5  | nappe | 5.02 | 0.167 | 0.222 | 37.9 s | **38.0 s** | 1.96 | 7.20 | 4436 |
| 10 | nappe | 3.16 | 0.105 | 0.176 | 50.2 s | **50.3 s** | 2.65 | 4.81 | 2218 |
| 15 | nappe | 2.41 | 0.080 | 0.154 | 59.8 s | **59.9 s** | 3.18 | 3.85 | 1479 |

Transit time insensitive to roughness: at W = 5 cm, n = 0.011–0.013 gives
37.3–38.5 s (±0.6 s).

### 4.2 High flow — 50 kg/min = 3000 kg/h

| W [cm] | regime | yc/h | V [m/s] | **transit** | **cycle** | hold-up [kg] | envelope [cm] |
|---|---|---|---|---|---|---|---|
| 5  | **skimming** | 1.02 | 0.81 | 8.2 s | **8.4 s** | 11.8 (sheet 6.9 + cavity 4.9) | 5.4 |
| 10 | high nappe | 0.64 | 0.43 | 16.9 s | **17.1 s** | 12.0 | 6.3 |
| 15 | nappe | 0.49 | 0.38 | 19.6 s | **19.8 s** | 14.2 | 5.5 |

15× more flow → transit drops ~5×, hold-up rises ~5×, and the 5 cm channel
crosses into skimming flow. See `fig_scenarios.png`.

### 4.3 Cycle time of one water parcel (headline)

"One cycle" = top (1st) stair → down all 33 steps → fall into the bottom tank.
It equals the staircase descent plus a ~0.14 s exit fall:

| W [cm] | cycle @ 200 kg/h | cycle @ 50 kg/min |
|---|---|---|
| 5  | **38.0 s** | **8.4 s** |
| 10 | **50.3 s** | **17.1 s** |
| 15 | **59.9 s** | **19.8 s** |

The pump-up leg (tank → +1 m → top stair) adds only ~1–2 s of pipe transport,
so the **gravity descent dominates the per-parcel cycle**. The interval between
successive passes of a given parcel is governed instead by tank turnover
(§6): 2 h/pass at design flow, 8 min/pass at high flow.

---

## 5. Hold-up — definition and location

Reported hold-up = liquid physically supported by the staircase at any instant:
the tread films in nappe flow, or the skimming sheet plus step-cavity vortices
in skimming flow. It is **not** the tank inventory.

For nappe flow there is also a small quantity of water in the falling jets
between steps:

```
jet hold-up = ṁ·N·t_fall
```

At 200 kg/h this is `0.0556·33·0.0782 = 0.14 kg`; at 50 kg/min it is
`0.833·33·0.0782 = 2.15 kg`. The reported table hold-up excludes this airborne
jet mass because it is not retained on the stairs; cycle time includes it
through the explicit free-fall term.

- **200 kg/h (nappe):** 2–3 kg on the treads as the 3.9–7.2 mm film flowing
  across each of the 33 treads (drawing to critical depth at every brink),
  ~60–95 g/step, plus ~0.14 kg in falling jets.
- **50 kg/min:** 12–14 kg on the staircase. Nappe widths → a thicker
  1.8–2.3 cm tread film plus ~2.15 kg in falling jets. Skimming (5 cm) → splits
  into the fast sheet over the step edges (~6.9 kg) plus recirculating vortices
  in the triangular step cavities (~4.9 kg), with no discrete per-step falling
  jets.

Hold-up is ≤3.5 % of the tank even at high flow, so the staircase cannot
accumulate a meaningful fraction of the inventory; the tank always holds
≥386 kg. There is no flooding risk; the roof is sized for splash, not bulk
water (§7).

---

## 6. Pump selection and recirculation

### 6.1 Duty and recommendation

| | 200 kg/h | 50 kg/min |
|---|---|---|
| Flow | 0.20 m³/h (3.3 L/min) | 3.0 m³/h (50 L/min) |
| Static head | 1.0 m | 1.0 m |
| Add pipe/fitting/entrance losses | +~0.5–1 m | +~0.5–1.5 m |
| **Total dynamic head (design to)** | ~1.5–2 m | ~1.5–2.5 m |
| Hydraulic power | 0.55 W | 8.2 W |
| Shaft power (η ≈ 0.35) | ~1.6 W | ~23 W |
| Tank turnover | 2.0 h/pass | 8.0 min/pass |

**At 200 kg/h — use a positive-displacement (peristaltic or small gear) pump.**
A catalogue centrifugal pump at 0.2 m³/h runs near shut-off, where it suffers
internal recirculation, vibration and self-heating, and the delivered flow is
hypersensitive to head — unacceptable for a process that needs steady metered
flow. A PD pump gives near-constant, head-insensitive, accurately metered flow
and self-primes. Size the motor ~10–20 W (efficiency + margin).

**At 50 kg/min — a small centrifugal circulator is acceptable** (3 m³/h at low
head is a normal small-pump duty); the earlier low-flow objection disappears.
A PD pump still gives the steadiest per-pass volume if exact dosing matters.
Size the motor ~40–80 W.

**Common requirements:** 316 stainless steel or food-grade polymer wetted parts
(corrosion + hygiene; mandatory if this is a UV water-treatment loop — see §7);
variable speed (VFD for centrifugal, speed control for PD) to trim flow and set
the number of passes; a strainer at the tank suction; and, if a centrifugal is
used at the low-flow case, a min-flow recirculation/bypass line. NPSH is a
non-issue at 1 m lift. Pump heat into the water is negligible (≤0.007 °C/pass).

### 6.2 "Pump all the water through several times"

One pass turns the whole tank over: 2.0 h (design) or 8.0 min (high flow). N
passes ≈ N × turnover (e.g. 10 passes = 20 h design / 1.3 h high flow). The
staircase descent (8–60 s) is <1 % of a pass, so **the tank+pump set the cycle
rate, not the stairs**. To get more passes in a given time, raise the pump flow.

**Regime headroom (design flow):** skimming would not begin until ~2,870 kg/h
(5 cm), 5,740 kg/h (10 cm) or 8,610 kg/h (15 cm) — 14×–43× the design flow — so
there is large room to increase flow while staying in nappe regime (until the
50 kg/min case, where 5 cm does reach skimming).

### 6.3 Piping calculations and sizing

Pipe velocity is calculated from

```
A = πD²/4
Vpipe = Q/A
Re = ρVpipeD/μ
```

Straight-pipe friction head is estimated with Darcy-Weisbach:

```
hf/L = f·Vpipe²/(2gD)
```

using `f = 64/Re` if laminar and the smooth-turbulent Blasius estimate
`f = 0.3164/Re^0.25` otherwise. Fittings, elbows, valves, strainers, distributors
and the top inlet add extra minor losses, so pump TDH is specified higher than
straight-pipe loss alone.

#### 200 kg/h piping

`Q = 5.566×10⁻⁵ m³/s = 3.34 L/min`

| Pipe ID | Velocity | Re | straight head loss |
|---|---:|---:|---:|
| 6 mm | 1.97 m/s | 11,800 | 1.00 m/m |
| 8 mm | 1.11 m/s | 8,800 | 0.25 m/m |
| 10 mm | 0.71 m/s | 7,100 | 0.09 m/m |
| 12 mm | 0.49 m/s | 5,900 | 0.04 m/m |

Recommended: **8–10 mm ID** for a peristaltic/gear pump. 6 mm works but wastes
head; 12 mm is fine if the run is long but is bulkier and can trap air.

#### 50 kg/min piping

`Q = 8.348×10⁻⁴ m³/s = 50.1 L/min`

| Pipe ID | Velocity | Re | straight head loss |
|---|---:|---:|---:|
| 15 mm | 4.72 m/s | 70,600 | 1.47 m/m |
| 20 mm | 2.66 m/s | 52,900 | 0.38 m/m |
| 25 mm | 1.70 m/s | 42,400 | 0.13 m/m |
| 32 mm | 1.04 m/s | 33,100 | 0.04 m/m |
| 40 mm | 0.66 m/s | 26,500 | 0.01 m/m |

Recommended: **DN32 preferred** (DN25 minimum) for the high-flow case. DN15 and
DN20 are too fast/noisy and add avoidable head. Use a flooded suction, short
suction line, strainer before the pump, vents at high points, drains at low
points, and a stilling box/distribution manifold at the top stair so the inlet
arrives as a uniform sheet instead of a jet.

---

## 7. Roof / wall clearance and the overhead lamp

Two concerns: (a) contain splash, and (b) fit a 3–4 cm-diameter lamp running the
length of each straight flight (≈1.0 m, 5 steps → ~7 lamps). A 3–4 cm quartz
lamp dosing recirculated water that passes "many times" is consistent with **UV
disinfection** (dose accumulates over passes). The lamp must stay **dry and
un-fouled** — water on a hot quartz UV lamp causes thermal-shock cracking, and
scale fouling destroys UV transmission.

Required interior clear height (perpendicular to the treads):

```
clear height = water+spray envelope + 3 cm (dry gap) + 4 cm (lamp) + 3 cm (mount/cooling)
```

| Condition | envelope | WITHOUT lamp | WITH 4 cm lamp |
|---|---|---|---|
| 200 kg/h (W 5–15 cm) | 4–6 cm | ~6–9 cm | **13–15 cm** |
| 50 kg/min, W 5 cm (skimming) | 5.4 cm | 8 cm | **15.4 cm** |
| 50 kg/min, W 10 cm (worst) | 6.3 cm | 9 cm | **16.3 cm** |
| 50 kg/min, W 15 cm | 5.5 cm | 8 cm | **15.5 cm** |

**Recommendation:** ≈16–18 cm interior clear height along the flights (up from
~10 cm without a lamp), roof parallel to the 8.53° slope; ≈20–22 cm at the 6
corners to cover superelevation/splash (§8). Centre the lamp above the channel,
keep ≥3 cm air gap above the splash crest, and fit a transparent splash guard if
the envelope is uncertain. See cross-section `fig_roof_lamp.png`.

---

## 8. Wrapping the staircase around a 1 m³ cube

| Quantity | Value |
|---|---|
| Steps per face | 5 (1.0 m flight) |
| Flights | 7 |
| Corner (90°) turns | 6 |
| Loops around cube | 1.75 |
| Total drop | 0.99 m (fits the 1 m cube height) |
| Footprint | 6.6 m run folded into 1 m × 1 m |

- **Straight-run hydraulics are unchanged** — transit, depth and hold-up depend
  on slope and flow per unit width, not on the plan layout.
- **The 6 corners change behaviour locally.** Outer-wall superelevation
  `Δz = V²W/(g·r)`:

  | W [cm] | V [m/s] | Δz @ r = 5 cm | Δz @ r = 10 cm |
  |---|---|---|---|
  | 5  | 0.81 | 6.7 cm | 3.4 cm |
  | 10 | 0.43 | 3.8 cm | 1.9 cm |
  | 15 | 0.38 | 4.4 cm | 2.2 cm |

  → put a small **stilling landing/pool at each corner** to kill momentum and
  redistribute the sheet, raise the wall/roof locally, and accept a few percent
  extra transit and **lamp/UV dead-bands** at the corners.

---

## 9. Pitfalls and operating envelope

1. **Full-width wetting at low flow.** Γ = ṁ/W = 1.11 / 0.56 / 0.37 kg/(m·s) at
   design flow; minimum wetting for water on metal ≈ 0.07–0.30 kg/(m·s). The
   **15 cm channel sits just above this band** → most prone to rivulets / dry
   stripes → uneven residence and reduced contact area. 5 and 10 cm wet
   robustly. At 50 kg/min all widths are fully wetted (Γ ≥ 5.6).
2. **Turbulent / wavy films** (film Re 1.5k–67k): surface waves, splashing,
   strong aeration. Good for aeration/UV mixing; it means real transit scatters
   ~±20 % around the model and validates the turbulent (Manning / Darcy) closure.
3. **Regime shift with flow.** 5 cm goes skimming at ~2,870 kg/h; above that,
   step-by-step contact (residence per pass) drops sharply. Choose flow/width to
   stay in the intended regime.
4. **Centrifugal min-flow** instability at the low-flow case (§6).
5. **Evaporation:** exposed film area 0.33–0.99 m²; loss ~0.02–0.30 kg/h (still
   → mild airflow). Over long multi-pass runs this removes water, cools it and
   concentrates solutes — provide make-up/level control.
6. **Inlet free-fall:** delivering the 1 m lift directly above the top step
   impacts at 4.4 m/s and splashes ~1 m. Deliver tangentially at top-step level
   or via a stilling box; only the head needs extra height.
7. **Nappe ventilation:** keep air under each nappe or the jet clings to the
   riser (Coandă) and the regime/timing shift.
8. **Long-run water quality / lamp fouling:** recirculating warm aerated water
   invites biofilm; keep the loop drainable and the lamp sleeves cleanable.

**Width recommendation.** At 200 kg/h, 5–10 cm for robust wetting (5 cm fastest,
least hold-up). At 50 kg/min, **10 cm** keeps step-by-step nappe contact and
full wetting (best for residence/UV dose per pass); 5 cm maximises throughput
but goes skimming (less contact per pass, but more passes per hour).

---

## 10. Assumptions and limitations

- Steady, 1-D depth-averaged open-channel flow; brink = critical-depth control
  (nappe); Manning n = 0.012.
- Nappe launched at brink ≈ critical velocity (true brink velocity ~1.4×;
  shortens wetted tread ~1 cm — negligible on transit).
- Skimming model: uniform normal depth, equivalent Darcy f = 0.2 (range
  0.17–1.0; transit ∝ √f), air concentration 0.4 for bulking; step-cavity
  hold-up taken as the full triangular volume (upper bound). Uniform-velocity
  transit ignores the short crest acceleration → the 8 s figure is a slight
  under-estimate.
- The 10 cm / 50 kg/min "nappe" film (23 mm on a 30 mm riser) is close to pool
  merging → treat as transition, ±25 %.
- Corner superelevation from the standard `V²W/(g·r)` bend formula — a sizing
  guide, not an exact value (depends on corner detail).
- 3-D rivulets, surface tension and air entrainment set the ±~20 % band. A
  free-surface CFD (VOF) run is the next fidelity step if sub-10 % accuracy is
  required.

---

## 11. Nomenclature

| Symbol | Meaning | Unit |
|---|---|---|
| ṁ, Q, q | mass / volumetric / unit-width flow | kg/s, m³/s, m²/s |
| h, l, W, N | riser, tread, width, step count | m, m, m, – |
| yc, y, d | critical depth, tread depth, normal depth | m |
| Vc, V | critical / mean velocity | m/s |
| Fr, Re | Froude, Reynolds number | – |
| Sf, S0, θ | friction slope, bed slope, pseudo-bottom angle | –, –, ° |
| n, f | Manning roughness, Darcy friction factor | –, – |
| Γ | wetting rate ṁ/W | kg/(m·s) |
| T, Δz | transit time, bend superelevation | s, m |

## 12. References

1. Chanson, H. (1994/2002). *Hydraulics of Stepped Chutes and Spillways.*
   Nappe/skimming onset `(yc/h) = 1.057 − 0.465 h/l`; skimming friction factor.
2. Chow, V.T. (1959). *Open-Channel Hydraulics.* Gradually-varied flow, free
   overfall, critical-depth control.
3. Hartley, D.E. & Murgatroyd, W. (1964). Minimum-wetting-rate criteria for thin
   liquid films.

---

## Appendix — Reproduce

```
pip install -r StairCascade/requirements.txt
python StairCascade/stair_cascade_sim.py     # prints both flow cases + wrap report; writes results.json
python StairCascade/make_plots.py            # writes the 5 figures
```

Figures: `fig_profiles.png` (tread depth profile), `fig_summary.png`
(transit/hold-up vs width), `fig_cascade.png` (cascade side view),
`fig_scenarios.png` (transit vs width, both flows), `fig_roof_lamp.png` (roof +
lamp cross-section).
