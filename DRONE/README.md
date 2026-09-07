# Speed5 — 5" 6S digital speed quad

Ultra-fast quad in the 3–5" class. Video and control are both digital. Analog 5.8 GHz video is not used.

| | |
|---|---|
| Frame | AOS 5 EVO (5" aero truss, deadcat camera path) |
| Motors | iFlight XING2 2207 **2550 KV** |
| Props | Gemfan Hurricane **5.1 × 5.26 × 2** (high pitch, 2-blade) |
| Battery | Tattu R-Line V5 **6S 1050 mAh 150C** LiHV, XT60 |
| Video | **HDZero 1W** + Nano V3 16:9 |
| Control | ELRS **Gemini Xrossband** (2.4 GHz + 900 MHz) |
| AUW | ~543 g with pack |
| Modelled top speed | **290 km/h** (`python3 speed_model.py`) |

The catalogue in `speed_model.py` is calibrated against the published iFlight Mach R5 Ultra 6S figure (2207 2100 KV, ~245 km/h class). Ranking, not a GPS guarantee.

Run:

```
python3 speed_model.py
python3 -m unittest test_speed_model.py
```

---

## Size: 5", not 3", not 5.5"

Thrust scales roughly with D⁴, prop power with D⁵. A 3" disc is 34% of a 5.1" disc area. Absolute top speed in this class is a 5" problem.

| Config | Modelled top | Why it loses |
|---|---|---|
| 3" 6S 4200 KV | 237 km/h | Reynolds + stator thermal + electronics that do not shrink with the frame |
| 3.5" 6S 3400 KV | 259 km/h | Same, less severe |
| **5" 6S 2550 KV 2-blade** | **290 km/h** | Disc area + pitch speed + 2207 copper |
| 5.5" 6S 1800 KV | 216 km/h | Correct KV for 5.5" is low; fights "ultra high KV" |
| 5" 4S 3800 KV | 281 km/h | Close on paper, loses on current (below) |

Going **bigger than 5"** (5.5 / 6 / 7) raises top speed only if KV is dropped so tip Mach stays sane. That is a long-range / cine airframe, not this brief. Going **smaller than 5"** raises RPM and pitch speed on paper and then dumps it as heat and transonic tips.

Stay at **5.1" props on a 5" wheelbase**.

---

## 6S, not 4S

RPM ≈ KV × V. Matching a 6S 2550 KV motor on 4S needs ~3800 KV. Power at the pack is V × I. Same shaft power on 4S is ~1.5× the current.

From the model, 5" high-pitch 2-blade:

| | 6S 2550 KV LiHV | 4S 3800 KV LiHV |
|---|---|---|
| Loaded RPM | 48.3k | 45.6k |
| Modelled top | 290 km/h | 281 km/h |
| Pack current (punch) | 61 A | 74 A |
| I²R proxy (4S / 6S) | 1.00 | **1.49** |

6S wins three times:

1. **Current.** ESC FETs, XT60, pigtail, and cell IR all heat as I²R. 4S burns ~50% more lead/FET loss at the same power. That loss is voltage sag, which is lost RPM, which is lost speed.
2. **Pack.** A 6S 1050 150C is the standard high-C race cell. 4S 1300–1500 150C exists, but the extra copper in the 3800 KV winding and the 4-in-1 is the weaker link.
3. **Headroom.** 6S 2550 KV is already at tip Mach ~0.96. 4S 3800 KV is not slower because of voltage; it is slower because sag and heat pull RPM down.

4S is the right answer for a 3–4" 2004 build, or a 5" 2700–3100 KV freestyle quad that must stay light on a budget. It is not the speed answer.

**8S** (1500 KV 2207) is the actual next step for record speed. Rejected here: the brief asked 4S vs 6S, HDZero 1W and most 4-in-1 BECs are happier on 6S, and 8S 1050 packs are scarce. Revisit 8S only after 6S motors are thermally maxed.

Charge **4.35 V/cell for speed runs**. 6 × 4.35 = 26.1 V full. That is free RPM. Practice at 4.20 V.

---

## Motors and props

**XING2 2207 2550 KV** is past the usual 6S 5" race window (1750–2100 KV). Mach R5 Ultra uses 2100 KV and already claims 245 km/h. 2550 KV exists to push pitch speed after the prop unloads in a straight.

No-load at 26.1 V: 2550 × 26.1 = 66.6k RPM. The model assumes 78% of that at speed (48.3k) because the 5.26" pitch still loads the motor. Tip Mach 0.96. Above ~0.92 the code applies a compressibility penalty; spinning yet faster is how you cook magnets, not how you gain km/h.

**2-blade 5.1×5.26** is the high-pitch choice. Two blades: less profile drag at Vmax, less current, less lock-in. Three blades (51466-3) if the 2-blade yaws in gusts or the camera is unreadable.

If winding temp > 90 °C after a 20 s run, drop to **2150 KV** or to the 4.8" pitch 2-blade (`speed5_6s_2100kv_safer`, 260 km/h modelled, 42 A). That is the thermal fallback, not the primary.

---

## Signal: HDZero, not DJI, not analog

At 70 m/s, 25 ms of extra video delay is 1.8 m of ground. The pilot is late on every correction.

| System | Glass-to-glass | Failure mode | Mass | Verdict |
|---|---|---|---|---|
| Analog 5.8 | ~9–15 ms | snow | low | forbidden by brief |
| **HDZero 1W** | **~14–20 ms, fixed** | sparkles, then break | ~20 g VTX | **chosen** |
| Walksnail Avatar | ~22–30 ms, variable | smear / freeze | mid | wrong for Vmax |
| DJI O4 (race mode) | ~15–28 ms, variable | freeze / retransmit | higher | better picture, worse when the link is ugly |
| OpenIPC | ~35–70 ms | software | DIY | latency disqualifies |

HDZero is the unconventional pick versus the 2026 default (DJI O4). It is the correct pick for speed: one-way, line-by-line, latency does not balloon as RSSI falls. Image breaks like analog (sparkles) instead of holding a stale frame while the quad covers another 10 m.

**1 W, not the Race 25–400 mW VTX.** A 1 km straight at 250 km/h is 14 s. The aircraft is hundreds of metres out, often tail-on, carbon edge-on to the pilot. 25 mW is a park flyer number.

Camera: **Nano V3 16:9**, 15–25° uptilt. No GoPro. DVR in the VTX/goggles. Action-cam mass is speed given back.

Goggles: HDZero Goggles / G2. There is no adapter that makes this feed DJI Goggles 3.

---

## ELRS: Gemini Xrossband, not a single 2.4 ceramic

Default race build: 2.4 GHz ELRS, one T-antenna, 500 Hz. That is fine until the airframe is a carbon knife flying away at 70 m/s and the 2.4 Fresnel zone is the aircraft itself.

**RadioMaster XR4** (dual LR1121) in **Gemini Xrossband**:

- Antenna 1: 900 MHz (penetration, diffraction around carbon and terrain)
- Antenna 2: 2.4 GHz (low latency, 500–1000 Hz)

Both packets, every frame. LQ holds when one band is in a null.

Packet rate:

- Close-in speed run: 2.4 Gemini **D500** or **F1000**
- Long straight / receding tail: **X100 / X150 Full Res** (true dual-band)

TX: Nomad or Bandit on ELRS 3.5+. Single-band 2.4 receivers will not sync in X-modes.

Antenna placement is the unconventional part. Pagodas and mushrooms on the canopy are drag and they sit in motor noise. **Dual-band T-antennas in TPU along the rear arms**, 90° polarization, ≥80 mm apart, 15 mm off carbon. 5.8 HDZero antenna on the tail, 45°, ≥40 mm from the 2.4 element.

---

## Frame: AOS 5 EVO, not a 70 g race brick

Drag is ½ ρ CdA V². At 70 m/s that term dominates mass. 40 g extra frame is irrelevant; arm section is not.

AOS 5 EVO (Chris Rosser): topology-optimised **truss arms**, 5" props, **deadcat** so the digital camera never sees a prop disc. Digital compression (and HDZero sparkle) treats a prop in frame as noise. Clean glass is a speed aid because the pilot can keep the nose down.

Why not:

- **Mach R5 190 mm brick** — lighter, proven 245 km/h analog. Worse camera path, bluff arms. Use it only if the AOS cannot be sourced.
- **AOS 5.5 EVO** — 235 mm, 5.5" props, wants 1800 KV 2208. Out of the 3–5" bound and the high-KV bound.
- **True-X 5" freestyle** — props in the digital image, no arm aero.
- **3" whoop duct** — ducts are speed brakes.

Hardware: 2.5 mm plates as designed, steel motor screws, threadlock. No aluminium camera side plates next to the 5.8 SMA.

---

## Weight budget

| Item | g |
|---|---|
| AOS 5 EVO + hardware | 125 |
| XING2 2207 × 4 | 128 |
| F7 + 60–80 A 4-in-1 | 22 |
| HDZero 1W + Nano V3 + 5.8 antenna | 28 |
| XR4 + dual-band Ts | 10 |
| M10 GPS | 9 |
| Cap, 12 AWG pigtail, TPU, strap | 35 |
| **Dry** | **357** |
| Tattu 6S 1050 150C | 186 |
| **AUW** | **543** |

GPS is on purpose. The only number that matters after the build is a GPS max-speed log, not a bench KV argument.

---

## Electrical

See `wiring.md`. Short 12 AWG XT60. 1000 µF 50 V on the ESC. HDZero on the 8–16 V pad, never the 5 V rail. Bidirectional DShot300, AM32, RPM filter. Starting CLI: `betaflight/speed5_cli.txt`.

Pack C: 1050 mAh × 150C = 157 A paper. Modelled punch is 61 A average pack; phase current is higher. If punch sag < 21.5 V, the pack is the limiter — step to R-Line 6S 1300 150C and accept ~20 g.

---

## How to fly it

Speed runs are a straight, a clock, and a battery. 20–40 s full send. Land. Check motor temp with a finger / IR. If the magnets smell, the KV is too high for that pitch.

Betaflight GPS Rescue is the recovery, not a cruise mode. Do not fly this over people, roads, or unlicensed spectrum. AUW 543 g is over the 250 g remote-id threshold in several jurisdictions; that is the operator's problem to satisfy.

Expected real GPS, level, no dive: **230–270 km/h** if the 2550 KV motors survive the pitch. Modelled 290 km/h is the drag-limited ceiling with the stated load factor. Dive + tailwind is how marketing numbers get to 245 on 2100 KV; this build is meant to beat that on a straight.

---

## Files

| File | |
|---|---|
| `speed_model.py` | 4S/6S, 3"/5"/5.5" ranking |
| `test_speed_model.py` | unit tests |
| `bom.json` | parts |
| `wiring.md` | UART / power / antennas |
| `betaflight/speed5_cli.txt` | starting dump |
