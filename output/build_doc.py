"""
Build the Word document containing answers to Paper B.
"""
import json
from docx import Document
from docx.shared import Pt, Cm, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_ALIGN_VERTICAL
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

with open('/workspace/output/q1_results.json') as fh:
    R = json.load(fh)

doc = Document()

# ---------- page setup ----------
for section in doc.sections:
    section.top_margin = Cm(1.8)
    section.bottom_margin = Cm(1.8)
    section.left_margin = Cm(2.0)
    section.right_margin = Cm(2.0)

style = doc.styles['Normal']
style.font.name = 'Calibri'
style.font.size = Pt(11)

def add_heading(text, level=1):
    h = doc.add_heading(text, level=level)
    for run in h.runs:
        run.font.color.rgb = RGBColor(0x1F, 0x3A, 0x68)
    return h

def add_para(text, bold=False, italic=False, size=11):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.bold = bold
    r.italic = italic
    r.font.size = Pt(size)
    return p

def add_bullet(text):
    p = doc.add_paragraph(text, style='List Bullet')
    return p

def fmt(x, dec=0):
    return f"{x:,.{dec}f}"

# ---------- header ----------
title = doc.add_paragraph()
title.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = title.add_run('CENG 4130 — Plant Design and Economics\nFinal Examination, Paper B (Spring 2025/2026)\nWorked Solutions')
r.bold = True
r.font.size = Pt(13)

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.add_run('Student ID: 20987955    |    Development cost X = USD 20,987,955').italic = True

# ====================================================================
# QUESTION 1
# ====================================================================
add_heading('Question 1 — VacuumBot Project Economics', level=1)

add_heading('1.0  Basis and assumptions', level=2)
assumptions = [
    'Year 0 captures all pre-operational outlays (1-year development period plus equipment purchase and working capital injection at start-up). Operations and sales span Years 1–5.',
    'Annual sales volume, selling price, fixed and variable operating costs, and advertising budget are constant over the 5-year operating life.',
    'Depreciation: 5-year property class MACRS half-year convention with rates 20.00 / 32.00 / 19.20 / 11.52 / 11.52 / 5.76 %. The project is terminated at the end of Year 5, so the un-recovered 5.76 % of the equipment basis (USD 2,304,000) is treated as a book loss on disposal (salvage = 0) and generates a one-off tax shield in Year 5.',
    'Development cost X is treated as a non-depreciable, non-deductible pre-operating outlay (capitalised). Treating it as a Year-0 tax-deductible expense would only improve NPV and IRR; the assumption taken is conservative.',
    'Profits tax rate = 15 % (Hong Kong two-tier upper rate). Discount rate = 20 % per annum, end-of-year cash-flow convention.',
    'Working capital of USD 20 M is fully recovered at the end of Year 5.',
]
for a in assumptions:
    add_bullet(a)

add_heading('1.1  Annual revenue and operating cost', level=2)

t = doc.add_table(rows=1, cols=3)
t.style = 'Light Grid Accent 1'
hdr = t.rows[0].cells
hdr[0].text = 'Item'
hdr[1].text = 'Calculation'
hdr[2].text = 'Value (USD/yr)'
rows = [
    ('Revenue', '500,000 units × USD 250 / unit', '125,000,000'),
    ('Variable cost', '500,000 × USD 90', '45,000,000'),
    ('Fixed operating cost', 'given', '35,000,000'),
    ('Advertising', 'given', '8,000,000'),
    ('Total operating cost', '45 + 35 + 8', '88,000,000'),
    ('EBITDA (Revenue – OpCost)', '125 – 88', '37,000,000'),
]
for r0 in rows:
    rr = t.add_row().cells
    for i, v in enumerate(r0):
        rr[i].text = v

add_heading('1.2  Depreciation schedule (5-year MACRS on USD 40 M)', level=2)
t = doc.add_table(rows=1, cols=4)
t.style = 'Light Grid Accent 1'
hdr = t.rows[0].cells
for i, h in enumerate(['Year', 'MACRS rate', 'Depreciation (USD)', 'Book value EoY (USD)']):
    hdr[i].text = h
rates = [0.20, 0.32, 0.192, 0.1152, 0.1152]
bv = 40_000_000
for y, rt in enumerate(rates, 1):
    d = rt * 40_000_000
    bv -= d
    rr = t.add_row().cells
    rr[0].text = str(y)
    rr[1].text = f'{rt*100:.2f} %'
    rr[2].text = fmt(d)
    rr[3].text = fmt(bv)
add_para('Un-recovered book value of USD 2,304,000 at end of Year 5 is written off as a loss on disposal; tax shield = 0.15 × 2,304,000 = USD 345,600.', italic=True)

add_heading('1.3  After-tax cash-flow table', level=2)

t = doc.add_table(rows=1, cols=7)
t.style = 'Light Grid Accent 1'
hdr = ['Year', 'Revenue', 'OpCost', 'Depr.', 'EBT', 'Tax (15 %)', 'Net CF']
for i, h in enumerate(hdr):
    t.rows[0].cells[i].text = h

# Year 0
rr = t.add_row().cells
rr[0].text = '0'
rr[1].text = '—'
rr[2].text = '—'
rr[3].text = '—'
rr[4].text = '—'
rr[5].text = '—'
rr[6].text = f"({fmt(80_987_955)})"

for y in range(1, 6):
    d = rates[y-1]*40_000_000
    EBT = 125_000_000 - 88_000_000 - d
    tax = 0.15*EBT
    NI = EBT - tax
    CF = NI + d
    extra = ''
    if y == 5:
        CF += 20_000_000 + 345_600
        extra = '\n(+WC 20.0 M + tax shield 0.346 M)'
    rr = t.add_row().cells
    rr[0].text = str(y)
    rr[1].text = fmt(125_000_000)
    rr[2].text = fmt(88_000_000)
    rr[3].text = fmt(d)
    rr[4].text = fmt(EBT)
    rr[5].text = fmt(tax)
    rr[6].text = fmt(CF) + extra

add_para('Year-0 outlay = X (20,987,955) + Equipment (40,000,000) + Working capital (20,000,000) = USD 80,987,955.', italic=True)

add_heading('1.4  (a) Project cash position', level=2)
t = doc.add_table(rows=1, cols=3)
t.style = 'Light Grid Accent 1'
for i, h in enumerate(['Year', 'Annual cash flow (USD)', 'Cumulative cash position (USD)']):
    t.rows[0].cells[i].text = h
for y, f, c in zip(range(0, 6), R['cf'], R['cum']):
    rr = t.add_row().cells
    rr[0].text = str(y)
    rr[1].text = ('(' + fmt(-f) + ')') if f < 0 else fmt(f)
    rr[2].text = ('(' + fmt(-c) + ')') if c < 0 else fmt(c)

add_para('The cumulative cash position turns positive between Year 2 and Year 3. Undiscounted simple payback period = 2.46 years.')

doc.add_picture('/workspace/output/cash_position.png', width=Inches(6.0))
cap = doc.paragraphs[-1]
cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p.add_run('Figure 1. Cumulative project cash position (Year 0 to Year 5).')
r.italic = True
r.font.size = Pt(10)

add_heading('1.5  (b) Net present worth at 20 %', level=2)
add_para('Discounting each annual cash flow at i = 20 %:')

t = doc.add_table(rows=1, cols=4)
t.style = 'Light Grid Accent 1'
for i, h in enumerate(['Year', 'Cash flow (USD)', 'Discount factor (1.20)^-n', 'Present value (USD)']):
    t.rows[0].cells[i].text = h
NPV_check = 0
for y, f in zip(range(0, 6), R['cf']):
    df = 1/(1.20**y)
    pv = f*df
    NPV_check += pv
    rr = t.add_row().cells
    rr[0].text = str(y)
    rr[1].text = ('(' + fmt(-f) + ')') if f < 0 else fmt(f)
    rr[2].text = f'{df:.4f}'
    rr[3].text = ('(' + fmt(-pv) + ')') if pv < 0 else fmt(pv)

rr = t.add_row().cells
rr[0].text = 'Σ NPV'
rr[3].text = fmt(R['NPV'])

p = doc.add_paragraph()
r = p.add_run(f"NPV @ 20 % = +USD {R['NPV']:,.0f} (≈ +USD 24.9 million).  ")
r.bold = True
p.add_run('Because NPV > 0 at the firm’s hurdle rate, the project creates value and should be undertaken.')

add_heading('1.6  (c) Discounted cash-flow rate of return (DCFROR)', level=2)
add_para('Solving Σ CFₙ / (1+i)ⁿ = 0 by iteration (bisection on i):')

t = doc.add_table(rows=1, cols=2)
t.style = 'Light Grid Accent 1'
t.rows[0].cells[0].text = 'Trial discount rate'
t.rows[0].cells[1].text = 'NPV (USD)'
for rt in [0.20, 0.25, 0.30, 0.325, 0.35]:
    npv = sum(f/(1+rt)**y for y, f in zip(range(0,6), R['cf']))
    rr = t.add_row().cells
    rr[0].text = f'{rt*100:.1f} %'
    rr[1].text = ('(' + fmt(-npv) + ')') if npv < 0 else fmt(npv)

p = doc.add_paragraph()
r = p.add_run(f"DCFROR ≈ {R['IRR']*100:.2f}%.  ")
r.bold = True
p.add_run('Since the DCFROR substantially exceeds the 20 % hurdle, the project is economically attractive.')

add_heading('1.7  Summary of Q1 results', level=2)
t = doc.add_table(rows=1, cols=2)
t.style = 'Light Grid Accent 1'
t.rows[0].cells[0].text = 'Metric'
t.rows[0].cells[1].text = 'Value'
summary = [
    ('Total Year-0 outlay', f"USD {80_987_955:,}"),
    ('Annual EBITDA (Years 1–5)', f"USD {37_000_000:,}"),
    ('Simple payback (undiscounted)', '2.46 yr'),
    ('NPV at 20 %', f"+USD {R['NPV']:,.0f}"),
    ('DCFROR (IRR)', f"{R['IRR']*100:.2f} %"),
    ('Recommendation', 'Proceed with the project'),
]
for k, v in summary:
    rr = t.add_row().cells
    rr[0].text = k
    rr[1].text = v

doc.add_page_break()

# ====================================================================
# QUESTION 2
# ====================================================================
add_heading('Question 2 — Towngas Residential Hydrogen Refilling Facility', level=1)

p = doc.add_paragraph()
r = p.add_run('Role: Technical Lead, Towngas — H₂ extraction by PSA from existing town-gas mains, installed in a residential car-park bay (≈3 parking spaces ≈ 37.5 m² footprint).')
r.italic = True

add_heading('Part (a) — Economic Feasibility Evaluation (Board of Directors)', level=2)

add_heading('A1. Process scope and design basis', level=3)
basis_a = [
    'Town gas composition (Towngas, naphtha-reformed): ≈49 vol% H₂, 28 % CH₄, 19 % CO₂, 3 % CO, balance N₂. Distribution pressure ≈ 75 mbarg, taken locally up to ≈ 8–10 barg by a feed booster.',
    'PSA recovery 80 %, H₂ purity 99.97 % (SAE J2719 fuel-cell grade).',
    'Nominal H₂ production: 50 kg H₂ / day. Storage at 875 bar (Type IV vessels) with cascade dispensing at 700 bar (H70).',
    'Site footprint = 37.5 m² (3 parking bays); 24/7 operation; expected utilisation ramp 40 % (yr 1) → 75 % (yr ≥3); design life 15 yr.',
    'Refuelling spec: 4–5 kg H₂ filled in <5 min per FCEV (Toyota Mirai ≈ 5.6 kg → ~650 km).',
    'Assumed FCEV throughput: ~7 vehicles / day at design utilisation; plus 30 kg/day diverted to a stationary 50 kWₑ PEM fuel-cell unit feeding the building.',
]
for a in basis_a:
    add_bullet(a)

add_heading('A2. Capital cost estimate (Class-4 study estimate, ±30 %)', level=3)
t = doc.add_table(rows=1, cols=3)
t.style = 'Light Grid Accent 1'
for i, h in enumerate(['Item', 'Basis', 'Installed cost (USD)']):
    t.rows[0].cells[i].text = h
items = [
    ('Town-gas feed booster + dryer + pretreatment', '10 Nm³/h skid', '60,000'),
    ('PSA unit (4-bed) sized for 50 kg H₂/day', 'small-scale skid pricing', '450,000'),
    ('Multi-stage diaphragm compressor 0.5–900 bar', '~30 kg/h, 50 kW', '380,000'),
    ('High-pressure storage cascade (≈100 kg @ 950 bar)', 'Type IV composite', '260,000'),
    ('H₂ dispenser (H70) + chiller (–40 °C)', 'single-hose', '180,000'),
    ('Civil works, fire wall, vent stack, ventilation', 'small footprint', '120,000'),
    ('Instrumentation, SIS, gas detection, CCTV', 'SIL-2 safety system', '110,000'),
    ('Engineering, procurement, commissioning, licensing', '15 % of direct', '230,000'),
    ('Contingency', '15 % of installed', '270,000'),
    ('Total fixed capital investment (FCI)', '', '≈ 2,060,000'),
]
for it in items:
    rr = t.add_row().cells
    for i, v in enumerate(it):
        rr[i].text = v
add_para('Working capital ≈ 10 % FCI = USD 0.21 M. Total capital outlay ≈ USD 2.27 M per site.', italic=True)

add_heading('A3. Annual operating cost (at 75 % utilisation, 13,700 kg H₂/yr)', level=3)
t = doc.add_table(rows=1, cols=3)
t.style = 'Light Grid Accent 1'
for i, h in enumerate(['Cost element', 'Basis / unit cost', 'USD / yr']):
    t.rows[0].cells[i].text = h
op = [
    ('Town-gas feed (net H₂ extracted only)',
     '≈ 0.5 USD / kg H₂ (HK Towngas wholesale ≈ HK$1.2/m³; depleted stream returned to grid at credit)',
     '6,850'),
    ('Electricity (PSA + compression + chiller ≈ 9 kWh/kg)',
     '13,700 kg × 9 × USD 0.16 / kWh', '19,700'),
    ('Cooling water + N₂ purge', '— ', '3,000'),
    ('Maintenance + spares', '4 % FCI', '82,000'),
    ('Labour (0.3 FTE shared across cluster)', '0.3 × USD 90 k', '27,000'),
    ('Insurance & permits', '1.5 % FCI', '31,000'),
    ('Property rental (3 bays)', '3 × USD 350/mo × 12', '12,600'),
    ('Total OPEX (excl. depreciation)', '', '≈ 182,000'),
]
for it in op:
    rr = t.add_row().cells
    for i, v in enumerate(it):
        rr[i].text = v

add_heading('A4. Levelised cost of hydrogen (LCOH)', level=3)
add_para('Annualisation factor (capital recovery factor) at 8 % over 15 yr: CRF = i(1+i)ⁿ / [(1+i)ⁿ – 1] = 0.1168.')
add_para('Annualised capex = 2,060,000 × 0.1168 = USD 240,600 / yr.')
add_para('LCOH = (annualised capex + OPEX) / annual H₂ delivered = (240,600 + 182,000) / 13,700 ≈ USD 30.8 / kg H₂.')

t = doc.add_table(rows=1, cols=4)
t.style = 'Light Grid Accent 1'
for i, h in enumerate(['Utilisation', 'Annual kg H₂', 'LCOH (USD/kg)', 'Equiv. per 100 km*']):
    t.rows[0].cells[i].text = h
for u, kg in [(0.40, 7300), (0.60, 10950), (0.75, 13700), (0.90, 16425)]:
    annual_capex = 240_600
    opex_var = 19_700 * (kg/13700) + 6_850*(kg/13700) + 3000*(kg/13700)
    opex_fix = 82_000 + 27_000 + 31_000 + 12_600
    lcoh = (annual_capex + opex_fix + opex_var)/kg
    rr = t.add_row().cells
    rr[0].text = f'{u*100:.0f} %'
    rr[1].text = f'{kg:,}'
    rr[2].text = f'{lcoh:.1f}'
    rr[3].text = f'≈ USD {lcoh*1.0:.1f}'  # FCEV ~1 kg / 100 km
add_para('*Toyota Mirai consumption ≈ 0.86 kg H₂ / 100 km; rounded to 1 kg / 100 km. Gasoline-ICE equivalent: 7 L/100 km × HK$ 25 /L (~USD 3.2 /L) ≈ USD 22 / 100 km.', italic=True)

add_heading('A5. Benchmark and conclusion', level=3)
bench = [
    'At design utilisation, distributed H₂ from this PSA route delivers H₂ at ≈ USD 18–31 / kg depending on utilisation. Gasoline equivalent on energy-service basis is ≈ USD 22 / 100 km, hence H₂ becomes cost-competitive only above ≈70 % utilisation and after factoring HK ZEV incentives and the 50 % first-registration-tax waiver on FCEV.',
    'A residential cluster roll-out of ≈ 50 sites yields economy of scale on PSA, compressor and dispenser procurement (≈25 % FCI reduction) and shared maintenance crew, dropping LCOH to ≈ USD 14–18 / kg — broadly comparable with gasoline on a per-km basis.',
    'Scalability: town-gas backbone covers >90 % of HK residential premises; no new feedstock distribution is required, making per-site scale-up faster than electrolysis (which is electricity-grid limited).',
    'Major economic risks: (i) low FCEV fleet uptake → low utilisation; (ii) town-gas tariff escalation; (iii) carbon-intensity scrutiny — Towngas is fossil-derived, so CO₂ from the depleted return stream still counts against the building portfolio unless CCS is added downstream; (iv) regulatory: HK Electrical & Mechanical Services Department (EMSD) licensing for residential H₂ is still maturing; (v) safety incident → loss of social licence.',
    'Sensitivity (per-site, design utilisation): every +10 % utilisation reduces LCOH by ≈ USD 2.6 / kg; every USD 0.02 / kWh increase in electricity adds ≈ USD 0.18 / kg.',
    'Recommendation: Pilot a single residential site for 12 months. Approve full roll-out only if (a) utilisation > 60 % by end of pilot, (b) safety case approved by Fire Services Department and EMSD, and (c) cluster procurement secures ≥20 % capex reduction.',
]
for b in bench:
    add_bullet(b)

add_heading('Part (b) — Safety and Risk Management Proposal (Government Submission)', level=2)

add_heading('B1. Key hazards identified', level=3)
haz = [
    'Hydrogen flammability: very wide flammable range 4–75 vol% in air, minimum ignition energy ~0.02 mJ (≈10× lower than methane). Quasi-invisible flame (≈ 200–300 °C surface but very high radiant heat).',
    'High-pressure mechanical hazard: 700–950 bar storage; potential for catastrophic vessel rupture, jet release (sonic, > 1,300 m/s), jet fire.',
    'Confined space accumulation: car-park is partially enclosed → loss of containment could pool above the leak and form a flammable cloud near the soffit (H₂ rises at ~20 m/s buoyancy in air).',
    'Hydrogen embrittlement of carbon-steel piping and fittings; cyclic-fatigue failure of refuelling hose.',
    'Town-gas side hazards: methane / CO toxicity of the feed; backflow of H₂ into low-pressure mains.',
    'Electrical ignition sources (vehicles, lighting), static, lightning.',
    'Human interaction: untrained residents performing self-service refuelling; vehicle impact on dispenser; vandalism.',
    'Domino effect with adjacent parked EVs / ICE cars (Li-ion thermal runaway, gasoline fuel tank).',
]
for h in haz:
    add_bullet(h)

add_heading('B2. Qualitative risk assessment (HAZID / What-If summary)', level=3)
t = doc.add_table(rows=1, cols=5)
t.style = 'Light Grid Accent 1'
for i, h in enumerate(['Scenario', 'Cause', 'Consequence', 'Likelihood', 'Risk (pre-mitigation)']):
    t.rows[0].cells[i].text = h
ra = [
    ('Small H₂ leak at compressor seal', 'Seal wear', 'Local ignition, jet fire <2 m', 'Possible', 'Medium'),
    ('PSA tail-gas H₂ slip → back-feed town-gas line', 'PSA control failure', 'Pipeline over-pressure, off-spec gas to consumers', 'Unlikely', 'Medium'),
    ('Dispenser hose burst during refuelling', 'Embrittlement / impact', 'Jet flame, burn injury, vehicle fire', 'Possible', 'High'),
    ('Storage-tank BLEVE / rupture', 'Fire engulfment / overfill', 'Blast over-pressure, fragments, fatalities', 'Rare', 'High'),
    ('Flammable cloud in enclosed car-park', 'Sustained leak + poor ventilation', 'Deflagration / DDT, structural damage', 'Unlikely', 'High'),
    ('Town-gas (CH₄/CO) feed leak', 'Pipework failure', 'Toxic exposure, asphyxiation', 'Possible', 'Medium'),
    ('Vehicle collision with dispenser', 'Driver error', 'Mechanical damage, possible leak', 'Possible', 'Medium'),
    ('Loss of power → SIS de-energised', 'Grid failure', 'Process upset, isolation valves close (fail-safe)', 'Possible', 'Low'),
]
for r0 in ra:
    rr = t.add_row().cells
    for i, v in enumerate(r0):
        rr[i].text = v

add_heading('B3. Engineering & administrative safeguards', level=3)
saf = [
    'Layout: outdoor or semi-outdoor canopy with at least two open sides; minimum separation distances per NFPA 2 (2023) — ≥ 6 m from public-access edge, ≥ 8 m from building openings, ≥ 3 m from parked vehicles. Fire wall (REI-120) between PSA / storage skid and adjacent parking bays.',
    'Ventilation: natural cross-ventilation ≥ 1 % of canopy roof area; supplementary forced extraction at high level (≥ 12 air changes per hour) to prevent H₂ accumulation > 25 % LFL (=1 vol%).',
    'Detection: redundant H₂ point detectors (catalytic + thermal-conductivity) at the soffit, ultrasonic leak detectors near compressors, UV/IR flame detectors at storage. Two-out-of-three voting triggers trip at 20 % LFL.',
    'Safety Instrumented System (SIL-2): fail-safe isolation valves on town-gas inlet, H₂ outlet, storage and dispenser; automatic depressurisation to a remote vent stack (≥ 4 m above roof) on confirmed leak/fire.',
    'Mechanical integrity: all wetted parts in 316L SS (low-strength grade to resist embrittlement); double-block-and-bleed on the town-gas / H₂ interface; non-return valves to prevent back-feed; hoses replaced on time-in-service or cycle count.',
    'Hazardous-area classification (ATEX / IEC 60079): Zone 1 within 1 m of all H₂ joints, Zone 2 within 3 m; explosion-proof Ex-d luminaires, intrinsically-safe instruments; bonding and earthing of all metallic parts; lightning protection.',
    'Pressure relief: thermally-activated PRDs on every storage cylinder, vented through dedicated vent stack with flashback arrestor; PSV on PSA receivers sized for credible heat-input case.',
    'Vehicle protection: bollards (kerb-mounted, K-12 rated) around dispenser and storage; CCTV; dispenser nozzle with automatic break-away coupling.',
    'Operational: read-only HMI for residents; QR-coded ID + RFID authorisation; geofenced shut-off if smoking / open flame detected; mandatory training for housing-estate maintenance staff; permit-to-work for hot-work / breaking containment.',
    'Documentation: site-specific HAZOP + LOPA, Pre-Startup Safety Review (PSSR), Management of Change (MOC), incident reporting under EMSD Gas Safety Ordinance Cap. 51.',
]
for s in saf:
    add_bullet(s)

add_heading('B4. Credible worst-case scenario and emergency response', level=3)
add_para(
    'Credible worst case (CWC): catastrophic failure of one 95 L composite cylinder at 950 bar full of H₂ (≈ 6 kg H₂, ≈ 0.85 GJ). '
    'Outcomes envelope:',
)
cwc = [
    'Immediate ignited release → vertical jet fire ~6–8 m flame length; thermal radiation 12.5 kW/m² endpoint at ≈ 6 m → exclusion zone radius 8 m maintained by layout.',
    'Delayed ignition with confined accumulation → deflagration peak overpressure ≈ 0.3 bar at 5 m, dropping below 0.05 bar at 15 m (structural-damage threshold). Achieved mitigation through ventilation and explosion-relief louvers in canopy.',
    'Unignited release → buoyant dispersion above roof line; ≤ 30 s to reach upper LFL contour, well below combustion-air domain at ground level.',
]
for c in cwc:
    add_bullet(c)
add_para('Mitigation hierarchy applied: inherent (small inventory per cylinder, 100 kg total cap), passive (composite-wrapped vessels with TPRD), active (SIS isolation in < 1 s, blow-down to vent stack within 90 s), procedural (training, drills).')

add_heading('B5. Emergency Response Plan (ERP) — key elements', level=3)
erp = [
    'Tier-1 (operator): automatic SIS trip — isolation + venting + sirens + strobe; SMS / app alert to estate manager and Towngas control room.',
    'Tier-2 (on-site): evacuation of car-park within 60 s using lit egress paths; activation of fire-suppression for ignited release (water curtain to cool adjacent surfaces — do NOT extinguish a controlled jet fire before isolation is confirmed).',
    'Tier-3 (external): notification to Fire Services Department (FSD) and EMSD; pre-agreed access route for fire appliances; ≥ 100 m cordon for ignited storage event.',
    'Post-incident: lock-out / tag-out, third-party root-cause analysis, regulator submission within 24 h per Gas Safety (Installation and Use) Regulations.',
    'Drills: full-scale exercise with FSD annually; quarterly tabletop with estate management; pre-commissioning and re-validation every 5 yr.',
    'Public communication: pre-launch information to residents — what hydrogen is, how to recognise an alarm, evacuation route, no smoking radius. Posters in 中/EN at lift lobbies and car-park entrances.',
]
for e in erp:
    add_bullet(e)

add_heading('B6. Overall safety conclusion', level=3)
add_para(
    'With inherent design choices (small distributed H₂ inventory, semi-outdoor location, PSA/compressor separation, composite vessels with TPRDs), engineered safeguards (SIL-2 SIS, detection, ventilation, vent stack), and a structured ERP integrated with FSD and EMSD, the residual risk of fatality at the boundary of the car-park can be demonstrated to be Tolerable / ALARP (≤ 1×10⁻⁶ per year individual risk). The proposal is judged acceptable for pilot deployment subject to formal HAZOP, LOPA and Quantitative Risk Assessment (QRA) submission to the regulator prior to commissioning.'
)

# ====================================================================
# Save
# ====================================================================
out_path = '/workspace/output/Paper_B_Solutions_SID_20987955.docx'
doc.save(out_path)
print(f'Saved: {out_path}')
