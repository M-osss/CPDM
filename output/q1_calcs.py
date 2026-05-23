"""
Question 1 cash flow analysis for CENG 4130 Paper B.
SID = 20987955 -> development cost X = USD 20,987,955.
"""
import numpy as np
import matplotlib.pyplot as plt

# ---------- Inputs ----------
SID = 20_987_955
X = SID                       # development cost (USD)
units = 500_000               # units sold per year
price = 250                   # USD/unit
revenue = units * price       # USD/yr = 125,000,000

advertising = 8_000_000
equip = 40_000_000
WC = 20_000_000
fixed_op = 35_000_000
var_unit = 90
variable_op = var_unit * units   # 45,000,000
tax_rate = 0.15
disc = 0.20

# 5-year MACRS depreciation rates (IRS half-year convention)
macrs = [0.20, 0.32, 0.192, 0.1152, 0.1152, 0.0576]
dep = [r * equip for r in macrs]      # USD per year (yrs 1..6)
# project life = 5 yrs -> the un-recovered 5.76% becomes a book loss at disposal (salvage = 0)
dep_y1_y5 = dep[:5]
book_value_end = equip - sum(dep_y1_y5)   # = 0.0576 * 40M

# ---------- Annual operating cash flow ----------
op_cost = fixed_op + variable_op + advertising      # 88,000,000
EBITDA = revenue - op_cost                          # 37,000,000

print(f"X = development cost = USD {X:,}")
print(f"Annual revenue      = USD {revenue:,}")
print(f"Annual op cost      = USD {op_cost:,}")
print(f"Annual EBITDA       = USD {EBITDA:,}")
print(f"Annual depreciation (Y1..Y5): {[f'{d:,.0f}' for d in dep_y1_y5]}")
print(f"Remaining book value at end Y5 = USD {book_value_end:,.0f}\n")

cf = {}
for yr in range(1, 6):
    d = dep_y1_y5[yr-1]
    EBT = EBITDA - d
    tax = tax_rate * EBT
    NI = EBT - tax
    operating_cf = NI + d
    cf[yr] = operating_cf
    print(f"Yr {yr}: dep={d:,.0f}, EBT={EBT:,.0f}, tax={tax:,.0f}, NI={NI:,.0f}, CF={operating_cf:,.0f}")

# Terminal (end of year 5)
loss_on_disposal = book_value_end            # salvage = 0 -> entire BV is a tax loss
tax_shield_disposal = tax_rate * loss_on_disposal
terminal = WC + tax_shield_disposal
cf[5] += terminal
print(f"\nTerminal Y5: WC recovery {WC:,} + tax shield on disposal "
      f"{tax_shield_disposal:,.0f}  -> Y5 total CF = {cf[5]:,.0f}")

# Year 0 (initial outlay = development + equipment + working capital)
cf[0] = -(X + equip + WC)
print(f"\nYr 0 outlay = -(X + Equip + WC) = USD {cf[0]:,.0f}")

# ---------- Cumulative cash position ----------
years = list(range(0, 6))
flows = [cf[y] for y in years]
cum = np.cumsum(flows)
print("\nYear |  Cash flow (USD)  | Cumulative position (USD)")
for y, f, c in zip(years, flows, cum):
    print(f"  {y}  | {f:>16,.0f} | {c:>16,.0f}")

# Payback (simple, undiscounted)
for i in range(1, len(years)):
    if cum[i-1] < 0 <= cum[i]:
        pb = (i-1) + (-cum[i-1]) / flows[i]
        print(f"\nUndiscounted payback period ≈ {pb:.2f} yr")
        break

# ---------- NPV ----------
NPV = sum(f / (1+disc)**y for y, f in zip(years, flows))
print(f"\nNPV at {disc*100:.0f}% = USD {NPV:,.0f}  -> "
      f"{'ACCEPT' if NPV > 0 else 'REJECT'} project")

# ---------- IRR (DCFROR) ----------
def npv(rate):
    return sum(f / (1+rate)**y for y, f in zip(years, flows))

lo, hi = 0.05, 1.0
for _ in range(80):
    mid = 0.5*(lo+hi)
    if npv(mid) > 0:
        lo = mid
    else:
        hi = mid
IRR = 0.5*(lo+hi)
print(f"\nDCFROR (IRR) ≈ {IRR*100:.2f}%")

# Check a few discount rates
print("\nNPV vs rate:")
for r in [0.10, 0.15, 0.20, 0.25, 0.30, 0.325, 0.35]:
    print(f"  r = {r*100:5.1f}%  ->  NPV = {npv(r):>15,.0f}")

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(years, cum/1e6, marker='o', linewidth=2, color='#1f4e79')
ax.axhline(0, color='k', linewidth=0.8)
ax.fill_between(years, cum/1e6, 0,
                where=(cum < 0), color='#fde0dc', alpha=0.7, label='Cash deficit')
ax.fill_between(years, cum/1e6, 0,
                where=(cum >= 0), color='#dcf0dc', alpha=0.7, label='Cash surplus')
ax.set_xlabel('Project year')
ax.set_ylabel('Cumulative cash position (USD million)')
ax.set_title('Project Cash Position over the 5-Year Operating Life')
ax.grid(True, linestyle=':')
for y, c in zip(years, cum):
    ax.annotate(f'{c/1e6:,.1f}', (y, c/1e6),
                textcoords='offset points', xytext=(0, 10),
                ha='center', fontsize=9)
ax.legend(loc='lower right')
plt.tight_layout()
plt.savefig('/workspace/output/cash_position.png', dpi=160)
print('\nSaved plot -> /workspace/output/cash_position.png')

# Save numerical table for the doc
import json
out = {
    'X': X,
    'revenue': revenue,
    'op_cost': op_cost,
    'EBITDA': EBITDA,
    'depreciation': dep_y1_y5,
    'book_value_end': book_value_end,
    'cf': flows,
    'cum': cum.tolist(),
    'NPV': NPV,
    'IRR': IRR,
}
with open('/workspace/output/q1_results.json', 'w') as fh:
    json.dump(out, fh, indent=2)
