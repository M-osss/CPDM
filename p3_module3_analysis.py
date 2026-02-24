import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy import stats

df = pd.read_csv("2026_Module3_P3.csv")
df.columns = df.columns.str.strip()

df.rename(columns={"X_gL": "X", "Y_gL": "Y", "K_log": "K"}, inplace=True)
df["T"] = (df["Thickener"] == "B").astype(int)

print("=" * 70)
print("DATA SUMMARY")
print("=" * 70)
print(df.to_string(index=False))
print(f"\nN = {len(df)}")
print(df[["X", "Y", "K"]].describe())

# =====================================================================
# PART A: Quadratic regression K ~ X, Y (ignore T)
# =====================================================================
print("\n" + "=" * 70)
print("PART A: QUADRATIC REGRESSION (ignoring thickener type)")
print("=" * 70)

# Model 1: Quadratic without interaction
print("\n--- Model 1: K = b0 + b1*X + b2*Y + b3*X^2 + b4*Y^2 ---")
model1 = ols("K ~ X + Y + I(X**2) + I(Y**2)", data=df).fit()
print(model1.summary())

# Model 2: Quadratic with interaction
print("\n--- Model 2: K = b0 + b1*X + b2*Y + b3*X^2 + b4*Y^2 + b5*X*Y ---")
model2 = ols("K ~ X + Y + I(X**2) + I(Y**2) + X:Y", data=df).fit()
print(model2.summary())

# Partial F-test for the interaction term
rss_reduced = model1.ssr
rss_full = model2.ssr
df_num = model2.df_model - model1.df_model
df_den = model2.df_resid
F_stat = ((rss_reduced - rss_full) / df_num) / (rss_full / df_den)
p_value = 1 - stats.f.cdf(F_stat, df_num, df_den)

print("\n--- Partial F-test for interaction term X*Y ---")
print(f"  RSS (without XY): {rss_reduced:.4f}")
print(f"  RSS (with XY):    {rss_full:.4f}")
print(f"  F-statistic:      {F_stat:.4f}")
print(f"  p-value:          {p_value:.4f}")
print(f"  Conclusion:       {'Interaction IS significant (synergistic/antagonistic effect present)' if p_value < 0.05 else 'Interaction is NOT significant at alpha=0.05'}")

print("\n--- Model Comparison ---")
print(f"  Model 1 (no interaction): R² = {model1.rsquared:.4f}, Adj R² = {model1.rsquared_adj:.4f}, AIC = {model1.aic:.2f}, BIC = {model1.bic:.2f}")
print(f"  Model 2 (with interaction): R² = {model2.rsquared:.4f}, Adj R² = {model2.rsquared_adj:.4f}, AIC = {model2.aic:.2f}, BIC = {model2.bic:.2f}")

# Also show a simple linear model for comparison
print("\n--- Model 0 (linear only): K = b0 + b1*X + b2*Y ---")
model0 = ols("K ~ X + Y", data=df).fit()
print(model0.summary())

# Interpretation of interaction coefficient
if model2.pvalues["X:Y"] < 0.05:
    coeff_xy = model2.params["X:Y"]
    if coeff_xy < 0:
        print(f"\nThe interaction coefficient b5 = {coeff_xy:.4f} is negative and significant.")
        print("This indicates an ANTAGONISTIC effect: increasing both X and Y together")
        print("reduces K more than would be expected from their individual effects alone.")
    else:
        print(f"\nThe interaction coefficient b5 = {coeff_xy:.4f} is positive and significant.")
        print("This indicates a SYNERGISTIC effect: increasing both X and Y together")
        print("enhances K more than would be expected from their individual effects alone.")
else:
    print(f"\nThe interaction coefficient (p = {model2.pvalues['X:Y']:.4f}) is NOT significant.")
    print("No evidence of synergistic or antagonistic interaction between X and Y.")

# =====================================================================
# PART B: Include thickener type T as binary variable
# =====================================================================
print("\n" + "=" * 70)
print("PART B: INCLUDE THICKENER TYPE T (T=0 for A, T=1 for B)")
print("=" * 70)

# Full model: add T to the best Part A model (quadratic + interaction)
print("\n--- Full Model: K = b0 + b1*X + b2*Y + b3*X^2 + b4*Y^2 + b5*X*Y + b6*T ---")
model_full = ols("K ~ X + Y + I(X**2) + I(Y**2) + X:Y + T", data=df).fit()
print(model_full.summary())

# Partial F-test: does T improve the model?
rss_no_T = model2.ssr
rss_with_T = model_full.ssr
df_num_T = model_full.df_model - model2.df_model
df_den_T = model_full.df_resid
F_T = ((rss_no_T - rss_with_T) / df_num_T) / (rss_with_T / df_den_T)
p_T = 1 - stats.f.cdf(F_T, df_num_T, df_den_T)
print(f"\n--- Partial F-test for T ---")
print(f"  RSS (without T): {rss_no_T:.4f}")
print(f"  RSS (with T):    {rss_with_T:.4f}")
print(f"  F-statistic:     {F_T:.4f}")
print(f"  p-value:         {p_T:.4f}")
print(f"  Conclusion:      {'T IS significant — thickener impacts K' if p_T < 0.05 else 'T is NOT significant at alpha=0.05'}")

# =====================================================================
# MODEL SIMPLIFICATION: backward elimination
# =====================================================================
print("\n" + "=" * 70)
print("MODEL SIMPLIFICATION (backward elimination, alpha = 0.05)")
print("=" * 70)

terms_sequence = [
    ("K ~ X + Y + I(X**2) + I(Y**2) + X:Y + T", "Full model with T"),
]

current_formula = "K ~ X + Y + I(X**2) + I(Y**2) + X:Y + T"
step = 0
while True:
    step += 1
    m = ols(current_formula, data=df).fit()
    print(f"\n--- Step {step}: {current_formula} ---")
    print(m.summary())

    pvals = m.pvalues.drop("Intercept")
    worst_term = pvals.idxmax()
    worst_p = pvals.max()

    if worst_p > 0.05:
        print(f"\n  >> Removing '{worst_term}' (p = {worst_p:.4f}) — not significant at alpha=0.05")
        term_map = {
            "I(X ** 2)": "I(X**2)",
            "I(Y ** 2)": "I(Y**2)",
            "X:Y": "X:Y",
        }
        remove_str = term_map.get(worst_term, worst_term)
        parts = [t.strip() for t in current_formula.split("~")[1].split("+")]
        parts = [t for t in parts if t != remove_str]
        current_formula = "K ~ " + " + ".join(parts)
    else:
        print(f"\n  >> All remaining terms significant (worst p = {worst_p:.4f}). Stopping.")
        break

print("\n" + "=" * 70)
print("FINAL SIMPLIFIED MODEL")
print("=" * 70)
final_model = ols(current_formula, data=df).fit()
print(final_model.summary())
print(f"\nFinal formula: {current_formula}")
print(f"R² = {final_model.rsquared:.4f}, Adj R² = {final_model.rsquared_adj:.4f}")
print(f"AIC = {final_model.aic:.2f}, BIC = {final_model.bic:.2f}")

print("\n--- Coefficient Interpretations ---")
for name, val in final_model.params.items():
    p = final_model.pvalues[name]
    print(f"  {name:15s}: coeff = {val:+.4f}, p = {p:.4f}")

# =====================================================================
# PART C: Plot the model and find optimal formula
# =====================================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import minimize

b = model_full.params
def K_model(X, Y, T):
    return (b["Intercept"] + b["X"]*X + b["Y"]*Y
            + b["I(X ** 2)"]*X**2 + b["I(Y ** 2)"]*Y**2
            + b["X:Y"]*X*Y + b["T"]*T)

# --- Figure 1: Contour plots for both thickeners ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

x_range = np.linspace(0, 10, 300)
y_range = np.linspace(0, 10, 300)
Xg, Yg = np.meshgrid(x_range, y_range)
mask = Xg + Yg <= 10

for ax, T_val, label in zip(axes, [0, 1], ["Thickener A (T=0)", "Thickener B (T=1)"]):
    Kg = K_model(Xg, Yg, T_val)
    Kg_masked = np.where(mask, Kg, np.nan)

    contour = ax.contourf(Xg, Yg, Kg_masked, levels=20, cmap="RdYlGn")
    cs = ax.contour(Xg, Yg, Kg_masked, levels=20, colors="k", linewidths=0.3)
    ax.clabel(cs, inline=True, fontsize=7, fmt="%.1f")
    fig.colorbar(contour, ax=ax, label="K (log units)")

    ax.plot([0, 10], [10, 0], "k--", linewidth=1.5, label="X + Y = 10")

    sub = df[df["T"] == T_val]
    ax.scatter(sub["X"], sub["Y"], c="black", edgecolors="white", s=50, zorder=5, label="Data points")

    ax.set_xlabel("X (g/L)", fontsize=12)
    ax.set_ylabel("Y (g/L)", fontsize=12)
    ax.set_title(label, fontsize=13)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.legend(loc="upper right", fontsize=9)

fig.suptitle(
    r"$K = 1.671 + 0.635X + 0.955Y - 0.037X^2 - 0.073Y^2 - 0.279XY + 0.549T$",
    fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig("part_c_contour.png", dpi=200, bbox_inches="tight")
plt.close()
print("\nSaved: part_c_contour.png")

# --- Figure 2: 3D surface for thickener B ---
fig2 = plt.figure(figsize=(10, 7))
ax3d = fig2.add_subplot(111, projection="3d")

Kg_B = K_model(Xg, Yg, 1)
Kg_B_masked = np.where(mask, Kg_B, np.nan)

surf = ax3d.plot_surface(Xg, Yg, Kg_B_masked, cmap="RdYlGn",
                         edgecolor="none", alpha=0.9)
fig2.colorbar(surf, ax=ax3d, shrink=0.5, label="K (log units)")

sub_B = df[df["T"] == 1]
ax3d.scatter(sub_B["X"], sub_B["Y"], sub_B["K"], c="black", s=40, zorder=5,
             depthshade=False, label="Data (Thickener B)")

ax3d.set_xlabel("X (g/L)")
ax3d.set_ylabel("Y (g/L)")
ax3d.set_zlabel("K (log units)")
ax3d.set_title("Thickener B: K vs X, Y", fontsize=13)
ax3d.view_init(elev=30, azim=225)
ax3d.legend()
plt.savefig("part_c_3d_surface.png", dpi=200, bbox_inches="tight")
plt.close()
print("Saved: part_c_3d_surface.png")

# =====================================================================
# OPTIMIZATION: maximize K subject to X+Y <= 10, X >= 0, Y >= 0
# =====================================================================
print("\n" + "=" * 70)
print("OPTIMIZATION: Best hand sanitizer formula")
print("=" * 70)

for T_val, T_name in [(0, "A"), (1, "B")]:
    result = minimize(
        lambda xy: -K_model(xy[0], xy[1], T_val),
        x0=[2, 2],
        bounds=[(0, None), (0, None)],
        constraints=[{"type": "ineq", "fun": lambda xy: 10 - xy[0] - xy[1]}],
        method="SLSQP"
    )
    X_opt, Y_opt = result.x
    K_opt = -result.fun
    print(f"\n  Thickener {T_name} (T={T_val}):")
    print(f"    X* = {X_opt:.4f} g/L")
    print(f"    Y* = {Y_opt:.4f} g/L")
    print(f"    X* + Y* = {X_opt + Y_opt:.4f} g/L")
    print(f"    K* = {K_opt:.4f} log units")

# Also check boundary / corner points for thickener B
print("\n--- Grid search verification (Thickener B) ---")
best_K = -np.inf
best_xy = (0, 0)
for xi in np.linspace(0, 10, 10001):
    for yi in [0, 10 - xi]:
        if yi < 0:
            continue
        k_val = K_model(xi, yi, 1)
        if k_val > best_K:
            best_K = k_val
            best_xy = (xi, yi)

# Also check interior critical point via gradient = 0
# dK/dX = b1 + 2*b3*X + b5*Y = 0
# dK/dY = b2 + 2*b4*Y + b5*X = 0
A_mat = np.array([[2*b["I(X ** 2)"], b["X:Y"]],
                  [b["X:Y"], 2*b["I(Y ** 2)"]]])
b_vec = np.array([-b["X"], -b["Y"]])
try:
    crit = np.linalg.solve(A_mat, b_vec)
    X_crit, Y_crit = crit
    if X_crit >= 0 and Y_crit >= 0 and X_crit + Y_crit <= 10:
        K_crit = K_model(X_crit, Y_crit, 1)
        print(f"  Interior critical point: X={X_crit:.4f}, Y={Y_crit:.4f}, K={K_crit:.4f}")
    else:
        print(f"  Interior critical point ({X_crit:.4f}, {Y_crit:.4f}) outside feasible region")
except np.linalg.LinAlgError:
    print("  No interior critical point (singular Hessian)")

print(f"  Boundary grid best: X={best_xy[0]:.4f}, Y={best_xy[1]:.4f}, K={best_K:.4f}")

# Final optimal
T_opt = 1
res_final = minimize(
    lambda xy: -K_model(xy[0], xy[1], T_opt),
    x0=[2, 2],
    bounds=[(0, None), (0, None)],
    constraints=[{"type": "ineq", "fun": lambda xy: 10 - xy[0] - xy[1]}],
    method="SLSQP"
)
X_star, Y_star = res_final.x
K_star = -res_final.fun

print("\n" + "=" * 70)
print("OPTIMAL HAND SANITIZER FORMULA")
print("=" * 70)
print(f"  Compound X concentration: {X_star:.4f} g/L")
print(f"  Compound Y concentration: {Y_star:.4f} g/L")
print(f"  X + Y = {X_star + Y_star:.4f} g/L  (< 10 g/L constraint satisfied)")
print(f"  Thickener: B")
print(f"  Predicted K = {K_star:.4f} log units")
print(f"  ({10**(K_star):.1f}-fold germ reduction, i.e. {(1 - 10**(-K_star))*100:.4f}% kill rate)")

# --- Figure 3: Contour with optimum marked ---
fig3, ax3 = plt.subplots(figsize=(8, 7))
Kg_B = K_model(Xg, Yg, 1)
Kg_B_masked = np.where(mask, Kg_B, np.nan)
contour = ax3.contourf(Xg, Yg, Kg_B_masked, levels=25, cmap="RdYlGn")
cs = ax3.contour(Xg, Yg, Kg_B_masked, levels=25, colors="k", linewidths=0.3)
ax3.clabel(cs, inline=True, fontsize=7, fmt="%.1f")
fig3.colorbar(contour, ax=ax3, label="K (log units)")
ax3.plot([0, 10], [10, 0], "k--", linewidth=1.5, label="X + Y = 10")
ax3.plot(X_star, Y_star, "r*", markersize=20, markeredgecolor="black",
         markeredgewidth=1.2, zorder=10,
         label=f"Optimum ({X_star:.2f}, {Y_star:.2f})\nK = {K_star:.2f}")
sub_B = df[df["T"] == 1]
ax3.scatter(sub_B["X"], sub_B["Y"], c="black", edgecolors="white", s=50, zorder=5,
            label="Data (Thickener B)")
ax3.set_xlabel("X (g/L)", fontsize=12)
ax3.set_ylabel("Y (g/L)", fontsize=12)
ax3.set_title("Thickener B: Optimal Sanitizer Formula", fontsize=13)
ax3.set_xlim(0, 10)
ax3.set_ylim(0, 10)
ax3.legend(loc="upper right", fontsize=10)
plt.tight_layout()
plt.savefig("part_c_optimum.png", dpi=200, bbox_inches="tight")
plt.close()
print("\nSaved: part_c_optimum.png")
