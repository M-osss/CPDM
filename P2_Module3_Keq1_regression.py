import numpy as np
import pandas as pd
import statsmodels.api as sm
from itertools import combinations

df = pd.read_csv("2026_Module3_P2.csv")

T = df["T_K"].values
P = df["P_bar"].values
y_CH4 = df["y_CH4"].values
y_H2O = df["y_H2O"].values
y_CO = df["y_CO"].values
y_H2 = df["y_H2"].values

Keq1 = (P / 1.0)**2 * (y_CO * y_H2**3) / (y_CH4 * y_H2O)

df["Keq1"] = Keq1
print("=" * 70)
print("Computed K_eq,1 values")
print("=" * 70)
print(df[["T_K", "P_bar", "Keq1"]].to_string(index=False))
print()

# ---------------------------------------------------------------------------
# LINEAR MODEL: Keq1 = aT + bP + intercept
# ---------------------------------------------------------------------------
print("=" * 70)
print("LINEAR MODEL: Keq1 = a*T + b*P + intercept")
print("=" * 70)

X_lin = pd.DataFrame({"T": T, "P": P})
X_lin_const = sm.add_constant(X_lin)
model_lin = sm.OLS(Keq1, X_lin_const).fit()
print(model_lin.summary())
print()

# Linear with only T
print("-" * 70)
print("LINEAR MODEL (T only): Keq1 = a*T + intercept")
print("-" * 70)
X_T = sm.add_constant(pd.DataFrame({"T": T}))
model_lin_T = sm.OLS(Keq1, X_T).fit()
print(model_lin_T.summary())
print()

# Linear with only P
print("-" * 70)
print("LINEAR MODEL (P only): Keq1 = b*P + intercept")
print("-" * 70)
X_P = sm.add_constant(pd.DataFrame({"P": P}))
model_lin_P = sm.OLS(Keq1, X_P).fit()
print(model_lin_P.summary())
print()

# ---------------------------------------------------------------------------
# QUADRATIC MODEL: Keq1 = aT + bP + cT^2 + dP^2 + ePT + intercept
# ---------------------------------------------------------------------------
print("=" * 70)
print("QUADRATIC MODEL: Keq1 = a*T + b*P + c*T^2 + d*P^2 + e*P*T + intercept")
print("=" * 70)

X_quad = pd.DataFrame({
    "T": T,
    "P": P,
    "T2": T**2,
    "P2": P**2,
    "PT": P * T,
})
X_quad_const = sm.add_constant(X_quad)
model_quad = sm.OLS(Keq1, X_quad_const).fit()
print(model_quad.summary())
print()

# ---------------------------------------------------------------------------
# Quadratic model elimination: try all subsets, rank by adjusted R^2 and AIC
# ---------------------------------------------------------------------------
print("=" * 70)
print("QUADRATIC SUBSET SELECTION (all subsets, ranked by Adj. R^2)")
print("=" * 70)

features = ["T", "P", "T2", "P2", "PT"]
results = []

for k in range(1, len(features) + 1):
    for combo in combinations(features, k):
        X_sub = sm.add_constant(X_quad[list(combo)])
        m = sm.OLS(Keq1, X_sub).fit()
        results.append({
            "features": combo,
            "n_features": k,
            "R2": m.rsquared,
            "Adj_R2": m.rsquared_adj,
            "AIC": m.aic,
            "BIC": m.bic,
        })

results_df = pd.DataFrame(results).sort_values("Adj_R2", ascending=False)
print(results_df.head(15).to_string(index=False))
print()

best_row = results_df.iloc[0]
best_features = list(best_row["features"])
print("=" * 70)
print(f"BEST SUBSET: {best_features}")
print(f"  Adj R^2 = {best_row['Adj_R2']:.6f}")
print(f"  AIC     = {best_row['AIC']:.4f}")
print(f"  BIC     = {best_row['BIC']:.4f}")
print("=" * 70)
print()

# Refit and print summary for the best subset
X_best = sm.add_constant(X_quad[best_features])
model_best = sm.OLS(Keq1, X_best).fit()
print("BEST MODEL SUMMARY:")
print(model_best.summary())
print()

# ---------------------------------------------------------------------------
# Also compare linear vs quadratic full models
# ---------------------------------------------------------------------------
print("=" * 70)
print("MODEL COMPARISON SUMMARY")
print("=" * 70)
print(f"{'Model':<30} {'R^2':>8} {'Adj R^2':>10} {'AIC':>12} {'BIC':>12}")
print("-" * 72)
models = {
    "Linear (T+P)": model_lin,
    "Linear (T only)": model_lin_T,
    "Linear (P only)": model_lin_P,
    "Quadratic (full)": model_quad,
    f"Best subset {best_features}": model_best,
}
for name, m in models.items():
    print(f"{name:<30} {m.rsquared:>8.5f} {m.rsquared_adj:>10.5f} {m.aic:>12.4f} {m.bic:>12.4f}")
print()

all_features_set = set(features)
eliminated = all_features_set - set(best_features)
print(f"Factors eliminated from quadratic model: {eliminated if eliminated else 'None'}")
print(f"Metric(s) used for selection: Adjusted R-squared (primary), AIC/BIC (secondary)")
