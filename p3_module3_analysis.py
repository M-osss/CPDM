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
