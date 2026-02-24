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
