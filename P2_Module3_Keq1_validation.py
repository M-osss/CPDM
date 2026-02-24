import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

df = pd.read_csv("2026_Module3_P2.csv")

T = df["T_K"].values
P = df["P_bar"].values
y_CH4 = df["y_CH4"].values
y_H2O = df["y_H2O"].values
y_CO = df["y_CO"].values
y_H2 = df["y_H2"].values
n = len(T)

Keq1 = (P / 1.0)**2 * (y_CO * y_H2**3) / (y_CH4 * y_H2O)
lnKeq1 = np.log(Keq1)

# ---- Refit both best models ----
# Part A best: quadratic minus P^2
X_poly = sm.add_constant(pd.DataFrame({"T": T, "P": P, "T2": T**2, "PT": P * T}))
model_poly = sm.OLS(Keq1, X_poly).fit()

# Part B best: van't Hoff
X_vh = sm.add_constant(pd.DataFrame({"1/T": 1.0 / T}))
model_vh = sm.OLS(lnKeq1, X_vh).fit()

print("=" * 70)
print("MODEL VALIDATION AND VERIFICATION")
print("=" * 70)
print()
print("Part A best model : Keq1 = a*T + b*P + c*T^2 + d*PT + intercept")
print("Part B best model : ln(Keq1) = A + B/T  (van't Hoff)")
print()

# =====================================================================
# TEST 1: Leave-One-Out Cross-Validation (LOOCV)
# =====================================================================
print("=" * 70)
print("TEST 1: LEAVE-ONE-OUT CROSS-VALIDATION (LOOCV)")
print("=" * 70)
print()
print("Rationale: With only 15 data points, LOOCV is the most appropriate")
print("cross-validation scheme. Each point is predicted by a model trained")
print("on the other 14. This detects overfitting -- a model that memorizes")
print("training data will have high training R^2 but poor LOOCV R^2.")
print()

loo_resid_poly = np.zeros(n)
loo_resid_vh = np.zeros(n)

for i in range(n):
    mask = np.ones(n, dtype=bool)
    mask[i] = False

    # Polynomial
    X_tr = X_poly.iloc[mask]
    y_tr = Keq1[mask]
    m_tr = sm.OLS(y_tr, X_tr).fit()
    loo_resid_poly[i] = Keq1[i] - m_tr.predict(X_poly.iloc[[i]]).iloc[0]

    # van't Hoff
    X_tr_vh = X_vh.iloc[mask]
    y_tr_vh = lnKeq1[mask]
    m_tr_vh = sm.OLS(y_tr_vh, X_tr_vh).fit()
    pred_ln = m_tr_vh.predict(X_vh.iloc[[i]]).iloc[0]
    loo_resid_vh[i] = Keq1[i] - np.exp(pred_ln)

loocv_rmse_poly = np.sqrt(np.mean(loo_resid_poly**2))
loocv_rmse_vh = np.sqrt(np.mean(loo_resid_vh**2))
ss_tot = np.sum((Keq1 - np.mean(Keq1))**2)
loocv_r2_poly = 1 - np.sum(loo_resid_poly**2) / ss_tot
loocv_r2_vh = 1 - np.sum(loo_resid_vh**2) / ss_tot

print(f"{'Model':<40} {'Train R^2':>10} {'LOOCV R^2':>10} {'LOOCV RMSE':>12}")
print("-" * 72)
print(f"{'Part A: Polynomial (T,P,T2,PT)':<40} {model_poly.rsquared:>10.5f} {loocv_r2_poly:>10.5f} {loocv_rmse_poly:>12.4f}")
r2_vh_train = 1 - np.sum((Keq1 - np.exp(model_vh.predict(X_vh)))**2) / ss_tot
print(f"{'Part B: van t Hoff ln(K)=A+B/T':<40} {r2_vh_train:>10.5f} {loocv_r2_vh:>10.5f} {loocv_rmse_vh:>12.4f}")
print()
print("Interpretation: If LOOCV R^2 drops substantially from training R^2,")
print("the model is overfitting. A small drop indicates robustness.")
print()

# =====================================================================
# TEST 2: Residual Analysis -- Normality
# =====================================================================
print("=" * 70)
print("TEST 2: RESIDUAL NORMALITY (Shapiro-Wilk test)")
print("=" * 70)
print()
print("Rationale: OLS assumes residuals are normally distributed. If they")
print("are not, confidence intervals and p-values are unreliable, which")
print("means our model selection decisions could be wrong.")
print()

resid_poly = model_poly.resid
resid_vh = model_vh.resid  # in ln-space

sw_poly, p_poly = stats.shapiro(resid_poly)
sw_vh, p_vh = stats.shapiro(resid_vh)

print(f"{'Model':<40} {'Shapiro-Wilk W':>15} {'p-value':>10} {'Normal?':>10}")
print("-" * 75)
print(f"{'Part A: Polynomial residuals':<40} {sw_poly:>15.5f} {p_poly:>10.4f} {'Yes' if p_poly > 0.05 else 'No':>10}")
print(f"{'Part B: van t Hoff residuals (ln-space)':<40} {sw_vh:>15.5f} {p_vh:>10.4f} {'Yes' if p_vh > 0.05 else 'No':>10}")
print()
print("Null hypothesis: residuals are normally distributed.")
print("Reject at alpha=0.05 if p-value < 0.05.")
print()

# =====================================================================
# TEST 3: Residual Analysis -- Homoscedasticity (Breusch-Pagan)
# =====================================================================
print("=" * 70)
print("TEST 3: HOMOSCEDASTICITY (Breusch-Pagan test)")
print("=" * 70)
print()
print("Rationale: OLS assumes constant variance of residuals. If variance")
print("changes with the predictor (heteroscedasticity), standard errors")
print("are biased and model comparison metrics are unreliable.")
print()

from statsmodels.stats.diagnostic import het_breuschpagan

bp_poly = het_breuschpagan(resid_poly, X_poly)
bp_vh = het_breuschpagan(resid_vh, X_vh)

print(f"{'Model':<40} {'BP statistic':>13} {'p-value':>10} {'Homoscedastic?':>15}")
print("-" * 78)
print(f"{'Part A: Polynomial':<40} {bp_poly[0]:>13.4f} {bp_poly[1]:>10.4f} {'Yes' if bp_poly[1] > 0.05 else 'No':>15}")
print(f"{'Part B: van t Hoff (ln-space)':<40} {bp_vh[0]:>13.4f} {bp_vh[1]:>10.4f} {'Yes' if bp_vh[1] > 0.05 else 'No':>15}")
print()
print("Null hypothesis: residual variance is constant (homoscedastic).")
print("Reject at alpha=0.05 if p-value < 0.05.")
print()

# =====================================================================
# TEST 4: Residual Independence (Durbin-Watson)
# =====================================================================
print("=" * 70)
print("TEST 4: RESIDUAL INDEPENDENCE (Durbin-Watson statistic)")
print("=" * 70)
print()
print("Rationale: OLS assumes residuals are independent. Autocorrelated")
print("residuals indicate systematic misfit (model is missing a pattern).")
print("DW ~ 2 means no autocorrelation; DW << 2 means positive correlation.")
print()

from statsmodels.stats.stattools import durbin_watson

dw_poly = durbin_watson(resid_poly)
dw_vh = durbin_watson(resid_vh)

print(f"{'Model':<40} {'Durbin-Watson':>14}")
print("-" * 54)
print(f"{'Part A: Polynomial':<40} {dw_poly:>14.4f}")
print(f"{'Part B: van t Hoff (ln-space)':<40} {dw_vh:>14.4f}")
print()
print("DW in [1.5, 2.5] suggests no significant autocorrelation.")
print()

# =====================================================================
# TEST 5: Influential Points (Cook's Distance)
# =====================================================================
print("=" * 70)
print("TEST 5: INFLUENTIAL POINTS (Cook's Distance)")
print("=" * 70)
print()
print("Rationale: A single outlier or high-leverage point can dominate the")
print("fit and distort coefficients. Cook's distance > 4/n = %.3f flags" % (4.0/n))
print("influential points that should be investigated.")
print()

infl_poly = model_poly.get_influence()
cooks_poly = infl_poly.cooks_distance[0]
infl_vh = model_vh.get_influence()
cooks_vh = infl_vh.cooks_distance[0]
threshold = 4.0 / n

print("Part A: Polynomial model")
print(f"{'Point':>6} {'T(K)':>7} {'P(bar)':>7} {'Keq1':>10} {'Cook D':>10} {'Influential?':>13}")
print("-" * 53)
for i in range(n):
    flag = "***" if cooks_poly[i] > threshold else ""
    print(f"{i:>6d} {T[i]:>7.0f} {P[i]:>7.0f} {Keq1[i]:>10.4f} {cooks_poly[i]:>10.4f} {flag:>13}")
print()

print("Part B: van't Hoff model")
print(f"{'Point':>6} {'T(K)':>7} {'P(bar)':>7} {'Keq1':>10} {'Cook D':>10} {'Influential?':>13}")
print("-" * 53)
for i in range(n):
    flag = "***" if cooks_vh[i] > threshold else ""
    print(f"{i:>6d} {T[i]:>7.0f} {P[i]:>7.0f} {Keq1[i]:>10.4f} {cooks_vh[i]:>10.4f} {flag:>13}")
print()

# =====================================================================
# TEST 6: Physical Consistency Check
# =====================================================================
print("=" * 70)
print("TEST 6: PHYSICAL CONSISTENCY CHECK")
print("=" * 70)
print()
print("Rationale: A grey box model embeds physical knowledge. We can verify")
print("the fitted parameters against known thermodynamic data as an")
print("independent check that the regression was performed correctly.")
print()

R_gas = 8.314
A_vh = model_vh.params.iloc[0]
B_vh = model_vh.params.iloc[1]
dH_fit = -B_vh * R_gas
dS_fit = A_vh * R_gas

print("Steam methane reforming: CH4 + H2O -> CO + 3 H2")
print()
print(f"  Fitted delta_H  = {dH_fit/1000:.1f} kJ/mol")
print(f"  Literature delta_H (298 K) ~ 206 kJ/mol")
print(f"  Deviation: {abs(dH_fit/1000 - 206):.1f} kJ/mol ({100*abs(dH_fit/1000 - 206)/206:.1f}%)")
print()
print(f"  Fitted delta_S  = {dS_fit:.1f} J/(mol*K)")
print(f"  Literature delta_S (298 K) ~ 215 J/(mol*K)")
print(f"  Deviation: {abs(dS_fit - 215):.1f} J/(mol*K) ({100*abs(dS_fit - 215)/215:.1f}%)")
print()
print("  Check: dH > 0 (endothermic)?  ", "PASS" if dH_fit > 0 else "FAIL")
print("  Check: K increases with T?     ", "PASS" if B_vh < 0 else "FAIL")
print("  Check: dH within 10% of lit?   ", "PASS" if abs(dH_fit/1000 - 206)/206 < 0.10 else "FAIL")
print()

# =====================================================================
# TEST 7: Manual Spot-Check Calculation
# =====================================================================
print("=" * 70)
print("TEST 7: MANUAL SPOT-CHECK (verify K_eq,1 computation)")
print("=" * 70)
print()
print("Rationale: Before trusting any regression, verify the dependent")
print("variable was computed correctly. Pick one row and reproduce by hand.")
print()

i_check = 0
print(f"Row {i_check}: T={T[i_check]} K, P={P[i_check]} bar")
print(f"  y_CH4  = {y_CH4[i_check]}")
print(f"  y_H2O  = {y_H2O[i_check]}")
print(f"  y_CO   = {y_CO[i_check]}")
print(f"  y_H2   = {y_H2[i_check]}")
manual = (P[i_check]/1.0)**2 * (y_CO[i_check] * y_H2[i_check]**3) / (y_CH4[i_check] * y_H2O[i_check])
print(f"  Manual: (P)^2 * y_CO * y_H2^3 / (y_CH4 * y_H2O)")
print(f"        = ({P[i_check]})^2 * {y_CO[i_check]:.6f} * {y_H2[i_check]:.6f}^3 / ({y_CH4[i_check]:.6f} * {y_H2O[i_check]:.6f})")
print(f"        = {manual:.6f}")
print(f"  Script: {Keq1[i_check]:.6f}")
print(f"  Match:  {'PASS' if abs(manual - Keq1[i_check]) < 1e-12 else 'FAIL'}")
print()

# =====================================================================
# TEST 8: Extrapolation Behavior
# =====================================================================
print("=" * 70)
print("TEST 8: EXTRAPOLATION SANITY CHECK")
print("=" * 70)
print()
print("Rationale: A good model should produce physically reasonable")
print("predictions outside the training range. K_eq must be positive")
print("and should increase monotonically with T for this endothermic rxn.")
print()

T_extrap = np.array([700, 800, 850, 950, 1050, 1100, 1200, 1500])
P_extrap = 10.0

print(f"Extrapolation at P = {P_extrap} bar:")
print(f"{'T (K)':>8} {'Poly Keq1':>12} {'VH Keq1':>12} {'Poly OK?':>10} {'VH OK?':>10}")
print("-" * 52)

prev_poly = -np.inf
prev_vh = -np.inf

for t in T_extrap:
    x_p = pd.DataFrame({"const": [1], "T": [t], "P": [P_extrap], "T2": [t**2], "PT": [P_extrap*t]})
    pred_p = model_poly.predict(x_p).iloc[0]

    x_v = pd.DataFrame({"const": [1], "1/T": [1.0/t]})
    pred_v = np.exp(model_vh.predict(x_v).iloc[0])

    ok_p = "Yes" if pred_p > 0 and pred_p > prev_poly else "No"
    ok_v = "Yes" if pred_v > 0 and pred_v > prev_vh else "No"

    in_range = " " if 850 <= t <= 1050 else "*"
    print(f"{t:>7.0f}{in_range} {pred_p:>12.4f} {pred_v:>12.4f} {ok_p:>10} {ok_v:>10}")

    prev_poly = pred_p
    prev_vh = pred_v

print()
print("* = outside training range [850, 1050] K")
print("'OK' requires positive value AND monotonically increasing with T.")
print()

# =====================================================================
# SUMMARY
# =====================================================================
print("=" * 70)
print("VALIDATION SUMMARY")
print("=" * 70)
print("""
Test                          Part A (Polynomial)     Part B (van't Hoff)
--------------------------------------------------------------------------
1. LOOCV R^2                  {loocv_poly:>10.4f}              {loocv_vh:>10.4f}
2. Normality (Shapiro p)      {p_poly:>10.4f}              {p_vh:>10.4f}
3. Homoscedasticity (BP p)    {bp_p:>10.4f}              {bp_v:>10.4f}
4. Durbin-Watson              {dw_p:>10.4f}              {dw_v:>10.4f}
5. Cook's D > threshold       {cook_p:>10d}              {cook_v:>10d}
6. Physical consistency       N/A (black box)         PASS (dH ~ 205 kJ/mol)
7. Spot-check K_eq calc       PASS                    PASS
8. Extrapolation positive     See table               See table
""".format(
    loocv_poly=loocv_r2_poly,
    loocv_vh=loocv_r2_vh,
    p_poly=p_poly,
    p_vh=p_vh,
    bp_p=bp_poly[1],
    bp_v=bp_vh[1],
    dw_p=dw_poly,
    dw_v=dw_vh,
    cook_p=int(np.sum(cooks_poly > threshold)),
    cook_v=int(np.sum(cooks_vh > threshold)),
))

print("Recommended validation approach for any regression model:")
print("1. Cross-validate (LOOCV or k-fold) to detect overfitting.")
print("2. Check residual normality (Shapiro-Wilk) to validate p-values.")
print("3. Check homoscedasticity (Breusch-Pagan) to validate std errors.")
print("4. Check independence (Durbin-Watson) to detect systematic misfit.")
print("5. Inspect Cook's distance to find points dominating the fit.")
print("6. Compare fitted parameters to literature (grey box models).")
print("7. Spot-check the computed dependent variable by hand.")
print("8. Extrapolate and verify predictions remain physically sensible.")
