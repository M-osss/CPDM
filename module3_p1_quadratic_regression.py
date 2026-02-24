#!/usr/bin/env python3
"""Fit full quadratic regression for Module 3 Problem 1 data.

Model form:
Y = b0 + b1*A + b2*B + b3*C + b4*A^2 + b5*B^2 + b6*C^2 + b7*A*B + b8*B*C + b9*A*C
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


def load_module3_p1(csv_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    a_vals: list[float] = []
    b_vals: list[float] = []
    c_vals: list[float] = []
    y_vals: list[float] = []

    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            a_vals.append(float(row["A"]))
            b_vals.append(float(row["B"]))
            c_vals.append(float(row["C"]))
            y_vals.append(float(row["Y_percent_reduction"]))

    return (
        np.asarray(a_vals, dtype=float),
        np.asarray(b_vals, dtype=float),
        np.asarray(c_vals, dtype=float),
        np.asarray(y_vals, dtype=float),
    )


def build_design_matrix(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    return np.column_stack(
        [
            np.ones_like(a),
            a,
            b,
            c,
            a * a,
            b * b,
            c * c,
            a * b,
            b * c,
            a * c,
        ]
    )


def fit_quadratic_model(x: np.ndarray, y: np.ndarray) -> dict[str, np.ndarray | float | int]:
    beta, _, rank, _ = np.linalg.lstsq(x, y, rcond=None)
    y_hat = x @ beta
    residuals = y - y_hat

    n_obs = y.shape[0]
    n_params_without_intercept = x.shape[1] - 1
    sse = float(np.sum(residuals**2))
    sst = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = 1.0 - (sse / sst)
    adj_r_squared = 1.0 - (1.0 - r_squared) * (n_obs - 1) / (n_obs - n_params_without_intercept - 1)
    rmse = float(np.sqrt(sse / (n_obs - n_params_without_intercept - 1)))

    return {
        "beta": beta,
        "rank": rank,
        "y_hat": y_hat,
        "residuals": residuals,
        "sse": sse,
        "r_squared": r_squared,
        "adj_r_squared": adj_r_squared,
        "rmse": rmse,
    }


def main() -> None:
    csv_path = Path(__file__).resolve().parent / "2026_Module3_P1.csv"
    a, b, c, y = load_module3_p1(csv_path)
    x = build_design_matrix(a, b, c)
    result = fit_quadratic_model(x, y)

    beta = result["beta"]
    names = ["b0", "b1_A", "b2_B", "b3_C", "b4_A2", "b5_B2", "b6_C2", "b7_AB", "b8_BC", "b9_AC"]

    print("Fitted coefficients")
    for name, value in zip(names, beta):
        print(f"{name} = {value:.10f}")

    print("\nFitted equation")
    print(
        "Y = "
        f"{beta[0]:.10f} + ({beta[1]:.10f})A + ({beta[2]:.10f})B + ({beta[3]:.10f})C + "
        f"({beta[4]:.10f})A^2 + ({beta[5]:.10f})B^2 + ({beta[6]:.10f})C^2 + "
        f"({beta[7]:.10f})AB + ({beta[8]:.10f})BC + ({beta[9]:.10f})AC"
    )

    print("\nModel diagnostics")
    print(f"rank = {result['rank']}")
    print(f"SSE = {result['sse']:.10f}")
    print(f"R^2 = {result['r_squared']:.10f}")
    print(f"Adjusted R^2 = {result['adj_r_squared']:.10f}")
    print(f"RMSE = {result['rmse']:.10f}")


if __name__ == "__main__":
    main()
