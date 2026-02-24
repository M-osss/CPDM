"""
Three-factor Box-Behnken design for 2026_Module3_P1.csv data.
Generates design matrix and fits second-order response surface model.
"""

import csv
import numpy as np
from pathlib import Path


def generate_box_behnken_3factor(n_center: int = 3) -> np.ndarray:
    """
    Generate 3-factor Box-Behnken design matrix (coded -1, 0, 1).
    12 edge runs + n_center center replicates.
    """
    # Edge midpoints: two factors at ±1, one at 0
    # (A,B) plane, C=0
    ab_plane = np.array([
        [-1, -1, 0], [-1, 1, 0], [1, -1, 0], [1, 1, 0]
    ])
    # (A,C) plane, B=0
    ac_plane = np.array([
        [-1, 0, -1], [-1, 0, 1], [1, 0, -1], [1, 0, 1]
    ])
    # (B,C) plane, A=0
    bc_plane = np.array([
        [0, -1, -1], [0, -1, 1], [0, 1, -1], [0, 1, 1]
    ])
    edges = np.vstack([ab_plane, ac_plane, bc_plane])
    centers = np.zeros((n_center, 3))
    return np.vstack([edges, centers])


def fit_rsm_second_order(X: np.ndarray, y: np.ndarray) -> dict:
    """
    Fit second-order response surface: y = b0 + sum(bi*xi) + sum(bii*xi^2) + sum(bij*xi*xj)
    """
    n = X.shape[0]
    # Build design matrix: [1, x1, x2, x3, x1^2, x2^2, x3^2, x1*x2, x1*x3, x2*x3]
    X_rsm = np.column_stack([
        np.ones(n),
        X[:, 0], X[:, 1], X[:, 2],
        X[:, 0]**2, X[:, 1]**2, X[:, 2]**2,
        X[:, 0]*X[:, 1], X[:, 0]*X[:, 2], X[:, 1]*X[:, 2]
    ])
    coeffs, residuals, rank, s = np.linalg.lstsq(X_rsm, y, rcond=None)
    y_pred = X_rsm @ coeffs
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    return {
        'coefficients': coeffs,
        'r_squared': r_squared,
        'y_pred': y_pred,
        'terms': ['intercept', 'A', 'B', 'C', 'A2', 'B2', 'C2', 'AB', 'AC', 'BC']
    }


def main():
    workspace = Path(__file__).parent
    csv_path = workspace / '2026_Module3_P1.csv'
    if not csv_path.exists():
        csv_path = Path('2026_Module3_P1.csv')

    with open(csv_path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    factors = ['A', 'B', 'C']
    response = 'Y_percent_reduction'
    design = np.array([[float(r[f]) for f in factors] for r in rows])
    y = np.array([float(r[response]) for r in rows])

    # Generate reference design
    ref_design = generate_box_behnken_3factor(n_center=3)
    print("Generated 3-factor Box-Behnken design (12 edges + 3 center):")
    for i, row in enumerate(ref_design):
        print(f"  {i+1:2d}: A={row[0]:4.1f} B={row[1]:4.1f} C={row[2]:4.1f}")

    # Fit RSM
    rsm = fit_rsm_second_order(design, y)
    print(f"\nSecond-order RSM fit: R² = {rsm['r_squared']:.4f}")
    print("\nCoefficients:")
    for term, coef in zip(rsm['terms'], rsm['coefficients']):
        print(f"  {term}: {coef:.4f}")

    # Output design to CSV (design only, no response)
    out_path = workspace / 'box_behnken_design_3factor.csv'
    with open(out_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(factors)
        w.writerows(ref_design.tolist())
    print(f"\nDesign written to {out_path}")


if __name__ == '__main__':
    main()
