"""Sensitivity analysis for Stage-2 PFR."""
from __future__ import annotations

import csv
import itertools
import json
import sys
import time
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

try:
    from .pfr_heattransfer import Stage2State, solve_stage2
    from .pfr_ideal import get_species_index, make_diluted_feed, FEED_FORMULA
except ImportError:
    from pfr_heattransfer import Stage2State, solve_stage2
    from pfr_ideal import get_species_index, make_diluted_feed, FEED_FORMULA

ROOT = Path(__file__).resolve().parent


def run_single_case(
    T_in: float,
    P_in: float,
    T_furnace: float,
    v_z0: float,
    n2_fraction: float,
    D_inner: float = 0.05,
    L: float = 10.0,
    verbose: bool = False,
) -> Dict[str, Any]:
    feed = make_diluted_feed("C2H6", n2_fraction)
    
    state = Stage2State(
        D_inner=D_inner,
        L=L,
        T_in=T_in,
        P_in=P_in,
        v_z0=v_z0,
        T_furnace=T_furnace,
        Nz=50,
        Nr_wall=3,
    )
    
    try:
        result = solve_stage2(
            state, 
            feed=feed, 
            max_outer_iter=30,
            verbose=False
        )
        
        spec_idx = get_species_index()
        Y_out = result["Y"][:, -1]
        Y_out = np.clip(Y_out, 0, None)
        Y_out /= Y_out.sum()
        
        Y_in = result["Y"][:, 0]
        Y_in = np.clip(Y_in, 0, None)
        Y_in /= Y_in.sum()
        
        i_eth = spec_idx.get("C2H6", -1)
        i_ene = spec_idx.get("C2H4", -1)
        i_ace = spec_idx.get("C2H2", -1)
        i_h2 = spec_idx.get("H2", -1)
        i_n2 = spec_idx.get("N2", -1)
        
        Y_C2H6_out = Y_out[i_eth] if i_eth >= 0 else 0
        Y_C2H4_out = Y_out[i_ene] if i_ene >= 0 else 0
        Y_C2H2_out = Y_out[i_ace] if i_ace >= 0 else 0
        Y_H2_out = Y_out[i_h2] if i_h2 >= 0 else 0
        Y_N2_out = Y_out[i_n2] if i_n2 >= 0 else 0
        
        Y_C2H6_in = Y_in[i_eth] if i_eth >= 0 else 0
        
        conversion = 1.0 - Y_C2H6_out / Y_C2H6_in if Y_C2H6_in > 1e-10 else 0
        ethane_reacted = Y_C2H6_in - Y_C2H6_out
        selectivity_C2H4 = Y_C2H4_out / ethane_reacted if ethane_reacted > 1e-10 else 0
        selectivity_C2H2 = Y_C2H2_out / ethane_reacted if ethane_reacted > 1e-10 else 0
        
        T_out = result["T_gas"][-1]
        T_wall_out = result["T_wall_inner"][-1]
        q_avg = np.mean(result["q_inner"]) / 1e3
        
        return {
            "T_in_K": T_in,
            "P_in_bar": P_in / 1e5,
            "T_furnace_K": T_furnace,
            "v_z0_m_s": v_z0,
            "n2_fraction": n2_fraction,
            "D_inner_mm": D_inner * 1000,
            "L_m": L,
            "converged": result["converged"],
            "n_iter": result["n_iter"],
            "T_out_K": T_out,
            "T_wall_out_K": T_wall_out,
            "q_avg_kW_m2": q_avg,
            "Y_C2H6_out": Y_C2H6_out,
            "Y_C2H4_out": Y_C2H4_out,
            "Y_C2H2_out": Y_C2H2_out,
            "Y_H2_out": Y_H2_out,
            "Y_N2_out": Y_N2_out,
            "conversion": conversion,
            "selectivity_C2H4": selectivity_C2H4,
            "selectivity_C2H2": selectivity_C2H2,
            "error": None,
        }
        
    except Exception as e:
        return {
            "T_in_K": T_in,
            "P_in_bar": P_in / 1e5,
            "T_furnace_K": T_furnace,
            "v_z0_m_s": v_z0,
            "n2_fraction": n2_fraction,
            "D_inner_mm": D_inner * 1000,
            "L_m": L,
            "converged": False,
            "n_iter": 0,
            "T_out_K": None,
            "T_wall_out_K": None,
            "q_avg_kW_m2": None,
            "Y_C2H6_out": None,
            "Y_C2H4_out": None,
            "Y_C2H2_out": None,
            "Y_H2_out": None,
            "Y_N2_out": None,
            "conversion": None,
            "selectivity_C2H4": None,
            "selectivity_C2H2": None,
            "error": str(e),
        }


def run_sensitivity_analysis(
    T_in_values: List[float] = None,
    P_in_values: List[float] = None,
    T_furnace_values: List[float] = None,
    v_z0_values: List[float] = None,
    n2_fraction_values: List[float] = None,
    output_csv: str = "sensitivity_results.csv",
    output_json: str = "sensitivity_results.json",
) -> List[Dict[str, Any]]:
    if T_in_values is None:
        T_in_values = [900.0, 950.0, 1000.0, 1050.0]
    if P_in_values is None:
        P_in_values = [5e5, 10e5, 20e5]
    if T_furnace_values is None:
        T_furnace_values = [1300.0, 1400.0, 1500.0]
    if v_z0_values is None:
        v_z0_values = [3.0, 5.0, 10.0]
    if n2_fraction_values is None:
        n2_fraction_values = [0.0, 0.3, 0.5]
    
    param_grid = list(itertools.product(
        T_in_values,
        P_in_values,
        T_furnace_values,
        v_z0_values,
        n2_fraction_values,
    ))
    
    total_cases = len(param_grid)
    print(f"Sensitivity Analysis: {total_cases} cases")
    print("=" * 60)
    print(f"T_in:       {T_in_values} K")
    print(f"P_in:       {[p/1e5 for p in P_in_values]} bar")
    print(f"T_furnace:  {T_furnace_values} K")
    print(f"v_z0:       {v_z0_values} m/s")
    print(f"N2 frac:    {n2_fraction_values}")
    print("=" * 60)
    
    results = []
    start_time = time.time()
    
    for i, (T_in, P_in, T_furnace, v_z0, n2_frac) in enumerate(param_grid):
        elapsed = time.time() - start_time
        if i > 0:
            eta = elapsed / i * (total_cases - i)
            eta_str = f"ETA: {eta/60:.1f} min"
        else:
            eta_str = ""
        
        print(f"[{i+1}/{total_cases}] T_in={T_in:.0f}K, P={P_in/1e5:.0f}bar, "
              f"T_furn={T_furnace:.0f}K, v={v_z0:.0f}m/s, N2={n2_frac*100:.0f}% {eta_str}")
        
        result = run_single_case(
            T_in=T_in,
            P_in=P_in,
            T_furnace=T_furnace,
            v_z0=v_z0,
            n2_fraction=n2_frac,
        )
        results.append(result)
        
        if result["converged"]:
            print(f"    -> Conv={result['conversion']*100:.1f}%, "
                  f"Y_C2H4={result['Y_C2H4_out']*100:.2f}%, "
                  f"Sel_C2H4={result['selectivity_C2H4']*100:.1f}%")
        else:
            print(f"    -> FAILED: {result.get('error', 'did not converge')}")
    
    csv_path = ROOT / output_csv
    fieldnames = list(results[0].keys())
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    print(f"\nResults saved to: {csv_path}")
    
    json_path = ROOT / output_json
    with open(json_path, "w") as f:
        json.dump({
            "metadata": {
                "generated_at": datetime.now().isoformat(),
                "total_cases": total_cases,
                "elapsed_time_s": time.time() - start_time,
                "parameters": {
                    "T_in_K": T_in_values,
                    "P_in_bar": [p/1e5 for p in P_in_values],
                    "T_furnace_K": T_furnace_values,
                    "v_z0_m_s": v_z0_values,
                    "n2_fraction": n2_fraction_values,
                }
            },
            "results": results,
        }, f, indent=2)
    print(f"JSON saved to: {json_path}")
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    converged_results = [r for r in results if r["converged"]]
    failed_results = [r for r in results if not r["converged"]]
    
    print(f"Converged: {len(converged_results)}/{total_cases}")
    print(f"Failed: {len(failed_results)}/{total_cases}")
    
    if converged_results:
        Y_C2H4_values = [r["Y_C2H4_out"] for r in converged_results if r["Y_C2H4_out"] is not None]
        conv_values = [r["conversion"] for r in converged_results if r["conversion"] is not None]
        sel_values = [r["selectivity_C2H4"] for r in converged_results if r["selectivity_C2H4"] is not None]
        
        print(f"\nY_C2H4_out: min={min(Y_C2H4_values)*100:.2f}%, "
              f"max={max(Y_C2H4_values)*100:.2f}%, "
              f"mean={np.mean(Y_C2H4_values)*100:.2f}%")
        print(f"Conversion: min={min(conv_values)*100:.1f}%, "
              f"max={max(conv_values)*100:.1f}%, "
              f"mean={np.mean(conv_values)*100:.1f}%")
        print(f"Selectivity C2H4: min={min(sel_values)*100:.1f}%, "
              f"max={max(sel_values)*100:.1f}%, "
              f"mean={np.mean(sel_values)*100:.1f}%")
        
        yields = [(r, r["conversion"] * r["selectivity_C2H4"]) 
                  for r in converged_results 
                  if r["conversion"] is not None and r["selectivity_C2H4"] is not None]
        if yields:
            best = max(yields, key=lambda x: x[1])
            print(f"\nBest C2H4 yield ({best[1]*100:.1f}%):")
            print(f"  T_in={best[0]['T_in_K']:.0f}K, P={best[0]['P_in_bar']:.0f}bar, "
                  f"T_furn={best[0]['T_furnace_K']:.0f}K, v={best[0]['v_z0_m_s']:.0f}m/s, "
                  f"N2={best[0]['n2_fraction']*100:.0f}%")
    
    return results


if __name__ == "__main__":
    results = run_sensitivity_analysis(
        T_in_values=[950.0, 1000.0, 1100, 1200, 1300],
        P_in_values=[2e5, 5e5, 10e5, 20e5, 40e5],
        T_furnace_values=[1500.0],
        v_z0_values=[5.0, 10.0],
        n2_fraction_values=[0.3],
    )
