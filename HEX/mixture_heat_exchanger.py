#!/usr/bin/env python
"""Mixture heat exchanger using thermo library."""
import numpy as np

# =============================================================================
# INPUT PARAMETERS — edit this section only
# =============================================================================

# --- Hot side inlet ---
HOT_COMPONENTS       = ["water"]
HOT_TEMPERATURE_IN   = 773.15    # K
HOT_PRESSURE_IN      = 101325    # Pa
HOT_FLOW_MOL         = 1000      # mol/s
HOT_COMPOSITION      = [1.0]     # mole fractions (order matches HOT_COMPONENTS)

# --- Cold side inlet ---
COLD_COMPONENTS      = ["benzene"]
COLD_TEMPERATURE_IN  = 323.15    # K
COLD_PRESSURE_IN     = 101325    # Pa
COLD_FLOW_MOL        = 150       # mol/s
COLD_COMPOSITION     = [1.0]     # mole fractions (order matches COLD_COMPONENTS)

# --- Operation mode ---
# "T_hot_out": fix hot outlet temperature; "T_cold_out": fix cold outlet temperature
OPERATION_MODE       = "T_cold_out"
TARGET_VALUE         = 600.15    # K (outlet temperature to fix)

# --- Design (optional; set to None for design mode — area calculated from duty) ---
AREA                = None      # m^2
OVERALL_U            = 800    # W/m^2/K (estimated if None: 800 liq-liq, 200 liq-vap, 30 vap-vap, very rough estimate, do not trust. 
                                # Very strongly recommended to set this to a value you know is correct or taken from studies)

# =============================================================================
# END INPUT PARAMETERS
# =============================================================================

import math
from scipy.optimize import brentq
import warnings
import sys
try:
    from thermo import Chemical
    THERMO_AVAILABLE = True
except ImportError:
    Chemical = None
    THERMO_AVAILABLE = False
    warnings.warn("The 'thermo' library is not installed. Please install it using 'pip install thermo'.")

class MixtureHeatExchanger:
    def __init__(self, hot_components=None, cold_components=None):
        if not THERMO_AVAILABLE:
            raise ImportError("The 'thermo' library is essential for this model to run.")
        print("="*80)
        print("INDUSTRIAL HEAT EXCHANGER MODEL (DATABASE-DRIVEN) - CORRECTED")
        print("="*80)
        self.hot_components = [self._load_component(c) for c in hot_components or ['water']]
        self.cold_components = [self._load_component(c) for c in cold_components or ['toluene']]
        print(f"Hot side: {[c['name'] for c in self.hot_components]}")
        print(f"Cold side: {[c['name'] for c in self.cold_components]}")
        self.area = None
        self.U = None
    
    def _load_component(self, name):
        try:
            chem = Chemical(name)
            if not (chem.MW and math.isfinite(chem.MW)):
                raise ValueError("Component found but has invalid properties.")
            
            Psat_func = None
            try:
                if hasattr(chem, 'Psat') and callable(chem.Psat):
                    Psat_func = chem.Psat
            except:
                pass
            
            return {
                'name': name,
                'MW': chem.MW,
                'Tc': chem.Tc,
                'Pc': chem.Pc,
                'omega': chem.omega,
                'Tb': chem.Tb,
                'Hvap': chem.Hvap(chem.Tb) if callable(chem.Hvap) else chem.Hvap,
                'obj': chem,
                'Psat': Psat_func
            }
        except Exception as e:
            print(f"Error loading component '{name}': {e}", file=sys.stderr)
            raise ValueError(f"Could not load properties for component: {name}")
    
    def _cp_component(self, T, P, comp):
        try:
            is_vapor = False
            if comp['Tc'] and not math.isnan(comp['Tc']):
                if T > comp['Tc']:
                    is_vapor = True
                else:
                    if comp['Psat']:
                        try:
                            Psat = comp['Psat'](T)
                            if Psat is not None and not math.isnan(Psat) and P < Psat:
                                is_vapor = True
                        except:
                            if T > comp['Tb']:
                                is_vapor = True
                    else:
                        if T > comp['Tb']:
                            is_vapor = True
            
            if is_vapor:
                cp_attr = comp['obj'].Cpg
            else:
                cp_attr = comp['obj'].Cpl
                
            if callable(cp_attr):
                cp = cp_attr(T)
            else:
                cp = cp_attr
                
            if cp is not None and math.isfinite(cp):
                return cp
                
        except Exception as e:
            print(f"Warning: Error calculating Cp for {comp['name']} at T={T}K, P={P}Pa: {str(e)}")
        
        if comp['name'].lower() == 'water':
            return 33.6 if T > 373.15 else 75.3
        if T > 500:
            return 30.0
        return 75.0

    def _cp_mixture(self, T, P, z, components):
        cp_mix = sum(zi * self._cp_component(T, P, comp) for zi, comp in zip(z, components))
        return cp_mix
    
    def setup_conditions(self, T_hot_in, T_cold_in, molar_flow_hot, molar_flow_cold, 
                        hot_composition=None, cold_composition=None,
                        P_hot=101325, P_cold=101325,
                        area=None, U=None,
                        mode='T_hot_out', target_value=None):
        self.T_hot_in = T_hot_in
        self.T_cold_in = T_cold_in
        self.molar_flow_hot = molar_flow_hot
        self.molar_flow_cold = molar_flow_cold
        self.P_hot = P_hot
        self.P_cold = P_cold
        self.area = area
        self.U = U
        self.operation_mode = mode
        self.target_value = target_value
        self.hot_composition = np.array(hot_composition or [1.0])
        self.cold_composition = np.array(cold_composition or [1.0])
        print(f"\nHot inlet: {T_hot_in:.1f}K, {P_hot/1000:.1f} kPa")
        print(f"Cold inlet: {T_cold_in:.1f}K, {P_cold/1000:.1f} kPa")
        if self.area is not None and self.U is not None:
            print(f'Heat exchanger: {self.area:.1f} m2, U = {self.U:.0f} W/(m2.K)')
        print(f'Operation mode: {mode} = {target_value}')
    
    def solve(self):
        if self.operation_mode == 'T_hot_out':
            Q = self._calculate_heat_duty('hot', self.T_hot_in, self.target_value, self.P_hot)
            T_hot_out = self.target_value
            T_cold_out = self._solve_outlet_temp('cold', Q, self.P_cold)
        elif self.operation_mode == 'T_cold_out':
            Q = self._calculate_heat_duty('cold', self.T_cold_in, self.target_value, self.P_cold)
            T_cold_out = self.target_value
            T_hot_out = self._solve_outlet_temp('hot', Q, self.P_hot)
        else:
            raise ValueError(f"Unsupported mode: {self.operation_mode}")
        
        if T_cold_out > self.T_hot_in or T_hot_out < self.T_cold_in:
            print("\n--- THERMODYNAMIC VIOLATION WARNING ---")
            print(f"The calculated outlet temperatures (Th_out={T_hot_out:.1f}K, Tc_out={T_cold_out:.1f}K) are physically impossible.")
            print("This requires the cold stream to get hotter than the hot stream's inlet, which violates the 2nd Law of Thermodynamics.")
            print("The process as specified is not feasible.")
            print("-----------------------------------------\n")
        
        if self.area is None or self.U is None:
            self.U = self.U or self._estimate_U()
            lmtd = self._calculate_lmtd(T_hot_out, T_cold_out)
            if lmtd > 0:
                self.area = Q / (self.U * lmtd)
                print(f"DESIGN MODE: Calculated Area = {self.area:.2f} m^2 for U = {self.U} W/m^2.K")
        else:
            lmtd = self._calculate_lmtd(T_hot_out, T_cold_out)
            if lmtd > 0:
                Q_spec = self.U * self.area * lmtd
                if Q_spec < Q:
                    print(f"\n--- EQUIPMENT WARNING: Undersized ---")
                    print(f"Required Duty: {Q/1e3:.1f} kW | Equipment Capacity: {Q_spec/1e3:.1f} kW")
                    print("---------------------------------------\n")
        
        self.print_results(T_hot_out, T_cold_out, Q)
        return True
    
    def _estimate_U(self):
        hot_vapor = any(self._is_vapor(self.T_hot_in, self.P_hot, comp) 
                       for comp in self.hot_components)
        cold_vapor = any(self._is_vapor(self.T_cold_in, self.P_cold, comp) 
                        for comp in self.cold_components)
        
        if hot_vapor and cold_vapor:
            return 30
        elif hot_vapor or cold_vapor:
            return 200
        return 800

    def _is_vapor(self, T, P, comp):
        try:
            if comp['Tc'] and T > comp['Tc']:
                return True
            
            if comp['Psat']:
                Psat = comp['Psat'](T)
                if Psat is not None and P < Psat:
                    return True
            else:
                if T > comp['Tb']:
                    return True
        except:
            pass
        return False
    
    def _calculate_heat_duty(self, side, T_in, T_out, P):
        if side == 'hot':
            components = self.hot_components
            composition = self.hot_composition
            molar_flow = self.molar_flow_hot
        else:
            components = self.cold_components
            composition = self.cold_composition
            molar_flow = self.molar_flow_cold
        
        delta_H = self._calculate_enthalpy_change_rigorous(T_in, T_out, composition, components, P)
        return abs(molar_flow * delta_H)
    
    def _solve_outlet_temp(self, side, Q, P):
        if side == 'hot':
            T_in = self.T_hot_in
            Q_sign = -1
            composition = self.hot_composition
            components = self.hot_components
            molar_flow = self.molar_flow_hot
            T_upper = min(self.T_hot_in + 50, 2000)
            T_lower = max(self.T_cold_in - 50, 200)
        else:
            T_in = self.T_cold_in
            Q_sign = 1
            composition = self.cold_composition
            components = self.cold_components
            molar_flow = self.molar_flow_cold
            T_upper = self.T_hot_in
            T_lower = max(self.T_cold_in - 50, 200)

        def energy_balance(T_out):
            delta_H = self._calculate_enthalpy_change_rigorous(T_in, T_out, composition, components, P)
            if molar_flow == 0:
                return float(-Q)
            return float(molar_flow * delta_H - (Q * Q_sign))

        try:
            return brentq(energy_balance, T_lower, T_upper, xtol=0.1)
        except ValueError:
            eb_lo = energy_balance(T_lower)
            eb_hi = energy_balance(T_upper)
            if side == 'cold' and eb_hi < 0:
                Q_max = molar_flow * self._calculate_enthalpy_change_rigorous(
                    T_in, T_upper, composition, components, P)
                print("\n--- PHYSICS INFEASIBILITY ---")
                print(f"Cold stream cannot absorb required duty {Q/1e6:.1f} MW.")
                print(f"Max absorption (heated to {T_upper:.0f} K): {Q_max/1e6:.1f} MW.")
                print("Remedies: increase cold flow, reduce hot flow, or relax hot outlet temperature.")
                print("--------------------------------\n")
                raise ValueError("Infeasible: cold stream cannot absorb heat duty.") from None
            elif side == 'hot' and eb_lo > 0:
                print("\n--- PHYSICS INFEASIBILITY ---")
                print(f"Hot stream cannot supply required duty {Q/1e6:.1f} MW.")
                print("Remedies: increase hot flow, reduce cold flow, or relax cold outlet temperature.")
                print("--------------------------------\n")
                raise ValueError("Infeasible: hot stream cannot supply heat duty.") from None
            else:
                warnings.warn(f"Outlet temperature solver failed for the {side} side.")
            return T_in
    
    def _calculate_max_heat_transfer(self, P_hot, P_cold):
        Q_hot_potential = self._calculate_heat_duty('hot', self.T_hot_in, self.T_cold_in, P_hot)
        Q_cold_potential = self._calculate_heat_duty('cold', self.T_cold_in, self.T_hot_in, P_cold)
        return min(Q_hot_potential, Q_cold_potential)
    
    def _calculate_enthalpy_change_rigorous(self, T_initial, T_final, z, components, P):
        h_total = 0.0
        if len(components) == 1:
            comp = components[0]
            phase_change_temp = None
            
            try:
                if comp['Psat']:
                    low, high = max(0.5 * comp['Tb'], 1), min(1.5 * comp['Tb'], comp['Tc'])
                    for _ in range(20):
                        mid = (low + high) / 2
                        try:
                            Psat_mid = comp['Psat'](mid)
                            if Psat_mid < P:
                                low = mid
                            else:
                                high = mid
                        except:
                            break
                    phase_change_temp = (low + high) / 2
                else:
                    phase_change_temp = comp['Tb']
            except:
                phase_change_temp = comp['Tb']
            
            if phase_change_temp and not math.isnan(phase_change_temp):
                if T_initial < phase_change_temp and T_final > phase_change_temp:
                    h_liquid = self._mixture_enthalpy_integral(T_initial, phase_change_temp, z, components, P)
                    h_vapor = self._mixture_enthalpy_integral(phase_change_temp, T_final, z, components, P)
                    h_total = h_liquid + comp['Hvap'] + h_vapor
                elif T_initial > phase_change_temp and T_final < phase_change_temp:
                    h_vapor = self._mixture_enthalpy_integral(T_initial, phase_change_temp, z, components, P)
                    h_liquid = self._mixture_enthalpy_integral(phase_change_temp, T_final, z, components, P)
                    h_total = h_vapor - comp['Hvap'] + h_liquid
                else:
                    h_total = self._mixture_enthalpy_integral(T_initial, T_final, z, components, P)
            else:
                h_total = self._mixture_enthalpy_integral(T_initial, T_final, z, components, P)
        else:
            h_total = self._mixture_enthalpy_integral(T_initial, T_final, z, components, P)
            
        return h_total
    
    def _mixture_enthalpy_integral(self, T_ref, T, z, components, P):
        if abs(T - T_ref) < 1e-6:
            return 0.0
            
        num_points = max(20, int(abs(T - T_ref)/5) + 1)
        T_points = np.linspace(T_ref, T, num_points)
        cp_vals = [self._cp_mixture(t, P, z, components) for t in T_points]
        return float((np.trapezoid if hasattr(np, 'trapezoid') else np.trapz)(cp_vals, T_points))
    
    def _calculate_lmtd(self, T_hot_out, T_cold_out):
        delta_T1 = self.T_hot_in - T_cold_out
        delta_T2 = T_hot_out - self.T_cold_in
        
        if delta_T1 > 1e-6 and delta_T2 > 1e-6:
            if abs(delta_T1 - delta_T2) < 1e-6:
                return delta_T1
            return (delta_T1 - delta_T2) / math.log(delta_T1 / delta_T2)
        elif delta_T1 <= 1e-6 or delta_T2 <= 1e-6:
            print("\nWARNING: Temperature crossover detected")
            return 0
        return 0
    
    def print_results(self, T_hot_out, T_cold_out, Q):
        print("\n--- SIMULATION RESULTS ---")
        print(f"Heat Duty: {Q/1000:.2f} kW")
        print(f"Hot Stream Outlet: {T_hot_out:.2f} K ({T_hot_out-273.15:.2f} C)")
        print(f"Cold Stream Outlet: {T_cold_out:.2f} K ({T_cold_out-273.15:.2f} C)")
        
        for i, comp in enumerate(self.hot_components):
            if self._is_vapor(self.T_hot_in, self.P_hot, comp) != self._is_vapor(T_hot_out, self.P_hot, comp):
                phase_in = "vapor" if self._is_vapor(self.T_hot_in, self.P_hot, comp) else "liquid"
                phase_out = "vapor" if self._is_vapor(T_hot_out, self.P_hot, comp) else "liquid"
                print(f"  Hot stream: {comp['name']} changed phase from {phase_in} to {phase_out}")
        
        for i, comp in enumerate(self.cold_components):
            if self._is_vapor(self.T_cold_in, self.P_cold, comp) != self._is_vapor(T_cold_out, self.P_cold, comp):
                phase_in = "vapor" if self._is_vapor(self.T_cold_in, self.P_cold, comp) else "liquid"
                phase_out = "vapor" if self._is_vapor(T_cold_out, self.P_cold, comp) else "liquid"
                print(f"  Cold stream: {comp['name']} changed phase from {phase_in} to {phase_out}")
        
        if self.area:
            print(f"Heat Exchanger Area: {self.area:.2f} m^2")
        print("--------------------------")

def main():
    hex_model = MixtureHeatExchanger(
        hot_components=HOT_COMPONENTS,
        cold_components=COLD_COMPONENTS,
    )
    hex_model.setup_conditions(
        T_hot_in=HOT_TEMPERATURE_IN,
        T_cold_in=COLD_TEMPERATURE_IN,
        molar_flow_hot=HOT_FLOW_MOL,
        molar_flow_cold=COLD_FLOW_MOL,
        hot_composition=HOT_COMPOSITION,
        cold_composition=COLD_COMPOSITION,
        P_hot=HOT_PRESSURE_IN,
        P_cold=COLD_PRESSURE_IN,
        area=AREA,
        U=OVERALL_U,
        mode=OPERATION_MODE,
        target_value=TARGET_VALUE,
    )
    hex_model.solve()

if __name__ == '__main__':
    main()