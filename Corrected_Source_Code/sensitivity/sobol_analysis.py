"""
Stage 6: Global Sensitivity Analysis via Sobol Variance Decomposition (SALib).
Computes first-order (S1) and total-order (ST) Sobol sensitivity indices for:
1. Mean Clinical Prevalence
2. Seasonal Peak Transmission Timing (Day of Year)
3. Cumulative Cases Averted by SMC Intervention
Evaluated across all three representative climate zones (Zone 1 Sahel, Zone 2 Guinea Savanna, Zone 3 Humid Forest).
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
from config.config import FIGURES_DIR, TABLES_DIR, SEED
from src.models.pde_solver import MalariaPDESolver
from src.data_ingestion.nasa_power_downloader import get_or_download_climate

# Define the 8 key epidemiological parameters and their physically plausible prior bounds
SOBOL_PROBLEM = {
    'num_vars': 8,
    'names': ['D_M', 'a_biting', 'b_H', 'b_M', 'gamma_H', 'omega_H', 'K_rain_scale', 'rho_SMC'],
    'bounds': [
        [0.01, 0.25],     # D_M (km^2/day)
        [0.20, 0.50],     # a_biting (bites/day)
        [0.25, 0.80],     # b_H (transmission to human)
        [0.20, 0.70],     # b_M (transmission to mosquito)
        [0.005, 0.030],   # gamma_H (recovery rate day^-1)
        [0.001, 0.010],   # omega_H (immunity waning day^-1)
        [150.0, 650.0],   # K_rain_scale (carrying capacity saturation mm)
        [0.50, 0.95]      # rho_SMC (SMC protection coverage)
    ]
}

def evaluate_model_for_sobol(param_values, climate_df, zone_name):
    """
    Evaluates the fast PDE-ODE simulation for a matrix of parameter sets.
    Returns:
    - Y_prev: Mean human prevalence
    - Y_peak: Peak timing (day of year)
    - Y_averted: Cases averted by SMC
    """
    n_samples = param_values.shape[0]
    Y_prev = np.zeros(n_samples)
    Y_peak = np.zeros(n_samples)
    Y_averted = np.zeros(n_samples)
    
    # 240-day high-transmission seasonal cycle
    t_span = (0.0, 240.0)
    
    # Baseline interventions for SMC (cycles at days 90, 118, 146, 174)
    interventions_smc = {
        "SMC_cycles": [90.0, 118.0, 146.0, 174.0],
        "SMC_coverage": 0.80
    }
    
    # Representative grid resolution for UQ
    solver = MalariaPDESolver(L=5.0, Nx=7)
    
    print(f"[Sobol GSA] Evaluating {n_samples} parameter evaluations for {zone_name}...")
    
    for i in range(n_samples):
        p_row = param_values[i]
        custom_params = {
            "D_M": p_row[0],
            "a_biting_rate": p_row[1],
            "b_H": p_row[2],
            "b_M": p_row[3],
            "gamma_H": p_row[4],
            "omega_H": p_row[5],
            "K_rain_scale": p_row[6],
            "rho_SMC_efficacy": p_row[7]
        }
        solver.p.update(custom_params)
        
        # Run counterfactual without SMC
        sol_no_smc = solver.simulate_fast(t_span, climate_df, interventions=None, dt=0.6)
        # Run scenario with SMC
        sol_smc = solver.simulate_fast(t_span, climate_df, interventions=interventions_smc, dt=0.6)
        
        Nx = solver.Nx
        # Human infectious profile
        IH_no_smc = np.mean(sol_no_smc.y[6*Nx : 7*Nx, :], axis=0)
        IH_smc = np.mean(sol_smc.y[6*Nx : 7*Nx, :], axis=0)
        NH = np.mean(solver.N_H)
        
        # Prevalence time series
        prev_series = IH_smc / max(1.0, NH)
        
        Y_prev[i] = np.mean(prev_series)
        # Peak timing (day of max prevalence)
        peak_idx = np.argmax(prev_series)
        Y_peak[i] = sol_smc.t[peak_idx] if len(sol_smc.t) > peak_idx else 220.0
        # Cases averted
        Y_averted[i] = np.maximum(0.0, trapezoid(IH_no_smc - IH_smc, sol_no_smc.t))
        
    return Y_prev, Y_peak, Y_averted

def run_sobol_sensitivity_analysis():
    np.random.seed(SEED)
    # Generate Sobol quasi-random sequence (Saltelli / Sobol scheme)
    # N=16 -> 16 * (8 + 2) = 160 evaluations per zone
    N = 32
    param_values = sobol_sample.sample(SOBOL_PROBLEM, N, calc_second_order=False, seed=SEED)
    
    zones = {
        "Zone 1 (Sahel)": "Jigawa",
        "Zone 2 (Guinea Savanna)": "FCT",
        "Zone 3 (Humid Forest)": "Abia"
    }
    
    all_results = []
    
    for zone_label, rep_site in zones.items():
        clim = get_or_download_climate(rep_site)
        Y_prev, Y_peak, Y_averted = evaluate_model_for_sobol(param_values, clim, zone_label)
        
        # Analyze Sobol indices
        Si_prev = sobol_analyze.analyze(SOBOL_PROBLEM, Y_prev, calc_second_order=False)
        Si_peak = sobol_analyze.analyze(SOBOL_PROBLEM, Y_peak, calc_second_order=False)
        Si_avert = sobol_analyze.analyze(SOBOL_PROBLEM, Y_averted, calc_second_order=False)
        
        for idx, p_name in enumerate(SOBOL_PROBLEM['names']):
            all_results.append({
                "Zone": zone_label,
                "Parameter": p_name,
                "S1_Prevalence": round(max(0.0, float(Si_prev['S1'][idx])), 4),
                "ST_Prevalence": round(max(0.0, float(Si_prev['ST'][idx])), 4),
                "S1_PeakTiming": round(max(0.0, float(Si_peak['S1'][idx])), 4),
                "ST_PeakTiming": round(max(0.0, float(Si_peak['ST'][idx])), 4),
                "S1_CasesAverted": round(max(0.0, float(Si_avert['S1'][idx])), 4),
                "ST_CasesAverted": round(max(0.0, float(Si_avert['ST'][idx])), 4)
            })
            
    df_sobol = pd.DataFrame(all_results)
    csv_path = os.path.join(TABLES_DIR, "table4_sobol_indices.csv")
    df_sobol.to_csv(csv_path, index=False)
    
    # LaTeX export
    tex_path = os.path.join(TABLES_DIR, "table4_sobol_indices.tex")
    with open(tex_path, "w") as f:
        f.write("% Global Sensitivity Analysis: Total Sobol Indices (ST)\n")
        f.write("\\begin{tabular}{lcccccc}\n\\hline\n")
        f.write("Parameter & Zone 1 Prev & Zone 1 Avert & Zone 2 Prev & Zone 2 Avert & Zone 3 Prev & Zone 3 Avert \\\\\n\\hline\n")
        for p_name in SOBOL_PROBLEM['names']:
            z1 = df_sobol[(df_sobol['Zone'].str.contains('Zone 1')) & (df_sobol['Parameter'] == p_name)].iloc[0]
            z2 = df_sobol[(df_sobol['Zone'].str.contains('Zone 2')) & (df_sobol['Parameter'] == p_name)].iloc[0]
            z3 = df_sobol[(df_sobol['Zone'].str.contains('Zone 3')) & (df_sobol['Parameter'] == p_name)].iloc[0]
            f.write(f"{p_name} & {z1['ST_Prevalence']:.3f} & {z1['ST_CasesAverted']:.3f} & {z2['ST_Prevalence']:.3f} & {z2['ST_CasesAverted']:.3f} & {z3['ST_Prevalence']:.3f} & {z3['ST_CasesAverted']:.3f} \\\\\n")
        f.write("\\hline\n\\end{tabular}\n")
        
    # Generate Publication Figure
    plt.rcParams.update({"font.size": 9.5, "font.family": "serif"})
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=300, sharey=True)
    
    outcomes = [
        ("ST_Prevalence", "Total Sobol Index (ST): Mean Prevalence"),
        ("ST_PeakTiming", "Total Sobol Index (ST): Peak Timing"),
        ("ST_CasesAverted", "Total Sobol Index (ST): Cases Averted")
    ]
    
    zone_colors = {"Zone 1 (Sahel)": "#d95f02", "Zone 2 (Guinea Savanna)": "#7570b3", "Zone 3 (Humid Forest)": "#1b9e77"}
    x_indices = np.arange(len(SOBOL_PROBLEM['names']))
    bar_width = 0.25
    
    for ax_idx, (col_metric, title_str) in enumerate(outcomes):
        for z_i, (zone_label, _) in enumerate(zones.items()):
            sub = df_sobol[df_sobol['Zone'] == zone_label]
            vals = sub[col_metric].values
            axes[ax_idx].bar(x_indices + (z_i - 1) * bar_width, vals, width=bar_width, label=zone_label, color=zone_colors[zone_label], edgecolor="black", linewidth=0.6)
            
        axes[ax_idx].set_title(title_str, fontweight="bold", fontsize=10)
        axes[ax_idx].set_xticks(x_indices)
        axes[ax_idx].set_xticklabels(SOBOL_PROBLEM['names'], rotation=45, ha="right")
        axes[ax_idx].set_ylim(0, 1.0)
        axes[ax_idx].grid(True, linestyle="--", alpha=0.5, axis="y")
        if ax_idx == 0:
            axes[ax_idx].set_ylabel("Total-Order Sobol Index (ST)")
            axes[ax_idx].legend(loc="upper right", fontsize=8)
            
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "fig3_sobol_sensitivity.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
    
    print(f"\n[Sobol GSA] Complete! Saved to:\n  Table: {csv_path}\n  Figure: {fig_path}")

if __name__ == "__main__":
    run_sobol_sensitivity_analysis()
