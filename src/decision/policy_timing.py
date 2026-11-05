"""
Stage 9 & 10: Decision Analysis (Optimal Intervention Timing) & Structural Uncertainty.
1. Evaluates posterior ensemble over SMC intervention start dates (May 15 to Sept 15).
2. Computes:
   - Probability that each start window is optimal P(tau* = tau_k)
   - Expected loss / regret from mistiming relative to empirical practice (July/August)
3. Structural Uncertainty: Re-evaluates timing and disease burden across 3 alternative
   larval carrying capacity formulations:
   - Form A: Michael & Parham Saturation
   - Form B: Exponential Saturating
   - Form C: Piecewise Linear Capped
Exports:
- tables/table7_decision_policy_timing.csv
- tables/table8_structural_uncertainty.csv
- figures/fig6_decision_timing_regret.png
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
from config.config import FIGURES_DIR, TABLES_DIR
from src.models.pde_solver import MalariaPDESolver
from src.data_ingestion.nasa_power_downloader import get_or_download_climate

# Candidate SMC Start Dates (Day of Year)
CANDIDATE_TIMINGS = {
    "May 15 (DOY 135)": 135.0,
    "June 15 (DOY 166)": 166.0,
    "July 15 (DOY 196) [Empirical North]": 196.0,
    "August 15 (DOY 227) [Empirical Middle]": 227.0,
    "September 15 (DOY 258)": 258.0
}

STRUCTURAL_FORMS = ["standard", "exponential", "linear_capped"]

def run_decision_and_structural_analysis():
    print("\n[Decision & Structural Analysis] Running policy timing evaluations...")
    
    # Load calibrated parameters
    df_calib = pd.read_csv(os.path.join(TABLES_DIR, "table5_mcmc_calibration_diagnostics.csv"))
    
    zones = {
        "Zone 1 (Sahel)": "Jigawa",
        "Zone 2 (Guinea Savanna)": "FCT",
        "Zone 3 (Humid Forest)": "Abia"
    }
    
    decision_records = []
    structural_records = []
    
    solver = MalariaPDESolver(L=5.0, Nx=7)
    
    # Timing evaluation results dictionary for plotting
    timing_curves = {}
    
    for zone_label, rep_site in zones.items():
        clim = get_or_download_climate(rep_site)
        clim["date"] = pd.to_datetime(clim["date"])
        # Use representative evaluation year (e.g., 2021)
        clim_eval = clim[(clim["date"] >= "2021-01-01") & (clim["date"] <= "2021-12-31")].reset_index(drop=True)
        
        # Propagate posterior uncertainty: load posterior samples
        post_file = os.path.join(TABLES_DIR, f"posterior_samples_{rep_site.lower()}.csv")
        if os.path.exists(post_file):
            df_post = pd.read_csv(post_file)
            sample_indices = np.linspace(0, len(df_post) - 1, min(25, len(df_post))).astype(int)
            post_draws = df_post.iloc[sample_indices].to_dict('records')
        else:
            sub_cal = df_calib[df_calib["Zone"] == zone_label]
            p_mean = {r["Parameter"]: r["Posterior_Mean"] for _, r in sub_cal.iterrows()}
            post_draws = [p_mean]
            
        t_keys = list(CANDIDATE_TIMINGS.keys())
        averted_ensemble = np.zeros((len(post_draws), len(t_keys)))
        
        for s_idx, theta_draw in enumerate(post_draws):
            solver.p.update(theta_draw)
            sol_b = solver.simulate_fast((0.0, 365.0), clim_eval, interventions=None, dt=1.0)
            gamma_val = solver.p.get("gamma_H", 0.0167)
            # Clinical episodes = integral of gamma_H * I_H(x, t) dx dt
            cases_b = np.sum(sol_b.y[6*solver.Nx : 7*solver.Nx, :]) * gamma_val * solver.dx * 1.0
            
            for k_idx, (t_label, t_start) in enumerate(CANDIDATE_TIMINGS.items()):
                interventions = {
                    "SMC_cycles": [t_start, t_start + 28.0, t_start + 56.0, t_start + 84.0],
                    "SMC_coverage": 0.80
                }
                sol_i = solver.simulate_fast((0.0, 365.0), clim_eval, interventions=interventions, dt=1.0)
                cases_i = np.sum(sol_i.y[6*solver.Nx : 7*solver.Nx, :]) * gamma_val * solver.dx * 1.0
                averted_ensemble[s_idx, k_idx] = max(0.0, cases_b - cases_i)
                
        # Posterior probability that each window is optimal: argmax per posterior sample
        opt_indices_per_sample = np.argmax(averted_ensemble, axis=1)
        prob_optimal = np.array([np.mean(opt_indices_per_sample == k) for k in range(len(t_keys))])
        
        # Expected cases averted across posterior ensemble
        mean_averted = np.mean(averted_ensemble, axis=0)
        
        # Expected Regret / Loss = mean over posterior draws of (max_averted_draw - chosen_averted_draw)
        max_per_sample = np.max(averted_ensemble, axis=1, keepdims=True)
        regret_ensemble = max_per_sample - averted_ensemble
        expected_regret = np.mean(regret_ensemble, axis=0)
        
        # Empirical practice comparison index
        emp_idx = 2 if "Zone 1" in zone_label else 3
        loss_vs_empirical = np.mean(averted_ensemble[:, emp_idx:emp_idx+1] - averted_ensemble, axis=0)
        
        timing_curves[zone_label] = {t_keys[k]: mean_averted[k] for k in range(len(t_keys))}
        best_k = np.argmax(prob_optimal)
        
        for k, t_label in enumerate(t_keys):
            decision_records.append({
                "Zone": zone_label,
                "Representative_Site": rep_site,
                "Timing_Window": t_label,
                "Cases_Averted_Index": round(mean_averted[k], 1),
                "Prob_Optimal": round(prob_optimal[k], 3),
                "Expected_Loss_Regret": round(expected_regret[k], 1),
                "Loss_vs_Empirical_Practice": round(loss_vs_empirical[k], 1),
                "Is_Optimal": (k == best_k)
            })
            
        # Structural Uncertainty Check across all 5 candidate start dates
        for s_form in STRUCTURAL_FORMS:
            sol_base_struct = solver.simulate_fast((0.0, 365.0), clim_eval, interventions=None, carrying_capacity_form=s_form, dt=0.5)
            cases_base_struct = np.sum(sol_base_struct.y[6*solver.Nx : 7*solver.Nx, :])
            for t_label_full, t_start in CANDIDATE_TIMINGS.items():
                t_short = t_label_full.split(' (')[0]
                interventions = {
                    "SMC_cycles": [t_start, t_start + 28.0, t_start + 56.0, t_start + 84.0],
                    "SMC_coverage": 0.80
                }
                sol_struct = solver.simulate_fast((0.0, 365.0), clim_eval, interventions=interventions, carrying_capacity_form=s_form, dt=0.5)
                cases_struct = np.sum(sol_struct.y[6*solver.Nx : 7*solver.Nx, :])
                av_struct = max(0.0, cases_base_struct - cases_struct)
                
                structural_records.append({
                    "Zone": zone_label,
                    "Carrying_Capacity_Form": s_form,
                    "Intervention_Timing": t_short,
                    "Cases_Averted": round(av_struct, 1)
                })
                
    df_dec = pd.DataFrame(decision_records)
    csv_dec = os.path.join(TABLES_DIR, "table7_decision_policy_timing.csv")
    df_dec.to_csv(csv_dec, index=False)
    
    df_struct = pd.DataFrame(structural_records)
    csv_struct = os.path.join(TABLES_DIR, "table8_structural_uncertainty.csv")
    df_struct.to_csv(csv_struct, index=False)
    
    # LaTeX Tables Export
    tex_dec = os.path.join(TABLES_DIR, "table7_decision_policy_timing.tex")
    with open(tex_dec, "w") as f:
        f.write("% Decision Analysis: Intervention Timing Optimization & Expected Loss\n")
        f.write("\\begin{tabular}{llcccc}\n\\hline\n")
        f.write("Zone & Start Window & Averted Burden & $P(\\tau^* = \\tau_k)$ & Regret (Loss) & Relative to Current \\\\\n\\hline\n")
        for _, r in df_dec.iterrows():
            f.write(f"{r['Zone']} & {r['Timing_Window'].split('[')[0].strip()} & {r['Cases_Averted_Index']} & {r['Prob_Optimal']:.2f} & {r['Expected_Loss_Regret']} & {r['Loss_vs_Empirical_Practice']} \\\\\n")
        f.write("\\hline\n\\end{tabular}\n")
        
    tex_struct = os.path.join(TABLES_DIR, "table8_structural_uncertainty.tex")
    with open(tex_struct, "w") as f:
        f.write("% Structural Uncertainty: Model Form Sensitivity Across All Candidate Windows\n")
        f.write("\\begin{tabular}{llccccc}\n\\hline\n")
        f.write("Zone & Habitat Formulation & May 15 & June 15 & July 15 & August 15 & September 15 \\\\\n\\hline\n")
        short_dates = ["May 15", "June 15", "July 15", "August 15", "September 15"]
        for z_name in zones.keys():
            for s_form in STRUCTURAL_FORMS:
                row_vals = []
                for d_str in short_dates:
                    v = df_struct[(df_struct["Zone"] == z_name) & (df_struct["Carrying_Capacity_Form"] == s_form) & (df_struct["Intervention_Timing"] == d_str)]["Cases_Averted"].values[0]
                    row_vals.append(f"{v:.1f}")
                f.write(f"{z_name} & {s_form.capitalize()} & {' & '.join(row_vals)} \\\\\n")
        f.write("\\hline\n\\end{tabular}\n")
        
    # Generate Publication Figure 6 (Timing Curve & Expected Loss Comparison)
    plt.rcParams.update({"font.size": 9.5, "font.family": "serif"})
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), dpi=300)
    
    timing_labels_short = ["May 15", "Jun 15", "Jul 15", "Aug 15", "Sep 15"]
    zone_colors = {"Zone 1 (Sahel)": "#d95f02", "Zone 2 (Guinea Savanna)": "#7570b3", "Zone 3 (Humid Forest)": "#1b9e77"}
    
    for z_name, av_dict in timing_curves.items():
        vals = list(av_dict.values())
        norm_vals = np.array(vals) / max(1.0, max(vals))
        axes[0].plot(timing_labels_short, norm_vals, 'o-', color=zone_colors[z_name], linewidth=2.0, markersize=7, label=z_name)
        # Optimal point
        opt_idx = np.argmax(norm_vals)
        axes[0].scatter(timing_labels_short[opt_idx], 1.0, color="gold", s=160, marker="*", edgecolors="black", linewidth=1.2, zorder=5)
        
    axes[0].set_title("(a) Normalized Intervention Efficacy vs Start Date", fontweight="bold")
    axes[0].set_ylabel("Relative Burden Averted (Normalized)")
    axes[0].set_xlabel("SMC 1st Cycle Start Window")
    axes[0].grid(True, linestyle="--", alpha=0.5)
    axes[0].legend(loc="lower left", fontsize=8.5)
    
    # Plot B: Expected Regret / Loss from Mistiming by Zone
    x_pos = np.arange(len(timing_labels_short))
    b_width = 0.25
    for z_i, z_name in enumerate(zones.keys()):
        sub_d = df_dec[df_dec["Zone"] == z_name]
        regrets = sub_d["Expected_Loss_Regret"].values
        axes[1].bar(x_pos + (z_i - 1)*b_width, regrets, width=b_width, color=zone_colors[z_name], label=z_name, edgecolor="black", linewidth=0.6)
        
    axes[1].set_title("(b) Expected Regret (Loss) from Suboptimal Timing", fontweight="bold")
    axes[1].set_ylabel("Expected Loss (Averted Burden Deficit)")
    axes[1].set_xlabel("SMC Start Date")
    axes[1].set_xticks(x_pos)
    axes[1].set_xticklabels(timing_labels_short)
    axes[1].grid(True, linestyle="--", alpha=0.5, axis="y")
    axes[1].legend(loc="upper right", fontsize=8.5)
    
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "fig6_decision_timing_regret.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
    
    print(f"\n[Decision & Structural Suite] Complete! Saved:\n  Tables: {csv_dec}, {csv_struct}\n  Figure: {fig_path}")

if __name__ == "__main__":
    run_decision_and_structural_analysis()
