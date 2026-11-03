#!/usr/bin/env python3
"""
================================================================================
END-TO-END RESEARCH PIPELINE RUNNER
Title: Statistical validation and uncertainty quantification of PDE simulation
       for malaria intervention timing across three representative climate zones in Nigeria
Target Journal: Infectious Disease Modelling / Epidemics / Mathematical Biosciences (Elsevier)
================================================================================
This script executes the complete research pipeline sequentially:
  Stage 2: Objective Agro-Climatic Zone Clustering (37 units)
  Stage 3: Climate Forcing Ingestion & Synthetic Truth Data Engine
  Stage 4 & 5: Numerical Verification & Convergence Analysis (6 tests)
  Stage 6: Global Sensitivity Analysis via Sobol Decomposition
  Stage 7: Bayesian MCMC Calibration on Zone Centroids
  Stage 8: Spatiotemporal Validation & Spatial Transferability (37 units)
  Stage 9 & 10: Decision Analysis, Optimal Timing & Structural Uncertainty

Usage:
    python run_all.py [--fast] [--synthetic-only]
================================================================================
"""

import os
import sys
import time
import argparse
import subprocess

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PROJECT_ROOT)

def run_stage(stage_num: int, stage_name: str, script_relpath: str, extra_args=None):
    print("=" * 80)
    print(f"[*] EXECUTING STAGE {stage_num}: {stage_name}")
    print(f"[*] Script: {script_relpath}")
    print("=" * 80)
    
    script_full_path = os.path.join(PROJECT_ROOT, script_relpath)
    cmd = [sys.executable, script_full_path]
    if extra_args:
        cmd.extend(extra_args)
        
    t0 = time.time()
    res = subprocess.run(cmd, cwd=PROJECT_ROOT)
    elapsed = time.time() - t0
    
    if res.returncode != 0:
        print(f"[!] ERROR: Stage {stage_num} failed with return code {res.returncode}")
        sys.exit(res.returncode)
    else:
        print(f"[✓] Stage {stage_num} completed successfully in {elapsed:.2f} seconds.\n")

def main():
    parser = argparse.ArgumentParser(description="End-to-end execution of the malaria PDE pipeline.")
    parser.add_argument("--fast", action="store_true", help="Run fast verification and sampling routines where supported.")
    args = parser.parse_args()

    print("\n" + "=" * 80)
    print(" NIGERIA MALARIA PDE SIMULATION, UQ & VALIDATION PIPELINE ")
    print(" Target: Elsevier (Infectious Disease Modelling / Epidemics)")
    print("=" * 80 + "\n")
    
    start_total = time.time()
    
    # Stage 2: Clustering
    run_stage(2, "Agro-Climatic Zone Clustering (37 admin-1 units)", 
              os.path.join("src", "clustering", "zone_classifier.py"))
    
    # Stage 3: Synthetic Truth & Climate Ingestion
    run_stage(3, "Climate Ingestion & Synthetic Truth Generator", 
              os.path.join("src", "data_ingestion", "synthetic_generator.py"))
    
    # Stage 4 & 5: Numerical Verification
    run_stage(5, "PDE Solver Numerical Verification & Order Convergence", 
              os.path.join("src", "verification", "run_verification.py"))
    
    # Stage 6: Global Sensitivity Analysis
    run_stage(6, "Sobol Global Variance-Based Sensitivity Decomposition", 
              os.path.join("src", "sensitivity", "sobol_analysis.py"))
    
    # Stage 7: Bayesian MCMC Calibration
    run_stage(7, "Bayesian MCMC Calibration on Zone Centroid Sites", 
              os.path.join("src", "calibration", "bayesian_mcmc.py"))
    
    # Stage 8: Spatiotemporal Holdout & Transferability
    run_stage(8, "Spatiotemporal Holdout Validation & Spatial Transferability", 
              os.path.join("src", "validation", "uq_validation.py"))
    
    # Stage 9 & 10: Decision Analysis & Structural Uncertainty
    run_stage(9, "Decision Policy Timing Analysis & Structural Model Comparison", 
              os.path.join("src", "decision", "policy_timing.py"))

    total_time = time.time() - start_total
    print("=" * 80)
    print(f"[✓] ALL STAGES COMPLETED SUCCESSFULLY IN {total_time/60:.2f} MINUTES!")
    print(f"[✓] Tables written to: {os.path.join(PROJECT_ROOT, 'tables')}")
    print(f"[✓] Figures written to: {os.path.join(PROJECT_ROOT, 'figures')}")
    print("=" * 80 + "\n")

if __name__ == "__main__":
    main()
