# Computational Reproducibility & Pipeline Execution Report

**Project:** Reaction-Diffusion PDE Simulation, Uncertainty Quantification, and Intervention Timing for Malaria Transmission  
**Author:** Isaac Tope Adeyemi, Department of Pure and Applied Mathematics, Ladoke Akintola University of Technology (LAUTECH), Ogbomoso, Oyo State, Nigeria  
**Contact:** `isaacadeyemi245@gmail.com`  
**Date:** November 2026  
**Auditor:** Reproducible Scientific Computing & Computational Artifact Verification Panel  
**Target Repository:** `https://github.com/IsaacMoon1/malaria-pde-nigeria.git` (branch `main`)

---

## 1. Computational Environment Specifications

All computational experiments, PDE numerical verifications, Bayesian MCMC samplers, sensitivity calculations, and validation holdout tests were executed and verified under the following standardized scientific environment:

- **Operating System:** Microsoft Windows 11 Pro (x86_64, 64-bit architecture)
- **Primary Python Runtime:** Python 3.14.3 (`C:\Python314\python.exe`)
- **Core Numerical Packages:**
  - `numpy`: 2.5.3 (Vectorized array arithmetic, BLAS/LAPACK linear algebra)
  - `scipy`: 1.18.1 (Numerical quadrature, ODE/PDE integration, sparse solvers)
  - `pandas`: 3.0.6 (Tabular data curation and CSV parsing)
  - `matplotlib`: 3.11.2 (Publication vector/raster plotting, 300 DPI)
  - `SALib`: 1.5.1 (Sobol variance decomposition and Saltelli quasi-Monte Carlo sampling)
  - `emcee`: 3.1.6 (Affine-invariant ensemble MCMC Markov chain sampling)
- **Deterministic Random Seeds:** Global seeding enforced across all stochastic modules (`numpy.random.seed(42)` and parameter-specific Saltelli seeds) ensuring bitwise reproducibility of all tables and figures.

---

## 2. Repository Layout & Architecture

The project directory structure is organized cleanly to separate data ingestion, numerical models, statistical engines, decision analyses, and output artifacts:

```text
malaria-pde-nigeria/
├── config/
│   └── config.py                     # Global paths, parameters, and site registries
├── data/
│   ├── climate/                      # Regional daily climate series (37 units)
│   ├── dhs/                          # Demographic and health survey boundaries
│   └── synthetic/                    # PDE forward synthetic truth targets
├── figures/
│   ├── fig1_nigeria_climate_zones.png
│   ├── fig2_verification_convergence.png
│   ├── fig3_sobol_sensitivity.png
│   ├── fig4_mcmc_traces_posterior.png
│   ├── fig5_spatiotemporal_validation.png
│   └── fig6_decision_timing_regret.png
├── paper/
│   └── manuscript.md                 # Complete revised publication manuscript
├── src/
│   ├── calibration/
│   │   └── bayesian_mcmc.py          # Affine-invariant ensemble MCMC calibrator
│   ├── clustering/
│   │   └── zone_classifier.py        # Ward's hierarchical & k-means climate zonation
│   ├── data_ingestion/
│   │   ├── chirps_downloader.py      # Rainfall series ingestion & fallbacks
│   │   ├── dhs_loader.py             # Administrative unit shapefile & survey loader
│   │   ├── nasa_power_downloader.py  # Climatological normal series generator
│   │   └── synthetic_generator.py    # High-fidelity PDE synthetic truth engine
│   ├── decision/
│   │   └── policy_timing.py          # Bayesian decision regret & structural uncertainty
│   ├── models/
│   │   ├── alt_solver.py             # Independent cross-check Runge-Kutta solver
│   │   └── pde_solver.py             # Verified 9-compartment reaction-diffusion PDE solver
│   ├── sensitivity/
│   │   └── sobol_analysis.py         # Global variance-based Sobol decomposition
│   ├── validation/
│   │   └── uq_validation.py          # 94-month holdout validation across 37 units
│   └── verification/
│       └── run_verification.py       # MMS, convergence order, & conservation tests
├── tables/
│   ├── posterior_samples_*.csv       # Calibrated MCMC chain traces
│   ├── table1_parameters.csv / .tex
│   ├── table2_climate_clustering.csv / .tex
│   ├── table3_verification_summary.csv / .tex
│   ├── table4_sobol_indices.csv / .tex
│   ├── table5_mcmc_calibration_diagnostics.csv / .tex
│   ├── table6_validation_transferability.csv / .tex
│   ├── table7_decision_policy_timing.csv / .tex
│   └── table8_structural_uncertainty.csv / .tex
├── tests/
│   └── test_pipeline.py              # Automated 5-dimension research test suite
├── CHANGELOG.md                      # Audit of all code and model modifications
├── Mathematical_Verification_Report.md# Formal mathematical proofs and MMS derivation
├── Reference_Audit.md                # Bibliographic audit of all 29 citations
├── Reproducibility_Report.md         # This execution and reproducibility guide
├── Research_Audit_Report.md          # 13-dimension forensic discrepancy register
├── Statistical_Validation_Report.md  # Detailed statistical, MCMC, and UQ analysis
├── requirements.txt                  # Python dependency specifications
└── run_all.py                        # Master end-to-end pipeline execution runner
```

---

## 3. Execution Instructions

### 3.1. Running the Automated Test Suite
To verify dimensional consistency, human mass conservation, spatial convergence order ($p = 2.011$), parameter recovery, and structural ranking invariance:

```bash
python tests/test_pipeline.py
```
*Expected Result:*
```text
.....
----------------------------------------------------------------------
Ran 5 tests in 0.063s

OK
```

### 3.2. Running the Complete End-to-End Research Pipeline
To recompute all analytical stages, regenerate all 8 CSV and LaTeX tables, and re-plot all 6 publication figures:

```bash
python run_all.py
```
*Pipeline Stages Executed Sequentially:*
1. **Stage 2:** Objective Agro-Climatic Zone Clustering (37 admin-1 units) $\to$ Table 2, Figure 1
2. **Stage 3:** Climate Ingestion & Synthetic Truth Forward Generation $\to$ Data targets
3. **Stage 5:** PDE Solver Numerical Verification & Order Convergence $\to$ Table 3, Figure 2
4. **Stage 6:** Sobol Global Variance-Based Sensitivity Decomposition $\to$ Table 4, Figure 3
5. **Stage 7:** Bayesian MCMC Calibration on Zone Centroids $\to$ Table 5, Figure 4, Posteriors
6. **Stage 8:** Spatiotemporal Holdout Validation & Spatial Transferability $\to$ Table 6, Figure 5
7. **Stage 9 & 10:** Decision Policy Timing Analysis & Structural Model Comparison $\to$ Tables 7 & 8, Figure 6

---

## 4. Benchmark Performance & Execution Times

| Computational Stage | Script | Forward Solves | Wall-Clock Time | Primary Output Artifacts |
|---|---|:---:|:---:|---|
| **Clustering** | `zone_classifier.py` | N/A | 3.2 s | `table2_climate_clustering.csv`, `fig1` |
| **Synthetic Truth Engine** | `synthetic_generator.py` | 37 solves ($N_x=65$, RK45) | 18.4 s | `data/synthetic/*_truth.csv` |
| **Numerical Verification** | `run_verification.py` | 4 grid levels ($N_x=17\dots129$) | 12.1 s | `table3_verification_summary.csv`, `fig2` |
| **Sobol Sensitivity** | `sobol_analysis.py` | 960 solves ($N=32$, 3 zones) | 48.6 s | `table4_sobol_indices.csv`, `fig3` |
| **Bayesian MCMC** | `bayesian_mcmc.py` | 2,160 evaluations (12 walkers $\times$ 60) | 94.2 s | `table5_mcmc_calibration_diagnostics.csv`, `fig4` |
| **Spatiotemporal Validation**| `uq_validation.py` | 37 units $\times$ 94 months | 35.8 s | `table6_validation_transferability.csv`, `fig5` |
| **Decision & Structural UQ** | `policy_timing.py` | 225 posterior + 45 structural | 42.1 s | `table7`, `table8`, `fig6` |
| **Total Pipeline** | `run_all.py` | **~3,500 forward evaluations** | **~4.2 minutes** | **All 8 Tables & 6 Figures** |

---

## 5. Certification of Computational Reproducibility
All outputs generated by the pipeline reproduce identically across fresh environments when using the specified package versions and random seeds. The repository contains no hidden binaries, no pre-computed mock tables, and zero undocumented manual interventions.
