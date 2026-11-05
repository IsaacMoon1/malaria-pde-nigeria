"""
Automated Research Verification Test Suite
Covers:
1. Dimensional consistency, non-negativity, and human mass conservation.
2. Numerical convergence order verification (MMS).
3. Calibration diagnostics and parameter recovery validity.
4. Sobol sensitivity index bounds (0 <= S_i, S_Ti <= 1.05).
5. Decision-theoretic ranking invariance across structural forms.
"""

import os
import sys
import unittest
import numpy as np
import pandas as pd

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, PROJECT_ROOT)

from src.models.pde_solver import MalariaPDESolver
from config.config import TABLES_DIR

class TestMalariaPDEPipeline(unittest.TestCase):

    def setUp(self):
        self.solver = MalariaPDESolver(L=5.0, Nx=7)
        # Mock climate array for 30 days
        dates = pd.date_range("2021-01-01", periods=30, freq="D")
        self.mock_clim = pd.DataFrame({
            "date": dates,
            "temp_mean": np.full(30, 26.5),
            "rainfall": np.full(30, 10.0)
        })

    def test_mass_conservation_and_non_negativity(self):
        """Verify human population is conserved to machine precision and all states remain non-negative."""
        sol = self.solver.simulate_fast((0.0, 30.0), self.mock_clim, dt=0.5)
        # Check non-negativity
        self.assertTrue(np.all(sol.y >= -1e-12), "Negative densities detected in solution!")
        
        # Check human mass conservation across all spatial nodes and time points
        Nx = self.solver.Nx
        S_H = sol.y[4*Nx : 5*Nx, :]
        E_H = sol.y[5*Nx : 6*Nx, :]
        I_H = sol.y[6*Nx : 7*Nx, :]
        R_H = sol.y[7*Nx : 8*Nx, :]
        P_H = sol.y[8*Nx : 9*Nx, :]
        N_H = S_H + E_H + I_H + R_H + P_H
        
        # Max deviation from initial N_H (which is 1000.0 per node)
        expected_NH = self.solver.N_H[:, None]
        max_deviation = np.max(np.abs(N_H - expected_NH))
        self.assertLess(max_deviation, 1e-10, f"Human mass conservation violated! Max leak = {max_deviation}")

    def test_convergence_order(self):
        """Verify spatial finite difference achieves theoretical 2nd order convergence (p >= 1.9)."""
        table3_path = os.path.join(TABLES_DIR, "table3_verification_summary.csv")
        self.assertTrue(os.path.exists(table3_path), "Table 3 verification summary missing!")
        df = pd.read_csv(table3_path)
        row = df[df["Metric"].str.contains("Asymptotic Order", case=False, na=False)]
        p_val = float(row["Value"].values[0])
        self.assertGreaterEqual(p_val, 1.90, f"Observed convergence order {p_val} is below second-order threshold!")

    def test_mcmc_diagnostics_and_recovery(self):
        """Verify parameter recovery success rate >= 80% (at least 10 of 12 parameters within 95% HPDI)."""
        table5_path = os.path.join(TABLES_DIR, "table5_mcmc_calibration_diagnostics.csv")
        self.assertTrue(os.path.exists(table5_path), "Table 5 calibration diagnostics missing!")
        df = pd.read_csv(table5_path)
        
        # Check R_hat values are strictly finite and non-negative
        r_hats = df["Gelman_Rubin_Rhat"].values
        for r_hat in r_hats:
            self.assertTrue(np.isfinite(r_hat) and r_hat >= 1.0, f"Invalid R_hat: {r_hat}")
            
        # Check parameter recovery rate
        recoveries = df["Recovery_Success"].values
        success_count = np.sum(recoveries == True)
        self.assertGreaterEqual(success_count, 10, f"Recovery rate too low: {success_count}/12 recovered")

    def test_sobol_index_properties(self):
        """Verify Sobol sensitivity indices obey mathematical bounds 0 <= S_i, S_Ti <= 1.05."""
        table4_path = os.path.join(TABLES_DIR, "table4_sobol_indices.csv")
        self.assertTrue(os.path.exists(table4_path), "Table 4 Sobol indices missing!")
        df = pd.read_csv(table4_path)
        for _, row in df.iterrows():
            for target in ["Prevalence", "PeakTiming", "CasesAverted"]:
                s1 = float(row[f"S1_{target}"])
                st = float(row[f"ST_{target}"])
                self.assertGreaterEqual(s1, -0.05, f"First order index significantly negative: {s1}")
                self.assertLessEqual(s1, 1.05, f"First order index {s1} > 1.05")
                self.assertGreaterEqual(st, -0.05, f"Total order index significantly negative: {st}")
                self.assertLessEqual(st, 1.05, f"Total order index {st} > 1.05")

    def test_decision_structural_ranking_invariance(self):
        """Verify that August 15 achieves maximal cases averted across all 3 structural forms."""
        table8_path = os.path.join(TABLES_DIR, "table8_structural_uncertainty.csv")
        self.assertTrue(os.path.exists(table8_path), "Table 8 structural uncertainty missing!")
        df = pd.read_csv(table8_path)
        for zone in df["Zone"].unique():
            for form in df["Carrying_Capacity_Form"].unique():
                sub = df[(df["Zone"] == zone) & (df["Carrying_Capacity_Form"] == form)]
                best_timing = sub.sort_values(by="Cases_Averted", ascending=False).iloc[0]["Intervention_Timing"]
                self.assertEqual(best_timing, "August 15", 
                                 f"Ranking invariance violated for {zone} under {form}: optimal is {best_timing}")

if __name__ == "__main__":
    unittest.main()
