"""
Stage 4 & 5: Independent Cross-Check Solver.
Implements an independent 2nd-order Runge-Kutta (Explicit Midpoint / RK2) finite-difference solver.
Used for verification and cross-solver concordance against the primary adaptive SciPy MOL engine.
Spatial order of accuracy: O(dx^2). Temporal order: O(dt^2).
"""

import os
import sys
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
from src.models.pde_solver import MalariaPDESolver

class IndependentPDESolver:
    def __init__(self, L=5.0, Nx=21, dt=0.01, params=None):
        self.L = float(L)
        self.Nx = int(Nx)
        self.dt = float(dt)
        self.dx = self.L / (self.Nx - 1)
        self.x = np.linspace(0, self.L, self.Nx)
        self.mol_solver = MalariaPDESolver(L=L, Nx=Nx, params=params)
        self.p = self.mol_solver.p

    def solve(self, t_span, climate_df, interventions=None, y0=None):
        """
        Integrates forward in time using independent fixed-step 2nd-Order Midpoint RK2 Scheme:
        y_{n+1/2} = y_n + (dt/2) * f(t_n, y_n)
        y_{n+1} = y_n + dt * f(t_n + dt/2, y_{n+1/2})
        """
        t0, tf = t_span
        n_steps = int(np.ceil((tf - t0) / self.dt))
        t_eval = np.linspace(t0, tf, n_steps + 1)
        actual_dt = (tf - t0) / n_steps
        self.dt = actual_dt
        
        Nx = self.Nx
        if y0 is None:
            L0 = np.ones(Nx) * 500.0 * self.mol_solver.breeding_mask
            MS0 = np.ones(Nx) * 300.0
            ME0 = np.ones(Nx) * 20.0
            MI0 = np.ones(Nx) * 10.0
            SH0 = self.mol_solver.N_H * 0.70
            EH0 = self.mol_solver.N_H * 0.05
            IH0 = self.mol_solver.N_H * 0.15
            RH0 = self.mol_solver.N_H * 0.10
            PH0 = np.zeros(Nx)
            y = np.concatenate([L0, MS0, ME0, MI0, SH0, EH0, IH0, RH0, PH0])
        else:
            y = y0.copy()
            
        temp_arr = climate_df["temp_mean"].values
        rain_arr = climate_df["rainfall"].values
        
        def climate_interpolator(t_val):
            idx = int(np.clip(t_val, 0, len(temp_arr) - 1))
            return temp_arr[idx], rain_arr[idx]
            
        y_history = [y.copy()]
        
        for step in range(n_steps):
            t_curr = t_eval[step]
            
            # Predictor (Euler half-step)
            k1 = self.mol_solver.rhs(t_curr, y, climate_interpolator, interventions)
            y_mid = np.maximum(0.0, y + 0.5 * self.dt * k1)
            
            # Corrector (Midpoint full step)
            k2 = self.mol_solver.rhs(t_curr + 0.5 * self.dt, y_mid, climate_interpolator, interventions)
            y = np.maximum(0.0, y + self.dt * k2)
            
            y_history.append(y.copy())
            
        return t_eval, np.array(y_history)
