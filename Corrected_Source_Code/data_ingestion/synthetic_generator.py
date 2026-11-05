"""
Stage 3: Synthetic Ground-Truth Generator.
Generates synthetic observations (monthly clinical incidence and survey prevalence)
for all 37 Nigerian administrative units from 2000-01-01 to 2026-10-31 using
an independent numerical simulation with explicit, known true parameters.
Every output is stamped with 'SYNTHETIC_TRUTH'.
"""

import os
import sys
import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
from config.config import DATA_SYNTHETIC_DIR, NIGERIA_ADMIN1, START_YEAR, SEED
from src.data_ingestion.nasa_power_downloader import get_or_download_climate
from src.clustering.zone_classifier import CLIMATE_NORMALS_DATA

np.random.seed(SEED)

# Known ground-truth parameters for synthetic generation
TRUE_SYNTHETIC_PARAMETERS = {
    "D_M": 0.085,             # Mosquito diffusion (km^2/day)
    "a_biting_rate": 0.33,    # Bites / mosquito / day
    "b_H": 0.50,              # Human infection prob
    "b_M": 0.45,              # Mosquito infection prob
    "gamma_H": 0.0167,        # Human recovery (1/60 days)
    "omega_H": 0.0027,        # Immunity loss (1/365 days)
    "K_rain_scale": 450.0,    # Larval habitat carrying capacity scaling
    "larvicide_eff": 0.75,
    "smc_cycle_efficacy": 0.88
}

def generate_synthetic_truth_dataset(end_date="2026-10-31"):
    """
    Generate complete monthly synthetic malaria burden records across all 37 units.
    """
    all_records = []
    
    # Process each state
    for item in CLIMATE_NORMALS_DATA:
        state_name = item["state"]
        climate_df = get_or_download_climate(state_name)
        climate_df["date"] = pd.to_datetime(climate_df["date"])
        climate_df = climate_df[climate_df["date"] <= pd.to_datetime(end_date)].copy()
        
        # Resample climate to monthly
        climate_df["year_month"] = climate_df["date"].dt.to_period("M")
        monthly_clim = climate_df.groupby("year_month").agg({
            "temp_mean": "mean",
            "rainfall": "sum",
            "humidity": "mean"
        }).reset_index()
        
        # Transmission suitability function (temperature-dependent Briére & degree days)
        # T_opt ~ 27-28C, zero below 16C, zero above 35C
        temps = monthly_clim["temp_mean"].values
        rain = monthly_clim["rainfall"].values
        
        # Mordecai et al. (2013) parabolic/thermal suitability factor
        t_suit = np.maximum(0.0, (temps - 16.0) * (35.0 - temps) / ((28.0 - 16.0) * (35.0 - 28.0)))
        # Rainfall habitat saturation (log-logistic response)
        rain_suit = rain / (rain + TRUE_SYNTHETIC_PARAMETERS["K_rain_scale"])
        
        # Baseline vector abundance index
        vector_index = t_suit * rain_suit
        
        # Smooth with antecedent lag (1 month development delay)
        vector_lag = np.zeros_like(vector_index)
        for t in range(len(vector_index)):
            if t == 0:
                vector_lag[t] = vector_index[t]
            else:
                vector_lag[t] = 0.65 * vector_lag[t-1] + 0.35 * vector_index[t]
                
        # Simulate annual SMC intervention in Sahel / Northern Guinea Savanna starting around 2014
        # Standard calendar practice: July-October cycles
        zone_str = item["koppen"]
        has_smc = (zone_str in ["BSh", "Aw"])
        
        # Prevalence calculation
        base_prev = vector_lag * 0.55
        
        # Apply SMC impact for months 7,8,9,10 post-2014
        prev_series = base_prev.copy()
        dates_period = monthly_clim["year_month"].dt.to_timestamp()
        
        for t, dt in enumerate(dates_period):
            if has_smc and dt.year >= 2014 and dt.month in [7, 8, 9, 10]:
                prev_series[t] *= (1.0 - 0.55 * TRUE_SYNTHETIC_PARAMETERS["smc_cycle_efficacy"])
                
        # Add realistic observation / binomial sampling noise (N = 1200 children per survey)
        sample_size = 1200
        sim_positives = np.random.binomial(sample_size, np.clip(prev_series, 0.01, 0.85))
        noisy_prev = sim_positives / sample_size
        
        # Monthly clinical incidence per 1000 person-months
        clinical_inc = noisy_prev * 1000 * 0.12 * (vector_lag + 0.1)
        
        for t, row in monthly_clim.iterrows():
            all_records.append({
                "date": str(row["year_month"]),
                "state": state_name,
                "temp_mean": round(row["temp_mean"], 2),
                "rainfall_sum": round(row["rainfall"], 1),
                "humidity_mean": round(row["humidity"], 1),
                "vector_density_index": round(vector_lag[t], 4),
                "true_prevalence": round(prev_series[t], 4),
                "observed_prevalence": round(noisy_prev[t], 4),
                "clinical_incidence_per_1000": round(clinical_inc[t], 2),
                "sample_size": sample_size,
                "DATA_STATUS": "SYNTHETIC_TRUTH"
            })
            
    df_truth = pd.DataFrame(all_records)
    out_path = os.path.join(DATA_SYNTHETIC_DIR, "synthetic_truth_malaria.csv")
    df_truth.to_csv(out_path, index=False)
    print(f"[Synthetic Generator] Successfully generated {len(df_truth)} monthly records from 2000 to {end_date}.")
    print(f"[Synthetic Generator] Saved to: {out_path}")
    return df_truth

if __name__ == "__main__":
    generate_synthetic_truth_dataset()
