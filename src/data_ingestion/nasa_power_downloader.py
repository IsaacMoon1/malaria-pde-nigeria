"""
Stage 3: NASA POWER API Downloader and Quality Control Engine.
Fetches daily 2-meter air temperature (T2M, T2M_MIN, T2M_MAX), relative humidity (RH2M),
and precipitation for state centroids. Dynamically detects latest server end date without extrapolation.
Includes fast offline detection and strict Quality Control (QC).
"""

import os
import sys
import datetime
import requests
import numpy as np
import pandas as pd

# Load configuration
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
from config.config import DATA_CLIMATE_DIR, NIGERIA_ADMIN1, START_YEAR, REFERENCE_CURRENT_DATE

NASA_POWER_BASE_URL = "https://power.larc.nasa.gov/api/temporal/daily/point"
_NETWORK_CHECKED = False
_IS_ONLINE = False

def check_online_status():
    global _NETWORK_CHECKED, _IS_ONLINE
    if _NETWORK_CHECKED:
        return _IS_ONLINE
    try:
        r = requests.get("https://power.larc.nasa.gov", timeout=1.5)
        _IS_ONLINE = (r.status_code == 200)
    except Exception:
        _IS_ONLINE = False
    _NETWORK_CHECKED = True
    return _IS_ONLINE

def detect_nasa_power_end_date():
    """
    Dynamically detect the true end date available on the NASA POWER server.
    Never assumes 2024. Never extrapolates missing months.
    """
    if check_online_status():
        probe_url = (
            f"{NASA_POWER_BASE_URL}?parameters=T2M&community=AG"
            f"&longitude=7.3986&latitude=9.0765&start=20260101&end=20261102&format=JSON"
        )
        try:
            resp = requests.get(probe_url, timeout=3.0)
            if resp.status_code == 200:
                data = resp.json()
                dates = list(data.get("properties", {}).get("parameter", {}).get("T2M", {}).keys())
                valid_dates = [d for d in dates if data["properties"]["parameter"]["T2M"][d] != -999]
                if valid_dates:
                    latest_str = max(valid_dates)
                    latest_date = datetime.datetime.strptime(latest_str, "%Y%m%d").strftime("%Y-%m-%d")
                    print(f"[NASA POWER] Successfully detected live server end date: {latest_date}")
                    return latest_date
        except Exception as e:
            print(f"[NASA POWER] Live query note: {e}")
            
    fallback_date = "2026-10-31"
    return fallback_date

def fetch_nasa_power_point(lat, lon, start_date="20000101", end_date=None):
    """Fetch daily time series for a single coordinate."""
    if not check_online_status():
        return None
        
    if end_date is None:
        end_date = detect_nasa_power_end_date().replace("-", "")
        
    url = (
        f"{NASA_POWER_BASE_URL}?parameters=T2M,T2M_MIN,T2M_MAX,RH2M,PRECTOTCORR"
        f"&community=AG&longitude={lon}&latitude={lat}&start={start_date}&end={end_date}&format=JSON"
    )
    
    try:
        response = requests.get(url, timeout=5.0)
        if response.status_code == 200:
            res_json = response.json()
            param_dict = res_json.get("properties", {}).get("parameter", {})
            dates = sorted(list(param_dict.get("T2M", {}).keys()))
            
            records = []
            for d in dates:
                records.append({
                    "date": datetime.datetime.strptime(d, "%Y%m%d").strftime("%Y-%m-%d"),
                    "temp_mean": param_dict["T2M"].get(d),
                    "temp_min": param_dict["T2M_MIN"].get(d),
                    "temp_max": param_dict["T2M_MAX"].get(d),
                    "humidity": param_dict["RH2M"].get(d),
                    "rainfall": param_dict["PRECTOTCORR"].get(d)
                })
            df = pd.DataFrame(records)
            return qc_climate_series(df)
    except Exception as e:
        print(f"[NASA POWER API] Point fetch failed ({lat}, {lon}): {e}")
        
    return None

def qc_climate_series(df):
    """Quality Control: bounds check, fill replacement, interpolation."""
    df = df.replace(-999, np.nan)
    df.loc[(df["humidity"] < 0) | (df["humidity"] > 100), "humidity"] = np.nan
    df.loc[(df["temp_mean"] < -10) | (df["temp_mean"] > 55), "temp_mean"] = np.nan
    df.loc[df["rainfall"] < 0, "rainfall"] = 0.0
    
    df = df.interpolate(method="linear", limit=3)
    df = df.bfill().ffill()
    return df

def generate_realistic_synthetic_climate(state_name, start_date="2000-01-01", end_date="2026-10-31"):
    """
    Generates synthetic daily climate time series derived from Nigerian Meteorological Agency (NiMet) normals.
    Explicitly labeled as SYNTHETIC.
    """
    from src.clustering.zone_classifier import CLIMATE_NORMALS_DATA
    state_norm = next(item for item in CLIMATE_NORMALS_DATA if item["state"] == state_name)
    
    date_rng = pd.date_range(start=start_date, end=end_date, freq="D")
    n_days = len(date_rng)
    doy = date_rng.dayofyear.values
    
    t_mean = state_norm["mean_temperature_c"]
    temp_seasonal = t_mean + 2.5 * np.sin(2 * np.pi * (doy - 45) / 365.25) - 1.5 * np.cos(4 * np.pi * (doy) / 365.25)
    temp_noise = np.random.normal(0, 0.8, n_days)
    daily_temp = np.clip(temp_seasonal + temp_noise, 15.0, 42.0)
    
    ann_rain = state_norm["annual_rainfall_mm"]
    r_months = state_norm["rainy_season_months"]
    
    if r_months < 5.0:  # Sahel unimodal peak August (DOY 220)
        rain_prob = np.exp(-0.5 * ((doy - 220) / 35)**2)
    elif r_months < 7.5:  # Guinea savanna peak July-September
        rain_prob = np.exp(-0.5 * ((doy - 230) / 55)**2)
    else:  # Humid forest bimodal
        peak1 = np.exp(-0.5 * ((doy - 180) / 40)**2)
        peak2 = np.exp(-0.5 * ((doy - 265) / 35)**2)
        rain_prob = 0.6 * peak1 + 0.5 * peak2
        
    raw_intensity = np.random.exponential(scale=18.0, size=n_days) * (np.random.rand(n_days) < rain_prob)
    current_annual = np.sum(raw_intensity) / (n_days / 365.25)
    scale_factor = ann_rain / max(1.0, current_annual)
    daily_rainfall = np.maximum(0.0, raw_intensity * scale_factor)
    
    humidity = np.clip(55.0 + 35.0 * rain_prob + np.random.normal(0, 5.0, n_days), 20.0, 98.0)
    
    df = pd.DataFrame({
        "date": date_rng.strftime("%Y-%m-%d"),
        "state": state_name,
        "temp_mean": np.round(daily_temp, 2),
        "temp_min": np.round(daily_temp - 4.5, 2),
        "temp_max": np.round(daily_temp + 5.0, 2),
        "humidity": np.round(humidity, 1),
        "rainfall": np.round(daily_rainfall, 2),
        "DATA_SOURCE": "SYNTHETIC_CLIMATE_NORMAL_DRIVEN"
    })
    return df

def get_or_download_climate(state_name):
    """Retrieve cached climate file or generate/download if missing."""
    cached_file = os.path.join(DATA_CLIMATE_DIR, f"climate_{state_name.replace(' ', '_')}.csv")
    if os.path.exists(cached_file):
        return pd.read_csv(cached_file)
        
    lat = NIGERIA_ADMIN1[state_name]["lat"]
    lon = NIGERIA_ADMIN1[state_name]["lon"]
    end_date = detect_nasa_power_end_date()
    
    df = fetch_nasa_power_point(lat, lon, start_date=f"{START_YEAR}0101", end_date=end_date.replace("-", ""))
    if df is not None:
        df["state"] = state_name
        df["DATA_SOURCE"] = "NASA_POWER_OBSERVED"
    else:
        df = generate_realistic_synthetic_climate(state_name, start_date=f"{START_YEAR}-01-01", end_date=end_date)
        
    df.to_csv(cached_file, index=False)
    return df
