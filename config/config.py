"""
Global configuration, random seed locks, and administrative metadata for Nigeria Malaria PDE Pipeline.
All random number generators must derive deterministically from SEED = 42.
Target: Elsevier Journal (Infectious Disease Modelling / Epidemics / Mathematical Biosciences)
"""

import os
import random
import numpy as np

# Deterministic Seed Enforcement
SEED = 42
random.seed(SEED)
np.random.seed(SEED)

# Workspace Directories
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(BASE_DIR, "data")
DATA_BOUNDARIES_DIR = os.path.join(DATA_DIR, "boundaries")
DATA_CLIMATE_DIR = os.path.join(DATA_DIR, "climate")
DATA_MALARIA_DIR = os.path.join(DATA_DIR, "malaria")
DATA_SYNTHETIC_DIR = os.path.join(DATA_DIR, "synthetic")

OUTPUT_DIR = os.path.join(BASE_DIR, "output")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")
TABLES_DIR = os.path.join(BASE_DIR, "tables")
PAPER_DIR = os.path.join(BASE_DIR, "paper")

for d in [DATA_DIR, DATA_BOUNDARIES_DIR, DATA_CLIMATE_DIR, DATA_MALARIA_DIR, 
          DATA_SYNTHETIC_DIR, OUTPUT_DIR, FIGURES_DIR, TABLES_DIR, PAPER_DIR]:
    os.makedirs(d, exist_ok=True)

# Data switch: "SYNTHETIC" or "REAL"
DATA_MODE = "SYNTHETIC"

# Benchmark Reference Date: Current local time reference
REFERENCE_CURRENT_DATE = "2026-11-02"
START_YEAR = 2000

# Geographic metadata: 36 Nigerian States + FCT (37 administrative units)
# Centroid coordinates (Latitude, Longitude) and initial regional grouping
NIGERIA_ADMIN1 = {
    "Abia": {"lat": 5.4527, "lon": 7.5248, "initial_zone": "Zone 3"},
    "Adamawa": {"lat": 9.3265, "lon": 12.3984, "initial_zone": "Zone 2"},
    "Akwa Ibom": {"lat": 5.0079, "lon": 7.8497, "initial_zone": "Zone 3"},
    "Anambra": {"lat": 6.2209, "lon": 7.0723, "initial_zone": "Zone 3"},
    "Bauchi": {"lat": 10.3158, "lon": 9.8442, "initial_zone": "Zone 1"},
    "Bayelsa": {"lat": 4.7719, "lon": 6.0699, "initial_zone": "Zone 3"},
    "Benue": {"lat": 7.3369, "lon": 8.7404, "initial_zone": "Zone 2"},
    "Borno": {"lat": 11.8333, "lon": 13.1500, "initial_zone": "Zone 1"},
    "Cross River": {"lat": 5.8702, "lon": 8.5988, "initial_zone": "Zone 3"},
    "Delta": {"lat": 5.7040, "lon": 5.9339, "initial_zone": "Zone 3"},
    "Ebonyi": {"lat": 6.2649, "lon": 8.0137, "initial_zone": "Zone 3"},
    "Edo": {"lat": 6.5438, "lon": 5.8987, "initial_zone": "Zone 3"},
    "Ekiti": {"lat": 7.7190, "lon": 5.3110, "initial_zone": "Zone 3"},
    "Enugu": {"lat": 6.5364, "lon": 7.4383, "initial_zone": "Zone 3"},
    "FCT": {"lat": 9.0765, "lon": 7.3986, "initial_zone": "Zone 2"},
    "Gombe": {"lat": 10.2897, "lon": 11.1673, "initial_zone": "Zone 1"},
    "Imo": {"lat": 5.5720, "lon": 7.0588, "initial_zone": "Zone 3"},
    "Jigawa": {"lat": 12.2280, "lon": 9.5616, "initial_zone": "Zone 1"},
    "Kaduna": {"lat": 10.5105, "lon": 7.4165, "initial_zone": "Zone 2"},
    "Kano": {"lat": 12.0022, "lon": 8.5919, "initial_zone": "Zone 1"},
    "Katsina": {"lat": 12.9908, "lon": 7.6018, "initial_zone": "Zone 1"},
    "Kebbi": {"lat": 12.4504, "lon": 4.1999, "initial_zone": "Zone 1"},
    "Kogi": {"lat": 7.7337, "lon": 6.6906, "initial_zone": "Zone 2"},
    "Kwara": {"lat": 8.9669, "lon": 4.6025, "initial_zone": "Zone 2"},
    "Lagos": {"lat": 6.5244, "lon": 3.3792, "initial_zone": "Zone 3"},
    "Nasarawa": {"lat": 8.5378, "lon": 8.5204, "initial_zone": "Zone 2"},
    "Niger": {"lat": 9.9309, "lon": 5.5983, "initial_zone": "Zone 2"},
    "Ogun": {"lat": 6.9980, "lon": 3.4737, "initial_zone": "Zone 3"},
    "Ondo": {"lat": 7.2571, "lon": 5.2058, "initial_zone": "Zone 3"},
    "Osun": {"lat": 7.5629, "lon": 4.5200, "initial_zone": "Zone 3"},
    "Oyo": {"lat": 7.8430, "lon": 3.9368, "initial_zone": "Zone 3"},  # Ecotone / border
    "Plateau": {"lat": 9.2182, "lon": 9.5179, "initial_zone": "Zone 2"},
    "Rivers": {"lat": 4.8156, "lon": 7.0498, "initial_zone": "Zone 3"},
    "Sokoto": {"lat": 13.0609, "lon": 5.2415, "initial_zone": "Zone 1"},
    "Taraba": {"lat": 7.8704, "lon": 9.7800, "initial_zone": "Zone 2"},
    "Yobe": {"lat": 12.2939, "lon": 11.4390, "initial_zone": "Zone 1"},
    "Zamfara": {"lat": 12.1222, "lon": 6.2236, "initial_zone": "Zone 1"},
}

ZONE_NAMES = {
    "Zone 1": "Sahel / Semi-Arid North",
    "Zone 2": "Guinea Savanna / Middle Belt",
    "Zone 3": "Humid Forest / Southern Belt"
}

# Proposed Representative Sites (to be confirmed dynamically via centroid clustering in Stage 2)
CANDIDATE_REPRESENTATIVE_SITES = {
    "Zone 1": "Sokoto",
    "Zone 2": "Kwara",     # Ilorin
    "Zone 3": "Oyo"        # Ibadan
}
