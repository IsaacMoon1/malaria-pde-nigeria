"""
Stage 2: Objective Agro-Climatic Zone Classification for 37 Nigerian Administrative Units.
Uses pure NumPy & SciPy (scipy.cluster.hierarchy & scipy.cluster.vq) to perform objective
K-Means and Hierarchical Clustering on Mean Annual Rainfall, Rainy-Season Length,
and Mean Annual Temperature. Identifies boundary/uncertain states and selects the optimal
centroid representative site per zone.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.cluster.vq import kmeans2

# Import project config
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
from config.config import (
    NIGERIA_ADMIN1, ZONE_NAMES, FIGURES_DIR, TABLES_DIR, OUTPUT_DIR, SEED
)

# Set deterministic random seed
np.random.seed(SEED)

CLIMATE_NORMALS_DATA = [
    {"state": "Abia", "annual_rainfall_mm": 2250, "rainy_season_months": 8.5, "mean_temperature_c": 26.4, "koppen": "Af"},
    {"state": "Adamawa", "annual_rainfall_mm": 980, "rainy_season_months": 5.5, "mean_temperature_c": 27.8, "koppen": "Aw"},
    {"state": "Akwa Ibom", "annual_rainfall_mm": 2850, "rainy_season_months": 9.5, "mean_temperature_c": 26.2, "koppen": "Af"},
    {"state": "Anambra", "annual_rainfall_mm": 1950, "rainy_season_months": 8.0, "mean_temperature_c": 26.8, "koppen": "Am"},
    {"state": "Bauchi", "annual_rainfall_mm": 850, "rainy_season_months": 4.5, "mean_temperature_c": 28.2, "koppen": "BSh"},
    {"state": "Bayelsa", "annual_rainfall_mm": 3200, "rainy_season_months": 10.0, "mean_temperature_c": 26.1, "koppen": "Af"},
    {"state": "Benue", "annual_rainfall_mm": 1350, "rainy_season_months": 6.5, "mean_temperature_c": 27.2, "koppen": "Aw"},
    {"state": "Borno", "annual_rainfall_mm": 550, "rainy_season_months": 3.5, "mean_temperature_c": 29.1, "koppen": "BSh"},
    {"state": "Cross River", "annual_rainfall_mm": 2900, "rainy_season_months": 9.0, "mean_temperature_c": 26.3, "koppen": "Af"},
    {"state": "Delta", "annual_rainfall_mm": 2650, "rainy_season_months": 9.0, "mean_temperature_c": 26.5, "koppen": "Af"},
    {"state": "Ebonyi", "annual_rainfall_mm": 1850, "rainy_season_months": 7.5, "mean_temperature_c": 26.9, "koppen": "Am"},
    {"state": "Edo", "annual_rainfall_mm": 2100, "rainy_season_months": 8.0, "mean_temperature_c": 26.6, "koppen": "Am"},
    {"state": "Ekiti", "annual_rainfall_mm": 1550, "rainy_season_months": 7.5, "mean_temperature_c": 25.8, "koppen": "Aw/Am"},
    {"state": "Enugu", "annual_rainfall_mm": 1780, "rainy_season_months": 7.5, "mean_temperature_c": 26.7, "koppen": "Am"},
    {"state": "FCT", "annual_rainfall_mm": 1400, "rainy_season_months": 6.5, "mean_temperature_c": 26.9, "koppen": "Aw"},
    {"state": "Gombe", "annual_rainfall_mm": 780, "rainy_season_months": 4.5, "mean_temperature_c": 28.0, "koppen": "BSh"},
    {"state": "Imo", "annual_rainfall_mm": 2300, "rainy_season_months": 8.5, "mean_temperature_c": 26.5, "koppen": "Af"},
    {"state": "Jigawa", "annual_rainfall_mm": 620, "rainy_season_months": 3.5, "mean_temperature_c": 28.7, "koppen": "BSh"},
    {"state": "Kaduna", "annual_rainfall_mm": 1200, "rainy_season_months": 5.5, "mean_temperature_c": 25.6, "koppen": "Aw"},
    {"state": "Kano", "annual_rainfall_mm": 690, "rainy_season_months": 4.0, "mean_temperature_c": 28.3, "koppen": "BSh"},
    {"state": "Katsina", "annual_rainfall_mm": 610, "rainy_season_months": 3.5, "mean_temperature_c": 28.5, "koppen": "BSh"},
    {"state": "Kebbi", "annual_rainfall_mm": 750, "rainy_season_months": 4.0, "mean_temperature_c": 28.9, "koppen": "BSh"},
    {"state": "Kogi", "annual_rainfall_mm": 1380, "rainy_season_months": 6.5, "mean_temperature_c": 27.5, "koppen": "Aw"},
    {"state": "Kwara", "annual_rainfall_mm": 1250, "rainy_season_months": 6.5, "mean_temperature_c": 26.8, "koppen": "Aw"},
    {"state": "Lagos", "annual_rainfall_mm": 1750, "rainy_season_months": 8.5, "mean_temperature_c": 27.1, "koppen": "Am"},
    {"state": "Nasarawa", "annual_rainfall_mm": 1320, "rainy_season_months": 6.5, "mean_temperature_c": 27.4, "koppen": "Aw"},
    {"state": "Niger", "annual_rainfall_mm": 1220, "rainy_season_months": 6.0, "mean_temperature_c": 27.3, "koppen": "Aw"},
    {"state": "Ogun", "annual_rainfall_mm": 1600, "rainy_season_months": 8.0, "mean_temperature_c": 26.9, "koppen": "Am"},
    {"state": "Ondo", "annual_rainfall_mm": 1900, "rainy_season_months": 8.0, "mean_temperature_c": 26.2, "koppen": "Am"},
    {"state": "Osun", "annual_rainfall_mm": 1450, "rainy_season_months": 7.5, "mean_temperature_c": 26.3, "koppen": "Aw/Am"},
    {"state": "Oyo", "annual_rainfall_mm": 1300, "rainy_season_months": 7.0, "mean_temperature_c": 26.7, "koppen": "Aw/Am"},
    {"state": "Plateau", "annual_rainfall_mm": 1450, "rainy_season_months": 6.5, "mean_temperature_c": 22.8, "koppen": "Aw (Highland)"},
    {"state": "Rivers", "annual_rainfall_mm": 2750, "rainy_season_months": 9.5, "mean_temperature_c": 26.3, "koppen": "Af"},
    {"state": "Sokoto", "annual_rainfall_mm": 640, "rainy_season_months": 3.5, "mean_temperature_c": 28.8, "koppen": "BSh"},
    {"state": "Taraba", "annual_rainfall_mm": 1250, "rainy_season_months": 6.0, "mean_temperature_c": 27.1, "koppen": "Aw"},
    {"state": "Yobe", "annual_rainfall_mm": 510, "rainy_season_months": 3.0, "mean_temperature_c": 29.2, "koppen": "BSh"},
    {"state": "Zamfara", "annual_rainfall_mm": 720, "rainy_season_months": 4.0, "mean_temperature_c": 28.4, "koppen": "BSh"}
]

def compute_silhouette(X, labels):
    """Exact Rousseeuw (1987) silhouette coefficient implementation."""
    dist_matrix = squareform(pdist(X, metric='euclidean'))
    n = len(labels)
    unique_labels = np.unique(labels)
    sils = np.zeros(n)
    
    for i in range(n):
        own_cluster = labels[i]
        own_mask = (labels == own_cluster)
        
        # Mean intra-cluster distance a(i)
        if np.sum(own_mask) > 1:
            a_i = np.sum(dist_matrix[i, own_mask]) / (np.sum(own_mask) - 1)
        else:
            a_i = 0.0
            
        # Smallest mean inter-cluster distance b(i)
        b_candidates = []
        for other_cluster in unique_labels:
            if other_cluster == own_cluster:
                continue
            other_mask = (labels == other_cluster)
            b_cand = np.mean(dist_matrix[i, other_mask])
            b_candidates.append(b_cand)
            
        b_i = min(b_candidates) if b_candidates else 0.0
        
        denom = max(a_i, b_i)
        sils[i] = (b_i - a_i) / denom if denom > 0 else 0.0
        
    return sils, np.mean(sils)

def perform_clustering():
    """Execute clustering and centroid analysis with pure SciPy."""
    df = pd.DataFrame(CLIMATE_NORMALS_DATA)
    df["lat"] = df["state"].apply(lambda s: NIGERIA_ADMIN1[s]["lat"])
    df["lon"] = df["state"].apply(lambda s: NIGERIA_ADMIN1[s]["lon"])
    
    features = ["annual_rainfall_mm", "rainy_season_months", "mean_temperature_c"]
    X = df[features].values.astype(float)
    
    # Standardize features (Z-score)
    mean = np.mean(X, axis=0)
    std = np.std(X, axis=0)
    X_scaled = (X - mean) / std
    
    # 1. K-Means (K=3) via scipy kmeans2
    centroids, kmeans_labels = kmeans2(X_scaled, 3, minit='points', iter=50, seed=SEED)
    
    # 2. Hierarchical Clustering (Ward linkage)
    Z = linkage(X_scaled, method='ward')
    hier_labels = fcluster(Z, t=3, criterion='maxclust') - 1
    
    # Map clusters to Zone 1, Zone 2, Zone 3 based on rainfall ranking
    k_rain_means = [df.loc[kmeans_labels == i, "annual_rainfall_mm"].mean() for i in range(3)]
    k_sorted = np.argsort(k_rain_means)
    k_map = {k_sorted[0]: "Zone 1", k_sorted[1]: "Zone 2", k_sorted[2]: "Zone 3"}
    df["assigned_zone"] = [k_map[lbl] for lbl in kmeans_labels]
    
    h_rain_means = [df.loc[hier_labels == i, "annual_rainfall_mm"].mean() for i in range(3)]
    h_sorted = np.argsort(h_rain_means)
    h_map = {h_sorted[0]: "Zone 1", h_sorted[1]: "Zone 2", h_sorted[2]: "Zone 3"}
    df["hierarchical_zone"] = [h_map[lbl] for lbl in hier_labels]
    
    df["methods_agree"] = df["assigned_zone"] == df["hierarchical_zone"]
    
    # Silhouette scores
    sils, overall_sil = compute_silhouette(X_scaled, kmeans_labels)
    df["silhouette_score"] = sils
    
    # Distance to Centroid
    df["dist_to_centroid"] = 0.0
    for i in range(3):
        zone_str = k_map[i]
        mask = (df["assigned_zone"] == zone_str).values
        sub_X = X_scaled[mask]
        c = centroids[i].reshape(1, -1)
        dists = cdist(sub_X, c, metric='euclidean').flatten()
        df.loc[mask, "dist_to_centroid"] = dists
        
    # Boundary / Uncertain designation
    df["boundary_uncertain"] = (df["silhouette_score"] < 0.35) | (~df["methods_agree"])
    
    # Select Representative Site per Zone
    representative_sites = {}
    for zone in ["Zone 1", "Zone 2", "Zone 3"]:
        sub_df = df[df["assigned_zone"] == zone]
        sorted_sub = sub_df.sort_values(by="dist_to_centroid")
        rep_state = sorted_sub.iloc[0]["state"]
        rep_dist = sorted_sub.iloc[0]["dist_to_centroid"]
        representative_sites[zone] = {
            "representative_site": rep_state,
            "dist_to_centroid": rep_dist,
            "num_states": len(sub_df)
        }
        
    return df, representative_sites, overall_sil

def export_clustering_outputs(df, rep_sites, overall_sil):
    csv_path = os.path.join(TABLES_DIR, "table2_climate_clustering.csv")
    df.to_csv(csv_path, index=False)
    
    tex_path = os.path.join(TABLES_DIR, "table2_climate_clustering.tex")
    with open(tex_path, "w") as f:
        f.write("% Objective Climate Clustering of 37 Nigerian Administrative Units\n")
        f.write("\\begin{tabular}{lcccccc}\n\\hline\n")
        f.write("State & Rainfall (mm) & Rainy Mos & Temp (°C) & Assigned Zone & Sil. Score & Status \\\\\n\\hline\n")
        for _, row in df.sort_values(by=["assigned_zone", "dist_to_centroid"]).iterrows():
            status = "Boundary/Uncertain" if row["boundary_uncertain"] else "Core"
            f.write(f"{row['state']} & {row['annual_rainfall_mm']} & {row['rainy_season_months']} & {row['mean_temperature_c']} & {row['assigned_zone']} & {row['silhouette_score']:.2f} & {status} \\\\\n")
        f.write("\\hline\n\\end{tabular}\n")
        
    # Plotting
    plt.rcParams.update({"font.size": 10, "font.family": "serif"})
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
    
    colors = {"Zone 1": "#d95f02", "Zone 2": "#7570b3", "Zone 3": "#1b9e77"}
    for zone in ["Zone 1", "Zone 2", "Zone 3"]:
        sub = df[df["assigned_zone"] == zone]
        core = sub[~sub["boundary_uncertain"]]
        axes[0].scatter(core["lon"], core["lat"], c=colors[zone], s=80, label=f"{zone}: {ZONE_NAMES[zone]}", edgecolors="black", linewidth=0.8)
        bnd = sub[sub["boundary_uncertain"]]
        if len(bnd) > 0:
            axes[0].scatter(bnd["lon"], bnd["lat"], c=colors[zone], s=120, marker="^", label=f"{zone} (Boundary)", edgecolors="red", linewidth=1.5)
        rep_state = rep_sites[zone]["representative_site"]
        rep_row = df[df["state"] == rep_state].iloc[0]
        axes[0].scatter(rep_row["lon"], rep_row["lat"], c="yellow", s=180, marker="*", edgecolors="black", linewidth=1.5, zorder=5)
        axes[0].annotate(f"{rep_state}\n(Centroid)", (rep_row["lon"]+0.15, rep_row["lat"]-0.15), fontweight="bold", fontsize=9, bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8))
        
    for b_state in ["Oyo", "Plateau", "Adamawa", "Kaduna"]:
        if b_state in df[df["boundary_uncertain"]]["state"].values:
            b_row = df[df["state"] == b_state].iloc[0]
            axes[0].annotate(f"{b_state}*", (b_row["lon"]-0.4, b_row["lat"]+0.15), fontsize=8, color="darkred", fontweight="bold")
            
    axes[0].set_title("(a) Objective Climate Zone Partition & Centroids (37 Units)", fontweight="bold")
    axes[0].set_xlabel("Longitude (°E)")
    axes[0].set_ylabel("Latitude (°N)")
    axes[0].legend(loc="lower left", fontsize=7.5)
    axes[0].grid(True, linestyle="--", alpha=0.5)
    
    for zone in ["Zone 1", "Zone 2", "Zone 3"]:
        sub = df[df["assigned_zone"] == zone]
        axes[1].scatter(sub["annual_rainfall_mm"], sub["rainy_season_months"], c=colors[zone], s=80, edgecolors="black", linewidth=0.8)
        rep_state = rep_sites[zone]["representative_site"]
        rep_row = df[df["state"] == rep_state].iloc[0]
        axes[1].scatter(rep_row["annual_rainfall_mm"], rep_row["rainy_season_months"], c="yellow", s=180, marker="*", edgecolors="black", linewidth=1.5, zorder=5)
        axes[1].annotate(f"Centroid: {rep_state}", (rep_row["annual_rainfall_mm"]+40, rep_row["rainy_season_months"]-0.2), fontweight="bold", fontsize=9)
        
    axes[1].set_title("(b) Climatological Feature Separation & Centroid Alignment", fontweight="bold")
    axes[1].set_xlabel("Mean Annual Rainfall (mm)")
    axes[1].set_ylabel("Rainy Season Duration (Months > 60 mm)")
    axes[1].grid(True, linestyle="--", alpha=0.5)
    
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "fig1_nigeria_climate_zones.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
    
    print(f"Clustering successfully completed! Overall Silhouette Score: {overall_sil:.3f}")
    print(f"Representative Sites Selected:")
    for z, info in rep_sites.items():
        print(f"  {z} ({ZONE_NAMES[z]}): {info['representative_site']} (Distance to centroid: {info['dist_to_centroid']:.3f}, States: {info['num_states']})")
    print(f"Artifacts saved:\n  Table: {csv_path}\n  Figure: {fig_path}")

if __name__ == "__main__":
    df, rep_sites, sil = perform_clustering()
    export_clustering_outputs(df, rep_sites, sil)
