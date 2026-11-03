# Statistical Validation and Uncertainty Quantification of a Reaction-Diffusion PDE Model for Malaria Intervention Timing: A Multi-Zone Synthetic Benchmark

**Isaac Tope Adeyemi**<sup>1,*</sup>

<sup>1</sup> *Department of Pure and Applied Mathematics, Faculty of Pure and Applied Sciences, Ladoke Akintola University of Technology (LAUTECH), P.M.B. 4000, Ogbomoso, Oyo State, Nigeria*  
\* *Corresponding author. Email: isaacadeyemi245@gmail.com*

---

## Abstract

Spatial reaction-diffusion partial differential equation (PDE) models capture the localized transmission gradients occurring between vector breeding habitats and residential settlements. However, their operational translation requires formal numerical verification and statistical uncertainty quantification (UQ). Here, we formulate, verify, and validate a hybrid reaction-diffusion PDE coupled human-mosquito transmission model ($S_H, E_H, I_H, R_H, P_H$) across a 1D spatial transect. To reflect diverse seasonal regimes, we establish a synthetic benchmark spanning 37 administrative units classified into three eco-climatic clusters: semi-arid (Zone 1, $n=10$), savanna (Zone 2, $n=12$), and humid forest (Zone 3, $n=15$), with an overall mean silhouette width of $0.4590$. The numerical method of lines (MOL) achieves an observed spatial convergence order of $p = 2.011$ against an exact diffusion solution and conserves human mass to machine precision ($\le 3.020 \times 10^{-15}$). Global Sobol sensitivity analysis demonstrates that transmission efficiencies ($b_H$, total-order $S_T = 0.3094–0.4273$; $b_M$, $S_T = 0.1718–0.3392$) dominate prevalence variance. In an out-of-sample temporal holdout (94 months), Markov Chain Monte Carlo (MCMC) parameter ensembles calibrated on centroid sites generalize across uncalibrated transfer states within the same climate cluster (Zone 1 transfer RMSE $0.0829 \pm 0.0031$ vs centroid $0.0833$), though empirical coverage exhibits under-coverage in humid regimes ($65.1\%$). In the semi-arid benchmark, decision-theoretic regret analysis illustrates that intervention timing synchronized with peak vector abundance averts $4,711.5$ case units compared to $4,082.7$ units for premature initiation. This benchmark provides a verified computational foundation for spatial transmission modeling.

**Keywords:** Mathematical epidemiology; Reaction-diffusion equations; Seasonal Malaria Chemoprevention; Uncertainty quantification; Bayesian calibration; Decision analysis

---

## 1. Introduction

Malaria caused by *Plasmodium falciparum* remains a critical global public health challenge, with sub-Saharan Africa bearing the overwhelming majority of illness and mortality (World Health Organization [WHO], 2023). Seasonal Malaria Chemoprevention (SMC)—the intermittent administration of full therapeutic courses of sulfadoxine-pyrimethamine plus amodiaquine (SPAQ) to children during high-transmission windows—has proven highly effective at preventing clinical episodes when chemoprotective blood titers coincide with seasonal vector transmission (Cairns et al., 2012; Greenwood, 2010; Meremikwu et al., 2012). 

However, the efficacy of chemoprotective interventions depends substantially on operational deployment timing. In seasonal environments, monsoonal precipitation drives larval habitat expansion, followed by temperature-dependent aquatic development, vector emergence, and extrinsic incubation of sporozoites (Bhatt et al., 2015; Sultan & Janicot, 2003; Weiss et al., 2019). When intervention rounds are initiated too early or too late relative to peak infectious vector inoculation, prophylactic drug concentrations wane prematurely or miss the main wave of transmission, generating substantial opportunity loss (Cairns et al., 2012).

Classical ordinary differential equation (ODE) compartmental models of the Ross-Macdonald lineage assume spatially homogeneous populations (Griffin et al., 2010; Smith et al., 2007). Such well-mixed frameworks cannot resolve the localized spatial boundary dynamics that occur between aquatic vector breeding habitats (e.g., river floodplains, borrow pits, wetlands) and residential human settlements (Depinay et al., 2004; Service, 1993). Female *Anopheles* mosquitoes disperse actively via spatial diffusion in search of blood meals, creating spatial gradients of infectious biting density (Lou & Zhao, 2010; Okuneye & Gumel, 2017). Reaction-diffusion partial differential equations (PDEs) provide an expressive mathematical formulation for coupling vector spatial dispersal with host transmission dynamics (Lou & Zhao, 2010).

Despite their mathematical realism, PDE epidemiological models are rarely subjected to end-to-end numerical verification and statistical uncertainty quantification (UQ). Computational disease models frequently suffer from three limitations:
1. **Unverified Numerics:** Discretization schemes are rarely benchmarked for empirical spatial and temporal convergence orders, manufactured solution accuracy, or closed-system mass conservation.
2. **Weak Identifiability & Diagnostics:** Multi-parameter spatial systems are often simulated without assessing parameter identifiability, Markov chain Monte Carlo (MCMC) mixing, or out-of-sample predictive uncertainty (Gelman et al., 2013; Saltelli et al., 2008).
3. **Unvalidated Spatial Transferability:** Calibrated model parameters are seldom tested for their capacity to generalize across diverse climatological settings without refitting (Gneiting & Raftery, 2007).

To address these challenges within an in silico setting, this paper develops a verified, coupled reaction-diffusion PDE-ODE modeling framework and evaluates its numerical stability, sensitivity, calibration identifiability, and spatial transferability across a multi-zone synthetic benchmark reflecting West African bioclimatic gradients. Specifically, we:
1. Establish a multi-zone synthetic benchmark across 37 units partitioned into three distinct seasonal regimes via unsupervised clustering, with formal silhouette diagnostics.
2. Formulate a 1D reaction-diffusion PDE system coupling four mosquito compartments ($L, M_S, M_E, M_I$) with a five-compartment human epidemiological framework ($S_H, E_H, I_H, R_H, P_H$).
3. Subject the numerical solver to six formal verification tests, confirming spatial convergence order ($p = 2.011$), manufactured solution fidelity, cross-solver concordance, and mass conservation.
4. Execute variance-based global sensitivity analysis (Sobol indices) to identify dominant biological drivers.
5. Perform Bayesian MCMC calibration on zone centroid sites, diagnosing parameter recovery, autocorrelation time, and identifiability limits.
6. Evaluate out-of-sample temporal holdout and spatial transferability across 34 uncalibrated transfer units.
7. Demonstrate a decision-theoretic regret analysis framework to quantify opportunity loss from intervention mistiming.

---

## 2. Climatological Zonation & Benchmark Architecture

### 2.1 Multi-Zone Synthetic Benchmark Design
To provide an in silico testing ground representative of sub-Saharan malaria transmission regimes, we constructed a 37-unit synthetic benchmark parameterized by long-term bioclimatic normals derived from West African meteorological patterns (NiMet, 2023; Sultan & Janicot, 2003). Each administrative unit $i \in \{1, \dots, 37\}$ is defined by three standardized features:
1. **Mean Annual Precipitation ($P_{\text{ann}}$, mm/year):** Cumulative annual precipitation.
2. **Rainy Season Duration ($L_{\text{rain}}$, months):** Months with mean precipitation $\ge 60\text{ mm}$.
3. **Mean Ambient Temperature ($\bar{T}$, $^\circ\text{C}$):** Mean temperature influencing vector development and extrinsic incubation rates.

### 2.2 Unsupervised Clustering and Silhouette Diagnostics
Features were standardized to unit variance. We evaluated cluster configurations for $k \in \{2, 3, 4, 5\}$ using $k$-means and Ward agglomerative hierarchical clustering (Ward, 1963). Exact silhouette coefficients (Rousseeuw, 1987) were computed for each configuration:
$$s(i) = \frac{b(i) - a(i)}{\max(a(i), b(i))}$$
where $a(i)$ is the mean intra-cluster distance and $b(i)$ is the minimum mean distance to any alternative cluster. 

The evaluation yielded mean silhouette widths of $0.5573$ for $k=2$, $0.4590$ for $k=3$, $0.3880$ for $k=4$, and $0.3366$ for $k=5$. While $k=2$ achieved higher mathematical compactness, $k=3$ was retained to preserve the three distinct eco-epidemiological regimes well-documented in tropical disease ecology (semi-arid, savanna, and humid forest).

```
       Cluster 1: Semi-Arid / Sahelian Regime (n=10 units)
       Rainfall: 510–850 mm | Rainy Season: 3.0–4.5 months | Temp: 28.0–29.2 °C
       Mean Silhouette: 0.7376 | Centroid Site: Jigawa (dist = 0.183282)
                                  |
                                  v
       Cluster 2: Savanna / Middle Belt Regime (n=12 units)
       Rainfall: 980–1,600 mm | Rainy Season: 5.5–8.0 months | Temp: 25.6–27.8 °C
       Mean Silhouette: 0.5179 | Centroid Site: FCT (dist = 0.135474)
                                  |
                                  v
       Cluster 3: Humid Forest / Southern Regime (n=15 units)
       Rainfall: 1,450–3,200 mm | Rainy Season: 6.5–10.0 months | Temp: 22.8–27.1 °C
       Mean Silhouette: 0.2261 | Centroid Site: Abia (dist = 0.175118)
```

The 37 units partitioned into:
- **Zone 1 (Semi-Arid, $n=10$):** High silhouette cohesion (mean $s = 0.7376$, zero negative silhouettes). Representative centroid site: **Jigawa** ($s = 0.819235$, distance to centroid $= 0.183282$).
- **Zone 2 (Savanna, $n=12$):** Moderate cohesion (mean $s = 0.5179$). Representative centroid site: **FCT** ($s = 0.672892$, distance to centroid $= 0.135474$).
- **Zone 3 (Humid Forest, $n=15$):** Substantial cluster overlap (mean $s = 0.2261$). Four states exhibited negative silhouette widths (Ebonyi $-0.164$, Ekiti $-0.192$, Enugu $-0.197$, Lagos $-0.049$), reflecting transitional boundaries between coastal maritime and interior savanna climates. Across the 37 units, 11 disagreements occurred between $k$-means and Ward linkage, illustrating significant boundary uncertainty that was tracked in downstream validation.

Table 2 reports the complete clustering results and boundary classifications.

---

## 3. Mathematical Model Formulation

### 3.1 Spatial Settlement-Breeding Domain
The model is formulated across a 1D spatial transect $\Omega = [0, L_x]$ with $L_x = 5.0\text{ km}$, representing the boundary interface between a localized vector breeding site and an adjacent residential settlement. The sub-domain $x \in [0, 1.5]\text{ km}$ represents the primary larval habitat zone ($\chi_{\text{breed}}(x) = 1.0$), while $x \in [2.0, 5.0]\text{ km}$ represents the human residential community.

Adult mosquito dispersal is governed by spatial diffusion with isotropic diffusivity $D_M$ ($\text{km}^2/\text{day}$). Homogeneous Neumann (zero-flux) boundary conditions are enforced on the boundaries $\partial \Omega$:
$$\left. \frac{\partial M_S}{\partial x} \right|_{x=0, L_x} = \left. \frac{\partial M_E}{\partial x} \right|_{x=0, L_x} = \left. \frac{\partial M_I}{\partial x} \right|_{x=0, L_x} = 0$$

### 3.2 Mosquito Dynamics (PDE System)
The mosquito population is structured into aquatic immatures (eggs, larvae, pupae) $L(x,t)$, susceptible adult females $M_S(x,t)$, latently infected (exposed) females $M_E(x,t)$, and infectious adult females $M_I(x,t)$:

$$\frac{\partial L}{\partial t} = \beta_E \cdot \left( M_S + M_E + M_I \right) \left( 1 - \frac{L}{K(R, x)} \right) - \left( \gamma_L(T) + d_L(T) \right) L$$

$$\frac{\partial M_S}{\partial t} = D_M \frac{\partial^2 M_S}{\partial x^2} + \frac{1}{2} \gamma_L(T) L - \lambda_M(x,t) M_S - \mu_M(T) M_S$$

$$\frac{\partial M_E}{\partial t} = D_M \frac{\partial^2 M_E}{\partial x^2} + \lambda_M(x,t) M_S - \left( \gamma_M(T) + \mu_M(T) \right) M_E$$

$$\frac{\partial M_I}{\partial t} = D_M \frac{\partial^2 M_I}{\partial x^2} + \gamma_M(T) M_E - \mu_M(T) M_I$$

where $\lambda_M(x,t) = a(T) b_M \frac{I_H(x,t)}{N_H(x,t)}$ is the vector force of infection, $a(T)$ is the temperature-dependent human-biting frequency, $b_M$ is the transmission probability per bite to mosquito, and the factor $1/2$ reflects a 1:1 sex ratio at adult emergence.

Thermal response functions for larval development $\gamma_L(T)$, extrinsic incubation rate $\gamma_M(T)$, and adult vector mortality $\mu_M(T)$ are parameterized following non-linear Briere and Arrhenius formulations (Mordecai et al., 2013, 2019):
$$\gamma_M(T) = \frac{T - T_{\text{min,sp}}}{\text{DD}_{\text{sp}}} \quad \text{for } T > T_{\text{min,sp}}$$

Hydrological carrying capacity $K(R, x)$ couples local daily rainfall $R(t)$ ($\text{mm}/\text{day}$) to breeding surface area (Parham & Michael, 2010):
$$K(R, x) = K_0 \cdot \chi_{\text{breed}}(x) \cdot \frac{R(t)}{R(t) + K_{\text{half}}}$$

### 3.3 Human Dynamics and Intervention Regimes
The human population is formulated as a five-compartment SEIRS system with chemoprotective prophylaxis ($S_H, E_H, I_H, R_H, P_H$), where total host density $N_H(x)$ is conserved:

$$\frac{\partial S_H}{\partial t} = \mu_H N_H + \omega_H R_H + \psi_{\text{SMC}} P_H - \lambda_H(x,t) S_H - \mu_H S_H - \rho_{\text{SMC}}(t) S_H$$

$$\frac{\partial E_H}{\partial t} = \lambda_H(x,t) S_H - \left( \sigma_H + \mu_H \right) E_H - \rho_{\text{SMC}}(t) E_H$$

$$\frac{\partial I_H}{\partial t} = \sigma_H E_H - \left( \gamma_H + \mu_H \right) I_H - \rho_{\text{SMC}}(t) \eta_{\text{clear}} I_H$$

$$\frac{\partial R_H}{\partial t} = \gamma_H I_H - \left( \omega_H + \mu_H \right) R_H - \rho_{\text{SMC}}(t) R_H$$

$$\frac{\partial P_H}{\partial t} = \rho_{\text{SMC}}(t) \left( S_H + E_H + \eta_{\text{clear}} I_H + R_H \right) - \left( \psi_{\text{SMC}} + \mu_H \right) P_H$$

Summing the human equations yields:
$$\frac{\partial}{\partial t}\left( S_H + E_H + I_H + R_H + P_H \right) = \mu_H N_H - \mu_H \left( S_H + E_H + I_H + R_H + P_H \right) = 0$$
guaranteeing total population mass conservation. The human force of infection is:
$$\lambda_H(x,t) = a(T) b_H \frac{M_I(x,t)}{N_H(x,t)}$$

Interventions are implemented as operational pulses:
- **Seasonal Malaria Chemoprevention (SMC):** Four discrete monthly cycles deployed at 28-day intervals ($t_k = t_{\text{start}} + 28k\text{ days}$). Each round clears existing parasitemia with clearance efficacy $\eta_{\text{clear}} = 0.95$ and transitions recipients into the protected compartment $P_H$. Prophylactic protection decays at rate $\psi_{\text{SMC}} = 1/28\text{ day}^{-1}$ (Cairns et al., 2012).
- **Indoor Residual Spraying (IRS):** Implemented as an exponential pulse in adult vector mortality: $d_{\text{IRS}}(t) = d_0 \exp(-\kappa_{\text{IRS}} (t - t_{\text{IRS}}))$.
- **Larviciding:** Modulates aquatic carrying capacity: $K_{\text{eff}}(t) = K(R, x) \cdot (1 - \epsilon_{\text{larv}}(t))$.

Table 1 summarizes all 20 literature-derived parameters and their associated DOIs.

---

## 4. Numerical Verification & Algorithmic Concordance

### 4.1 Discretization and Convergence Testing
The 1D spatial domain is discretized into $N_x$ nodes with uniform spacing $\Delta x = L_x / (N_x - 1)$. Spatial derivatives are approximated using 2nd-order central finite differences:
$$\left. \frac{\partial^2 u}{\partial x^2} \right|_{x_i} = \frac{u_{i+1} - 2u_i + u_{i-1}}{\Delta x^2} + \mathcal{O}(\Delta x^2)$$

Temporal integration is executed via a 2nd-order explicit Runge-Kutta midpoint integrator with adaptive sub-stepping constrained by the parabolic stability condition $\Delta t \le \Delta x^2 / (2 D_M)$.

### 4.2 Verification Benchmarks
The numerical solver was evaluated across six formal verification tests (Table 3):
1. **Analytic Diffusion Accuracy:** Tested against the exact fundamental solution $u(x,t) = \exp(-\pi^2 D_M t / L_x^2) \cos(\pi x / L_x)$. At $N_x = 129$, the observed spatial $L_2$ error was $5.055 \times 10^{-6}$ (pre-specified tolerance: $< 1.0 \times 10^{-4}$).
2. **Observed Spatial Convergence Order:** Refining the grid resolution across $N_x \in \{17, 33, 65, 129\}$ yielded an empirical convergence order of:
   $$p = \frac{\ln(E_{65} / E_{129})}{\ln(2)} = 2.011$$
   consistent with theoretical second-order accuracy ($p \approx 2.0 \pm 0.1$).
3. **Temporal Integration Stability:** Verified under adaptive Runge-Kutta-Fehlberg (RK45, $\text{rtol}=10^{-4}, \text{atol}=10^{-6}$) and fixed-step midpoint integration.
4. **Independent Solver Cross-Check:** Numerical trajectories were compared against an independently implemented finite-difference code (`IndependentPDESolver`). Cross-solver correlation reached Pearson $r = 1.00000$ with Normalized Root Mean Square Error $\text{NRMSE} = 0.01\%$.
5. **Human Mass Conservation:** Total population density was tracked over 365 days. The maximum relative mass leak was $3.020 \times 10^{-15}$, confirming conservation to machine precision.
6. **Biological Invariance:** All state variables strictly satisfied non-negativity: $\min_{i,t} y_i(x,t) \ge 0.0$.

---

## 5. Global Sensitivity Analysis

A variance-based global sensitivity analysis was conducted using the Sobol (2001) decomposition method via the `SALib` library (Herman & Usher, 2017). Eight parameters were sampled over biological intervals: $D_M \in [0.01, 0.25]$, biting rate $a \in [0.2, 0.5]$, human transmission efficiency $b_H \in [0.25, 0.8]$, mosquito susceptibility $b_M \in [0.2, 0.7]$, human clearance rate $\gamma_H \in [0.005, 0.03]$, immunity waning $\omega_H \in [0.001, 0.01]$, rainfall carrying capacity scale $K_{\text{scale}} \in [0.5, 2.0]$, and SMC coverage $\rho_{\text{SMC}} \in [0.5, 0.95]$.

Saltelli sampling was deployed with base sample size $N = 256$, generating $N(2k + 2) = 4,608$ model evaluations per ecological zone. First-order ($S_1$) and total-order ($S_T$) indices were computed for community infection prevalence (Table 4).

Transmission efficiencies ($b_H, b_M$) accounted for the majority of variance across all three ecological regimes:
- In Zone 1 (Semi-Arid), $S_{T, b_H} = 0.4273$ and $S_{T, b_M} = 0.2478$.
- In Zone 2 (Savanna), $S_{T, b_H} = 0.3557$ and $S_{T, b_M} = 0.3392$.
- In Zone 3 (Humid Forest), $S_{T, b_H} = 0.3094$ and $S_{T, b_M} = 0.1718$.

Adult mosquito spatial diffusivity ($D_M$) exhibited total-order sensitivity in Zone 2 ($S_T = 0.1131$), reflecting vector dispersal across transitional savanna settlement perimeters, whereas immunity waning ($\omega_H$) contributed $S_T = 0.1478$ in Zone 1.

---

## 6. Bayesian Calibration & Parameter Recovery

### 6.1 Calibration Framework
Calibration was performed using the affine-invariant ensemble MCMC sampler `emcee` (Foreman-Mackey et al., 2013) on representative centroid sites: Jigawa (Zone 1), FCT (Zone 2), and Abia (Zone 3).

The parameter vector $\boldsymbol{\theta} = \{D_M, b_H, b_M, \gamma_H\}$ was assigned bounded uniform prior distributions matching physiological constraints:
$$\pi(\boldsymbol{\theta}) \propto \mathcal{U}(0.01, 0.25) \times \mathcal{U}(0.25, 0.8) \times \mathcal{U}(0.2, 0.7) \times \mathcal{U}(0.005, 0.03)$$

Under the synthetic ground truth benchmark (generated directly with known parameters $\boldsymbol{\theta}^* = \{0.085, 0.500, 0.450, 0.0167\}$), the Gaussian pseudo-likelihood is:
$$\ln \mathcal{L}(\boldsymbol{y}_{\text{obs}} | \boldsymbol{\theta}) = -\frac{1}{2} \sum_{t=1}^{T} \left[ \frac{(y_{\text{obs}}(t) - y_{\text{sim}}(t; \boldsymbol{\theta}))^2}{\sigma_{\text{obs}}^2} + \ln(2\pi \sigma_{\text{obs}}^2) \right]$$
with observational standard deviation $\sigma_{\text{obs}} = 0.05$.

### 6.2 Convergence Diagnostics and Parameter Recovery
Ensemble sampling was conducted with 16 parallel walkers over 600 iterations per walker, discarding the initial 200 iterations as burn-in (Table 5):
- In Zone 1 (Jigawa), posterior means and 95% highest posterior density intervals (HPDI) bounded the true parameters: $D_M = 0.0965$ [0.0260, 0.1886] (true $0.0850$), $b_H = 0.4519$ [0.2683, 0.7228] (true $0.5000$), and $b_M = 0.2916$ [0.2080, 0.5251] (true $0.4500$).
- In Zone 2 (FCT), recovered values were $D_M = 0.0999$ [0.0327, 0.1952], $b_H = 0.5180$ [0.3392, 0.7452], and $b_M = 0.3677$ [0.2012, 0.6924].
- In Zone 3 (Abia), recovered values were $D_M = 0.1236$ [0.0249, 0.1976], $b_H = 0.5464$ [0.3282, 0.7124], and $b_M = 0.5142$ [0.2813, 0.6831].
- Human clearance rate ($\gamma_H$) exhibited structural shrinkage toward the upper boundary ($0.0270–0.0284$), illustrating parameter correlation between clinical duration and waning immunity.

---

## 7. Spatiotemporal Validation & Spatial Transferability

### 7.1 Temporal Holdout Scheme
To evaluate predictive generalizability, time series were partitioned into a historical calibration window (January 2000 to December 2018; 228 monthly observations) and an out-of-sample temporal holdout (January 2019 to October 2026; 94 holdout months). Predictions were scored via Root Mean Square Error (RMSE), Continuous Ranked Probability Score (CRPS) (Gneiting & Raftery, 2007), and 95% credible interval empirical coverage.

### 7.2 Spatial Transferability Across All 37 Units
To test spatial transferability, posterior parameter distributions calibrated exclusively on zone centroid sites were transferred—without recalibration—to simulate the remaining states in each respective zone under local climate forcing (Table 6):
- **Zone 1 (Semi-Arid):** The mean generalization performance across the 9 uncalibrated transfer states was $\text{RMSE} = 0.0829 \pm 0.0031$, $\text{CRPS} = 0.0502$, and 95% coverage $= 79.7\%$, closely matching the calibrated Jigawa centroid ($\text{RMSE} = 0.0833$, $\text{CRPS} = 0.0514$, coverage $= 79.8\%$).
- **Zone 2 (Savanna):** The 11 transfer states achieved a mean $\text{RMSE} = 0.0726 \pm 0.0029$, $\text{CRPS} = 0.0395$, and coverage $= 84.3\%$ (centroid FCT: $\text{RMSE} = 0.0710$, $\text{CRPS} = 0.0378$, coverage $= 86.2\%$).
- **Zone 3 (Humid Forest):** Transfer states yielded a mean $\text{RMSE} = 0.0870 \pm 0.0092$, $\text{CRPS} = 0.0526$, and coverage $= 65.1\%$ (centroid Abia: $\text{RMSE} = 0.0865$, $\text{CRPS} = 0.0547$, coverage $= 59.6\%$).

The consistency between centroid and spatial transfer performance demonstrates that agro-climatic clustering captures the primary drivers of seasonal transmission variation.

---

## 8. Decision Analysis & Intervention Timing Policy

### 8.1 Bayesian Decision Formulation
We formulate SMC deployment timing as a Bayesian decision problem under climate and parameter uncertainty (Berger, 1985). Five candidate initiation dates were evaluated:
$$\mathcal{A} = \{\text{May 15 (DOY 135)}, \text{June 15 (DOY 166)}, \text{July 15 (DOY 196)}, \text{August 15 (DOY 227)}, \text{September 15 (DOY 258)}\}$$

For each candidate start date $a \in \mathcal{A}$ and posterior parameter sample $\boldsymbol{\theta}^{(m)}$, the model evaluates cumulative clinical cases averted over the transmission season:
$$U(a, \boldsymbol{\theta}^{(m)}) = \int_{0}^{T_{\text{seas}}} \left( I_{H, \text{baseline}}(x, t) - I_{H, \text{SMC}}(x, t; a) \right) dt$$

The expected opportunity loss (regret) relative to the optimal schedule is:
$$\mathcal{R}(a) = \max_{a' \in \mathcal{A}} \bar{U}(a') - \bar{U}(a)$$

### 8.2 Timing Regret Across Climate Regimes
As documented in Table 7:
- In **Zone 1 (Semi-Arid)**, the optimal start date is **August 15**, achieving $4,711.5$ cases averted. Initiating SMC on July 15 (premature initiation) averts $4,082.7$ cases, incurring an expected loss/regret of $628.8$ units. In the semi-arid zone, peak monsoon rains occur in late August, meaning that a July initiation depletes drug titers before peak infectious biting.
- In **Zone 2 (Savanna)** and **Zone 3 (Humid Forest)**, longer wet seasons shift the maximum prevention window to **September 15** (cases averted $= 2,077.7$ in Zone 2 and $1,393.3$ in Zone 3). Early initiation in May or June incurs severe opportunity loss ($>1,180$ units in Zone 2).

---

## 9. Structural Uncertainty Quantification

To evaluate whether intervention timing recommendations depend on the functional form of larval habitat capacity, three alternative carrying-capacity models $K(R)$ were evaluated:
1. **Standard Saturation (Parham & Michael, 2010):** $K(R) = K_0 \frac{R}{R + K_{\text{half}}}$.
2. **Exponential Formulation:** $K(R) = K_0 (1 - \exp(-R / R_0))$.
3. **Linear Capped Formulation:** $K(R) = K_0 \min(1.0, R / R_{\text{max}})$.

Across all three functional forms (Table 8), the ranking of intervention start dates remained invariant: August 15 averted $4,711.5$ case units versus $4,082.7$ for July 15 in Zone 1 across all three specifications. This demonstrates that intervention timing optimality is governed by the temporal phasing of monsoonal precipitation rather than the mathematical formulation of larval density dependence.

---

## 10. Discussion

This study demonstrates the feasibility of combining formal numerical verification with Bayesian uncertainty quantification in spatial vector-host epidemiology:
- **Verified Numerics:** We observed an empirical spatial convergence order of $p = 2.011$ for the diffusion operator and confirmed human mass conservation to machine precision ($\le 3.020 \times 10^{-15}$).
- **Objective Zonation:** Unsupervised clustering with silhouette diagnostics partitioned the 37 units into three seasonal regimes (mean $s = 0.4590$), while explicitly flagging 11 boundary states exhibiting cluster disagreement.
- **Identifiability & Transferability:** Bayesian calibration recovered transmission efficiency parameters ($b_H, b_M$) within 95% HPDIs, while identifying structural shrinkage in human clearance rates ($\gamma_H$). Centroid calibrations transferred to uncalibrated states with consistent RMSE ($0.0829 \pm 0.0031$ in Zone 1), demonstrating spatial stability across homogenous climate clusters.

---

## 11. Limitations

1. **Synthetic Benchmark:** Calibration and validation were conducted on an in silico benchmark. While synthetic benchmarks allow rigorous assessment of parameter recovery, translational policy deployment requires final calibration against empirical health facility incidence (DHIS2) and survey micro-data (DHS/MIS).
2. **Spatial Dimension:** The model employs a 1D settlement-breeding transect. While this captures boundary dispersal gradients, it does not represent complex 2D hydrology or urban microclimates.
3. **Age Structure & Resistance:** The model assumes a single homogeneous host population and constant chemoprotective efficacy ($\eta_{\text{clear}} = 0.95$). Emergence of resistance mutations or age-dependent immunity could alter protective longevity.

---

## 12. Data and Code Availability

The complete computational pipeline—including downloaders, verification suites, MCMC calibration, validation scripts, and figure/table generators—is open-source under the MIT License at:  
`https://github.com/IsaacMoon1/malaria-pde-nigeria`.

The workflow executes via a single command:
```bash
python run_all.py
```

---

## Acknowledgements
The author acknowledges the high-performance computational facilities and open-source scientific software packages (NumPy, SciPy, pandas, Matplotlib, SALib, emcee) that supported this research.

---

## References

Bhatt, S., Weiss, D. J., Cameron, E., Bisanzio, D., Mappin, B., Dalrymple, U., Nguyen, T. D., Andolina, C., Smith, D. L., Horby, P., Hay, S. I., & Gething, P. W. (2015). The effect of malaria control on *Plasmodium falciparum* in Africa between 2000 and 2015. *Nature*, 526(7572), 207–211. https://doi.org/10.1038/nature15535

Berger, J. O. (1985). *Statistical Decision Theory and Bayesian Analysis* (2nd ed.). Springer-Verlag. https://doi.org/10.1007/978-1-4757-4286-2

Cairns, M., Roca-Feltrer, A., Garske, T., Wilson, A. L., Diallo, D., Milligan, P. J., Ghani, A. C., & Greenwood, B. M. (2012). Estimating the potential public health impact of seasonal malaria chemoprevention in African children. *Nature Communications*, 3, 881. https://doi.org/10.1038/ncomms1879

Depinay, J. M., Mbogo, C. M., Killeen, G. F., Knols, B. G., Beier, J. C., Carlson, J., Dushoff, J., Billingsley, P., & Mwambi, H. (2004). A simulation model of African *Anopheles* ecology and population dynamics for the analysis of malaria transmission. *Malaria Journal*, 3, 29. https://doi.org/10.1186/1475-2875-3-29

Foreman-Mackey, D., Hogg, D. W., Lang, D., & Goodman, J. (2013). emcee: The MCMC hammer. *Publications of the Astronomical Society of the Pacific*, 125(925), 306–312. https://doi.org/10.1086/670067

Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. *Statistical Science*, 7(4), 457–472. https://doi.org/10.1214/ss/1177011136

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). Chapman and Hall/CRC. https://doi.org/10.1201/b16018

Gneiting, T., & Raftery, A. E. (2007). Strictly proper scoring rules, prediction, and estimation. *Journal of the American Statistical Association*, 102(477), 359–378. https://doi.org/10.1198/016214506000001437

Greenwood, B. (2010). Anti-malarial drugs and the prevention of malaria in the population of malaria endemic areas. *Malaria Journal*, 9(Suppl 3), S2. https://doi.org/10.1186/1475-2875-9-S3-S2

Griffin, J. T., Hollingsworth, T. D., Okell, L. C., Churcher, T. S., White, M., Hinsley, W., Bousema, T., Drakeley, C. J., Esteva, L., Basáñez, M. G., & Ghani, A. C. (2010). Reducing *Plasmodium falciparum* malaria transmission in Africa: A model-based evaluation of intervention strategies. *PLoS Medicine*, 7(8), e1000324. https://doi.org/10.1371/journal.pmed.1000324

Herman, J., & Usher, W. (2017). SALib: An open-source Python library for sensitivity analysis. *Journal of Open Source Software*, 2(9), 97. https://doi.org/10.21105/joss.00097

Lou, Y., & Zhao, X. Q. (2010). A reaction-diffusion malaria model with incubation period and seasonal development. *Journal of Mathematical Biology*, 62(4), 543–568. https://doi.org/10.1007/s00285-010-0346-7

Meremikwu, M. M., Donegan, S., Esu, E., & Oringanje, C. (2012). Seasonal malaria chemoprevention for preventing malaria in children in areas of highly seasonal transmission. *Cochrane Database of Systematic Reviews*, 2, CD003756. https://doi.org/10.1002/14651858.CD003756.pub4

Mordecai, E. A., Paaijmans, K. P., Johnson, L. R., Balzer, C., Ben-Horin, T., de Moor, E., McNally, A., Pawar, S., Ryan, S. J., Smith, T. C., & Lafferty, K. D. (2013). Optimal temperature for malaria transmission is lower than previously predicted. *Ecology Letters*, 16(1), 22–30. https://doi.org/10.1111/ele.12015

Mordecai, E. A., Caldwell, J. M., Grossman, M. K., Lippi, C. A., Johnson, L. R., Neira, M., Rohr, J. R., Ryan, S. J., Savage, V., Shocket, M. S., Sippy, R., Stewart Ibarra, A. M., Thomas, M. B., & Villena, O. (2019). Thermal biology of mosquito-borne disease. *Ecology Letters*, 22(10), 1690–1708. https://doi.org/10.1111/ele.13335

NiMet. (2023). *2023 Seasonal Climate Prediction (SCP)*. Nigerian Meteorological Agency, Federal Ministry of Aviation, Abuja, Nigeria.

Okuneye, K., & Gumel, A. B. (2017). Analysis of a temperature- and rainfall-dependent model for malaria transmission dynamics. *Mathematical Biosciences*, 287, 72–92. https://doi.org/10.1016/j.mbs.2017.02.007

Parham, P. E., & Michael, E. (2010). Modeling the effects of weather and climate change on malaria transmission. *Environmental Health Perspectives*, 118(5), 620–626. https://doi.org/10.1289/ehp.0901256

Rousseeuw, P. J. (1987). Silhouettes: A graphical aid to the interpretation and validation of cluster analysis. *Journal of Computational and Applied Mathematics*, 20, 53–65. https://doi.org/10.1016/0377-0427(87)90125-7

Saltelli, A., Ratto, M., Andres, T., Campolongo, F., Cariboni, J., Gatelli, D., Saisana, M., & Tarantola, S. (2008). *Global Sensitivity Analysis: The Primer*. John Wiley & Sons. https://doi.org/10.1002/9780470725184

Service, M. W. (1993). *Mosquito Ecology: Field Sampling Methods* (2nd ed.). Chapman & Hall. https://doi.org/10.1007/978-94-015-8113-4

Smith, D. L., McKenzie, F. E., Snow, R. W., & Hay, S. I. (2007). Revisiting the basic reproductive number for malaria and its implications for malaria control. *PLoS Biology*, 5(3), e42. https://doi.org/10.1371/journal.pbio.0050042

Sobol, I. M. (2001). Global sensitivity indices for nonlinear mathematical models and their Monte Carlo estimates. *Mathematics and Computers in Simulation*, 55(1–3), 271–280. https://doi.org/10.1016/S0378-4754(00)00270-6

Sultan, B., & Janicot, S. (2003). The West African monsoon dynamics. Part II: The "pre-onset" and "onset" of the summer monsoon. *Journal of Climate*, 16(21), 3407–3427. https://doi.org/10.1175/1520-0442(2003)016<3407:TWAMDP>2.0.CO;2

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P. C. (2021). Rank-normalization, folding, and localization: An improved $\widehat{R}$ for assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667–718. https://doi.org/10.1214/20-BA1221

Ward, J. H. (1963). Hierarchical grouping to optimize an objective function. *Journal of the American Statistical Association*, 58(301), 236–244. https://doi.org/10.1080/01621459.1963.10500845

Weiss, D. J., Lucas, T. C., Nguyen, M., Nandi, A. K., Bisanzio, D., Battle, K. E., Cameron, E., Twohig, K. A., Pfeffer, D. A., Rozier, J. A., Gibson, H. S., Rao, P. C., Hedley, E., Collins, E. L., Simpson, N., Fleming, K. A., Howes, R. E., Gething, P. W., & Bhatt, S. (2019). Mapping the global prevalence, incidence, and mortality of *Plasmodium falciparum*, 2000–17: A spatial and temporal modelling study. *The Lancet*, 394(10195), 322–331. https://doi.org/10.1016/S0140-6736(19)31097-9

World Health Organization. (2022). *WHO Guidelines for Malaria - 3 June 2022*. World Health Organization, Geneva. https://apps.who.int/iris/handle/10665/354781

World Health Organization. (2023). *World Malaria Report 2023*. World Health Organization, Geneva. https://iris.who.int/handle/10665/374526

---

## Tables and Figures Reference

### Table 1: Biophysical and Epidemiological Parameter Registry
*Summary of 20 model parameters, biological descriptions, baseline values, units, feasible exploration intervals, and literature citations with authentic DOIs. (Saved at `tables/table1_parameters.csv` and `tables/table1_parameters.tex`).*

### Table 2: Climatological Clustering and Representative Centroid Metrics
*Summary of 37 administrative units partitioned into Zone 1 (Semi-Arid), Zone 2 (Savanna), and Zone 3 (Humid Forest), including mean annual precipitation, rainy season duration, mean temperature, silhouette scores, and centroid distances. (Saved at `tables/table2_climate_clustering.csv` and `tables/table2_climate_clustering.tex`).*

### Table 3: Summary of Numerical Verification Tests
*Formal verification suite results demonstrating analytic diffusion accuracy ($5.055 \times 10^{-6}$), empirical spatial convergence ($p = 2.011$), cross-solver concordance ($r = 1.0000$, $\text{NRMSE} = 0.01\%$), and machine-precision mass conservation ($\le 3.020 \times 10^{-15}$). (Saved at `tables/table3_verification_summary.csv` and `tables/table3_verification_summary.tex`).*

### Table 4: Global Sensitivity Analysis (Sobol Variance Decomposition)
*First-order ($S_1$) and total-order ($S_T$) Sobol sensitivity indices for infection prevalence across the three ecological zones. (Saved at `tables/table4_sobol_indices.csv` and `tables/table4_sobol_indices.tex`).*

### Table 5: Bayesian MCMC Calibration and Parameter Recovery Diagnostics
*Prior bounds, true data-generating values, posterior means, 95% highest posterior density intervals (HPDI), Gelman-Rubin $\hat{R}$, and effective sample sizes for calibrated centroid sites. (Saved at `tables/table5_mcmc_calibration_diagnostics.csv` and `tables/table5_mcmc_calibration_diagnostics.tex`).*

### Table 6: Spatiotemporal Validation and Spatial Transferability Performance
*Out-of-sample temporal holdout (94 months) metrics (RMSE, CRPS, 95% CI coverage) comparing calibrated centroid representative sites against 34 spatial transfer states. (Saved at `tables/table6_validation_transferability.csv` and `tables/table6_validation_transferability.tex`).*

### Table 7: Bayesian Decision Analysis and Intervention Timing Regret
*Evaluation of candidate SMC initiation windows, clinical cases averted index, and expected timing regret across seasonal regimes. (Saved at `tables/table7_decision_policy_timing.csv` and `tables/table7_decision_policy_timing.tex`).*

### Table 8: Structural Model Uncertainty Under Alternative Carrying Capacities
*Comparison of cases averted across standard saturation, exponential, and linear-capped carrying-capacity formulations, demonstrating structural invariance of optimal timing windows. (Saved at `tables/table8_structural_uncertainty.csv` and `tables/table8_structural_uncertainty.tex`).*

---

### Figure Captions

- **Figure 1: Climatological Classification and Spatial Distribution of the 37 Benchmark Units.** Map displaying objective partitioning into Zone 1 (Semi-Arid, gold), Zone 2 (Savanna, green), and Zone 3 (Humid Forest, blue), with centroid representative sites (Jigawa, FCT, Abia) highlighted alongside transitional boundary states. (File at `figures/fig1_nigeria_climate_zones.png`).
- **Figure 2: Numerical Verification, Spatial Order of Accuracy, and Cross-Solver Concordance.** Log-log plot of spatial discretization error versus grid spacing establishing $p = 2.011$ convergence order, alongside temporal trajectory alignment between the primary method-of-lines solver and the independent Runge-Kutta solver ($r = 1.00000$). (File at `figures/fig2_verification_convergence.png`).
- **Figure 3: Global Variance-Based Sensitivity Analysis (Sobol Indices).** Bar plots of first-order ($S_1$) and total-order ($S_T$) Sobol indices illustrating the predominant influence of transmission efficiencies ($b_H, b_M$) and vector diffusion ($D_M$) across all three eco-epidemiological zones. (File at `figures/fig3_sobol_sensitivity.png`).
- **Figure 4: Bayesian MCMC Calibration Traces and Marginal Posterior Distributions.** Ensemble MCMC Markov chains and posterior density marginals for $D_M, b_H, b_M, \gamma_H$, showing parameter recovery against known ground truth. (File at `figures/fig4_mcmc_traces_posterior.png`).
- **Figure 5: Out-of-Sample Spatiotemporal Validation and Spatial Transferability.** Multi-panel time-series validation curves comparing model posterior predictive credible bands against observed longitudinal prevalence over the temporal holdout, contrasting calibrated centroid sites with spatial transfer states. (File at `figures/fig5_spatiotemporal_validation.png`).
- **Figure 6: Bayesian Decision Analysis: Optimal Intervention Timing and Opportunity Regret.** Expected cases averted curves and opportunity regret distributions across candidate SMC deployment windows, illustrating the superiority of August 15 initiation in Zone 1. (File at `figures/fig6_decision_timing_regret.png`).
