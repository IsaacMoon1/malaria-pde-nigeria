"""
Epidemiological, entomological, and intervention parameters compiled from literature.
Every value extracted from published literature is explicitly tagged as:
'LITERATURE VALUE: verify before citing' along with its formal citation and DOI.
Target: Elsevier Journal standards.
"""

import pandas as pd

# Comprehensive Parameter Registry
PARAMETERS_TABLE = [
    {
        "symbol": "D_M",
        "description": "Mosquito spatial diffusion coefficient",
        "value": 0.085,
        "units": "km^2 / day",
        "bounds": [0.01, 0.25],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Service, M. W. (1997). Mosquito ecology: field sampling methods. Chapman & Hall; Okuneye, K., & Gumel, A. B. (2017). Analysis of a temperature- and rainfall-dependent model for malaria transmission dynamics. Mathematical Biosciences, 287, 72-92.",
        "doi": "10.1016/j.mbs.2017.02.007"
    },
    {
        "symbol": "T_min_larva",
        "description": "Minimum threshold temperature for larval development",
        "value": 16.0,
        "units": "°C",
        "bounds": [14.0, 18.0],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Mordecai, E. A. et al. (2013). Optimal temperature for malaria transmission is lower than previously predicted. Ecology Letters, 16(1), 22-30.",
        "doi": "10.1111/ele.12015"
    },
    {
        "symbol": "T_opt_larva",
        "description": "Optimum temperature for larval development rate",
        "value": 28.0,
        "units": "°C",
        "bounds": [26.0, 31.0],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Mordecai, E. A. et al. (2019). Thermal biology of mosquito-borne disease. Ecology Letters, 22(10), 1690-1708.",
        "doi": "10.1111/ele.13335"
    },
    {
        "symbol": "T_max_larva",
        "description": "Maximum thermal threshold for larval survival",
        "value": 35.0,
        "units": "°C",
        "bounds": [33.0, 38.0],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Mordecai, E. A. et al. (2013). Optimal temperature for malaria transmission is lower than previously predicted. Ecology Letters, 16(1), 22-30.",
        "doi": "10.1111/ele.12015"
    },
    {
        "symbol": "d_L_baseline",
        "description": "Baseline aquatic larval mortality rate",
        "value": 0.10,
        "units": "day^-1",
        "bounds": [0.05, 0.25],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Depinay, J. M. et al. (2004). A simulation model of African Anopheles ecology and population dynamics for the analysis of malaria transmission. Acta Tropica, 89(1), 1-19.",
        "doi": "10.1016/j.actatropica.2003.09.006"
    },
    {
        "symbol": "gamma_L_baseline",
        "description": "Larval pupation and emergence rate at 25°C",
        "value": 0.08,
        "units": "day^-1",
        "bounds": [0.04, 0.15],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Parham, P. E., & Michael, E. (2010). Modeling the effects of weather and climate change on malaria transmission. Environmental Health Perspectives, 118(5), 620-626.",
        "doi": "10.1289/ehp.0901256"
    },
    {
        "symbol": "beta_E",
        "description": "Average oviposition rate per adult female",
        "value": 30.0,
        "units": "eggs / mosquito / day",
        "bounds": [15.0, 60.0],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Okuneye, K., & Gumel, A. B. (2017). Analysis of a temperature- and rainfall-dependent model for malaria transmission dynamics. Mathematical Biosciences, 287, 72-92.",
        "doi": "10.1016/j.mbs.2017.02.007"
    },
    {
        "symbol": "a_biting_rate",
        "description": "Human-biting frequency per adult mosquito",
        "value": 0.33,
        "units": "bites / mosquito / day",
        "bounds": [0.20, 0.50],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Griffin, J. T. et al. (2010). Reducing Plasmodium falciparum malaria transmission in Africa: a model-based evaluation of intervention strategies. PLoS Medicine, 7(8), e1000324.",
        "doi": "10.1371/journal.pmed.1000324"
    },
    {
        "symbol": "b_H",
        "description": "Transmission probability per infectious bite to human",
        "value": 0.50,
        "units": "dimensionless",
        "bounds": [0.25, 0.80],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Griffin, J. T. et al. (2014). Potential for elimination of Plasmodium falciparum malaria with contemporary interventions: a modeling study. Lancet Global Health, 2(4), e222-e233.",
        "doi": "10.1016/S2214-109X(14)70198-4"
    },
    {
        "symbol": "b_M",
        "description": "Transmission probability from infectious human to susceptible mosquito",
        "value": 0.45,
        "units": "dimensionless",
        "bounds": [0.20, 0.70],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Griffin, J. T. et al. (2010). Reducing Plasmodium falciparum malaria transmission in Africa: a model-based evaluation of intervention strategies. PLoS Medicine, 7(8), e1000324.",
        "doi": "10.1371/journal.pmed.1000324"
    },
    {
        "symbol": "mu_M_baseline",
        "description": "Adult mosquito baseline natural mortality rate at 25°C",
        "value": 0.10,
        "units": "day^-1",
        "bounds": [0.06, 0.20],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Mordecai, E. A. et al. (2019). Thermal biology of mosquito-borne disease. Ecology Letters, 22(10), 1690-1708.",
        "doi": "10.1111/ele.13335"
    },
    {
        "symbol": "DD_parasite",
        "description": "Degree-days required for parasite sporogony (EIP)",
        "value": 111.0,
        "units": "degree-days",
        "bounds": [100.0, 125.0],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Detinova, T. S. (1962). Age-grouping methods in Diptera of medical importance. WHO Monograph Series, 47; Mordecai et al. (2013).",
        "doi": "10.1111/ele.12015"
    },
    {
        "symbol": "T_min_sporogony",
        "description": "Minimum thermal threshold for P. falciparum development",
        "value": 16.0,
        "units": "°C",
        "bounds": [14.5, 17.5],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Mordecai, E. A. et al. (2013). Optimal temperature for malaria transmission is lower than previously predicted. Ecology Letters, 16(1), 22-30.",
        "doi": "10.1111/ele.12015"
    },
    {
        "symbol": "alpha_H",
        "description": "Human intrinsic incubation rate (1 / 14 days)",
        "value": 0.0714,
        "units": "day^-1",
        "bounds": [0.05, 0.10],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "White, M. T. et al. (2015). Modelling the overt and hidden burden of malaria. Parasites & Vectors, 8, 28.",
        "doi": "10.1186/s13071-015-0638-4"
    },
    {
        "symbol": "gamma_H",
        "description": "Human recovery rate from clinical malaria",
        "value": 0.0167,
        "units": "day^-1",
        "bounds": [0.008, 0.035],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Griffin, J. T. et al. (2010). Reducing Plasmodium falciparum malaria transmission in Africa: a model-based evaluation of intervention strategies. PLoS Medicine, 7(8), e1000324.",
        "doi": "10.1371/journal.pmed.1000324"
    },
    {
        "symbol": "omega_H",
        "description": "Loss of naturally acquired clinical immunity",
        "value": 0.0027,
        "units": "day^-1",
        "bounds": [0.001, 0.005],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Griffin, J. T. et al. (2014). Potential for elimination of Plasmodium falciparum malaria with contemporary interventions: a modeling study. Lancet Global Health, 2(4), e222-e233.",
        "doi": "10.1016/S2214-109X(14)70198-4"
    },
    {
        "symbol": "mu_H",
        "description": "Natural human daily mortality rate (approx. 55-year life expectancy)",
        "value": 0.00005,
        "units": "day^-1",
        "bounds": [0.00004, 0.00007],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "World Health Organization. (2024). World Malaria Report 2024. Geneva: WHO.",
        "doi": "ISBN:978-92-4-009999-0"
    },
    {
        "symbol": "psi_SMC",
        "description": "Decay rate of chemoprotection per SMC cycle (~28 days efficacy)",
        "value": 0.0357,
        "units": "day^-1",
        "bounds": [0.025, 0.050],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Cairns, M. et al. (2012). Estimating the potential public health impact of seasonal malaria chemoprevention in African children. Nature Communications, 3, 881.",
        "doi": "10.1038/ncomms1879"
    },
    {
        "symbol": "rho_SMC_efficacy",
        "description": "Clinical protective efficacy of SMC (SPAQ regimen)",
        "value": 0.88,
        "units": "dimensionless",
        "bounds": [0.75, 0.95],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "World Health Organization. (2022). WHO Guidelines for malaria. Geneva: World Health Organization.",
        "doi": "WHO/UCN/GMP/2022.01"
    },
    {
        "symbol": "delta_IRS_mortality",
        "description": "IRS-induced adult vector mortality increase",
        "value": 0.25,
        "units": "day^-1",
        "bounds": [0.10, 0.45],
        "tag": "LITERATURE VALUE: verify before citing",
        "citation": "Pluess, B. et al. (2010). Indoor residual spraying for preventing malaria. Cochrane Database of Systematic Reviews, 2010(4), CD006657.",
        "doi": "10.1002/14651858.CD006657.pub2"
    }
]

def get_parameters_dataframe():
    """Return parameters formatted as pandas DataFrame for validation and paper export."""
    return pd.DataFrame(PARAMETERS_TABLE)

def get_default_parameter_dict():
    """Return key-value dict of default parameter values."""
    return {p["symbol"]: p["value"] for p in PARAMETERS_TABLE}
