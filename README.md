# Ganga River Mainstem: Environmental Science & Pressure Analysis

> **Course / Project Context:** Environmental Science (ES) Midsem Research Project  
> **Core Research Question:** *How have environmental pressures around culturally or religiously important rivers changed over time?*  
> A spatial, temporal, demographic, and hydrological investigation of the Ganga River mainstem (2012–2024) across water quality observations, population dynamics, hydrological assimilative capacity, and mega religious gathering events.

---

## 📌 Executive Summary & Current Project Status

This repository contains the complete analytical pipeline, processed datasets, statistical tests, and publication-ready figures for our Environmental Science study on the Ganga River mainstem. 

We examine the longitudinal environmental gradient of the river across **5 strategically selected monitoring stations** stretching from the Himalayan foothills to the lower Gangetic plain:
1. **Haridwar (`Har-Ki-Pauri Ghat`)** – Upper baseline stretch / foothill entry
2. **Kanpur (`Jajmau Bridge`)** – Major industrial and tannery cluster
3. **Prayagraj (`Sangam`)** – Yamuna-Ganga confluence & Maha Kumbh Mela mega-gathering hub
4. **Varanasi (`Malviya Bridge`)** – Dense cultural, cremation, and pilgrimage center
5. **Patna (`Gandhi Ghat`)** – High-density lower alluvial plain

---

## 🚀 What Has Been Done So Far

### 1. Data Collection & Multi-Source Inventory (`docs/` & `data/raw/`)
- **CPCB NWMP Water Quality Reports (2012–2024):** Acquired annual national water quality data published by the Central Pollution Control Board (CPCB) covering DO, BOD, pH, conductivity, nitrate, and fecal/total coliform.
- **Hydrography & Basin Architecture (HydroSHEDS / HydroRIVERS):** Obtained high-resolution HydroBASINS (Level 1–12) and HydroRIVERS Asia vector reach datasets, extracting reach geometry and long-term modeled discharge (`flow_avg_cms`).
  - *Data Policy Note:* Observed gauge discharge on the India-WRIS portal is classified/restricted by the Central Water Commission (CWC) due to transboundary river sensitivities. Modeled discharge from HydroRIVERS was utilized as an objective hydrodynamic baseline.
- **Demographic Rasters (WorldPop Project):** Downloaded 1 km resolution gridded human population count and density GeoTIFFs for India across 5-year intervals (2000, 2005, 2010, 2015, 2020).
- **Religious Events Chronicle (`events.csv`):** Manually compiled metadata, estimated footfall (up to 120 million bathers), dates, and official source links for mega mass bathing events (Maha Kumbh 2013, Ardh Kumbh 2019, Haridwar Kumbh 2021, Chhath Puja, Ganga Dussehra, and Kartik Purnima).
- **Monitoring Station Gazetteer (`docs/stations.csv`):** Standardized GPS coordinates, town names, administrative states, and parameter lists across all 5 key monitoring points.

### 2. Data Cleaning, Parsing & Spatial Processing (`notebooks/` & `data/processed/`)
- **PDF Extraction & Format Harmonization:** Built custom parsers and regex patterns to normalize shifting annual report structures:
  - Harmonized 2012–2014 triplets (Min, Max, Mean) with 2015+ pairs (Min, Max).
  - Addressed multi-line line-wrapping, missing parameter columns (e.g., Haridwar & Patna 2017/2018 nitrate records), and eliminated unanchored regex false positives.
- **Demographic Linear Interpolation (2012–2024):** Because WorldPop data is released on a 5-year discrete interval (and official Census 2021 was delayed), we implemented linear interpolation to generate continuous annual population density estimates within a 1 km station buffer.
- **GIS Reach Extraction:** Clipped and aligned the Ganga river mainstem network (`ganga_rivers_clipped.shp`) and paired each station with upstream catchment area and average river discharge (`station_avg_flow.csv`).
- **Master Dataset Integration:** Merged all four domains into a unified analytical matrix:
  👉 `data/processed/wq_es_integrated_master.csv`

### 3. Statistical Testing & Environmental Hypothesis Evaluation
- **Spearman Rank Correlation Analyses:**
  - *Population Density vs BOD:* Evaluated anthropogenic organic loading vs downstream tributary dilution dynamics ($\rho = -0.52, p < 0.001$).
  - *Population Density vs Fecal Coliform:* Strong positive correlation linking municipal raw sewage discharge with bacterial contamination ($\rho = +0.32, p = 0.027$).
  - *River Discharge vs BOD Concentration:* Evaluated hydrodynamic assimilative capacity and tributary buffering.
- **Non-Parametric Event Testing (Mann-Whitney U Test):**
  - Compared Prayagraj Sangam Kumbh event years (2013, 2019) against baseline non-event years to examine the signature of acute episodic gathering shocks within annual monitoring envelopes.
- **CPCB Class B (Outdoor Bathing) Standard Compliance:**
  - Systematically evaluated each station against official Indian regulatory thresholds:
    - Dissolved Oxygen (DO) $\ge 5.0\text{ mg/L}$
    - Biochemical Oxygen Demand (BOD) $\le 3.0\text{ mg/L}$
    - Fecal Coliform $\le 2500\text{ MPN/100 mL}$

### 4. Cause-and-Effect Matrix & Engineering Solutions
- Formulated an Environmental Science synthesis connecting observed water quality degradation to:
  - **Streeter-Phelps deoxygenation kinetics** (Kanpur/Sangam DO sag)
  - **Pathogen proliferation** (untreated nala outfalls in Varanasi and Patna)
  - **Hydrodynamic assimilative capacity** (Yamuna, Ghaghara, Gandak, Son tributary inflows)
  - **Policy & engineering interventions** (Namami Gange STP commissioning, Sisamau Nala interception & diversion in 2018, mandatory environmental flow / e-flow policies).

### 5. Automated End-to-End Pipeline (`src/midsem_es_analysis.py`)
- Engineered a modular, self-contained Python pipeline that:
  - Ingests raw data extracts
  - Cleans water quality records
  - Interpolates demographics
  - Merges hydrology and event metadata
  - Computes all statistical tests and outputs CSV tables
  - Generates 5 publication-ready charts (300 DPI)

### 6. Enhanced Analysis Module (`src/enhanced_analysis.py`)
- **Water Quality Index (WQI):** Implemented NSF-WQI adapted methodology with weighted sub-indices (DO 31%, Fecal Coliform 30%, BOD 23%, pH 16%) to produce a single composite score (0–100) per station per year.
- **Multi-Dimensional Pressure Profiling:** Constructed radar/spider chart overlays comparing 6 normalized pressure axes (BOD Load, DO Deficit, Pathogen Risk, Population Pressure, Dilution Deficit, Nutrient Load) across all 5 stations.
- **Namami Gange Policy Impact Assessment:** Pre-vs-Post (2012–2015 vs 2017–2024) paired Mann-Whitney U tests with percentage change analysis for BOD, DO, and Fecal Coliform at every station.
- **Coliform Exceedance Ratio Analysis:** Computed the ratio of observed fecal coliform to CPCB safe limit (2500 MPN/100mL) showing each station's departure from safety over time.
- **Regulatory Compliance Dashboard:** Traffic-light heatmap showing % compliance across DO, BOD, Fecal Coliform, and pH parameters per station.
- Generates 5 additional publication-ready figures (Fig 6–10) and 2 new statistical tables (Table 4–5).

---

## 📊 Summary of Key Findings

| Monitoring Station | Town | Modeled Flow ($m^3/s$) | Mean Pop. Density ($persons/km^2$) | DO Min ($mg/L$) | BOD Max ($mg/L$) | Peak Fecal Coliform ($MPN/100mL$) | CPCB Class B Bathing Compliance |
|---|---|---|---|---|---|---|---|
| **Har-Ki-Pauri** | Haridwar | 566.1 | 4,676 | 7.1 | 1.04 | 170 | **Complies** |
| **Jajmau Bridge** | Kanpur | 998.9 | 2,156 | 2.8 | 6.30 | 17,000 | **Fails** *(DO Deficit, BOD Excess, Coliform)* |
| **Sangam** | Prayagraj | 3,061.9 | 4,551 | 4.7 | 9.20 | 33,000 | **Fails** *(DO Deficit, BOD Excess, Coliform)* |
| **Malviya Bridge** | Varanasi | 3,266.8 | 3,850 | 6.1 | 6.10 | 49,000 | **Fails** *(BOD Excess, Severe Coliform)* |
| **Gandhi Ghat** | Patna | 10,040.1 | 5,937 | 6.5 | 3.00 | 1,700,000 | **Fails** *(Coliform Extreme, up to 680× limit)* |

### Key Scientific Takeaways:
1. **The Upper Stretch Buffer:** Haridwar remains the only consistently compliant station due to high foothill reaeration, low water temperatures, and minimal upstream industrial discharge.
2. **Industrial & Confluence Sags:** Kanpur and Prayagraj experience severe DO dips below $3.0\text{ mg/L}$ and $5.0\text{ mg/L}$, driven by high carbonaceous/nitrogenous organic waste discharge.
3. **The Bacterial Crisis:** While BOD shows stabilization post-2018 (assisted by the Sisamau Nala interception in Kanpur and STP capacity increases under Namami Gange), fecal coliform levels in Varanasi and Patna exceed safe bathing limits by 20× to 680×, revealing an acute municipal sanitation and septage deficit.
4. **Hydrological Assimilative Paradox:** Discharge expands 18-fold from Haridwar ($566\text{ m}^3\text{/s}$) to Patna ($10,040\text{ m}^3\text{/s}$). High tributary volume dilutes BOD concentration at Patna despite peak population density, but pathogen persistence keeps bacterial contamination dangerously elevated.
5. **Episodic Festival Shock Loading:** Millions of bathers during Kumbh Mela cause localized acute surges; however, annual minimum/maximum monitoring envelopes tend to smooth out episodic spikes, highlighting the critical need for continuous real-time water quality monitoring during festivals.

---

## 📁 Repository Directory Structure

```text
Ganga-River-stem_ES-Project/
│
├── README.md                                  <- Project overview, status, and instructions (this file)
├── test_parse.py                              <- Regex validation & parser test suite
├── test_parse_complete.py                     <- Parser verification for all 5 monitoring stations
├── test_parse_details.py                      <- Parameter verification across all annual records
│
├── src/
│   └── midsem_es_analysis.py                  <- Master executable pipeline (clean -> merge -> stats -> plots)
│
├── docs/
│   ├── stations.csv                           <- GPS coordinates and metadata for 5 key monitoring stations
│   ├── data_inventory.csv                     <- Detailed metadata inventory & data limitations log
│   ├── cleaning_log.md                        <- Chronological log of data transformations
│   └── Download_log.txt                       <- Data sourcing and policy restriction notes
│
├── data/
│   ├── raw/
│   │   ├── cpcb/                              <- Raw CPCB water quality reports / tables
│   │   ├── hydrosheds/                        <- HydroBASINS & HydroRIVERS Asia vector files
│   │   ├── worldpop/                          <- 1 km population density rasters (2000-2020)
│   │   ├── events/                            <- events.csv (Kumbh Melas, Chhath Puja, etc.)
│   │   ├── indiawris/                         <- India-WRIS metadata notes
│   │   └── osm/                               <- OpenStreetMap extract configs
│   │
│   └── processed/
│       ├── wq_annual_clean.csv                <- Normalized CPCB water quality records (2012-2024)
│       ├── station_population_trends.csv      <- Extracted WorldPop densities at station buffers
│       ├── station_population_interpolated.csv<- Annual interpolated population densities (2012-2024)
│       ├── station_avg_flow.csv               <- HydroRIVERS reach discharge mapped to stations
│       ├── ganga_rivers_clipped.*             <- Clipped Ganga basin river network shapefiles
│       └── wq_es_integrated_master.csv        <- Master unified dataset for all analyses
│
├── outputs/
│   ├── figures/                               <- Publication-ready figures (300 DPI)
│   │   ├── fig1_spatial_water_quality_gradient.png  <- DO, BOD, Coliform transects
│   │   ├── fig2_temporal_trends_policy.png          <- 2012-2024 longitudinal trends & policy markers
│   │   ├── fig3_demographics_vs_pollution.png       <- Population density vs BOD & Fecal Coliform
│   │   ├── fig4_hydrology_assimilative_capacity.png  <- Flow dilution & organic mass flux
│   │   └── fig5_religious_events_impact.png         <- Kumbh Mela vs baseline comparisons
│   │
│   └── tables/                                <- Analytical summary & statistical CSV tables
│       ├── table1_station_baseline.csv        <- Baseline values & CPCB Class B compliance
│       ├── table2_statistical_tests.csv       <- Spearman & Mann-Whitney test results
│       └── table3_environmental_cause_effect_mitigation.csv <- ES synthesis matrix
│
└── notebooks/
    ├── 01_clean.ipynb                         <- Exploratory cleaning notebook
    ├── clean_1.ipynb                          <- Secondary cleaning workflow
    ├── data_clean_process.ipynb               <- Detailed data parsing & processing notebook
    └── ES Project.ipynb                       <- Initial project formulation
```

---

## 🖼️ Deliverables Overview

### Generated Figures (`outputs/figures/`)
- **Figure 1 (`fig1_spatial_water_quality_gradient.png`):** Spatial boxplots showing upstream-to-downstream gradients for DO minima, BOD maxima, and Fecal Coliform (logarithmic scale) against CPCB Class B thresholds.
- **Figure 2 (`fig2_temporal_trends_policy.png`):** Multi-year longitudinal trajectories (2012–2024) displaying BOD reductions and DO improvements following the Namami Gange launch (2014) and Sisamau Nala diversion (2018).
- **Figure 3 (`fig3_demographics_vs_pollution.png`):** Scatter plots and Spearman correlation fits linking riparian population density with organic load (BOD) and pathogen counts (Fecal Coliform).
- **Figure 4 (`fig4_hydrology_assimilative_capacity.png`):** Modeled river discharge vs BOD concentration highlighting dilution capacity, alongside estimated total organic mass flux ($BOD \times Q$).
- **Figure 5 (`fig5_religious_events_impact.png`):** Boxplot and strip plot distributions comparing Kumbh Mela event years against non-event baseline years at Prayagraj Sangam and Haridwar.
- **Figure 6 (`fig6_water_quality_index_heatmap.png`):** Station × Year heatmap showing composite Water Quality Index (WQI) scores using NSF-WQI adapted methodology with weighted sub-indices for DO, Fecal Coliform, BOD, and pH.
- **Figure 7 (`fig7_environmental_pressure_radar.png`):** Radar/spider chart comparing multi-dimensional environmental pressure profiles (BOD Load, DO Deficit, Pathogen Risk, Population Pressure, Dilution Deficit, Nutrient Load) across all 5 stations.
- **Figure 8 (`fig8_namami_gange_policy_impact.png`):** Pre-vs-Post Namami Gange (2014) paired boxplot comparisons for BOD, DO, and Fecal Coliform at each station.
- **Figure 9 (`fig9_coliform_exceedance_analysis.png`):** Fecal coliform exceedance ratio timeline and average exceedance bar chart showing departure from CPCB 2500 MPN safe limit.
- **Figure 10 (`fig10_compliance_dashboard.png`):** Traffic-light compliance dashboard showing percentage of monitoring years meeting each CPCB Class B parameter threshold.

### Analytical Tables (`outputs/tables/`)
- **Table 1 (`table1_station_baseline.csv`):** Comprehensive baseline metrics, discharge values, mean population densities, and formal regulatory compliance verdicts.
- **Table 2 (`table2_statistical_tests.csv`):** Formal hypothesis tests, Spearman $\rho$, Mann-Whitney $U$, sample sizes, and p-values.
- **Table 3 (`table3_environmental_cause_effect_mitigation.csv`):** Environmental science cause-and-effect matrix linking observed parameters to microbial processes, ecological/public health impacts, and targeted engineering mitigations.
- **Table 4 (`table4_water_quality_index.csv`):** Composite WQI scores (mean, min, max, latest) and quality category classification per station.
- **Table 5 (`table5_namami_gange_impact.csv`):** Pre-vs-Post Namami Gange statistical comparison with Mann-Whitney U test results, percentage change, and significance verdicts per station per parameter.

---

## 🛠️ How to Run the Pipeline

### Prerequisites
Ensure you have Python 3.8+ installed with the following packages:
```bash
pip install numpy pandas scipy matplotlib seaborn
```

### Reproduce Analysis and Figures
From the repository root directory, run:
```bash
# Step 1: Run core pipeline (Figures 1–5, Tables 1–3)
python3 src/midsem_es_analysis.py

# Step 2: Run enhanced analysis (Figures 6–10, Tables 4–5)
python3 src/enhanced_analysis.py
```

**Core Pipeline (Step 1)** executes all 5 stages:
1. Cleans CPCB annual water quality data $\to$ `data/processed/wq_annual_clean.csv`
2. Generates annual population interpolations $\to$ `data/processed/station_population_interpolated.csv`
3. Integrates all domains $\to$ `data/processed/wq_es_integrated_master.csv`
4. Runs statistical tests and writes summary tables $\to$ `outputs/tables/`
5. Plots and saves Figures 1–5 at 300 DPI $\to$ `outputs/figures/`

**Enhanced Analysis (Step 2)** generates:
6. Water Quality Index (WQI) heatmap $\to$ `fig6`
7. Multi-dimensional pressure radar chart $\to$ `fig7`
8. Namami Gange policy impact assessment $\to$ `fig8`
9. Coliform exceedance ratio analysis $\to$ `fig9`
10. Regulatory compliance dashboard $\to$ `fig10`

---

## 👥 Next Steps & Collaboration Roadmap for Teammates

Here are the remaining tasks and areas where teammates can jump in:
- **Presentation Deck Assembly:** Incorporate the 5 generated charts from `outputs/figures/` and tables from `outputs/tables/` into our midsem presentation slides.
- **OpenStreetMap (OSM) Riparian Mapping:** Use the Overpass API to extract and count cremation grounds, ghat geometries, and stormwater nala mouths in a 2 km buffer around each station.
- **Midsem Report Drafting:** Use the cause-and-effect matrix (`outputs/tables/table3_environmental_cause_effect_mitigation.csv`) to draft the environmental science narrative and engineering policy recommendations.
- **Seasonal Analysis:** If daily or seasonal pre-/post-monsoon data is obtained, analyze wet-season dilution vs dry-season lean flow concentration effects.

---

*Maintained for the Environmental Science Course Project. Feel free to open issues or submit PRs!*
