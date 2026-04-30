# US Transportation Noise Health Impact Assessment

Analysis code for: **"Health burden of transportation noise in the United States: a national assessment with equity analysis"**

Published in *Environmental Research* (ER-26-4325).

**Author:** David Rojas-Rueda, Colorado State University

## Overview

This repository contains the Monte Carlo uncertainty analysis code and input data used to generate the results reported in the manuscript. The simulation propagates uncertainty through 12 parameters simultaneously across both YLD (years lived with disability) and YLL (years of life lost) components.

## Requirements

- Python 3.8+
- numpy
- scipy

```bash
pip install numpy scipy
```

## Usage

```bash
# Default: 10,000 iterations, seed 42
python us_noise_hia_monte_carlo.py

# Custom iterations and seed
python us_noise_hia_monte_carlo.py 10000 42
```

Results are saved to `mc_results.json` and printed to stdout.

## Files

| File | Description |
|------|-------------|
| `us_noise_hia_monte_carlo.py` | 12-parameter Monte Carlo simulation script |
| `us_national_hia_input.csv` | National noise exposure distribution (6 bands, population counts) |
| `LICENSE` | MIT License |

## Parameters Varied

| # | Parameter | Distribution | Central value | Range |
|---|-----------|-------------|---------------|-------|
| 1 | RR IHD | Normal (log) | 1.08 | 1.01-1.15 |
| 2 | RR ischemic stroke | Normal (log) | 1.025 | 1.009-1.041 |
| 3 | RR cardiomyopathy | Normal (log) | 1.04 | 1.02-1.07 |
| 4 | RR CVD mortality | Normal (log) | 1.05 | 1.02-1.07 |
| 5 | LAeq-to-Lden offset | PERT | +3 dB | 0-5 dB |
| 6 | Lden-to-Lnight offset | PERT | -9 dB | -11 to -7 dB |
| 7 | DW annoyance | PERT | 0.011 | 0.005-0.025 |
| 8 | DW sleep disturbance | PERT | 0.010 | 0.006-0.07 |
| 9 | DW IHD | PERT | 0.08 | 0.05-0.12 |
| 10 | IHD episode duration | PERT | 5 yr | 3-8 yr |
| 11 | Road noise fraction | Uniform | 0.8 (mean) | 0.6-1.0 |
| 12 | Remaining life expectancy | PERT | 14.5 yr | 12-17 yr |

## Expected Output

With default parameters (10,000 iterations, seed 42):

| Metric | Value |
|--------|-------|
| Total DALYs (deterministic) | 272,578 |
| Total DALYs (MC median) | 245,486 |
| 95% UI | 144,849 - 393,815 |
| MCSE (median) | 815 |
| CVM | 0.041% |

## Data Sources

Input data are publicly available from:

- **Noise exposure:** National Transportation Noise Exposure Map ([download](https://deohs.washington.edu/national-transportation-noise-exposure-map-download))
- **Baseline disease rates:** Global Burden of Disease Study 2023 ([GBD Results](https://vizhub.healthdata.org/gbd-results/))
- **Social vulnerability:** CDC/ATSDR Social Vulnerability Index 2022 ([SVI](https://www.atsdr.cdc.gov/placeandhealth/svi/))

## Exposure-Response Functions

| Outcome | Source |
|---------|--------|
| IHD incidence | van Kempen et al. 2018 (WHO Environmental Noise Guidelines) |
| Ischemic stroke | Pershagen et al. 2025 (Environ Epidemiol) |
| Cardiomyopathy | Sorensen et al. 2024 (Redox Biol) |
| CVD mortality | Sorensen et al. 2024 (Redox Biol) |
| Type 2 diabetes | Vienneau et al. 2024 (Environ Health) |
| High annoyance | Guski et al. 2017 (WHO systematic review) |
| High sleep disturbance | Smith et al. 2022 (Environ Health Perspect) |
| Disability weights | Charalampous et al. 2024 (BMJ Public Health) |

## Citation

Rojas-Rueda D. Health burden of transportation noise in the United States: a national assessment with equity analysis. *Environmental Research.* 2026.

## License

MIT License
