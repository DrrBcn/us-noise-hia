# US Transportation Noise Health Impact Assessment — Monte Carlo Code

Analysis code for: **"Health burden of transportation noise in the United States: a national assessment with equity analysis"**

Published in *Environmental Research* (ER-26-4325).

## Overview

This repository contains the Monte Carlo uncertainty analysis code used to generate the 95% uncertainty intervals reported in the manuscript. The simulation propagates uncertainty through 12 parameters simultaneously across both YLD (years lived with disability) and YLL (years of life lost) components.

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

## Parameters Varied

| # | Parameter | Distribution | Values |
|---|-----------|-------------|--------|
| 1 | RR IHD | Normal(log) | 1.08 (1.01–1.15) |
| 2 | RR ischemic stroke | Normal(log) | 1.025 (1.009–1.041) |
| 3 | RR cardiomyopathy / heart failure | Normal(log) | 1.04 (1.02–1.07) |
| 4 | RR CVD mortality | Normal(log) | 1.05 (1.02–1.07) |
| 5 | LAeq-to-Lden offset | PERT | (0, 3, 5) dB |
| 6 | Lden-to-Lnight offset | PERT | (-11, -9, -7) dB |
| 7 | DW annoyance | PERT | (0.005, 0.011, 0.025) |
| 8 | DW sleep disturbance | PERT | (0.006, 0.010, 0.07) |
| 9 | DW IHD | PERT | (0.05, 0.08, 0.12) |
| 10 | IHD episode duration | PERT | (3, 5, 8) years |
| 11 | Road noise fraction | Uniform | (0.6, 1.0) |
| 12 | Remaining life expectancy | PERT | (12, 14.5, 17) years |

## Data Sources

Input data are publicly available from:

- **Noise exposure:** National Transportation Noise Exposure Map ([download](https://deohs.washington.edu/national-transportation-noise-exposure-map-download))
- **Baseline disease rates:** Global Burden of Disease Study 2023 ([GBD Results](https://vizhub.healthdata.org/gbd-results/))
- **Social vulnerability:** CDC/ATSDR Social Vulnerability Index 2022 ([SVI](https://www.atsdr.cdc.gov/placeandhealth/svi/))

## Exposure-Response Functions

- **IHD:** van Kempen et al. 2018 (WHO Environmental Noise Guidelines), Int J Environ Res Public Health 15(2):379
- **Stroke:** Pershagen et al. 2025, Environ Epidemiol 9(3):e400
- **Cardiomyopathy / heart failure:** Engelmann et al. 2023, ETC/HE Report 2023/11
- **CVD mortality:** Munzel et al. 2024, Circ Res 134:e62–e82
- **T2DM:** Vienneau et al. 2024, Environ Health 23(1):46
- **Annoyance:** Guski et al. 2017, Int J Environ Res Public Health 14(12):1539
- **Sleep disturbance:** Smith et al. 2022, Environ Health Perspect 130(7):076001
- **Disability weights:** Charalampous et al. 2024, BMJ Public Health 2:e000470

## Citation

Rojas-Rueda D. Health burden of transportation noise in the United States: a national assessment with equity analysis. *Environmental Research.* 2026.

## License

MIT License
