"""
Monte Carlo Uncertainty Analysis for US Transportation Noise Health Impact Assessment

Manuscript: "Health burden of transportation noise in the United States: 
a national assessment with equity analysis"
Journal: Environmental Research (ER-26-4325)
Author: David Rojas-Rueda, Colorado State University

This script implements the 12-parameter Monte Carlo simulation described 
in the manuscript. It propagates uncertainty through both YLD (years lived 
with disability) and YLL (years of life lost) components.

Data sources (not included; download separately):
- National Transportation Noise Exposure Map: https://deohs.washington.edu/national-transportation-noise-exposure-map-download
- Global Burden of Disease 2023: https://vizhub.healthdata.org/gbd-results/
- CDC Social Vulnerability Index 2022: https://www.atsdr.cdc.gov/placeandhealth/svi/

Requirements: numpy, scipy
"""

import numpy as np
from scipy import stats
import json
import sys

# ============================================================
# EXPOSURE DATA (from NTNE, +3 dB Lden offset)
# ============================================================

TOTAL_POPULATION = 324_412_244
FRAC_30PLUS = 0.647
POP_30PLUS = TOTAL_POPULATION * FRAC_30PLUS

# Exposure bands: LAeq midpoint, Lden midpoint (+3), population
BANDS = [
    {"laeq_mid": 47.5, "lden_mid": 50.5, "lnight_mid": 41.5, "pop": 43_644_807},
    {"laeq_mid": 55.0, "lden_mid": 58.0, "lnight_mid": 49.0, "pop": 39_031_553},
    {"laeq_mid": 65.0, "lden_mid": 68.0, "lnight_mid": 59.0, "pop":  8_817_826},
    {"laeq_mid": 75.0, "lden_mid": 78.0, "lnight_mid": 69.0, "pop":  2_409_053},
    {"laeq_mid": 85.0, "lden_mid": 88.0, "lnight_mid": 79.0, "pop":    535_416},
    {"laeq_mid": 95.0, "lden_mid": 98.0, "lnight_mid": 89.0, "pop":     22_000},
]

# Baseline disease rates (per 100,000 adults 30+, GBD 2023 national averages)
RATE_IHD = 540
RATE_STROKE = 179
RATE_CM = 25
RATE_T2DM = 735
RATE_CVD_MORT = 480

# Default disability weights and durations
DEFAULTS = {
    "rr_ihd": 1.08,
    "rr_stroke": 1.025,
    "rr_cm": 1.04,
    "rr_t2dm": 1.07,
    "rr_cvd_mort": 1.05,
    "dw_annoyance": 0.011,
    "dw_sleep": 0.010,
    "dw_ihd": 0.08,
    "dur_ihd": 5,
    "dw_stroke": 0.07,
    "dur_stroke": 4,
    "dw_cm": 0.05,
    "dur_cm": 6,
    "dw_t2dm": 0.05,
    "dur_t2dm": 10,
    "remaining_le": 14.5,
    "lden_offset": 3,
    "lnight_offset": -9,
    "road_fraction": 1.0,
    "counterfactual": 53,
}


def pert_params(a, b, c):
    """Convert PERT (min, mode, max) to Beta distribution parameters."""
    if b == a or c == a:
        return a, 0, 1, 1
    lam = 4
    mu = (a + lam * b + c) / (lam + 2)
    if mu == a:
        alpha = 1
    else:
        alpha = ((mu - a) * (2 * b - a - c)) / ((b - mu) * (c - a))
    beta_param = alpha * (c - mu) / (mu - a) if (mu - a) != 0 else 1
    alpha = max(alpha, 0.5)
    beta_param = max(beta_param, 0.5)
    return a, c, alpha, beta_param


def sample_pert(rng, a, b, c, n):
    """Sample from a PERT distribution."""
    lo, hi, alpha, beta_param = pert_params(a, b, c)
    samples = stats.beta.rvs(alpha, beta_param, size=n, random_state=rng)
    return lo + samples * (hi - lo)


def sample_lognormal_rr(rng, rr_central, rr_lo, rr_hi, n):
    """Sample RR from Normal distribution on log scale."""
    log_rr = np.log(rr_central)
    se = (np.log(rr_hi) - np.log(rr_lo)) / (2 * 1.96)
    log_samples = rng.normal(log_rr, se, n)
    return np.exp(log_samples)


def calc_burden(params, bands=BANDS, total_pop=TOTAL_POPULATION, pop_30plus=POP_30PLUS):
    """Calculate total DALYs for a single parameter set."""
    lden_offset = params["lden_offset"]
    lnight_offset = params["lnight_offset"]
    road_frac = params["road_fraction"]
    counterfactual = params["counterfactual"]

    # Annoyance (Guski 2017, capped at 75 dB Lden)
    n_ha = 0
    for b in bands:
        lden = b["laeq_mid"] + lden_offset
        lden_cap = min(lden, 75)
        if lden_cap >= 45:
            pct = max(0, 78.9270 - 3.1162 * lden_cap + 0.0342 * lden_cap ** 2)
            n_ha += b["pop"] * pct / 100

    # Sleep disturbance (Smith 2022, capped at 65 dB Lnight)
    n_hsd = 0
    for b in bands:
        lden = b["laeq_mid"] + lden_offset
        lnight = lden + lnight_offset
        ln_cap = min(lnight, 65)
        if ln_cap >= 40:
            pct = max(0, 31.18323 - 1.47351 * ln_cap + 0.01851 * ln_cap ** 2)
            n_hsd += b["pop"] * pct / 100

    # PAF-based outcomes (no Lden cap for log-linear RRs)
    def calc_paf(rr_per10):
        numerator = 0
        for b in bands:
            lden = b["laeq_mid"] + lden_offset
            if lden > counterfactual:
                delta = lden - counterfactual
                rr_i = rr_per10 ** (delta / 10)
                p_i = b["pop"] / total_pop  # Corrected denominator
                numerator += p_i * (rr_i - 1)
        return numerator / (numerator + 1) if numerator > 0 else 0

    ihd_cases = calc_paf(params["rr_ihd"]) * RATE_IHD * pop_30plus / 100_000
    stroke_cases = calc_paf(params["rr_stroke"]) * RATE_STROKE * pop_30plus / 100_000
    cm_cases = calc_paf(params["rr_cm"]) * RATE_CM * pop_30plus / 100_000
    t2dm_cases = calc_paf(params["rr_t2dm"]) * RATE_T2DM * pop_30plus / 100_000
    cvd_deaths = calc_paf(params["rr_cvd_mort"]) * RATE_CVD_MORT * pop_30plus / 100_000

    # YLD
    yld = (
        ihd_cases * params["dw_ihd"] * params["dur_ihd"]
        + stroke_cases * params["dw_stroke"] * DEFAULTS["dur_stroke"]
        + cm_cases * params["dw_cm"] * DEFAULTS["dur_cm"]
        + t2dm_cases * params["dw_t2dm"] * DEFAULTS["dur_t2dm"]
        + n_ha * params["dw_annoyance"]
        + n_hsd * params["dw_sleep"]
    )

    # YLL
    yll = cvd_deaths * params["remaining_le"]

    # Apply road fraction
    yld *= road_frac
    yll *= road_frac
    total_dalys = yld + yll

    return {
        "total_dalys": total_dalys,
        "yld": yld,
        "yll": yll,
        "ihd_cases": ihd_cases * road_frac,
        "n_ha": n_ha * road_frac,
        "n_hsd": n_hsd * road_frac,
    }


def run_monte_carlo(n_iter=10_000, seed=42):
    """Run the 12-parameter Monte Carlo simulation."""
    rng = np.random.default_rng(seed)

    # Sample all 12 parameters
    rr_ihd = sample_lognormal_rr(rng, 1.08, 1.01, 1.15, n_iter)
    rr_stroke = sample_lognormal_rr(rng, 1.025, 1.009, 1.041, n_iter)
    rr_cm = sample_lognormal_rr(rng, 1.04, 1.02, 1.07, n_iter)
    rr_cvd_mort = sample_lognormal_rr(rng, 1.05, 1.02, 1.07, n_iter)
    lden_offset = sample_pert(rng, 0, 3, 5, n_iter)
    lnight_offset = sample_pert(rng, -11, -9, -7, n_iter)
    dw_ann = sample_pert(rng, 0.005, 0.011, 0.025, n_iter)
    dw_sleep = sample_pert(rng, 0.006, 0.010, 0.07, n_iter)
    dw_ihd = sample_pert(rng, 0.05, 0.08, 0.12, n_iter)
    dur_ihd = sample_pert(rng, 3, 5, 8, n_iter)
    road_frac = rng.uniform(0.6, 1.0, n_iter)
    remaining_le = sample_pert(rng, 12, 14.5, 17, n_iter)

    # Storage
    results = {k: np.zeros(n_iter) for k in ["total_dalys", "yld", "yll", "ihd_cases", "n_ha", "n_hsd"]}

    for i in range(n_iter):
        params = dict(DEFAULTS)
        params.update({
            "rr_ihd": rr_ihd[i],
            "rr_stroke": rr_stroke[i],
            "rr_cm": rr_cm[i],
            "rr_cvd_mort": rr_cvd_mort[i],
            "rr_t2dm": DEFAULTS["rr_t2dm"],  # Not varied
            "lden_offset": lden_offset[i],
            "lnight_offset": lnight_offset[i],
            "dw_annoyance": dw_ann[i],
            "dw_sleep": dw_sleep[i],
            "dw_ihd": dw_ihd[i],
            "dur_ihd": dur_ihd[i],
            "road_fraction": road_frac[i],
            "remaining_le": remaining_le[i],
        })
        out = calc_burden(params)
        for k in results:
            results[k][i] = out[k]

    # Convergence check (CVM)
    running_mean = np.cumsum(results["total_dalys"]) / np.arange(1, n_iter + 1)
    last_1000 = running_mean[-1000:]
    cvm = np.std(last_1000) / np.mean(last_1000) * 100

    # MCSE
    mcse_median = 1.253 * np.std(results["total_dalys"]) / np.sqrt(n_iter)

    # Variance decomposition (squared Spearman correlations)
    param_arrays = {
        "Road fraction": road_frac,
        "DW annoyance": dw_ann,
        "DW sleep": dw_sleep,
        "LAeq-Lden offset": lden_offset,
        "RR CVD mortality": rr_cvd_mort,
        "Remaining LE": remaining_le,
        "Lden-Lnight offset": lnight_offset,
        "RR IHD": rr_ihd,
        "IHD duration": dur_ihd,
        "DW IHD": dw_ihd,
        "RR stroke": rr_stroke,
        "RR cardiomyopathy": rr_cm,
    }

    variance_decomp = {}
    for name, arr in param_arrays.items():
        rho = stats.spearmanr(arr, results["total_dalys"])[0]
        variance_decomp[name] = round(rho ** 2 * 100, 1)

    summary = {
        "n_iterations": n_iter,
        "n_parameters": 12,
        "seed": seed,
        "total_dalys": {
            "median": float(np.median(results["total_dalys"])),
            "p2.5": float(np.percentile(results["total_dalys"], 2.5)),
            "p97.5": float(np.percentile(results["total_dalys"], 97.5)),
            "mean": float(np.mean(results["total_dalys"])),
        },
        "total_yld": {
            "median": float(np.median(results["yld"])),
            "p2.5": float(np.percentile(results["yld"], 2.5)),
            "p97.5": float(np.percentile(results["yld"], 97.5)),
            "mean": float(np.mean(results["yld"])),
        },
        "yll": {
            "median": float(np.median(results["yll"])),
            "p2.5": float(np.percentile(results["yll"], 2.5)),
            "p97.5": float(np.percentile(results["yll"], 97.5)),
            "mean": float(np.mean(results["yll"])),
        },
        "rate_per_100k": {
            "median": float(np.median(results["total_dalys"]) / TOTAL_POPULATION * 100_000),
            "p2.5": float(np.percentile(results["total_dalys"], 2.5) / TOTAL_POPULATION * 100_000),
            "p97.5": float(np.percentile(results["total_dalys"], 97.5) / TOTAL_POPULATION * 100_000),
        },
        "mcse_median_dalys": float(mcse_median),
        "effective_sample_size": n_iter,
        "convergence_cvm_pct": float(cvm),
        "variance_decomposition": variance_decomp,
    }

    return summary


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 42

    print(f"Running Monte Carlo with {n} iterations, seed {seed}...")
    result = run_monte_carlo(n_iter=n, seed=seed)

    print(f"\nTotal DALYs: median {result['total_dalys']['median']:,.0f}")
    print(f"  95% UI: {result['total_dalys']['p2.5']:,.0f} - {result['total_dalys']['p97.5']:,.0f}")
    print(f"  MCSE: {result['mcse_median_dalys']:,.0f}")
    print(f"  CVM: {result['convergence_cvm_pct']:.3f}%")
    print(f"\nRate per 100K: {result['rate_per_100k']['median']:.1f}")
    print(f"  95% UI: {result['rate_per_100k']['p2.5']:.1f} - {result['rate_per_100k']['p97.5']:.1f}")
    print(f"\nVariance decomposition:")
    for k, v in sorted(result["variance_decomposition"].items(), key=lambda x: -x[1]):
        print(f"  {k}: {v}%")

    with open("mc_results.json", "w") as f:
        json.dump(result, f, indent=2)
    print("\nResults saved to mc_results.json")
