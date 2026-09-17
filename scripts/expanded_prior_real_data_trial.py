"""Exercise corrected UltraNest expanded-prior warm starts on saved real light curves.

This is a manual validation helper. It fits each data set with a deliberately
narrow Rp/R* prior, expands that prior, and compares the corrected warm-started
fit with an independent cold fit over the same expanded prior.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path

import numpy as np

from exotic.api.elca import lc_fitter


DATASETS = (
    {
        "name": "WASP-43b full transit",
        "path": Path(r"D:\WASP43b_codex_current_fix_run\temp\NormalizedFlux_WASP-43 b_2026-03-20.txt"),
        "format": "normalized",
        "prior": {
            "rprs": 0.1594,
            "tmid": 2461120.68976,
            "ars": 4.86,
            "per": 0.813475,
            "inc": 82.6,
            "ecc": 0.0,
            "omega": 90.0,
        },
        "initial_rprs_bounds": [0.155, 0.165],
        "expanded_rprs_bounds": [0.05, 0.30],
        "tmid_bounds": [2461120.686, 2461120.694],
    },
    {
        "name": "TrES-5b full transit",
        "path": Path(
            r"D:\TrES5b_20260716_baron_rp\TrES5b_20260716_baron_rp"
            r"\20260718_112527\Diagnostics\comp7\working_artifacts"
            r"\FinalLightCurve_TrES-5b_2026-07-16.csv"
        ),
        "format": "final_lightcurve",
        "prior": {
            "rprs": 0.143,
            "tmid": 2461238.83817,
            "ars": 6.1,
            "per": 1.48224686,
            "inc": 84.27,
            "ecc": 0.0,
            "omega": 0.0,
        },
        "initial_rprs_bounds": [0.140, 0.146],
        "expanded_rprs_bounds": [0.05, 0.30],
        "tmid_bounds": [2461238.834, 2461238.842],
    },
    {
        "name": "KELT-20b one-sided partial transit",
        "path": Path(
            r"D:\KELT-20\20260718_003832_codex_full_test"
            r"\working_artifacts\NormalizedFlux_KELT-20b_2026-07-15.txt"
        ),
        "format": "normalized",
        "prior": {
            "rprs": 0.1144,
            "tmid": 2461237.776,
            "ars": 7.42,
            "per": 3.4741085,
            "inc": 86.12,
            "ecc": 0.0,
            "omega": 0.0,
        },
        "initial_rprs_bounds": [0.110, 0.120],
        "expanded_rprs_bounds": [0.02, 0.30],
        "tmid_bounds": [2461237.736, 2461237.816],
    },
)


def load_normalized(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    time_values = np.asarray([float(row["BJD"]) for row in rows], dtype=float)
    flux = np.asarray([float(row["Norm Flux"]) for row in rows], dtype=float)
    error = np.asarray([float(row["Norm Err"]) for row in rows], dtype=float)
    airmass = np.asarray([float(row["AM"]) for row in rows], dtype=float)
    return time_values, flux, error, airmass


def load_final_lightcurve(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    values = np.loadtxt(path, delimiter=",", comments="#")
    return values[:, 0], values[:, 2], values[:, 3], values[:, 5]


def load_dataset(config: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if config["format"] == "normalized":
        return load_normalized(config["path"])
    return load_final_lightcurve(config["path"])


def complete_prior(values: dict) -> dict:
    return {
        **values,
        "u0": 0.0,
        "u1": 0.0,
        "u2": 0.0,
        "u3": 0.0,
        "a0": 1.0,
        "a1": 1.0,
        "a2": 0.0,
    }


def summarize_fit(fit: lc_fitter, elapsed_seconds: float) -> dict:
    samples = np.asarray(fit.results["samples"], dtype=float)
    rprs_index = list(fit.sampled_keys).index("rprs")
    rprs_samples = samples[:, rprs_index]
    q16, median, q84 = np.quantile(rprs_samples, [0.16, 0.5, 0.84])
    results = fit.results
    return {
        "elapsed_seconds": elapsed_seconds,
        "ncall": int(results.get("ncall", 0)),
        "posterior_sample_count": int(samples.shape[0]),
        "rprs_q16": float(q16),
        "rprs_median": float(median),
        "rprs_q84": float(q84),
        "rprs_stdev": float(np.std(rprs_samples, ddof=1)),
        "rprs_maximum_likelihood": float(fit.parameters["rprs"]),
        "logz": float(results.get("logz", math.nan)),
        "logzerr": float(results.get("logzerr", math.nan)),
        "warmstart_attempted": bool(
            getattr(fit, "ultranest_expanded_prior_warmstart_attempted", False)
        ),
        "warmstart_applied": bool(
            getattr(fit, "ultranest_expanded_prior_warmstart_applied", False)
        ),
        "warmstart_note": getattr(
            fit, "ultranest_expanded_prior_warmstart_note", None
        ),
        "warmstart_source_sample_count": int(
            getattr(
                fit,
                "ultranest_expanded_prior_warmstart_source_sample_count",
                0,
            )
        ),
        "warmstart_source_effective_sample_size": float(
            getattr(
                fit,
                "ultranest_expanded_prior_warmstart_effective_sample_size",
                0.0,
            )
        ),
        "warmstart_expanded_keys": list(
            getattr(
                fit,
                "ultranest_expanded_prior_warmstart_expanded_keys",
                [],
            )
        ),
        "warmstart_full_prior_fraction": float(
            getattr(
                fit,
                "ultranest_expanded_prior_warmstart_full_prior_fraction",
                math.nan,
            )
        ),
    }


def run_fit(
    *,
    time_values: np.ndarray,
    flux: np.ndarray,
    error: np.ndarray,
    airmass: np.ndarray,
    prior: dict,
    bounds: list[float],
    tmid_bounds: list[float],
    seed: int,
    live_points: int,
    warmstart_source: lc_fitter | None = None,
) -> tuple[lc_fitter, dict]:
    np.random.seed(seed)
    started = time.perf_counter()
    fit = lc_fitter(
        time_values,
        flux,
        error,
        airmass,
        prior,
        {"rprs": list(bounds), "tmid": list(tmid_bounds)},
        mode="ns",
        jd_times=time_values,
        verbose=False,
        use_impactparameter_rather_than_inclination_to_fit=False,
        ultranest_min_num_live_points=live_points,
        ultranest_warmstart_source=warmstart_source,
    )
    elapsed = time.perf_counter() - started
    return fit, summarize_fit(fit, elapsed)


def run_dataset(config: dict, live_points: int, seed: int) -> dict:
    time_values, flux, error, airmass = load_dataset(config)
    finite = (
        np.isfinite(time_values)
        & np.isfinite(flux)
        & np.isfinite(error)
        & np.isfinite(airmass)
        & (error > 0)
    )
    time_values = time_values[finite]
    flux = flux[finite]
    error = error[finite]
    airmass = airmass[finite]
    prior = complete_prior(config["prior"])

    initial_fit, initial_summary = run_fit(
        time_values=time_values,
        flux=flux,
        error=error,
        airmass=airmass,
        prior=prior,
        bounds=config["initial_rprs_bounds"],
        tmid_bounds=config["tmid_bounds"],
        seed=seed,
        live_points=live_points,
    )
    warm_fit, warm_summary = run_fit(
        time_values=time_values,
        flux=flux,
        error=error,
        airmass=airmass,
        prior=prior,
        bounds=config["expanded_rprs_bounds"],
        tmid_bounds=config["tmid_bounds"],
        seed=seed + 1,
        live_points=live_points,
        warmstart_source=initial_fit,
    )
    _, cold_summary = run_fit(
        time_values=time_values,
        flux=flux,
        error=error,
        airmass=airmass,
        prior=prior,
        bounds=config["expanded_rprs_bounds"],
        tmid_bounds=config["tmid_bounds"],
        seed=seed + 2,
        live_points=live_points,
    )

    combined_sigma = math.hypot(
        warm_summary["rprs_stdev"], cold_summary["rprs_stdev"]
    )
    posterior_z = (
        abs(warm_summary["rprs_median"] - cold_summary["rprs_median"])
        / combined_sigma
        if combined_sigma > 0
        else math.inf
    )
    logz_sigma = math.hypot(
        warm_summary["logzerr"], cold_summary["logzerr"]
    )
    logz_z = (
        abs(warm_summary["logz"] - cold_summary["logz"]) / logz_sigma
        if np.isfinite(logz_sigma) and logz_sigma > 0
        else math.nan
    )
    old_lower, old_upper = config["initial_rprs_bounds"]
    warm_outside_old = (
        warm_summary["rprs_median"] < old_lower
        or warm_summary["rprs_median"] > old_upper
    )
    cold_outside_old = (
        cold_summary["rprs_median"] < old_lower
        or cold_summary["rprs_median"] > old_upper
    )
    passed = (
        warm_summary["warmstart_applied"]
        and warm_summary["warmstart_expanded_keys"] == ["rprs"]
        and posterior_z <= 1.0
        and (not np.isfinite(logz_z) or logz_z <= 3.0)
        and warm_outside_old == cold_outside_old
    )
    return {
        "name": config["name"],
        "source_path": str(config["path"]),
        "point_count": int(time_values.size),
        "time_min_bjd_tdb": float(np.min(time_values)),
        "time_max_bjd_tdb": float(np.max(time_values)),
        "initial_rprs_bounds": list(config["initial_rprs_bounds"]),
        "expanded_rprs_bounds": list(config["expanded_rprs_bounds"]),
        "tmid_bounds": list(config["tmid_bounds"]),
        "initial": initial_summary,
        "expanded_warm": warm_summary,
        "expanded_cold": cold_summary,
        "comparison": {
            "posterior_median_difference_sigma": float(posterior_z),
            "logz_difference_sigma": float(logz_z),
            "warm_median_outside_initial_bounds": bool(warm_outside_old),
            "cold_median_outside_initial_bounds": bool(cold_outside_old),
        },
        "passed": bool(passed),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--live-points", type=int, default=100)
    parser.add_argument("--seed", type=int, default=4729)
    parser.add_argument(
        "--dataset-index",
        type=int,
        action="append",
        help="Zero-based dataset index to run; repeat to select more than one.",
    )
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    trial_started = time.perf_counter()
    results = []
    selected_indices = (
        list(range(len(DATASETS)))
        if args.dataset_index is None
        else args.dataset_index
    )
    for index in selected_indices:
        if index < 0 or index >= len(DATASETS):
            parser.error(f"--dataset-index must be between 0 and {len(DATASETS) - 1}")
        config = DATASETS[index]
        print(f"TRIAL START: {config['name']}", flush=True)
        result = run_dataset(
            config,
            live_points=max(40, args.live_points),
            seed=args.seed + index * 100,
        )
        results.append(result)
        print(
            "TRIAL DONE: "
            f"{config['name']} | pass={result['passed']} | "
            f"warm={result['expanded_warm']['warmstart_applied']} | "
            f"posterior_z={result['comparison']['posterior_median_difference_sigma']:.3f} | "
            f"logz_z={result['comparison']['logz_difference_sigma']:.3f}",
            flush=True,
        )

    payload = {
        "live_points": max(40, args.live_points),
        "seed": args.seed,
        "elapsed_seconds": time.perf_counter() - trial_started,
        "all_passed": all(result["passed"] for result in results),
        "datasets": results,
    }
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"RESULTS: {args.output}", flush=True)
    print(f"ALL PASSED: {payload['all_passed']}", flush=True)
    return 0 if payload["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
