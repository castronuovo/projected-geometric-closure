# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Vitantonio Castronuovo

"""Controlled synthetic injection--recovery calibration.

This executable reuses the frozen response-space contract defined in
``benchmark_projected_spectral.py``.  It does not ingest observational data
and is not a survey forecast.  The mocks inject the declared signal, nuisance
columns, and correlated Gaussian noise before applying the same whitening,
nuisance projection, non-negative amplitude fit, and mass scan used by the
deterministic benchmark.

Two predeclared diagnostics are calibrated on statistically independent mock
sets: (i) lack of fit of the best protected single-mode template and (ii) a
fixed dual-cone witness obtained from the deterministic signed target.  A
window-mismatch stress test changes the injection operator while leaving the
recovery operator fixed.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import benchmark_projected_spectral as benchmark


ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = ROOT / "injection_recovery_config.json"
OUTPUT_CSV = ROOT / "injection_recovery_summary.csv"
OUTPUT_JSON = ROOT / "injection_recovery_summary.json"
OUTPUT_FIGURE = ROOT / "figs" / "figS3_injection_recovery.png"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="JSON analysis contract (default: %(default)s)",
    )
    return parser.parse_args()


def load_config(path: Path) -> dict:
    with path.resolve().open(encoding="utf-8") as stream:
        config = json.load(stream)
    required = {
        "calibration_realizations",
        "decision_quantile",
        "evaluation_realizations",
        "model_mismatch",
        "nuisance_coefficient_sigma",
        "seed",
        "target_norms",
    }
    if set(config) != required:
        raise ValueError("Injection--recovery configuration keys do not match the contract.")
    if int(config["calibration_realizations"]) < 1000:
        raise ValueError("At least 1000 calibration realizations are required.")
    if int(config["evaluation_realizations"]) < 1000:
        raise ValueError("At least 1000 evaluation realizations are required.")
    if not 0.5 < float(config["decision_quantile"]) < 1.0:
        raise ValueError("decision_quantile must lie strictly between 0.5 and 1.")
    target_norms = np.asarray(config["target_norms"], dtype=float)
    if target_norms.ndim != 1 or target_norms.size < 2:
        raise ValueError("target_norms must be a one-dimensional grid.")
    if target_norms[0] != 0.0 or np.any(np.diff(target_norms) <= 0.0):
        raise ValueError("target_norms must start at zero and be strictly increasing.")
    return config


def safe_matmul(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Evaluate a matrix product and reject non-finite numerical output."""
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        result = np.matmul(left, right)
    if not np.all(np.isfinite(result)):
        raise RuntimeError("Non-finite value in injection--recovery matrix product.")
    return result


WHITENED_NUISANCE = benchmark.whiten(
    benchmark.nuisance_observed, benchmark.covariance_cholesky
)
NUISANCE_BASIS = benchmark.nuisance_orthobasis
RECOVERY_PROJECTOR = (
    np.eye(benchmark.n_data) - safe_matmul(NUISANCE_BASIS, NUISANCE_BASIS.T)
)
RECOVERY_TEMPLATES = benchmark.projected_cone_design(benchmark.mass_bank)
TEMPLATE_NORMS_SQUARED = np.sum(RECOVERY_TEMPLATES**2, axis=0)
if np.any(TEMPLATE_NORMS_SQUARED <= 0.0):
    raise RuntimeError("The recovery bank contains a null projected template.")


def harden(data: np.ndarray) -> np.ndarray:
    """Whiten and nuisance-project one vector or a column-major mock matrix."""
    whitened = benchmark.whiten(data, benchmark.covariance_cholesky)
    return safe_matmul(RECOVERY_PROJECTOR, whitened)


def normalize_signal(signal: np.ndarray, target_norm: float) -> np.ndarray:
    """Scale a raw signal to a declared nuisance-hardened target norm."""
    if target_norm == 0.0:
        return np.zeros_like(signal)
    norm = float(np.linalg.norm(harden(signal)))
    if norm <= 0.0:
        raise RuntimeError("Cannot normalize a signal removed by nuisance projection.")
    return signal * (target_norm / norm)


def fit_single_mode(hardened_data: np.ndarray) -> dict[str, np.ndarray]:
    """Profile non-negative amplitude and scan the fixed single-mode bank."""
    correlations = safe_matmul(RECOVERY_TEMPLATES.T, hardened_data)
    positive_correlations = np.maximum(correlations, 0.0)
    improvements = positive_correlations**2 / TEMPLATE_NORMS_SQUARED[:, None]
    best_indices = np.argmax(improvements, axis=0)
    columns = np.arange(hardened_data.shape[1])
    data_norms_squared = np.sum(hardened_data**2, axis=0)
    best_improvements = improvements[best_indices, columns]
    best_amplitudes = (
        positive_correlations[best_indices, columns]
        / TEMPLATE_NORMS_SQUARED[best_indices]
    )
    return {
        "residual_chi2": np.maximum(data_norms_squared - best_improvements, 0.0),
        "mass": benchmark.mass_bank[best_indices],
        "amplitude": best_amplitudes,
        "improvements": improvements,
    }


def draw_hardened_mocks(
    rng: np.random.Generator,
    signal: np.ndarray,
    realization_count: int,
    nuisance_sigma: float,
) -> np.ndarray:
    """Inject signal, declared nuisances, and correlated Gaussian noise."""
    standard_noise = rng.standard_normal((benchmark.n_data, realization_count))
    raw_noise = safe_matmul(benchmark.covariance_cholesky, standard_noise)
    nuisance_coefficients = nuisance_sigma * rng.standard_normal(
        (benchmark.nuisance_observed.shape[1], realization_count)
    )
    raw_mocks = (
        signal[:, None]
        + safe_matmul(benchmark.nuisance_observed, nuisance_coefficients)
        + raw_noise
    )
    return harden(raw_mocks)


def mismatched_window_signal(width: float) -> np.ndarray:
    """Generate the protected single mode with a shifted injection window."""
    injection_window = benchmark.block_diagonal(
        [benchmark.window_matrix(benchmark.log_k_bins, width)]
        * benchmark.n_redshift
    )
    raw_response = benchmark.protected_response(
        benchmark.k_bins,
        benchmark.redshifts,
        benchmark.single_mass,
        np.asarray([1.0]),
        0.040,
    )
    return safe_matmul(injection_window, raw_response)


def fixed_dual_witness() -> tuple[np.ndarray, dict[str, float]]:
    """Return the normalized predeclared witness for the signed target."""
    signed_target = benchmark.hardened_signal(benchmark.outside_signal)
    _, residual, iterations = benchmark.nonnegative_least_squares(
        benchmark.cone_design, signed_target
    )
    witness = -residual
    witness_norm = float(np.linalg.norm(witness))
    if witness_norm <= 0.0:
        raise RuntimeError("The signed target did not produce a dual witness.")
    witness /= witness_norm
    generator_products = safe_matmul(benchmark.cone_design.T, witness)
    if float(np.min(generator_products)) < -1.0e-10:
        raise RuntimeError("The fixed witness is not dual feasible.")
    signed_target_product = float(safe_matmul(witness, signed_target))
    if signed_target_product >= 0.0:
        raise RuntimeError("The fixed witness does not separate the signed target.")
    return witness, {
        "construction_iterations": int(iterations),
        "minimum_generator_product": float(np.min(generator_products)),
        "signed_target_product": signed_target_product,
    }


def quantile(values: np.ndarray, probability: float) -> float:
    return float(np.quantile(values, probability, method="higher"))


def summarize(
    scenario_name: str,
    target_norm: float,
    fitted: dict[str, np.ndarray],
    single_threshold: float,
    witness_scores: np.ndarray,
    witness_threshold: float,
    true_mass: float | None,
) -> dict[str, float | str | None]:
    residuals = fitted["residual_chi2"]
    masses = fitted["mass"]
    amplitudes = fitted["amplitude"]
    row: dict[str, float | str | None] = {
        "scenario": scenario_name,
        "target_norm": target_norm,
        "single_mode_threshold": single_threshold,
        "single_mode_rejection_rate": float(np.mean(residuals > single_threshold)),
        "mean_single_mode_residual_chi2": float(np.mean(residuals)),
        "median_recovered_mass_h_Mpc": float(np.median(masses)),
        "recovered_mass_q16_h_Mpc": float(np.quantile(masses, 0.16)),
        "recovered_mass_q84_h_Mpc": float(np.quantile(masses, 0.84)),
        "median_recovered_amplitude": float(np.median(amplitudes)),
        "dual_witness_threshold": witness_threshold,
        "dual_witness_rejection_rate": float(np.mean(witness_scores > witness_threshold)),
        "mean_dual_witness_score": float(np.mean(witness_scores)),
        "true_mass_h_Mpc": true_mass,
        "nominal_profile_delta_chi2_1_coverage": None,
    }
    if true_mass is not None:
        true_index = int(np.argmin(np.abs(benchmark.mass_bank - true_mass)))
        best_improvement = np.max(fitted["improvements"], axis=0)
        true_delta = best_improvement - fitted["improvements"][true_index]
        row["nominal_profile_delta_chi2_1_coverage"] = float(
            np.mean(true_delta <= 1.0)
        )
    return row


def write_csv(rows: list[dict]) -> None:
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_figure(rows: list[dict], decision_quantile: float) -> None:
    labels = {
        "single_mode": "matched single mode",
        "single_mode_window_mismatch": "window-mismatched single mode",
        "positive_two_mode": "positive two-mode",
        "signed_out_of_class": "signed out-of-class",
    }
    colors = {
        "single_mode": "#1f4e79",
        "single_mode_window_mismatch": "#6a7d8f",
        "positive_two_mode": "#c55a11",
        "signed_out_of_class": "#a61c3c",
    }
    styles = {
        "single_mode": "-",
        "single_mode_window_mismatch": ":",
        "positive_two_mode": "--",
        "signed_out_of_class": "-.",
    }
    figure, axes = plt.subplots(1, 2, figsize=(7.2, 3.15), sharex=True, sharey=True)
    for scenario_name in labels:
        selected = [row for row in rows if row["scenario"] == scenario_name]
        norms = [float(row["target_norm"]) for row in selected]
        axes[0].plot(
            norms,
            [float(row["single_mode_rejection_rate"]) for row in selected],
            marker="o",
            markersize=3.4,
            linewidth=1.5,
            linestyle=styles[scenario_name],
            color=colors[scenario_name],
            label=labels[scenario_name],
        )
        axes[1].plot(
            norms,
            [float(row["dual_witness_rejection_rate"]) for row in selected],
            marker="o",
            markersize=3.4,
            linewidth=1.5,
            linestyle=styles[scenario_name],
            color=colors[scenario_name],
            label=labels[scenario_name],
        )
    nominal_size = 1.0 - decision_quantile
    for axis in axes:
        axis.axhline(nominal_size, color="0.45", linewidth=0.9, linestyle="--")
        axis.set_ylim(0.0, 1.02)
        axis.set_xlabel(r"injected nuisance-hardened target norm $\mathcal{N}$")
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("empirical rejection probability")
    axes[0].set_title("Best-single-mode lack of fit")
    axes[1].set_title("Fixed dual-cone witness")
    axes[0].legend(frameon=False, fontsize=7.0, loc="upper left")
    figure.tight_layout()
    figure.savefig(OUTPUT_FIGURE, dpi=300, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    arguments = parse_arguments()
    config = load_config(arguments.config)
    seed = int(config["seed"])
    calibration_count = int(config["calibration_realizations"])
    evaluation_count = int(config["evaluation_realizations"])
    decision_quantile = float(config["decision_quantile"])
    nuisance_sigma = float(config["nuisance_coefficient_sigma"])
    target_norms = [float(value) for value in config["target_norms"]]
    mismatch_width = float(config["model_mismatch"]["window_injection_logk_width"])

    witness, witness_validation = fixed_dual_witness()
    calibration_rng = np.random.default_rng(seed)
    witness_null_mocks = draw_hardened_mocks(
        calibration_rng,
        np.zeros(benchmark.n_data),
        calibration_count,
        nuisance_sigma,
    )
    witness_threshold = quantile(
        -safe_matmul(witness, witness_null_mocks), decision_quantile
    )

    single_thresholds = {}
    for target_norm in target_norms:
        signal = normalize_signal(benchmark.single_signal, target_norm)
        calibration_mocks = draw_hardened_mocks(
            calibration_rng, signal, calibration_count, nuisance_sigma
        )
        single_thresholds[str(target_norm)] = quantile(
            fit_single_mode(calibration_mocks)["residual_chi2"],
            decision_quantile,
        )

    scenario_signals = {
        "single_mode": benchmark.single_signal,
        "single_mode_window_mismatch": mismatched_window_signal(mismatch_width),
        "positive_two_mode": benchmark.two_mode_signal,
        "signed_out_of_class": benchmark.outside_signal,
    }
    true_masses = {
        "single_mode": float(benchmark.single_mass[0]),
        "single_mode_window_mismatch": float(benchmark.single_mass[0]),
        "positive_two_mode": None,
        "signed_out_of_class": None,
    }
    rows = []
    evaluation_rng = np.random.default_rng(seed + 1)
    for scenario_name, base_signal in scenario_signals.items():
        for target_norm in target_norms:
            signal = normalize_signal(base_signal, target_norm)
            mocks = draw_hardened_mocks(
                evaluation_rng, signal, evaluation_count, nuisance_sigma
            )
            fitted = fit_single_mode(mocks)
            witness_scores = -safe_matmul(witness, mocks)
            rows.append(
                summarize(
                    scenario_name,
                    target_norm,
                    fitted,
                    single_thresholds[str(target_norm)],
                    witness_scores,
                    witness_threshold,
                    true_masses[scenario_name],
                )
            )

    nuisance_test_rng = np.random.default_rng(seed + 2)
    nuisance_test = safe_matmul(
        benchmark.nuisance_observed,
        nuisance_test_rng.standard_normal(
            (benchmark.nuisance_observed.shape[1], 128)
        ),
    )
    nuisance_projection_maximum = float(np.max(np.abs(harden(nuisance_test))))
    if nuisance_projection_maximum > 1.0e-10:
        raise RuntimeError("Declared nuisance injections were not projected out.")

    matched_rows = [row for row in rows if row["scenario"] == "single_mode"]
    nominal_size = 1.0 - decision_quantile
    maximum_size_error = max(
        abs(float(row["single_mode_rejection_rate"]) - nominal_size)
        for row in matched_rows
    )
    null_witness_size = float(
        np.mean(
            -safe_matmul(
                witness,
                draw_hardened_mocks(
                    evaluation_rng,
                    np.zeros(benchmark.n_data),
                    evaluation_count,
                    nuisance_sigma,
                ),
            )
            > witness_threshold
        )
    )
    if maximum_size_error > 0.025:
        raise RuntimeError("Matched single-mode test failed its empirical-size audit.")
    if abs(null_witness_size - nominal_size) > 0.025:
        raise RuntimeError("Dual-witness test failed its empirical-size audit.")

    payload = {
        "contract": {
            **config,
            "observational_data_ingested": False,
            "interpretation": "conditional null-point response-space calibration, not a survey forecast",
            "noise_model": "correlated Gaussian draws from the fixed synthetic covariance",
            "recovery": "fixed nuisance projection, non-negative amplitude, 321-point mass bank",
            "dual_test": "fixed deterministic signed-target witness; no witness scan",
        },
        "calibration": {
            "single_mode_statistic": "q1 absolute best-single-mode lack of fit",
            "single_mode_reference_null_mass_h_Mpc": float(
                benchmark.single_mass[0]
            ),
            "single_mode_lack_of_fit_thresholds": single_thresholds,
            "uniform_composite_H1_calibration": False,
            "q_spec_likelihood_ratio_calibrated": False,
            "dual_witness_threshold": witness_threshold,
            "independent_calibration_and_evaluation_seeds": True,
        },
        "validation": {
            "nuisance_projection_maximum_absolute_residual": nuisance_projection_maximum,
            "maximum_matched_single_mode_size_error": maximum_size_error,
            "evaluation_null_dual_witness_rejection_rate": null_witness_size,
            "dual_witness": witness_validation,
        },
        "results": rows,
    }
    write_csv(rows)
    with OUTPUT_JSON.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
    make_figure(rows, decision_quantile)
    print(
        "Controlled injection--recovery calibration complete: "
        f"{len(rows)} scenario--norm rows, "
        f"{calibration_count} calibration and {evaluation_count} evaluation mocks per row."
    )


if __name__ == "__main__":
    main()
