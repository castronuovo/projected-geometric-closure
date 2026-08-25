# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Vitantonio Castronuovo

"""Translate structural residual fractions into conditional information needs.

This script is not a survey forecast. It reads the frozen projected benchmark
and evaluates the exact scaling

    Delta chi2_miss = f_res|X * Delta chi2_target|X.

For a user-supplied quadratic decision threshold T, the nuisance-hardened
target must therefore satisfy Delta chi2_target|X = T / f_res|X. The null
distribution of T is deliberately not assigned here; it must be calibrated
with the end-to-end mocks of a concrete analysis.
"""

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
INPUT_JSON = ROOT / "survey_projected_benchmark.json"
OUTPUT_CSV = ROOT / "conditional_information_requirement.csv"
OUTPUT_FIGURE = ROOT / "figs" / "figS2_conditional_information_requirement.png"
DEFAULT_THRESHOLDS = (1.0, 4.0, 9.0, 16.0, 25.0)
SCENARIOS = ("positive_two_mode", "signed_out_of_class")


def parse_thresholds(raw: str) -> tuple[float, ...]:
    """Parse a comma-separated list of strictly positive thresholds."""
    thresholds = tuple(float(value.strip()) for value in raw.split(","))
    if not thresholds or any(not np.isfinite(value) or value <= 0 for value in thresholds):
        raise argparse.ArgumentTypeError("thresholds must be finite and positive")
    return thresholds


def load_residual_fractions(path: Path) -> dict[str, float]:
    """Load and validate the structural squared-norm residual fractions."""
    with path.open(encoding="utf-8") as stream:
        payload = json.load(stream)
    fractions = {
        scenario: float(payload["results"][scenario]["hardened_residual_fraction"])
        for scenario in SCENARIOS
    }
    for scenario, fraction in fractions.items():
        if not np.isfinite(fraction) or not 0.0 < fraction <= 1.0:
            raise ValueError(f"invalid residual fraction for {scenario}: {fraction}")
    return fractions


def conditional_requirement(fraction: float, threshold: float) -> tuple[float, float]:
    """Return required target Delta chi2 and its covariance-weighted norm."""
    target_delta_chi2 = threshold / fraction
    return target_delta_chi2, float(np.sqrt(target_delta_chi2))


def write_csv(
    path: Path,
    fractions: dict[str, float],
    thresholds: tuple[float, ...],
) -> None:
    """Write machine-readable conditional requirements."""
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(
            [
                "scenario",
                "hardened_residual_fraction",
                "quadratic_threshold_T",
                "required_target_delta_chi2",
                "required_target_norm",
            ]
        )
        for scenario in SCENARIOS:
            fraction = fractions[scenario]
            for threshold in thresholds:
                target_delta_chi2, target_norm = conditional_requirement(
                    fraction, threshold
                )
                writer.writerow(
                    [scenario, fraction, threshold, target_delta_chi2, target_norm]
                )


def write_figure(
    path: Path,
    fractions: dict[str, float],
    thresholds: tuple[float, ...],
) -> None:
    """Plot the conditional target-information requirement."""
    threshold_grid = np.geomspace(min(thresholds) / 2.0, max(thresholds), 240)
    labels = {
        "positive_two_mode": "positive two-mode",
        "signed_out_of_class": "signed out-of-class",
    }
    colors = {
        "positive_two_mode": "#c55a11",
        "signed_out_of_class": "#a61c3c",
    }
    figure, axis = plt.subplots(figsize=(6.2, 3.8))
    for scenario in SCENARIOS:
        fraction = fractions[scenario]
        required_norm = np.sqrt(threshold_grid / fraction)
        axis.plot(
            threshold_grid,
            required_norm,
            linewidth=2.2,
            color=colors[scenario],
            label=rf"{labels[scenario]} ($f_{{\rm res|X}}={fraction:.3f}$)",
        )
        marker_norm = np.sqrt(np.asarray(thresholds) / fraction)
        axis.scatter(thresholds, marker_norm, s=24, color=colors[scenario], zorder=3)

    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel(r"generic calibrated quadratic threshold $T$")
    axis.set_ylabel(r"required target norm $(T/f_{\rm res|X})^{1/2}$")
    axis.set_title("Conditional information requirement (not a survey forecast)")
    axis.grid(which="both", alpha=0.22)
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(path, dpi=300, bbox_inches="tight")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--thresholds",
        type=parse_thresholds,
        default=DEFAULT_THRESHOLDS,
        help="comma-separated positive quadratic thresholds (default: 1,4,9,16,25)",
    )
    arguments = parser.parse_args()
    fractions = load_residual_fractions(INPUT_JSON)
    write_csv(OUTPUT_CSV, fractions, arguments.thresholds)
    write_figure(OUTPUT_FIGURE, fractions, arguments.thresholds)
    print(f"Wrote {OUTPUT_CSV.name} and {OUTPUT_FIGURE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
