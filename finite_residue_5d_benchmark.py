# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Vitantonio Castronuovo

"""Numerical audit of the exactly soluble finite-residue 5D scalar benchmark.

This script validates the local interval model used in the manuscript.  It is
an analytic benchmark, not a gravitational-master-sector calculation, a
survey likelihood, or a naturalness analysis.
"""

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "finite_residue_5d_spectrum.csv"
JSON_PATH = ROOT / "finite_residue_5d_benchmark.json"
FIGURE_PATH = ROOT / "figs" / "figS4_finite_residue_5d_benchmark.png"

LENGTH = 1.0
BOUNDARY_KINETIC_LENGTH = 12.0
DELTA_H = 4.0e-3
MODE_COUNT = 512
MOMENTUM_COUNT = 500
BISECTION_STEPS = 120
VALIDATION_TOLERANCE = 5.0e-11


def spectral_equation(mass: float, boundary_coefficient: float) -> float:
    """Return -m cot(mL) + r_b m^2 - h."""
    return (
        -mass / np.tan(mass * LENGTH)
        + BOUNDARY_KINETIC_LENGTH * mass**2
        - boundary_coefficient
    )


def bisect_root(left: float, right: float, boundary_coefficient: float) -> float:
    """Bisection for one simple root with a declared sign change."""
    f_left = spectral_equation(left, boundary_coefficient)
    f_right = spectral_equation(right, boundary_coefficient)
    if not (f_left < 0.0 < f_right):
        raise RuntimeError(
            f"Root is not bracketed: [{left:.8g}, {right:.8g}] "
            f"with values [{f_left:.8g}, {f_right:.8g}]."
        )
    for _ in range(BISECTION_STEPS):
        midpoint = 0.5 * (left + right)
        f_midpoint = spectral_equation(midpoint, boundary_coefficient)
        if f_midpoint > 0.0:
            right = midpoint
        else:
            left = midpoint
    return 0.5 * (left + right)


def positive_mass_roots(boundary_coefficient: float, count: int) -> np.ndarray:
    """Return the light root and the first count-1 heavy roots."""
    branch_offset = 1.0e-10 / LENGTH
    roots = [
        bisect_root(
            branch_offset,
            0.5 * np.pi / LENGTH - branch_offset,
            boundary_coefficient,
        )
    ]
    for branch in range(1, count):
        roots.append(
            bisect_root(
                branch * np.pi / LENGTH + branch_offset,
                (branch + 0.5) * np.pi / LENGTH - branch_offset,
                boundary_coefficient,
            )
        )
    return np.asarray(roots)


def heavy_gap_at_criticality() -> float:
    """Return the first nonzero mass at h_c=-1/L."""
    branch_offset = 1.0e-10 / LENGTH
    return bisect_root(
        np.pi / LENGTH + branch_offset,
        1.5 * np.pi / LENGTH - branch_offset,
        -1.0 / LENGTH,
    )


def boundary_residues(masses: np.ndarray) -> np.ndarray:
    """Return Z_n=|u_n(L)|^2 for modes normalized in the r_b inner product."""
    sine = np.sin(masses * LENGTH)
    bulk_norm = LENGTH / 2.0 - np.sin(2.0 * masses * LENGTH) / (4.0 * masses)
    norm = bulk_norm + BOUNDARY_KINETIC_LENGTH * sine**2
    return sine**2 / norm


def exact_resolvent(momentum_squared: np.ndarray, boundary_coefficient: float) -> np.ndarray:
    root_q = np.sqrt(momentum_squared)
    denominator = (
        root_q / np.tanh(root_q * LENGTH)
        + boundary_coefficient
        + BOUNDARY_KINETIC_LENGTH * momentum_squared
    )
    return 1.0 / denominator


def main() -> None:
    critical_h = -1.0 / LENGTH
    boundary_coefficient = critical_h + DELTA_H
    a_b = BOUNDARY_KINETIC_LENGTH + LENGTH / 3.0

    masses = positive_mass_roots(boundary_coefficient, MODE_COUNT)
    masses_squared = masses**2
    residues = boundary_residues(masses)
    weights = BOUNDARY_KINETIC_LENGTH * residues
    cumulative_weight = np.cumsum(weights)

    light_mass_squared_leading = DELTA_H / a_b
    light_mass_squared_second_order = (
        light_mass_squared_leading
        - LENGTH**3 * DELTA_H**2 / (45.0 * a_b**3)
    )
    critical_light_weight = BOUNDARY_KINETIC_LENGTH / a_b
    critical_heavy_weight = LENGTH / (LENGTH + 3.0 * BOUNDARY_KINETIC_LENGTH)

    critical_gap = heavy_gap_at_criticality()
    rho_b = BOUNDARY_KINETIC_LENGTH / LENGTH
    x_gap = critical_gap * LENGTH
    gap_equation_residual = np.tan(x_gap) - x_gap / (1.0 + rho_b * x_gap**2)

    q_max = 0.2 * critical_gap**2
    momentum_squared = np.logspace(
        np.log10(max(1.0e-3 * masses_squared[0], 1.0e-8)),
        np.log10(q_max),
        MOMENTUM_COUNT,
    )
    resolvent_exact = exact_resolvent(momentum_squared, boundary_coefficient)
    response_exact = (
        BOUNDARY_KINETIC_LENGTH * momentum_squared * resolvent_exact
    )
    spectral_resolvent = np.sum(
        residues[None, :] / (momentum_squared[:, None] + masses_squared[None, :]),
        axis=1,
    )
    light_response = weights[0] * momentum_squared / (
        momentum_squared + masses_squared[0]
    )
    heavy_response = response_exact - light_response
    heavy_weight = 1.0 - weights[0]
    resolved_gap = masses[1]
    heavy_bound = heavy_weight * momentum_squared / (
        momentum_squared + resolved_gap**2
    )

    maximum_relative_resolvent_error = float(
        np.max(np.abs(spectral_resolvent / resolvent_exact - 1.0))
    )
    maximum_heavy_bound_violation = float(
        max(0.0, np.max(heavy_response - heavy_bound))
    )
    minimum_heavy_response = float(np.min(heavy_response))
    truncated_weight = float(cumulative_weight[-1])
    exact_weight_sum = 1.0
    weight_shortfall = exact_weight_sum - truncated_weight

    validations = {
        "all_residues_positive": bool(np.all(residues > 0.0)),
        "partial_weight_below_exact_sum": bool(
            truncated_weight <= exact_weight_sum + VALIDATION_TOLERANCE
        ),
        "spectral_resolvent_relative_error": maximum_relative_resolvent_error,
        "heavy_response_minimum": minimum_heavy_response,
        "heavy_bound_maximum_violation": maximum_heavy_bound_violation,
        "critical_gap_equation_residual": float(abs(gap_equation_residual)),
        "critical_gap_in_declared_interval": bool(
            np.pi < x_gap <= 4.493409
        ),
    }

    if not validations["all_residues_positive"]:
        raise RuntimeError("A non-positive boundary residue was found.")
    if not validations["partial_weight_below_exact_sum"]:
        raise RuntimeError("The truncated residue sum exceeds the exact sum rule.")
    if minimum_heavy_response < -VALIDATION_TOLERANCE:
        raise RuntimeError("The reconstructed heavy response is negative.")
    if maximum_heavy_bound_violation > VALIDATION_TOLERANCE:
        raise RuntimeError("The heavy-sector remainder bound was violated.")
    if maximum_relative_resolvent_error > VALIDATION_TOLERANCE:
        raise RuntimeError("The truncated spectral resolvent has not converged.")
    if abs(gap_equation_residual) > VALIDATION_TOLERANCE:
        raise RuntimeError("The critical heavy root fails the analytic gap equation.")
    if not validations["critical_gap_in_declared_interval"]:
        raise RuntimeError("The critical heavy root lies outside the declared interval.")

    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "mode",
                "mass_squared",
                "boundary_residue",
                "normalized_weight",
                "cumulative_normalized_weight",
            ]
        )
        for index in range(MODE_COUNT):
            writer.writerow(
                [
                    index,
                    f"{masses_squared[index]:.17g}",
                    f"{residues[index]:.17g}",
                    f"{weights[index]:.17g}",
                    f"{cumulative_weight[index]:.17g}",
                ]
            )

    output = {
        "scope": (
            "Exactly soluble local 5D scalar interval benchmark; not the full "
            "coupled gravitational master sector and not a survey calculation."
        ),
        "contract": {
            "L": LENGTH,
            "r_b": BOUNDARY_KINETIC_LENGTH,
            "h_c": critical_h,
            "delta_h": DELTA_H,
            "h": boundary_coefficient,
            "mode_count": MODE_COUNT,
            "momentum_count": MOMENTUM_COUNT,
            "q_max_over_critical_gap_squared": 0.2,
        },
        "light_pole": {
            "mass_squared_exact": float(masses_squared[0]),
            "mass_squared_leading": light_mass_squared_leading,
            "mass_squared_second_order": light_mass_squared_second_order,
            "leading_relative_error": float(
                abs(light_mass_squared_leading / masses_squared[0] - 1.0)
            ),
            "second_order_relative_error": float(
                abs(light_mass_squared_second_order / masses_squared[0] - 1.0)
            ),
            "weight_exact": float(weights[0]),
            "weight_critical_limit": critical_light_weight,
        },
        "heavy_sector": {
            "first_mass_squared_at_h": float(masses_squared[1]),
            "critical_gap_squared": float(critical_gap**2),
            "critical_gap_x": float(x_gap),
            "heavy_weight_exact": float(heavy_weight),
            "heavy_weight_critical_limit": critical_heavy_weight,
        },
        "sum_rule": {
            "exact_normalized_weight": exact_weight_sum,
            "truncated_normalized_weight": truncated_weight,
            "truncation_shortfall": float(weight_shortfall),
            "exact_boundary_residue_sum": 1.0 / BOUNDARY_KINETIC_LENGTH,
            "truncated_boundary_residue_sum": float(np.sum(residues)),
        },
        "validation": validations,
    }
    with JSON_PATH.open("w", encoding="utf-8") as stream:
        json.dump(output, stream, indent=2, sort_keys=True)
        stream.write("\n")

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
        }
    )
    figure, axes = plt.subplots(1, 3, figsize=(9.4, 3.0))

    axes[0].semilogx(momentum_squared, response_exact, lw=2.0, label="exact response")
    axes[0].semilogx(momentum_squared, light_response, "--", lw=1.8, label="light pole")
    axes[0].set_xlabel(r"physical momentum squared $q L^2$")
    axes[0].set_ylabel(r"normalized response $S_\Xi$")
    axes[0].set_title("Light-pole approximation")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.2)

    axes[1].loglog(momentum_squared, heavy_response, lw=2.0, label="exact heavy remainder")
    axes[1].loglog(momentum_squared, heavy_bound, ":", lw=2.0, label="gap bound")
    axes[1].set_xlabel(r"physical momentum squared $q L^2$")
    axes[1].set_ylabel(r"heavy contribution $R_H$")
    axes[1].set_title("Controlled heavy sector")
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=0.2)

    mode_index = np.arange(1, MODE_COUNT + 1)
    axes[2].semilogx(mode_index, cumulative_weight, lw=2.0)
    axes[2].axhline(1.0, color="black", ls=":", lw=1.3, label=r"exact $\sum_n w_n=1$")
    axes[2].set_xlabel("number of retained modes")
    axes[2].set_ylabel(r"cumulative normalized residue")
    axes[2].set_title("Finite-residue sum rule")
    axes[2].set_ylim(max(0.0, weights[0] - 0.01), 1.002)
    axes[2].legend(frameon=False)
    axes[2].grid(alpha=0.2)

    figure.suptitle(
        r"Local 5D scalar benchmark: $r_b/L=12$, "
        r"$\delta h\,L=4\times10^{-3}$",
        y=1.01,
    )
    figure.tight_layout()
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(FIGURE_PATH, dpi=300, bbox_inches="tight")
    plt.close(figure)

    print(
        "Validated the finite-residue 5D scalar benchmark: "
        f"m_*^2={masses_squared[0]:.8g}, w_*={weights[0]:.8g}, "
        f"eta^2={critical_gap**2:.8g}, "
        f"partial sum={truncated_weight:.10f}."
    )


if __name__ == "__main__":
    main()
