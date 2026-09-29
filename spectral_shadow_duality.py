#!/usr/bin/env python3
"""Controlled illustration of finite-data degeneracy and continuum rigidity.

This is not a proof of the analytic uniqueness theorem.  It constructs two
distinct non-negative atomic measures that agree under a three-dimensional
sampling map and verifies that their continuous Stieltjes responses differ
away from the retained samples.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = Path(__file__).with_name("spectral_shadow_duality.json")
FIGURE_PATHS = (
    ROOT / "figs" / "fig4_spectral_shadow_duality.png",
)


def kernel(u: np.ndarray, support: np.ndarray) -> np.ndarray:
    """Stieltjes sampling matrix K_ij = 1/(u_i+t_j)."""
    return 1.0 / (u[:, None] + support[None, :])


def main() -> None:
    support = np.array([0.10, 0.28, 0.75, 2.0, 5.4, 14.0])
    retained_u = np.array([0.16, 1.05, 7.0])
    finite_map = kernel(retained_u, support)

    # A null direction exists because six positive spectral weights are
    # projected into three retained response values.
    _, _, vh = np.linalg.svd(finite_map, full_matrices=True)
    null_direction = vh[-1]
    null_direction /= np.max(np.abs(null_direction))

    central_weights = np.array([0.80, 0.62, 0.48, 0.36, 0.28, 0.22])
    admissible = central_weights[np.abs(null_direction) > 1.0e-14] / np.abs(
        null_direction[np.abs(null_direction) > 1.0e-14]
    )
    step = 0.78 * np.min(admissible)
    weights_a = central_weights + step * null_direction
    weights_b = central_weights - step * null_direction

    finite_a = np.sum(finite_map * weights_a[None, :], axis=1)
    finite_b = np.sum(finite_map * weights_b[None, :], axis=1)
    finite_residual = np.linalg.norm(finite_a - finite_b)

    dense_u = np.geomspace(0.06, 30.0, 700)
    dense_map = kernel(dense_u, support)
    response_a = np.sum(dense_map * weights_a[None, :], axis=1)
    response_b = np.sum(dense_map * weights_b[None, :], axis=1)
    dense_difference = np.abs(response_a - response_b)

    dense_rank = int(np.linalg.matrix_rank(dense_map, tol=1.0e-12))
    finite_rank = int(np.linalg.matrix_rank(finite_map, tol=1.0e-12))
    max_dense_difference = float(np.max(dense_difference))

    assert np.all(weights_a >= 0.0) and np.all(weights_b >= 0.0)
    assert finite_rank == retained_u.size
    assert finite_residual < 1.0e-12
    assert dense_rank == support.size
    assert max_dense_difference > 1.0e-4

    summary = {
        "support": support.tolist(),
        "retained_u": retained_u.tolist(),
        "weights_a": weights_a.tolist(),
        "weights_b": weights_b.tolist(),
        "finite_map_rank": finite_rank,
        "dense_fixed_support_rank": dense_rank,
        "finite_response_residual_l2": float(finite_residual),
        "maximum_continuum_response_difference": max_dense_difference,
        "interpretation": (
            "Two distinct positive measures are exactly degenerate under the "
            "declared three-sample map but differ as continuous Stieltjes functions."
        ),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8})
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.15), constrained_layout=True)

    x = np.arange(support.size)
    width = 0.38
    axes[0].bar(x - width / 2, weights_a, width, label="measure A", color="#2468a2")
    axes[0].bar(x + width / 2, weights_b, width, label="measure B", color="#d66a2c")
    axes[0].set_xticks(x, [f"{value:g}" for value in support], rotation=35)
    axes[0].set_xlabel(r"spectral support $t_j$")
    axes[0].set_ylabel(r"positive weight $w_j$")
    axes[0].set_title("(a) Distinct positive spectra")
    axes[0].legend(frameon=False)

    axes[1].plot(dense_u, response_a, color="#2468a2", lw=2.0, label=r"$F_A(u)$")
    axes[1].plot(dense_u, response_b, color="#d66a2c", lw=1.8, ls="--", label=r"$F_B(u)$")
    axes[1].scatter(retained_u, finite_a, color="black", marker="o", s=25, zorder=4,
                    label="retained data")
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"resolved scale variable $u$")
    axes[1].set_ylabel(r"Stieltjes response $F(u)$")
    axes[1].set_title("(b) Exact finite-data degeneracy")
    axes[1].legend(frameon=False)

    axes[2].plot(dense_u, dense_difference, color="#6f3c8f", lw=2.0)
    displayed_finite_residual = np.maximum(np.abs(finite_a - finite_b), 2.0e-16)
    axes[2].scatter(retained_u, displayed_finite_residual, color="black", s=25, zorder=4,
                    label="retained equalities")
    axes[2].axhline(1.0e-12, color="0.55", lw=1.0, ls=":", label="numerical tolerance")
    axes[2].set_xscale("log")
    axes[2].set_yscale("log")
    axes[2].set_ylim(1.0e-17, max_dense_difference * 2.0)
    axes[2].set_xlabel(r"resolved scale variable $u$")
    axes[2].set_ylabel(r"$|F_A(u)-F_B(u)|$")
    axes[2].set_title("(c) Continuum response separates")
    axes[2].legend(frameon=False, loc="upper right")

    for axis in axes:
        axis.grid(alpha=0.20, which="both")

    for path in FIGURE_PATHS:
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
