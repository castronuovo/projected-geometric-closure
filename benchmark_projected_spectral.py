# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Vitantonio Castronuovo

"""Reproducible analytic and survey-inspired benchmarks for Paper III.

This is not a survey likelihood or forecast.  It tests finite-scale spectral
diagnostics and covariance-profiled identifiability under a fixed,
survey-inspired linear analysis contract.
"""

import json
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.ticker import NullFormatter


ROOT = Path(__file__).resolve().parent
FIG_SPECTRAL = ROOT / "figs" / "fig1_geometric_spectral_test.png"
FIG_IDENTIFIABILITY = ROOT / "figs" / "fig2_identifiability_benchmark.png"
FIG_FINITE_SCALE = ROOT / "figs" / "figS1_finite_scale_spectral_diagnostics.png"
FIG_FLOW = ROOT / "figs" / "fig0_conceptual_flow.png"
FIG_SURVEY = ROOT / "figs" / "fig3_survey_projected_spectral_benchmark.png"
CSV = ROOT / "benchmark_identifiability.csv"
SURVEY_CSV = ROOT / "survey_projected_benchmark.csv"
SURVEY_JSON = ROOT / "survey_projected_benchmark.json"


def orthogonal_residual(signal: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """Return (I-P_basis) signal using a rank-revealing SVD."""
    if basis.size == 0:
        return signal.copy()
    u, singular, _ = np.linalg.svd(basis, full_matrices=False)
    tolerance = np.finfo(float).eps * max(basis.shape) * singular[0]
    rank = int(np.sum(singular > tolerance))
    return signal - u[:, :rank] @ (u[:, :rank].T @ signal)


def polynomial_basis(coordinate: np.ndarray, degree: int) -> np.ndarray:
    """Centered polynomial basis with columns 1,t,...,t**degree."""
    scaled = (coordinate - coordinate.mean()) / coordinate.std()
    return np.column_stack([scaled**order for order in range(degree + 1)])


centers = np.linspace(-4.0, 4.0, 161)
local = np.linspace(-0.6, 0.6, 31)
degrees = (0, 1, 2)
profiled_norm = {degree: [] for degree in degrees}
eta_squared = {degree: [] for degree in degrees}

for center in centers:
    log_x = center + local
    x = np.exp(log_x)
    signal = x**2 / (1.0 + x**2)
    signal_norm = np.linalg.norm(signal)
    for degree in degrees:
        basis = polynomial_basis(log_x, degree)
        residual = orthogonal_residual(signal, basis)
        profiled_norm[degree].append(np.linalg.norm(residual))
        eta_squared[degree].append(
            np.dot(residual, residual) / np.dot(signal, signal)
        )

normalization = max(max(values) for values in profiled_norm.values())
rows = []
for index, center in enumerate(centers):
    row = [center]
    for degree in degrees:
        row.extend(
            [
                profiled_norm[degree][index] / normalization,
                eta_squared[degree][index],
            ]
        )
    rows.append(row)

header = ["ln_window_center"]
for degree in degrees:
    header.extend([f"profiled_norm_p{degree}", f"eta2_p{degree}"])
np.savetxt(CSV, np.asarray(rows), delimiter=",", header=",".join(header), comments="")

# Cross-regime information benchmark at a window centered on the turnover.
x_resolved = np.exp(local)
signal_resolved = x_resolved**2 / (1.0 + x_resolved**2)
background_growth = np.ones_like(signal_resolved)
g_ss = signal_resolved @ signal_resolved
g_sx = signal_resolved @ background_growth
g_xx = background_growth @ background_growth
external_ratio = np.logspace(-3, 3, 250)
information_fraction = 1.0 - (g_sx**2 / g_ss) / (
    g_xx * (1.0 + external_ratio)
)

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

# Positive-spectrum benchmark: one isolated mode versus two positive modes.
u_spectral = np.logspace(-2.5, 2.5, 400)


def spectral_moments(
    u: np.ndarray, weights: np.ndarray, masses_squared: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return A1, A2, A3 for a discrete positive spectral measure."""
    denominator = u[:, None] + masses_squared[None, :]
    return tuple(
        np.sum(weights[None, :] / denominator**order, axis=1)
        for order in (1, 2, 3)
    )


single_moments = spectral_moments(
    u_spectral, np.asarray([1.0]), np.asarray([1.0])
)
mixture_moments = spectral_moments(
    u_spectral, np.asarray([0.65, 0.35]), np.asarray([0.4, 4.0])
)

spectral_figure, spectral_axes = plt.subplots(1, 3, figsize=(7.1, 2.45))
spectral_styles = (
    ("single mode", single_moments, "#1f4e79"),
    ("two-mode mixture", mixture_moments, "#c55a11"),
)
for label, moments, color in spectral_styles:
    a1, a2, a3 = moments
    response = u_spectral * a1
    effective_turnover = a1 / a2 - u_spectral
    spectral_curvature = a1 * a3 / a2**2 - 1.0
    spectral_axes[0].semilogx(
        np.sqrt(u_spectral), response, color=color, linewidth=1.8, label=label
    )
    spectral_axes[1].semilogx(
        np.sqrt(u_spectral), effective_turnover, color=color, linewidth=1.8
    )
    spectral_axes[2].semilogx(
        np.sqrt(u_spectral), spectral_curvature, color=color, linewidth=1.8
    )

spectral_axes[0].set_ylabel(
    r"$\mathcal{R}=-\mathcal{D}_{\rm app}/(\rho_m\Delta_m)$"
)
spectral_axes[0].legend(frameon=False, loc="lower right")
spectral_axes[1].set_ylabel(r"$t_{\rm eff}$")
spectral_axes[2].set_ylabel(r"$A_1A_3/A_2^2-1$")
spectral_axes[2].set_ylim(bottom=-0.002)
for spectral_axis, panel in zip(spectral_axes, ("a", "b", "c")):
    spectral_axis.set_xlabel(r"scale $k$ (reference units)")
    spectral_axis.grid(alpha=0.18, linewidth=0.6)
    spectral_axis.text(
        0.04,
        0.92,
        f"({panel})",
        transform=spectral_axis.transAxes,
        fontweight="bold",
        va="top",
    )

spectral_figure.tight_layout()
spectral_figure.savefig(FIG_SPECTRAL, dpi=300, bbox_inches="tight")

# Finite-scale chord diagnostic used only in the Supporting Information.
u_left, u_right = 0.05, 20.0
u_middle = np.logspace(
    np.log10(u_left) + 0.015, np.log10(u_right) - 0.015, 300
)


def inverse_response(
    u: np.ndarray, weights: np.ndarray, masses_squared: np.ndarray
) -> np.ndarray:
    """Return G=1/F for a positive discrete spectral measure."""
    denominator = u[:, None] + masses_squared[None, :]
    f_value = np.sum(weights[None, :] / denominator, axis=1)
    return 1.0 / f_value


single_weights = np.asarray([1.0])
single_masses = np.asarray([1.0])
mixture_weights = np.asarray([0.65, 0.35])
mixture_masses = np.asarray([0.4, 4.0])

finite_figure, finite_axes = plt.subplots(1, 2, figsize=(7.1, 2.8))
for label, weights, masses, color in (
    ("single mode", single_weights, single_masses, "#1f4e79"),
    ("two-mode mixture", mixture_weights, mixture_masses, "#c55a11"),
):
    u_full = np.logspace(np.log10(u_left), np.log10(u_right), 400)
    g_full = inverse_response(u_full, weights, masses)
    g_end = inverse_response(np.asarray([u_left, u_right]), weights, masses)
    chord_full = (
        (u_right - u_full) * g_end[0]
        + (u_full - u_left) * g_end[1]
    ) / (u_right - u_left)
    finite_axes[0].semilogx(
        u_full, g_full, color=color, linewidth=1.9, label=label
    )
    finite_axes[0].semilogx(
        u_full, chord_full, color=color, linestyle="--", linewidth=1.0
    )

    g_middle = inverse_response(u_middle, weights, masses)
    chord_middle = (
        (u_right - u_middle) * g_end[0]
        + (u_middle - u_left) * g_end[1]
    ) / (u_right - u_left)
    chord_gap = (g_middle - chord_middle) / g_middle
    finite_axes[1].semilogx(
        u_middle, chord_gap, color=color, linewidth=1.9, label=label
    )

finite_axes[0].set_xlabel(r"$u=k^2$ (reference units)")
finite_axes[0].set_ylabel(r"inverse response $G(u)=1/F(u)$")
finite_axes[0].set_title("Response and endpoint chords")
finite_axes[0].legend(frameon=False)
finite_axes[1].set_xlabel(r"middle scale $u_2$")
finite_axes[1].set_ylabel(r"normalized chord gap")
finite_axes[1].set_title("Finite-scale spectral discriminator")
finite_axes[1].axhline(0.0, color="0.35", linestyle=":", linewidth=1.0)
finite_axes[1].legend(frameon=False)
for finite_axis, panel in zip(finite_axes, ("a", "b")):
    finite_axis.grid(alpha=0.18, linewidth=0.6)
    finite_axis.text(
        0.04,
        0.92,
        f"({panel})",
        transform=finite_axis.transAxes,
        fontweight="bold",
        va="top",
    )

finite_figure.tight_layout()
finite_figure.savefig(FIG_FINITE_SCALE, dpi=300, bbox_inches="tight")

figure, axes = plt.subplots(1, 2, figsize=(7.1, 2.9))

colors = ("#1f4e79", "#c55a11", "#548235")
for degree, color in zip(degrees, colors):
    axes[0].plot(
        centers,
        np.asarray(profiled_norm[degree]) / normalization,
        color=color,
        linewidth=1.8,
        label=rf"$p={degree}$",
    )
axes[0].axvline(0.0, color="0.35", linestyle="--", linewidth=1.0)
axes[0].set_xlabel(r"window center $\ln(k/k_{\rm turn})$")
axes[0].set_ylabel("profiled response norm (normalized)")
axes[0].set_title("Resolved-turnover information")
axes[0].set_ylim(bottom=0.0)
axes[0].legend(frameon=False)

axes[1].semilogx(
    external_ratio,
    information_fraction,
    color="#7030a0",
    linewidth=2.0,
)
axes[1].axhline(1.0, color="0.35", linestyle=":", linewidth=1.0)
axes[1].set_xlabel(r"external background information $F_b/F_g$")
axes[1].set_ylabel(r"profiled fraction $F_{A|\alpha}/F_{AA}$")
axes[1].set_title("Cross-regime strengthening")
axes[1].set_ylim(0.0, 1.04)

for axis in axes:
    axis.grid(alpha=0.18, linewidth=0.6)

figure.tight_layout()
figure.savefig(FIG_IDENTIFIABILITY, dpi=300, bbox_inches="tight")

# ---------------------------------------------------------------------------
# Conceptual map: every arrow carries the assumption used at that step.
# ---------------------------------------------------------------------------

flow_figure, flow_axis = plt.subplots(figsize=(7.0, 3.05))
flow_axis.set_xlim(0.0, 1.0)
flow_axis.set_ylim(0.0, 1.0)
flow_axis.axis("off")

flow_boxes = (
    ("5D projected\nresponse", 0.04, 0.64, "#d9eaf7"),
    (r"$\{\mu,\Sigma\}$", 0.375, 0.64, "#e8f1f8"),
    (
        "rank-one apparent sources"
        + "\n"
        + r"$\mathcal{D}_{\rm app}=\mathcal{A}_{\rm app}$",
        0.71,
        0.64,
        "#e2f0d9",
    ),
    ("Stieltjes hierarchy\nand saturation", 0.71, 0.16, "#fff2cc"),
    (
        r"survey map"
        + "\n"
        + r"$\{\mathbf{W},\mathbf{C},\mathbf{X}\}$",
        0.375,
        0.16,
        "#fce4d6",
    ),
    (
        "profiled information\n" + r"$\mathbf{F}_{S|X}$",
        0.04,
        0.16,
        "#eadcf8",
    ),
)
box_width = 0.25
box_height = 0.22
for label, box_x, box_y, color in flow_boxes:
    box = FancyBboxPatch(
        (box_x, box_y),
        box_width,
        box_height,
        boxstyle="round,pad=0.012,rounding_size=0.015",
        linewidth=0.9,
        edgecolor="0.25",
        facecolor=color,
    )
    flow_axis.add_patch(box)
    flow_axis.text(
        box_x + box_width / 2,
        box_y + box_height / 2,
        label,
        ha="center",
        va="center",
        fontsize=7.3,
    )

arrow_labels = (
    "Gauss–Codazzi",
    r"protected $\Sigma=1$",
    r"$d\nu\geq0$",
    "forward model",
    "profiling",
)
arrow_endpoints = (
    ((0.04 + box_width, 0.64 + box_height / 2),
     (0.375, 0.64 + box_height / 2)),
    ((0.375 + box_width, 0.64 + box_height / 2),
     (0.71, 0.64 + box_height / 2)),
    ((0.71 + box_width / 2, 0.64),
     (0.71 + box_width / 2, 0.16 + box_height)),
    ((0.71, 0.16 + box_height / 2),
     (0.375 + box_width, 0.16 + box_height / 2)),
    ((0.375, 0.16 + box_height / 2),
     (0.04 + box_width, 0.16 + box_height / 2)),
)
arrow_label_positions = (
    (0.3325, 0.93, "center"),
    (0.6675, 0.93, "center"),
    (0.855, 0.51, "left"),
    (0.6675, 0.08, "center"),
    (0.3325, 0.08, "center"),
)
for arrow_label, endpoints, label_position in zip(
    arrow_labels, arrow_endpoints, arrow_label_positions
):
    arrow = FancyArrowPatch(
        endpoints[0],
        endpoints[1],
        arrowstyle="-|>",
        mutation_scale=10,
        linewidth=1.0,
        color="0.25",
        shrinkA=3,
        shrinkB=3,
    )
    flow_axis.add_patch(arrow)
    flow_axis.text(
        label_position[0],
        label_position[1],
        arrow_label,
        ha=label_position[2],
        va="center",
        fontsize=7.0,
        color="0.2",
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.4},
    )

flow_figure.tight_layout(pad=0.25)
flow_figure.savefig(FIG_FLOW, dpi=300, bbox_inches="tight")

# ---------------------------------------------------------------------------
# Fixed survey-inspired contract.
#
# Three redshift slices and logarithmic scale bins approximate the resolved
# structure of a full-shape analysis without using a catalogue or likelihood.
# Window mixing, correlated covariance, and all nuisance responses are applied
# before profiling.  Absolute significances are illustrative because the
# covariance normalization is synthetic.
# ---------------------------------------------------------------------------


def block_diagonal(blocks: list[np.ndarray]) -> np.ndarray:
    """Construct a dense block-diagonal matrix."""
    rows = sum(block.shape[0] for block in blocks)
    columns = sum(block.shape[1] for block in blocks)
    output = np.zeros((rows, columns))
    row_offset = 0
    column_offset = 0
    for block in blocks:
        row_count, column_count = block.shape
        output[
            row_offset : row_offset + row_count,
            column_offset : column_offset + column_count,
        ] = block
        row_offset += row_count
        column_offset += column_count
    return output


def window_matrix(log_k: np.ndarray, width: float = 0.105) -> np.ndarray:
    """Return a row-normalized non-diagonal mixing matrix in log k."""
    separation = log_k[:, None] - log_k[None, :]
    core = np.exp(-0.5 * (separation / width) ** 2)
    shoulder = 0.08 * np.exp(
        -0.5 * ((separation - 0.12) / (1.6 * width)) ** 2
    )
    matrix = core + shoulder
    return matrix / matrix.sum(axis=1, keepdims=True)


def protected_response(
    k_values: np.ndarray,
    redshifts: np.ndarray,
    masses: np.ndarray,
    weights: np.ndarray,
    beta: float,
) -> np.ndarray:
    """Return a redshift-major protected response vector."""
    slices = []
    for redshift in redshifts:
        scale_factor = 1.0 / (1.0 + redshift)
        beta_at_epoch = beta * scale_factor**1.2
        denominator = (
            k_values[:, None] ** 2
            + (scale_factor * masses[None, :]) ** 2
        )
        shape = np.sum(
            weights[None, :] * k_values[:, None] ** 2 / denominator,
            axis=1,
        )
        slices.append(beta_at_epoch * shape)
    return np.concatenate(slices)


def whiten(array: np.ndarray, covariance_cholesky: np.ndarray) -> np.ndarray:
    """Apply the fixed lower-Cholesky whitening operator."""
    return np.linalg.solve(covariance_cholesky, array)


def profiled_residual(
    target: np.ndarray,
    design: np.ndarray,
    covariance_cholesky: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Profile a target against design columns in covariance metric."""
    whitened_target = whiten(target, covariance_cholesky)
    whitened_design = whiten(design, covariance_cholesky)
    coefficients, _, _, _ = np.linalg.lstsq(
        whitened_design, whitened_target, rcond=None
    )
    residual = whitened_target - whitened_design @ coefficients
    return residual, coefficients


redshifts = np.asarray([0.50, 0.80, 1.10])
k_bins = np.geomspace(0.02, 0.20, 18)
log_k_bins = np.log(k_bins)
n_redshift = redshifts.size
n_scale = k_bins.size
n_data = n_redshift * n_scale

single_window = window_matrix(log_k_bins)
survey_window = block_diagonal([single_window] * n_redshift)
if not np.allclose(survey_window.sum(axis=1), 1.0, atol=1.0e-12):
    raise RuntimeError("Window rows are not normalized.")

scale_index = np.arange(n_scale)
scale_correlation = np.exp(
    -np.abs(scale_index[:, None] - scale_index[None, :]) / 2.3
)
redshift_index = np.arange(n_redshift)
redshift_correlation = np.exp(
    -np.abs(redshift_index[:, None] - redshift_index[None, :]) / 1.1
)
correlation = (
    0.86 * np.kron(redshift_correlation, scale_correlation)
    + 0.14 * np.eye(n_data)
)
sigma_slices = []
for redshift in redshifts:
    sigma_slices.append(
        (0.0065 + 0.010 * (k_bins / k_bins.max()) ** 1.35)
        * (1.0 + 0.18 * (redshift - redshifts.mean()))
    )
sigma_vector = np.concatenate(sigma_slices)
covariance = sigma_vector[:, None] * correlation * sigma_vector[None, :]
minimum_covariance_eigenvalue = float(np.linalg.eigvalsh(covariance).min())
if minimum_covariance_eigenvalue <= 0.0:
    raise RuntimeError("Survey-inspired covariance is not positive definite.")
covariance_cholesky = np.linalg.cholesky(covariance)

# Five columns spanning three nuisance families: independent slice amplitudes,
# a common smooth tilt, and a common smooth curvature.
scaled_log_k = (log_k_bins - log_k_bins.mean()) / log_k_bins.std()
nuisance_columns = []
for slice_index in range(n_redshift):
    column = np.zeros(n_data)
    begin = slice_index * n_scale
    column[begin : begin + n_scale] = 1.0
    nuisance_columns.append(column)
nuisance_columns.append(np.tile(scaled_log_k, n_redshift))
nuisance_columns.append(np.tile(scaled_log_k**2 - 1.0, n_redshift))
nuisance_raw = np.column_stack(nuisance_columns)
nuisance_observed = np.einsum(
    "ij,jk->ik", survey_window, nuisance_raw, optimize=False
)

single_mass = np.asarray([0.150])
two_masses = np.asarray([0.100, 0.400])
outside_masses = np.asarray([0.100, 0.400])
single_signal = np.einsum(
    "ij,j->i",
    survey_window,
    protected_response(
        k_bins, redshifts, single_mass, np.asarray([1.0]), 0.040
    ),
    optimize=False,
)
two_mode_signal = np.einsum(
    "ij,j->i",
    survey_window,
    protected_response(
        k_bins,
        redshifts,
        two_masses,
        np.asarray([0.50, 0.50]),
        0.040,
    ),
    optimize=False,
)
outside_signal = np.einsum(
    "ij,j->i",
    survey_window,
    protected_response(
        k_bins,
        redshifts,
        outside_masses,
        np.asarray([-0.50, 1.50]),
        0.040,
    ),
    optimize=False,
)

mass_bank = np.unique(
    np.concatenate([np.geomspace(0.10, 0.80, 321), single_mass])
)


def best_single_mode_fit(target: np.ndarray) -> dict[str, float]:
    """Fit a non-negative single-mode template after nuisance profiling."""
    whitened_target = whiten(target, covariance_cholesky)
    target_chi2 = float(np.dot(whitened_target, whitened_target))
    whitened_nuisance = whiten(nuisance_observed, covariance_cholesky)
    nuisance_u, nuisance_singular, _ = np.linalg.svd(
        whitened_nuisance, full_matrices=False
    )
    nuisance_rank_local = int(
        np.sum(
            nuisance_singular
            > np.finfo(float).eps
            * max(whitened_nuisance.shape)
            * nuisance_singular[0]
        )
    )
    nuisance_basis = nuisance_u[:, :nuisance_rank_local]
    hardened_target = whitened_target - nuisance_basis @ (
        nuisance_basis.T @ whitened_target
    )
    nuisance_hardened_chi2 = float(
        np.dot(hardened_target, hardened_target)
    )
    best = {
        "mass": np.nan,
        "amplitude": np.nan,
        "residual_chi2": np.inf,
        "raw_residual_fraction": np.inf,
        "hardened_residual_fraction": np.inf,
    }
    for mass in mass_bank:
        unit_template = np.einsum(
            "ij,j->i",
            survey_window,
            protected_response(
                k_bins,
                redshifts,
                np.asarray([mass]),
                np.asarray([1.0]),
                1.0,
            ),
            optimize=False,
        )
        whitened_template = whiten(unit_template, covariance_cholesky)
        hardened_template = whitened_template - nuisance_basis @ (
            nuisance_basis.T @ whitened_template
        )
        template_norm = float(
            np.dot(hardened_template, hardened_template)
        )
        amplitude = max(
            0.0,
            float(np.dot(hardened_template, hardened_target) / template_norm),
        )
        residual = hardened_target - amplitude * hardened_template
        residual_chi2 = float(np.dot(residual, residual))
        if residual_chi2 < best["residual_chi2"]:
            best = {
                "mass": float(mass),
                "amplitude": amplitude,
                "residual_chi2": residual_chi2,
                "raw_residual_fraction": residual_chi2 / target_chi2,
                "hardened_residual_fraction": (
                    residual_chi2 / nuisance_hardened_chi2
                ),
            }
    return best


def nuisance_fraction(target: np.ndarray) -> float:
    """Return the signal fraction that survives nuisance-only profiling."""
    residual, _ = profiled_residual(
        target, nuisance_observed, covariance_cholesky
    )
    whitened_target = whiten(target, covariance_cholesky)
    return float(
        np.dot(residual, residual)
        / np.dot(whitened_target, whitened_target)
    )


survey_scenarios = {
    "single_mode": single_signal,
    "positive_two_mode": two_mode_signal,
    "signed_out_of_class": outside_signal,
}
scenario_masses = {
    "single_mode": (float(single_mass[0]), np.nan),
    "positive_two_mode": (float(two_masses[0]), float(two_masses[1])),
    "signed_out_of_class": (
        float(outside_masses[0]),
        float(outside_masses[1]),
    ),
}
survey_results = {}
for scenario_name, scenario_signal in survey_scenarios.items():
    fit_result = best_single_mode_fit(scenario_signal)
    fit_result["nuisance_eta2"] = nuisance_fraction(scenario_signal)
    survey_results[scenario_name] = fit_result

whitened_nuisance = whiten(nuisance_observed, covariance_cholesky)
u_nuisance, singular_nuisance, _ = np.linalg.svd(
    whitened_nuisance, full_matrices=False
)
nuisance_rank = int(
    np.sum(
        singular_nuisance
        > np.finfo(float).eps
        * max(whitened_nuisance.shape)
        * singular_nuisance[0]
    )
)
nuisance_projector = (
    u_nuisance[:, :nuisance_rank] @ u_nuisance[:, :nuisance_rank].T
)
response_matrix = np.column_stack([single_signal, two_mode_signal])
whitened_response = whiten(response_matrix, covariance_cholesky)
profiled_information = whitened_response.T @ (
    np.eye(n_data) - nuisance_projector
) @ whitened_response
profiled_eigenvalues = np.linalg.eigvalsh(profiled_information)

metadata = {
    "contract": {
        "redshifts": redshifts.tolist(),
        "k_min_h_Mpc": float(k_bins.min()),
        "k_max_h_Mpc": float(k_bins.max()),
        "bins_per_slice": int(n_scale),
        "window_logk_width": 0.105,
        "nuisance_families": [
            "independent_slice_amplitude",
            "common_logk_tilt",
            "common_logk_curvature",
        ],
        "nuisance_columns": int(nuisance_observed.shape[1]),
        "observational_data_ingested": False,
        "interpretation": "survey-inspired structural benchmark, not forecast",
        "target_amplitude_law": "beta(a)=0.040*a**1.2",
        "targets": {
            "single_mode": {
                "masses_h_Mpc": single_mass.tolist(),
                "weights": [1.0],
            },
            "positive_two_mode": {
                "masses_h_Mpc": two_masses.tolist(),
                "weights": [0.50, 0.50],
            },
            "signed_out_of_class": {
                "masses_h_Mpc": outside_masses.tolist(),
                "weights": [-0.50, 1.50],
            },
        },
        "single_mode_bank": {
            "logarithmic_grid_points": 321,
            "minimum_mass_h_Mpc": 0.100,
            "maximum_mass_h_Mpc": 0.800,
            "exact_single_target_mass_added": True,
            "amplitude_profiled": True,
        },
    },
    "validation": {
        "window_row_sum_max_error": float(
            np.max(np.abs(survey_window.sum(axis=1) - 1.0))
        ),
        "covariance_minimum_eigenvalue": minimum_covariance_eigenvalue,
        "nuisance_rank": nuisance_rank,
        "profiled_information_eigenvalues": profiled_eigenvalues.tolist(),
    },
    "results": survey_results,
}
with SURVEY_JSON.open("w", encoding="utf-8") as output_file:
    json.dump(metadata, output_file, indent=2, sort_keys=True)
    output_file.write("\n")

with SURVEY_CSV.open("w", newline="", encoding="utf-8") as output_file:
    writer = csv.writer(output_file)
    writer.writerow(
        [
            "scenario",
            "target_mass_1",
            "target_mass_2",
            "nuisance_eta2",
            "best_single_mass",
            "best_single_amplitude",
            "single_fit_residual_chi2",
            "single_fit_raw_residual_fraction",
            "single_fit_hardened_residual_fraction",
        ]
    )
    for scenario_name in survey_scenarios:
        result = survey_results[scenario_name]
        writer.writerow(
            [
                scenario_name,
                scenario_masses[scenario_name][0],
                scenario_masses[scenario_name][1],
                result["nuisance_eta2"],
                result["mass"],
                result["amplitude"],
                result["residual_chi2"],
                result["raw_residual_fraction"],
                result["hardened_residual_fraction"],
            ]
        )

survey_figure, survey_axes = plt.subplots(2, 2, figsize=(7.1, 5.25))
image_window = survey_axes[0, 0].imshow(
    single_window,
    origin="lower",
    aspect="auto",
    cmap="Blues",
    extent=[k_bins.min(), k_bins.max(), k_bins.min(), k_bins.max()],
)
survey_axes[0, 0].set_xscale("log")
survey_axes[0, 0].set_yscale("log")
window_ticks = [0.02, 0.05, 0.10, 0.20]
window_tick_labels = ["0.02", "0.05", "0.10", "0.20"]
survey_axes[0, 0].set_xticks(window_ticks, labels=window_tick_labels)
survey_axes[0, 0].set_yticks(window_ticks, labels=window_tick_labels)
survey_axes[0, 0].xaxis.set_minor_formatter(NullFormatter())
survey_axes[0, 0].yaxis.set_minor_formatter(NullFormatter())
survey_axes[0, 0].set_xlabel(r"input $k\,[h\,{\rm Mpc}^{-1}]$")
survey_axes[0, 0].set_ylabel(r"observed $k\,[h\,{\rm Mpc}^{-1}]$")
survey_axes[0, 0].set_title("Non-diagonal window")
survey_figure.colorbar(
    image_window, ax=survey_axes[0, 0], fraction=0.046, pad=0.04
)

image_correlation = survey_axes[0, 1].imshow(
    correlation,
    origin="lower",
    aspect="auto",
    cmap="RdBu_r",
    vmin=-1.0,
    vmax=1.0,
)
survey_axes[0, 1].set_xlabel("data-vector bin")
survey_axes[0, 1].set_ylabel("data-vector bin")
survey_axes[0, 1].set_title("Correlated covariance")
survey_figure.colorbar(
    image_correlation, ax=survey_axes[0, 1], fraction=0.046, pad=0.04
)

middle_begin = n_scale
middle_end = 2 * n_scale
for label, signal, color, style in (
    ("single mode", single_signal, "#1f4e79", "-"),
    ("positive two-mode", two_mode_signal, "#c55a11", "--"),
    ("signed out-of-class", outside_signal, "#a61c3c", ":"),
):
    survey_axes[1, 0].plot(
        k_bins,
        signal[middle_begin:middle_end],
        color=color,
        linestyle=style,
        linewidth=1.9,
        label=label,
    )
survey_axes[1, 0].set_xscale("log")
survey_axes[1, 0].set_xlabel(r"$k\,[h\,{\rm Mpc}^{-1}]$")
survey_axes[1, 0].set_ylabel(r"windowed response at $z=0.8$")
survey_axes[1, 0].set_title("Projected spectral shapes")
survey_axes[1, 0].legend(frameon=False, fontsize=7.2)

scenario_labels = ("single\nmode", "positive\ntwo-mode", "signed\nout-of-class")
residual_fractions = [
    survey_results[name]["hardened_residual_fraction"]
    for name in survey_scenarios
]
display_fractions = [max(value, 1.0e-4) for value in residual_fractions]
bars = survey_axes[1, 1].bar(
    scenario_labels,
    display_fractions,
    color=("#1f4e79", "#c55a11", "#a61c3c"),
    alpha=0.88,
)
survey_axes[1, 1].set_yscale("log")
survey_axes[1, 1].set_ylim(1.0e-4, 1.0)
survey_axes[1, 1].set_ylabel(r"fraction of nuisance-hardened signal")
survey_axes[1, 1].set_title("After window, covariance, and nuisances")
survey_axes[1, 1].grid(axis="y", alpha=0.2)
bar_labels = tuple(
    "numerical zero" if value < 1.0e-8 else f"{value:.3f}"
    for value in residual_fractions
)
for bar, bar_label in zip(bars, bar_labels):
    survey_axes[1, 1].text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() * 1.12,
        bar_label,
        ha="center",
        va="bottom",
        fontsize=7.0,
    )

for panel_axis, panel_label in zip(
    survey_axes.flat, ("a", "b", "c", "d")
):
    panel_axis.text(
        0.03,
        0.95,
        f"({panel_label})",
        transform=panel_axis.transAxes,
        va="top",
        ha="left",
        fontweight="bold",
    )

survey_figure.tight_layout()
survey_figure.savefig(FIG_SURVEY, dpi=300, bbox_inches="tight")
