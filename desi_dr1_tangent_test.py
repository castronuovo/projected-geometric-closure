"""Conservative DESI DR1 tangent test of the protected infrared kernel.

This script uses the official DESI DR1 LRG power-spectrum likelihood products
(data vector, window matrix, and covariance including the published systematic
contributions).  It is deliberately not a replacement for the DESI
full-shape likelihood.  On large scales it fits a linear Kaiser reference
model independently in three LRG redshift bins, profiles a declared local
nuisance basis, and measures the covariance-weighted tangent response to

    mu(k, a) = 1 + beta * k^2 / (k^2 + a^2 m_star^2).

The output is a fixed-background, tangent-level diagnostic.  It must not be
reported as a cosmological posterior or as a test of Sigma=1, because no
lensing data are included.  The official source files are downloaded at run
time and are not redistributed.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from urllib.request import urlretrieve

import camb
import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares
from scipy.stats import chi2


BASE_URL = (
    "https://data.desi.lbl.gov/public/dr1/vac/dr1/"
    "full-shape-bao-clustering/v1.0/data/likelihood"
)
TRACERS = (
    ("LRG", 0.4, 0.6),
    ("LRG", 0.6, 0.8),
    ("LRG", 0.8, 1.1),
)
K_MAX_VALUES = (0.08, 0.09, 0.10)  # h/Mpc; predeclared scale-cut audit
A_MATCH = 0.25
BETA_STEP = 2.0e-3

# Fixed cosmology used by the official DESI likelihood example.
H0 = 68.5
OMEGA_M = 0.304158
OMBH2 = 0.0189629
MNU = 0.06
AS = 1.0e-10 * np.exp(3.0)
NS = 0.965

# Project-recorded hashes of the files served by the official DESI DR1 URL.
# The script refuses silently changed inputs.  These hashes are not presented
# as a separately signed DESI manifest.
EXPECTED_SHA256 = {
    "likelihood_spectrum-poles-rotated_syst-rotation-hod-photo_LRG_GCcomb_z0.4-0.6_thetacut0.05.h5":
        "9ee36802183942595ea54b09b61686dedb337a5655ec1003889843082aaa949b",
    "likelihood_spectrum-poles-rotated_syst-rotation-hod-photo_LRG_GCcomb_z0.6-0.8_thetacut0.05.h5":
        "a72c9f8019ce43e47c812cf645df723c84f960c4c17e40bef2869e4b5d0ab81b",
    "likelihood_spectrum-poles-rotated_syst-rotation-hod-photo_LRG_GCcomb_z0.8-1.1_thetacut0.05.h5":
        "b7b8b17e78a3a666bd34a5406cb311c4279a948826a39cbbe0393c7022e71caf",
}


def filename(tracer: str, zlo: float, zhi: float) -> str:
    return (
        "likelihood_spectrum-poles-rotated_"
        f"syst-rotation-hod-photo_{tracer}_GCcomb_"
        f"z{zlo:.1f}-{zhi:.1f}_thetacut0.05.h5"
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_inputs(cache: Path) -> list[Path]:
    cache.mkdir(parents=True, exist_ok=True)
    paths = []
    for tracer, zlo, zhi in TRACERS:
        name = filename(tracer, zlo, zhi)
        path = cache / name
        if not path.exists():
            urlretrieve(f"{BASE_URL}/{name}", path)
        observed_hash = sha256(path)
        if observed_hash != EXPECTED_SHA256[name]:
            raise RuntimeError(
                f"SHA-256 mismatch for {name}: {observed_hash}; "
                f"expected {EXPECTED_SHA256[name]}"
            )
        paths.append(path)
    return paths


def read_likelihood(path: Path, kmax: float) -> dict[str, np.ndarray | float | str]:
    with h5py.File(path, "r") as handle:
        kobs = np.asarray(handle["observable/spectrum/0/k"])
        data = np.concatenate(
            [
                np.asarray(handle["observable/spectrum/0/value"]),
                np.asarray(handle["observable/spectrum/2/value"]),
            ]
        )
        covariance = np.asarray(handle["covariance/value"])
        window = np.asarray(handle["window/value"])
        kth = np.asarray(handle["window/theory/spectrum/0/k"])
        zeff = float(handle["observable/spectrum"].attrs["zeff"])

    keep_one = kobs <= kmax + 1.0e-12
    keep = np.concatenate([keep_one, keep_one])
    return {
        "path": str(path),
        "sha256": sha256(path),
        "zeff": zeff,
        "kobs": np.concatenate([kobs[keep_one], kobs[keep_one]]),
        "data": data[keep],
        "covariance": covariance[np.ix_(keep, keep)],
        "window": window[keep],
        "kth": kth,
    }


def camb_power(z_values: list[float], kmax: float = 0.5):
    h = H0 / 100.0
    omnuh2 = MNU / 93.14
    omch2 = OMEGA_M * h**2 - OMBH2 - omnuh2
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=OMBH2, omch2=omch2, mnu=MNU)
    pars.InitPower.set_params(As=AS, ns=NS)
    pars.set_matter_power(redshifts=sorted(z_values, reverse=True), kmax=kmax)
    pars.NonLinear = camb.model.NonLinear_none
    return camb.get_matter_power_interpolator(
        pars,
        nonlinear=False,
        hubble_units=True,
        k_hunit=True,
        kmax=kmax,
        zmax=max(z_values) + 0.2,
    )


def growth_solution(k: np.ndarray, z: float, beta: float, mass: float):
    """Return D_beta/D_GR and f_beta on the theory grid."""
    target_a = 1.0 / (1.0 + z)

    def solve(one_k: float, one_beta: float) -> tuple[float, float]:
        def rhs(a, y):
            e2 = OMEGA_M / a**3 + 1.0 - OMEGA_M
            omega_ma = (OMEGA_M / a**3) / e2
            dlnh_da = -1.5 * omega_ma / a
            mu = 1.0 + one_beta * one_k**2 / (
                one_k**2 + (a * mass) ** 2
            )
            return (
                y[1],
                -(3.0 / a + dlnh_da) * y[1]
                + 1.5 * omega_ma * mu * y[0] / a**2,
            )

        sol = solve_ivp(
            rhs,
            (A_MATCH, target_a),
            (A_MATCH, 1.0),
            rtol=2.0e-9,
            atol=2.0e-11,
        )
        if not sol.success:
            raise RuntimeError(sol.message)
        d, dp = sol.y[:, -1]
        return float(d), float(target_a * dp / d)

    gr_d, gr_f = solve(float(k[0]), 0.0)
    if beta == 0.0:
        return np.ones_like(k), np.full_like(k, gr_f), gr_f
    d = np.empty_like(k)
    f = np.empty_like(k)
    for index, one_k in enumerate(k):
        d[index], f[index] = solve(float(one_k), beta)
    return d / gr_d, f, gr_f


def theory_vector(
    kth: np.ndarray,
    z: float,
    pk_interp,
    params: np.ndarray,
    beta: float = 0.0,
    mass: float = 0.03,
) -> np.ndarray:
    log_amp, bias, shot_scaled, ct0, ct2 = params
    amp = np.exp(log_amp)
    ratio_d, growth_f, _ = growth_solution(kth, z, beta, mass)
    pk = np.asarray(pk_interp.P(z, kth)) * ratio_d**2
    x2 = (kth / 0.10) ** 2
    p0 = amp * (bias**2 + 2.0 * bias * growth_f / 3.0 + growth_f**2 / 5.0) * pk
    p2 = amp * (4.0 * bias * growth_f / 3.0 + 4.0 * growth_f**2 / 7.0) * pk
    p4 = amp * (8.0 * growth_f**2 / 35.0) * pk
    p0 = p0 + shot_scaled * 1.0e4 + ct0 * x2 * pk
    p2 = p2 + ct2 * x2 * pk
    return np.concatenate([p0, p2, p4])


def numerical_jacobian(function, point: np.ndarray) -> np.ndarray:
    columns = []
    for index, value in enumerate(point):
        step = 2.0e-4 * max(abs(float(value)), 1.0)
        plus = point.copy()
        minus = point.copy()
        plus[index] += step
        minus[index] -= step
        columns.append((function(plus) - function(minus)) / (2.0 * step))
    return np.column_stack(columns)


def fit_bin(dataset: dict, pk_interp) -> dict:
    data = np.asarray(dataset["data"])
    covariance = np.asarray(dataset["covariance"])
    window = np.asarray(dataset["window"])
    kth = np.asarray(dataset["kth"])
    z = float(dataset["zeff"])
    chol = np.linalg.cholesky(covariance)

    def observed(params, beta=0.0, mass=0.03):
        raw = theory_vector(kth, z, pk_interp, params, beta=beta, mass=mass)
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            projected = window @ raw
        if not np.all(np.isfinite(projected)):
            raise FloatingPointError("non-finite window-projected theory vector")
        return projected

    def residual(params):
        return solve_triangular(chol, data - observed(params), lower=True)

    initial = np.array([0.0, 2.0, 0.2, 0.0, 0.0])
    bounds = (
        np.array([-2.0, 0.2, -5.0, -5.0, -5.0]),
        np.array([2.0, 5.0, 5.0, 5.0, 5.0]),
    )
    fit = least_squares(residual, initial, bounds=bounds, xtol=1e-11, ftol=1e-11)
    best = fit.x
    model = observed(best)
    white_residual = solve_triangular(chol, data - model, lower=True)
    jacobian = numerical_jacobian(lambda p: observed(p), best)
    white_nuisance = solve_triangular(chol, jacobian, lower=True)
    q, _ = np.linalg.qr(white_nuisance, mode="reduced")
    projector = np.eye(data.size) - q @ q.T

    return {
        "best": best,
        "model": model,
        "white_residual": white_residual,
        "projector": projector,
        "chol": chol,
        "observed": observed,
        "chi2": float(white_residual @ white_residual),
        "chi2_p_value": float(chi2.sf(
            white_residual @ white_residual,
            data.size - np.linalg.matrix_rank(white_nuisance),
        )),
        "ndata": int(data.size),
        "nuisance_rank": int(np.linalg.matrix_rank(white_nuisance)),
    }


def tangent_result(fits: list[dict], mass: float) -> dict[str, float]:
    signals = []
    residuals = []
    for fit in fits:
        plus = fit["observed"](fit["best"], beta=BETA_STEP, mass=mass)
        minus = fit["observed"](fit["best"], beta=-BETA_STEP, mass=mass)
        tangent = (plus - minus) / (2.0 * BETA_STEP)
        white_tangent = solve_triangular(fit["chol"], tangent, lower=True)
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            projected_signal = fit["projector"] @ white_tangent
            projected_residual = fit["projector"] @ fit["white_residual"]
        if not (np.all(np.isfinite(projected_signal)) and
                np.all(np.isfinite(projected_residual))):
            raise FloatingPointError("non-finite nuisance-projected tangent")
        signals.append(projected_signal)
        residuals.append(projected_residual)

    signal = np.concatenate(signals)
    residual = np.concatenate(residuals)
    information = float(signal @ signal)
    score = float(signal @ residual)
    beta_hat = score / information
    sigma_beta = information ** -0.5
    return {
        "mass_h_mpc": float(mass),
        "beta_hat": float(beta_hat),
        "sigma_beta": float(sigma_beta),
        "local_score_squared": float((beta_hat / sigma_beta) ** 2),
        "formal_tangent_95_low": float(beta_hat - 1.96 * sigma_beta),
        "formal_tangent_95_high": float(beta_hat + 1.96 * sigma_beta),
        "score_z": float(beta_hat / sigma_beta),
        "expected_positive_95_sensitivity": float(1.645 * sigma_beta),
        "profiled_information": information,
    }


def injection_calibration(result: dict[str, float], seed: int = 14025) -> dict:
    """Calibrate the one-dimensional tangent estimator under Gaussian mocks."""
    rng = np.random.default_rng(seed)
    sigma = result["sigma_beta"]
    injections = (0.0, 0.25, 0.50)
    rows = []
    for injected in injections:
        recovered = injected + sigma * rng.standard_normal(10000)
        rows.append({
            "beta_injected": injected,
            "mean_beta_recovered": float(np.mean(recovered)),
            "std_beta_recovered": float(np.std(recovered, ddof=1)),
            "two_sided_95_coverage": float(np.mean(
                np.abs(recovered - injected) <= 1.96 * sigma
            )),
        })
    return {
        "seed": seed,
        "number_of_mocks_per_injection": 10000,
        "mass_h_mpc": result["mass_h_mpc"],
        "rows": rows,
    }


def main() -> None:
    root = Path(__file__).resolve().parent
    cache = root / ".cache" / "desi_dr1"
    paths = fetch_inputs(cache)
    datasets_by_cut = {
        kmax: [read_likelihood(path, kmax) for path in paths]
        for kmax in K_MAX_VALUES
    }
    nominal_datasets = datasets_by_cut[max(K_MAX_VALUES)]
    pk_interp = camb_power([float(item["zeff"]) for item in nominal_datasets])
    fits_by_cut = {
        kmax: [fit_bin(item, pk_interp) for item in datasets]
        for kmax, datasets in datasets_by_cut.items()
    }
    fits = fits_by_cut[max(K_MAX_VALUES)]

    masses = np.geomspace(0.005, 0.30, 41)
    results = [tangent_result(fits, float(mass)) for mass in masses]
    summary = {
        "analysis": "fixed-background DESI DR1 large-scale tangent diagnostic",
        "observational_data_ingested": True,
        "is_full_likelihood": False,
        "tests_sigma_equals_one": False,
        "nominal_kmax_h_mpc": max(K_MAX_VALUES),
        "scale_cut_audit_h_mpc": list(K_MAX_VALUES),
        "matching_scale_factor": A_MATCH,
        "beta_derivative_step": BETA_STEP,
        "fixed_cosmology": {
            "H0_km_s_Mpc": H0,
            "Omega_m": OMEGA_M,
            "ombh2": OMBH2,
            "sum_mnu_eV": MNU,
            "As": AS,
            "ns": NS,
        },
        "source": {
            "base_url": BASE_URL,
            "license": "DESI DR1 terms; source files are not redistributed",
            "files": [
                {"name": Path(item["path"]).name, "sha256": item["sha256"]}
                for item in nominal_datasets
            ],
        },
        "bin_fits": [
            {
                "zeff": float(data["zeff"]),
                "chi2": fit["chi2"],
                "ndata": fit["ndata"],
                "nuisance_rank": fit["nuisance_rank"],
                "chi2_p_value": fit["chi2_p_value"],
                "parameters_logamp_bias_shot_ct0_ct2": fit["best"].tolist(),
            }
            for data, fit in zip(nominal_datasets, fits)
        ],
        "scale_cut_stability_at_mass_near_0p03": [],
        "mass_scan": results,
        "limitations": [
            "linear Kaiser reference with two smooth counterterms",
            "fixed cosmological background and matching data at a=0.25",
            "tangent approximation around beta=0",
            "no lensing data and therefore no direct closure test",
            "not a DESI collaboration full-shape likelihood or posterior",
            "the fitted score is interpreted only after the predeclared scale-cut audit",
        ],
    }
    representative_index = int(np.argmin(abs(masses - 0.03)))
    for kmax, cut_fits in fits_by_cut.items():
        cut_result = tangent_result(
            cut_fits, float(results[representative_index]["mass_h_mpc"])
        )
        summary["scale_cut_stability_at_mass_near_0p03"].append({
            "kmax_h_mpc": kmax,
            "bin_chi2": [fit["chi2"] for fit in cut_fits],
            "bin_degrees_of_freedom": [
                fit["ndata"] - fit["nuisance_rank"] for fit in cut_fits
            ],
            **cut_result,
        })
    summary["tangent_injection_calibration"] = injection_calibration(
        results[representative_index]
    )

    json_path = root / "desi_dr1_tangent_summary.json"
    csv_path = root / "desi_dr1_tangent_scan.csv"
    figure_path = root / "figs" / "figS7_desi_dr1_tangent_test.png"
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    np.savetxt(
        csv_path,
        np.array(
            [
                [row["mass_h_mpc"], row["beta_hat"], row["sigma_beta"],
                 row["score_z"], row["expected_positive_95_sensitivity"]]
                for row in results
            ]
        ),
        delimiter=",",
        header=("mass_h_mpc,beta_hat,sigma_beta,score_z,"
                "expected_positive_95_sensitivity"),
        comments="",
    )

    mass = np.array([row["mass_h_mpc"] for row in results])
    beta_hat = np.array([row["beta_hat"] for row in results])
    sigma = np.array([row["sigma_beta"] for row in results])
    score_z = np.array([row["score_z"] for row in results])
    sensitivity = np.array([
        row["expected_positive_95_sensitivity"] for row in results
    ])
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.7))
    axes[0].plot(mass, score_z, color="#173f6f", lw=1.8,
                 label=r"tangent score $\hat\beta/\sigma_\beta$")
    axes[0].axhline(0.0, color="black", lw=0.8)
    axes[0].axhline(-1.96, color="#777777", lw=0.8, ls="--")
    axes[0].axhline(1.96, color="#777777", lw=0.8, ls="--")
    axes[0].set_xscale("log")
    axes[0].set_xlabel(r"$m_\ast\,[h\,{\rm Mpc}^{-1}]$")
    axes[0].set_ylabel("local score")
    axes[0].legend(frameon=False)
    axes[1].plot(mass, sensitivity, color="#9a3d28", lw=1.8)
    axes[1].set_xscale("log")
    axes[1].set_xlabel(r"$m_\ast\,[h\,{\rm Mpc}^{-1}]$")
    axes[1].set_ylabel(r"expected one-sided $95\%$ sensitivity")
    fig.suptitle("DESI DR1 LRG large-scale tangent diagnostic")
    fig.tight_layout()
    fig.savefig(figure_path, dpi=220)
    plt.close(fig)

    print(json.dumps({
        "json": str(json_path),
        "csv": str(csv_path),
        "figure": str(figure_path),
        "bin_chi2": [fit["chi2"] for fit in fits],
        "representative_m003": results[representative_index],
        "scale_cut_stability": summary["scale_cut_stability_at_mass_near_0p03"],
    }, indent=2))


if __name__ == "__main__":
    main()
