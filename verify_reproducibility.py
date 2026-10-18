# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Vitantonio Castronuovo

"""Compare regenerated scientific outputs with frozen reference files."""

import csv
from functools import lru_cache
import json
import math
import sys
from pathlib import Path


RELATIVE_TOLERANCE = 1.0e-11
ABSOLUTE_TOLERANCE = 1.0e-14
OPTIMIZATION_RELATIVE_TOLERANCE = 1.0e-7
OPTIMIZATION_ABSOLUTE_TOLERANCE = 1.0e-12
# The DESI fit optimizes nearly degenerate nuisance directions. Independent
# runs retain the scores but need a separate, explicit numerical tolerance.
DESI_RELATIVE_TOLERANCE = 3.0e-5
DESI_ABSOLUTE_TOLERANCE = 1.0e-6
# Fixed-witness reconstruction and attained maximum residual retain the
# independent-integrator absolute tolerance. Comparisons between separately
# optimized witnesses allow 1e-4 of their approximation error in individual
# endpoints/states, without relaxing the maximum-residual comparison.
JOINT_HISTORY_ABSOLUTE_TOLERANCE = 5.0e-11
JOINT_HISTORY_FIT_ERROR_FRACTION = 1.0e-4
NUMERICAL_ZERO_TOLERANCE = 1.0e-10
NUMERICAL_ZERO_DIAGNOSTICS = (
    "minimum_generator_product",
    "nuisance_projection_maximum_absolute_residual",
    "generator_step_halving_error",
    "nesting_max_absolute_error",
    "rk4_step_halving_error",
)
CSV_FILES = (
    "benchmark_identifiability.csv",
    "cone_grid_convergence.csv",
    "finite_residue_5d_spectrum.csv",
    "survey_projected_benchmark.csv",
    "conditional_information_requirement.csv",
    "injection_recovery_summary.csv",
    "cross_epoch_transport_benchmark.csv",
    "causal_growth_transport_benchmark.csv",
    "desi_dr1_tangent_scan.csv",
    "causal_second_order_certified_range.csv",
)
JSON_FILES = (
    "finite_residue_5d_benchmark.json",
    "survey_projected_benchmark.json",
    "injection_recovery_summary.json",
    "causal_growth_transport_summary.json",
    "causal_transport_support_results.json",
    "causal_transport_free_amplitude_results.json",
    "spectral_shadow_duality.json",
    "desi_dr1_tangent_summary.json",
    "causal_second_order_allowance_results.json",
    "causal_amplitude_reachability_results.json",
    "gapped_residue_inverse_diagnostic_results.json",
    "bounded_joint_history_results.json",
)
FILES = CSV_FILES + JSON_FILES
V15_FILES_NOT_IN_LEGACY_CI = frozenset((
    "desi_dr1_tangent_scan.csv",
    "spectral_shadow_duality.json",
    "desi_dr1_tangent_summary.json",
))
V16_FILES_NOT_IN_V15_CI = frozenset((
    "causal_second_order_certified_range.csv",
    "causal_second_order_allowance_results.json",
))
V17_FILES_NOT_IN_V16_CI = frozenset(("causal_amplitude_reachability_results.json",
                                                 "gapped_residue_inverse_diagnostic_results.json",
                                                 "bounded_joint_history_results.json"))
PLATFORM_SENSITIVE_JSON_FILES = (
    "causal_transport_support_results.json",
    "causal_transport_free_amplitude_results.json",
    "causal_second_order_allowance_results.json",
)


def close(reference: float, regenerated: float, *, optimization: bool = False,
          desi: bool = False, absolute_tolerance=None) -> bool:
    if math.isnan(reference) or math.isnan(regenerated):
        return math.isnan(reference) and math.isnan(regenerated)
    return math.isclose(
        reference,
        regenerated,
        rel_tol=(DESI_RELATIVE_TOLERANCE if desi else
                 OPTIMIZATION_RELATIVE_TOLERANCE if optimization else
                 RELATIVE_TOLERANCE),
        abs_tol=(absolute_tolerance if absolute_tolerance is not None else
                 DESI_ABSOLUTE_TOLERANCE if desi else
                 OPTIMIZATION_ABSOLUTE_TOLERANCE if optimization else
                 ABSOLUTE_TOLERANCE),
    )


def compare_csv(reference_path: Path, regenerated_path: Path, *, optimization=False,
                desi=False) -> None:
    with reference_path.open(newline="", encoding="utf-8") as stream:
        reference_rows = list(csv.reader(stream))
    with regenerated_path.open(newline="", encoding="utf-8") as stream:
        regenerated_rows = list(csv.reader(stream))

    if not reference_rows or reference_rows[0] != regenerated_rows[0]:
        raise AssertionError(f"CSV header mismatch: {regenerated_path.name}")
    if len(reference_rows) != len(regenerated_rows):
        raise AssertionError(f"CSV row-count mismatch: {regenerated_path.name}")

    for row_index, (reference_row, regenerated_row) in enumerate(
        zip(reference_rows[1:], regenerated_rows[1:]), start=2
    ):
        if len(reference_row) != len(regenerated_row):
            raise AssertionError(
                f"CSV column-count mismatch: {regenerated_path.name}:{row_index}"
            )
        for column_index, (reference_value, regenerated_value) in enumerate(
            zip(reference_row, regenerated_row), start=1
        ):
            try:
                reference_number = float(reference_value)
                regenerated_number = float(regenerated_value)
            except ValueError:
                if reference_value != regenerated_value:
                    raise AssertionError(
                        f"CSV value mismatch: {regenerated_path.name}:"
                        f"{row_index}:{column_index}"
                    )
            else:
                diagnostic = reference_rows[0][column_index-1]
                if diagnostic in NUMERICAL_ZERO_DIAGNOSTICS and (
                    abs(reference_number) <= NUMERICAL_ZERO_TOLERANCE
                    and abs(regenerated_number) <= NUMERICAL_ZERO_TOLERANCE
                ):
                    continue
                if not close(reference_number, regenerated_number,
                             optimization=optimization, desi=desi):
                    raise AssertionError(
                        f"CSV numerical mismatch: {regenerated_path.name}:"
                        f"{row_index}:{column_index}"
                    )


def compare_json(reference, regenerated, path: str = "root", *, optimization=False,
                 desi=False, absolute_tolerance=None) -> None:
    if isinstance(reference, bool) or reference is None or isinstance(reference, str):
        if reference != regenerated:
            raise AssertionError(f"JSON value mismatch at {path}")
        return
    if isinstance(reference, int):
        if not isinstance(regenerated, int) or isinstance(regenerated, bool):
            raise AssertionError(f"JSON integer mismatch at {path}")
        if reference != regenerated:
            raise AssertionError(f"JSON integer mismatch at {path}")
        return
    if isinstance(reference, float):
        if not isinstance(regenerated, (int, float)) or isinstance(regenerated, bool):
            raise AssertionError(f"JSON type mismatch at {path}")
        if path.endswith(NUMERICAL_ZERO_DIAGNOSTICS):
            if (
                abs(float(reference)) <= NUMERICAL_ZERO_TOLERANCE
                and abs(float(regenerated)) <= NUMERICAL_ZERO_TOLERANCE
            ):
                return
        if not close(float(reference), float(regenerated), optimization=optimization,
                     desi=desi, absolute_tolerance=absolute_tolerance):
            raise AssertionError(f"JSON numerical mismatch at {path}: "
                                 f"{reference!r} != {regenerated!r}")
        return
    if isinstance(reference, list):
        if not isinstance(regenerated, list) or len(reference) != len(regenerated):
            raise AssertionError(f"JSON list mismatch at {path}")
        for index, (reference_item, regenerated_item) in enumerate(
            zip(reference, regenerated)
        ):
            compare_json(reference_item, regenerated_item, f"{path}[{index}]",
                         optimization=optimization, desi=desi,
                         absolute_tolerance=absolute_tolerance)
        return
    if isinstance(reference, dict):
        if not isinstance(regenerated, dict) or reference.keys() != regenerated.keys():
            raise AssertionError(f"JSON object mismatch at {path}")
        for key in reference:
            compare_json(reference[key], regenerated[key], f"{path}.{key}",
                         optimization=optimization, desi=desi,
                         absolute_tolerance=absolute_tolerance)
        return
    raise TypeError(f"Unsupported JSON value at {path}")


def compare_causal_amplitude(reference, regenerated):
    """Compare physical outputs; bound method diagnostics, report provenance."""
    reference, regenerated = dict(reference), dict(regenerated)
    if reference.keys() != regenerated.keys():
        raise AssertionError("Causal-amplitude output schema changed")
    for key in ("numpy", "scipy"):
        old, new = reference.pop(key), regenerated.pop(key)
        if old != new:
            print(f"Causal-amplitude runtime provenance: {key} {old} -> {new}")
    diagnostic_limits = {
        "relative_tangent_finite_difference_error": 1e-6,
        "max_fitted_power_residual": 2e-12,
        "independent_Radau_max_power_residual": 2e-11,
        "DOP853_vs_Radau_fitted_difference": 2e-11,
        "DOP853_vs_Radau_target_difference": 2e-11,
    }
    for key, limit in diagnostic_limits.items():
        for output in (reference, regenerated):
            value = output.pop(key)
            if not isinstance(value, (int, float)) or not 0 <= value < limit:
                raise AssertionError(f"Causal-amplitude diagnostic failed: {key}")
    # Iteration residuals can fluctuate at roundoff, but failed fits cannot.
    for output in (reference, regenerated):
        iterations = output.pop("iterations")
        if not iterations or not 0 <= iterations[-1]["max_power_residual"] < 2e-12:
            raise AssertionError("Causal-amplitude fit did not converge")
    compare_json(reference, regenerated, optimization=True)


def compare_inverse_diagnostic(reference, regenerated):
    """Compare inverse bounds; report hashes of regenerated forward inputs."""
    reference, regenerated = dict(reference), dict(regenerated)
    old, new = reference.pop("input_sha256"), regenerated.pop("input_sha256")
    if old.keys() != new.keys():
        raise AssertionError("Inverse-diagnostic input schema changed")
    if old != new:
        print("Inverse-diagnostic input hashes changed after forward regeneration")
    compare_json(reference, regenerated)


@lru_cache(maxsize=64)
def _joint_history_prediction(coefficients):
    # Reuse the declared forward model to check the recorded witness. This
    # replays a fixed history; it does not run or alter the optimizer.
    from bounded_joint_history import integrate
    prediction, _, states = integrate(coefficients)
    return prediction.tolist(), states


def _check_joint_history_fit(case, name):
    from bounded_joint_history import restrictions
    import numpy as np

    fit = case[name]
    coefficients = fit["coefficients"]
    if (len(coefficients) != 8 or any(isinstance(x, bool) or
            not isinstance(x, (int, float)) or not math.isfinite(x)
            for x in coefficients)):
        raise AssertionError("Bounded joint-history coefficients invalid")
    if fit["parameter_count"] != (2 if name == "H2" else 8):
        raise AssertionError("Bounded joint-history parameter count changed")
    if name == "H2" and any(x != block[0] for block in
            (coefficients[:4], coefficients[4:]) for x in block):
        raise AssertionError("Bounded joint-history tied coefficients changed")
    expected = restrictions(np.asarray(coefficients))
    for key, value in expected.items():
        compare_json(value, fit[key], f"{name}.{key}")
    if not expected["admissible"] or not fit["optimizer_success"]:
        raise AssertionError("Bounded joint-history admissibility/optimizer failed")
    lo, hi = fit["amplitude_bounds"]
    if not (0.015 <= lo <= hi <= 0.025 and
            0 <= fit["conservative_log_slope_bound"] <= 2+1e-12):
        raise AssertionError("Bounded joint-history prior violated")

    prediction, states = _joint_history_prediction(tuple(coefficients))
    compare_json(prediction, fit["prediction"], f"{name}.forward_prediction",
                 absolute_tolerance=JOINT_HISTORY_ABSOLUTE_TOLERANCE)
    compare_json(states, fit["states"], f"{name}.forward_states",
                 absolute_tolerance=JOINT_HISTORY_ABSOLUTE_TOLERANCE)
    residual = [p-t for p, t in zip(fit["prediction"], case["target"])]
    if len(residual) != 8:
        raise AssertionError("Bounded joint-history snapshot count changed")
    for key, indices in (("max_joint_residual", range(8)),
                         ("max_power_residual", (0, 1, 4, 5)),
                         ("max_growthrate_residual", (2, 3, 6, 7))):
        compare_json(max(abs(residual[i]) for i in indices), fit[key],
                     f"{name}.{key}")
    if name == "H8" and abs(fit["max_joint_residual"] -
            fit["independent_Radau_max_joint_residual"]) > JOINT_HISTORY_ABSOLUTE_TOLERANCE:
        raise AssertionError("Bounded joint-history independent residual failed")


def _joint_fit_comparison_tolerance(reference, regenerated):
    return max(JOINT_HISTORY_ABSOLUTE_TOLERANCE,
               JOINT_HISTORY_FIT_ERROR_FRACTION * max(
                   reference["max_joint_residual"], regenerated["max_joint_residual"]))


def compare_bounded_joint_history(reference, regenerated):
    """Check each history witness, then compare its physical outputs.

    Optimizer coefficients and their amplitude/slope envelopes are solution
    dependent. They must be internally consistent and admissible, rather than
    identical between platforms. Maximum residuals retain a 5e-11 absolute
    comparison. Individual endpoints/states may differ by at most 1e-4 of the
    attained approximation error (with a 5e-11 numerical floor). The contract
    and structural metadata retain strict checks.
    """
    reference = json.loads(json.dumps(reference))
    regenerated = json.loads(json.dumps(regenerated))
    if reference.keys() != regenerated.keys():
        raise AssertionError("Bounded joint-history output schema changed")
    for key in ("numpy", "scipy"):
        old, new = reference.pop(key), regenerated.pop(key)
        if old != new:
            print(f"Bounded joint-history runtime provenance: {key} {old} -> {new}")
    compare_json(reference.pop("contract"), regenerated.pop("contract"), "contract")
    reference_cases = reference.pop("cases")
    regenerated_cases = regenerated.pop("cases")
    if len(reference_cases) != len(regenerated_cases):
        raise AssertionError("Bounded joint-history case count changed")
    compare_json(reference, regenerated, optimization=True)
    for old_case, new_case in zip(reference_cases, regenerated_cases):
        for output in (old_case, new_case):
            for name in ("H2", "H8"):
                _check_joint_history_fit(output, name)
        for name in ("H2", "H8"):
            old_fit, new_fit = old_case.pop(name), new_case.pop(name)
            if old_fit.keys() != new_fit.keys():
                raise AssertionError("Bounded joint-history fit schema changed")
            for fit in (old_fit, new_fit):
                for key in ("coefficients", "evaluations", "amplitude_bounds",
                            "conservative_log_slope_bound"):
                    fit.pop(key)
            tolerance = _joint_fit_comparison_tolerance(old_fit, new_fit)
            for key in ("max_joint_residual", "independent_Radau_max_joint_residual"):
                if key in old_fit:
                    compare_json(old_fit.pop(key), new_fit.pop(key), f"{name}.{key}",
                                 absolute_tolerance=JOINT_HISTORY_ABSOLUTE_TOLERANCE)
            compare_json(old_fit, new_fit, name, absolute_tolerance=tolerance)
        compare_json(old_case, new_case)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: verify_reproducibility.py REFERENCE_DIRECTORY")
    reference_directory = Path(sys.argv[1]).resolve()
    regenerated_directory = Path(__file__).resolve().parent
    missing = {name for name in FILES if not (reference_directory / name).is_file()}
    if missing and missing not in (
        V17_FILES_NOT_IN_V16_CI,
        V17_FILES_NOT_IN_V16_CI | V16_FILES_NOT_IN_V15_CI,
        V17_FILES_NOT_IN_V16_CI | V16_FILES_NOT_IN_V15_CI | V15_FILES_NOT_IN_LEGACY_CI,
    ):
        raise FileNotFoundError(
            "Incomplete frozen-output set: " + ", ".join(sorted(missing))
        )

    for filename in CSV_FILES:
        if filename in missing:
            continue
        compare_csv(reference_directory / filename, regenerated_directory / filename,
                    optimization=filename == "causal_second_order_certified_range.csv",
                    desi=filename == "desi_dr1_tangent_scan.csv")

    for json_filename in JSON_FILES:
        if json_filename in missing:
            continue
        with (reference_directory / json_filename).open(encoding="utf-8") as stream:
            reference_json = json.load(stream)
        with (regenerated_directory / json_filename).open(encoding="utf-8") as stream:
            regenerated_json = json.load(stream)
        if json_filename == "causal_amplitude_reachability_results.json":
            compare_causal_amplitude(reference_json, regenerated_json)
            continue
        if json_filename == "bounded_joint_history_results.json":
            compare_bounded_joint_history(reference_json, regenerated_json)
            continue
        if json_filename == "gapped_residue_inverse_diagnostic_results.json":
            compare_inverse_diagnostic(reference_json, regenerated_json)
            continue
        compare_json(
            reference_json,
            regenerated_json,
            optimization=json_filename in PLATFORM_SENSITIVE_JSON_FILES,
            desi=json_filename == "desi_dr1_tangent_summary.json",
        )
    count = len(FILES) - len(missing)
    if missing:
        print(f"Legacy CI subset: {count} of {len(FILES)} scientific outputs "
              "agree within the declared numerical tolerances; "
              "this is partial coverage. Missing: " + ", ".join(sorted(missing)))
    else:
        print(f"All {len(FILES)} scientific outputs agree within the declared "
              "numerical tolerances.")



if __name__ == "__main__":
    main()
