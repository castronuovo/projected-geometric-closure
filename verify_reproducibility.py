# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Vitantonio Castronuovo

"""Compare regenerated scientific outputs with frozen reference files."""

import csv
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
)
FILES = CSV_FILES + JSON_FILES
V15_FILES_NOT_IN_LEGACY_CI = frozenset((
    "desi_dr1_tangent_scan.csv",
    "spectral_shadow_duality.json",
    "desi_dr1_tangent_summary.json",
))
PLATFORM_SENSITIVE_JSON_FILES = (
    "causal_transport_support_results.json",
    "causal_transport_free_amplitude_results.json",
)


def close(reference: float, regenerated: float, *, optimization: bool = False,
          desi: bool = False) -> bool:
    if math.isnan(reference) or math.isnan(regenerated):
        return math.isnan(reference) and math.isnan(regenerated)
    return math.isclose(
        reference,
        regenerated,
        rel_tol=(DESI_RELATIVE_TOLERANCE if desi else
                 OPTIMIZATION_RELATIVE_TOLERANCE if optimization else
                 RELATIVE_TOLERANCE),
        abs_tol=(DESI_ABSOLUTE_TOLERANCE if desi else
                 OPTIMIZATION_ABSOLUTE_TOLERANCE if optimization else
                 ABSOLUTE_TOLERANCE),
    )


def compare_csv(reference_path: Path, regenerated_path: Path, *, desi=False) -> None:
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
                if not close(reference_number, regenerated_number, desi=desi):
                    raise AssertionError(
                        f"CSV numerical mismatch: {regenerated_path.name}:"
                        f"{row_index}:{column_index}"
                    )


def compare_json(reference, regenerated, path: str = "root", *, optimization=False,
                 desi=False) -> None:
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
                     desi=desi):
            raise AssertionError(f"JSON numerical mismatch at {path}")
        return
    if isinstance(reference, list):
        if not isinstance(regenerated, list) or len(reference) != len(regenerated):
            raise AssertionError(f"JSON list mismatch at {path}")
        for index, (reference_item, regenerated_item) in enumerate(
            zip(reference, regenerated)
        ):
            compare_json(reference_item, regenerated_item, f"{path}[{index}]",
                         optimization=optimization, desi=desi)
        return
    if isinstance(reference, dict):
        if not isinstance(regenerated, dict) or reference.keys() != regenerated.keys():
            raise AssertionError(f"JSON object mismatch at {path}")
        for key in reference:
            compare_json(reference[key], regenerated[key], f"{path}.{key}",
                         optimization=optimization, desi=desi)
        return
    raise TypeError(f"Unsupported JSON value at {path}")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: verify_reproducibility.py REFERENCE_DIRECTORY")
    reference_directory = Path(sys.argv[1]).resolve()
    regenerated_directory = Path(__file__).resolve().parent
    missing = {name for name in FILES if not (reference_directory / name).is_file()}
    if missing and missing != V15_FILES_NOT_IN_LEGACY_CI:
        raise FileNotFoundError(
            "Incomplete frozen-output set: " + ", ".join(sorted(missing))
        )

    for filename in CSV_FILES:
        if filename in missing:
            continue
        compare_csv(reference_directory / filename, regenerated_directory / filename,
                    desi=filename == "desi_dr1_tangent_scan.csv")

    for json_filename in JSON_FILES:
        if json_filename in missing:
            continue
        with (reference_directory / json_filename).open(encoding="utf-8") as stream:
            reference_json = json.load(stream)
        with (regenerated_directory / json_filename).open(encoding="utf-8") as stream:
            regenerated_json = json.load(stream)
        compare_json(
            reference_json,
            regenerated_json,
            optimization=json_filename in PLATFORM_SENSITIVE_JSON_FILES,
            desi=json_filename == "desi_dr1_tangent_summary.json",
        )
    if missing:
        print("Legacy CI subset: 14 of 17 scientific outputs agree within the "
              "declared numerical tolerances; three v1.5.0 outputs were "
              "not regenerated by this workflow.")
    else:
        print("All 17 scientific outputs agree within the declared "
              "numerical tolerances.")


if __name__ == "__main__":
    main()
