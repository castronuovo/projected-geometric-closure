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
FILES = (
    "benchmark_identifiability.csv",
    "survey_projected_benchmark.csv",
    "survey_projected_benchmark.json",
)


def close(reference: float, regenerated: float) -> bool:
    if math.isnan(reference) or math.isnan(regenerated):
        return math.isnan(reference) and math.isnan(regenerated)
    return math.isclose(
        reference,
        regenerated,
        rel_tol=RELATIVE_TOLERANCE,
        abs_tol=ABSOLUTE_TOLERANCE,
    )


def compare_csv(reference_path: Path, regenerated_path: Path) -> None:
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
                if not close(reference_number, regenerated_number):
                    raise AssertionError(
                        f"CSV numerical mismatch: {regenerated_path.name}:"
                        f"{row_index}:{column_index}"
                    )


def compare_json(reference, regenerated, path: str = "root") -> None:
    if isinstance(reference, bool) or reference is None or isinstance(reference, str):
        if reference != regenerated:
            raise AssertionError(f"JSON value mismatch at {path}")
        return
    if isinstance(reference, (int, float)):
        if not isinstance(regenerated, (int, float)) or isinstance(regenerated, bool):
            raise AssertionError(f"JSON type mismatch at {path}")
        if not close(float(reference), float(regenerated)):
            raise AssertionError(f"JSON numerical mismatch at {path}")
        return
    if isinstance(reference, list):
        if not isinstance(regenerated, list) or len(reference) != len(regenerated):
            raise AssertionError(f"JSON list mismatch at {path}")
        for index, (reference_item, regenerated_item) in enumerate(
            zip(reference, regenerated)
        ):
            compare_json(reference_item, regenerated_item, f"{path}[{index}]")
        return
    if isinstance(reference, dict):
        if not isinstance(regenerated, dict) or reference.keys() != regenerated.keys():
            raise AssertionError(f"JSON object mismatch at {path}")
        for key in reference:
            compare_json(reference[key], regenerated[key], f"{path}.{key}")
        return
    raise TypeError(f"Unsupported JSON value at {path}")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: verify_reproducibility.py REFERENCE_DIRECTORY")
    reference_directory = Path(sys.argv[1]).resolve()
    regenerated_directory = Path(__file__).resolve().parent

    for filename in FILES[:2]:
        compare_csv(reference_directory / filename, regenerated_directory / filename)

    json_filename = FILES[2]
    with (reference_directory / json_filename).open(encoding="utf-8") as stream:
        reference_json = json.load(stream)
    with (regenerated_directory / json_filename).open(encoding="utf-8") as stream:
        regenerated_json = json.load(stream)
    compare_json(reference_json, regenerated_json)
    print("Scientific outputs agree within the declared numerical tolerances.")


if __name__ == "__main__":
    main()
