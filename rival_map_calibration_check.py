# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Vitantonio Castronuovo

"""Synthetic check of the ordered rival-family decision contract.

This executable validates the decision states and illustrates why a mixture
quantile need not control each member of a composite family. It is not a
survey likelihood or a rival-family exclusion.
"""

import json

import numpy as np


SEED = 260926
CALIBRATION_COUNT = 100_000
EVALUATION_COUNT = 100_000
FAMILYWISE_ALPHA = 0.05
NUMBER_OF_FAMILIES = 2
ROW_ALPHA = FAMILYWISE_ALPHA / NUMBER_OF_FAMILIES
SHADOW_ALPHA = 0.05
FAMILY_MEANS = (-4.0, 4.0)
FAMILY_SIGMAS = (1.0, 1.6)
MINIMUM_POWER = 0.2


def scores(y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Profile the two-point rival family and evaluate the shadow score."""
    q_family = np.minimum((y - FAMILY_MEANS[0]) ** 2,
                          (y - FAMILY_MEANS[1]) ** 2)
    return q_family, y ** 2


def upper_tail_cutoff(values: np.ndarray, tail_probability: float) -> float:
    """Conservative empirical upper-tail quantile for a continuous toy null."""
    return float(np.quantile(values, 1.0 - tail_probability, method="higher"))


def decide(q_family: float, q_shadow: float, family_cutoff: float,
           shadow_cutoff: float, *, available: bool,
           power: float) -> str:
    """Apply the ordered N/I/S/E/J rule in the article."""
    if not available:
        return "N"
    if power < MINIMUM_POWER:
        return "I"
    if q_family <= family_cutoff:
        return "S"
    if q_shadow <= shadow_cutoff:
        return "E"
    return "J"


def run() -> dict:
    rng = np.random.default_rng(SEED)
    calibration = [scores(mean + sigma * rng.standard_normal(CALIBRATION_COUNT))[0]
                   for mean, sigma in zip(FAMILY_MEANS, FAMILY_SIGMAS)]
    pointwise_cutoffs = [upper_tail_cutoff(q, ROW_ALPHA) for q in calibration]
    uniform_cutoff = max(pointwise_cutoffs)
    mixture_cutoff = upper_tail_cutoff(np.concatenate(calibration), ROW_ALPHA)
    shadow_calibration = scores(rng.standard_normal(CALIBRATION_COUNT))[1]
    shadow_cutoff = upper_tail_cutoff(shadow_calibration, SHADOW_ALPHA)

    evaluation = []
    for mean, sigma in zip(FAMILY_MEANS, FAMILY_SIGMAS):
        q_family, q_shadow = scores(mean + sigma * rng.standard_normal(
            EVALUATION_COUNT))
        evaluation.append({
            "mean": mean,
            "sigma": sigma,
            "family_rejection_uniform": float(np.mean(q_family > uniform_cutoff)),
            "family_rejection_mixture": float(np.mean(q_family > mixture_cutoff)),
            "false_E_uniform": float(np.mean((q_family > uniform_cutoff) &
                                             (q_shadow <= shadow_cutoff))),
        })

    q_family_alt, q_shadow_alt = scores(rng.standard_normal(EVALUATION_COUNT))
    power = float(np.mean((q_family_alt > uniform_cutoff) &
                          (q_shadow_alt <= shadow_cutoff)))
    assert power >= MINIMUM_POWER
    truth_table = {
        "N": decide(0, 0, uniform_cutoff, shadow_cutoff,
                    available=False, power=power),
        "I": decide(0, 0, uniform_cutoff, shadow_cutoff,
                    available=True, power=0.0),
        "S": decide(0, 0, uniform_cutoff, shadow_cutoff,
                    available=True, power=power),
        "E": decide(uniform_cutoff + 1, 0, uniform_cutoff,
                    shadow_cutoff, available=True, power=power),
        "J": decide(uniform_cutoff + 1, shadow_cutoff + 1,
                    uniform_cutoff, shadow_cutoff,
                    available=True, power=power),
    }
    assert all(expected == actual for expected, actual in truth_table.items())
    assert all(row["family_rejection_uniform"] < ROW_ALPHA + 0.003
               for row in evaluation)
    assert evaluation[1]["family_rejection_mixture"] > ROW_ALPHA + 0.003

    return {
        "status": "synthetic contract check; not an observational exclusion",
        "seed": SEED,
        "calibration_count_per_null": CALIBRATION_COUNT,
        "evaluation_count_per_null": EVALUATION_COUNT,
        "familywise_alpha": FAMILYWISE_ALPHA,
        "number_of_families": NUMBER_OF_FAMILIES,
        "bonferroni_row_alpha": ROW_ALPHA,
        "shadow_alpha": SHADOW_ALPHA,
        "uniform_cutoff": uniform_cutoff,
        "prior_mixture_cutoff": mixture_cutoff,
        "evaluation": evaluation,
        "protected_alternative_E_power": power,
        "truth_table": truth_table,
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
