# Projected geometric closure

This directory contains the reproducibility materials for
*Projected Geometric Dark-Sector Closure: Positive-Spectrum Saturation and
Observable Identifiability*.

Repository: <https://github.com/castronuovo/projected-geometric-closure>

The calculation has two layers:

1. analytic tests of positive spectral measures, including the Stieltjes
   hierarchy, the single-mode saturation diagnostics, and the finite-scale
   chord test;
2. a fixed survey-inspired projection benchmark with three redshift slices,
   logarithmic scale bins, a non-diagonal window matrix, correlated covariance,
   and five nuisance columns spanning three smooth nuisance families.

The late-time scale and redshift coverage is DESI-like at the level of
binning. No released DESI tracer selection, window, covariance, catalogue, or
likelihood is used.

The second layer compares a protected single mode, a positive two-mode
mixture, and a signed response outside the positive-spectrum class. It is a
structural identifiability benchmark. It does not ingest DESI or KiDS data, run
a Boltzmann solver, evaluate a likelihood, or provide a survey forecast.
The single-mode fit profiles a non-negative amplitude over 321 masses in
`0.100--0.800 h Mpc^-1`, within the inherited quasi-static benchmark range.

## Environment

The stored output was generated with:

- Python 3.9.6
- NumPy 2.0.2
- Matplotlib 3.9.4

Install the recorded dependencies with:

```bash
python3 -m pip install -r requirements.txt
```

## Reproduction

Run from this directory:

```bash
python3 benchmark_projected_spectral.py
```

The command regenerates:

- `benchmark_identifiability.csv`
- `survey_projected_benchmark.csv`
- `survey_projected_benchmark.json`
- `figs/fig0_conceptual_flow.png`
- `figs/fig1_geometric_spectral_test.png`
- `figs/fig2_identifiability_benchmark.png`
- `figs/fig3_survey_projected_spectral_benchmark.png`
- `figs/figS1_finite_scale_spectral_diagnostics.png`

The script stops if the window rows are not normalized or if the covariance is
not positive definite. The JSON records the fixed analysis contract, the
minimum covariance eigenvalue, the nuisance rank, the profiled-information
eigenvalues, and the single-mode residuals.

## Interpretation

The stored benchmark gives the following fractions of the nuisance-hardened
signal left by the best joint fit of a positive single-mode template and the
declared nuisance columns:

- protected single mode: numerical zero;
- positive two-mode mixture: `0.030`;
- signed out-of-class response: `0.174`.

These values are conditional on the fixed synthetic window, covariance,
redshift slices, scale range, and nuisance basis. They demonstrate
class-conditional projected separation in this benchmark only. They are not
detection probabilities or statements about current survey sensitivity.

## Integrity

Run:

```bash
shasum -a 256 -c SHA256SUMS
```

to verify the frozen source, results, and figures. Rendering hashes can change
with the Matplotlib or font-stack version; the JSON and CSV files are the
machine-readable scientific outputs.

## License

- `benchmark_projected_spectral.py` is licensed under the
  [BSD 3-Clause License](LICENSE).
- The documentation, CSV and JSON outputs, and figures are licensed under
  [Creative Commons Attribution 4.0 International](LICENSE-CONTENT).

Copyright (c) 2026 Vitantonio Castronuovo.
