# Changelog

## 1.1.1 — 2026-08-13

- Replaced byte-exact cross-platform output comparison with an explicit
  numerical-tolerance audit while preserving exact structural checks.
- Added a standalone reproducibility verifier for the frozen CSV and JSON
  scientific outputs.
- Updated the workflow actions to their Node 24 generations.

This patch does not alter the scientific model, benchmarks, or manuscript
claims.

## 1.1.0 — 2026-08-13

- Added machine-readable checks of the normalized Jensen envelopes.
- Added a direct numerical check of the exact spectral-variance identity.
- Added bounds on the differential width `0 <= W <= 1`.
- Distinguished and tested the slope- and ratio-defined turnover estimators.
- Added dominant-pole, heavy-support, and narrow-spectrum error-bound checks.
- Renamed the spectral-figure turnover axis from `t_eff` to `t_slope`.
- Added continuous-integration reproduction checks.

This version extends the analytic audit only. It does not introduce
observational data, a Boltzmann solver, a survey likelihood, or a forecast.

## 1.0.0 — 2026-08-01

- Frozen initial public reproducibility release for the projected geometric
  closure manuscript.
- Added the analytic spectral figures and the survey-inspired projected
  identifiability benchmark.
