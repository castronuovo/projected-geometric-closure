# Changelog

## 1.2.1 — 2026-09-05

- Corrected reproducibility metadata to distinguish the 321-point logarithmic mass grid from the exact reference mass, for 322 distinct templates.
- Updated the injection--recovery contract and SHA-256 manifest accordingly.
- Clarified package documentation without changing stored benchmark outputs or numerical results.

## 1.2.0 — 2026-08-25

- Added an exactly soluble local five-dimensional scalar benchmark with
  numerical validation of the generalized Robin spectrum, positive boundary
  residues, critical light-pole expansion, finite-residue sum rule, heavy gap,
  and bounded tower remainder.
- Added sharp gap-certified lower and upper bounds on the total unweighted
  heavy residue, their simultaneous-band and uncertain-light-pole envelopes,
  the jointly sharp generalized-moment programs with finite atomic
  certificates, and the conditional critical inversion to \(r_b/L\).
- Added machine-readable 5D spectrum and validation outputs together with a
  diagnostic figure and continuous-integration reproduction checks.
- Added controlled synthetic injection--recovery mocks with independent
  calibration and evaluation ensembles, exact nuisance injection and
  profiling, a fixed dual-cone witness, and a window-mismatch stress test.
- Added machine-readable conditional null-point size, fixed-alternative
  rejection-probability, mass-recovery, and nominal profile-interval-coverage
  summaries over six injected target norms.
- Added the injection--recovery power figure and continuous-integration
  reproduction checks for its CSV and JSON outputs.
- Added a Lipschitz mesh-to-continuum error bound and a sufficient robust
  margin for promoting finite-grid dual witnesses to continuum certificates.
- Added a four-level positive-cone mesh-convergence audit without injected
  support nodes and stored its results in a machine-readable CSV and JSON.
- Added a nuisance-hardened positive spectral cone, non-negative least-squares
  class distances, and stored dual certificates of cone exclusion.
- Added the finite-data conic Carathéodory identifiability obstruction,
  clarifying that finite linear
  observables cannot uniquely establish a continuum or higher-dimensional
  spectral ontology.
- Distinguished the exact nuisance-orthogonal ambient dimension and its
  Carathéodory atom bound from the tolerance-dependent effective cone rank.
- Added positive-semidefinite Loewner checks, the exact Cauchy-factorization
  rank theorem, deterministic eigenvalue perturbation margins, and numerical
  rank diagnostics for ideal deconvolved one- and two-mode responses.
- Added an exact normalized spectral-transport identity and validation for a
  fixed physical spectral measure across epochs.
- Added an independent executable that maps each stored structural residual
  fraction to the nuisance-hardened target information required for a generic
  calibrated quadratic threshold.
- Added machine-readable conditional requirements for selectable thresholds.
- Added a conditional-information figure and continuous-integration coverage.
- Redesigned the conceptual flow figure with a clearer numbered hierarchy and
  explicit assumption labels while preserving the original scientific logic.
- Added manuscript-facing advance and decision maps while preserving the
  distinction between a structural scaling and a survey forecast.

## 1.1.2 — 2026-08-13

- Corrected the benchmark description to identify the stored residuals as
  fractions of the nuisance-hardened squared covariance-weighted norm.
- Added the corresponding residual norm ratios to the documentation.
- Updated the survey-benchmark figure label without changing its numerical
  inputs or stored scientific outputs.

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
