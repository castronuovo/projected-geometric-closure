# Changelog

## 1.4.1 — 2026-09-11

- Made the declared logarithmic mass grids platform-independent by replacing near-coincident anchor nodes rather than relying on bitwise deduplication.
- Preserved exact comparison of integer metadata while applying explicit cross-platform tolerances only to the two SciPy optimization outputs introduced in v1.4.0.
- This portability patch does not alter the model, analysis contract, qualitative conclusions, or manuscript claims.

## 1.4.0 — 2026-09-11

- Implemented a fixed-background growth integration with temporal spectral weights applied before causal propagation, and a tangent map to the fractional matter-power response.
- Added a uniform power-response remainder and a conservative distance bound for the amplitude-restricted stationary class on a two-mass dictionary. Stored both informative and inconclusive amplitude regimes.
- Extended the stationary comparison to arbitrary positive mass support with continuum-grid and unresolved-tail allowances, and verified that the controlled small-amplitude separation survives a smooth weight transition.
- Added an adversarial stationary-shape fit with independent amplitudes in two source intervals. The resulting residual is below the applicable uniform error allowance, so this broader comparison is explicitly reported as inconclusive.
- Added exact-versus-tangent checks, source-aligned RK4 step halving, nested-cone and NNLS KKT checks, CSV/JSON outputs, and the causal-growth figure.
- Included both transport CSV files in the numerical verifier and frozen CI outputs, and added CFF validation.
- Corrected `cff-version` to the metadata schema version `1.2.0` and advanced the software package version to `1.4.0`.
- The new calculation is an illustrative matter-power benchmark, not a survey forecast or a validated tracer/lensing pipeline.

## 1.3.0 — 2026-09-06

- Added a nested-cone cross-epoch transport diagnostic that separates a shared stationary positive spectrum from independently positive epoch responses under a fixed amplitude history.
- Added stationary and positive weight-drift benchmark targets; the latter remains inside the epoch-separated cone while leaving a fractional stationary-cone gap of `0.155119`.
- Added machine-readable cross-epoch benchmark output and the corresponding diagnostic figure.
- Extended the reproducibility verifier, frozen analysis contract, documentation, and SHA-256 manifest to cover the new benchmark.
- Preserved the synthetic response-space scope: this release does not ingest observational data, run a Boltzmann solver, evaluate a survey likelihood, or provide a forecast.

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
