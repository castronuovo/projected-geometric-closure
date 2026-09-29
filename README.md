# Spectral identifiability of a protected cosmological response

This directory contains the reproducibility materials for
*Spectral Identifiability of a Protected Cosmological Response:
Conditional Exclusion and the Limits of a Two-Time Shadow*.

Repository: <https://github.com/castronuovo/projected-geometric-closure>

Current public release: **v1.5.0**. The release tag fixes the source,
machine-readable outputs, figures, pinned dependencies, and SHA-256 manifest
for this reproducibility package. Earlier releases remain available.

## Scientific architecture

The package implements the Paper-III inference bridge in the companion
program:

1. Paper I defines the protected late-time single-pole kernel.
2. Paper II supplies its observational entry and falsification gates.
3. Paper III applies progressive exclusion to test compatibility with a
   source-visible positive response class. Exact continuum agreement fixes its
   cyclic minimal operator up to unitary equivalence; finite projected data
   do not establish that exact identification.
4. A separate two-time theory develops an explicit (4+2)-dimensional parent
   realization of the surviving class.

The claimed uniqueness is conditional on the predeclared rival bank and the
retained observable contract. A theory reproducing the complete operator and
cross-channel signature belongs to the same observable equivalence class; the
package does not convert that equivalence into an unconditional proof of the
parent ontology.

The numerical and symbolic support has six layers:

1. an exactly soluble local five-dimensional scalar interval benchmark that
   validates the generalized Robin spectrum, positive boundary residues,
   finite-residue sum rule, critical light-pole expansion, heavy gap, and
   bounded tower remainder;
2. analytic tests of positive spectral measures, including the Stieltjes
   hierarchy, normalized Jensen envelopes, the variance identity, distinct
   slope- and ratio-defined turnover estimators, controlled single-pole
   bounds, sharp gap-certified inversion of the total heavy residue, the
   finite-scale chord test, Loewner positivity, exact
   spectral-complexity rank, deterministic perturbation margins, and normalized
   fixed-physical-spectrum transport;
3. a fixed survey-inspired response-space benchmark with three redshift slices,
   logarithmic scale bins, a non-diagonal window matrix, correlated covariance,
   five nuisance columns spanning three smooth nuisance families, and a
   finite-grid positive spectral cone with stored dual certificates, a
   four-level mesh-convergence audit, and nested stationary versus
   epoch-separated cones testing cross-epoch spectral transport;
4. a controlled synthetic injection--recovery layer with independent
   calibration and evaluation ensembles, a fixed dual-cone witness, and a
   predeclared window-mismatch stress test; and
5. an executable conditional-information calculator that translates the
   structural residual fractions into the target norm required for a generic,
   externally calibrated quadratic threshold; and
6. a finite-data duality construction showing that distinct positive spectral
   measures can agree on the retained samples while differing as continuum
   Stieltjes responses.

## Two-time-parent inference scope

The manuscript additionally derives an analytic result that does not require a
numerical pipeline. Every single-field Einstein-frame potential admits an
exact homogeneous two-scalar Weyl lift. The lift is therefore compatible with
the inherited infrared kernel but is not identifiable from that kernel,
mediator mass, or reduced potential alone. The associated matter-coupling
analysis yields three distinct branches: universal screened coupling,
species-dependent dark coupling, and locally constrained weak universal
coupling. None is silently identified with the universal kernel.

Independent Wolfram checks and the maintained Lean theorem remain in the
author's local development workspace. The graphical-abstract source and
generator are also kept with the manuscript assets, outside this package.
The two legacy scripts under the manuscript's `script/` directory deliberately
raise an obsolete-benchmark error and are not part of this package.

The late-time scale and redshift coverage is DESI-like at the level of
binning. No released DESI tracer selection, window, covariance, catalogue, or
likelihood is used.

That statement applies to the original synthetic benchmark. The separate
`desi_dr1_tangent_test.py` development demonstrator now ingests official DESI
DR1 LRG power-spectrum data vectors, windows, and covariance products. It is
not used to replace or retrospectively reinterpret the frozen synthetic
results.

## Restricted DESI DR1 tangent demonstrator

Run `python desi_dr1_tangent_test.py` to download the three official DESI DR1
combined-footprint LRG power-spectrum likelihood containers, verify their
project-recorded SHA-256 hashes, and regenerate
`desi_dr1_tangent_summary.json`, `desi_dr1_tangent_scan.csv`, and
`figs/figS7_desi_dr1_tangent_test.png`.

The calculation uses a fixed-background CAMB linear spectrum, integrates the
scale-dependent growth equation from matched data at `a=0.25`, propagates
linear Kaiser multipoles through the released windows, and profiles five local
nuisance directions per redshift bin. The nominal reference fits have
chi-squared values `18.72`, `24.77`, and `27.55` for 27 nominal degrees of
freedom. The score near `m_star = 0.03 h/Mpc` changes from `-1.69` to about
`-2.47` across the predeclared `k_max = 0.08, 0.09, 0.10 h/Mpc` audit, and the
unconstrained estimate lies outside the validated local tangent neighborhood.
Consequently the stored output is an inapplicability diagnostic, not a
detection, exclusion, posterior, or DESI Collaboration result. It does not
test `Sigma = 1` because no lensing likelihood is included.

The second layer compares a protected single mode, a positive two-mode
mixture, and a signed response outside the positive-spectrum class. It is a
structural identifiability benchmark. It does not ingest DESI or KiDS data, run
a Boltzmann solver, integrate a growth equation, map the kernel response to a
named survey observable, evaluate a likelihood, or provide a survey forecast.
The single-mode fit profiles a non-negative amplitude over a 321-point logarithmic grid plus the exact single-target mass (322 distinct templates) in
`0.100--0.800 h Mpc^-1`, within the inherited quasi-static benchmark range.

## Causal matter-power extension

Run `python causal_growth_transport.py` to reproduce
`causal_growth_transport_benchmark.csv`,
`causal_growth_transport_summary.json`, and
`figs/figS6_causal_growth_transport.png`.

This additional example integrates the linear-density growth equation on a
fixed flat matter–Lambda background and computes `P_m/P_m_GR - 1` with
identical matched initial data. Two non-negative temporal source blocks act
**before** growth propagation. A two-mass target produces nested
stationary and independent-history tangent cones; one scale-independent
amplitude nuisance is removed per output epoch. The covariance and window
are synthetic. The script does not compute Kaiser multipoles or lensing.

A uniform Neumann remainder controls the transformation from growth to
power for a declared bound on total response amplitude. Distance to the
unbounded stationary tangent cone, minus that error allowance, lower-bounds
distance to the fully integrated amplitude-restricted stationary class.
At amplitude `0.0003` the margin is positive; at `0.003` the conservative
bound is inconclusive. These are approximation diagnostics, not detection
significances.

Run `python causal_transport_support_audit.py` to extend the stationary
comparison to arbitrary positive mass support, including a continuum, with
explicit grid and unresolved-tail allowances. For amplitude `0.0003`, the
certified lower margin is `0.00062252` for the step target and `0.000573451`
for a smooth weight transition.

Run `python causal_transport_free_amplitude.py` for the adversarial comparison
in which one stationary normalized spectrum has a separate non-negative total
amplitude in each source interval. This freedom reduces the tangent residual
to `0.000382366`; the applicable uniform error allowance is `0.000406502`.
The comparison is therefore inconclusive. It does not establish compatibility
or robust discrimination with a free amplitude history.
The JSON records the complete contract and all scan points.

Validation includes source-aligned fourth-order Runge–Kutta step halving,
stationary tangent recovery, independent-history recovery, cone inclusion,
NNLS KKT conditions, and comparison of actual errors with the analytic bound.
Floating-point convergence is not interval-arithmetic certification.

## Environment and execution

The additional `rival_map_calibration_check.py` is a deterministic synthetic
contract check for the ordered N/I/S/E/J decision rule. It compares a
least-favorable composite-family cutoff with a prior-mixture cutoff and
evaluates the former on independent draws. It does not use DESI data, test a
physical rival family, or produce an observational exclusion.

The stored output was generated with:

- Python 3.9.6
- NumPy 1.26.4
- Matplotlib 3.9.4
- SciPy 1.13.1 (robustness audits only)
- h5py 3.12.1 and CAMB 1.6.0 (DESI tangent demonstrator only)

Install the recorded dependencies with:

```bash
python3 -m pip install -r requirements.txt
```

## Reproduction

Run from this directory:

```bash
python3 benchmark_projected_spectral.py
python3 spectral_shadow_duality.py
python3 finite_residue_5d_benchmark.py
python3 injection_recovery.py
python3 conditional_information_requirement.py
python3 causal_growth_transport.py
python3 causal_transport_support_audit.py
python3 causal_transport_free_amplitude.py
python3 desi_dr1_tangent_test.py
python3 rival_map_calibration_check.py
```

The DESI command downloads official released likelihood containers into a
local cache and verifies their recorded hashes. These input files are not
redistributed with this package.

To compare a complete regenerated run with the frozen release outputs, copy
the frozen CSV and JSON files to a separate directory before running the nine
benchmark commands above, then run
`python3 verify_reproducibility.py FROZEN_DIRECTORY`. The synthetic contract
check is self-validating and has no frozen scientific output.
The verifier covers all nine CSV and eight JSON scientific outputs, including
the DESI diagnostic. Its DESI-specific tolerances are `3e-5` relative and
`1e-6` absolute because profiled nuisance optima vary slightly across
independent numerical runs; structural metadata and integer counts still
match exactly. This comparison does not validate the linear Kaiser forward
model against survey mocks or turn the DESI score into a constraint.
The public GitHub workflow still regenerates the 14 inherited outputs. When
given that exact legacy reference set, the verifier reports the partial scope
explicitly. A complete v1.5.0 replay requires running the three additional
scripts listed above and retaining all 17 frozen outputs before comparison.

Optional thresholds can be supplied explicitly, for example:

```bash
python3 conditional_information_requirement.py --thresholds 1,4,9,16,25
```

These commands regenerate:

- `benchmark_identifiability.csv`
- `spectral_shadow_duality.json`
- `cone_grid_convergence.csv`
- `cross_epoch_transport_benchmark.csv`
- `finite_residue_5d_spectrum.csv`
- `finite_residue_5d_benchmark.json`
- `survey_projected_benchmark.csv`
- `survey_projected_benchmark.json`
- `conditional_information_requirement.csv`
- `injection_recovery_summary.csv`
- `injection_recovery_summary.json`
- `causal_growth_transport_benchmark.csv`
- `causal_growth_transport_summary.json`
- `causal_transport_support_results.json`
- `causal_transport_free_amplitude_results.json`
- `figs/fig1_geometric_spectral_test.png`
- `figs/fig2_identifiability_benchmark.png`
- `figs/fig3_survey_projected_spectral_benchmark.png`
- `figs/fig4_spectral_shadow_duality.png`
- `figs/figS1_finite_scale_spectral_diagnostics.png`
- `figs/figS2_conditional_information_requirement.png`
- `figs/figS3_injection_recovery.png`
- `figs/figS4_finite_residue_5d_benchmark.png`
- `figs/figS5_cross_epoch_transport.png`
- `figs/figS6_causal_growth_transport.png`
- `desi_dr1_tangent_scan.csv`
- `desi_dr1_tangent_summary.json`
- `figs/figS7_desi_dr1_tangent_test.png`

The benchmark script stops if the window rows are not normalized, if the covariance is
not positive definite, or if a declared spectral inequality is violated
beyond the recorded numerical tolerance. The JSON records the fixed analysis
contract, the minimum covariance eigenvalue, the nuisance rank, the
profiled-information eigenvalues, the single-mode residuals, and
machine-readable validation of the normalized spectral identities and
approximation bounds, Loewner matrices, spectral transport, positive-cone
separation, finite-grid convergence, and the causal-growth tangent remainder.

The dedicated five-dimensional benchmark independently records the exact
light mass and residue, the critical heavy gap, the finite-residue sum-rule
convergence, the truncated spectral reconstruction of the closed-form
resolvent, and the heavy-remainder inequality. Its scope is the explicit local
scalar interval model in the manuscript; it does not identify that mediator
with the full coupled gravitational master sector.

## Interpretation

The stored benchmark gives the following fractions of the nuisance-hardened
squared covariance-weighted norm (quadratic information) left by the best
joint fit of a positive single-mode template and the declared nuisance
columns:

- protected single mode: numerical zero;
- positive two-mode mixture: `0.030`;
- signed out-of-class response: `0.174`.

The corresponding residual norm ratios are approximately `0.172` and `0.417`
for the positive two-mode and signed targets. These values are conditional on the fixed synthetic window, covariance,
redshift slices, scale range, and nuisance basis. They demonstrate
class-conditional projected separation in this benchmark only. They are not
detection probabilities or statements about current survey sensitivity.

The same three targets are tested against a 323-generator positive spectral
cone after identical whitening and nuisance hardening. The single-mode and
positive two-mode targets are compatible to numerical tolerance. The signed
target retains a cone-distance fraction `0.174090`; its residual provides a
dual witness whose product with every cone generator is non-negative to the
stored tolerance while its product with the target is negative. This is a
finite-grid convex-cone validation, not evidence that the data resolve a
continuous spectrum.

The nuisance-orthogonal ambient data space has exact dimension
\(54-5=49\); this gives the exact conic-Caratheodory bound of at most 49
atoms for an exactly compatible projected cone point. The sampled cone design
has effective numerical rank 16 at the declared absolute singular-value
tolerance \(10^{-10}\). That tolerance-dependent value is retained only as a
compression diagnostic and is not used as an exact atom bound.

An independent mesh audit repeats the cone projection on 81, 161, 321, and
641 logarithmically spaced masses without inserting the target support. The
signed-target cone-distance fraction converges from `0.1740943` to
`0.1740901`, while the positive-target residuals tend to numerical zero. The
stored covering radii and residual fractions diagnose finite-grid stability;
they are not, by themselves, a continuum dual certificate. The latter also
requires the Lipschitz margin derived in the Supplemental Material.

For a generic quadratic threshold `T`, the second executable evaluates

```text
Delta chi2_target|X = T / f_res|X
required target norm = sqrt(T / f_res|X).
```

For example, `T=4` requires target norms `11.615` and `4.793` for the positive
two-mode and signed targets, respectively. `T` is not interpreted as a
Gaussian significance. A concrete analysis must calibrate it with its own
boundary conditions, mass scan, look-elsewhere prescription, and end-to-end
mocks.

## Controlled injection--recovery calibration

`injection_recovery.py` draws correlated Gaussian noise from the fixed
synthetic covariance, injects random amplitudes of all five declared nuisance
columns, and adds one of four predeclared response-space targets. Every mock is
then passed through the same whitening, nuisance projection, non-negative
amplitude fit, and 321-point logarithmic grid plus exact reference mass as the deterministic benchmark. The
configuration in `injection_recovery_config.json` fixes the random seed, 5000
calibration realizations, 5000 statistically independent evaluation
realizations, a 95-percent decision quantile, and target norms
`0, 4, 8, 12, 20, 30`.

The best-single-mode lack-of-fit threshold is calibrated independently at
each injected target norm using a fixed reference null mass
\(m=0.150\,h\,{\rm Mpc}^{-1}\). Its conditional null-point rejection rate stays
between `0.051` and `0.057`; the nominal profile-Delta-chi-squared-one mass
interval covers the injected mass in `0.669--0.692` of the nonzero-signal
realizations. Changing the injection window width from `0.105` to `0.125`
while retaining the nominal recovery window leaves the rejection rate below
`0.063` over the stored norm grid. This is a limited response-space mismatch
test, not a validation of an observational forward model.

This executable calibrates the absolute lack-of-fit statistic \(q_1\) at the
stated null point. It does not establish uniform size over the composite
one-mode class and does not evaluate or calibrate the composite likelihood
ratio \(q_{\rm spec}\) defined as the target statistic for a survey-level
analysis.

For the positive two-mode target, the single-mode rejection probability is
`0.130`, `0.311`, and `0.758` at target norms `12`, `20`, and `30`. A fixed
dual-cone witness, constructed before the Monte Carlo draws and calibrated on
an independent noise-only ensemble, has evaluation null rejection rate
`0.0448`. For the signed out-of-class target its rejection probabilities are
`0.495`, `0.949`, and `0.999` at norms `4`, `8`, and `12`. It remains
conservative for the two compatible positive targets because they lie inside,
rather than on the least-favorable boundary of, the tested half-space.
With 5000 evaluation realizations, the binomial Monte Carlo standard error is
at most `0.0071` for every reported rejection probability.

These numbers calibrate the statistical behavior of the frozen synthetic
response-space contract. They do not ingest a catalogue, propagate a named
survey selection, or establish observational sensitivity. A survey analysis
must rebuild and recalibrate the mocks with its complete growth--Weyl forward
model, windows, covariance, nuisance basis, scale cuts, and likelihood.

## Normalized spectral validation

The `validation.normalized_spectral_diagnostics` object in
`survey_projected_benchmark.json` records, for the one-mode and positive
two-mode benchmarks:

- the finite total residue and the harmonic/arithmetic endpoint scales;
- maximum violations of the two Jensen bounds;
- the residual in the exact variance identity;
- the bounds `0 <= W <= 1`;
- monotonicity checks for `t_slope` and `t_ratio`;
- equality of the two estimators for a single pole;
- dominant-pole, heavy-support, and narrow-spectrum error-bound checks;
- positive-semidefinite Loewner matrices, their exact finite-atom rank
  interpretation, and numerical ranks;
- the fixed normalized physical-spectrum transport residual.

The separate `validation.positive_spectral_cone` object records the projected
cone rank, mass grid, non-negative least-squares distances, active weights,
and dual-certificate checks.

These are deterministic analytic-grid checks, not observational constraints.
The mathematical results remain conditional on the hypotheses stated in the
article and Supplemental Material.

## Finite-residue five-dimensional benchmark

`finite_residue_5d_benchmark.py` solves the generalized Robin eigenproblem
for the declared dimensionless contract \(L=1\), \(r_b/L=12\), and
\(\delta h\,L=4\times10^{-3}\). It checks:

- positivity of every retained boundary residue;
- the exact normalized sum rule \(\sum_n w_n=1\);
- the first- and second-order critical light-mass expansions;
- the critical heavy-gap equation and interval;
- reconstruction of the closed-form boundary resolvent with 512 modes; and
- positivity and the analytic gap bound for the heavy-sector remainder.

The stored benchmark gives
\(m_\ast^2L^2=3.24324135\times10^{-4}\),
\(w_\ast=0.97297184\), and
\(\eta^2L^2=10.0341804\). The 512-mode spectral reconstruction agrees with
the exact resolvent to better than \(10^{-11}\) over the declared momentum
window. These numbers validate the analytic benchmark and are not fitted
physical parameters.

For any declared decomposition into one identified light pole and a positive
heavy measure supported above the certified gap, the normalized response also
provides sharp one-scale lower and upper bounds on the total heavy residue.
Their intersection across scales is either a conservative allowed interval or
an empty-set falsification of that subclass. In the critical limit only, this
interval maps monotonically to a conditional bound on \(r_b/L\). The proof and
the envelope for an uncertain light-pole location are given in the
Supplemental Material; they do not require an upper cutoff on the tower. The
jointly sharp multiscale endpoints are generalized linear moment programs and
admit finite atomic certificates with at most one more atom than the number of
sampled scales.

## Integrity

Run:

```bash
shasum -a 256 -c SHA256SUMS
```

to verify the frozen source, results, and figures. Rendering hashes can change
with the Matplotlib or font-stack version; the JSON and CSV files are the
machine-readable scientific outputs.

The current frozen public package is release
[`v1.5.0`](https://github.com/castronuovo/projected-geometric-closure/releases/tag/v1.5.0),
published on 29 September 2026. It adds the DESI tangent demonstrator,
spectral-duality benchmark, rival-map calibration check, and their applicable
machine-readable outputs. Earlier tags, including
[`v1.4.1`](https://github.com/castronuovo/projected-geometric-closure/releases/tag/v1.4.1),
remain immutable and are not overwritten.

## License

- `benchmark_projected_spectral.py` is licensed under the
  [BSD 3-Clause License](LICENSE).
- `injection_recovery.py` and `conditional_information_requirement.py` are
  licensed under the [BSD 3-Clause License](LICENSE).
- `finite_residue_5d_benchmark.py` is licensed under the
  [BSD 3-Clause License](LICENSE).
- The documentation, CSV and JSON outputs, and figures are licensed under
  [Creative Commons Attribution 4.0 International](LICENSE-CONTENT).

Copyright (c) 2026 Vitantonio Castronuovo.
