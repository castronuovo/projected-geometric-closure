# Projected geometric dark-sector closure

This directory contains the reproducibility materials for
*Projected Geometric Dark-Sector Closure: Five-Dimensional Boundary-Scalar
Benchmark, Single-Pole Stability, and Observable Spectral Complexity*.

Repository: <https://github.com/castronuovo/projected-geometric-closure>

The calculation has five layers:

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
   finite-grid positive spectral cone with stored dual certificates and a
   four-level mesh-convergence audit;
4. a controlled synthetic injection--recovery layer with independent
   calibration and evaluation ensembles, a fixed dual-cone witness, and a
   predeclared window-mismatch stress test; and
5. an executable conditional-information calculator that translates the
   structural residual fractions into the target norm required for a generic,
   externally calibrated quadratic threshold.

The late-time scale and redshift coverage is DESI-like at the level of
binning. No released DESI tracer selection, window, covariance, catalogue, or
likelihood is used.

The second layer compares a protected single mode, a positive two-mode
mixture, and a signed response outside the positive-spectrum class. It is a
structural identifiability benchmark. It does not ingest DESI or KiDS data, run
a Boltzmann solver, integrate a growth equation, map the kernel response to a
named survey observable, evaluate a likelihood, or provide a survey forecast.
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
python3 finite_residue_5d_benchmark.py
python3 injection_recovery.py
python3 conditional_information_requirement.py
```

Optional thresholds can be supplied explicitly, for example:

```bash
python3 conditional_information_requirement.py --thresholds 1,4,9,16,25
```

These commands regenerate:

- `benchmark_identifiability.csv`
- `cone_grid_convergence.csv`
- `finite_residue_5d_spectrum.csv`
- `finite_residue_5d_benchmark.json`
- `survey_projected_benchmark.csv`
- `survey_projected_benchmark.json`
- `conditional_information_requirement.csv`
- `injection_recovery_summary.csv`
- `injection_recovery_summary.json`
- `figs/fig0_conceptual_flow.png`
- `figs/fig1_geometric_spectral_test.png`
- `figs/fig2_identifiability_benchmark.png`
- `figs/fig3_survey_projected_spectral_benchmark.png`
- `figs/figS1_finite_scale_spectral_diagnostics.png`
- `figs/figS2_conditional_information_requirement.png`
- `figs/figS3_injection_recovery.png`
- `figs/figS4_finite_residue_5d_benchmark.png`

The benchmark script stops if the window rows are not normalized, if the covariance is
not positive definite, or if a declared spectral inequality is violated
beyond the recorded numerical tolerance. The JSON records the fixed analysis
contract, the minimum covariance eigenvalue, the nuisance rank, the
profiled-information eigenvalues, the single-mode residuals, and
machine-readable validation of the normalized spectral identities and
approximation bounds, Loewner matrices, spectral transport, and positive-cone
separation and finite-grid convergence.

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
amplitude fit, and 321-point mass scan as the deterministic benchmark. The
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

The release target corresponding to the revised manuscript is
[`v1.2.0`](https://github.com/castronuovo/projected-geometric-closure/releases/tag/v1.2.0).
It was published on 25 August 2026 and is the frozen reproducibility package
supporting the submitted manuscript.
The earlier `v1.0.0` and `v1.1.0` tags remain immutable and are not
overwritten.

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
