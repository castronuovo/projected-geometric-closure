# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Vitantonio Castronuovo
"""Causal, fixed-background matter-power transport benchmark; not a forecast.

Only NumPy and Matplotlib are required. The finite mass bank is deliberately
two-dimensional: this is a declared finite dictionary, not a continuum test.
Initial growth data are fixed at a=1/4; no primordial matching is inferred.
"""

import csv
import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
K = np.geomspace(0.03, 0.15, 12)  # h/Mpc, comoving
MASSES = np.array([0.1, 0.4])  # h/Mpc, physical inverse lengths
TIMES = np.linspace(np.log(0.25), 0.0, 9)[1:]
SPLIT = np.log(0.5)
AMPLITUDES = (0.0003, 0.001, 0.003, 0.01, 0.03)
SIGMA = 0.01  # arbitrary fractional-power normalization, not survey errors


def background(n):
    om = 0.3 * np.exp(-3 * n) / (0.3 * np.exp(-3 * n) + 0.7)
    return 2 - 1.5 * om, 1.5 * om


def evolve(weights, steps=2048):
    """RK4, aligning all output times and the source discontinuity with steps.

    State: GR D,D'; exact E=D/D_GR,E'; first-order u,u'. Each spectral
    weight is dimensionless and enters mu-1 directly, not through beta^2.
    """
    if steps < 16 or steps % 8:
        raise ValueError('steps must be a multiple of eight and at least 16')
    y = np.zeros((6, len(K)))
    y[0] = 0.25
    y[1] = 0.25  # fixed matching data; not an exact LambdaCDM growing mode
    y[2] = 1.0
    n0 = np.log(0.25)
    dt = -n0 / steps
    outputs = []

    def rhs(n, state, epoch):
        acoef, bcoef = background(n)
        d, dp, e, ep, u, up = state
        kernel = K[:, None] ** 2 / (
            K[:, None] ** 2 + np.exp(2*n) * MASSES[None, :] ** 2
        )
        h = kernel @ weights[epoch]
        damping = 2 * dp / d + acoef
        return np.array([dp, bcoef*d-acoef*dp, ep,
                         bcoef*h*e-damping*ep, up, bcoef*h-damping*up])

    for j in range(steps):
        n = n0 + j * dt
        # One-sided source on each step, including its endpoint evaluations.
        epoch = 0 if j < steps // 2 else 1
        k1 = rhs(n, y, epoch)
        k2 = rhs(n+dt/2, y+dt*k1/2, epoch)
        k3 = rhs(n+dt/2, y+dt*k2/2, epoch)
        k4 = rhs(n+dt, y+dt*k3, epoch)
        y += dt * (k1+2*k2+2*k3+k4)/6
        if (j+1) % (steps//8) == 0:
            outputs.append(y.copy())
    states = np.array(outputs)
    return (states[:, 2]**2-1).ravel(), (2*states[:, 4]).ravel()


def project(vector, basis):
    return (vector - basis @ (basis.T @ vector)) / SIGMA


def nnls_small(design, target):
    """Enumerate all faces (at most 16); verify primal/dual KKT conditions."""
    best = (float(target @ target), np.zeros(design.shape[1]))
    for size in range(1, design.shape[1]+1):
        for active in itertools.combinations(range(design.shape[1]), size):
            w = np.linalg.lstsq(design[:, active], target, rcond=None)[0]
            if np.min(w) < -1e-11:
                continue
            full = np.zeros(design.shape[1])
            full[list(active)] = np.maximum(w, 0)
            residual = target - design @ full
            q = float(residual @ residual)
            if q < best[0]:
                best = q, full
    residual = target - design @ best[1]
    scale = max(1.0, np.linalg.norm(design)*np.linalg.norm(target))
    assert np.max(design.T @ residual) < 1e-9 * scale
    assert np.max(np.abs(best[1] * (design.T @ residual))) < 1e-9 * scale
    return best[0], best[1]


def run(steps=2048):
    ndata = len(TIMES)*len(K)
    # One independent scale-independent amplitude nuisance per output epoch.
    nuisance = np.kron(np.eye(len(TIMES)), np.ones((len(K), 1)))
    basis = np.linalg.qr(nuisance, mode='reduced')[0]
    columns = []
    coarse_error = 0.0
    memory_amplitude = 0.0
    for epoch in range(2):
        for mass in range(2):
            w = np.zeros((2, 2)); w[epoch, mass] = 1
            _, tangent = evolve(w, steps)
            if epoch == 0 and mass == 0:
                memory_amplitude = float(np.max(tangent.reshape(8, 12)[-1]))
                assert memory_amplitude > 0
            _, coarse = evolve(w, steps//2)
            coarse_error = max(coarse_error, float(np.max(np.abs(tangent-coarse))))
            columns.append(project(tangent, basis))
    epoch_design = np.column_stack(columns)
    stat_design = epoch_design[:, :2] + epoch_design[:, 2:]
    _, stationary_unit = evolve(np.array([[1., 0.], [1., 0.]]), steps)
    nesting_error = float(np.max(np.abs(stat_design[:, 0]-project(stationary_unit, basis))))
    assert nesting_error < 1e-9
    assert coarse_error < 1e-8

    # A>=1/2, B<=3/2, f_GR>0 imply damping>=1/2 and positive Green kernel.
    # Uniform sup-norm bound for |h|<=beta, all masses and all retained times.
    duration = -np.log(0.25)
    norm_coefficient = 1.5*(2*duration-4*(1-np.exp(-duration/2)))
    rows = []
    sample = None
    for beta in AMPLITUDES:
        for scenario in ('stationary', 'weight_drift'):
            weights = beta*np.array([[0.8, 0.2],
                                     [0.8, 0.2] if scenario=='stationary' else [0.2, 0.8]])
            exact, tangent = evolve(weights, steps)
            exact_coarse, _ = evolve(weights, steps//2)
            integration_error = float(np.max(np.abs(exact-exact_coarse)))
            assert integration_error < 1e-9
            y, yt = project(exact, basis), project(tangent, basis)
            qs, _ = nnls_small(stat_design, y)
            qe, _ = nnls_small(epoch_design, y)
            qst, _ = nnls_small(stat_design, yt)
            qet, _ = nnls_small(epoch_design, yt)
            assert qs >= qe - 1e-10
            assert qet < 1e-15
            if scenario == 'stationary':
                assert qst < 1e-15
            r = beta*norm_coefficient
            assert r < 1
            power_bound = 2*r*r/(1-r) + (r/(1-r))**2
            epsilon = np.sqrt(ndata)*power_bound/SIGMA
            raw_error = float(np.max(np.abs(exact-tangent)))
            assert raw_error <= power_bound+1e-10
            actual_error = float(np.linalg.norm(y-yt))
            delta_bound = 2*epsilon*(np.sqrt(qs)+np.sqrt(qe))+2*epsilon**2
            assert abs((qs-qe)-(qst-qet)) <= delta_bound+1e-10
            # Distance to the unbounded tangent cone is a conservative lower
            # bound for the tangent class with total amplitude <= beta.
            margin = np.sqrt(qs)-epsilon
            rows.append(dict(scenario=scenario, amplitude=beta,
                             q_stationary=qs, q_epoch=qe, delta_q=qs-qe,
                             tangent_q_stationary=qst, tangent_q_epoch=qet,
                             actual_projected_error=actual_error,
                             uniform_projected_error_bound=epsilon,
                             stationary_exclusion_margin=margin,
                             raw_power_error=raw_error, raw_power_error_bound=power_bound,
                             operator_norm_bound=r, rk4_step_halving_error=integration_error))
            if beta == 0.003 and scenario == 'weight_drift':
                sample = (exact.reshape(8, 12), tangent.reshape(8, 12))
    summary = dict(
        status='illustrative finite-dictionary matter-power benchmark; not a survey forecast',
        observable='P_m/P_m_GR - 1 at identical initial power and fixed background',
        background=dict(Omega_m0=0.3, Omega_Lambda0=0.7),
        matching=dict(a_initial=0.25, D_initial=0.25, Dprime_initial=0.25,
                      tangent_initial=0.0, primordial_matching_claimed=False),
        k_h_Mpc=K.tolist(), masses_h_Mpc=MASSES.tolist(),
        output_scale_factors=np.exp(TIMES).tolist(), temporal_split_a=0.5,
        amplitudes=list(AMPLITUDES), fractional_power_sigma=SIGMA,
        covariance='sigma^2 I; synthetic', window='identity',
        nuisance='one scale-independent amplitude per output epoch',
        source='piecewise constant non-negative mass weights before growth integration',
        certification_domain='same two-mass bank, same temporal basis, total weight <= amplitude in each epoch; floating-point diagnostics, not interval arithmetic',
        integration_steps=steps, generator_step_halving_error=coarse_error,
        early_source_response_at_final_epoch=memory_amplitude,
        nesting_max_absolute_error=nesting_error, rows=rows)
    (ROOT/'causal_growth_transport_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    with (ROOT/'causal_growth_transport_benchmark.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.5), constrained_layout=True)
    for epoch in (1, 4, 7):
        line, = axes[0].plot(K, sample[0][epoch], label=f'a={np.exp(TIMES[epoch]):.2f}')
        axes[0].plot(K, sample[1][epoch], '--', color=line.get_color())
    axes[0].set(xscale='log', xlabel=r'$k\ [h\,\mathrm{Mpc}^{-1}]$',
                ylabel=r'$P_m/P_{m,\mathrm{GR}}-1$', title=r'Causal weight drift, amplitude $0.003$')
    axes[0].legend(fontsize=8)
    drift = [v for v in rows if v['scenario']=='weight_drift']
    axes[1].loglog(AMPLITUDES, [np.sqrt(v['q_stationary']) for v in drift], 'o-', label='Stationary-cone distance')
    axes[1].loglog(AMPLITUDES, [v['uniform_projected_error_bound'] for v in drift], 's--', label='Uniform linearization bound')
    axes[1].loglog(AMPLITUDES, [v['actual_projected_error'] for v in drift], '^:', label='Actual tangent error')
    axes[1].set(xlabel='Total response amplitude', ylabel='Whitened projected norm', title='Fixed contract; no detection significance')
    axes[1].legend(fontsize=7)
    (ROOT/'figs').mkdir(exist_ok=True)
    fig.savefig(ROOT/'figs/figS6_causal_growth_transport.png', dpi=220)
    plt.close(fig)
    print('Causal growth benchmark passed: nesting, KKT, stationary recovery, bounds, step halving.')
    for row in rows:
        if row['scenario']=='weight_drift':
            print('amplitude', row['amplitude'], 'margin', row['stationary_exclusion_margin'])


if __name__ == '__main__':
    run()
