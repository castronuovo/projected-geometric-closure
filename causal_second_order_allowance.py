# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Vitantonio Castronuovo
"""Second-order linearization allowance for the causal cone test; not a forecast.

Extends causal_growth_transport.py and causal_transport_support_audit.py without
modifying them. Floating-point diagnostics on a declared mass grid, not
validated interval arithmetic. Synthetic fixed-background contract only.

Exact structure. At fixed background and for a non-negative measure nu with
h(k,N)=int q dnu in [0,alpha], E=D/D_GR=sum_n T^n 1 with T positive and linear
in nu. Hence, for R_P=E^2-1,

    R_P = 2 e1 + Q + T3,   e1=T1,   Q = 2 T^2 1 + (T1)^2,
    T3 = 2R3 + 2 e1 e2 + 2 e1 R3 + (e2+R3)^2 >= 0,

with e2=T^2 1, R3=sum_{n>=3}T^n 1. Q is bilinear in nu:
Q = int int S(m,m') dnu dnu, S = 2 T_m T_m' 1 + (T_m 1)(T_m' 1). The
projected second-order norm is bounded by alpha^2 n2, n2 = sup ||P S_sym||/sigma,
and T3 by the Neumann tail with r = alpha*rho0, rho0 = sup_N (T_{h=1} 1)(N).
"""

import csv
import importlib.util
import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import nnls

HERE = Path(__file__).resolve().parent


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


original = _load('original_growth', 'causal_growth_transport.py')
audit = _load('support_audit', 'causal_transport_support_audit.py')

K = original.K
SIGMA = original.SIGMA
NK, NT = len(K), 8
NDATA = NK * NT
STEPS = 512
LOG_GRID_POINTS = 241
MASS_RANGE = (1.0e-3, 1.0e2)
ALPHAS = (0.0003, 0.001, 0.003, 0.01, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1)
N0 = np.log(0.25)


def tail_bound_t3(r):
    """Sup-norm bound on T3 for operator norm r<1 (all terms non-negative)."""
    return (2*r**3/(1-r) + 2*r**3 + 2*r**4/(1-r)
            + (r**2 + r**3/(1-r))**2)


def old_allowance(alpha):
    r = alpha*1.5*(2*(-N0) - 4*(1-np.exp(N0/2)))
    return np.sqrt(NDATA)*(2*r*r/(1-r) + (r/(1-r))**2)/SIGMA


def project_epoch_mean(vectors):
    """Projection orthogonal to one scale-independent amplitude per epoch."""
    v = vectors.reshape(NT, NK, -1)
    return (v - v.mean(axis=1, keepdims=True)).reshape(NDATA, -1)/SIGMA


def pair_norms(steps):
    """Return n2 grid data, rho0 and the validation against stored errors."""
    dt = -N0/steps
    hh = dt/2
    grid = np.arange(2*steps + 1)
    times = N0 + hh*grid

    def gr_rhs(n, y):
        a, b = original.background(n)
        return np.array([y[1], b*y[0] - a*y[1]])

    y = np.array([0.25, 0.25])
    gr = [y.copy()]
    for i in range(2*steps):
        n = N0 + i*hh
        k1 = gr_rhs(n, y)
        k2 = gr_rhs(n + hh/2, y + hh*k1/2)
        k3 = gr_rhs(n + hh/2, y + hh*k2/2)
        k4 = gr_rhs(n + hh, y + hh*k3)
        y = y + hh*(k1 + 2*k2 + 2*k3 + k4)/6
        gr.append(y.copy())
    gr = np.array(gr)
    a_coef, b_coef = original.background(times)
    damping = 2*gr[:, 1]/gr[:, 0] + a_coef

    logs = np.geomspace(MASS_RANGE[0], MASS_RANGE[1], LOG_GRID_POINTS)
    for anchor in (0.1, 0.4):
        if not np.any(np.isclose(logs, anchor, rtol=1e-12)):
            logs = np.sort(np.r_[logs, anchor])
    masses = np.r_[0.0, logs]   # m=0 gives q=1: exactly nuisance-degenerate
    nm = len(masses)

    def q_of(i):
        return K[None, :]**2/(K[None, :]**2 + np.exp(2*times[i])*masses[:, None]**2)

    def rk4(state, rhs, i0):
        k1 = rhs(i0, state)
        k2 = rhs(i0 + 1, state + dt*k1/2)
        k3 = rhs(i0 + 1, state + dt*k2/2)
        k4 = rhs(i0 + 2, state + dt*k3)
        return state + dt*(k1 + 2*k2 + 2*k3 + k4)/6

    def rhs1(i, s):
        return np.array([s[1], -damping[i]*s[1] + b_coef[i]*q_of(i)])

    s = np.zeros((2, nm, NK))
    u_all, v_all, e1_out = [s[0].copy()], [s[1].copy()], []
    for j in range(steps):
        s = rk4(s, rhs1, 2*j)
        u_all.append(s[0].copy())
        v_all.append(s[1].copy())
        if (j + 1) % (steps//8) == 0:
            e1_out.append(s[0].copy())
    u_all, v_all, e1_out = np.array(u_all), np.array(v_all), np.array(e1_out)

    def e1_at(i):
        if i % 2 == 0:
            return u_all[i//2]
        j = i//2   # cubic Hermite midpoint
        return (u_all[j] + u_all[j+1])/2 + dt/8*(v_all[j] - v_all[j+1])

    def rhs2(i, s):
        source = b_coef[i]*q_of(i)[:, None, :]*e1_at(i)[None, :, :]
        return np.array([s[1], -damping[i]*s[1] + source])

    s = np.zeros((2, nm, nm, NK))
    c_out = []
    for j in range(steps):
        s = rk4(s, rhs2, 2*j)
        if (j + 1) % (steps//8) == 0:
            c_out.append(s[0].copy())
    c_out = np.array(c_out)                       # (epoch, m, m', k)

    s_full = 2*c_out + e1_out[:, :, None, :]*e1_out[:, None, :, :]
    s_sym = (s_full + s_full.transpose(0, 2, 1, 3))/2
    vec = s_sym.transpose(0, 3, 1, 2).reshape(NDATA, nm*nm)
    norms = np.linalg.norm(project_epoch_mean(vec), axis=0).reshape(nm, nm)
    index = np.unravel_index(int(norms.argmax()), norms.shape)
    rho0 = float(u_all[:, 0, 0].max())

    # Validation: Q of the stored stationary two-mass target vs stored error.
    stored = {}
    with (HERE/'causal_growth_transport_benchmark.csv').open(newline='') as f:
        for row in csv.DictReader(f):
            if row['scenario'] == 'stationary':
                stored[float(row['amplitude'])] = float(row['actual_projected_error'])
    i1 = int(np.argmin(abs(masses - 0.1)))
    i2 = int(np.argmin(abs(masses - 0.4)))
    cube = vec.reshape(NDATA, nm, nm)
    validation = []
    for alpha, actual in sorted(stored.items()):
        w = np.zeros(nm)
        w[i1], w[i2] = 0.8*alpha, 0.2*alpha
        q_vector = (cube*w[None, :, None]*w[None, None, :]).sum(axis=(1, 2))
        predicted = float(np.linalg.norm(project_epoch_mean(q_vector[:, None])))
        validation.append(dict(amplitude=alpha, second_order=predicted,
                               stored_actual_error=actual,
                               relative_difference=(actual - predicted)/actual))
    return dict(n2_grid=float(norms.max()),
                n2_argmax_masses=[float(masses[index[0]]), float(masses[index[1]])],
                rho0=rho0, validation=validation, n_masses=nm)


def certified_range(n2_cover, rho0):
    masses = audit.declared_mass_grid(8193)
    raw = audit.bank(masses)
    basis = np.linalg.qr(np.kron(np.eye(NT), np.ones((NK, 1))), mode='reduced')[0]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)   # Accelerate matmul noise
        design = original.project(raw, basis)
    assert np.all(np.isfinite(design))
    norms = np.linalg.norm(design, axis=0)
    durations = original.TIMES - np.log(0.25)
    coefficients = 1.5*(2*durations - 4*(1 - np.exp(-durations/2)))
    lipschitz = float(np.linalg.norm(np.repeat(coefficients, NK))/SIGMA)
    radius = float(np.max(np.diff(np.log(masses)))/2)
    rows = []
    for alpha in ALPHAS:
        exact, _ = original.evolve(alpha*np.array([[0.8, 0.2], [0.2, 0.8]]))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            y = original.project(exact, basis)
        indices = np.unique(np.r_[np.arange(0, len(masses), 16), len(masses) - 1,
                                  np.flatnonzero(np.isin(masses, [0.1, 0.4]))])
        indices = indices[norms[indices] > 1e-14]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            fit, _ = nnls(design[:, indices]/norms[indices], y, maxiter=10000)
            residual = y - (design[:, indices]/norms[indices])@fit
        distance = float(np.linalg.norm(residual))
        witness = residual/distance
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            maximum = float(np.max(witness@design))
            target = float(witness@y)
        penalty = alpha*max(0.0, maximum + lipschitz*radius)
        eps_old = old_allowance(alpha)
        r = alpha*rho0
        eps_second = alpha**2*n2_cover
        eps_third = np.sqrt(NDATA)*tail_bound_t3(r)/SIGMA
        rows.append(dict(
            amplitude=alpha, dual_target=target, support_penalty=penalty,
            allowance_previous=eps_old, allowance_second_order=eps_second,
            allowance_third_order=eps_third,
            allowance_new=eps_second + eps_third,
            margin_previous=target - penalty - eps_old,
            margin_new=target - penalty - eps_second - eps_third))
    return rows, lipschitz, radius


def crossing(rows, key):
    """Interpolate sampled margins; this does not certify a root or interval."""
    for a, b in zip(rows[:-1], rows[1:]):
        if a[key] > 0 >= b[key]:
            fraction = a[key]/(a[key] - b[key])
            return a['amplitude'] + fraction*(b['amplitude'] - a['amplitude'])
    return None


def plot_results(rows, validation):
    """Render stored diagnostics without repeating the scientific calculation."""
    alphas = np.array([r['amplitude'] for r in rows])
    signal = np.array([r['dual_target'] - r['support_penalty'] for r in rows])
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), constrained_layout=True)
    axes[0].loglog(alphas, [r['allowance_previous'] for r in rows], 's--',
                   label='Previous uniform allowance')
    axes[0].loglog(alphas, [r['allowance_new'] for r in rows], 'o-',
                   label='Second-order allowance')
    actual = [(v['amplitude'], v['stored_actual_error'])
              for v in validation]
    axes[0].loglog(*zip(*actual), '^:', label='Actual stationary error')
    axes[0].loglog(alphas, signal, 'k-', lw=1, label='Dual signal minus support penalty')
    axes[0].set(xlabel='Total response amplitude', ylabel='Whitened projected norm',
                title='Linearization allowance')
    axes[0].legend(fontsize=7)
    for key, label, style in (('margin_previous', 'Previous', 's--'),
                              ('margin_new', 'Second-order', 'o-')):
        axes[1].semilogx(alphas, [r[key]/s for r, s in zip(rows, signal)], style,
                         label=label)
    axes[1].axhline(0.0, color='k', lw=0.8)
    axes[1].set(xlabel='Total response amplitude', ylim=(-1.0, 1.0),
                ylabel='Estimated margin / (dual signal minus penalty)',
                title='Estimated exclusion margin (clipped at $-1$)')
    axes[1].legend(fontsize=7)
    (HERE/'figs').mkdir(exist_ok=True)
    fig.savefig(HERE/'figs/figS8_second_order_allowance.png', dpi=220)
    plt.close(fig)


def run():
    main = pair_norms(STEPS)
    half = pair_norms(STEPS//2)
    logs = np.geomspace(MASS_RANGE[0], MASS_RANGE[1], LOG_GRID_POINTS)
    spacing = float(np.log(logs[1]/logs[0]))
    dq = 0.5*spacing/2                       # |dq/dln m| <= 1/2 at half spacing
    lip = 3*main['rho0']**2*np.sqrt(NDATA)/SIGMA
    low_tail = (MASS_RANGE[0]/K.min())**2
    high_tail = (K.max()/(0.25*MASS_RANGE[1]))**2
    n2_cover = main['n2_grid'] + lip*(2*dq) + lip*(low_tail + high_tail)
    rows, lipschitz, radius = certified_range(n2_cover, main['rho0'])
    # Grid-level variant (no Lipschitz covering of the mass continuum).
    for row in rows:
        r = row['amplitude']*main['rho0']
        row['allowance_new_grid_level'] = (
            row['amplitude']**2*main['n2_grid']
            + np.sqrt(NDATA)*tail_bound_t3(r)/SIGMA)
        row['margin_new_grid_level'] = (row['dual_target'] - row['support_penalty']
                                        - row['allowance_new_grid_level'])
    stored_bound_holds = all(
        v['stored_actual_error'] <= v['amplitude']**2*n2_cover
        + np.sqrt(NDATA)*tail_bound_t3(v['amplitude']*main['rho0'])/SIGMA
        for v in main['validation'])
    summary = dict(
        status='audit extension; floating point on a declared mass grid; synthetic contract; no significance',
        steps=STEPS, n_masses_pair_grid=main['n_masses'],
        n2_grid=main['n2_grid'], n2_argmax_masses_h_Mpc=main['n2_argmax_masses'],
        n2_step_halving_difference=abs(main['n2_grid'] - half['n2_grid']),
        n2_with_lipschitz_covering=n2_cover,
        lipschitz_constant_per_unit_q=lip,
        low_mass_tail_q_shift=low_tail, high_mass_tail_q_shift=high_tail,
        rho0_unit_operator_norm=main['rho0'],
        previous_operator_norm_coefficient=1.5*(2*(-N0) - 4*(1 - np.exp(N0/2))),
        validation_second_order_vs_stored_stationary_error=main['validation'],
        new_allowance_bounds_all_stored_stationary_errors=bool(stored_bound_holds),
        first_order_support_penalty_lipschitz=lipschitz,
        first_order_covering_radius=radius,
        alpha_star_previous=crossing(rows, 'margin_previous'),
        alpha_star_new=crossing(rows, 'margin_new'),
        alpha_star_new_grid_level=crossing(rows, 'margin_new_grid_level'),
        numpy_version=np.__version__)
    assert stored_bound_holds
    (HERE/'causal_second_order_allowance_results.json').write_text(
        json.dumps(summary, indent=2) + '\n')
    fields = list(rows[0])
    with (HERE/'causal_second_order_certified_range.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    plot_results(rows, main['validation'])
    print('n2 grid %.4f; with covering %.4f; rho0 %.5f' % (
        main['n2_grid'], n2_cover, main['rho0']))
    print('alpha* previous %s, new %s, new (grid level) %s' % (
        summary['alpha_star_previous'], summary['alpha_star_new'],
        summary['alpha_star_new_grid_level']))


if __name__ == '__main__':
    run()
