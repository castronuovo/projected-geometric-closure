"""Differential audit only. Does not modify the release or manuscript outputs.

Run with NumPy, SciPy and Matplotlib installed. Floating-point diagnostics,
not validated interval arithmetic. The continuum inequality is derived in
the accompanying report; no statistical significance is computed.
"""
import importlib.util
import json
from pathlib import Path
import numpy as np
from scipy.optimize import nnls
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
SOURCE = HERE / 'causal_growth_transport.py'
spec = importlib.util.spec_from_file_location('original_growth', SOURCE)
original = importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)


def declared_mass_grid(points):
    """Build the declared grid with platform-independent anchor insertion."""
    masses = np.geomspace(.0001, 100., points)
    for anchor in (.1, .4):
        nearest = int(np.argmin(np.abs(np.log(masses/anchor))))
        if np.isclose(masses[nearest], anchor, rtol=1e-13, atol=0.):
            masses[nearest] = anchor
        else:
            masses = np.sort(np.r_[masses, anchor])
    return masses


def bank(masses, steps=2048):
    """Vectorized independent RK4 tangent propagation for stationary masses."""
    k = original.K[:, None]
    masses = masses[None, :]
    d, dp = .25, .25
    u = np.zeros((12, masses.size))
    up = u.copy()
    n0 = np.log(.25)
    dt = -n0 / steps
    output = []

    def rhs(n, state):
        d, dp, u, up = state
        a, b = original.background(n)
        q = k*k/(k*k+np.exp(2*n)*masses*masses)
        return dp, b*d-a*dp, up, b*q-(2*dp/d+a)*up

    def add(state, deriv, h):
        return tuple(x+h*v for x, v in zip(state, deriv))

    state = (d, dp, u, up)
    for j in range(steps):
        n = n0+j*dt
        v1 = rhs(n, state)
        v2 = rhs(n+dt/2, add(state, v1, dt/2))
        v3 = rhs(n+dt/2, add(state, v2, dt/2))
        v4 = rhs(n+dt, add(state, v3, dt))
        state = tuple(x+dt*(a+2*b+2*c+d)/6
                      for x, a, b, c, d in zip(state, v1, v2, v3, v4))
        if (j+1) % (steps//8) == 0:
            output.append(2*state[2])
    return np.array(output).reshape(96, masses.size)


def smooth_target(alpha, tolerance):
    """C-infinity tanh drift of width 0.1 in ln(a), independently integrated."""
    initial = np.zeros((4, 12)); initial[0:2] = .25; initial[2] = 1.

    def rhs(n, vector):
        d, dp, e, ep = vector.reshape(4, 12)
        a, b = original.background(n)
        transition = .5*(1+np.tanh((n-np.log(.5))/.1))
        weights = alpha*np.array([.8-.6*transition, .2+.6*transition])
        q = original.K[:, None]**2/(original.K[:, None]**2+
                                   np.exp(2*n)*original.MASSES[None, :]**2)
        h = q@weights
        return np.array([dp, b*d-a*dp, ep, b*h*e-(2*dp/d+a)*ep]).ravel()

    solution = solve_ivp(rhs, (np.log(.25), 0.), initial.ravel(),
                         t_eval=original.TIMES, method='DOP853',
                         rtol=tolerance, atol=tolerance*.01)
    assert solution.success
    return (solution.y.T.reshape(8, 4, 12)[:, 2]**2-1).ravel()


def run():
    # Deliberately widen support before seeing the fit; retain all scan points.
    masses = declared_mass_grid(8193)
    raw = bank(masses)
    coarse = bank(masses, 1024)
    nuisance = np.kron(np.eye(8), np.ones((12, 1)))
    basis = np.linalg.qr(nuisance, mode='reduced')[0]
    design = original.project(raw, basis)
    original_error = 0.
    for m in (.1, .4):
        w = np.zeros((2, 2)); w[:, 0 if m == .1 else 1] = 1.
        _, tangent = original.evolve(w)
        original_error = max(original_error, float(np.max(np.abs(
            raw[:, np.flatnonzero(masses == m)[0]]-tangent))))
    # |dq/d log(m)| <= 1/2. Propagation supplies a per-output bound C(T).
    durations = original.TIMES-np.log(.25)
    coefficients = 1.5*(2*durations-4*(1-np.exp(-durations/2)))
    lipschitz = float(np.linalg.norm(np.repeat(coefficients, 12))/original.SIGMA)
    radius = float(np.max(np.diff(np.log(masses)))/2)
    # Q removes q=1 exactly (one amplitude per epoch). The two unresolved
    # tails therefore approach zero after projection, with explicit bounds.
    low_tail = 2*lipschitz*(masses[0]/original.K.min())**2
    high_tail = 2*lipschitz*(original.K.max()/(.25*masses[-1]))**2
    rows = []
    for scenario, alpha in [(s, a) for s in ('step', 'smooth') for a in original.AMPLITUDES]:
        if scenario == 'step':
            exact, tangent = original.evolve(alpha*np.array([[.8, .2], [.2, .8]]))
            exact_coarse, _ = original.evolve(alpha*np.array([[.8, .2], [.2, .8]]), 1024)
        else:
            exact = smooth_target(alpha, 1e-13)
            exact_coarse = smooth_target(alpha, 1e-11)
        y = original.project(exact, basis)
        norms = np.linalg.norm(design, axis=0)
        # Use a moderate subgrid to construct a witness; then validate it on
        # the full grid. A nonoptimal witness still gives a valid lower bound.
        indices = np.unique(np.r_[np.arange(0, len(masses), 16), len(masses)-1,
                                 np.flatnonzero(np.isin(masses, [.1, .4]))])
        indices = indices[norms[indices] > 1e-14]
        fit, _ = nnls(design[:, indices]/norms[indices], y, maxiter=10000)
        residual = y-(design[:, indices]/norms[indices])@fit
        distance = float(np.linalg.norm(residual))
        witness = residual/distance
        maximum = float(np.max(witness@design))
        r = alpha*coefficients[-1]
        epsilon = np.sqrt(96)*(2*r*r/(1-r)+(r/(1-r))**2)/original.SIGMA
        support_penalty = alpha*max(0., maximum+lipschitz*radius)
        margin = float(witness@y-support_penalty-epsilon)
        global_penalty = alpha*max(0., maximum+lipschitz*radius, low_tail, high_tail)
        rows.append(dict(scenario=scenario, amplitude=alpha, subgrid_fit_distance=distance,
                         target_integration_difference=float(np.max(np.abs(exact-exact_coarse))),
                         dual_target=float(witness@y), grid_dual_supremum=maximum,
                         bounded_continuum_support_penalty=support_penalty,
                         uniform_tangent_allowance=float(epsilon),
                         continuum_lower_margin=margin,
                         all_nonnegative_masses_lower_margin=float(witness@y-global_penalty-epsilon),
                         fitted_total_weight=float(np.sum(fit/norms[indices]))))
        print(rows[-1], flush=True)
    output = dict(status='Exploratory synthetic audit; floating point, no significance',
                  support_mass_h_Mpc=[.0001, 100.], grid_size=len(masses),
                  unresolved_low_mass_generator_bound=low_tail,
                  unresolved_high_mass_generator_bound=high_tail,
                  numpy_version=np.__version__,
                  log_mass_covering_radius=radius, projected_lipschitz=lipschitz,
                  max_generator_step_halving=float(np.max(np.abs(raw-coarse))),
                  agreement_with_original_generator=original_error, rows=rows)
    (HERE/'causal_transport_support_results.json').write_text(json.dumps(output, indent=2)+'\n')


if __name__ == '__main__':
    run()
