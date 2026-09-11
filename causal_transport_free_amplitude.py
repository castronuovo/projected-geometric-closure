"""Small adversarial test: stationary spectral shape with free epoch amplitudes.

The competitor has weights w[e,m] = A[e] p[m], with A[e] >= 0,
p[m] >= 0 and sum(p)=1. Thus its relative spectrum is stationary while its
total amplitude can change independently in the two source intervals.
This is a nonconvex rank-one cone; alternating NNLS is run from many starts.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar, nnls

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "causal_growth_transport.py"
spec = importlib.util.spec_from_file_location("original_growth", SOURCE)
original = importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)


def epoch_bank(masses, steps=2048):
    """Return projected tangent generators D[e][:,m] before whitening."""
    k = original.K[:, None]
    masses = masses[None, :]
    n0 = np.log(.25)
    dt = -n0 / steps
    outputs = [[], []]

    for source_epoch in (0, 1):
        d, dp = .25, .25
        u = np.zeros((12, masses.size)); up = u.copy()

        def rhs(n, state, active):
            d0, dp0, u0, up0 = state
            a, b = original.background(n)
            q = k*k/(k*k+np.exp(2*n)*masses*masses) if active else 0.
            return dp0, b*d0-a*dp0, up0, b*q-(2*dp0/d0+a)*up0

        def add(state, deriv, h):
            return tuple(x+h*v for x, v in zip(state, deriv))

        state = (d, dp, u, up)
        for j in range(steps):
            n = n0+j*dt
            active = (j < steps//2) == (source_epoch == 0)
            v1 = rhs(n, state, active)
            v2 = rhs(n+dt/2, add(state, v1, dt/2), active)
            v3 = rhs(n+dt/2, add(state, v2, dt/2), active)
            v4 = rhs(n+dt, add(state, v3, dt), active)
            state = tuple(x+dt*(a+2*b+2*c+d1)/6
                          for x, a, b, c, d1 in zip(state, v1, v2, v3, v4))
            if (j+1) % (steps//8) == 0:
                outputs[source_epoch].append(2*state[2])
    return [np.array(x).reshape(96, masses.size) for x in outputs]


def rank_one_fit(designs, target, seed):
    """Alternating NNLS; returns a stationary shape and two amplitudes."""
    rng = np.random.default_rng(seed)
    amplitudes = 10**rng.uniform(-4.5, -2.5, 2)
    previous = np.inf
    for _ in range(1000):
        combined = amplitudes[0]*designs[0]+amplitudes[1]*designs[1]
        shape, _ = nnls(combined, target, maxiter=20000)
        total = shape.sum()
        if total == 0:
            return np.inf, amplitudes, shape
        shape /= total
        amplitudes *= total
        epoch_vectors = np.column_stack([d@shape for d in designs])
        amplitudes, residual = nnls(epoch_vectors, target)
        if abs(previous-residual) <= 1e-13*max(1., previous):
            break
        previous = residual
    return residual, amplitudes, shape


def ratio_profile_fit(designs, target):
    """Profile the only nonlinear parameter lambda=A1/A0 on log scale."""
    cache = {}

    def objective(log_ratio):
        key = float(log_ratio)
        if key not in cache:
            ratio = np.exp(key)
            weights, residual = nnls(designs[0]+ratio*designs[1], target,
                                     maxiter=30000)
            cache[key] = (float(residual), weights)
        return cache[key][0]

    grid = np.linspace(-12., 12., 121)
    values = np.array([objective(x) for x in grid])
    candidates = []
    for i in range(1, len(grid)-1):
        if values[i] <= values[i-1] and values[i] <= values[i+1]:
            candidates.append(minimize_scalar(objective,
                              bounds=(grid[i-1], grid[i+1]), method="bounded",
                              options={"xatol": 1e-12}))
    best = min(candidates, key=lambda x: x.fun)
    residual, weights = cache[min(cache, key=lambda x: abs(x-best.x))]
    ratio = float(np.exp(best.x))
    total = float(weights.sum())
    shape = weights/total
    amplitudes = np.array([total, ratio*total])
    return residual, amplitudes, shape, float(values.min()), int(len(candidates))


def run():
    # Wide, predeclared grid plus the injected masses. End tails are nulled by
    # the per-output scale-independent nuisance and checked by enlargement.
    masses = np.unique(np.r_[np.geomspace(.0001, 100., 1025), .1, .4])
    nuisance = np.kron(np.eye(8), np.ones((12, 1)))
    basis = np.linalg.qr(nuisance, mode="reduced")[0]
    designs = [original.project(x, basis) for x in epoch_bank(masses)]
    rows = []
    for alpha in original.AMPLITUDES:
        exact, tangent = original.evolve(alpha*np.array([[.8, .2], [.2, .8]]))
        y_exact = original.project(exact, basis)
        y_tangent = original.project(tangent, basis)
        fits = [rank_one_fit(designs, y_tangent, seed) for seed in range(32)]
        local_residual, _, _ = min(fits, key=lambda x: x[0])
        residual, amplitudes, shape, scan_minimum, minima_count = ratio_profile_fit(
            designs, y_tangent)
        candidate = sum(a*(d@shape) for a, d in zip(amplitudes, designs))
        exact_residual = float(np.linalg.norm(y_exact-candidate))
        duration = -np.log(.25)
        norm_coefficient = 1.5*(2*duration-4*(1-np.exp(-duration/2)))
        response_cap = float(np.max(amplitudes))
        operator_bound = response_cap*norm_coefficient
        power_bound = (2*operator_bound**2/(1-operator_bound)
                       +(operator_bound/(1-operator_bound))**2)
        uniform_allowance = float(np.sqrt(96)*power_bound/original.SIGMA)
        unrestricted = np.column_stack(designs)
        _, epoch_residual = nnls(unrestricted, y_tangent, maxiter=30000)
        active = np.flatnonzero(shape > 1e-8*shape.max())
        rows.append(dict(
            amplitude=alpha,
            stationary_shape_free_amplitude_tangent_residual=float(residual),
            residual_against_exact_target=exact_residual,
            fitted_response_cap=response_cap,
            uniform_error_allowance_at_fitted_cap=uniform_allowance,
            exclusion_certified=bool(residual > uniform_allowance),
            unrestricted_epoch_cone_tangent_residual=float(epoch_residual),
            fitted_epoch_amplitudes=amplitudes.tolist(),
            active_mass_count=int(active.size),
            active_mass_range_h_Mpc=[float(masses[active].min()),
                                     float(masses[active].max())] if active.size else None,
            best_residual_spread=[float(min(x[0] for x in fits)),
                                  float(max(x[0] for x in fits))],
            profiled_amplitude_ratio=float(amplitudes[1]/amplitudes[0]),
            coarse_log_ratio_scan_minimum=scan_minimum,
            number_of_profile_minima=minima_count,
            agreement_with_best_alternating_fit=float(residual-local_residual)))
        print(rows[-1], flush=True)
    result = dict(
        status="Exploratory finite-grid adversarial test; no significance or global-optimum certificate",
        competitor="one nonnegative spectral shape shared across two source intervals, with two independent nonnegative total amplitudes",
        grid_mass_h_Mpc=[float(masses[0]), float(masses[-1])],
        grid_size=int(masses.size), starts=32,
        numpy=np.__version__, rows=rows)
    (HERE/"causal_transport_free_amplitude_results.json").write_text(
        json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    run()
