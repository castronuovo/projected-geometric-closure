"""Synthetic four-sample counterfit with a positive temporal amplitude.

Compare a two-mass target with a fixed-mass model on a fixed background.
This response-model example does not provide a survey constraint.
"""
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
K = np.array([0.035, 0.14])
MASS = 0.15
HEAVY = 0.30
BASE = 0.02
MIX = 0.001
START = float(np.log(0.25))
EPOCHS = np.array([np.log(0.5), 0.0])
EDGES = np.r_[START, EPOCHS]
CENTERS = np.concatenate([a + (b-a)*np.array([1/3, 2/3])
                          for a, b in zip(EDGES[:-1], EDGES[1:])])
WIDTHS = np.repeat(np.diff(EDGES)/8, 2)
NDATA = K.size * EPOCHS.size


def bumps(n):
    x = (n-CENTERS)/WIDTHS
    out = np.zeros_like(x)
    inside = np.abs(x) < 1
    out[inside] = np.exp(1-1/(1-x[inside]**2))
    return out


def background(n):
    om = 0.3*np.exp(-3*n)/(0.3*np.exp(-3*n)+0.7)
    return 2-1.5*om, 1.5*om


def drift(n):
    return MIX*(1+np.tanh((n-float(EPOCHS[0]))/0.1))/2


def integrate(theta=None, target=False, method="DOP853", rtol=2e-11,
              atol=2e-13, max_step=0.012):
    theta = np.zeros(NDATA) if theta is None else np.asarray(theta)
    state = np.zeros(2 + 2*K.size + 2*K.size*NDATA)
    state[:2] = 0.25
    state[2:2+K.size] = 1
    outputs, jacobians = [], []

    def rhs(n, y):
        a, b = background(n)
        e = y[2:2+K.size]
        ep = y[2+K.size:2+2*K.size]
        s = y[2+2*K.size:].reshape(2, K.size, NDATA)
        psi = bumps(n)
        q = K*K/(K*K+np.exp(2*n)*MASS*MASS)
        amp = BASE+psi@theta
        h = amp*q
        if target:
            qh = K*K/(K*K+np.exp(2*n)*HEAVY*HEAVY)
            w = drift(n)
            h = BASE*((1-w)*q+w*qh)
        c = 2*y[1]/y[0]+a
        ds = np.stack([s[1], b*(h[:, None]*s[0]
                               + q[:, None]*psi[None, :]*e[:, None])
                              - c*s[1]])
        return np.r_[y[1], b*y[0]-a*y[1], ep, b*h*e-c*ep, ds.ravel()]

    knots = np.unique(np.r_[START, CENTERS-WIDTHS, CENTERS,
                            CENTERS+WIDTHS, EPOCHS])
    for left, right in zip(knots[:-1], knots[1:]):
        sol = solve_ivp(rhs, (left, right), state, method=method,
                        rtol=rtol, atol=atol, max_step=max_step)
        if not sol.success:
            raise RuntimeError(sol.message)
        state = sol.y[:, -1]
        if np.any(np.isclose(right, EPOCHS, atol=1e-14, rtol=0)):
            e = state[2:2+K.size]
            outputs.append(e*e-1)
            sens = state[2+2*K.size:].reshape(2, K.size, NDATA)[0]
            jacobians.append(2*e[:, None]*sens)
    return np.concatenate(outputs), np.vstack(jacobians)


def main():
    target, _ = integrate(target=True)
    baseline, jac = integrate()
    theta = np.linalg.solve(jac, target-baseline)
    iterations = []
    for step in range(8):
        pred, current_jac = integrate(theta)
        residual = pred-target
        iterations.append(dict(step=step, max_power_residual=float(max(abs(residual)))))
        if max(abs(residual)) < 2e-13:
            break
        theta -= np.linalg.solve(current_jac, residual)
    # Independent integration method and tighter temporal/accuracy settings.
    check, _ = integrate(theta, method="Radau", rtol=8e-13,
                         atol=8e-15, max_step=0.006)
    target_check, _ = integrate(target=True, method="Radau", rtol=8e-13,
                                atol=8e-15, max_step=0.006)
    grid = np.linspace(START, 0, 4001)
    amplitudes = np.array([BASE+bumps(n)@theta for n in grid])
    # Supports are disjoint and each bump ranges over [0,1]: exact bounds.
    amp_lower = BASE+min(0., float(theta.min()))
    amp_upper = BASE+max(0., float(theta.max()))
    instant_gap = max(float(max(abs(
        (BASE+bumps(n)@theta)*K*K/(K*K+np.exp(2*n)*MASS*MASS)
        - BASE*((1-drift(n))*K*K/(K*K+np.exp(2*n)*MASS*MASS)
                + drift(n)*K*K/(K*K+np.exp(2*n)*HEAVY*HEAVY))))) for n in grid)
    fd = np.column_stack([(integrate(np.eye(NDATA)[i]*1e-5)[0]
                           - integrate(-np.eye(NDATA)[i]*1e-5)[0])/2e-5
                          for i in range(NDATA)])
    result = dict(status="numerically_verified synthetic finite-power counterfit; no survey or parent inference",
                  k_h_Mpc=K.tolist(), mass_h_Mpc=MASS, second_mass_h_Mpc=HEAVY,
                  baseline_amplitude=BASE, maximum_second_mass_fraction=MIX,
                  start_N=START, output_N=EPOCHS.tolist(), centers_N=CENTERS.tolist(),
                  widths_N=WIDTHS.tolist(), fitted_coefficients=theta.tolist(),
                  target_power_response=target.tolist(), baseline_power_response=baseline.tolist(),
                  fitted_power_response=pred.tolist(),
                  baseline_jacobian=jac.tolist(), jacobian_singular_values=np.linalg.svd(jac,compute_uv=False).tolist(),
                  jacobian_condition=float(np.linalg.cond(jac)),
                  relative_tangent_finite_difference_error=float(np.linalg.norm(fd-jac)/np.linalg.norm(jac)),
                  exact_amplitude_bounds=[amp_lower, amp_upper], sampled_amplitude_range=[float(amplitudes.min()),float(amplitudes.max())],
                  max_instantaneous_force_response_difference=instant_gap,
                  max_fitted_power_residual=float(max(abs(pred-target))),
                  independent_Radau_max_power_residual=float(max(abs(check-target_check))),
                  DOP853_vs_Radau_fitted_difference=float(max(abs(pred-check))),
                  DOP853_vs_Radau_target_difference=float(max(abs(target-target_check))),
                  iterations=iterations, numpy=np.__version__, scipy=scipy.__version__,
                  limitations=["four ideal matter-power samples", "fixed background and matched data", "phenomenological amplitude freedom; no covariant completion", "no uniform reachable radius or survey significance", "floating point, not interval arithmetic"])
    assert amp_lower > 0
    assert result["max_fitted_power_residual"] < 2e-12
    assert result["independent_Radau_max_power_residual"] < 2e-11
    assert result["relative_tangent_finite_difference_error"] < 1e-6
    (HERE/"causal_amplitude_reachability_results.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:result[k] for k in ("fitted_coefficients","exact_amplitude_bounds","jacobian_singular_values","jacobian_condition","max_fitted_power_residual","independent_Radau_max_power_residual","relative_tangent_finite_difference_error","max_instantaneous_force_response_difference")},indent=2))


if __name__ == "__main__":
    main()
