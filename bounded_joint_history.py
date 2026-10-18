# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Vitantonio Castronuovo

"""Joint growth/velocity comparison with predeclared smooth history restrictions.

Fixed synthetic response model, not a survey likelihood or covariant parent.
The calculation reuses the reviewed growth/sensitivity equations, not their
old numerical outputs. All three predeclared two-mass targets are reported.
"""
from pathlib import Path
import json
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.special import hyp2f1

HERE=Path(__file__).resolve().parent
K=np.array([.035,.14]); MASS=.15; HEAVY=.30; BASE=.02
START=float(np.log(.25)); EPOCHS=np.array([np.log(.5),0.])
EDGES=np.r_[START,EPOCHS]
CENTERS=np.concatenate([a+(b-a)*np.arange(1,5)/5 for a,b in zip(EDGES[:-1],EDGES[1:])])
WIDTHS=np.repeat(.09*np.diff(EDGES),4)
TARGET_FRACTIONS=(.0001,.001,.01)
COARSE=np.zeros((8,2));COARSE[:4,0]=1;COARSE[4:,1]=1
NORM_SCALE=.01


def bumps(n):
    x=(n-CENTERS)/WIDTHS
    out=np.zeros_like(x);inside=np.abs(x)<1
    out[inside]=np.exp(1-1/(1-x[inside]**2))
    return out


def background(n):
    om=.3*np.exp(-3*n)/(.3*np.exp(-3*n)+.7)
    return 2-1.5*om,1.5*om


def initial_growth_rate():
    a=np.exp(START);om=.3/(.3+.7*a**3)
    g=hyp2f1(1/3,1,11/6,-.7/.3*a**3)
    # Hamilton astro-ph/0006089 equation(4), source-reused directly.
    return -1-om/2+(1-om)+5*om/(2*g)


def integrate(theta=None,target_fraction=None,method="DOP853",rtol=2e-12,atol=2e-14,max_step=.01):
    theta=np.zeros(8) if theta is None else np.asarray(theta)
    y=np.zeros(2+4+32);y[0]=1;y[1]=initial_growth_rate();y[2:4]=1
    values=[];jac=[];states=[]
    def rhs(n,y):
        aa,bb=background(n);e=y[2:4];ep=y[4:6]
        sensitivity=y[6:].reshape(2,2,8)
        psi=bumps(n);q=K*K/(K*K+np.exp(2*n)*MASS*MASS)
        response=(BASE+psi@theta)*q
        if target_fraction is not None:
            qh=K*K/(K*K+np.exp(2*n)*HEAVY*HEAVY)
            response=BASE*((1-target_fraction)*q+target_fraction*qh)
        cc=2*y[1]/y[0]+aa
        ds=np.stack((sensitivity[1],bb*(response[:,None]*sensitivity[0]+q[:,None]*psi[None,:]*e[:,None])-cc*sensitivity[1]))
        return np.r_[y[1],bb*y[0]-aa*y[1],ep,bb*response*e-cc*ep,ds.ravel()]
    knots=np.unique(np.r_[START,CENTERS-WIDTHS,CENTERS,CENTERS+WIDTHS,EPOCHS])
    for left,right in zip(knots[:-1],knots[1:]):
        sol=solve_ivp(rhs,(left,right),y,method=method,rtol=rtol,atol=atol,max_step=max_step)
        if not sol.success:raise RuntimeError(sol.message)
        y=sol.y[:,-1]
        if np.any(np.isclose(right,EPOCHS,atol=1e-14,rtol=0)):
            e=y[2:4];ep=y[4:6];ss=y[6:].reshape(2,2,8)
            values.append(np.r_[e*e-1,ep/e])
            jac.append(np.vstack((2*e[:,None]*ss[0],(ss[1]*e[:,None]-ep[:,None]*ss[0])/e[:,None]**2)))
            states.append(dict(E=e.tolist(),Eprime=ep.tolist(),reference_growth=y[0],reference_growth_rate=y[1]/y[0]))
    return np.concatenate(values),np.vstack(jac),states


def restrictions(theta):
    # Disjoint support and0<=psi<=1 give global amplitude bounds.
    lo=BASE+min(0.,float(theta.min()));hi=BASE+max(0.,float(theta.max()))
    # For the chosen bump |psi_x|<=8/e<3. This is an analytic conservative
    # bound; evaluate one active bump at a time, not a finite-grid claim.
    log_bound=float(np.max(3*np.abs(theta)/(WIDTHS*(BASE+np.minimum(0,theta)))))
    return dict(amplitude_bounds=[lo,hi],conservative_log_slope_bound=log_bound,
                amplitude_prior_pass=bool(lo>=.015-1e-14 and hi<=.025+1e-14),
                slope_prior_pass=bool(log_bound<=2+1e-12),
                admissible=bool(lo>=.015-1e-14 and hi<=.025+1e-14 and log_bound<=2+1e-12))


def coefficient_bounds(mapping):
    # Sufficient analytic slope box: 3|theta_i| <= 2 w_i (b0+min(theta_i,0)).
    lower=-2*WIDTHS/(3+2*WIDTHS)
    upper=2*WIDTHS/3
    return (np.array([max(lower[mapping[:,j]!=0]) for j in range(mapping.shape[1])]),
            np.array([min(upper[mapping[:,j]!=0]) for j in range(mapping.shape[1])]))


def fit(target,mapping):
    cached={}
    def evaluate(x):
        key=tuple(x)
        if key not in cached:cached[key]=integrate(BASE*mapping@x)
        return cached[key]
    def fun(x):return (evaluate(x)[0]-target)/NORM_SCALE
    def jac(x):return (evaluate(x)[1]@mapping)*BASE/NORM_SCALE
    answer=least_squares(fun,np.zeros(mapping.shape[1]),jac=jac,bounds=coefficient_bounds(mapping),
                         xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60)
    theta=BASE*mapping@answer.x;prediction,jacobian,states=integrate(theta)
    residual=(prediction-target).reshape(2,4)
    out=dict(coefficients=theta.tolist(),parameter_count=mapping.shape[1],
             max_power_residual=float(np.max(np.abs(residual[:,:2]))),
             max_growthrate_residual=float(np.max(np.abs(residual[:,2:]))),
             max_joint_residual=float(np.max(np.abs(residual))),
             evaluations=int(answer.nfev),optimizer_success=bool(answer.success),
             prediction=prediction.tolist(),states=states,**restrictions(theta))
    return out


def main():
    base,jac,_=integrate();cases=[]
    for fraction in TARGET_FRACTIONS:
        target,_,target_states=integrate(target_fraction=fraction)
        fine=fit(target,np.eye(8));coarse=fit(target,COARSE)
        if fine['admissible']:
            test,_,_=integrate(fine['coefficients'],method="Radau",rtol=4e-13,atol=4e-15,max_step=.005)
            exact_target,_,_=integrate(target_fraction=fraction,method="Radau",rtol=4e-13,atol=4e-15,max_step=.005)
            fine['independent_Radau_max_joint_residual']=float(np.max(np.abs(test-exact_target)))
        cases.append(dict(heavy_fraction=fraction,target=target.tolist(),target_states=target_states,H8=fine,H2=coarse))
    result=dict(status="numerical joint-snapshot comparison under declared temporal priors; not a survey or complete two-time parent inference",contract=dict(Omega_m0=.3,Omega_Lambda0=.7,k_h_Mpc=K.tolist(),masses_h_Mpc=[MASS,HEAVY],baseline=BASE,start_a=.25,output_a=[.5,1],heavy_fractions=list(TARGET_FRACTIONS),amplitude_prior=[.015,.025],log_slope_prior=2,centers_N=CENTERS.tolist(),halfwidths_N=WIDTHS.tolist(),initial_f=initial_growth_rate(),reference="source pure growing matter+Lambda, normalized D0(start)=1",observable_order="per epoch: P(k1),P(k2),f(k1)-f0,f(k2)-f0",H2="four coefficients tied within each source interval",coefficient_restriction="3 abs(theta_i) <= 2 width_i (baseline+min(theta_i,0)); sufficient slope box, not all slope-bounded histories",computational_tolerance=1e-10),baseline_joint_jacobian=jac.tolist(),singular_values=np.linalg.svd(jac,compute_uv=False).tolist(),condition=float(np.linalg.cond(jac)),cases=cases,numpy=np.__version__,scipy=scipy.__version__,limitations=["phenomenological temporal priors, not derived parent dynamics","ideal finite mode/time samples with identity projection","attained optimizer residuals do not certify a global minimum or exclusion of all admissible histories","floating point; no uniform physical or statistical radius","prior power-only counterfit not rerun"])
    (HERE/'bounded_joint_history_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'condition':result['condition'],'cases':[dict(heavy_fraction=x['heavy_fraction'],H8={k:v for k,v in x['H8'].items() if k in ('max_joint_residual','amplitude_bounds','conservative_log_slope_bound','admissible','independent_Radau_max_joint_residual')},H2={k:v for k,v in x['H2'].items() if k in ('max_joint_residual','amplitude_bounds','conservative_log_slope_bound','admissible')}) for x in cases]},indent=2),flush=True)


if __name__=='__main__':main()
