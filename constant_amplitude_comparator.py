# SPDX-License-Identifier: BSD-3-Clause
"""Constant-amplitude one-mass comparator for the frozen JH1 contract.

Postprocessing-level addition: reuses background, initial data, targets and
observable definitions of bounded_joint_history.py (imported,
not modified; its stored results are not overwritten). Fits the same two-mass
targets with a one-mass model with CONSTANT amplitude c (mass fixed to the
light mass) and, as a secondary comparator, constant c with free mass.
Floating point; no validated ODE enclosure.
"""
import sys, json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import bounded_joint_history as J

def observables(c,mass,fraction=None,method="DOP853",rtol=2e-12,atol=2e-14):
    y=np.array([1,J.initial_growth_rate(),1,1,0,0],float)
    def rhs(n,y):
        aa,bb=J.background(n);e=y[2:4];ep=y[4:6]
        q=J.K**2/(J.K**2+np.exp(2*n)*mass**2)
        if fraction is None: resp=c*q
        else:
            qh=J.K**2/(J.K**2+np.exp(2*n)*J.HEAVY**2)
            resp=J.BASE*((1-fraction)*q+fraction*qh)
        cc=2*y[1]/y[0]+aa
        return np.r_[y[1],bb*y[0]-aa*y[1],ep,bb*resp*e-cc*ep]
    out=[];n0=J.START
    for n1 in J.EPOCHS:
        s=solve_ivp(rhs,(n0,n1),y,method=method,rtol=rtol,atol=atol,max_step=.01)
        y=s.y[:,-1];n0=n1;e=y[2:4];ep=y[4:6]
        out.append(np.r_[e*e-1,ep/e])
    return np.concatenate(out)

def fit(target,free_mass,bounded):
    lo=[.015 if bounded else 0.]; hi=[.025 if bounded else 1.]
    x0=[J.BASE]
    if free_mass: lo.append(.01);hi.append(1.);x0.append(J.MASS)
    f=lambda x:(observables(x[0],x[1] if free_mass else J.MASS)-target)/J.NORM_SCALE
    r=least_squares(f,x0,bounds=(lo,hi),xtol=1e-14,ftol=1e-14,gtol=1e-14)
    res=observables(r.x[0],r.x[1] if free_mass else J.MASS)-target
    return dict(params=r.x.tolist(),max_joint_residual=float(np.abs(res).max()),
                max_power_residual=float(np.abs(res.reshape(2,4)[:,:2]).max()),
                max_growthrate_residual=float(np.abs(res.reshape(2,4)[:,2:]).max()),success=bool(r.success))

def main():
    ref=json.load(open(HERE/'bounded_joint_history_results.json'))
    cases=[]
    for case in ref['cases']:
        fr=case['heavy_fraction'];target=np.array(case['target'])
        d0=float(np.abs(observables(0,J.MASS,fr)-target).max());assert d0<1e-9,d0   # reproduces frozen target
        null=observables(J.BASE,J.MASS)-target
        c=dict(heavy_fraction=fr,target_reproduction_max_diff=d0,target_max_abs=float(np.abs(target).max()),
               unfit_light_mass_const_BASE_max_residual=float(np.abs(null).max()),
               const_amp_fixed_mass_unbounded=fit(target,False,False),
               const_amp_fixed_mass_prior_box=fit(target,False,True),
               const_amp_free_mass_unbounded=fit(target,True,False),
               H8_frozen=case['H8']['max_joint_residual'],H2_frozen=case['H2']['max_joint_residual'])
        c['ratio_const_fixed_over_H8']=c['const_amp_fixed_mass_unbounded']['max_joint_residual']/c['H8_frozen']
        cases.append(c)
    json.dump(dict(status='postprocessing comparator under frozen JH1 contract; not survey or parent inference',cases=cases),open(HERE/'constant_amplitude_comparator_results.json','w'),indent=2)
    for c in cases:
        print(c['heavy_fraction'],'target',round(c['target_max_abs'],5),'unfit',c['unfit_light_mass_const_BASE_max_residual'],'cfix',c['const_amp_fixed_mass_unbounded']['max_joint_residual'],c['const_amp_fixed_mass_unbounded']['params'],'cbox',c['const_amp_fixed_mass_prior_box']['max_joint_residual'],'cfree',c['const_amp_free_mass_unbounded']['max_joint_residual'],c['const_amp_free_mass_unbounded']['params'],'H8',c['H8_frozen'],'H2',c['H2_frozen'])
main()
