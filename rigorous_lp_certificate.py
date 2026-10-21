# SPDX-License-Identifier: BSD-3-Clause
"""Analytic lower-bound construction evaluated in floating point (not a validated enclosure).

Class C = {b : |b-0.02| <= eta=0.005, |b'| <= Lb=0.05}  (contained in the JH1 amplitude
prior [0.015,0.025]; the slope bound is implied by |dlnb/dN|<=2 only as an upper
envelope, so every history of the original prior with that slope bound lies in C).
For a competitor of mass m we linearize about the constant reference b_ref=0.02 at
the SAME mass: E = sum_j S^j E_ref with S the (nonnegative-kernel) integral operator of
the perturbation dg=db*q, so
   O(b,m) = O(b_ref,m) + J db + Rem,   |Rem| <= explicit bound from Neumann series.
Lower bound:  max|O(b,m)-O_t| >= t_LP(box,Lipschitz) - max_j Rem_j - interpolation allowance.
The constants (rho, rho', E_ref, J) are evaluated by floating-point ODE integration,
not by validated enclosures. Target O_t: constant b=0.02 at m0=0.15 (frozen JH1 contract).
"""
import sys, json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linprog
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import bounded_joint_history as J
N0=J.START;ETA=.005;LB=.05;BREF=.02
M=int(sys.argv[1]) if __name__ == '__main__' and len(sys.argv)>1 else 400
TAG=sys.argv[2] if __name__ == '__main__' and len(sys.argv)>2 else ''
M0=float(sys.argv[3]) if __name__ == '__main__' and len(sys.argv)>3 else .15
FRAC=float(sys.argv[4]) if __name__ == '__main__' and len(sys.argv)>4 else 0.
RLIST=[float(x) for x in sys.argv[5].split(',')] if __name__ == '__main__' and len(sys.argv)>5 else [.9,.95,.98,1.02,1.05,1.1,1.2,1.5]
nodes=np.linspace(N0,0,M+1);dN=nodes[1]-nodes[0]
ep_idx=[int(np.argmin(np.abs(nodes-e))) for e in J.EPOCHS]
assert all(abs(nodes[i]-e)<1e-12 for i,e in zip(ep_idx,J.EPOCHS))
def reference_and_jacobian(mass):
    K=J.K;nm=M+1
    # states: D0,D0',  E[2],E'[2],  rho[2],rho'[2],  dE[nm,2],dE'[nm,2]
    y0=np.zeros(2+4+4+4*nm);y0[0]=1;y0[1]=J.initial_growth_rate();y0[2:4]=1
    def rhs(n,y):
        aa,bb=J.background(n)
        D0,D0p=y[0],y[1];E=y[2:4];Ep=y[4:6];rho=y[6:8];rhop=y[8:10]
        q=K*K/(K*K+np.exp(2*n)*mass*mass);g=BREF*q;cc=2*D0p/D0+aa
        dE=y[10:10+2*nm].reshape(nm,2);dEp=y[10+2*nm:].reshape(nm,2)
        hat=np.maximum(0,1-np.abs(n-nodes)/dN)[:,None]
        d_dE=dEp
        d_dEp=-cc*dEp+bb*g[None,:]*dE+bb*q[None,:]*hat*E[None,:]
        return np.r_[D0p,bb*D0-aa*D0p,Ep,bb*g*E-cc*Ep,rhop,bb*g*rho-cc*rhop+bb*q,d_dE.ravel(),d_dEp.ravel()]
    y=y0;out={}
    for a,b in zip(range(M),range(1,M+1)):
        s=solve_ivp(rhs,(nodes[a],nodes[b]),y,method='DOP853',rtol=1e-11,atol=1e-13,max_step=dN)
        if not s.success: raise RuntimeError(s.message)
        y=s.y[:,-1]
        if b in ep_idx:out[b]=y.copy()
    return out,y
def build(mass):
    states,yend=reference_and_jacobian(mass)
    nm=M+1;Jm=[];O=[]
    for i in ep_idx:
        y=states[i];E=y[2:4];Ep=y[4:6]
        dE=y[10:10+2*nm].reshape(nm,2);dEp=y[10+2*nm:].reshape(nm,2)
        O+= list(E*E-1)+list(Ep/E)
        P=2*E[None,:]*dE;F=(dEp-(Ep/E)[None,:]*dE)/E[None,:]
        Jm+= [P[:,0],P[:,1],F[:,0],F[:,1]]
    Jm=np.array(Jm)                                    # (8,nm) response per unit hat height
    rho=yend[6:8];rhop=yend[8:10];Emax=yend[2:4];
    # epoch values for u=E', v=E
    ep_vals=[(states[i][2:4],states[i][4:6],states[i][8:10]) for i in ep_idx]
    return np.array(O),Jm,rho,rhop,Emax,ep_vals
def remainder(rho,rhop,Emax,ep_vals):
    sig=ETA*rho;sigp=ETA*rhop
    rem=[]
    for (v,u,epoch_rhop) in ep_vals:
        for typ in ('P','F'):
            for k in range(2):
                E=Emax[k];s=sig[k];sp=ETA*epoch_rhop[k]
                if not 0 <= s < 1:
                    raise ValueError("Neumann contraction condition failed")
                r2=E*s*s/(1-s);r2p=sp*E*s/(1-s)
                if typ=='P':
                    rem.append((s*E)**2+2*E*r2+2*s*E*r2+r2**2)
                else:
                    dv=s*E+r2;du=sp*E+r2p;umax=u[k]
                    if dv >= 1:
                        raise ValueError("Quotient denominator bound failed")
                    rem.append(r2p+umax*r2+(du*dv+umax*dv*dv)/(1-dv))
    # reorder to [P1,P2,F1,F2]*2 epochs
    order=[]
    for e in range(2):
        order+=[e*4+0,e*4+1,e*4+2,e*4+3]
    return np.array([rem[i] for i in order])
def target():
    J.MASS=M0;J.BASE=BREF
    if FRAC>0:
        J.MASS=.15  # frozen JH1 two-mass target: light .15, heavy .30
        return J.integrate(np.zeros(8),target_fraction=FRAC)[0]
    return J.integrate(np.zeros(8))[0]
def lp_box(c,Jm):
    nv=M+1
    cvec=np.r_[np.zeros(nv),1.]
    A=[];b=[]
    for j in range(8):
        A.append(np.r_[Jm[j],-1.]);b.append(-c[j])
        A.append(np.r_[-Jm[j],-1.]);b.append(c[j])
    for i in range(M):
        for sg in (1,-1):
            r=np.zeros(nv+1);r[i+1]=sg;r[i]=-sg;A.append(r);b.append(LB*dN)
    res=linprog(cvec,A_ub=np.array(A),b_ub=b,bounds=[(-ETA,ETA)]*nv+[(0,None)],method='highs')
    if not res.success:
        raise RuntimeError(res.message)
    # Weak duality with a box correction: no stationarity equality is assumed.
    # y <= 0 for A x <= b; ensure the coefficient of the unbounded t >= 0
    # is nonnegative, and minimize remaining residual coefficients over the box.
    y=np.minimum(res.ineqlin.marginals,0.)
    weight=-y[:16].sum()
    if weight > 0:
        y *= min(1.,(1.-1e-12)/weight)
    reduced=cvec-np.asarray(A).T@y
    if not np.all(np.isfinite(reduced)):
        raise ArithmeticError("Nonfinite dual stationarity residual")
    if reduced[-1] < 0:
        raise ArithmeticError("Dual t coefficient is negative")
    dual=float(np.asarray(b)@y-ETA*np.abs(reduced[:nv]).sum())
    if not np.isfinite(dual):
        raise ArithmeticError("Nonfinite dual objective")
    gap=float(res.fun-dual)
    if gap < -1e-10:
        raise ArithmeticError("Inconsistent primal/dual objectives")
    return dual,res.x[:nv],float(res.fun),gap
def nonlinear(db,mass):
    J.MASS=mass;J.BASE=BREF
    bn=BREF+db
    def bf(n):return float(np.interp(n,nodes,bn))
    y=np.array([1,J.initial_growth_rate(),1,1,0,0],float)
    def f(n,y):
        a_,b_=J.background(n);e=y[2:4];ep=y[4:6]
        q=J.K**2/(J.K**2+np.exp(2*n)*mass**2);cc=2*y[1]/y[0]+a_
        return np.r_[y[1],b_*y[0]-a_*y[1],ep,b_*bf(n)*q*e-cc*ep]
    # Restart at each interpolation knot: adaptive stages must not skip
    # slope discontinuities and underestimate their integration error.
    out=[];n0=N0
    for n1 in nodes[1:]:
        so=solve_ivp(f,(n0,n1),y,method='DOP853',rtol=1e-11,atol=1e-13)
        if not so.success: raise RuntimeError(so.message)
        y=so.y[:,-1];n0=n1
        if any(abs(n1-epoch)<1e-12 for epoch in J.EPOCHS):
            e=y[2:4];ep=y[4:6];out.append(np.r_[e*e-1,ep/e])
    if len(out) != len(J.EPOCHS):
        raise ValueError("History mesh must contain every output epoch")
    return np.concatenate(out)

def main():
    Ot=target();rows=[]
    for r in RLIST:
        m=r*M0;O,Jm,rho,rhop,Emax,epv=build(m)
        c=O-Ot;t_lp,dbopt,primal,dual_gap=lp_box(c,Jm)
        rem=remainder(rho,rhop,Emax,epv)
        # Uniform interpolation bound LB*dN/2, valid with variable positive
        # response kernels; an unweighted mean error is insufficient here.
        assert Jm.min()>=-1e-14
        allow=(LB*dN/2)*Jm.sum(axis=1)
        cert=t_lp-rem.max()-allow.max()
        nl=float(np.abs(nonlinear(dbopt,m)-Ot).max())
        rows.append(dict(mass_ratio=r,LP_box=float(primal),LP_dual_lower=float(t_lp),LP_duality_gap=dual_gap,remainder=float(rem.max()),interp_allowance=float(allow.max()),estimated_lower_bound=float(max(0.,cert)),nonlinear_residual_at_box_LP_optimum=nl,J_min=float(Jm.min()),sigma=(ETA*rho).tolist()))
        print(rows[-1],flush=True)
    json.dump(dict(status=f'target m0={M0} heavy_fraction={FRAC}; analytic lower-bound formula evaluated numerically; no validated numerical enclosure; class |b-0.02|<=0.005,|b-prime|<=0.05',M=M,eta=ETA,Lb=LB,rows=rows),open(HERE/f'rigorous_lp_certificate_M{M}{TAG}.json','w'),indent=2)
if __name__ == '__main__':
    main()
