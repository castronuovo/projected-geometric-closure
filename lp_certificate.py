# SPDX-License-Identifier: BSD-3-Clause
"""First-order LP lower bound for mass distinguishability (frozen JH1 contract).

Necessity record: nonlinear fits (JH1-JH3) give only upper estimates. To first
order in the amplitude b(N), E=1+E1 with E1''+cc E1'=bb*q*b (E=1,E'=0 at b=0);
observables P~2*E1, f-f0~E1'. Slope-bounded positive histories (0.01<=b<=0.04,
|b'|<=2 min b) are a convex polytope for piecewise-linear b, so the minimal
max-residual over the class is an LP; its optimum is a LOWER BOUND within the
first-order model. A second-order remainder (not rigorously bounded here) is
reported by comparing the linear model with the nonlinear forward map.
Uses imported frozen background/initial data; floating point.
"""
import sys, json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linprog
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import bounded_joint_history as J
M0=.15;N0=J.START;M=400
nodes=np.linspace(N0,0,M+1);dN=nodes[1]-nodes[0]
# reference growth D0 (pure growing) on a fine grid
fine=np.linspace(N0,0,8001)
def rhs(n,y):
    aa,bb=J.background(n);return [y[1],bb*y[0]-aa*y[1]]
s=solve_ivp(rhs,(N0,0),[1,J.initial_growth_rate()],t_eval=fine,rtol=1e-12,atol=1e-14,method='DOP853')
D0,D0p=s.y
aa=np.array([J.background(n)[0] for n in fine]);bbv=np.array([J.background(n)[1] for n in fine])
cc=2*D0p/D0+aa
Phi=np.concatenate([[0],np.cumsum((cc[1:]+cc[:-1])/2*np.diff(fine))])
def hat(n,i):  # piecewise-linear hat on nodes
    return np.maximum(0,1-np.abs(n-nodes[i])/dN)
def operator(mass):
    q=J.K[:,None]**2/(J.K[:,None]**2+np.exp(2*fine)[None,:]*mass**2)  # (2,F)
    rows=[]
    H=np.array([hat(fine,i) for i in range(M+1)])    # (M+1,F)
    for T in J.EPOCHS:
        iT=np.searchsorted(fine,T-1e-12)
        sel=slice(0,iT+1)
        w=np.exp(-(Phi[iT]-Phi[sel]))*bbv[sel]               # kernel for E1'
        # E1'(T)=int w q b ds ; E1(T)=int_{N0}^{T} E1'(u)du = int (int_s^T exp(-(Phi(u)-Phi(s))) du) bb q b ds
        u=fine[sel]
        W=np.zeros(iT+1)
        for j in range(iT+1):
            e=np.exp(-(Phi[j:iT+1]-Phi[j]))
            W[j]=np.trapz(e,u[j:iT+1])*bbv[j]
        dfn=lambda wv,qk:np.trapz((wv*qk[sel])[None,:]*H[:,sel],u,axis=1)
        Pk=[2*dfn(W,q[k]) for k in range(2)];Fk=[dfn(w,q[k]) for k in range(2)]
        rows+= [Pk[0],Pk[1],Fk[0],Fk[1]]
    return np.array(rows)   # (8,M+1)
def lp(mass,target_lin,tol_nodes=True):
    A=operator(mass)
    nv=M+1
    c=np.r_[np.zeros(nv),1.]
    Aub=[];bub=[]
    for j in range(8):
        Aub.append(np.r_[A[j],-1.]);bub.append(target_lin[j])
        Aub.append(np.r_[-A[j],-1.]);bub.append(-target_lin[j])
    for i in range(M):
        for sgn in (1,-1):
            for k in (i,i+1):
                r=np.zeros(nv+1);r[i+1]+=sgn;r[i]-=sgn;r[k]-=2*dN;Aub.append(r);bub.append(0.)
    res=linprog(c,A_ub=np.array(Aub),b_ub=bub,bounds=[(.01,.04)]*nv+[(0,None)],method='highs')
    return res.fun,res.x[:nv]
def nonlinear(bnodes,mass):
    def bf(n):return float(np.interp(n,nodes,bnodes))
    y=np.array([1,J.initial_growth_rate(),1,1,0,0],float)
    def f(n,y):
        a_,b_=J.background(n);e=y[2:4];ep=y[4:6]
        q=J.K**2/(J.K**2+np.exp(2*n)*mass**2);cc_=2*y[1]/y[0]+a_
        return np.r_[y[1],b_*y[0]-a_*y[1],ep,b_*bf(n)*q*e-cc_*ep]
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
    bt=np.full(M+1,.02)
    Lt=operator(M0)@bt;Nt=nonlinear(bt,M0)
    rows=[dict(note='first-order vs nonlinear target',max_abs_target=float(np.abs(Nt).max()),max_second_order_gap=float(np.abs(Nt-Lt).max()))]
    print(rows[0],flush=True)
    for r in (.9,.95,.98,1.05,1.1,1.2,1.5):
        m=r*M0;tmin,bopt=lp(m,Lt)
        nl=nonlinear(bopt,m)-Nt
        lin_unfit=float(np.abs(operator(m)@bt-Lt).max())
        rows.append(dict(mass_ratio=r,LP_lower_bound_first_order=float(tmin),nonlinear_residual_at_LP_optimum=float(np.abs(nl).max()),first_order_unfit=lin_unfit,b_range=[float(bopt.min()),float(bopt.max())]))
        print(rows[-1],flush=True)
    json.dump(dict(status='first-order LP lower bound within linear model; second-order remainder not rigorously bounded; floating point',M=M,rows=rows),open(HERE/'lp_certificate_results.json','w'),indent=2)
main()
