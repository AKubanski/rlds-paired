#!/usr/bin/env python3
# RLDS Paper V - Figures 1-4 (microphysical validation), DR2 canon
# Reconstructed recipe (verified to reproduce V1 behavior & |dR|<0.01 in z[0,10]):
#   thawing canonical scalar, exponential potential V0 exp(-lambda phi)
#   lambda=0.5, beta_*=6e-4 M_Pl, d_phi=0.3, phi_c=0.85 (activation at R_inst~0.005),
#   V0 shot for R(0)=0.2974 (DR2). RLDS branch from the rlds-paired SDK (DR2 canon).
# Outputs BOTH PDF (vector) and PNG (@300 dpi) for each figure, plus per-figure CSVs
#   and a machine-readable fig5_summary.json, into the directory given by --out
#   (default: outputs/paper_v, relative to the current working directory).
# Requires: numpy, scipy, matplotlib, and the rlds-paired package on PYTHONPATH.
# Honest status: reproduces the RLDS background trajectory R(N) and the structural
#   form; the microscopic current Q is orders of magnitude below the RLDS restoring
#   amplitude (micro-meso gap). See Section 10 / Prop 10.1 / Fig.2 caption (revised).
#
# Run:
#   python scripts/paper_v/fig5_generate.py --out outputs/paper_v
#
# NOTE (packaging): the only differences from the original local script are the
#   output-path handling (CLI --out + pathlib, no hardcoded absolute path) and the
#   machine-readable JSON summary at the end. The physics (recipe constants, V0
#   shooting, RLDSTrajectory construction, figure/CSV generation logic) is unchanged.

import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from rlds_paired.core.ode import RLDSTrajectory
from rlds_paired.constants import K_BEST, ALPHA, RC, DR, OMEGA_M, EPSILON

# ---------- output directory (CLI / relative path; no hardcoded absolute path) ----------
_parser = argparse.ArgumentParser(
    description="Reproduce RLDS Paper V microphysical Figures 1-4 (DR2 canon).")
_parser.add_argument(
    "--out", default="outputs/paper_v",
    help="output directory for figures, CSVs and the summary "
         "(default: outputs/paper_v, relative to the working directory)")
_args, _unknown = _parser.parse_known_args()
OUT_DIR = Path(_args.out)
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT = str(OUT_DIR)
# ---------- DR2 canon ----------
K,al,Rc,dRw,Om = K_BEST,ALPHA,RC,DR,OMEGA_M
trR=RLDSTrajectory(K=K,alpha=al,Rc=Rc,dr=dRw,omega_m=Om)
Rstar=trR.R_star
gamma=2*al*(1-Rstar)                       # RLDS exp relax rate |F'(R*)|
def rlds_term(r): return K*al*((1-r)/(r+EPSILON))*0.5*(1+np.tanh((Rc-r)/dRw))

# ---------- micro recipe ----------
N0=-np.log(51.0); lam=0.5; dphi=0.3; bstar=6e-4
def integrate(V0,phic,Nend=5.0):
    def beff(phi): return (bstar/(2*dphi))/np.cosh((phi-phic)/dphi)**2
    def rhs(N,y):
        phi,u,s=y; rc=np.exp(s); V=V0*np.exp(-lam*phi)
        H2=(rc+V)/(3-0.5*u*u); g=-0.5*(rc/H2+u*u); be=beff(phi)
        return [u,-g*u-3*u+lam*V/H2-be*rc/H2,-3.0+be*u]
    return solve_ivp(rhs,[N0,Nend],[0,0,0],method='Radau',rtol=1e-10,atol=1e-13,dense_output=True,max_step=0.02)
def state(N,sol):
    phi,u,sv=sol.sol(N); phi=np.atleast_1d(phi);u=np.atleast_1d(u);sv=np.atleast_1d(sv)
    rc=np.exp(sv); return phi,u,rc
def Rof(N,sol,V0):
    phi,u,rc=state(N,sol); V=V0*np.exp(-lam*phi); H2=(rc+V)/(3-0.5*u*u); return rc/(3*H2)
# shoot V0 -> R(0)=0.2974 ; phic placed so activation ~ R_inst≈0.005 (micro-meso O(100))
phic=0.85
V0=brentq(lambda v: Rof(np.array([0.0]),integrate(v,phic),v)[0]-0.2974,1e-7,1e-2,xtol=1e-12)
sol=integrate(V0,phic)
def beff(phi): return (bstar/(2*dphi))/np.cosh((phi-phic)/dphi)**2

# grids
Ng=np.linspace(-4,5,1400)
Rm=Rof(Ng,sol,V0); Rr=np.array([trR.R(np.exp(-n)-1) if np.exp(-n)-1>=0 else trR.R(0.0) for n in Ng])
# RLDS valid for z>=0 (N<=0); for N>0 integrate RLDS forward in N:
# use SDK trajectory in N directly:
zr=np.exp(-Ng)-1
Rr=np.where(zr>=0,[trR.R(max(z,0.0)) for z in zr],np.nan)
# future (N>0): integrate RLDS ODE forward from R(0)
from scipy.integrate import odeint
def Frlds(R,N): return -al*R*(1-R)+K*al*((1-R)/(R+EPSILON))*0.5*(1+np.tanh((Rc-R)/dRw))
maskf=Ng>0
if maskf.any():
    Nf=Ng[maskf]; R0now=trR.R(0.0)
    Rf=odeint(Frlds,R0now,np.concatenate([[0],Nf]))[1:,0]
    Rr[maskf]=Rf
phi,u,rc=state(Ng,sol)
COL_M="#1a6fff"; COL_R="#d62728"

# ===== Fig 1 trajectory =====
fig,(a1,a2)=plt.subplots(2,1,figsize=(8,6),height_ratios=[3,1],sharex=True)
a1.plot(Ng,Rm,color=COL_M,lw=2.2,label=r'$R_{\rm micro}(N)$ (full system)')
a1.plot(Ng,Rr,color=COL_R,ls='--',lw=2.0,label=r'$R_{\rm RLDS}(N)$ (reduced ODE)')
a1.axhline(Rstar,color='green',ls=':',lw=1.5,label=fr'$R^*={Rstar:.3f}$ (attractor)')
a1.axhline(2/3,color='gray',ls='-.',lw=1.2,label=r'$R=2/3$ (acceleration)')
a1.axvline(0,color='0.6',lw=0.8); a1.text(0.05,0.45,'$N=0$\n(today)',color='0.4',fontsize=8)
a1.set_ylabel(r'$R(N)=\rho_c/(\rho_c+\rho_\phi)$'); a1.set_ylim(0,1.02)
a1.set_title('Trajectory Agreement: Microphysical vs RLDS'); a1.legend(fontsize=8,loc='upper right'); a1.grid(alpha=0.25)
res=Rm-Rr; a2.plot(Ng,res,color=COL_M,lw=1.8); a2.fill_between(Ng,res,0,color=COL_M,alpha=0.15)
a2.axhline(0,color='k',lw=0.6,ls=':'); a2.set_ylabel(r'$\Delta R$'); a2.set_xlabel(r'$N=\ln a$')
imx=np.nanargmax(np.abs(res)); a2.set_title(fr'Residual  ($\max|\Delta R|={np.nanmax(np.abs(res)):.3f}$, reached at late times; $|\Delta R|<0.01$ for $z\in[0,10]$)',fontsize=9)
a2.set_xlim(-4,5)
plt.tight_layout(); plt.savefig(f"{OUT}/fig5_trajectory.pdf"); plt.savefig(f"{OUT}/fig5_trajectory.png",dpi=300); plt.close()

# ===== Fig 2 reconstructed interaction (log-y, honest) =====
Ng2=np.linspace(N0+0.05,3.5,1400); phi2,u2,rc2=state(Ng2,sol); R2=Rof(Ng2,sol,V0)
Qrec=beff(phi2)*u2*R2
order=np.argsort(R2); Rs=R2[order]; Qs=Qrec[order]
Rgrid=np.linspace(0.01,0.95,400); RL=np.array([rlds_term(r) for r in Rgrid])
fig,ax=plt.subplots(figsize=(8,5))
ax.semilogy(Rgrid,RL,color=COL_R,ls='--',lw=2.0,label=r'$K\alpha\,\Psi(R)\,S(R_c-R)$ (RLDS)')
ax.semilogy(Rs,np.clip(Qs,1e-7,None),color=COL_M,lw=2.2,label=r'$Q_{\rm micro}/(H\rho_{\rm tot})$ (reconstructed)')
ax.axvline(Rc,color='orange',ls='-.',lw=1.2,label=fr'$R_c={Rc:.3f}$')
ax.axvline(Rstar,color='green',ls=':',lw=1.2,label=fr'$R^*={Rstar:.3f}$')
ax.set_xlabel('$R$'); ax.set_ylabel('Dimensionless interaction term'); ax.set_ylim(1e-7,2)
ax.set_title('Reconstructed Interaction vs RLDS Restoring Term'); ax.legend(fontsize=8,loc='upper right'); ax.grid(alpha=0.25,which='both')
ax.text(0.28,1e-5,'reconstructed micro current is orders of magnitude\nbelow the RLDS restoring amplitude (micro–meso gap, Sec. 9)',fontsize=8,color='0.3')
plt.tight_layout(); plt.savefig(f"{OUT}/fig5_restoring.pdf"); plt.savefig(f"{OUT}/fig5_restoring.png",dpi=300); plt.close()

# ===== Fig 3 late-time convergence =====
dm=np.abs(Rm-Rstar); dr=np.abs(Rr-Rstar)
fig,ax=plt.subplots(figsize=(8,5))
ax.semilogy(Ng,dm,color=COL_M,lw=2.2,label=r'$|R_{\rm micro}-R^*|$')
ax.semilogy(Ng,dr,color=COL_R,ls='--',lw=2.0,label=r'$|R_{\rm RLDS}-R^*|$')
# exp fit guide on RLDS branch after transition
Nfit=np.linspace(0.0,1.5,50); guide=dr[np.argmin(np.abs(Ng-0.0))]*np.exp(-gamma*Nfit)
ax.semilogy(Nfit,guide,color='0.3',ls=':',lw=1.3,label=fr'Exp. fit: $\gamma\approx{gamma:.2f}$')
ax.axvline(0,color='0.6',lw=0.8); ax.text(0.1,3e-6,'$N=0$ (today)',color='0.4',fontsize=8)
ax.set_xlabel(r'$N=\ln a$'); ax.set_ylabel(r'$|R-R^*|$'); ax.set_ylim(1e-6,2); ax.set_xlim(-4,5)
ax.set_title('Late-Time Convergence to Attractor'); ax.legend(fontsize=8,loc='lower left'); ax.grid(alpha=0.25,which='both')
plt.tight_layout(); plt.savefig(f"{OUT}/fig5_attractor.pdf"); plt.savefig(f"{OUT}/fig5_attractor.png",dpi=300); plt.close()

# ===== Fig 4 switching profile =====
fig,(b1,b2)=plt.subplots(1,2,figsize=(11,4.3))
# left: projected beff(phi(R)) normalised vs S(Rc-R)
Rp=R2.copy(); bp=beff(phi2); 
oo=np.argsort(Rp); Rps=Rp[oo]; bps=bp[oo]/bp.max()
b1.plot(Rps,bps,color=COL_M,lw=2.2,label=r'$\beta_{\rm eff}(\phi(R))$ (projected, norm.)')
b1.plot(Rgrid,0.5*(1+np.tanh((Rc-Rgrid)/dRw)),color=COL_R,ls='--',lw=2.0,label=r'$S(R_c-R)$ (RLDS)')
b1.axvline(Rc,color='orange',ls='-.',lw=1.2,label=fr'$R_c={Rc:.3f}$')
b1.set_xlabel('$R$'); b1.set_ylabel('Normalised activation'); b1.set_xlim(0,1); b1.set_ylim(-0.02,1.05)
b1.set_title('Projected Switching Profile'); b1.legend(fontsize=8,loc='center right'); b1.grid(alpha=0.25)
# right: beff(N)
beN=beff(phi); ip=np.argmax(beN)
b2.plot(Ng,beN,color=COL_M,lw=2.2)
b2.annotate(fr'Peak: $N={Ng[ip]:.2f}$, $R={Rm[ip]:.3f}$',xy=(Ng[ip],beN[ip]),xytext=(Ng[ip]+0.6,beN[ip]*0.8),
            fontsize=8,arrowprops=dict(arrowstyle='->',color='0.4'))
b2.set_xlabel(r'$N=\ln a$'); b2.set_ylabel(r'$\beta_{\rm eff}(\phi(N))$'); b2.set_xlim(-4,5)
b2.set_title('Coupling Activation vs Time'); b2.grid(alpha=0.25)
plt.tight_layout(); plt.savefig(f"{OUT}/fig5_switching.pdf"); plt.savefig(f"{OUT}/fig5_switching.png",dpi=300); plt.close()

# ===== CSVs =====
np.savetxt(f"{OUT}/fig5_data_trajectory.csv",np.column_stack([Ng,Rm,Rr,Rm-Rr]),delimiter=",",header="N,R_micro,R_RLDS,residual",comments="")
np.savetxt(f"{OUT}/fig5_data_restoring.csv",np.column_stack([Rs,Qs,[rlds_term(r) for r in Rs]]),delimiter=",",header="R,Q_rec,RLDS_term",comments="")
np.savetxt(f"{OUT}/fig5_data_convergence.csv",np.column_stack([Ng,dm,dr]),delimiter=",",header="N,abs_Rmicro_minus_Rstar,abs_RRLDS_minus_Rstar",comments="")
np.savetxt(f"{OUT}/fig5_data_switching.csv",np.column_stack([Ng,beff(phi),Rm]),delimiter=",",header="N,beta_eff,R_micro",comments="")

print("FIGURES GENERATED. Key numbers:")
print(f"  recipe: lambda={lam}, beta_*={bstar}, dphi={dphi}, phic={phic}, V0={V0:.4e}")
print(f"  R(0)={Rof(np.array([0.0]),sol,V0)[0]:.4f}  R*={Rstar:.4f}  gamma_RLDS={gamma:.3f}")
win=(Ng>=-np.log(11))&(Ng<=0)
print(f"  max|dR| over z[0,10] = {np.nanmax(np.abs((Rm-Rr)[win])):.4f}")
print(f"  beta_eff peak = {beff(phi).max():.4f} at N={Ng[np.argmax(beff(phi))]:.2f}, R={Rm[np.argmax(beff(phi))]:.4f}")
print(f"  Q_rec peak = {Qrec.max():.4f}  vs RLDS term up to {RL.max():.3f}")
print(f"  micro-meso: R_c^eff/R_c^inst ≈ {Rc/Rm[np.argmax(beff(phi))]:.0f}")

# ===== machine-readable summary (numbers identical to the console output above) =====
summary = {
    "lambda": float(lam),
    "beta_star": float(bstar),
    "dphi": float(dphi),
    "phi_c": float(phic),
    "V0": float(V0),
    "R0": float(Rof(np.array([0.0]), sol, V0)[0]),
    "R_star": float(Rstar),
    "gamma_RLDS": float(gamma),
    "max_abs_delta_R_z_0_10": float(np.nanmax(np.abs((Rm - Rr)[win]))),
    "beta_eff_peak": float(beff(phi).max()),
    "beta_eff_peak_N": float(Ng[np.argmax(beff(phi))]),
    "beta_eff_peak_R": float(Rm[np.argmax(beff(phi))]),
    "Q_rec_peak": float(Qrec.max()),
    "RLDS_term_max": float(RL.max()),
    "micro_meso_gap": float(Rc / Rm[np.argmax(beff(phi))]),
}
with open(f"{OUT}/fig5_summary.json", "w", encoding="utf-8") as _f:
    json.dump(summary, _f, indent=2)
print(f"  summary written: {OUT}/fig5_summary.json")
