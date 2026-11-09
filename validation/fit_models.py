"""Stage 1.9c / Stage 2 preview: best-fit chi2 for all models with likelihood_v2."""
import json, sys, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
from rlds_paired import likelihood_v2 as L

MODELS = {
 "LCDM": (L.lcdm, ["Om","H0","wb"], [(0.2,0.45),(55,80),(0.019,0.026)], [[0.31,67.5,0.0224],[0.30,68.5,0.0223]]),
 "wCDM": (L.wcdm, ["Om","H0","wb","w"], [(0.2,0.45),(55,80),(0.019,0.026),(-1.5,-0.5)], [[0.31,67.5,0.0224,-1.0],[0.31,66.5,0.0224,-0.95]]),
 "CPL":  (L.cpl,  ["Om","H0","wb","w0","wa"], [(0.2,0.45),(55,80),(0.019,0.026),(-2,0),(-3,2)], [[0.31,67.5,0.0224,-1.0,0.0],[0.31,67.,0.0224,-0.85,-0.6]]),
 "RLDS-free": (L.rlds_free, ["Om","H0","wb","K","Rc","dR","alpha"],
      [(0.2,0.45),(55,80),(0.019,0.026),(0,0.25),(0.3,0.6),(0.03,0.2),(2.5,3.5)],
      [[0.31,67.5,0.0224,0.01,0.45,0.12,3.0],[0.31,67.,0.0224,0.01,0.565,0.17,2.9],[0.31,67.5,0.0224,0.0,0.45,0.12,3.0]]),
 "RLDS-shooting": (L.rlds_shooting, ["K","H0","wb","Rc","dR","alpha"],
      [(0.005,0.25),(55,80),(0.019,0.026),(0.3,0.6),(0.03,0.2),(3.01,3.3)],
      [[0.0096,67.8,0.0224,0.599,0.199,3.0101]]),
}
def fit(args):
    name, use = args
    fn, pn, bounds, starts = MODELS[name]
    lo=np.array([b[0] for b in bounds]); hi=np.array([b[1] for b in bounds])
    sn = L.load_sn() if "sn" in use else None
    def f(x):
        if np.any(x<lo) or np.any(x>hi): return 1e10
        if name=="CPL" and x[3]+x[4]>=0: return 1e10
        try: return L.chi2_total(fn(*x), sn, use)
        except Exception: return 1e10
    best=None; t=time.time()
    mf = 400 if name=="RLDS-shooting" else 4000
    for s in starts:
        r=minimize(f,s,method="Powell",bounds=list(zip(lo,hi)),options=dict(xtol=1e-5,ftol=1e-6,maxfev=mf))
        r=minimize(f,r.x,method="Nelder-Mead",options=dict(xatol=1e-6,fatol=1e-6,maxfev=mf))
        if best is None or r.fun<best.fun: best=r
    bg=fn(*best.x)
    parts=dict(cmb=L.chi2_cmb(bg),bao=L.chi2_bao(bg))
    if sn is not None: parts["sn"]=L.chi2_sn(bg,sn)
    out=dict(model=name,data="+".join(use),chi2=float(best.fun),params=dict(zip(pn,map(float,best.x))),
             Om_today=bg.Omega_m_today, omega_m_early=bg.omega_m_early, r_d=float(bg.r_d), parts=parts, sec=round(time.time()-t))
    print(json.dumps(out),flush=True); return out
if __name__=="__main__":
    which = sys.argv[1].split(",")
    combos=[("cmb","bao"),("cmb","bao","sn")]
    jobs=[(m,c) for m in which for c in combos]
    with Pool(2) as p: res=p.map(fit,jobs)
    json.dump(res,open(f"validation/fit_{'_'.join(which)}.json","w"),indent=1)
