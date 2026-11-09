"""Refine the shooting best fit from the ridge point found by shooting_scan.py."""
import json, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
from rlds_paired import likelihood_v2 as L
lo=np.array([0.005,55,0.019,0.3,0.03,3.01]); hi=np.array([0.25,80,0.026,0.6,0.2,3.3])
def job(use):
    sn=L.load_sn() if "sn" in use else None
    f=lambda x: 1e10 if (np.any(x<lo) or np.any(x>hi)) else L.chi2_total(L.rlds_shooting(*x),sn,use)
    x0=np.array([0.01,68.82,0.02263,0.599,0.199,3.0101])
    simplex=np.vstack([x0]+[x0+d for d in np.diag([0.001,0.3,0.0001,-0.01,-0.005,0.0005])])
    r=minimize(f,x0,method="Nelder-Mead",options=dict(initial_simplex=simplex,xatol=1e-6,fatol=1e-4,maxfev=900))
    bg=L.rlds_shooting(*r.x)
    out=dict(model="RLDS-shooting",data="+".join(use),chi2=float(r.fun),
             params=dict(zip(["K","H0","wb","Rc","dR","alpha"],map(float,r.x))),
             Om_today=bg.Omega_m_today,omega_m_early=bg.omega_m_early)
    print(json.dumps(out),flush=True); return out
if __name__=="__main__":
    with Pool(2) as p: res=p.map(job,[("cmb","bao"),("cmb","bao","sn")])
    json.dump(res,open("validation/fit_shooting_refined.json","w"),indent=1)
