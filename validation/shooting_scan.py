"""Is the shooting model's huge chi2 real or an optimizer failure? Targeted scan."""
import itertools, json, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
from rlds_paired import likelihood_v2 as L
def job(p):
    K,a,Rc,dR=p
    f=lambda x: L.chi2_total(L.rlds_shooting(K,x[0],x[1],Rc,dR,a),None,("cmb","bao")) if (55<x[0]<80 and 0.019<x[1]<0.026) else 1e10
    r=minimize(f,[68.0,0.0224],method="Nelder-Mead",options=dict(xatol=1e-4,fatol=1e-3,maxfev=150))
    return dict(K=K,alpha=a,Rc=Rc,dR=dR,chi2=float(r.fun),H0=float(r.x[0]),wb=float(r.x[1]))
if __name__=="__main__":
    grid=list(itertools.product([0.005,0.01],[3.01,3.03,3.1],[(0.6,0.2),(0.3,0.03)]))
    grid=[(K,a,rc,dr) for K,a,(rc,dr) in grid]
    with Pool(2) as p: res=p.map(job,grid)
    for r in sorted(res,key=lambda r:r["chi2"]): print(json.dumps(r))
    json.dump(res,open("validation/shooting_scan.json","w"),indent=1)
