import json, numpy as np
R=json.load(open("validation/fit_LCDM_wCDM_CPL.json"))+[r for r in json.load(open("validation/fit_RLDS-free_RLDS-shooting.json")) if r["model"]!="RLDS-shooting"]+json.load(open("validation/fit_shooting_refined.json"))
k={"LCDM":3,"wCDM":4,"CPL":5,"RLDS-shooting":6,"RLDS-free":7}
N={"cmb+bao":16,"cmb+bao+sn":1606}
for d in ("cmb+bao","cmb+bao+sn"):
    base=[r for r in R if r["model"]=="LCDM" and r["data"]==d][0]["chi2"]
    print(f"\n{d}  (chi2_LCDM = {base:.2f}, N = {N[d]})")
    print(f"{'model':15s} {'dchi2':>7s} {'dAIC':>7s} {'dBIC':>7s}  key params")
    for m in ["LCDM","wCDM","CPL","RLDS-shooting","RLDS-free"]:
        r=[r for r in R if r["model"]==m and r["data"]==d][0]; dk=k[m]-3; dc=r["chi2"]-base
        p={kk:round(v,4) for kk,v in r["params"].items() if kk not in ("H0","wb")}
        print(f"{m:15s} {dc:+7.2f} {dc+2*dk:+7.2f} {dc+dk*np.log(N[d]):+7.2f}  {p}  Om_today={r['Om_today']:.3f}")
