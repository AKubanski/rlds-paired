"""Extension of lA_convention.py (review of draft v4): stability of the l_A calibration
factor for late-time-modified expansion histories (wCDM, CPL), compared with CLASS."""
import numpy as np
from classy import Class
from rlds_paired import likelihood_v2 as L
pts = {"LCDM Planck mean": (0.02236, 0.1202, 67.27, -1.0, 0.0),
       "wCDM w=-0.9": (0.02236, 0.1202, 66.0, -0.9, 0.0),
       "CPL -0.85/-0.56": (0.02240, 0.1190, 67.5, -0.85, -0.56),
       "CPL -0.42/-1.75": (0.02240, 0.1190, 64.0, -0.42, -1.75),
       "phantom w=-1.1": (0.02236, 0.1202, 70.0, -1.1, 0.0)}
ratios = []
for name, (wb, wc, H0, w0, wa) in pts.items():
    h = H0 / 100; Om = (wb + wc + L.OMEGA_NU) / h ** 2
    bg = L.cpl(Om, H0, wb, w0, wa); R, lA, _ = bg.cmb_vector()
    c = Class()
    pars = {'h': h, 'omega_b': wb, 'omega_cdm': wc, 'N_ur': 2.0328, 'N_ncdm': 1, 'm_ncdm': 0.06}
    if (w0, wa) != (-1.0, 0.0):
        pars.update({'Omega_Lambda': 0.0, 'w0_fld': w0, 'wa_fld': wa, 'use_ppf': 'yes'})
    c.set(pars); c.compute()
    th = c.get_current_derived_parameters(['100*theta_star'])['100*theta_star']
    ratio = lA / (np.pi / (th / 100)); ratios.append(ratio)
    print(f"{name:20s} lA={lA:.3f}  pi/theta*={np.pi/(th/100):.3f}  ratio={ratio:.6f}")
    c.struct_cleanup(); c.empty()
print("spread of ratio: %.2e" % (max(ratios) - min(ratios)))
