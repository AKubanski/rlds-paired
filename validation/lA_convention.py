"""Stage 1.9b: size and stability of the l_A convention offset vs Chen+2019."""
import numpy as np
from classy import Class
from rlds_paired import likelihood_v2 as L
# Planck 2018 TT,TE,EE+lowE means (Planck 2018 VI, Table 2): wb, wc, H0
pts = {"Planck18 TTTEEE+lowE mean": (0.02236, 0.1202, 67.27),
       "shifted wc +2sig": (0.02236, 0.1230, 66.0),
       "shifted wb +2sig": (0.02266, 0.1202, 67.8),
       "wCDM-like low H0": (0.02236, 0.1202, 64.0)}
for name,(wb,wc,H0) in pts.items():
    h=H0/100; Om=(wb+wc+L.OMEGA_NU)/h**2
    bg=L.lcdm(Om,H0,wb); R,lA,_=bg.cmb_vector()
    c=Class(); c.set({'h':h,'omega_b':wb,'omega_cdm':wc,'N_ur':2.0328,'N_ncdm':1,'m_ncdm':0.06}); c.compute()
    th=c.get_current_derived_parameters(['100*theta_star'])['100*theta_star']
    print(f"{name:28s} R={R:.5f}  l_A={lA:.3f}  pi/theta*_CLASS={np.pi/(th/100):.3f}  ratio ours/(pi/theta*)={lA/(np.pi/(th/100)):.6f}")
