"""Stage 1.9a: background quantities of likelihood_v2 vs CLASS (Planck 2018 LCDM)."""
import numpy as np
from classy import Class
from rlds_paired import likelihood_v2 as L
h, ob, oc = 0.6736, 0.02237, 0.1200
c = Class(); c.set({'h':h,'omega_b':ob,'omega_cdm':oc,'N_ur':2.0328,'N_ncdm':1,'m_ncdm':0.06}); c.compute()
Om = c.Omega_m()
bg = L.lcdm(Om, 100*h, ob)
print(f"Omega_m = {Om:.5f}; omega_m_early = {bg.omega_m_early:.5f}")
rows=[]
for z in [0.1,0.5,1.0,2.33,10,100,1089.9]:
    dm_c = c.angular_distance(z)*(1+z); H_c = c.Hubble(z)*L.C_LIGHT/(100*h)
    rows.append((z, bg.DM(z)/dm_c-1, bg.E(z)/H_c-1))
    print(f"z={z:7.2f}  D_M rel.diff {rows[-1][1]:+.2e}   E rel.diff {rows[-1][2]:+.2e}")
d = c.get_current_derived_parameters(['rs_d','z_d','rs_star','z_star'])
print(f"r_d: Aubourg {bg.r_d:.3f}  CLASS {d['rs_d']:.3f}  rel {bg.r_d/d['rs_d']-1:+.2e}")
# r_s at the same redshift: CLASS comoving sound horizon at z
zs = bg.z_star
bk = c.get_background(); zz=bk['z']; rs=bk['comov.snd.hrz.']
rs_c = np.interp(zs, zz[::-1], rs[::-1])
print(f"z_* (Hu-Sugiyama) = {zs:.2f};  r_s(z_*): ours {bg.r_s(zs):.3f}  CLASS {rs_c:.3f}  rel {bg.r_s(zs)/rs_c-1:+.2e}")
v = bg.cmb_vector(); print("CMB vector (R, l_A, omega_b):", v, " prior means:", L.CMB_MEAN)
print("chi2_CMB at Planck best fit:", L.chi2_cmb(bg))
