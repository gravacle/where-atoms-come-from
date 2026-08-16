"""R2 -- close the two loopholes a defender could still reach for."""
import numpy as np, sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_R_MAPS_REFUTER")
from rmlib import *
def hr(t): print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)

hr("R2.1  THE CONTRACTION DOES NOT SMUGGLE A CONNECTION CHANGE (explicit gauge)")
print("""Contracting edge e is gauge-fixing a_e = 0.  Verified explicitly: build the gauge
theta_v that zeroes the contracted edge, apply a_e -> a_e + theta_tgt - theta_src,
and check both holonomies are unchanged to machine precision.""")
KU = K1()
rng = np.random.default_rng(271828182)          # SEED PUBLISHED
worst = 0.0
for _ in range(500):
    a = {e[0]: rng.uniform(0, 2*np.pi) for e in KU.edges}
    WF0, WC0 = KU.holon(a)
    for ename in ["e1", "e2", "e3", "e4", "e5", "e6"]:
        _, s, t = KU.edges[KU.ei[ename]]
        th = {v: 0.0 for v in KU.verts}
        th[t] = -a[ename]                        # zero out this one edge
        a2 = {en: a[en] + th[tg] - th[sc] for (en, sc, tg) in KU.edges}
        WF1, WC1 = KU.holon(a2)
        worst = max(worst, abs(WF1 - WF0), abs(WC1 - WC0), abs(a2[ename]))
print("  max over 500 random connections x 6 edges of")
print("  ( |W_F' - W_F| , |W_C' - W_C| , |a_e'| ) = %.3e" % worst)
print("  => the contraction is a legitimate move: holonomies EXACTLY preserved.")

hr("R2.2  THE lambda_A SHIFT IS NOT A FLUKE OF THE TEST POINT (2.0, 1.1)")
print("Sup over T^2 of |lambda_A(K1) - lambda_A(K1/e1)|, 2001x2001 midpoint grid")
print("(theta_j = 2pi(j+1/2)/2001), SENSE U and SENSE C:\n")
for tag, p in [("SENSE U (uniform)", None),
               ("SENSE C (S3's own p)", dict(v0=.4, v1=.15, v2=.15, v3=.15, v4=.15))]:
    K = K1(p=p); Kc = contract(K, "e1")
    n = 2001
    g = 2*np.pi*(np.arange(n)+0.5)/n
    F, C = np.meshgrid(g, g, indexing="ij")
    lo = np.log(np.abs(Z_from_pi(K.pi(), F, C, 1)))
    ln = np.log(np.abs(Z_from_pi(Kc.pi(), F, C, 1)))
    d = np.abs(lo - ln)
    print(f"  {tag:22s} sup|diff| = {d.max():.6f}   mean|diff| = {d.mean():.6f}   "
          f"frac{{diff>1e-9}} = {np.mean(d>1e-9):.4f}")
print("""
  lambda_A differs almost everywhere on the torus after ONE edge contraction, with a
  sup difference of order 1.  The (2.0,1.1) reading is typical, not selected.""")

hr("R2.3  FORMATION DIES WITH NO MAP AT ALL")
print("""K1 untouched -- 5 vertices, 6 edges, 1 face, Betti (1,1,0), chain map intact --
with the ready state moved onto the root alone (S3's own W-02 correction of record:
'four exhibited families never form, one of them on K1's own published ready state'):""")
for tag, p in [("p uniform          ", None),
               ("p all on v0 (root) ", dict(v0=1., v1=0., v2=0., v3=0., v4=0.)),
               ("p on v1,v2 only    ", dict(v0=0., v1=.5, v2=.5, v3=0., v4=0.))]:
    K = K1(p=p); pi = K.pi(); nS, rG, S = rank_G(pi)
    print(f"  {tag} pi={np.round(pi,4)} |S|={nS} rkG={rG} FORMS={forms(pi)!s:5s} "
          f"lambda_A(2.0,1.1)={lambda_A(pi,2.0,1.1):+.12f}")
print("""
  |S| 3 -> 1, rank G 2 -> 0, FORMS True -> False, lambda -> 0.000000000000 exactly:
  the CLAIM'S ENTIRE OUTCOME VECTOR, reproduced with the cell structure untouched and
  no map applied.  The tree collapse is not necessary for the outcome it is credited
  with, and therefore cannot be what causes it.""")
