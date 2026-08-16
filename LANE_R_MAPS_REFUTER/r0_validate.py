"""R0 -- rebuild the corpus from scratch and reproduce the numbers the CLAIM cites,
before attacking anything.  Every number below is computed here, none inherited."""
import numpy as np, sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_R_MAPS_REFUTER")
from rmlib import *

np.set_printoptions(precision=12, suppress=False)
def hr(t): print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)

hr("R0.1  K1 -- PUBLISHED INCIDENCE MATRICES")
K = K1()
print("vertices:", K.verts)
print("edges   :", [(e[0], e[1] + '->' + e[2]) for e in K.edges])
print("faces   :", K.faces)
print("loopF   :", K.loopF, "  loopC:", K.loopC)
print("\nd1 (rows=v0..v4, cols=e1..e6):\n", K.d1)
print("\nd2 (rows=e1..e6, cols=F):\n", K.d2.T, " (transposed for display)")
print("\nd1 @ d2  max|.| =", K.check_d2_zero(), "   (d^2 = 0)")
print("loopF closed? d1.chain max|.| =", K.loop_is_closed(K.loopF))
print("loopC closed? d1.chain max|.| =", K.loop_is_closed(K.loopC))
b, r = K.betti()
print("chi =", K.chi(), "  rank d1 =", r[0], " rank d2 =", r[1], "  betti =", b)

hr("R0.2  TRANSPORT BY EXPLICIT MATRICES vs THE CLASS FORMULA")
# S1 sec.6 worked instance and S3's headline connection
for tag, aedge in [
    ("S1 worked instance a=(pi/3 x3, pi/2 x3)",
     dict(e1=np.pi/3, e2=np.pi/3, e3=np.pi/3, e4=np.pi/2, e5=np.pi/2, e6=np.pi/2)),
    ("S3 headline f=2.0 c=1.1",
     dict(e1=2.0, e2=0.0, e3=0.0, e4=1.1, e5=0.0, e6=0.0)),
]:
    WF, WC = K.holon(aedge)
    f, c = np.angle(WF), np.angle(WC)
    dev = max(abs(K.Z_matrix(WF, WC, k) - Z_from_pi(K.pi(), f, c, k)) for k in range(1, 60))
    print(f"{tag}:  W_F={WF:.12f}  W_C={WC:.12f}")
    print(f"   max|Z_k(matrix) - Z_k(class formula)| over k<=59 = {dev:.3e}")

hr("R0.3  K1 UNDER SENSE U (uniform vertex weights) -- the CLAIM's normalisation")
piU = K.pi()
print("pi (00,10,01,11) =", piU, "   sum =", piU.sum())
nS, rG, S = rank_G(piU)
print("|S| =", nS, "  rank G =", rG, "  S =", S, "  FORMS =", forms(piU))
f, c = 2.0, 1.1
lA = lambda_A(piU, f, c)
print(f"lambda_A at (f,c)=(2.0,1.1) = {lA:.12f}")
print("   CLAIM cites                  -1.536645461686")
print("   check via matrices           %.12f" % float(np.log(abs(K.Z_matrix(np.exp(1j*f), np.exp(1j*c), 1)))))

print("\nlambda_B, four independent routes:")
gf, gc = 1.0 + np.sqrt(2.0), 1.0 + np.sqrt(3.0)     # generic (non-resonant) point
print("  generic point used for the direct Weyl average: (f,c) = (1+sqrt2, 1+sqrt3)")
print("  direct  k_n=n  N=4e6      : %.12f" % lambda_B_direct(piU, gf, gc, 4_000_000))
print("  grid    3000x3000 midpoint: %.12f" % mahler2(piU, 3000))
print("  grid    6001x6001 midpoint: %.12f" % mahler2(piU, 6001))
print("  MonteCarlo N=2e7 seed 20260816: %.12f" % mahler2_mc(piU))
print("  Cassaigne-Maillot (IMPORT) : %.12f" % cassaigne_maillot(piU[3], piU[1], piU[2]))
print("  CLAIM cites                  -0.756573585640")
print("  S4 sec.4.4 CONTROL 1 cites   -0.756573585634   (of record)")

hr("R0.4  K1 UNDER THE CORPUS'S OWN PUBLISHED READY STATE (SENSE C)")
Kc = K1(p=dict(v0=0.4, v1=0.15, v2=0.15, v3=0.15, v4=0.15), name="K1(S3 p)")
piC = Kc.pi()
print("p = (0.4,0.15,0.15,0.15,0.15)  ->  pi =", piC)
print("  Cassaigne-Maillot          : %.12f" % cassaigne_maillot(piC[3], piC[1], piC[2]))
print("  grid 6001                  : %.12f" % mahler2(piC, 6001))
print("  register erratum of record  -0.767507880")

hr("R0.5  THE CLAIM'S RING-TORUS NUMBER, IDENTIFIED")
x = -0.613104472886
print("claim's ring-torus lambda_B =", x, "   exp =", np.exp(x))
print("log(13/24) = %.12f" % np.log(13/24), "   <-- exact match")
print("S4 of-record ring tori: B0a -0.747659833   B0b -0.810930216 = log(4/9) = %.12f"
      % np.log(4/9))
print("=> the claim's ring torus is NEITHER B0a NOR B0b.  It is a THIRD object,")
print("   and log(max weight) is the Mahler measure precisely when one class")
print("   weight EXCEEDS 1/2 (no-triangle branch): 13/24 = %.6f > 1/2." % (13/24))
