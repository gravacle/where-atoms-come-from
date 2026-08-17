# W10-D leg 7 -- the register rows that are CLAIMS ABOUT THE CORPUS'S OWN BYTES (W-04 ERR-3, ERR-4)
# and the class-occupancy census of every carrier the corpus has ever named.
import subprocess, os, numpy as np
R = "/Users/bgm/MB Work/where-atoms-come-from"

print("="*100)
print("== 7A  W-04 ERR-3 RE-RUN AT THE BYTES: 'grep -in gravit returns EXACTLY ONE LINE' ==")
print("="*100)
arts = [f for f in sorted(os.listdir(R)) if f.endswith(".md")]
tot = 0
for f in arts:
    out = subprocess.run(["grep", "-in", "gravit", os.path.join(R, f)],
                         capture_output=True, text=True).stdout.strip()
    nl = len(out.splitlines()) if out else 0
    tot += nl
    if nl:
        for L in out.splitlines():
            print(f"    {f}:{L[:110]}")
print(f"  sealed .md artifacts scanned: {len(arts)}   total 'gravit' lines: {tot}")
print("  ERR-3's count is stated over the 11 artifacts sealed at W-04.  The register has grown")
print("  since (W-07, W-08 rows quote the word 'gravit'? -- printed above if so).  The ERR-3")
print("  finding is a CORPUS fact, CARRIER-INDEPENDENT, and no carrier can bear on it.")
for w in ("Gauss law", "backreaction", "edge mode", "coupling constant", "plaquette"):
    n = subprocess.run(f"grep -ril '{w}' '{R}'/*.md | wc -l", shell=True,
                       capture_output=True, text=True).stdout.strip()
    print(f"    artifacts containing '{w}': {n}")

print("\n"+"="*100)
print("== 7B  W-04 ERR-4: S1's 'SMALLEST COMPLEX' CLAIM, AND WHAT ITS COUNTEREXAMPLE'S CLASSES ARE ==")
print("="*100)
print("  ERR-4's counterexample: two triangles sharing an EDGE, one filled.  Built and checked:")
Vn = 4  # a=0 b=1 c=2 d=3
edges = [(0, 1), (1, 2), (2, 0), (1, 3), (3, 0)]
E = len(edges)
d1 = np.zeros((Vn, E))
for k, (s, t) in enumerate(edges):
    d1[t, k] += 1; d1[s, k] -= 1
d2 = np.zeros((E, 1)); d2[0, 0] = 1; d2[1, 0] = 1; d2[2, 0] = 1      # face abc = e_ab + e_bc + e_ca
r1, r2 = np.linalg.matrix_rank(d1), np.linalg.matrix_rank(d2)
print(f"    V={Vn} E={E} F=1  chi={Vn-E+1}  rank d1={r1} rank d2={r2}  b0={Vn-r1} "
      f"b1={E-r1-r2} b2={1-r2}  invariants={E-r1}  max|d1.d2|={np.abs(d1@d2).max():.1e}")
print("    K1's profile: V=5 E=6 F=1 chi=0 b0=1 b1=1 b2=0 invariants=2.")
print("    SAME chi, b1, b2 and invariant count on ONE FEWER VERTEX AND EDGE.  ERR-4 CONFIRMED;")
print("    S1:43-44's stated reason for the carrier is false at the bytes.")
FVs, CVs = {0, 1, 2}, {0, 1, 3}
from collections import Counter
cc = Counter(('1' if v in FVs else '0')+('1' if v in CVs else '0') for v in range(Vn))
print(f"    AND ITS CLASS MULTISET: {dict(cc)}  ->  THREE classes, class 00 EMPTY.")
print("    The carrier W-04 offers as the one S1 should have used is ALSO three-class.  Nothing")
print("    in this corpus, including its own counterexamples, has ever occupied all four.")

print("\n"+"="*100)
print("== 7C  CLASS-OCCUPANCY CENSUS OF EVERY CARRIER THE CORPUS HAS NAMED ==")
print("="*100)
CENSUS = [("B0a ring torus, disjoint", "{00:2, 01:3, 10:4}", 3),
          ("B0b ring torus, MEETING ", "{00:4, 01:1, 10:2, 11:2}", 4),
          ("B3  horn torus          ", "{01:2, 10:2, 11:1}", 3),
          ("B1  K1 (as handed)      ", "{01:2, 10:2, 11:1}", 3),
          ("B4  SPINDLE             ", "{00:1, 01:1, 10:1, 11:3}", 4),
          ("B5  double sphere       ", "-- no gamma_C, no datum --", 0),
          ("B2  K1 both filled      ", "{01:2, 10:2, 11:1}", 3),
          ("B1p K1-bridged          ", "{01:3, 10:3}", 2),
          ("B1q K1-bridged+spectator", "{00:1, 01:3, 10:3}", 3),
          ("B1s K1 subdivided       ", "{01:5, 10:5, 11:1}", 3),
          ("ERR-4's edge-shared pair", "{01:1, 10:1, 11:2}", 3)]
for nm, ms, k in CENSUS:
    tag = "**FOUR CLASSES**" if k == 4 else (f"{k} classes" if k else "no datum")
    print(f"  {nm}  {ms:28s}  {tag}")
print(f"\n  FOUR-CLASS CARRIERS IN THE ENTIRE CORPUS: 2 of 11.  Both appear in ONE table (S4:519)")
print("  and neither was run through W-01, W-02, W-05's N-list, W-06, W-07 or W-08 by any lane")
print("  before W-09 and this one.  EVERY OTHER ROW OF THIS PROGRAM WAS DERIVED ON A CARRIER")
print("  WITH AT MOST THREE OCCUPIED CLASSES.")

print("\n"+"="*100)
print("== 7D  W-01's 'CORRECT TRIVIAL LIMIT', ON FOUR CLASSES ==")
print("="*100)
CLS = ('00', '10', '01', '11'); EXP = {'00': (0, 0), '10': (1, 0), '01': (0, 1), '11': (1, 1)}
for nm, p in [("B1  K1", np.array([0, 2/5, 2/5, 1/5])), ("B0b", np.array([4/9, 2/9, 1/9, 2/9])),
              ("B4 ", np.array([1/6, 1/6, 1/6, 3/6]))]:
    k = np.arange(1, 1001)
    Z = sum(p[i]*np.exp(1j*k*(-EXP[CLS[i]][0]*0.0 + EXP[CLS[i]][1]*0.0)) for i in range(4))
    print(f"  {nm}  W_F = W_C = 1:  Z_k = {Z[0]:.15f} for all k;  min|Z_k| = {np.abs(Z).min():.15f}")
print("  At the trivial connection every character is 1, so Z_k = sum p = 1 on ANY carrier and")
print("  with ANY class occupancy.  CARRIER-INDEPENDENT, and it is an identity, not a test.")
