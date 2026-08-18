"""O-7 SWEEP: fracton / no-string-logical-operator region against R1-R3.
Every number below is computed by exact F_2 linear algebra on the stabiliser group.
SELF-CHECKS run first and must all PASS before any headline number is read."""
import sys, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O78_FRACTON_NONABELIAN")
from gf2 import *

def hr(t=""): print("\n" + "="*78); print(t); print("="*78)

from models import *

# ---------------------------------------------------------------- report
def report(name, gens, n, want_k=None, do_d=True, trials=1500):
    bad = all_commute(gens, n)
    k, r = code_k(gens, n)
    rx, rz, rtot, is_css = css_split(gens, n)
    ok = (bad == 0)
    line = f"{name:34s} n={n:5d}  gens={len(gens):5d}  rank={r:5d}  k={k:4d}  CSS={is_css!s:5s}  commute_violations={bad}"
    print(line)
    if want_k is not None:
        print(f"    SELF-CHECK k == {want_k} : {'PASS' if k==want_k else 'FAIL'}   (got {k})")
        ok = ok and (k == want_k)
    d = None
    if do_d and k > 0:
        d, m = min_logical_weight(gens, n, trials=trials)
        mode = "EXACT (exhaustive)" if (1<<m) <= (1<<22) else f"UPPER BOUND (randomised, {trials} restarts)"
        print(f"    min logical weight = {d}   [{mode}]   dim S-perp = {m}")
    return dict(name=name, n=n, k=k, css=is_css, d=d, commute_ok=(bad==0))

hr("PART 0 -- SELF-CHECKS AND POSITIVE CONTROLS (must all PASS)")
checks = []
checks.append(report("CONTROL 2D toric code L=3", *toric2d(3), want_k=2))
checks.append(report("CONTROL 2D toric code L=4", *toric2d(4), want_k=2))
checks.append(report("CONTROL 2D toric code L=5", *toric2d(5), want_k=2, trials=800))
checks.append(report("CONTROL Steane [[7,1,3]]",  *steane(),   want_k=1))
checks.append(report("CONTROL perfect [[5,1,3]]", *perfect5(), want_k=1))
checks.append(report("CONTROL 3D toric code L=2", *toric3d(2), want_k=3))
checks.append(report("CONTROL 3D toric code L=3", *toric3d(3), want_k=3, trials=600))
allpass = all(c["commute_ok"] for c in checks)
print(f"\n  ALL STABILISER GROUPS COMMUTE: {'PASS' if allpass else 'FAIL'}")
print("  d(toric2d,L) should equal L; d(Steane)=3; d([[5,1,3]])=3 -- these are the")
print("  independently-known answers that validate the distance routine.")

hr("PART 1 -- X-CUBE  (type-I fracton, Vijay-Haah-Fu)")
for L in (2,3,4):
    rr = report(f"X-cube L={L}", *xcube(L), trials=400)
    print(f"    predicted log2 GSD = 6L-3 = {6*L-3}   measured k = {rr['k']}   "
          f"{'MATCH' if rr['k']==6*L-3 else 'MISMATCH'}")

hr("PART 2 -- CHECKERBOARD MODEL")
for L in (2,4):
    rr = report(f"checkerboard L={L}", *checkerboard(L), trials=400)
    print(f"    12L-6 = {12*L-6} ; 6L-3 = {6*L-3} ; measured k = {rr['k']}")

hr("PART 3 -- HAAH CUBIC CODE 1  (type-II, NO string logical operators)")
haahk = {}
for L in (2,3,4,5,6,7,8):
    gens, n = haah(L)
    bad = all_commute(gens, n)
    k, r = code_k(gens, n)
    rx, rz, rtot, is_css = css_split(gens, n)
    haahk[L] = k
    print(f"Haah CC1 L={L:2d}   n={n:5d}  rank={r:5d}  k={k:4d}  CSS={is_css}  commute_violations={bad}")
print("\n  Haah's published closed form for L = 2^m is  k = 4L - 2 :")
for L in (2,4,8):
    print(f"    L={L}: 4L-2 = {4*L-2}   measured k = {haahk[L]}   "
          f"{'MATCH' if haahk[L]==4*L-2 else 'MISMATCH'}")
print("\n  distance probe on the smallest instances (upper bounds):")
for L in (2,3,4):
    gens, n = haah(L)
    d, m = min_logical_weight(gens, n, trials=400)
    print(f"    L={L}: min logical weight <= {d}   (dim S-perp = {m})")

hr("PART 4 -- CHAMON MODEL  (NON-CSS fracton)")
for L in (4,6):
    for sub in (False, True):
        gens, n = chamon(L, sublattice=sub)
        bad = all_commute(gens, n)
        k, r = code_k(gens, n)
        rx, rz, rtot, is_css = css_split(gens, n)
        tag = "fcc sublattice" if sub else "all sites"
        print(f"Chamon L={L} ({tag:14s})  n={n:5d}  gens={len(gens):5d}  rank={r:5d}  "
              f"k={k:4d}  CSS={is_css}  commute_violations={bad}")

hr("PART 5 -- THE CHAIN-COMPLEX QUESTION (G-12 coverage)")
print("""A CSS code IS, by construction, a length-2 F_2 chain complex:
    C_2 --d2=H_Z^T--> C_1 --d1=H_X--> C_0 ,   d1 d2 = H_X H_Z^T = 0,
    logical space = H_1 = ker H_X / im H_Z^T,   dim = k.
So the test 'is this candidate inside G-12's class?' is EXACTLY the test 'is it CSS?'
(up to local Clifford).  Below: the CSS column decided per candidate, and the
chain-complex identity d1 d2 = 0 verified numerically for each CSS candidate.""")
def chain_check(name, gens, n):
    HX = [xpart(g,n) for g in gens if zpart(g,n)==0 and xpart(g,n)!=0]
    HZ = [zpart(g,n) for g in gens if xpart(g,n)==0 and zpart(g,n)!=0]
    if not HX or not HZ:
        print(f"  {name:28s} : NOT CSS in this basis -- no d1,d2 to form"); return
    bad = 0
    for a in HX:
        for b in HZ:
            bad += popcount(a & b) & 1
    kk = code_k(gens,n)[0]
    # dim H_1 = dim ker(HX) - rank(HZ)
    rX,_ = rank_gf2(HX); rZ,_ = rank_gf2(HZ)
    dimH1 = (n - rX) - rZ
    print(f"  {name:28s} : d1 d2 = 0 ? {'PASS' if bad==0 else 'FAIL'}   "
          f"dim H_1 = {dimH1}   k = {kk}   {'MATCH' if dimH1==kk else 'MISMATCH'}")
chain_check("2D toric L=3", *toric2d(3))
chain_check("3D toric L=2", *toric3d(2))
chain_check("X-cube L=3", *xcube(3))
chain_check("checkerboard L=4", *checkerboard(4))
chain_check("Haah CC1 L=4", *haah(4))
chain_check("Chamon L=4", *chamon(4))
print("\nDONE.")
