"""LANE_SCALE_C_ACCUM  --  script 3: what does a BOUNDED REGION know about the records?

Gravity's defining feature is that a surface's flux grows with what is ENCLOSED.  The record
analogue that can actually be measured on these carriers is:

    reach(A)  =  the F_2 rank of  { v(P) : P a Pauli supported inside region A },
                 where v(P)_i = 1 iff P anticommutes with record R_i

-- how many independent record-bits an operation confined to A can disturb.  It is EXACT and
costs nothing: v is F_2-BILINEAR in P, so reach(A) is the rank of the 2|A| x k matrix of
sp(e_j, R_i) over the coordinates j living in A.  BASIS-INVARIANT: a change of record basis
right-multiplies by an element of GL(k,2) and does not move a rank.

A FIRST ATTEMPT AT THIS SCRIPT ASSUMED reach depended only on |A|.  A self-check refuted that
exhaustively (C(12,2), C(12,3), C(12,4) all give more than one value), so every number below is
reported as a RANGE over regions of that size -- min and max, exhaustive where the binomial fits,
sampled otherwise and labelled.

CONTROL: m independent [[4,2,2]] blocks at the same k.
"""
import sys, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_C_ACCUM")
import numpy as np
from math import comb
from record_model import symplectic_logicals
from s1_combinatorics import carrier_nn2, carrier_product, embed, sp2, f2_rank

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.append(s)

EXHAUST = 300000

def coord_images(n, R):
    k = len(R); rows = []
    for j in range(2 * n):
        e = [0] * (2 * n); e[j] = 1
        rows.append([sp2(e, R[i], n) for i in range(k)])
    return rows

def reach(rows, n, qubits):
    idx = list(qubits) + [n + q for q in qubits]
    return f2_rank([rows[j] for j in idx])

def reach_range(rows, n, r, rng, samples=4000):
    """(min, max, how) of reach over regions of size r -- exhaustive when C(n,r) fits."""
    if r > n:                                   # region bigger than the carrier: use the carrier
        v = reach(rows, n, range(n)); return (v, v, "capped")
    if comb(n, r) <= EXHAUST:
        vals = [reach(rows, n, s) for s in itertools.combinations(range(n), r)]
        return (min(vals), max(vals), "exact")
    vals = [reach(rows, n, range(r))]
    for _ in range(samples):
        vals.append(reach(rows, n, sorted(rng.choice(n, size=r, replace=False).tolist())))
    return (min(vals), max(vals), "sampled")

def cover_min(rows, n, k, rng):
    """smallest region size that can disturb EVERY record.  Exact where the search is
       exhaustive; otherwise an upper bound from sampling + a greedy run."""
    lo = -(-k // 2)
    for r in range(lo, n + 1):
        mn, mx, how = reach_range(rows, n, r, rng)
        if mx == k: return r, how
    return None, "none"

def build(kind, par):
    if kind == "family":
        n, S = carrier_nn2(par)
        pr = symplectic_logicals([s[:] for s in S], n)
        return n, S, [p[1] for p in pr]
    n, S = carrier_product(par)
    pr4 = symplectic_logicals([s[:] for s in carrier_product(1)[1]], 4)
    return n, S, [embed(p[1], 4, 4 * b, n) for b in range(par) for p in pr4]

def run():
    P("=" * 118)
    P("LANE_SCALE_C_ACCUM  script 3  --  reach(A): how many record-bits can a bounded region disturb?")
    P("=" * 118)
    rng = np.random.default_rng(0)
    CH = []
    NS = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22]
    F, C = {}, {}
    for n in NS:
        nn, S, R = build("family", n); F[n - 2] = (nn, R, coord_images(nn, R))
    for m in range(1, 11):
        nn, S, R = build("control", m); C[2 * m] = (nn, R, coord_images(nn, R))

    for tag, D in (("family", F), ("control", C)):
        for k, (nn, R, rows) in D.items():
            CH.append((tag, k, "full register reaches rank k", f2_rank(rows) == k))
            CH.append((tag, k, "reach(A) <= 2|A| for every contiguous A",
                       all(reach(rows, nn, range(r)) <= 2 * r for r in range(1, nn + 1))))
            CH.append((tag, k, "reach(empty region) = 0", reach(rows, nn, []) == 0))
    # the refuted assumption, kept as a standing check so it cannot creep back in
    for tag, D, kk in (("family", F, 10), ("control", C, 10)):
        nn, R, rows = D[kk]
        vals = set(reach(rows, nn, s) for s in itertools.combinations(range(nn), 3))
        CH.append((tag, kk, "REFUTED ASSUMPTION: reach is NOT a function of |A| alone "
                            "(this check passes when it varies)", len(vals) > 1))
    bad = [c for c in CH if not c[3]]
    P("")
    P("SELF-CHECKS")
    P("-" * 118)
    for c in bad: P(f"   FAIL {c[0]} k={c[1]} {c[2]}")
    P(f"   {len(CH)-len(bad)} / {len(CH)} pass" + ("   -- ALL PASS" if not bad else "   -- SOME FAILED"))
    if bad:
        P("   CONCLUSIONS VOID."); return

    RS = [1, 2, 3, 4, 6, 8]
    P("")
    P("TABLE R1   reach over ALL regions of size r, as min-max.  'e' = exhaustive over every")
    P("           region of that size, 's' = 4000 sampled regions plus the contiguous one.")
    P("-" * 118)
    P(f"{'k':>3} | {'n':>3} " + " ".join(f"{'r='+str(r):>9}" for r in RS)
      + f" || {'nC':>3} " + " ".join(f"{'r='+str(r):>9}" for r in RS))
    P("-" * 118)
    RR = {}
    for n in NS:
        k = n - 2
        def cells(D):
            nn, R, rows = D[k]
            out = f"{nn:>3} "
            for r in RS:
                mn, mx, how = reach_range(rows, nn, r, rng)
                RR[(id(D), k, r)] = (mn, mx)
                out += f" {str(mn)+'-'+str(mx)+how[0]:>9}"
            return out
        P(f"{k:>3} | {cells(F)} || {cells(C)}")
    P("-" * 118)

    P("")
    P("TABLE R2   THE GAUSS-LAW QUESTION.  Hold the region SIZE fixed and let the number of")
    P("           enclosed records grow.  Gravity-shaped behaviour would mean a fixed region")
    P("           learns MORE as k grows.  Reported: the BEST region of that size (max reach).")
    P("-" * 118)
    P(f"{'k':>3} | {'fam best r=2':>12} {'fam best r=4':>12} {'fam best r=6':>12} {'frac k (r=6)':>13}"
      f" || {'ctl best r=2':>12} {'ctl best r=4':>12} {'ctl best r=6':>12} {'frac k (r=6)':>13}"
      f" || {'min region for ALL k':>22}")
    for n in NS:
        k = n - 2
        nnf, _, rf = F[k]; nnc, _, rc = C[k]
        bf = [reach_range(rf, nnf, r, rng)[1] for r in (2, 4, 6)]
        bc = [reach_range(rc, nnc, r, rng)[1] for r in (2, 4, 6)]
        cf, hf = cover_min(rf, nnf, k, rng)
        cc, hc = cover_min(rc, nnc, k, rng)
        P(f"{k:>3} | {bf[0]:>12} {bf[1]:>12} {bf[2]:>12} {bf[2]/k:>13.4f}"
          f" || {bc[0]:>12} {bc[1]:>12} {bc[2]:>12} {bc[2]/k:>13.4f}"
          f" || fam {str(cf)+hf[0]:>6}  ctl {str(cc)+hc[0]:>6}")
    P("-" * 118)
    P("frac k = the fraction of the record content a size-6 region can touch")

    P("")
    P("READ  (filled from the numbers above)")
    P("-" * 118)
    ks = [n - 2 for n in NS]
    b6f = [reach_range(F[k][2], F[k][0], 6, rng)[1] for k in ks]
    b6c = [reach_range(C[k][2], C[k][0], 6, rng)[1] for k in ks]
    cvf = [cover_min(F[k][2], F[k][0], k, rng)[0] for k in ks]
    cvc = [cover_min(C[k][2], C[k][0], k, rng)[0] for k in ks]
    P(f"  best reach at FIXED region size r=6:  family {b6f}")
    P(f"                                        control {b6c}")
    P(f"  smallest region that touches ALL k:   family {cvf}   (hard lower bound {[-(-k//2) for k in ks]})")
    P(f"                                        control {cvc}")
    P("")
    P("Largest carrier: n = 22 (k = 20) family, n = 40 (k = 20) control.  Nothing stopped this")
    P("computation at that size; it is the F_2 rank of a 2n x k matrix.  The range was chosen to")
    P("match scripts 1 and 2 so the columns line up.")

if __name__ == "__main__":
    run()
    with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_C_ACCUM/s3_region_reach.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
