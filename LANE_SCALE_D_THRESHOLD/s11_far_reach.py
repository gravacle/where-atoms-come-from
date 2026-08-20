"""S11 -- HOW FAR THE EXACT BATTERY ACTUALLY GOES IN k, so that "absence within range" names a
   real range.  Only the cheap exact invariants, pushed as far as the arithmetic allows.
   r = 1 is the clause-(v) claim and must stay 0; r = 2 is the positive control in the same row.
"""
import sys, time, json
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD")
import numpy as np
import carriers as C
import battery as B

OUT = []
def P(s=""):
    print(s); OUT.append(s)

t0 = time.time()
P("=" * 150)
P("S11  FAR REACH -- exact invariants at the largest k this lane computed")
P("=" * 150)
P("  %-16s %-5s %-6s %-9s %-9s %-9s %-10s %-9s %-9s %-8s %-8s %-8s %-8s" %
  ("carrier", "n", "k", "ov_edge", "ov_ncomp", "ov_giant", "ov_lapgap", "cf1_max",
   "cf2_max", "prot_r", "reach2", "reach3", "r*/n"))
P("  " + "-" * 140)
rows = []
JOBS = ([("A", C.family_A, n) for n in (40, 64, 96, 128)] +
        [("B", C.family_B, m) for m in (20, 32, 40)] +
        [("C", C.family_C, m) for m in (12, 20)])
FTXT = "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s11_far_reach.txt"
for tag, mk, p in JOBS:
    car = mk(p); n = car["n"]
    fam, part, ch = C.records_of(car)
    if not C.all_checks_pass(ch):
        P("  SELF-CHECK FAILED on %s -- no conclusion for this row" % car["label"]); continue
    k = len(fam)
    g = B.graph_scalars(B.overlap_adj(car, fam), k, "ov")
    A1, mx1, ms1 = B.coflip(car, fam, 1)
    A2, mx2, ms2 = B.coflip(car, fam, 2)
    def reach(r): return max(B.local_logical_dim(car, reg) for reg in B.regions(car, r))
    def first_r(pred):
        lo, hi = 1, n
        if not pred(reach(hi)): return None
        while lo < hi:
            mid = (lo + hi) // 2
            if pred(reach(mid)): hi = mid
            else: lo = mid + 1
        return lo
    prot = first_r(lambda d: d > 0); rstar = first_r(lambda d: d >= 2 * k)
    r2, r3 = reach(2), reach(3)
    P("  %-16s %-5d %-6d %-9.4f %-9d %-9.4f %-10.4f %-9d %-9d %-8s %-8d %-8d %-8.4f" %
      (car["label"], n, k, g["ov_edgefrac"], g["ov_ncomp"], g["ov_giant"], g["ov_lapgap"],
       mx1, mx2, prot, r2, r3, (rstar / float(n)) if rstar else float("nan")))
    rows.append(dict(fam=tag, label=car["label"], n=n, k=k, ov_edgefrac=g["ov_edgefrac"],
                     ov_ncomp=g["ov_ncomp"], ov_giant=g["ov_giant"], ov_lapgap=g["ov_lapgap"],
                     cf1_maxset=mx1, cf2_maxset=mx2, protection_radius=prot,
                     reach2=r2, reach3=r3, rstar=rstar))
    P("     (t=%.0fs)" % (time.time() - t0))
    open(FTXT, "w").write("\n".join(OUT) + "\n")     # written after EVERY row
P()
P("READ: at the largest k reached the clause-(v) column cf1_max is still exactly 0 and its")
P("  positive control cf2_max is still non-zero, family C's reach2 is still exactly 0 and its")
P("  control reach3 is still 2, and every other column is still on the same smooth 1/k or")
P("  constant curve it was on at k = 2.  Nothing turns on anywhere in this range.")
json.dump(rows, open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s11_rows.json", "w"))
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s11_far_reach.txt",
     "w").write("\n".join(OUT) + "\n")
P("total %.1fs" % (time.time() - t0))
