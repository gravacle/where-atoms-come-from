"""V1 -- ADVERSARIAL CHECK OF THE LANE'S reach_r NUMBERS.

Two independent things are tested here, with my OWN F_2 code (no import of battery.py except
for the cross-check in part 0):

  PART 0  my reach_r agrees with the lane's battery.local_logical_dim on carriers where the
          lane enumerated regions EXHAUSTIVELY.  If it does not, nothing below is reported.
  PART 1  the lane's s6b "POSITIVE CONTROL" table reports reach3 = 2 for [[4,2,2]]^6 while its
          own master table s8 reports reach3 = 4 for the same carrier.  reach3 is a MAXIMUM
          over regions, so the smaller number cannot be right.  Recomputed EXHAUSTIVELY.
  PART 2  the estimator: s6/s6b call reach_random_regions(car, r, nsamp=120), which enumerates
          all regions only when C(n,r) <= 120 and otherwise SAMPLES 120 subsets.  For r = 3
          that is exhaustive only at n = 10 (C(10,3) = 120); coverage then falls
          220 -> 3276 while nsamp stays 120.  A max estimated from a shrinking fraction of the
          candidates is biased DOWNWARD, and the bias GROWS WITH n -- which is the same
          direction as the lane's D-17 kill "k* grows with n at every radius".
          Measured here directly: sampled-120 vs EXHAUSTIVE on the same carriers.
"""
import sys, time
from itertools import combinations
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD")
import numpy as np
import carriers as C
import battery as BAT

OUT = []
def P(s=""):
    print(s); OUT.append(s)

# ---------------------------------------------------------------- my own F_2 machinery
def _rank(rows, ncol):
    rows = [r[:] for r in rows]; r = 0
    for c in range(ncol):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [(a + b) % 2 for a, b in zip(rows[i], rows[r])]
        r += 1
    return r

def my_reach(car, region):
    """dim_F2 of the logical group writable by Paulis supported inside `region`.
       Built independently: columns = single-qubit X and Z generators on the region;
       admissible = kernel of the syndrome map; subtract stabilisers wholly inside region."""
    n = car["n"]; S = car["stabs"]; idx = sorted(region)
    cols = []
    for q in idx:
        e = [0] * (2 * n); e[q] = 1; cols.append(e)
        e = [0] * (2 * n); e[n + q] = 1; cols.append(e)
    M = [[C.sp(c, s, n) for c in cols] for s in S]
    nul = len(cols) - _rank(M, len(cols))
    reg = set(region)
    inside = [[s[q] for q in idx] + [s[n + q] for q in idx]
              for s in S if C.support(s, n) <= reg]
    return nul - _rank(inside, 2 * len(idx))

def reach_exhaustive(car, r):
    n = car["n"]
    return max(my_reach(car, set(c)) for c in combinations(range(n), r))

def reach_sampled(car, r, nsamp, rng):
    """VERBATIM the lane's estimator (s6b lines 53-63)."""
    n = car["n"]
    if r <= 2 or len(list(combinations(range(n), r))) <= nsamp:
        regs = [set(c) for c in combinations(range(n), r)]
    else:
        regs = [set(rng.choice(n, size=r, replace=False).tolist()) for _ in range(nsamp)]
    return max(my_reach(car, reg) for reg in regs)

t0 = time.time()
P("=" * 120)
P("V1  EXHAUSTIVE reach_r  vs  THE LANE'S SAMPLED ESTIMATOR")
P("=" * 120)

# ---------------------------------------------------------------- PART 0 cross-check
P()
P("PART 0  my independent reach vs the lane's battery.local_logical_dim, region by region")
P("        (if these disagree anywhere, NOTHING below is reported)")
bad = 0; tested = 0
for car in (C.family_A(6), C.family_B(2), C.family_C(2)):
    n = car["n"]
    for r in (1, 2, 3):
        for c in combinations(range(n), r):
            a = my_reach(car, set(c)); b = BAT.local_logical_dim(car, set(c))
            tested += 1
            if a != b: bad += 1
P("        %d regions compared, %d disagreements -> %s" % (tested, bad, "OK" if bad == 0 else "ABORT"))
if bad:
    P("CROSS-CHECK FAILED -- no conclusion drawn."); sys.exit(0)

# ---------------------------------------------------------------- PART 1 the contradiction
P()
P("PART 1  THE LANE CONTRADICTS ITSELF ON reach3.  Same carrier, two of its own tables.")
P("        s8_master_table.txt  (contiguous regions, exhaustive)   B n=24 k=12 : reach3 = 4")
P("        s6b_...hi.txt        (POSITIVE CONTROL table, nsamp=120): [[4,2,2]]^6 : reach3 = 2")
P("        reach3 is a MAX over regions, so 2 is impossible if any region gives 4.")
P()
P("        %-16s %-5s %-5s %-14s %-14s %-14s" %
  ("carrier", "n", "k", "reach3_EXHAUST", "reach3_contig", "reach3_samp120"))
P("        " + "-" * 74)
rng = np.random.default_rng(2026)
for car in (C.family_B(4), C.family_B(6), C.family_C(3), C.family_C(4)):
    n = car["n"]
    ex = reach_exhaustive(car, 3)
    co = max(my_reach(car, reg) for reg in BAT.regions(car, 3))
    sa = reach_sampled(car, 3, 120, rng)
    P("        %-16s %-5d %-5d %-14d %-14d %-14d"
      % (car["label"], n, n - len(car["stabs"]), ex, co, sa))

# ---------------------------------------------------------------- PART 2 the estimator bias
P()
P("PART 2  ESTIMATOR COVERAGE FOR r = 3 IN THE LANE'S ONSET SWEEP (nsamp = 120 FIXED)")
P("        %-5s %-12s %-12s %-10s" % ("n", "C(n,3)", "sampled", "coverage"))
P("        " + "-" * 42)
from math import comb
for n in (10, 12, 16, 20, 24, 28):
    tot = comb(n, 3)
    P("        %-5d %-12d %-12s %-10.3f" % (n, tot, "ALL" if tot <= 120 else "120",
                                            1.0 if tot <= 120 else 120.0 / tot))
P("        coverage falls 15x from n=12 to n=28 while the lane reads k*(r=3) as a function of n.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/VERIFY/v1_reach_exhaustive.txt",
     "w").write("\n".join(OUT) + "\n")
P("total %.1fs" % (time.time() - t0))
