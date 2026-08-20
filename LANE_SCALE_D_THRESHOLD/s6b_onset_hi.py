"""S6 -- THE ONE PLACE AN ONSET FROM EXACTLY ZERO CAN LIVE: PROTECTION vs k, AT FIXED n.

The structured families A, B, C hold k rigidly tied to n, so they cannot separate "a threshold
in the NUMBER of records" from "a threshold in the DENSITY of records".  This script breaks the
tie: it holds the carrier size n FIXED and sweeps k from 1 to n-1 over an ensemble of RANDOM
stabiliser carriers, and asks whether the quantity

    reach_r  =  dim_F2 of the record group that ONE region of r sites can write

turns on from EXACTLY ZERO at some k.  reach_r = 0 is clause (v) holding at radius r; reach_r > 0
is clause (v) failing at radius r.  That is a genuine qualitative change, not a trend.

D-17 IS THE WHOLE POINT HERE: the onset k* is measured at FIVE carrier sizes n.  If k* is the
same NUMBER at every n, it is a record-count threshold and the principal's hypothesis stands.
If k* scales with n, it is a DENSITY threshold and the number of records is not what matters.

The structured families are carried in the same table as the positive control: they have
reach_2 > 0 at every k, so a zero in the random ensemble is a real zero and not a dead probe.
"""
import sys, time, json
from itertools import combinations
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD")
import numpy as np
import carriers as C
import battery as B
from f2 import rank

OUT = []
def P(s=""):
    print(s); OUT.append(s)

def random_code(n, k, rng, ntrans=None):
    """A random stabiliser group of rank n-k.  Rejection sampling cannot reach n-k generators
       (a random Pauli commutes with j chosen ones only with probability 2^-j), so the group is
       built by TRANSPORT instead: start from the manifestly valid commuting independent set
       {Z_1,...,Z_(n-k)} and push it through a random element of Sp(2n,F_2), assembled from
       random transvections T_v(w) = w + <w,v> v.  Transvections preserve the symplectic form,
       so the image is still mutually commuting and still independent -- both re-checked below,
       and a failure is reported rather than used."""
    need = n - k
    S = []
    for i in range(need):
        v = [0] * (2 * n); v[n + i] = 1; S.append(v)
    if ntrans is None: ntrans = 3 * n
    for _ in range(ntrans):
        v = [int(x) for x in rng.integers(0, 2, size=2 * n)]
        if not any(v): continue
        S = [[(w[j] + C.sp(w, v, n) * v[j]) % 2 for j in range(2 * n)] for w in S]
    if any(C.sp(S[i], S[j], n) for i in range(need) for j in range(i + 1, need)): return None
    if rank(S, 2 * n) != need: return None
    return S

def reach_random_regions(car, r, nsamp, rng):
    """max over sampled size-r regions of local_logical_dim -- exact per region."""
    n = car["n"]
    best = 0
    if r <= 2 or len(list(combinations(range(n), r))) <= nsamp:
        regs = [set(c) for c in combinations(range(n), r)]
    else:
        regs = [set(rng.choice(n, size=r, replace=False).tolist()) for _ in range(nsamp)]
    for reg in regs:
        d = B.local_logical_dim(car, reg)
        if d > best: best = d
    return best

t0 = time.time()
P("=" * 150)
P("S6b HIGHER-STATISTICS ONSET OF LOCAL WRITABILITY vs k AT FIXED CARRIER SIZE n  (random stabiliser ensemble)")
P("=" * 150)
NENS = 32
rng = np.random.default_rng(2026)
NS = [10, 12, 16, 20, 24, 28]
RS = (1, 2, 3)
rows = []
spot = []
for n in NS:
    for k in range(1, n):
        agg = {r: [] for r in RS}
        made = 0
        for e in range(NENS):
            S = random_code(n, k, rng)
            if S is None: continue
            car = dict(name="R", label="rand[[%d,%d]]" % (n, k), n=n, dim=2 ** n, stabs=S)
            # D-18 SPOT CHECK: on a sample of these random carriers, actually construct the
            # record family and run the same symbolic clause self-checks S1 used.  reach_r
            # itself needs only the stabiliser group, so it is not run per carrier.
            if e == 0 and k in (1, n // 2, n - 1):
                fam, part, ch = C.records_of(car)
                spot.append((n, k, len(fam) == k, C.all_checks_pass(ch)))
            made += 1
            for r in RS:
                agg[r].append(reach_random_regions(car, r, 120, rng))
        if made == 0:
            P("  n=%d k=%d: no valid carrier built -- SKIPPED, no conclusion" % (n, k)); continue
        row = dict(n=n, k=k, made=made)
        for r in RS:
            v = agg[r]
            row["frac_pos_r%d" % r] = float(np.mean([1.0 if x > 0 else 0.0 for x in v]))
            row["mean_reach_r%d" % r] = float(np.mean(v))
        rows.append(row)
    P("  ... n=%d done  t=%.1fs" % (n, time.time() - t0))

P()
P("  D-18 SPOT CHECKS on the random ensemble (record family built and re-checked):")
for n, k, kok, chok in spot:
    P("     n=%-3d k=%-3d  family size correct=%-6s  symbolic self-checks=%s" % (n, k, kok, chok))
if not all(a and b for _, _, a, b in spot):
    P("  A SPOT CHECK FAILED -- the random-ensemble table below is NOT reported.")
    rows = []

P()
P("TABLE  fraction of random carriers in which ONE region of r sites can write ANY record")
P("       (0.0 = clause (v) holds at radius r for every carrier in the ensemble)")
P("  %-5s %-5s %-7s %-13s %-13s %-13s %-13s %-13s %-13s" %
  ("n", "k", "k/n", "P(reach1>0)", "P(reach2>0)", "P(reach3>0)", "mean_reach1",
   "mean_reach2", "mean_reach3"))
P("  " + "-" * 104)
for r in rows:
    P("  %-5d %-5d %-7.3f %-13.3f %-13.3f %-13.3f %-13.3f %-13.3f %-13.3f" %
      (r["n"], r["k"], r["k"] / r["n"], r["frac_pos_r1"], r["frac_pos_r2"], r["frac_pos_r3"],
       r["mean_reach_r1"], r["mean_reach_r2"], r["mean_reach_r3"]))

P()
P("POSITIVE CONTROL in the same table -- the structured families, where reach_2 > 0 always:")
P("  %-16s %-5s %-5s %-13s %-13s %-13s" % ("carrier", "n", "k", "reach1", "reach2", "reach3"))
P("  " + "-" * 70)
for mk, ps in ((C.family_A, (4, 8, 12, 16)), (C.family_B, (1, 2, 4, 6)), (C.family_C, (1, 2, 3, 4))):
    for p in ps:
        car = mk(p); fam, part, ch = C.records_of(car)
        vals = [reach_random_regions(car, r, 120, rng) for r in RS]
        P("  %-16s %-5d %-5d %-13d %-13d %-13d" % (car["label"], car["n"], len(fam), *vals))

P()
P("ONSET k*(n, r) = SMALLEST k at which more than half the ensemble has reach_r > 0.")
P("  D-17 READ: if k* is the SAME NUMBER at every n it is a RECORD-COUNT threshold;")
P("            if k*/n is the same FRACTION at every n it is a DENSITY threshold.")
P("  %-5s %-5s %-8s %-9s" % ("n", "r", "k*", "k*/n"))
P("  " + "-" * 32)
onsets = []
for n in NS:
    for r in RS:
        cand = [x["k"] for x in rows if x["n"] == n and x["frac_pos_r%d" % r] > 0.5]
        ks = min(cand) if cand else None
        onsets.append((n, r, ks))
        P("  %-5d %-5d %-8s %-9s" % (n, r, ks if ks else "none in range",
                                     ("%.3f" % (ks / float(n))) if ks else "-"))

json.dump(dict(rows=rows, onsets=onsets),
          open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s6b_rows.json", "w"))
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s6b_onset_random_carriers_hi.txt",
     "w").write("\n".join(OUT) + "\n")
P("total %.1fs" % (time.time() - t0))
