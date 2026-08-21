"""S4 -- IS ANY OF IT A PROPERTY OF THE CARRIER, OR OF THE CHOSEN GENERATORS?

The record family is fixed only up to a SYMPLECTIC CHANGE OF BASIS of the logical group:
any Sp(2k,F_2) transform of {(X_i,Z_i)} is an equally valid conjugate basis, and every record
in the new family still satisfies (i)-(iv) (the clauses are properties of the class, and the
transformed operators are still in N(S)\\S and still trace-balanced on every eigenspace).

So a "threshold at k" that appears in ONE basis and not others is a property of the
Gram-Schmidt output, not of the carrier.  This script re-measures the basis-DEPENDENT order
parameters over an ensemble of random symplectic bases (random transvections, which generate
Sp(2k,F_2)) and reports min / mean / max at every k, alongside the canonical value.

A quantity only counts as a candidate threshold if the SAME qualitative change survives the
whole ensemble.
"""
import sys, time, json
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD")
import numpy as np
import carriers as C
import battery as B

OUT = []
def P(s=""):
    print(s); OUT.append(s)

def logical_coords(k):
    """coordinates on the logical group: 2k bits, X_1..X_k then Z_1..Z_k"""
    return 2 * k

def sp_log(a, b, k):
    return sum(a[i] * b[k + i] + a[k + i] * b[i] for i in range(k)) % 2

def random_symplectic_family(fam, part, car, rng, ntrans=None):
    """Apply a random element of Sp(2k,F_2) (a product of transvections) to the logical
       generators and return the new commuting family, as Pauli (x|z) vectors."""
    n = car["n"]; k = len(fam)
    gens = [list(v) for v in fam] + [list(v) for v in part]          # Pauli reps
    coord = [[1 if j == i else 0 for j in range(2 * k)] for i in range(2 * k)]
    if ntrans is None: ntrans = 4 * k
    for _ in range(ntrans):
        v = [int(x) for x in rng.integers(0, 2, size=2 * k)]
        if not any(v): continue
        for i in range(2 * k):
            s = sp_log(coord[i], v, k)
            if s:
                coord[i] = [(coord[i][j] + v[j]) % 2 for j in range(2 * k)]
    # SELF-CHECK: the transform must preserve the symplectic form (pairing = identity)
    ok = all(sp_log(coord[i], coord[k + j], k) == (1 if i == j else 0)
             for i in range(k) for j in range(k)) and \
         all(sp_log(coord[i], coord[j], k) == 0 for i in range(k) for j in range(i + 1, k))
    newfam = []
    for i in range(k):
        v = [0] * (2 * n)
        for t in range(2 * k):
            if coord[i][t]:
                v = [(a + b) % 2 for a, b in zip(v, gens[t])]
        newfam.append(C.min_weight_rep(v, car))
    # SELF-CHECK: the new family must still be mutually commuting and in N(S)\S
    ok = ok and all(C.sp(newfam[i], newfam[j], n) == 0
                    for i in range(k) for j in range(i + 1, k))
    ok = ok and all(C.sp(r, s, n) == 0 for r in newfam for s in car["stabs"])
    ok = ok and all(any(r) for r in newfam)
    return newfam, ok

def scalars(car, fam):
    k = len(fam); n = car["n"]
    out = {}
    out.update(B.graph_scalars(B.overlap_adj(car, fam), k, "ov"))
    A, mx, ms = B.coflip(car, fam, 2)
    out.update(B.graph_scalars(A, k, "cf2"))
    out["cf2_maxset"] = mx; out["cf2_maxsingle"] = ms
    A1, mx1, ms1 = B.coflip(car, fam, 1)
    out["cf1_maxsingle"] = ms1
    out["wbar"] = sum(C.weight(v, n) for v in fam) / float(k)
    return out

KEYS = ["wbar", "ov_edgefrac", "ov_ncomp", "ov_giant", "ov_clus", "ov_lapgap",
        "cf2_edgefrac", "cf2_ncomp", "cf2_giant", "cf2_clus", "cf2_maxset", "cf2_maxsingle",
        "cf1_maxsingle"]

t0 = time.time()
P("=" * 150)
P("S4  BASIS SENSITIVITY OF THE BASIS-DEPENDENT ORDER PARAMETERS")
P("=" * 150)
P("  ensemble: 24 random symplectic bases per carrier (products of 4k random transvections),")
P("  every one re-checked for (a) symplectic form preserved, (b) family still commuting,")
P("  (c) family still inside N(S), (d) no member trivial.  A basis that fails is discarded and")
P("  reported, not silently used.")
NENS = 24
rng = np.random.default_rng(11)
rows = []
for tag, mk, params in (("A", C.family_A, [4, 6, 8, 10, 12, 16, 20, 24]),
                        ("B", C.family_B, [1, 2, 3, 4, 5, 6, 8, 10, 12]),
                        ("C", C.family_C, [1, 2, 3, 4, 5, 6, 8, 10, 12])):
    for p in params:
        car = mk(p)
        fam, part, ch = C.records_of(car)
        if not C.all_checks_pass(ch):
            P("  SELF-CHECK FAILED %s" % car["label"]); continue
        k = len(fam)
        canon = scalars(car, fam)
        ens = []; bad = 0
        for _ in range(NENS):
            nf, ok = random_symplectic_family(fam, part, car, rng)
            if not ok: bad += 1; continue
            ens.append(scalars(car, nf))
        row = dict(tag=tag, label=car["label"], k=k, nbad=bad, nens=len(ens))
        for key in KEYS:
            row[key + "_canon"] = canon[key]
            vals = [e[key] for e in ens]
            row[key + "_min"] = min(vals) if vals else None
            row[key + "_mean"] = float(np.mean(vals)) if vals else None
            row[key + "_max"] = max(vals) if vals else None
        rows.append(row)
        P("  measured %-14s k=%-3d ensemble=%d discarded=%d  t=%.1fs"
          % (car["label"], k, len(ens), bad, time.time() - t0))

for key in KEYS:
    P()
    P("ORDER PARAMETER  %s   -- canonical basis vs random-symplectic-basis ensemble" % key)
    P("  %-8s %-16s %-5s %-12s %-12s %-12s %-12s" %
      ("family", "carrier", "k", "canonical", "ens_min", "ens_mean", "ens_max"))
    P("  " + "-" * 84)
    for r in rows:
        def f(x):
            return ("%-12.4f" % x) if isinstance(x, float) else ("%-12s" % x)
        P("  %-8s %-16s %-5d %s %s %s %s" % (r["tag"], r["label"], r["k"],
                                             f(r[key + "_canon"]), f(r[key + "_min"]),
                                             f(r[key + "_mean"]), f(r[key + "_max"])))

json.dump(rows, open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s4_rows.json", "w"))
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s4_basis_sensitivity.txt",
     "w").write("\n".join(OUT) + "\n")
# T-35: post-write "total Ns" stdout line removed -- it was captured by reproduce.sh but never part of the sealed output, so it differed every run
