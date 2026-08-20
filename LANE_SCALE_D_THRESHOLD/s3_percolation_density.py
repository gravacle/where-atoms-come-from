"""S3 -- IS THERE A PERCOLATION THRESHOLD IN k, AND DOES IT MOVE WITH THE VENUE'S OWN SCALE?

D-17 in force.  A "giant component appears at k = k*" claim is worthless until k* is measured
at SEVERAL bath sizes.  If k* moves with the number of bath sites n_B, the onset is a DENSITY
threshold (records per site), not a threshold in the number of records.

The relation: two records INTERACT if their supports touch the SAME bath site.  Qubits are
assigned to bath sites two ways -- ROUND-ROBIN (q -> q mod n_B) and BLOCK (contiguous chunks)
-- because the assignment is venue data and a threshold that survives only one assignment is
an artifact of that assignment.

Both families appear in the SAME table (D-15).  Family B with BLOCK assignment and n_B >= m is
the carrier on which the effect CANNOT occur: it must return ncomp = k/2 and no giant
component, and it does -- that is the negative arm; family A in the same rows is the positive
arm that registers a giant component immediately.
"""
import sys, time, json
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD")
import numpy as np
import carriers as C
import battery as B

OUT = []
def P(s=""):
    print(s); OUT.append(s)

def site_map(n, nB, mode):
    if mode == "roundrobin": return [q % nB for q in range(n)]
    per = -(-n // nB)
    return [min(q // per, nB - 1) for q in range(n)]

def share_graph(car, fam, nB, mode):
    n = car["n"]; k = len(fam)
    sm = site_map(n, nB, mode)
    sites = [frozenset(sm[q] for q in C.support(v, n)) for v in fam]
    A = [[1 if (i != j and sites[i] & sites[j]) else 0 for j in range(k)] for i in range(k)]
    return A, sites

def run_family(carfn, ks, nBs, modes, tag):
    rows = []
    for kk in ks:
        car = carfn(kk)
        fam, part, ch = C.records_of(car)
        if not C.all_checks_pass(ch):
            P("  SELF-CHECK FAILED at %s -- no conclusion" % car["label"]); continue
        k = len(fam); n = car["n"]
        for mode in modes:
            for nB in nBs:
                if nB > n: continue
                A, sites = share_graph(car, fam, nB, mode)
                g = B.graph_scalars(A, k, "sh")
                occ = {}
                for s in sites:
                    for x in s: occ[x] = occ.get(x, 0) + 1
                rows.append(dict(tag=tag, label=car["label"], n=n, k=k, nB=nB, mode=mode,
                                 recs_per_site=(sum(occ.values()) / float(nB)),
                                 ncomp=g["sh_ncomp"], giant=g["sh_giant"],
                                 edgefrac=g["sh_edgefrac"], lapgap=g["sh_lapgap"],
                                 percolates=g["sh_percolates"]))
    return rows

t0 = time.time()
P("=" * 150)
P("S3  PERCOLATION OF THE RECORD-INTERACTION GRAPH vs k, AT SEVERAL BATH SIZES (D-17)")
P("=" * 150)

nBs = [2, 3, 4, 6, 8, 12, 16, 24]
ksA = [4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48]          # n for family A -> k = n-2
ksB = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24]          # m for family B -> k = 2m
rows = run_family(C.family_A, ksA, nBs, ("roundrobin", "block"), "A")
rows += run_family(C.family_B, ksB, nBs, ("roundrobin", "block"), "B")
ksC = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24]          # m for family C -> k = m
rows += run_family(C.family_C, ksC, nBs, ("roundrobin", "block"), "C")

P()
P("TABLE  interaction = two records touch the SAME bath site.  A and B in the same table.")
hdr = "  %-12s %-16s %-4s %-4s %-11s %-12s %-7s %-8s %-10s %-9s %s" % (
    "family", "carrier", "k", "nB", "assignment", "recs/site", "ncomp", "giant", "edgefrac",
    "lapgap", "percolates")
P(hdr); P("  " + "-" * (len(hdr) + 8))
for r in rows:
    P("  %-12s %-16s %-4d %-4d %-11s %-12.3f %-7d %-8.3f %-10.4f %-9.4f %s" %
      (r["tag"], r["label"], r["k"], r["nB"], r["mode"], r["recs_per_site"], r["ncomp"],
       r["giant"], r["edgefrac"], r["lapgap"], r["percolates"]))

P()
P("k*  = SMALLEST k at which the interaction graph percolates (giant component > k/2).")
P("     If k* MOVES with nB the onset is a DENSITY threshold, not a record-count threshold.")
P("  %-8s %-12s %-6s %-8s %s" % ("family", "assignment", "nB", "k*", "k*/nB"))
P("  " + "-" * 52)
summary = []
for tag in ("A", "B", "C"):
    for mode in ("roundrobin", "block"):
        for nB in nBs:
            cand = sorted([r["k"] for r in rows
                           if r["tag"] == tag and r["mode"] == mode and r["nB"] == nB
                           and r["percolates"]])
            ks_all = sorted(set(r["k"] for r in rows if r["tag"] == tag and r["mode"] == mode
                                and r["nB"] == nB))
            if not ks_all: continue
            kstar = cand[0] if cand else None
            summary.append((tag, mode, nB, kstar))
            P("  %-8s %-12s %-6d %-8s %s" % (tag, mode, nB, kstar if kstar else "none in range",
                                             ("%.2f" % (kstar / float(nB))) if kstar else "-"))

json.dump(rows, open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s3_rows.json", "w"))
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s3_percolation_density.txt",
     "w").write("\n".join(OUT) + "\n")
P("total %.1fs" % (time.time() - t0))
