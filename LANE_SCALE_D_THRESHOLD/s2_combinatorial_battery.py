"""S2 -- THE ORDER-PARAMETER BATTERY vs k, on TWO structurally different carrier families,
   in ONE table (D-15: the control column is not in a separate run).

   family A  [[n,n-2,2]]      one block, k = n-2   -- records SHARE the block
   family B  [[4,2,2]]^m      m blocks, k = 2m     -- records in different blocks share nothing
             B is the standing POSITIVE CONTROL for every quantity that A pins at a constant,
             and A is the positive control for every quantity B pins at a constant.

   Every number is exact.  Nothing is sampled.
"""
import sys, time, json
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD")
import numpy as np
import carriers as C
import battery as B

OUT = []
def P(s=""):
    print(s); OUT.append(s)

def measure(car, with_coflip=True, coflip_r=(1, 2, 3)):
    n = car["n"]
    fam, part, checks = C.records_of(car)
    k = len(fam)
    if not C.all_checks_pass(checks):
        return None
    row = dict(fam=car["name"], label=car["label"], n=n, k=k, log2dim=n)
    sup = [C.support(v, n) for v in fam]
    row["wbar"] = sum(len(s) for s in sup) / k
    row["wmax"] = max(len(s) for s in sup)
    row["wbar_over_n"] = row["wbar"] / n
    # ---- BASIS-DEPENDENT: support-overlap graph
    row.update(B.graph_scalars(B.overlap_adj(car, fam), k, "ov"))
    # ---- BASIS-DEPENDENT: co-flip graphs at region size r (r=1 is the clause-(v) claim)
    if with_coflip:
        for r in coflip_r:
            if r > n: continue
            A, maxset, maxsingle = B.coflip(car, fam, r)
            row.update(B.graph_scalars(A, k, "cf%d" % r))
            row["cf%d_maxset" % r] = maxset
            row["cf%d_maxsingle" % r] = maxsingle
    # ---- BASIS-INVARIANT: what one contiguous region can write
    for r in (1, 2, 3, 4):
        if r > n: continue
        d = max(B.local_logical_dim(car, reg) for reg in B.regions(car, r))
        row["reach%d" % r] = d
        row["reachfrac%d" % r] = d / (2.0 * k)
    # max over contiguous regions of size r is MONOTONE NON-DECREASING in r (a size-(r+1)
    # region contains a size-r one), so both thresholds are found by bisection -- exact.
    def reach(r):
        return max(B.local_logical_dim(car, reg) for reg in B.regions(car, r))
    def first_r(pred):
        lo, hi = 1, n
        if not pred(reach(hi)): return None
        while lo < hi:
            mid = (lo + hi) // 2
            if pred(reach(mid)): hi = mid
            else: lo = mid + 1
        return lo
    rstar = first_r(lambda d: d >= 2 * k)
    row["rstar"] = rstar
    row["rstar_over_n"] = (rstar / n) if rstar else None
    prot = first_r(lambda d: d > 0)
    row["protection_radius"] = prot
    # ---- BASIS-INVARIANT: capacity
    row["joint_block_dim"] = 2 ** (n - len(car["stabs"]) - k)
    row["k_over_log2dim"] = k / float(n)
    return row

COLS_A = ["label", "k", "log2dim", "wbar", "wbar_over_n",
          "ov_edgefrac", "ov_ncomp", "ov_giant", "ov_clus", "ov_lapgap", "ov_lammax_over_tr",
          "ov_percolates"]
COLS_B = ["label", "k", "cf1_edgefrac", "cf1_maxset", "cf1_maxsingle", "cf2_edgefrac", "cf2_ncomp", "cf2_giant",
          "cf2_clus", "cf2_lapgap", "cf2_percolates", "cf2_maxset", "cf2_maxsingle",
          "cf3_maxset", "cf3_maxsingle"]
COLS_C = ["label", "k", "protection_radius", "reach1", "reach2", "reach3", "reach4",
          "reachfrac2", "reachfrac4", "rstar", "rstar_over_n", "joint_block_dim",
          "k_over_log2dim"]

def table(rows, cols, title, note):
    P(); P(title); P(note)
    hdr = "  " + "".join("%-16s" % c for c in cols)
    P(hdr); P("  " + "-" * (16 * len(cols) - 2))
    for r in rows:
        cells = []
        for c in cols:
            v = r.get(c)
            if isinstance(v, float): cells.append("%-16.4f" % v)
            elif isinstance(v, bool): cells.append("%-16s" % v)
            else: cells.append("%-16s" % v)
        P("  " + "".join(cells))

t0 = time.time()
P("=" * 170)
P("S2  ORDER-PARAMETER BATTERY vs NUMBER OF RECORDS k -- exact, two families in one table")
P("=" * 170)

NA = [4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32]
MB = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16]
MC = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16]
rows = []
for n in NA:
    r = measure(C.family_A(n)); rows.append(r)
    P("... A n=%d done  t=%.1fs" % (n, time.time() - t0))
for m in MB:
    r = measure(C.family_B(m)); rows.append(r)
    P("... B m=%d done  t=%.1fs" % (m, time.time() - t0))
for m in MC:
    r = measure(C.family_C(m)); rows.append(r)
    P("... C m=%d done  t=%.1fs" % (m, time.time() - t0))

rows = [r for r in rows if r]
inter = []
for kk in sorted(set(r["k"] for r in rows)):
    for tag in ("A", "B", "C"):
        for r in rows:
            if r["k"] == kk and r["fam"] == tag: inter.append(r)

table(inter, COLS_A, "TABLE 1  SUPPORT-OVERLAP GRAPH  (BASIS-DEPENDENT)",
      "  families A, B and C are interleaved at EQUAL k -- each is the other's control")
table(inter, COLS_B, "TABLE 2  CO-FLIP GRAPH: records ONE admissible region-r operation flips TOGETHER (BASIS-DEPENDENT)",
      "  cf1_* is the clause-(v) claim and MUST be zero; cf2_*/cf3_* are the POSITIVE CONTROLS in the same table")
table(inter, COLS_C, "TABLE 3  WHAT ONE CONTIGUOUS REGION CAN WRITE  (BASIS-INVARIANT)",
      "  reach_r = dim_F2 of the logical group a single region of size r can write; 2k = dim of the whole record group")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s2_combinatorial_battery.txt",
     "w").write("\n".join(OUT) + "\n")
json.dump(rows, open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s2_rows.json", "w"), indent=1)
P("total %.1fs" % (time.time() - t0))
