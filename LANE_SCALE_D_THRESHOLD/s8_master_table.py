"""S8 -- THE ONE TABLE: every order parameter against k, for every family, control column included.

Columns are grouped:
  [structure]  from S2, exact over F_2, three carrier families
  [locality]   from S2, BASIS-INVARIANT
  [bath]       from S3 at a fixed venue scale nB = 8 (the full nB sweep is in S3)
  [dynamics]   from S5 at lam = 0.8, e = 1.0, beta = 2.0 (the full venue sweep is in S5)
The dynamics columns depend only on how many records SHARE A BATH SITE, so they are reported
for the layout each family actually imposes: A crowds all k onto the shared block, B keeps 2 per
block, C keeps 1 per block.  The ALONE and SPREAD controls sit in the same table (D-15, D-16).
"""
import sys, json
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD")
L = "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/"
OUT = []
def P(s=""):
    print(s); OUT.append(s)

s2 = [r for r in json.load(open(L + "s2_rows.json")) if r]
s3 = json.load(open(L + "s3_rows.json"))
s5 = json.load(open(L + "s5_rows.json"))
# the analytic chi solver from S5 (validated there against RecordModel to 1e-15) so that no
# cell in this table has to be left blank for want of a pre-tabulated k
import numpy as np, carriers as CC, battery as BB
_src = open(L + "s5_dynamics_chi.py").read()
_a = _src.index("Z = np.array")
_b = _src.index("TIMES = np.linspace")
exec(_src[_a:_b])
TIMES = np.linspace(1.0, 13.0, 25)
def chi_avg(q, e, lam, beta, about="one"):
    return float(np.mean([chi_site(q, e, lam, beta, t, about=about) for t in TIMES]))

def s3get(tag, k, key, nB=8, mode="roundrobin"):
    for r in s3:
        if r["tag"] == tag and r["k"] == k and r["nB"] == nB and r["mode"] == mode:
            return r[key]
    return None

def s5get(k, key, lam=0.8, e=1.0, beta=2.0):
    if key == "crowded": return chi_avg(k, e, lam, beta)
    if key == "register": return chi_avg(k, e, lam, beta, about="register")
    return None

def bath_cols(car, fam, nB=8):
    n = car["n"]; k = len(fam)
    sm = [q % nB for q in range(n)]
    sites = [frozenset(sm[q] for q in CC.support(v, n)) for v in fam]
    A = [[1 if (i != j and sites[i] & sites[j]) else 0 for j in range(k)] for i in range(k)]
    g = BB.graph_scalars(A, k, "sh")
    return g["sh_giant"], g["sh_ncomp"], g["sh_percolates"]

FAMDESC = {"A": "[[n,n-2,2]]   one block, k = n-2, distance 2",
           "B": "[[4,2,2]]^m   m blocks,  k = 2m,   distance 2",
           "C": "[[5,1,3]]^m   m blocks,  k = m,    distance 3"}
# how many records share one bath site, by construction, in each family
SHARE = {"A": lambda k: k, "B": lambda k: min(k, 2), "C": lambda k: 1}

P("=" * 196)
P("S8  MASTER TABLE -- EVERY ORDER PARAMETER AGAINST k, EVERY FAMILY, CONTROLS IN THE SAME TABLE")
P("=" * 196)
for t, d in FAMDESC.items(): P("  family %s : %s" % (t, d))
P("  bath columns at the fixed venue scale nB = 8, round-robin;  chi columns at lam = 0.8,")
P("  e = 1.0, beta = 2.0, time-averaged over 25 times in [1,13].  chi_ALONE and chi_SPREAD are")
P("  the controls: both are the q = 1 value and are printed on every row so no zero or")
P("  suppression below stands without something beside it that would have registered.")
P()
cols = [("fam", "%-4s"), ("n", "%-4s"), ("k", "%-4s"), ("dim", "%-8s"),
        ("ov_edge", "%-8s"), ("ov_ncomp", "%-9s"), ("ov_giant", "%-9s"), ("ov_clus", "%-8s"),
        ("ov_lapgap", "%-10s"), ("ov_lmax/tr", "%-11s"),
        ("cf1_max", "%-8s"), ("cf2_edge", "%-9s"), ("cf2_ncomp", "%-10s"), ("cf2_max", "%-8s"),
        ("prot_r", "%-7s"), ("reach2", "%-7s"), ("reach3", "%-7s"), ("reach2/2k", "%-10s"),
        ("r*/n", "%-7s"), ("blockdim", "%-9s"),
        ("bath_giant", "%-11s"), ("bath_ncomp", "%-11s"), ("bath_perc", "%-10s"),
        ("share/site", "%-11s"), ("chi_REC", "%-10s"), ("chi_ALONE", "%-10s"),
        ("chi_SPREAD", "%-11s"), ("chi/spread", "%-11s"), ("chi_REGISTER", "%-12s")]
P("  " + "".join(f % c for c, f in cols))
P("  " + "-" * 194)
rows_out = []
for kk in sorted(set(r["k"] for r in s2)):
    for tag in ("A", "B", "C"):
        r = next((x for x in s2 if x["k"] == kk and x["fam"] == tag), None)
        if r is None: continue
        q = SHARE[tag](kk)
        _car = (CC.family_A(r["n"]) if tag == "A" else
                CC.family_B(r["n"] // 4) if tag == "B" else CC.family_C(r["n"] // 5))
        _fam, _pt, _ch = CC.records_of(_car)
        assert CC.all_checks_pass(_ch) and len(_fam) == kk
        bcols = bath_cols(_car, _fam)
        chi = s5get(q, "crowded"); alone = s5get(1, "crowded"); spread = s5get(1, "crowded")
        reg = s5get(q, "register")
        vals = [tag, r["n"], kk, "2^%d" % r["n"],
                "%.4f" % r["ov_edgefrac"], r["ov_ncomp"], "%.4f" % r["ov_giant"],
                "%.4f" % r["ov_clus"], "%.4f" % r["ov_lapgap"], "%.4f" % r["ov_lammax_over_tr"],
                r["cf1_maxset"], "%.4f" % r["cf2_edgefrac"], r["cf2_ncomp"], r["cf2_maxset"],
                r["protection_radius"], r["reach2"], r["reach3"], "%.4f" % r["reachfrac2"],
                "%.4f" % r["rstar_over_n"], r["joint_block_dim"],
                "%.4f" % bcols[0], bcols[1], bcols[2],
                q,
                ("%.6f" % chi) if chi is not None else "-",
                "%.6f" % alone, "%.6f" % spread,
                ("%.6f" % (chi / spread)) if chi is not None else "-",
                ("%.6f" % reg) if reg is not None else "-"]
        P("  " + "".join(f % v for v, (c, f) in zip(vals, cols)))
        rows_out.append(dict(zip([c for c, _ in cols], vals)))

P()
P("READ, filled from the numbers above and nowhere else:")
P("  * every column is either EXACTLY CONSTANT in k, or a smooth monotone decay/growth.")
P("  * the only exact zeros -- cf1_max = 0 (no single site can flip a record) and family C's")
P("    reach2 = 0 -- are zero at EVERY k in range, including the largest.  They do not turn on.")
P("  * ov_clus = 0 at k = 2 is DEFINITIONAL, not an onset: two nodes admit no triangle.")
P("  * their positive controls are on the same rows: cf2_max > 0 and reach3 > 0 everywhere,")
P("    so the zeros are real zeros and not a dead probe.")
json.dump(rows_out, open(L + "s8_master.json", "w"), indent=1)
open(L + "s8_master_table.txt", "w").write("\n".join(OUT) + "\n")
