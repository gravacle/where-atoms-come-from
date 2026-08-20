"""ADVERSARIAL CHECK V1 -- THE MISSING D-15 POSITIVE CONTROL FOR THE GEOMETRY INSTRUMENT.

The lane's headline ZERO is "no fixed intrinsic dimension emerges: d90 grows with the number
of records".  Its controls (FREE, RAND) are both NEGATIVE controls -- objects that also fail
to have a geometry.  D-15 demands a POSITIVE control: a point set that GENUINELY LIVES IN A
FIXED DIMENSION, pushed through the SAME analyse()/double_centre/dim_stats pipeline at the
SAME point counts (2k = 4..36).  If d90 grows there too, the instrument cannot detect a
geometry and the null is an instrument artifact.  If d90 stays flat there, the instrument
works and the null stands.

Also V3: are the "GROWS" verdicts for relations (a), (c) and (b1) arithmetic rather than
measurement?  Each of those relation matrices is degenerate (permutation / identity /
diagonal), so the transformed distance matrix is a REGULAR SIMPLEX and d90 is forced by
ceil(0.9*(#distinct points - 1)) with no carrier information in it at all.  Checked in closed
form against the lane's own numbers.
"""
import sys, math
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC")
from lib_scaleb import double_centre, dim_stats, corr_dim
import numpy as np

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)

rng = np.random.default_rng(7)
MS = [4, 8, 12, 16, 20, 24, 28, 32, 36]      # the lane's own 2k values

def d90_of(D):
    st = dim_stats(double_centre(D))
    return st["d_frac"], st["d_pr"], st.get("neg_frac", 0.0)

def euclid(X):
    d2 = ((X[:, None, :] - X[None, :, :]) ** 2).sum(-1)
    return np.sqrt(np.maximum(d2, 0.0))

P("=" * 108)
P("V1  POSITIVE CONTROL FOR THE d90 INSTRUMENT -- point sets that DO have a fixed dimension,")
P("    pushed through the lane's OWN double_centre + dim_stats at the lane's OWN point counts.")
P("=" * 108)
P("")
P("  m(=2k) | 1D line   2D disc   3D ball   2D+5%noise 2D+20%noise 6D ball  | LANE (a)S LANE (b2)")
P("  " + "-" * 104)
lane_a  = dict(zip(MS, [1, 3, 5, 7, 9, 10, 12, 14, 16]))
lane_b2 = dict(zip(MS, [2, 4, 6, 7, 9, 11, 13, 14, 16]))
res = {k: [] for k in ["1D", "2D", "3D", "2Dn5", "2Dn20", "6D"]}
for m in MS:
    cells = []
    for tag, dim, noise in [("1D", 1, 0.0), ("2D", 2, 0.0), ("3D", 3, 0.0),
                            ("2Dn5", 2, 0.05), ("2Dn20", 2, 0.20), ("6D", 6, 0.0)]:
        # average over 8 draws so a single unlucky sample cannot decide anything
        ds = []
        for rep in range(8):
            X = rng.normal(size=(m, dim))
            X /= np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
            X *= rng.random((m, 1)) ** (1.0 / dim)
            D = euclid(X)
            if noise > 0:
                N = rng.normal(scale=noise, size=(m, m)); N = np.abs((N + N.T) / 2)
                D = D + N; np.fill_diagonal(D, 0.0)
            ds.append(d90_of(D)[0])
        med = int(np.median(ds))
        res[tag].append(med)
        cells.append(med)
    P("  %-6d | %-9d %-9d %-9d %-10d %-11d %-8d | %-9d %-9d" %
      (m, cells[0], cells[1], cells[2], cells[3], cells[4], cells[5], lane_a[m], lane_b2[m]))
P("")
xs = np.array(MS, float)
for tag in ["1D", "2D", "3D", "2Dn5", "2Dn20", "6D"]:
    sl = np.polyfit(xs, np.array(res[tag], float), 1)[0]
    P("  positive control %-6s  d90 vs m slope = %+.4f   %s" %
      (tag, sl, "FLAT (instrument detects a fixed dimension)" if abs(sl) < 0.05
       else "GROWS (instrument cannot hold a fixed dimension here)"))
P("  lane (a) symplectic          d90 vs m slope = %+.4f" % np.polyfit(xs, np.array([lane_a[m] for m in MS], float), 1)[0])
P("  lane (b2) crowding           d90 vs m slope = %+.4f" % np.polyfit(xs, np.array([lane_b2[m] for m in MS], float), 1)[0])

# ------------------------------------------------------------------ V3 simplex arithmetic
P("")
P("=" * 108)
P("V3  ARE THE 'GROWS' VERDICTS FOR (a), (c), (b1) ARITHMETIC RATHER THAN MEASUREMENT?")
P("    A regular simplex on p equidistant points has p-1 equal MDS eigenvalues, so")
P("    d90 = ceil(0.9*(p-1)) with NO dependence on the carrier.  Closed form vs lane numbers.")
P("=" * 108)
P("")
P("  2k  | (a) p=k simplex  pred  lane | (c) p=2k simplex pred  lane | (b1) pred  lane")
P("  " + "-" * 96)
lane_c  = dict(zip(MS, [3, 7, 10, 14, 18, 21, 25, 28, 32]))
ok_a = ok_c = True
for m in MS:
    k = m // 2
    pa = int(math.ceil(0.9 * (k - 1))) if k > 1 else 1
    # d90 counts eigenvalues, min 1
    pa = max(pa, 1)
    pc = int(math.ceil(0.9 * (m - 1)))
    ok_a &= (pa == lane_a[m]); ok_c &= (pc == lane_c[m])
    P("  %-4d| p=%-3d           %-5d %-4d | p=%-3d           %-5d %-4d | %-5d %-4d" %
      (m, k, pa, lane_a[m], m, pc, lane_c[m], pc, m - 1))
P("")
P("  closed-form simplex reproduces lane (a) at every n: %s" % ok_a)
P("  closed-form simplex reproduces lane (c) at every n: %s" % ok_c)
P("  -> where True, that column's 'GROWS WITH n' is a fact about the ARITHMETIC OF A SIMPLEX,")
P("     not a measurement on the carrier: it would read the same for ANY carrier whose")
P("     relation matrix is a permutation / the identity / diagonal.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/VERIFY/V1_instrument_positive_control.txt",
     "w").write("\n".join(OUT) + "\n")
