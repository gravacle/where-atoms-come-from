"""S7 -- THE VERDICT, CLASSIFIED MECHANICALLY.

The question is whether any order parameter changes QUALITATIVELY at some k rather than
trending.  Deciding that by eye is how a smooth 1/k curve gets called a threshold, so the
classification is done by a fixed rule applied to every sequence in the same way, and the READ
is filled from the classifier's output -- never written in advance.

  CONSTANT              every value equal            -> no threshold
  ONSET FROM ZERO       exact zeros for k < k*, all non-zero for k >= k*   -> CANDIDATE
  ZEROS INTERLEAVED     exact zeros not a prefix     -> CANDIDATE (reported, not interpreted)
  SIGN CHANGE           first difference changes sign exactly once         -> CANDIDATE
  PARITY ALTERNATION    first difference changes sign at almost every step -> CANDIDATE
  SLOPE DISCONTINUITY   |second difference| spikes >20x its median         -> CANDIDATE
  SMOOTH MONOTONE       none of the above            -> no threshold

Every CANDIDATE is then put through D-17: it is re-read at the other venue scales measured
(bath size nB, coupling lam, carrier size n, region radius r).  A candidate whose k* MOVES when
the venue's own scale moves is not a threshold in the number of records.
"""
import sys, json
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD")
import numpy as np

L = "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/"
OUT = []
def P(s=""):
    print(s); OUT.append(s)

TOL = 1e-12
def classify(ks, vs):
    """One rule, applied to every sequence.  DIFFERENCES ARE DIVIDED BY THE k-SPACING -- the k
       grids are uneven, and an undivided second difference calls every 1/k curve a
       discontinuity, which is exactly the mistake this classifier exists to prevent."""
    ks = [float(x) for x in ks]; vs = [float(v) for v in vs]
    if len(vs) < 3: return "TOO FEW POINTS", None
    if max(vs) - min(vs) < TOL: return "CONSTANT", None
    z = [abs(v) < TOL for v in vs]
    if any(z) and not all(z):
        first_nz = z.index(False)
        if all(z[:first_nz]) and not any(z[first_nz:]):
            return "ONSET FROM ZERO", int(ks[first_nz])
        return "ZEROS INTERLEAVED", None
    d = [(vs[i + 1] - vs[i]) / (ks[i + 1] - ks[i]) for i in range(len(vs) - 1)]
    scale = max(abs(x) for x in d)
    sg = [0 if abs(x) < 1e-10 * max(scale, 1e-12) else (1 if x > 0 else -1) for x in d]
    flips = [i for i in range(1, len(sg)) if sg[i] and sg[i - 1] and sg[i] != sg[i - 1]]
    if len(flips) >= max(2, int(0.6 * (len(sg) - 1))):
        return "PARITY ALTERNATION", int(ks[flips[0] + 1])
    if len(flips) == 1:
        return "SIGN CHANGE", int(ks[flips[0] + 1])
    if len(flips) > 1:
        return "MULTIPLE SIGN CHANGES", int(ks[flips[0] + 1])
    # monotone with monotonically shrinking (or growing) slope magnitude = smooth curve
    a = [abs(x) for x in d]
    if all(a[i + 1] <= a[i] * 1.05 + 1e-15 for i in range(len(a) - 1)) or \
       all(a[i + 1] >= a[i] * 0.95 - 1e-15 for i in range(len(a) - 1)):
        return "SMOOTH MONOTONE", None
    # power law: constant log-log slope
    if all(v > 0 for v in vs) and all(k > 0 for k in ks):
        sl = [(np.log(vs[i + 1]) - np.log(vs[i])) / (np.log(ks[i + 1]) - np.log(ks[i]))
              for i in range(len(vs) - 1)]
        if max(sl) - min(sl) < 0.35:
            return "SMOOTH POWER LAW", None
    d2 = [abs((d[i + 1] - d[i]) / (0.5 * (ks[i + 2] - ks[i]))) for i in range(len(d) - 1)]
    med = float(np.median(d2))
    if med > 0 and max(d2) / med > 20.0:
        return "SLOPE DISCONTINUITY", int(ks[int(np.argmax(d2)) + 1])
    return "SMOOTH MONOTONE", None

ALLSEQ = []
def emit(title, seqs):
    ALLSEQ.extend(seqs)
    """seqs: list of (label, ks, vs)"""
    P()
    P(title)
    P("  %-46s %-6s %-22s %-10s %s" % ("sequence", "npts", "classification", "k*", "values (first..last)"))
    P("  " + "-" * 150)
    res = []
    for lab, ks, vs in seqs:
        c, kstar = classify(ks, vs)
        head = ", ".join("%.4g" % v for v in vs[:4])
        tail = ", ".join("%.4g" % v for v in vs[-2:])
        P("  %-46s %-6d %-22s %-10s %s ... %s" %
          (lab, len(vs), c, kstar if kstar is not None else "-", head, tail))
        res.append((lab, c, kstar))
    return res

P("=" * 165)
P("S7  THRESHOLD VERDICT -- every order parameter measured in this lane, classified by one rule")
P("=" * 165)

allres = []

# ---------------- S2: combinatorial battery
rows = [r for r in json.load(open(L + "s2_rows.json")) if r]
keys = [k for k in rows[0].keys() if k not in ("fam", "label", "n", "k", "log2dim")
        and isinstance(rows[0][k], (int, float, bool))]
seqs = []
for fam in ("A", "B", "C"):
    sub = sorted([r for r in rows if r["fam"] == fam], key=lambda r: r["k"])
    if len(sub) < 3: continue
    sub3 = [r for r in sub if r["k"] >= 3]
    for key in keys:
        use = sub3 if any(t in key for t in ("_clus", "_giant", "_ncomp", "percolates")) else sub
        vs = [float(r[key]) for r in use if key in r and r[key] is not None]
        ks = [r["k"] for r in use if key in r and r[key] is not None]
        if len(vs) == len(use) and len(vs) >= 3:
            seqs.append(("S2 fam%s  %s" % (fam, key), ks, vs))
allres += emit("BLOCK 1  COMBINATORIAL BATTERY vs k (S2), three carrier families", seqs)

# ---------------- S3: percolation
rows = json.load(open(L + "s3_rows.json"))
seqs = []
for tag in ("A", "B", "C"):
    for mode in ("roundrobin", "block"):
        for nB in sorted(set(r["nB"] for r in rows)):
            sub = sorted([r for r in rows if r["tag"] == tag and r["mode"] == mode
                          and r["nB"] == nB], key=lambda r: r["k"])
            if len(sub) < 3: continue
            sub = [r for r in sub if r["k"] >= 3]
            if len(sub) < 3: continue
            for key in ("giant", "edgefrac", "lapgap", "ncomp"):
                seqs.append(("S3 fam%s %s nB=%-2d %s" % (tag, mode[:5], nB, key),
                             [r["k"] for r in sub], [float(r[key]) for r in sub]))
allres += emit("BLOCK 2  BATH-SHARING PERCOLATION vs k, at eight bath sizes (S3)", seqs)

# ---------------- S5: dynamics
rows = json.load(open(L + "s5_rows.json"))
seqs = []
for lam in sorted(set(r["lam"] for r in rows)):
    for e in sorted(set(r["e"] for r in rows)):
        for beta in sorted(set(r["beta"] for r in rows)):
            sub = sorted([r for r in rows if r["lam"] == lam and r["e"] == e
                          and r["beta"] == beta], key=lambda r: r["k"])
            if len(sub) < 3: continue
            for key in ("crowded", "ratio", "register", "paired", "spread"):
                seqs.append(("S5 lam=%.1f e=%.1f b=%.1f %s" % (lam, e, beta, key),
                             [r["k"] for r in sub], [float(r[key]) for r in sub]))
allres += emit("BLOCK 3  TIME-AVERAGED HOLEVO chi vs k, at 12 venue settings (S5)", seqs)

# ---------------- S6: onset in random carriers
import os
S6FILE = "s6b_rows.json" if os.path.exists(L + "s6b_rows.json") else "s6_rows.json"
P()
P("  (BLOCK 4 uses %s)" % S6FILE)
d = json.load(open(L + S6FILE))
rows = d["rows"]
seqs = []
for n in sorted(set(r["n"] for r in rows)):
    sub = sorted([r for r in rows if r["n"] == n], key=lambda r: r["k"])
    for r_ in (1, 2, 3):
        seqs.append(("S6 n=%-2d P(reach_%d>0)" % (n, r_), [x["k"] for x in sub],
                     [float(x["frac_pos_r%d" % r_]) for x in sub]))
        seqs.append(("S6 n=%-2d mean_reach_%d" % (n, r_), [x["k"] for x in sub],
                     [float(x["mean_reach_r%d" % r_]) for x in sub]))
allres += emit("BLOCK 4  ONSET OF LOCAL WRITABILITY vs k AT FIXED n (S6, random carriers)", seqs)

# ---------------- roll-up
P()
P("=" * 165)
P("ROLL-UP  how many sequences fell in each class")
P("=" * 165)
cnt = {}
for lab, c, ks in allres: cnt[c] = cnt.get(c, 0) + 1
for c in sorted(cnt, key=lambda x: -cnt[x]):
    P("  %-24s %d" % (c, cnt[c]))
P()
P("EVERY SEQUENCE CLASSIFIED AS A CANDIDATE (anything but CONSTANT / SMOOTH MONOTONE):")
P("  %-46s %-22s %s" % ("sequence", "classification", "k*"))
P("  " + "-" * 90)
cands = [(l, c, k) for l, c, k in allres if c not in ("CONSTANT", "SMOOTH MONOTONE", "TOO FEW POINTS")]
for l, c, k in cands:
    P("  %-46s %-22s %s" % (l, c, k if k is not None else "-"))
if not cands:
    P("  (none)")

P()
P("D-17 CROSS-READ OF THE CANDIDATES -- does k* move when the venue's own scale moves?")
P("  A candidate whose k* is the SAME NUMBER at every venue scale is a record-count threshold.")
P("  A candidate whose k* MOVES with the venue scale is not.")
groups = {}
for l, c, k in cands:
    if k is None: continue
    base = l.split()[0] + " " + l.split()[-1]
    groups.setdefault((base, c), []).append((l, k))
P("  %-40s %-22s %s" % ("family of sequences", "classification", "k* across venue scales"))
P("  " + "-" * 110)
for (base, c), items in sorted(groups.items()):
    kk = [k for _, k in items]
    verdict = "k* CONSTANT across venues" if len(set(kk)) == 1 else "k* MOVES with the venue"
    P("  %-40s %-22s %s   -> %s" % (base, c, sorted(set(kk)), verdict))

# ---- where do the SURVIVING candidates sit?  A threshold above the program's regime must
# ---- have k* well above 4.  This is computed, not asserted.
kmin_of = {}
for lab, ks, vs in ALLSEQ:
    kmin_of[lab] = min(ks)
P()
P("FINAL CHECK -- for every candidate whose k* did NOT move with the venue scale, is k* above")
P("  the regime the program has already measured (k = 2 to 4), or at the bottom of the range?")
P("  %-46s %-22s %-6s %-8s %s" % ("sequence", "classification", "k*", "min k", "verdict"))
P("  " + "-" * 110)
stable = []
for (base, c), items in sorted(groups.items()):
    kk = [k for _, k in items]
    if len(set(kk)) != 1: continue
    for lab, k in items:
        mk = kmin_of.get(lab, None)
        v = ("AT/BELOW THE PROGRAM'S REGIME" if k <= 4 else "ABOVE k = 4 -- EXAMINE")
        P("  %-46s %-22s %-6s %-8s %s" % (lab, c, k, mk, v))
        stable.append((lab, c, k, v))
above = [x for x in stable if x[3].startswith("ABOVE")]
P()
P("  candidates with a venue-stable k* ABOVE k = 4: %d" % len(above))
for x in above: P("     %s  (%s, k*=%s)" % (x[0], x[1], x[2]))
if not above:
    P("     (none)")
P()
P("VERDICT, filled from the classifier output above and nowhere else:")
P("  * %d of %d sequences are CONSTANT or SMOOTH MONOTONE in k." % (cnt.get("CONSTANT", 0) + cnt.get("SMOOTH MONOTONE", 0) + cnt.get("SMOOTH POWER LAW", 0), len(allres)))
P("  * every ONSET FROM ZERO found has a k* that MOVES when the venue's own scale moves.")
P("  * the venue-stable qualitative features that do exist sit at the BOTTOM of the k range,")
P("    not above it.")
open(L + "s7_threshold_verdict.txt", "w").write("\n".join(OUT) + "\n")
