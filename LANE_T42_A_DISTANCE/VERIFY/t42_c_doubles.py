"""T42-A part 3: d_W BETWEEN CONFIGURATIONS on quantum doubles with |G| = 4 = 2^2.  EXACT.

Two non-isomorphic groups of order 4, both powers of two:

  D(Z_2 x Z_2)  = two toric-code factors on the SAME lattice edges (one 4-dim qudit per
                  edge).  Weight counts QUDITS touched, so the two factors can SHARE support:
                  the naive product formula d = d1 + d2 must fail, and does (D-15 control).
                  Carriers L = 2 and 3; at L = 2 the full (x|z)-per-factor space is also
                  searched to re-verify the z=0 reduction on this carrier.

  D(Z_4)        = the Z_4 toric code on the 2 x 2 torus, oriented lattice, arithmetic mod 4
                  throughout (exact integers).  Probes TORSION: is d(2t) two steps or one?

All configuration distances are per-pair searched, all metric axioms checked over ALL pairs
and ALL triples, venue Aut groups enumerated exactly (16 points), every verdict gated.
"""
import sys, time
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_T42_A_DISTANCE/VERIFY"
sys.path.insert(0, LANE)
from t42_lib import (pc, solve_affine_f2, span_all, vec_to_mask, metric_axioms, aut_count,
                     toric_stabilizers, weight_xz)
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import symplectic_logicals

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s); print(s); sys.stdout.flush()

def gate(label, ok, detail=""):
    say(("PASS  " if ok else "FAIL  ") + label + (("  " + detail) if detail else ""))
    return ok

say("T42-A part 3  --  quantum doubles |G|=4, exact, started " + time.strftime("%F %T"))
say("")
ALL_OK = True

# ================================================================= D(Z_2 x Z_2)
def toric_cosets(L):
    """For each 2-bit class t: (x0, nullbasis) of the z=0 writer space, plus n."""
    n, stab, stars, plaqs = toric_stabilizers(L)
    pairs = symplectic_logicals(stab, n)
    logi = [(vec_to_mask(X), vec_to_mask(Z)) for X, Z in pairs]
    Zlogs = [logi[0][1], logi[1][1]]
    out = {}
    for t in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        rows = [p for p in plaqs] + [Zlogs[0] >> n, Zlogs[1] >> n]
        rhs = [0] * len(plaqs) + [t[0], t[1]]
        x0, nb, _ = solve_affine_f2(rows, rhs, n)
        out[t] = (x0, nb)
    return n, out

Z22_METRICS = {}
for L in (2, 3):
    t0 = time.time()
    n, cosets = toric_cosets(L)
    elems = {t: [x0 ^ v for v in span_all(nb)] for t, (x0, nb) in cosets.items()}
    labels4 = [(a, b) for a in range(2) for b in range(2)]
    configs = [(t1, t2) for t1 in labels4 for t2 in labels4]
    dmat = {}
    class_min = {}
    searched = 0
    for s in configs:
        for sp_ in configs:
            t1 = (s[0][0] ^ sp_[0][0], s[0][1] ^ sp_[0][1])
            t2 = (s[1][0] ^ sp_[1][0], s[1][1] ^ sp_[1][1])
            key = (t1, t2)
            if key not in class_min:
                best = min(pc(x1 | x2) for x1 in elems[t1] for x2 in elems[t2])
                class_min[key] = best
                searched += len(elems[t1]) * len(elems[t2])
            dmat[(s, sp_)] = class_min[key]
    say("D(Z2xZ2) L=%d  16 configurations, %d qudits; writer pairs searched: %d"
        % (L, n, searched))
    rep = metric_axioms(configs, dmat)
    ALL_OK &= gate("D(Z2xZ2) L=%d METRIC AXIOMS over ALL %d pairs, ALL %d triples"
                   % (L, rep["n_pairs"], rep["n_triples"]), rep["all"],
                   "" if rep["all"] else str(rep))
    # structure: d = L * popcount(t1 OR t2)  (joint-support crossing cost)
    or_ok = all(class_min[(t1, t2)] == L * pc(vec_to_mask(t1) | vec_to_mask(t2))
                for (t1, t2) in class_min)
    ALL_OK &= gate("D(Z2xZ2) L=%d d_W = L x |support(t1 OR t2)| -- crossing cost of the "
                   "JOINT homology support" % L, or_ok,
                   str({k: v for k, v in sorted(class_min.items())}) if not or_ok else "")
    # D-15 control: the imported product form d1+d2 must FAIL
    prod_fail = [k for k in class_min
                 if class_min[k] != class_min[(k[0], (0, 0))] + class_min[((0, 0), k[1])]]
    ALL_OK &= gate("D(Z2xZ2) L=%d D-15 CONTROL: naive product formula d1+d2 FAILS "
                   "(factors share qudit support)" % L, len(prod_fail) > 0,
                   "witness class %s: d=%d, d1+d2=%d" % (prod_fail[0],
                    class_min[prod_fail[0]],
                    class_min[(prod_fail[0][0], (0, 0))] + class_min[((0, 0), prod_fail[0][1])])
                   if prod_fail else "")
    kind, na = aut_count(configs, dmat)
    say("D(Z2xZ2) L=%d  D-22 venue Aut of the 16-point metric: |Aut| = %d (%s)" % (L, na, kind))
    ALL_OK &= gate("D(Z2xZ2) L=%d venue has geometry (Aut != S_16)" % L,
                   na != 20922789888000, "|Aut|=%d" % na)
    Z22_METRICS[L] = (configs, dmat, class_min)
    say("D(Z2xZ2) L=%d done in %.1fs" % (L, time.time() - t0))
    say("")

# z-freedom control at L=2 for the double: allow z-parts in BOTH factors
say("REDUCTION CONTROL for D(Z2xZ2) at L=2: full (x|z) search in both factors")
n, stab, stars, plaqs = toric_stabilizers(2)
pairs = symplectic_logicals(stab, n)
logi = [(vec_to_mask(X), vec_to_mask(Z)) for X, Z in pairs]
Zlogs = [logi[0][1], logi[1][1]]
mask_n = (1 << n) - 1
def full_row(g):
    return (g >> n) | ((g & mask_n) << n)
stab_masks = [vec_to_mask(sv) for sv in stab]
rows_full = [full_row(g) for g in stab_masks]
full_elems = {}
for t in [(0, 0), (1, 0), (0, 1), (1, 1)]:
    R = rows_full + [full_row(Zlogs[0]), full_row(Zlogs[1])]
    rhs = [0] * len(rows_full) + [t[0], t[1]]
    u0, nb, _ = solve_affine_f2(R, rhs, 2 * n)
    full_elems[t] = [u0 ^ v for v in span_all(nb)]
def touched(u):
    return (u & mask_n) | (u >> n)
ok_all = True
_, _, cm2 = Z22_METRICS[2]
for (t1, t2), want in sorted(cm2.items()):
    best = min(pc(touched(u1) | touched(u2))
               for u1 in full_elems[t1] for u2 in full_elems[t2])
    ok_all &= (best == want)
ALL_OK &= gate("D(Z2xZ2) L=2 allowing z-parts in both factors never lowers d_W "
               "(exhaustive, 2^8 x 2^8 per class pair)", ok_all)
say("")

# ================================================================= D(Z_4), 2x2 torus
say("D(Z_4) on the 2x2 torus, arithmetic mod 4, oriented lattice")
L = 4  # modulus, not lattice size -- lattice is 2x2 below
LL = 2
nE = 2 * LL * LL  # 8 edges
def h(i, j):
    return (i % LL) * LL + (j % LL)
def v(i, j):
    return LL * LL + (i % LL) * LL + (j % LL)
# oriented plaquette rows: circulation  +h(i,j) +v(i,j+1) -h(i+1,j) -v(i,j)
plaq_rows = []
for i in range(LL):
    for j in range(LL):
        r = [0] * nE
        r[h(i, j)] += 1
        r[v(i, j + 1)] += 1
        r[h(i + 1, j)] -= 1
        r[v(i, j)] -= 1
        plaq_rows.append([x % 4 for x in r])
# oriented star vectors: +h(i,j) -h(i,j-1) +v(i,j) -v(i-1,j)
star_vecs = []
for i in range(LL):
    for j in range(LL):
        r = [0] * nE
        r[h(i, j)] += 1
        r[h(i, j - 1)] -= 1
        r[v(i, j)] += 1
        r[v(i - 1, j)] -= 1
        star_vecs.append(tuple(x % 4 for x in r))

def dot4(a, b):
    return sum(x * y for x, y in zip(a, b)) % 4

# C = all x-strings commuting with every plaquette (brute force 4^8 = 65536, exact)
C = []
for m in range(4 ** nE):
    u = []
    mm = m
    for _ in range(nE):
        u.append(mm % 4)
        mm //= 4
    u = tuple(u)
    if all(dot4(u, p) == 0 for p in plaq_rows):
        C.append(u)
ALL_OK &= gate("D(Z4) |C| = 4^5 (kernel of oriented circulation, computed) ", len(C) == 4 ** 5,
               "|C|=%d" % len(C))
stars_ok = all(all(dot4(s, p) == 0 for p in plaq_rows) for s in star_vecs)
ALL_OK &= gate("D(Z4) every oriented star lies in C (stabilizers commute)", stars_ok)
# span of stars mod 4
span = {tuple([0] * nE)}
frontier = [tuple([0] * nE)]
while frontier:
    nxt = []
    for x in frontier:
        for s in star_vecs:
            y = tuple((a + b) % 4 for a, b in zip(x, s))
            if y not in span:
                span.add(y)
                nxt.append(y)
    frontier = nxt
ALL_OK &= gate("D(Z4) |span stars| = 4^3 (one relation: sum of stars = 0)", len(span) == 4 ** 3,
               "|B|=%d" % len(span))
# label functionals: pairing with the two uniform primal loops (row i0=0, column j0=0)
wA = [0] * nE
for j in range(LL):
    wA[h(0, j)] = 1          # primal row loop
wB = [0] * nE
for i in range(LL):
    wB[v(i, 0)] = 1          # primal column loop
vanA = all(dot4(wA, s) == 0 for s in span)
vanB = all(dot4(wB, s) == 0 for s in span)
ALL_OK &= gate("D(Z4) both Wilson functionals vanish on the whole star span (gauge invariant)",
               vanA and vanB)
from collections import Counter
labels = Counter()
class_elems = {}
for u in C:
    t = (dot4(u, wA), dot4(u, wB))
    labels[t] += 1
    class_elems.setdefault(t, []).append(u)
ALL_OK &= gate("D(Z4) the two functionals classify C into 16 equal classes of 64",
               len(labels) == 16 and all(c == 64 for c in labels.values()),
               str(dict(labels)) if not (len(labels) == 16) else "")
def wt4(u):
    return sum(1 for x in u if x != 0)
class_min4 = {t: min(wt4(u) for u in class_elems[t]) for t in class_elems}
say("D(Z4)  class minimal weights: " + str({t: class_min4[t] for t in sorted(class_min4)}))
# configurations = Z_4^2 labels; d(s,s') = min weight over the class of s'-s (searched above)
configs4 = [(a, b) for a in range(4) for b in range(4)]
dmat4 = {}
for s in configs4:
    for sp_ in configs4:
        t = ((sp_[0] - s[0]) % 4, (sp_[1] - s[1]) % 4)
        dmat4[(s, sp_)] = class_min4[t]
rep = metric_axioms(configs4, dmat4)
ALL_OK &= gate("D(Z4) METRIC AXIOMS over ALL %d pairs, ALL %d triples"
               % (rep["n_pairs"], rep["n_triples"]), rep["all"],
               "" if rep["all"] else str(rep))
# structure: support-Hamming vs Lee weight
supp_ok = all(class_min4[t] == LL * sum(1 for x in t if x != 0) for t in class_min4)
lee = lambda x: min(x, 4 - x)
lee_ok = all(class_min4[t] == LL * (lee(t[0]) + lee(t[1])) for t in class_min4)
ALL_OK &= gate("D(Z4) d_W = L x |support(t)|  (torsion is FREE: d(2t)=d(t) for t nonzero)",
               supp_ok)
ALL_OK &= gate("D(Z4) D-15 CONTROL: the imported Lee-weight form L*(lee(a)+lee(b)) FAILS",
               not lee_ok,
               "witness t=(2,0): d=%d, Lee form=%d" % (class_min4[(2, 0)], LL * 2))
kind, na4 = aut_count(configs4, dmat4)
say("D(Z4)  D-22 venue Aut of the 16-point metric: |Aut| = %d (%s)" % (na4, kind))
ALL_OK &= gate("D(Z4) venue has geometry (Aut != S_16)", na4 != 20922789888000)

# ---- are the two |G|=4 doubles the same 16-point geometry?  (computed, backtracking)
say("")
say("CROSS-DOUBLE COMPARISON: D(Z2xZ2) L=2 vs D(Z4) 2x2 -- same 16-point metric space?")
cfgA, dA, _ = Z22_METRICS[2]
spectA = sorted(dA[(p, q)] for p in cfgA for q in cfgA)
spectB = sorted(dmat4[(p, q)] for p in configs4 for q in configs4)
same_spect = (spectA == spectB)
say("distance multisets equal: %s" % same_spect)
iso_found = False
if same_spect:
    # backtracking isometry search
    A = list(cfgA); B = list(configs4)
    DA = [[dA[(p, q)] for q in A] for p in A]
    DB = [[dmat4[(p, q)] for q in B] for p in B]
    profA = [tuple(sorted(r)) for r in DA]
    profB = [tuple(sorted(r)) for r in DB]
    def bt(k, perm):
        if k == 16:
            return True
        for c in range(16):
            if c in perm[:k] or profB[c] != profA[k]:
                continue
            if all(DA[k][j] == DB[c][perm[j]] for j in range(k)):
                perm[k] = c
                if bt(k + 1, perm):
                    return True
        return False
    iso_found = bt(0, [0] * 16)
if iso_found:
    say("FINDING: an explicit isometry EXISTS -- the two non-isomorphic order-4 doubles "
        "carry THE SAME 16-point configuration geometry at this size (the metric does not "
        "see the group's torsion here)")
else:
    say("FINDING: no isometry -- the two order-4 doubles carry DIFFERENT 16-point "
        "geometries (the metric distinguishes the groups)")
say("")
say("ALL GATES: " + ("PASS" if ALL_OK else "AT LEAST ONE FAIL"))
with open(LANE + "/t42_c_doubles.OUT.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
