"""T42-A part 2: d_W BETWEEN CONFIGURATIONS on [[n,n-2,2]] carriers, n = 4..8.  EXACT, F_2.

Even n: the canonical two-generator code  S = <X^n, Z^n>.
Odd n:  X^n and Z^n anticommute, so the canonical pair does not exist; a two-generator
        distance-2 code is SEARCHED (deterministic seed), never nominated, and verified.

For every ordered configuration pair the FULL affine space of admissible Pauli writers with
that label action is enumerated over (x|z) -- no CSS reduction assumed -- and the minimum
weight (|supp x  U  supp z|) taken.

QUESTIONS ANSWERED WITH GATES:
  * metric axioms over ALL pairs and ALL triples,
  * group invariance of d_W (emerges from per-pair search, not assumed),
  * whether d_W = d_code x Hamming(labels) -- expected to FAIL here (this is the D-15
    control for the toric relation: the instrument must be able to say NO),
  * the counting obstruction: no linear relabeling can ever make it Hamming-proportional
    when #{classes at distance d_code} != #labels of Hamming weight 1,
  * basis dependence: two COMPUTED symplectic label bases give different d(label) tables,
    related by an explicit F_2-linear isometry -- the metric space is the invariant.
"""
import sys, time, random
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_T42_A_DISTANCE"
sys.path.insert(0, LANE)
from t42_lib import (pc, solve_affine_f2, span_all, sp_pair, weight_xz, vec_to_mask,
                     metric_axioms, aut_count)
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import symplectic_logicals

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s); print(s); sys.stdout.flush()

def gate(label, ok, detail=""):
    say(("PASS  " if ok else "FAIL  ") + label + (("  " + detail) if detail else ""))
    return ok

def full_row(g, n):
    mask = (1 << n) - 1
    return (g >> n) | ((g & mask) << n)

def mask_to_vec(m, ln):
    return [(m >> i) & 1 for i in range(ln)]

def code_distance_two_gen(g1, g2, n):
    """min weight over N(S) \\ span{g1,g2}; N(S) enumerated exactly."""
    rows = [full_row(g1, n), full_row(g2, n)]
    x0, nb, rank = solve_affine_f2(rows, [0, 0], 2 * n)
    S = {0, g1, g2, g1 ^ g2}
    best = None
    for v in span_all(nb):
        if v in S:
            continue
        w = weight_xz(v, n)
        if best is None or w < best:
            best = w
    return best

def find_odd_carrier(n, seed=42, attempts=50000):
    """SEARCH a two-generator [[n, n-2, 2]] carrier.  Deterministic seed; verified exactly."""
    rng = random.Random(seed)
    for _ in range(attempts):
        g1 = rng.randrange(1, 1 << (2 * n))
        g2 = rng.randrange(1, 1 << (2 * n))
        if g2 in (0, g1):
            continue
        if sp_pair(g1, g2, n) != 0:
            continue
        d = code_distance_two_gen(g1, g2, n)
        if d == 2:
            return g1, g2
    return None, None


def weight1_masks(n):
    out = []
    for i in range(n):
        out += [1 << i, 1 << (n + i), (1 << i) | (1 << (n + i))]
    return out


def has_weight1_logical(g1, g2, n):
    Sbar = {0, g1, g2, g1 ^ g2}
    for w in weight1_masks(n):
        if w in Sbar:
            continue
        if sp_pair(w, g1, n) == 0 and sp_pair(w, g2, n) == 0:
            return True
    return False


def odd_n_nonexistence(n):
    """THEOREM (computed): no two-generator [[n, n-2, 2]] exists for odd n.
       Chain, every link computed here or exhausted:
       (A) 16-case local lemma: at one site with local generator pair (a,b), all three
           weight-1 Paulis anticommute with a or b  IFF  a,b are distinct and non-identity
           IFF  a,b locally anticommute.
       (B) therefore a pair with NO weight-1 Pauli commuting with both generators must be
           locally anticommuting at EVERY site, so sp(g1,g2) = n mod 2: for odd n every
           commuting pair has a site that is not locally anticommuting.
       (C) the surviving weight-1 Pauli at such a site could only fail to be a logical by
           lying in S-bar = {I, g1, g2, g1g2}; a weight-1 member of S-bar means one of
           g1, g2, g1g2 has weight 1 -- that subfamily is EXHAUSTED directly.
       Returns (lemma_ok, subfamily_ok, n_subfamily_pairs)."""
    # (A) the 16-case local lemma, exact
    def local_anti(a, b):  # a,b in {0..3} encoding (x,z): 0=I,1=X,2=Z,3=Y
        ax, az = a & 1, (a >> 1) & 1
        bx, bz = b & 1, (b >> 1) & 1
        return (ax * bz + az * bx) % 2 == 1
    lemma_ok = True
    for a in range(4):
        for b in range(4):
            kills_all = all(
                ((px * ((a >> 1) & 1) + pz * (a & 1)) % 2 == 1) or
                ((px * ((b >> 1) & 1) + pz * (b & 1)) % 2 == 1)
                for (px, pz) in ((1, 0), (0, 1), (1, 1)))
            lemma_ok &= (kills_all == (a != 0 and b != 0 and a != b))
            lemma_ok &= ((a != 0 and b != 0 and a != b) == local_anti(a, b))
    # (C) exhaust the subfamily where some element of S-bar has weight 1
    sub_ok = True
    cnt = 0
    for w in weight1_masks(n):
        g1 = w
        for g2 in range(1, 1 << (2 * n)):
            if g2 == g1:
                continue
            if sp_pair(g1, g2, n) != 0:
                continue
            cnt += 1
            # d >= 2 would need no weight-1 logical
            if not has_weight1_logical(g1, g2, n):
                sub_ok = False
    return lemma_ok, sub_ok, cnt

say("T42-A part 2  --  [[n,n-2,2]] configuration distance, exact, started "
    + time.strftime("%F %T"))
say("")

ALL_OK = True
SUMMARY = {}

for n in (4, 5, 6, 7, 8):
    t0 = time.time()
    if n % 2 == 0:
        g1 = vec_to_mask([1] * n + [0] * n)
        g2 = vec_to_mask([0] * n + [1] * n)
        origin = "canonical <X^%d, Z^%d>" % (n, n)
    else:
        g1, g2 = find_odd_carrier(n)
        if g1 is None:
            say("n=%d  no two-generator [[%d,%d,2]] carrier found by search; "
                "testing NONEXISTENCE exactly" % (n, n, n - 2))
            lemma_ok, sub_ok, cnt = odd_n_nonexistence(n)
            ALL_OK &= gate("n=%d local 16-case lemma (kills-all == locally-anticommuting)"
                           % n, lemma_ok)
            ALL_OK &= gate("n=%d weight-1-generator subfamily EXHAUSTED, every pair has a "
                           "weight-1 logical (%d commuting pairs)" % (n, cnt), sub_ok)
            if n == 5:
                # full exhaustive confirmation over ALL commuting independent pairs
                W1 = weight1_masks(n)
                comm_mask = {}
                for g in range(1 << (2 * n)):
                    m = 0
                    for i, w in enumerate(W1):
                        if sp_pair(w, g, n) == 0:
                            m |= 1 << i
                    comm_mask[g] = m
                bad = 0
                tot = 0
                for g1_ in range(1, 1 << (2 * n)):
                    for g2_ in range(g1_ + 1, 1 << (2 * n)):
                        if sp_pair(g1_, g2_, n) != 0:
                            continue
                        tot += 1
                        both = comm_mask[g1_] & comm_mask[g2_]
                        Sbar = {0, g1_, g2_, g1_ ^ g2_}
                        found = False
                        m = both
                        while m:
                            i = (m & -m).bit_length() - 1
                            m &= m - 1
                            if W1[i] not in Sbar:
                                found = True
                                break
                        if not found:
                            bad += 1
                ALL_OK &= gate("n=5 EXHAUSTIVE over all %d commuting pairs: every one has a "
                               "weight-1 logical -> [[5,3,2]] (2-gen) DOES NOT EXIST" % tot,
                               bad == 0, "pairs with d>=2: %d" % bad)
            else:
                say("n=7  nonexistence: computed lemma + symplectic parity "
                    "(commuting pair on odd n cannot be locally anticommuting at every "
                    "site) + exhausted subfamily close every case")
            say("n=%d  FINDING: the [[n,n-2,2]] carrier EXISTS iff n is even; the "
                "obstruction IS the symplectic parity sp(g1,g2) = n mod 2" % n)
            say("")
            continue
        origin = "SEARCHED (seed 42), two commuting generators found and verified"
    dcode = code_distance_two_gen(g1, g2, n)
    ALL_OK &= gate("n=%d carrier is a distance-2 code (computed)" % n, dcode == 2,
                   origin + "  d_code=%d" % dcode)
    stab = [mask_to_vec(g1, 2 * n), mask_to_vec(g2, 2 * n)]
    pairs = symplectic_logicals(stab, n)
    k = n - 2
    ALL_OK &= gate("n=%d symplectic_logicals returns %d conjugate pairs (computed)" % (n, k),
                   len(pairs) == k)
    logi = [(vec_to_mask(X), vec_to_mask(Z)) for X, Z in pairs]
    Zlogs = [z for _, z in logi]
    configs = [tuple((m >> i) & 1 for i in range(k)) for m in range(1 << k)]
    srows = [full_row(g1, n), full_row(g2, n)]
    lrows = [full_row(z, n) for z in Zlogs]

    dmat = {}
    class_min = {}
    searched = 0
    space_dim = None
    for s in configs:
        for sp_ in configs:
            t = tuple(a ^ b for a, b in zip(s, sp_))
            x0, nb, rank = solve_affine_f2(srows + lrows, [0, 0] + list(t), 2 * n)
            assert x0 is not None
            if space_dim is None:
                space_dim = len(nb)
            best = min(weight_xz(x0 ^ v, n) for v in span_all(nb))
            searched += 1 << len(nb)
            dmat[(s, sp_)] = best
            class_min.setdefault(t, best)
    say("n=%d  configurations: %d   writer space per pair: 2^%d   writers searched: %d"
        % (n, len(configs), space_dim, searched))

    inv_ok = all(dmat[(s, sp_)] == class_min[tuple(a ^ b for a, b in zip(s, sp_))]
                 for s in configs for sp_ in configs)
    ALL_OK &= gate("n=%d d_W(s,s') depends only on s+s' (invariance emerged from search)" % n,
                   inv_ok)

    rep = metric_axioms(configs, dmat)
    ALL_OK &= gate("n=%d METRIC AXIOMS over ALL %d pairs and ALL %d triples"
                   % (n, rep["n_pairs"], rep["n_triples"]), rep["all"],
                   "" if rep["all"] else str({kk: rep[kk] for kk in
                                              ("identity_ok", "symmetry_ok",
                                               "indiscernible_ok", "triangle_ok")}))

    # d vs Hamming of the computed labels
    by_h = {}
    for t, dv in class_min.items():
        by_h.setdefault(sum(t), set()).add(dv)
    say("n=%d  d_W by Hamming weight of computed label: %s"
        % (n, {h: sorted(v) for h, v in sorted(by_h.items())}))
    ham_prop = all(v == {dcode * h} for h, v in by_h.items())
    ALL_OK &= gate("n=%d D-15 CONTROL: d_W DEVIATES from d_code x Hamming(labels) here "
                   "(the instrument can say no; the toric relation is not an artifact)" % n,
                   not ham_prop)
    # counting obstruction to ANY linear Hamming relabeling
    n_at_dcode = sum(1 for t, dv in class_min.items() if dv == dcode and any(t))
    obstruct = (n_at_dcode != k)
    say("n=%d  #classes at distance d_code = %d   #Hamming-weight-1 labels = %d  ->  %s"
        % (n, n_at_dcode, k,
           "NO linear relabeling can make d_W Hamming-proportional" if obstruct
           else "counting permits a Hamming relabeling"))

    # D-22: venue Aut
    N = len(configs)
    if N <= 16:
        kind, na = aut_count(configs, dmat)
        full_sym = 1
        for i in range(2, N + 1):
            full_sym *= i
        say("n=%d  D-22 venue Aut of the %d-point metric space: |Aut| = %d (%s); S_N = %d"
            % (n, N, na, kind, full_sym))
        if na == full_sym:
            say("n=%d  D-22 FINDING: the metric is UNIFORM -- this venue carries NO geometry "
                "to detect; no separation claim may be made on it" % n)
        else:
            say("n=%d  D-22: venue HAS geometry (Aut is a proper subgroup)" % n)
    else:
        # exact subgroup: translations x label-coordinate permutations.  The invariance gate
        # above already proves every translation is an isometry: d(s+g, s'+g) = f(s+s').
        perm_ok = True
        for i in range(k - 1):  # adjacent transpositions generate S_k
            def swap(t, i=i):
                t = list(t); t[i], t[i + 1] = t[i + 1], t[i]; return tuple(t)
            perm_ok &= all(class_min[swap(t)] == class_min[t] for t in class_min)
        sub = (1 << k)
        if perm_ok:
            fact = 1
            for i in range(2, k + 1):
                fact *= i
            sub *= fact
        vals = sorted(set(class_min.values()))
        nonuniform = len([v for v in vals if v != 0]) > 1
        say("n=%d  D-22 venue Aut: exact enumeration skipped at N=%d; VERIFIED isometry "
            "subgroup order >= %d (translations%s)" % (n, N, sub,
            " + label-coordinate S_k" if perm_ok else ""))
        ALL_OK &= gate("n=%d venue has geometry to detect (nonuniform distance spectrum %s "
                       "-> Aut != S_%d)" % (n, vals, N), nonuniform)

    SUMMARY[n] = {"dcode": dcode, "k": k, "classes_at_dcode": n_at_dcode,
                  "spectrum": sorted(set(class_min.values())),
                  "by_h": {h: sorted(v) for h, v in sorted(by_h.items())}}
    say("n=%d done in %.1fs" % (n, time.time() - t0))
    say("")

# ---------------- basis dependence demonstration on n=6
say("BASIS DEPENDENCE (n=6): two computed symplectic label bases, one metric space")
n = 6
g1 = vec_to_mask([1] * n + [0] * n)
g2 = vec_to_mask([0] * n + [1] * n)
stab = [mask_to_vec(g1, 2 * n), mask_to_vec(g2, 2 * n)]
pairs = symplectic_logicals(stab, n)
k = n - 2
logi = [(vec_to_mask(X), vec_to_mask(Z)) for X, Z in pairs]
# basis 2: X1' = X1^X2, Z2' = Z1^Z2 (a symplectic transvection pair), rest unchanged
logi2 = [row[:] if isinstance(row, list) else row for row in logi]
logi2 = list(logi)
X1, Z1 = logi[0]; X2, Z2 = logi[1]
logi2[0] = (X1 ^ X2, Z1)
logi2[1] = (X2, Z1 ^ Z2)
sym_ok = all(sp_pair(logi2[i][0], logi2[i][1], n) == 1 for i in range(k)) and \
         all(sp_pair(logi2[i][a], logi2[j][b], n) == 0
             for i in range(k) for j in range(k) for a in (0, 1) for b in (0, 1) if i != j)
ALL_OK &= gate("n=6 second basis is symplectic (computed check)", sym_ok)

def metric_for(Zlogs):
    srows = [full_row(g1, n), full_row(g2, n)]
    lrows = [full_row(z, n) for z in Zlogs]
    cm = {}
    for m in range(1 << k):
        t = tuple((m >> i) & 1 for i in range(k))
        x0, nb, _ = solve_affine_f2(srows + lrows, [0, 0] + list(t), 2 * n)
        cm[t] = min(weight_xz(x0 ^ v, n) for v in span_all(nb))
    return cm

cm1 = metric_for([z for _, z in logi])
cm2 = metric_for([z for _, z in logi2])
differs = any(cm1[t] != cm2[t] for t in cm1)
same_multiset = sorted(cm1.values()) == sorted(cm2.values())
ALL_OK &= gate("n=6 the per-label table CHANGES with the basis", differs,
               "example: " + str(next((t, cm1[t], cm2[t]) for t in cm1 if cm1[t] != cm2[t]))
               if differs else "")
ALL_OK &= gate("n=6 the distance multiset is basis-invariant", same_multiset)
# explicit linear isometry: label2(u) = pairings of u against Zlogs2; on classes this is a
# linear map M of labels; build M by evaluating on the X-logical generators of basis 1
M = []  # columns: image of unit label e_i
for i in range(k):
    u = logi[i][0]  # X_i of basis 1 flips label e_i in basis 1
    M.append(tuple(sp_pair(u, logi2[j][1], n) for j in range(k)))
def applyM(t):
    out = [0] * k
    for i, bit in enumerate(t):
        if bit:
            out = [a ^ b for a, b in zip(out, M[i])]
    return tuple(out)
iso_ok = all(cm2[applyM(t)] == cm1[t] for t in cm1)
ALL_OK &= gate("n=6 explicit F_2-linear relabeling is an isometry between the two tables",
               iso_ok, "M columns=" + str(M))
say("")
say("SUMMARY: " + str(SUMMARY))
say("ALL GATES: " + ("PASS" if ALL_OK else "AT LEAST ONE FAIL"))
with open(LANE + "/t42_b_codes.OUT.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
