"""S1b -- IS THE TOTAL MINIMUM WRITER WEIGHT BASIS-DEPENDENT?

S1 found Wflip_tot = 4, 8, 16, 24, 36, 48, 64 for k = 2..14 in the basis symplectic_logicals
happens to return -- roughly k^2/4, SUPER-LINEAR.  Before that can be called a growth law it
has to survive D-17: vary the venue's own scale.  Here the venue's own scale is THE CHOICE OF
SYMPLECTIC BASIS for the record family; a different conjugate basis is the same code with the
record set relabelled.

TWO INDEPENDENT PROBES OF BASIS DEPENDENCE.
  (A) RANDOM SYMPLECTIC BASIS CHANGES applied to the pairs symplectic_logicals returns, built
      from three moves that preserve the canonical Gram matrix:
          CNOT(i->j) : R_j += R_i, W_i += W_j     SWAP(i,j)      HAD(i) : R_i <-> W_i
      The Gram matrix is re-verified after every move; failures are discarded.
  (B) A DETERMINISTIC SEARCH FOR A MINIMUM-WEIGHT CONJUGATE BASIS.  No operator is nominated:
      we enumerate EVERY weight-2 element of N(S)\S from scratch, then greedily extract k
      conjugate pairs from that pool and VERIFY (a) the full canonical Gram matrix and (b) that
      the 2k operators together with the 2 stabilisers have F_2 rank 2n-2, i.e. they really do
      span N(S).  Only a set that passes both checks is used.

Lower bound: every writer is a non-identity logical, so its weight is >= d = 2; hence
Wflip_tot >= 2k for ANY basis.
"""
import sys, random, json, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import symplectic_logicals

OUT = []
def p(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)

POP = [bin(i).count("1") for i in range(1 << 16)]
def popcount(v):
    c = 0
    while v: c += POP[v & 0xFFFF]; v >>= 16
    return c

def to_int(vec, n):
    x = z = 0
    for i in range(n):
        if vec[i]: x |= 1 << i
        if vec[n + i]: z |= 1 << i
    return (x, z)

def wt(v): return popcount(v[0] | v[1])
def xor(a, b): return (a[0] ^ b[0], a[1] ^ b[1])
def spi(a, b): return (popcount(a[0] & b[1]) + popcount(a[1] & b[0])) % 2

def gram_ok(R, W, k):
    for i in range(k):
        for j in range(k):
            if spi(R[i], R[j]) != 0: return False
            if spi(W[i], W[j]) != 0: return False
            if spi(R[i], W[j]) != (1 if i == j else 0): return False
    return True

def f2_rank(vecs, n):
    rows = [(v[0] | (v[1] << n)) for v in vecs]
    r = 0
    for bit in range(2 * n):
        piv = next((i for i in range(r, len(rows)) if rows[i] >> bit & 1), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and (rows[i] >> bit & 1): rows[i] ^= rows[r]
        r += 1
    return r

def flip_only_total(R, W, k, stabg):
    """sum_i min weight over  W_i * S * <R_1..R_k>  -- writers that flip ONLY record i."""
    total, per = 0, []
    for i in range(k):
        best = None
        for g in stabg:
            cur = xor(W[i], g)
            w = wt(cur)
            if best is None or w < best: best = w
            m_prev = 0
            for m in range(1, 1 << k):
                gray = m ^ (m >> 1)
                b = (gray ^ m_prev).bit_length() - 1
                cur = xor(cur, R[b]); m_prev = gray
                w = wt(cur)
                if w < best: best = w
        per.append(best); total += best
    return total, per

def weight2_logicals(n, stab_i):
    """EVERY weight-2 element of N(S)\\S, enumerated -- nothing nominated."""
    out = []
    for a, b in itertools.combinations(range(n), 2):
        for ka in (1, 2, 3):
            for kb in (1, 2, 3):
                x = z = 0
                for s, kd in ((a, ka), (b, kb)):
                    if kd in (1, 3): x |= 1 << s
                    if kd in (2, 3): z |= 1 << s
                v = (x, z)
                if all(spi(v, s) == 0 for s in stab_i):
                    out.append(v)
    return out

def greedy_weight2_basis(pool, k, n, stab_i, rng, restarts=400):
    """Search the ENUMERATED weight-2 pool for k conjugate pairs.

    Structure used (derived, not assumed): within the pool, two elements commute automatically
    unless one is X-type and the other Z-type, in which case the symplectic form is the parity
    of their site-overlap.  So a candidate basis is fixed by choosing k pool elements as records
    and finding, for each, a pool element whose overlap-parity vector is the i-th unit vector.
    Every candidate is then put through the SAME two verifications: canonical Gram matrix, and
    F_2 rank 2k+2 together with the stabilisers.  Nothing is nominated and nothing is assumed."""
    def sites(v): return v[0] | v[1]
    Xtype = [v for v in pool if v[1] == 0]
    Ztype = [v for v in pool if v[0] == 0]
    if not Xtype or not Ztype: return None, None
    # candidate record sets: every "star" available in the pool, then random k-subsets
    cand_sets = []
    for c in range(n):
        st = [v for v in Xtype if (sites(v) >> c) & 1]
        if len(st) >= k: cand_sets.append(st[:k])
    for _ in range(restarts):
        cand_sets.append(rng.sample(Xtype, k))
    for R in cand_sets:
        Rs = [sites(v) for v in R]
        W = []
        ok = True
        for i in range(k):
            hit = None
            for q in Ztype:
                sq = sites(q)
                if all((popcount(Rs[j] & sq) % 2) == (1 if j == i else 0) for j in range(k)):
                    hit = q; break
            if hit is None: ok = False; break
            W.append(hit)
        if not ok: continue
        if gram_ok(R, W, k) and f2_rank(list(R) + W + list(stab_i), n) == 2 * k + 2:
            return list(R), W
    return None, None

res = {}
p("=" * 112)
p("S1b  TOTAL MINIMUM FLIP-ONLY WRITER WEIGHT vs CHOICE OF SYMPLECTIC BASIS  (D-17: vary the venue's own scale)")
p("=" * 112)
rng = random.Random(11)
NS = [4, 6, 8, 10, 12]
for n in NS:
    k = n - 2
    stab = [[1] * n + [0] * n, [0] * n + [1] * n]
    pairs = symplectic_logicals([s[:] for s in stab], n)
    if len(pairs) != k:
        p("SELF-CHECK FAILED at n=%d -- CONCLUDING NOTHING" % n); sys.exit(1)
    R = [to_int(a, n) for a, b in pairs]
    W = [to_int(b, n) for a, b in pairs]
    if not gram_ok(R, W, k):
        p("SELF-CHECK FAILED: Gram matrix wrong at n=%d -- CONCLUDING NOTHING" % n); sys.exit(1)
    S1i, S2i = to_int(stab[0], n), to_int(stab[1], n)
    stabg = [(0, 0), S1i, S2i, xor(S1i, S2i)]

    gs_tot, gs_per = flip_only_total(R, W, k, stabg)

    # (A) random symplectic basis changes
    best_tot, best_per = gs_tot, gs_per
    cur_R, cur_W = list(R), list(W)
    trials = 400 if k <= 8 else 120
    for _ in range(trials):
        RR, WW = list(cur_R), list(cur_W)
        for _m in range(rng.randint(1, 3 * k)):
            mv = rng.random(); i, j = rng.randrange(k), rng.randrange(k)
            if mv < 0.6 and i != j: RR[j] = xor(RR[j], RR[i]); WW[i] = xor(WW[i], WW[j])
            elif mv < 0.8 and i != j: RR[i], RR[j] = RR[j], RR[i]; WW[i], WW[j] = WW[j], WW[i]
            else: RR[i], WW[i] = WW[i], RR[i]
        if not gram_ok(RR, WW, k): continue
        t, per = flip_only_total(RR, WW, k, stabg)
        if t < best_tot: best_tot, best_per = t, per
        if t <= gs_tot: cur_R, cur_W = RR, WW

    # (B) deterministic construction from the enumerated weight-2 logicals
    pool = weight2_logicals(n, [S1i, S2i])
    R2, W2 = greedy_weight2_basis(pool, k, n, [S1i, S2i], rng)
    if R2 is None:
        w2_tot, w2_per, w2_ok = None, None, False
    else:
        w2_tot, w2_per = flip_only_total(R2, W2, k, stabg)
        w2_ok = True

    res[n] = dict(n=n, k=k, gram_schmidt_total=gs_tot, gram_schmidt_per=gs_per,
                  random_best_total=best_tot, random_best_per=best_per,
                  weight2_basis_found=w2_ok, weight2_pool_size=len(pool),
                  weight2_total=w2_tot, weight2_per=w2_per,
                  lower_bound_2k=2 * k, trials=trials)

p("")
p("  n   k | pool of weight-2 logicals | Wflip_tot Gram-Schmidt | Wflip_tot best random basis | Wflip_tot verified weight-2 basis | exact lower bound 2k")
p("-" * 112)
for n in NS:
    r = res[n]
    p("%3d %3d | %24d | %22d | %27d | %33s | %19d"
      % (r["n"], r["k"], r["weight2_pool_size"], r["gram_schmidt_total"], r["random_best_total"],
         str(r["weight2_total"]) if r["weight2_basis_found"] else "SEARCH FAILED", r["lower_bound_2k"]))
p("-" * 112)
p("")
attained = [n for n in NS if res[n]["weight2_basis_found"] and res[n]["weight2_total"] == 2 * (n - 2)]
p("READ (filled from the numbers above, not in advance):")
p("  the verified weight-2 conjugate basis attains the exact lower bound 2k at n = %s   (of %s tested)"
  % (attained, NS))
if len(attained) == len(NS):
    p("  the Gram-Schmidt column's ~k^2/4 rise is therefore an artefact of the BASIS, not a property of the code;")
    p("  it does not survive D-17.  The basis-independent minimum is Wflip_tot = 2k = 2(n-2), exactly LINEAR,")
    p("  and the per-record minimum is exactly 2 = d, CONSTANT at every N.")
else:
    p("  the weight-2 construction did NOT attain 2k everywhere; the basis-independent minimum is bracketed")
    p("  between 2k and the smallest value found, and no growth law is claimed.")
p("  NOTE: the random-basis column is a search, not a minimisation -- it is an upper bound only,")
p("  and it is superseded by the verified weight-2 column, which is a CONSTRUCTION.")

with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_A_ZERO_OR_GROWING/s1b_writer_weight_basis.json", "w") as f:
    json.dump({str(a): b for a, b in res.items()}, f, indent=1)
with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_A_ZERO_OR_GROWING/s1b_writer_weight_basis.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
