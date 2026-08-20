"""LANE_SCALE_C_ACCUM  --  script 1: EXACT combinatorics of the record structure at scale.

QUESTION: does anything ACCUMULATE with the number of records k -- the thing a COUNT could not?

CARRIER : [[n, n-2, 2]], n even, stabilisers X^(x)n and Z^(x)n, k = n-2 records, dim 2^n.
CONTROL : m independent NON-INTERACTING [[4,2,2]] blocks, k = 2m records, n_ctrl = 4m.
          Anything that grows the SAME way there is additive BY CONSTRUCTION, not collective.

Everything here is EXACT F_2 linear algebra on (x|z) Pauli vectors and exhaustive enumeration
of the admissible Pauli group -- no 2^n density matrix is ever built, so n runs far past what
script 2 (chi) can reach.

NEVER NOMINATE LOGICAL OPERATORS: every record comes from symplectic_logicals().

THE BASIS PROBLEM, and how it is handled.  "How many records does this operator disturb" and
"what does it cost to flip record i alone" are NOT invariants: a change of record basis
R -> A.R (A in GL(k,2)) permutes them.  The invariant object is the whole map

      w : F_2^k -> N,   w(v) = min weight of an ADMISSIBLE Pauli P with  sp(P,R_i) = v_i

("which records does P flip" -> "what does the cheapest such P weigh").  A change of record
basis PERMUTES w by the GL(k,2) action, so max_v w(v), mean_v w(v), and the weight at which
the realised patterns first SPAN F_2^k are basis-INVARIANT.  Those are the reported quantities;
basis-dependent ones are printed too but flagged.
"""
import sys, itertools, time
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
import numpy as np
from record_model import symplectic_logicals

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.append(s)

# ------------------------------------------------------------------ F_2 Pauli helpers
def sp2(a, b, n):
    """symplectic form: 1 iff the two Paulis ANTICOMMUTE"""
    return (sum(a[i] * b[n + i] for i in range(n)) + sum(a[n + i] * b[i] for i in range(n))) % 2

def f2_rank(rows):
    if not rows: return 0
    rows = [list(r) for r in rows]; r = 0; L = len(rows[0])
    for c in range(L):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [(x + y) % 2 for x, y in zip(rows[i], rows[r])]
        r += 1
    return r

def in_span(v, basis):
    return f2_rank(list(basis) + [list(v)]) == f2_rank(list(basis))

def masks(v, n):
    x = z = 0
    for i in range(n):
        if v[i]: x |= (1 << i)
        if v[n + i]: z |= (1 << i)
    return x, z

# ------------------------------------------------------------------ carriers
def carrier_nn2(n):
    return n, [[1] * n + [0] * n, [0] * n + [1] * n]

def carrier_product(m):
    n = 4 * m; S = []
    for b in range(m):
        x = [0] * (2 * n); z = [0] * (2 * n)
        for q in range(4 * b, 4 * b + 4):
            x[q] = 1; z[n + q] = 1
        S.append(x); S.append(z)
    return n, S

def embed(v, n_small, offset, n_big):
    out = [0] * (2 * n_big)
    for i in range(n_small):
        out[offset + i] = v[i]; out[n_big + offset + i] = v[n_small + i]
    return out

# ------------------------------------------------------------------ the invariant weight map
def weight_map(n, S, R, budget=2.0e8, verbose=False):
    """EXACT w(v) for EVERY v in F_2^k, by exhaustive enumeration of the admissible Pauli group.

       Admissible = commutes with every stabiliser (<=> [P,H]=0 since H = -sum of stabilisers).
       Both carriers here have stabilisers that are PURE X-type or PURE Z-type, so admissibility
       factorises: a pure-X stabiliser constrains only the z-part, a pure-Z one only the x-part.
       That is asserted, not assumed.  Returns (w array of length 2^k, n_admissible) or None if
       the enumeration would exceed `budget` Pauli operators."""
    k = len(R); full = (1 << n) - 1
    xcon, zcon = [], []            # masks: parity of popcount(x & mask) / popcount(z & mask) must be 0
    for s in S:
        sx, sz = masks(s, n)
        if sz == 0:   zcon.append(sx)       # pure X stabiliser: sp = popcount(z & sx)
        elif sx == 0: xcon.append(sz)       # pure Z stabiliser: sp = popcount(x & sz)
        else: return None, None, "stabiliser is neither pure X nor pure Z"
    # SIZE FIRST: the admissible group has index 2^(#independent stabilisers), so its order is
    # known before anything is allocated.  Bail before building arrays that cannot fit.
    est = 1 << (2 * n - f2_rank(S))
    if est > budget: return None, est, "over budget"
    def admissible_half(cons):
        v = np.arange(1 << n, dtype=np.uint64)
        keep = np.ones(v.size, dtype=bool)
        for msk in cons:
            keep &= (np.bitwise_count(np.bitwise_and(v, np.uint64(msk))) & 1) == 0
        return v[keep]
    xs = admissible_half(xcon); zs = admissible_half(zcon)
    total = int(xs.size) * int(zs.size)
    if total > budget: return None, total, "over budget"
    Rx = np.array([masks(r, n)[0] for r in R], dtype=np.uint64)
    Rz = np.array([masks(r, n)[1] for r in R], dtype=np.uint64)
    SHIFT = 5                                    # weight fits in 5 bits for n < 32
    seen = np.zeros((1 << k) << SHIFT, dtype=bool)
    per = max(1, int(2.0e6 // max(1, zs.size)))
    for a in range(0, xs.size, per):
        X = xs[a:a + per][:, None]; Z = zs[None, :]
        w = np.bitwise_count(np.bitwise_or(X, Z)).astype(np.uint64)
        pat = np.zeros(w.shape, dtype=np.uint64)
        for i in range(k):
            b = (np.bitwise_count(np.bitwise_and(X, Rz[i]))
                 + np.bitwise_count(np.bitwise_and(Z, Rx[i]))) & np.uint64(1)
            pat |= np.left_shift(b.astype(np.uint64), np.uint64(i))
        code = np.bitwise_or(np.left_shift(pat, np.uint64(SHIFT)), w)
        seen[np.unique(code)] = True
    seen = seen.reshape(1 << k, 1 << SHIFT)
    wv = np.full(1 << k, -1, dtype=np.int64)
    for wgt in range(1 << SHIFT):
        hit = seen[:, wgt] & (wv < 0)
        wv[hit] = wgt
    return wv, total, "exact"

def wmap_summary(wv, k):
    """basis-INVARIANT summaries of the weight map"""
    nz = wv[1:]                                   # drop the identity pattern v = 0
    assert wv[0] == 0, "the zero pattern must be realised at weight 0 by the identity"
    unreached = int((nz < 0).sum())
    good = nz[nz >= 0]
    # smallest weight at which the realised patterns SPAN F_2^k
    span_at = None
    for wgt in range(0, int(good.max()) + 1):
        pats = [i + 1 for i in range(len(nz)) if 0 <= nz[i] <= wgt]
        rows = [[(p >> j) & 1 for j in range(k)] for p in pats]
        if f2_rank(rows) == k: span_at = wgt; break
    return dict(wmax=int(good.max()), wmin=int(good.min()), wmean=float(good.mean()),
                wmean_all=float(wv.mean()), unreached=unreached, span_at=span_at,
                w_allones=int(wv[(1 << k) - 1]), exact=True)

# ------------------------------------------------------------------ cheap (large-n) quantities
def cheap(n, S, R, nbasis=25, seed=0):
    k = len(R)
    def ops_of_weight(w):
        for supp in itertools.combinations(range(n), w):
            for tags in itertools.product((1, 2, 3), repeat=w):
                v = [0] * (2 * n)
                for q, t in zip(supp, tags):
                    if t & 1: v[q] = 1
                    if t & 2: v[n + q] = 1
                yield v
    def stats(w):
        vecs = [[sp2(v, R[i], n) for i in range(k)] for v in ops_of_weight(w)]
        cnt = [sum(x) for x in vecs]
        nzv = [x for x in vecs if any(x)]
        return dict(max=max(cnt), mean=sum(cnt) / len(cnt),
                    rank=f2_rank(nzv), distinct=len(set(map(tuple, nzv))), n_ops=len(cnt))
    out = dict(dist1=stats(1), dist2=stats(2))
    rng = np.random.default_rng(seed); mx = []
    for _ in range(nbasis):
        while True:
            A = rng.integers(0, 2, size=(k, k))
            if f2_rank([list(map(int, r)) for r in A]) == k: break
        Rb = []
        for i in range(k):
            v = [0] * (2 * n)
            for j in range(k):
                if A[i, j]: v = [(a + b) % 2 for a, b in zip(v, R[j])]
            Rb.append(v)
        mx.append(max(sum(sp2(v, Rb[i], n) for i in range(k)) for v in ops_of_weight(1)))
    out["dist1_max_bases"] = (min(mx), max(mx))
    out["n_stab"] = f2_rank(S)
    out["code_log2"] = n - out["n_stab"]
    return out

def selfchecks(name, n, S, R, W=None):
    ch = []
    k = len(R)
    ch.append(("records commute with every stabiliser  [clause (ii): [R,H]=0]",
               all(sp2(r, s, n) == 0 for r in R for s in S), ""))
    ch.append(("no record lies in the stabiliser group  [clause (iii): non-trivial]",
               all(not in_span(r, S) for r in R), ""))
    ch.append(("records pairwise commute  [a joint family of bits]",
               all(sp2(R[i], R[j], n) == 0 for i in range(k) for j in range(k) if i != j), ""))
    ch.append(("records independent mod stabilisers",
               f2_rank(list(R) + list(S)) == k + f2_rank(S), ""))
    if W is not None:
        ch.append(("symplectic pairing matrix is the identity (non-degenerate)",
                   all(sp2(R[i], W[j], n) == (1 if i == j else 0)
                       for i in range(k) for j in range(k)), ""))
    return [(name, c, ok, note) for c, ok, note in ch]

# ------------------------------------------------------------------ run
def run():
    P("=" * 118)
    P("LANE_SCALE_C_ACCUM  script 1  --  EXACT combinatorics.  What accumulates with the number of records k?")
    P("=" * 118)
    CHECKS = []

    NS = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22]
    fam = {}
    for n in NS:
        nn, S = carrier_nn2(n)
        pairs = symplectic_logicals([s[:] for s in S], nn)
        Rr = [p[1] for p in pairs]; Wr = [p[0] for p in pairs]
        CHECKS += selfchecks(f"[[{n},{n-2},2]]", nn, S, Rr, Wr)
        CHECKS.append((f"[[{n},{n-2},2]]", "symplectic_logicals returned k = n-2 pairs",
                       len(pairs) == n - 2, f"{len(pairs)}"))
        fam[n] = dict(name=f"[[{n},{n-2},2]]", n=nn, S=S, R=Rr, W=Wr, k=len(Rr))
        fam[n].update(cheap(nn, S, Rr))

    ctl = {}
    n4, S4 = carrier_product(1)
    pairs4 = symplectic_logicals([s[:] for s in S4], 4)
    for m in range(1, 11):
        n, S = carrier_product(m)
        Rr = [embed(pr[1], 4, 4 * b, n) for b in range(m) for pr in pairs4]
        Wr = [embed(pr[0], 4, 4 * b, n) for b in range(m) for pr in pairs4]
        CHECKS += selfchecks(f"[[4,2,2]]^{m}", n, S, Rr, Wr)
        c = dict(name=f"[[4,2,2]]^{m}", n=n, S=S, R=Rr, W=Wr, k=len(Rr))
        c.update(cheap(n, S, Rr, nbasis=12))
        ctl[2 * m] = c

    # ---------------- the invariant weight map, exact where the enumeration fits
    P("")
    P("Building the EXACT weight map w(v) by exhaustive enumeration of the admissible Pauli group.")
    for n in NS:
        t0 = time.time()
        wv, tot, how = weight_map(fam[n]["n"], fam[n]["S"], fam[n]["R"])
        fam[n]["wmap_how"] = how; fam[n]["wmap_size"] = tot
        if how == "exact":
            fam[n]["wmap"] = wmap_summary(wv, fam[n]["k"])
            P(f"   {fam[n]['name']:14s} k={fam[n]['k']:2d}  {tot:>12,d} admissible Paulis"
              f"   {time.time()-t0:6.1f}s   EXACT")
        else:
            fam[n]["wmap"] = None
            P(f"   {fam[n]['name']:14s} k={fam[n]['k']:2d}  {tot:>12,d} admissible Paulis"
              f"   -- {how}, not enumerated")
    for kk in sorted(ctl):
        c = ctl[kk]; t0 = time.time()
        wv, tot, how = weight_map(c["n"], c["S"], c["R"])
        c["wmap_how"] = how; c["wmap_size"] = tot
        if how == "exact":
            c["wmap"] = wmap_summary(wv, c["k"])
            P(f"   {c['name']:14s} k={c['k']:2d}  {tot:>12,d} admissible Paulis"
              f"   {time.time()-t0:6.1f}s   EXACT")
        else:
            c["wmap"] = None
            P(f"   {c['name']:14s} k={c['k']:2d}  {tot:>12,d} admissible Paulis"
              f"   -- {how}, not enumerated")

    # SELF-CHECK / EXTRAPOLATION LICENCE: for a PRODUCT of disjoint blocks the weight map must be
    # exactly additive across blocks, since the supports are disjoint and the constraints decouple.
    # If that holds where both are computable, the control curve may be extended by additivity.
    add_ok, add_note = None, ""
    if ctl[2]["wmap"] and ctl[4]["wmap"] and ctl[6]["wmap"]:
        b = ctl[2]["wmap"]
        pred = [(m, m * b["wmax"], m * b["w_allones"]) for m in (2, 3)]
        got = [(m, ctl[2 * m]["wmap"]["wmax"], ctl[2 * m]["wmap"]["w_allones"]) for m in (2, 3)]
        add_ok = all(p == g for p, g in zip(pred, got))
        add_note = f"predicted {pred} vs measured {got}"
    CHECKS.append(("control", "block additivity of the weight map (licenses extrapolation)",
                   bool(add_ok), add_note))
    # EXTRAPOLATE the control by that verified additivity, and MARK every extrapolated cell.
    if add_ok:
        b = ctl[2]["wmap"]
        for kk in sorted(ctl):
            c = ctl[kk]
            if c.get("wmap") is None:
                m = kk // 2; K = kk
                mean_all = m * b["wmean_all"]
                c["wmap"] = dict(wmax=m * b["wmax"], wmin=b["wmin"],
                                 wmean=mean_all * (2 ** K) / (2 ** K - 1),
                                 wmean_all=mean_all, unreached=0, span_at=b["span_at"],
                                 w_allones=m * b["w_allones"], exact=False)

    # ---------------- self-checks
    P("")
    P("SELF-CHECKS  (a FAILING check voids every conclusion below it)")
    P("-" * 118)
    bad = [c for c in CHECKS if not c[2]]
    for c in bad: P(f"   FAIL  {c[0]:14s} {c[1]}   {c[3]}")
    P(f"   {len(CHECKS) - len(bad)} / {len(CHECKS)} checks pass"
      + ("   -- ALL PASS" if not bad else "   -- SOME FAILED, see above"))

    # ---------------- TABLE 1
    P("")
    P("TABLE 1   BASIS-INVARIANT SCALING.   Control columns carry the SAME k on a product of")
    P("          independent [[4,2,2]] blocks, where every growth is additive BY CONSTRUCTION.")
    P("-" * 118)
    P(f"{'k':>3} | {'n':>3} {'codeS':>5} {'nstab':>5} {'wmax':>4} {'wmean':>6} {'wall':>4} {'span':>4} {'d1rk':>4} {'d2rk':>4}"
      f" || {'nC':>3} {'codeS':>5} {'nstab':>5} {'wmax':>4} {'wmean':>6} {'wall':>4} {'span':>4} {'d1rk':>4} {'d2rk':>4}")
    P("-" * 118)
    def cells(r):
        if r is None:
            return f"{'-':>3} {'-':>5} {'-':>5} {'-':>4} {'-':>6} {'-':>4} {'-':>4} {'-':>4} {'-':>4}"
        wm = r.get("wmap")
        star = "" if (wm and wm.get("exact", True)) else "*"
        a = (f"{str(wm['wmax'])+star:>4} {wm['wmean']:>6.3f} {str(wm['w_allones'])+star:>4} {str(wm['span_at']):>4}"
             if wm else f"{'--':>4} {'--':>6} {'--':>4} {'--':>4}")
        return (f"{r['n']:>3} {r['code_log2']:>5} {r['n_stab']:>5} {a} "
                f"{r['dist1']['rank']:>4} {r['dist2']['rank']:>4}")
    for n in NS:
        r = fam[n]; P(f"{r['k']:>3} | {cells(r)} || {cells(ctl.get(r['k']))}")
    P("-" * 118)
    P("codeS = log2 dim code space (entropy of the maximally mixed code state, in bits)")
    P("nstab = number of independent stabiliser constraints")
    P("wmax  = max over all 2^k-1 non-zero disturbance patterns v of the MIN weight realising v  [INVARIANT]")
    P("wmean = mean of that min weight over all non-zero v                                       [INVARIANT]")
    P("wall  = min weight of one admissible operator that flips ALL k records at once")
    P("span  = smallest weight at which the realised patterns SPAN F_2^k                         [INVARIANT]")
    P("d1rk / d2rk = F_2 rank of the disturbance patterns of ALL weight-1 / weight-2 Paulis      [INVARIANT]")
    P("*     = control value EXTRAPOLATED by block additivity, which was verified exactly at k = 4 and 6")
    P("--    = the admissible Pauli group exceeded the 2e8-operator enumeration budget; no value claimed")

    # ---------------- TABLE 2
    P("")
    P("TABLE 2   BASIS-DEPENDENT quantities, printed so the basis dependence is VISIBLE and")
    P("          cannot be mistaken for a law.  d1max = max records disturbed by one 1-qubit Pauli.")
    P("-" * 118)
    P(f"{'k':>3} | {'fam d1max':>9} {'fam range over 25 bases':>24} {'fam d1mean':>10}"
      f" || {'ctl d1max':>9} {'ctl range over 12 bases':>24} {'ctl d1mean':>10}")
    for n in NS:
        r = fam[n]; c = ctl.get(r["k"])
        cs = (f"{c['dist1']['max']:>9} {str(c['dist1_max_bases']):>24} {c['dist1']['mean']:>10.3f}"
              if c else f"{'-':>9} {'-':>24} {'-':>10}")
        P(f"{r['k']:>3} | {r['dist1']['max']:>9} {str(r['dist1_max_bases']):>24} "
          f"{r['dist1']['mean']:>10.3f} || {cs}")

    # ---------------- SCALING READ
    P("")
    P("SCALING READ   (classifications computed FROM the tables above)")
    P("-" * 118)
    ks = [fam[n]["k"] for n in NS]
    def classify(ks, vals):
        pts = [(k, v) for k, v in zip(ks, vals) if v is not None]
        if len(pts) < 3: return f"too few points ({len(pts)})"
        K = np.array([p[0] for p in pts], float); V = np.array([p[1] for p in pts], float)
        if V.max() == V.min(): return f"CONSTANT at {V[0]:g}  over k = {int(K.min())}..{int(K.max())}"
        A = np.vstack([K, np.ones_like(K)]).T
        sol = np.linalg.lstsq(A, V, rcond=None)[0]
        rms_lin = float(np.sqrt(np.mean((V - A @ sol) ** 2)))
        A2 = np.vstack([K * K, K, np.ones_like(K)]).T
        sol2 = np.linalg.lstsq(A2, V, rcond=None)[0]
        rms_q = float(np.sqrt(np.mean((V - A2 @ sol2) ** 2)))
        A3 = np.vstack([np.log2(K), np.ones_like(K)]).T
        sol3 = np.linalg.lstsq(A3, V, rcond=None)[0]
        rms_log = float(np.sqrt(np.mean((V - A3 @ sol3) ** 2)))
        best = min([("LINEAR", rms_lin), ("QUADRATIC", rms_q), ("LOGARITHMIC", rms_log)],
                   key=lambda t: t[1])
        return (f"{best[0]:11s} slope_lin={sol[0]:+.4f} rms_lin={rms_lin:.4f} "
                f"rms_quad={rms_q:.4f} rms_log={rms_log:.4f}  "
                f"[k x{K[-1]/K[0]:.1f} -> value x{V[-1]/V[0]:.2f}]")
    def ser(d, key, sub=None, wm=None):
        out = []
        for n in NS:
            r = d[n] if d is fam else d.get(fam[n]["k"])
            if r is None: out.append(None); continue
            if wm: out.append(r["wmap"][wm] if r.get("wmap") else None)
            elif sub: out.append(r[key][sub])
            else: out.append(r[key])
        return out
    ROWS = [
        ("code-space entropy (bits)", ser(fam, "code_log2"), ser(ctl, "code_log2")),
        ("independent stabiliser constraints", ser(fam, "n_stab"), ser(ctl, "n_stab")),
        ("wmax  worst-case write cost [INVARIANT]", ser(fam, None, wm="wmax"), ser(ctl, None, wm="wmax")),
        ("wmean mean write cost [INVARIANT]", ser(fam, None, wm="wmean"), ser(ctl, None, wm="wmean")),
        ("wall  cost to flip every record at once", ser(fam, None, wm="w_allones"), ser(ctl, None, wm="w_allones")),
        ("span  weight at which locals span F_2^k", ser(fam, None, wm="span_at"), ser(ctl, None, wm="span_at")),
        ("rank of 1-qubit disturbance patterns", ser(fam, "dist1", "rank"), ser(ctl, "dist1", "rank")),
        ("rank of 2-qubit disturbance patterns", ser(fam, "dist2", "rank"), ser(ctl, "dist2", "rank")),
        ("distinct 1-qubit disturbance patterns", ser(fam, "dist1", "distinct"), ser(ctl, "dist1", "distinct")),
    ]
    for label, a, b in ROWS:
        P(f"  {label}")
        P(f"     family  {a}")
        P(f"             {classify(ks, a)}")
        P(f"     control {b}")
        P(f"             {classify(ks, b)}")
    P("")
    P(f"Largest carrier reached in script 1: n = {NS[-1]} qubits, k = {NS[-1]-2} records, "
      f"Hilbert dim 2^{NS[-1]} = {2**NS[-1]:,d} -- never built; only F_2 vectors of length {2*NS[-1]}.")
    P("Exact weight maps stop where the admissible Pauli group exceeds the enumeration budget "
      "(2e8 operators); that limit, not memory or time on the F_2 side, is what stopped the invariant columns.")
    return fam, ctl

if __name__ == "__main__":
    run()
    with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_C_ACCUM/s1_combinatorics.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
