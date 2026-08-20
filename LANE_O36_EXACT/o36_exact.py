"""
O-36 DECIDED EXACTLY: does a TRANSPORT-FIXED record exist on D(D_4)?

R transport-fixed record  <=>  R = R-dag, R^2 = I, [H,R] = 0, [A_h,R] = 0 for all h in G,
                               and Tr(P_E R) = 0 on EVERY H-eigenspace E.

Since [A_h,H] = 0, every H-eigenspace E is a G-representation under h -> A_h|_E.  R|_E must lie
in the COMMUTANT of that representation, which is exactly  (+)_i  I_{d_i} (x) M_{m_i}.
A self-adjoint unitary there has  Tr = sum_i d_i (2 p_i - m_i),  0 <= p_i <= m_i integer.
So the whole question is an exact integer feasibility problem ONCE (d_i, m_i) are exact.

(d_i, m_i) ARE OBTAINED BY CHARACTER THEORY, NOT BY SAMPLING.  record_model.commutant() is
sampling-based; an incomplete basis inflates multiplicities (the previous attempt returned
sum d_i m_i = 20 against dim E = 18).  Here:
  * the character table of G is computed from the CLASS ALGEBRA (Burnside/Dixon) -- exact up to
    linear algebra on integer matrices -- and then VERIFIED against both orthogonality relations,
    sum d_i^2 = |G|, integrality of d_i, and the regular character.
  * chi_E(h) = Tr(P_E A_h) directly.
  * m_i = (1/|G|) sum_h conj(chi_i(h)) chi_E(h), CHECKED to be a non-negative integer.
  * ASSERT sum_i d_i m_i = dim E.
"""
import sys, itertools
import numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import RecordModel, eigenspaces

np.set_printoptions(linewidth=200)
TOL = 1e-8

# ------------------------------------------------------------------ finite groups
def group_Zn(n):
    els = list(range(n))
    return els, (lambda a, b: (a + b) % n), 0, "Z_%d" % n

def group_ZmxZn(m, n):
    els = [(a, b) for a in range(m) for b in range(n)]
    return els, (lambda x, y: ((x[0] + y[0]) % m, (x[1] + y[1]) % n)), (0, 0), "Z_%dxZ_%d" % (m, n)

def group_dihedral(k):
    """D_k = <r,s | r^k = s^2 = e, s r s = r^-1>, order 2k.  Elements (a,b) = r^a s^b."""
    els = [(a, b) for a in range(k) for b in range(2)]
    def mul(x, y):
        a1, b1 = x; a2, b2 = y
        return ((a1 + (a2 if b1 == 0 else -a2)) % k, (b1 + b2) % 2)
    return els, mul, (0, 0), "D_%d" % k

def group_Q8():
    """Quaternion group as {+-1,+-i,+-j,+-k} encoded (sign, unit) with unit in 0..3."""
    els = [(s, u) for s in (0, 1) for u in range(4)]          # s=0 -> +, s=1 -> -
    # unit multiplication table for 1,i,j,k  -> (sign, unit)
    T = {}
    names = [0, 1, 2, 3]
    base = {(0, 0): (0, 0), (0, 1): (0, 1), (0, 2): (0, 2), (0, 3): (0, 3),
            (1, 0): (0, 1), (2, 0): (0, 2), (3, 0): (0, 3),
            (1, 1): (1, 0), (2, 2): (1, 0), (3, 3): (1, 0),
            (1, 2): (0, 3), (2, 1): (1, 3),
            (2, 3): (0, 1), (3, 2): (1, 1),
            (3, 1): (0, 2), (1, 3): (1, 2)}
    def mul(x, y):
        s1, u1 = x; s2, u2 = y
        s3, u3 = base[(u1, u2)]
        return ((s1 + s2 + s3) % 2, u3)
    return els, mul, (0, 0), "Q_8"

# ------------------------------------------------------------------ exact character table
def conjugacy_classes(els, mul):
    idx = {g: i for i, g in enumerate(els)}
    n = len(els)
    inv = [None] * n
    e = None
    for i, g in enumerate(els):
        for j, h in enumerate(els):
            if mul(g, h) == mul(h, g) and mul(g, h) == g and mul(g, h) == h:
                pass
    # identity: the unique x with x*x = x
    for i, g in enumerate(els):
        if mul(g, g) == g: e = i
    for i, g in enumerate(els):
        for j, h in enumerate(els):
            if idx[mul(g, h)] == e: inv[i] = j
    seen = [False] * n; classes = []
    for i in range(n):
        if seen[i]: continue
        cl = set()
        for j in range(n):
            cl.add(idx[mul(mul(els[j], els[i]), els[inv[j]])])
        for c in cl: seen[c] = True
        classes.append(sorted(cl))
    # put the identity class first
    classes.sort(key=lambda c: (e not in c, len(c), c))
    return classes, e, inv, idx

def character_table(els, mul):
    """Burnside/Dixon: simultaneously diagonalise the class-algebra matrices.
       M_i[j][k] = a_{i j k}  where  C_i C_j = sum_k a_{i j k} C_k."""
    n = len(els)
    classes, e, inv, idx = conjugacy_classes(els, mul)
    k = len(classes)
    cls_of = [0] * n
    for ci, c in enumerate(classes):
        for x in c: cls_of[x] = ci
    # structure constants: a_{i j k} = # of (x in C_i, y in C_j) with x*y = a fixed rep of C_k
    a = np.zeros((k, k, k), dtype=float)
    reps = [c[0] for c in classes]
    for i in range(k):
        for j in range(k):
            cnt = np.zeros(k)
            tally = {}
            for x in classes[i]:
                for y in classes[j]:
                    z = idx[mul(els[x], els[y])]
                    tally[z] = tally.get(z, 0) + 1
            for kk in range(k):
                a[i, j, kk] = tally.get(reps[kk], 0)
    M = [a[i] for i in range(k)]                       # M[i][j,kk] = a_{i j kk}
    rng = np.random.default_rng(7)
    for attempt in range(50):
        c = rng.normal(size=k)
        Msum = sum(c[i] * M[i] for i in range(k))
        w, V = np.linalg.eig(Msum)
        if len(w) == k and np.min(np.abs(w[:, None] - w[None, :] + np.eye(k) * 1e6)) > 1e-6:
            break
    else:
        raise RuntimeError("could not separate the class algebra")
    chars = []
    sizes = np.array([len(c) for c in classes], dtype=float)
    for r in range(k):
        v = V[:, r]
        omega = np.array([ (np.vdot(v, M[i] @ v) / np.vdot(v, v)) for i in range(k) ])
        d2 = len(els) / np.sum(np.abs(omega) ** 2 / sizes)
        d = np.sqrt(d2)
        chi = omega * d / sizes
        chars.append((d, chi))
    chars.sort(key=lambda t: (round(t[0].real), -np.sum(np.abs(t[1]))))
    dims = np.array([c[0] for c in chars])
    table = np.array([c[1] for c in chars])
    return classes, sizes, dims, table, cls_of, e

def verify_character_table(name, els, classes, sizes, dims, table, out):
    ok = True
    G = len(els)
    di = np.round(dims.real).astype(int)
    e1 = float(np.max(np.abs(dims - di)))
    ok &= e1 < 1e-6
    e2 = abs(int((di ** 2).sum()) - G)
    ok &= (e2 == 0)
    # row orthonormality  (1/|G|) sum_c |C| conj(chi_i) chi_j = delta_ij
    Grm = (table.conj() * sizes) @ table.T / G
    e3 = float(np.max(np.abs(Grm - np.eye(len(di)))))
    ok &= e3 < 1e-8
    # column orthogonality  sum_i conj(chi_i(c)) chi_i(c') = delta * |G|/|C|
    Col = table.conj().T @ table
    e4 = float(np.max(np.abs(Col - np.diag(G / sizes))))
    ok &= e4 < 1e-8
    # regular character decomposes with multiplicity d_i
    chi_reg = np.zeros(len(sizes)); chi_reg[0] = G
    mreg = (table.conj() * sizes) @ chi_reg / G
    e5 = float(np.max(np.abs(mreg - di)))
    ok &= e5 < 1e-8
    # chi_i(e) = d_i
    e6 = float(np.max(np.abs(table[:, 0] - dims)))
    ok &= e6 < 1e-9
    out.append("  char-table self-check [%s]: dims=%s  int_err=%.1e  sum d^2 - |G| = %d  "
               "row_orth_err=%.1e  col_orth_err=%.1e  regular_err=%.1e  chi(e)-d err=%.1e  -> %s"
               % (name, list(di), e1, e2, e3, e4, e5, e6, "PASS" if ok else "FAIL"))
    return ok, di

# ------------------------------------------------------------------ D(G) on the minimal torus
def build_DG(els, mul, idx, inv):
    n = len(els)
    N = n * n
    def A_of(h):
        M = np.zeros((N, N), dtype=complex)
        hi = idx[h]
        for i1, g1 in enumerate(els):
            for i2, g2 in enumerate(els):
                a = mul(mul(h, g1), els[inv[hi]])
                b = mul(mul(h, g2), els[inv[hi]])
                M[idx[a] * n + idx[b], i1 * n + i2] = 1.0
        return M
    As = [A_of(h) for h in els]
    A = sum(As) / n
    Bd = np.zeros(N)
    for i1, g1 in enumerate(els):
        for i2, g2 in enumerate(els):
            c = mul(mul(g1, g2), mul(els[inv[i1]], els[inv[i2]]))
            if idx[c] == idx[mul(g1, els[inv[i1]])]:      # c == e
                Bd[i1 * n + i2] = 1.0
    B = np.diag(Bd).astype(complex)
    H = -(A + B)
    return H, A, B, As

# ------------------------------------------------------------------ integer feasibility
def achievable_traces(blocks):
    """blocks = [(d_i, m_i)]. Achievable Tr = sum_i d_i (2 p_i - m_i). Returns set + a witness map."""
    reach = {0: []}
    for (d, m) in blocks:
        nxt = {}
        for s, path in reach.items():
            for p in range(m + 1):
                t = s + d * (2 * p - m)
                if t not in nxt: nxt[t] = path + [p]
        reach = nxt
    return reach

# ------------------------------------------------------------------ per-carrier pipeline
def analyse(gspec, out):
    els, mul, e_el, name = gspec
    classes, sizes, dims, table, cls_of, e = character_table(els, mul)
    ok, di = verify_character_table(name, els, classes, sizes, dims, table, out)
    if not ok:
        out.append("  ABORT %s: character table failed self-check; no conclusion drawn." % name)
        return None
    _, e2, inv, idx = conjugacy_classes(els, mul)
    H, A, B, As = build_DG(els, mul, idx, inv)
    n = len(els); N = n * n
    # sanity: A_h is a unitary rep of G
    rep_err = 0.0
    for i in range(min(n, 8)):
        for j in range(min(n, 8)):
            rep_err = max(rep_err, np.linalg.norm(As[i] @ As[j] - As[idx[mul(els[i], els[j])]]))
    comm_err = max(np.linalg.norm(As[i] @ H - H @ As[i]) for i in range(n))
    es = eigenspaces(H)
    out.append("  carrier D(%s): dim %d   rep-homomorphism err %.2e   max||[A_h,H]|| %.2e"
               % (name, N, rep_err, comm_err))
    out.append("  eigenvalues/multiplicities: %s" % [(round(float(v), 6), m) for v, _, m in es])
    per_E = []
    all_ok = True
    for (val, P, m) in es:
        chiE = np.array([np.trace(P @ As[i]) for i in range(n)])
        imag = float(np.max(np.abs(chiE.imag)))
        # multiplicities
        mult = []
        merr = 0.0
        for r in range(len(di)):
            s = 0j
            for i in range(n):
                s += np.conj(table[r, cls_of[i]]) * chiE[i]
            s /= n
            mr = float(s.real)
            merr = max(merr, abs(s.imag), abs(mr - round(mr)))
            mult.append(int(round(mr)))
        neg = any(x < 0 for x in mult)
        tot = sum(int(di[r]) * mult[r] for r in range(len(di)))
        consistent = (tot == m) and (not neg) and merr < 1e-6
        all_ok &= consistent
        blocks = [(int(di[r]), mult[r]) for r in range(len(di)) if mult[r] > 0]
        reach = achievable_traces(blocks) if blocks else {0: []}
        zero_ok = (0 in reach)
        per_E.append(dict(val=float(val), dim=int(m), mult=mult, blocks=blocks,
                          tot=tot, merr=merr, chiE_imag=imag, consistent=consistent,
                          zero=zero_ok, witness=reach.get(0), P=P,
                          tracespan=(min(reach), max(reach))))
    return dict(name=name, N=N, H=H, As=As, es=es, di=list(map(int, di)), per_E=per_E,
                table_ok=ok, blocks_ok=all_ok, abelian=all(len(c) == 1 for c in classes),
                idx=idx, els=els, mul=mul, inv=inv)

# ------------------------------------------------------------------ direct construction
def gauge_average(X, As):
    return sum(g @ X @ g.conj().T for g in As) / len(As)

def block_project(X, es):
    return sum(P @ X @ P for _, P, _ in es)

def minimal_projections(X, tol=1e-7):
    w, V = np.linalg.eigh(X)
    out, i = [], 0
    while i < len(w):
        j = i
        while j + 1 < len(w) and abs(w[j + 1] - w[i]) < tol: j += 1
        Q = V[:, i:j + 1]; out.append(Q @ Q.conj().T); i = j + 1
    return out


def eigenbasis_of(P, tol=1e-8):
    """Orthonormal columns spanning the range of the projector P -- exact, no big offsets."""
    w, V = np.linalg.eigh(P)
    return V[:, w > 0.5]

def minimal_projections_in(X, Q, tol=1e-7):
    """Minimal projections of X restricted to the subspace spanned by the orthonormal Q.
       Diagonalising the SMALL matrix Q-dag X Q avoids the 1e6-offset trick, whose dynamic range
       cost ~6 digits of eigenvector accuracy."""
    S = Q.conj().T @ X @ Q
    S = (S + S.conj().T) / 2
    w, V = np.linalg.eigh(S)
    out, i = [], 0
    while i < len(w):
        j = i
        while j + 1 < len(w) and abs(w[j + 1] - w[i]) < tol: j += 1
        W = Q @ V[:, i:j + 1]
        out.append(W @ W.conj().T); i = j + 1
    return out

def signs_to_zero(ranks):
    """choose s_a in {+1,-1} with sum s_a * rank_a = 0; DP. Returns list or None."""
    reach = {0: []}
    for r in ranks:
        nxt = {}
        for s, path in reach.items():
            for sg in (+1, -1):
                t = s + sg * r
                if t not in nxt: nxt[t] = path + [sg]
        reach = nxt
    return reach.get(0)

def construct_fixed_record(res, seed=0):
    """Project a random Hermitian onto commutant({H} u {A_h}) by EXACT group averaging, then
       look for a zero-trace self-adjoint unitary in it."""
    N = res["N"]; As = res["As"]; es = res["es"]
    rng = np.random.default_rng(seed)
    Y = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    Y = (Y + Y.conj().T) / 2
    X = gauge_average(block_project(Y, es), As)          # the two projections commute
    # verify X is really in the commutant
    err_H = np.linalg.norm(X @ res["H"] - res["H"] @ X)
    err_A = max(np.linalg.norm(X @ g - g @ X) for g in As)
    R = np.zeros((N, N), dtype=complex)
    detail = []
    feasible = True
    for (val, P, m) in es:
        Q = eigenbasis_of(P)
        projs = minimal_projections_in(X, Q)
        ranks = [int(round(float(np.trace(q).real))) for q in projs]
        span_ok = (sum(ranks) == m)
        sg = signs_to_zero(ranks)
        detail.append((float(val), int(m), sorted(ranks), sum(ranks), span_ok, sg is not None))
        if sg is None or not span_ok:
            feasible = False
            R = R + P            # placeholder; verdict already failed
        else:
            R = R + sum(s * q for s, q in zip(sg, projs))
    return dict(err_H=err_H, err_A=err_A, detail=detail, feasible=feasible, R=R)

def verify_record(R, H, es, As):
    N = R.shape[0]
    v = {}
    v["herm"] = float(np.linalg.norm(R - R.conj().T))
    v["invol"] = float(np.linalg.norm(R @ R - np.eye(N)))
    v["commH"] = float(np.linalg.norm(R @ H - H @ R))
    v["trace"] = float(max(abs(np.trace(P @ R)) for _, P, _ in es))
    v["commA"] = float(max(np.linalg.norm(R @ g - g @ R) for g in As))
    # clause (iii): non-constant on some eigenspace
    nc = 0.0
    for _, P, m in es:
        M = P @ R @ P
        nc = max(nc, float(np.linalg.norm(M - (np.trace(M) / m) * P)))
    v["nonconst"] = nc
    return v

def blind_record(res, seed=1):
    """POSITIVE CONTROL for D-15: a record built in the commutant of H ALONE (no transport
       constraint). It satisfies (i)-(iv) and should register NON-ZERO ||[A_h,R]|| wherever
       transport acts non-trivially."""
    N = res["N"]; es = res["es"]
    rng = np.random.default_rng(seed)
    Y = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    Y = (Y + Y.conj().T) / 2
    X = block_project(Y, es)
    R = np.zeros((N, N), dtype=complex)
    for (val, P, m) in es:
        projs = minimal_projections_in(X, eigenbasis_of(P))
        ranks = [int(round(float(np.trace(q).real))) for q in projs]
        sg = signs_to_zero(ranks)
        if sg is None: return None
        R = R + sum(s * q for s, q in zip(sg, projs))
    return R

# ------------------------------------------------------------------ main
def main():
    out = []
    out.append("=" * 118)
    out.append("O-36 EXACT: does a TRANSPORT-FIXED record exist on D(D_4)?")
    out.append("method: character theory (exact) for the block data (d_i,m_i); integer feasibility for Tr = 0;")
    out.append("        cross-checked against a DIRECT construction by exact finite-group averaging.")
    out.append("=" * 118)

    carriers = [
        group_dihedral(4),          # THE TARGET: D_4, non-abelian, order 8
        group_Zn(2),                # CONTROL 1 (D-15): C-43 measured 0/40 moved -> fixed record MUST exist
        group_ZmxZn(2, 2),          # CONTROL 2: abelian, order 4
        group_Zn(4),                # CONTROL 3: abelian, order 4
        group_Q8(),                 # D-17 scale: a DIFFERENT non-abelian group of the same order 8
        group_dihedral(8),          # D-17 scale: non-abelian, order 16 (dim 256)
        group_Zn(3),                # NEGATIVE CONTROL: |G| not a power of 2 -> C-41 says NO record at all
    ]

    out.append("")
    out.append("STEP 1 -- CHARACTER TABLES, VERIFIED (no sampling anywhere)")
    results = []
    for g in carriers:
        r = analyse(g, out)
        results.append(r)
        out.append("")

    out.append("=" * 118)
    out.append("STEP 2 -- PER-EIGENSPACE BLOCK DATA AND EXACT INTEGER FEASIBILITY OF Tr(P_E R) = 0")
    out.append("  columns: carrier | ab? | eig | dimE | irrep multiplicities m_i | blocks (d_i x m_i) | "
               "sum d_i m_i == dimE | Tr range | 0 reachable")
    out.append("-" * 118)
    for r in results:
        if r is None: continue
        for k, pe in enumerate(r["per_E"]):
            out.append("  %-10s %-4s %+6.2f  %5d  %-22s %-26s  %-14s %-14s %s"
                       % (("D(%s)" % r["name"]) if k == 0 else "",
                          ("yes" if r["abelian"] else "NO") if k == 0 else "",
                          pe["val"], pe["dim"], str(pe["mult"]),
                          " ".join("%dx%d" % (d, m) for d, m in pe["blocks"]),
                          "%d==%d %s" % (pe["tot"], pe["dim"], "OK" if pe["consistent"] else "**FAIL**"),
                          "[%d,%d]" % pe["tracespan"],
                          "YES" if pe["zero"] else "NO"))
        out.append("-" * 118)

    out.append("")
    out.append("=" * 118)
    out.append("STEP 3 -- VERDICT PER CARRIER (character-theoretic) vs DIRECT CONSTRUCTION (exact averaging)")
    out.append("  a transport-fixed record exists  <=>  0 is reachable on EVERY eigenspace")
    out.append("  columns: carrier | consistency | char verdict | constructed? | ||[A_h,R]|| for the CONSTRUCTED R"
               " | ||[A_h,R]|| for a TRANSPORT-BLIND record (positive control, D-15)")
    out.append("-" * 118)
    summary = []
    for r in results:
        if r is None: continue
        cons = all(pe["consistent"] for pe in r["per_E"])
        verdict = all(pe["zero"] for pe in r["per_E"]) if cons else None
        con = construct_fixed_record(r)
        if con["feasible"]:
            v = verify_record(con["R"], r["H"], r["es"], r["As"])
            built = "YES"
            cA = v["commA"]
            vs = ("herm %.1e invol %.1e [H,R] %.1e maxTr %.1e nonconst %.3f"
                  % (v["herm"], v["invol"], v["commH"], v["trace"], v["nonconst"]))
        else:
            built = "NO"; cA = float("nan"); vs = "-"
        Rb = blind_record(r)
        if Rb is None:
            cb = float("nan"); vb = "no zero-trace splitting exists even without transport"
        else:
            vv = verify_record(Rb, r["H"], r["es"], r["As"])
            cb = vv["commA"]
            vb = ("herm %.1e invol %.1e [H,R] %.1e maxTr %.1e" %
                  (vv["herm"], vv["invol"], vv["commH"], vv["trace"]))
        agree = ((verdict is True and built == "YES") or (verdict is False and built == "NO"))
        out.append("  %-12s %-11s %-14s %-12s %-14s %-14s  agree=%s"
                   % ("D(%s)" % r["name"], "OK" if cons else "**FAIL**",
                      ("EXISTS" if verdict else "NONE") if cons else "n/a",
                      built,
                      ("%.3e" % cA) if cA == cA else "n/a",
                      ("%.3e" % cb) if cb == cb else "n/a",
                      agree))
        out.append("       constructed-R checks : %s" % vs)
        out.append("       blind-R (control)    : %s" % vb)
        out.append("       commutant-projection residuals: ||[X,H]||=%.2e  max||[X,A_h]||=%.2e"
                   % (con["err_H"], con["err_A"]))
        for d in con["detail"]:
            out.append("       eig %+6.2f dimE %4d  minimal-projection ranks %s  sum=%d spans_E=%s zero-trace=%s"
                       % (d[0], d[1], d[2], d[3], d[4], d[5]))
        out.append("-" * 118)
        summary.append(dict(carrier="D(%s)" % r["name"], dim=r["N"], abelian=r["abelian"],
                            consistent=cons, verdict=verdict, built=(built == "YES"),
                            commA=cA, commA_blind=cb,
                            per_E=[(pe["val"], pe["dim"], pe["blocks"], pe["zero"]) for pe in r["per_E"]]))

    out.append("")
    out.append("=" * 118)
    out.append("READ (filled from the numbers above, not in advance)")
    tgt = [s for s in summary if s["carrier"] == "D(D_4)"]
    z2 = [s for s in summary if s["carrier"] == "D(Z_2)"]
    if tgt and z2:
        t = tgt[0]; c = z2[0]
        out.append("  D(D_4):  consistency %s   character verdict: transport-fixed record %s   direct construction: %s"
                   % ("PASSED" if t["consistent"] else "FAILED",
                      "EXISTS" if t["verdict"] else "DOES NOT EXIST",
                      "SUCCEEDED" if t["built"] else "FAILED"))
        out.append("  D(Z_2) control: verdict %s, construction %s  (C-43 measured 0/40 moved, so anything"
                   " but EXISTS here would mean the method is broken)"
                   % ("EXISTS" if c["verdict"] else "NONE", "SUCCEEDED" if c["built"] else "FAILED"))
    out.append("=" * 118)

    txt = "\n".join(out)
    print(txt)
    with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O36_EXACT/o36_exact.txt", "w") as f:
        f.write(txt + "\n")

    import json
    with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O36_EXACT/o36_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

if __name__ == "__main__":
    main()
