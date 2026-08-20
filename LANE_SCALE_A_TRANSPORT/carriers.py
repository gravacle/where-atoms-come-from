"""carriers.py -- the reusable construction shared by S1/S2/S3.  No top-level work is done
here, so importing it never re-runs an analysis."""
import sys, time, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_A_TRANSPORT")
import glib
def say(*a): print(*a); sys.stdout.flush()
TOL = 1e-8

def census(G):
    n = G.n
    cl = G.classes
    cent_sz = np.array([len(G.centralizer(g)) for g in range(n)], dtype=np.int64)
    # k(C_G(h)) : number of conjugacy classes of the centraliser, as a subgroup of G
    kcent = []
    for h in range(n):
        C = G.centralizer(h)
        Cs = set(C); seen = set(); kk = 0
        for g in C:
            if g in seen: continue
            orb = {G.conj(x, g) for x in C}
            seen |= orb; kk += 1
        kcent.append(kk)
    kcent = np.array(kcent, dtype=np.int64)
    FixX = cent_sz ** 2
    FixC = cent_sz * kcent
    a = int(FixX.sum()) // n
    m = int(FixC.sum()) // n
    Csz = n * G.k()
    assert int(FixC.sum()) % n == 0 and int(FixX.sum()) % n == 0, "Burnside not integral"
    d2 = m
    d1 = (a - m) + (Csz - m)
    d0 = n * n - Csz - (a - m)
    assert d2 + d1 + d0 == n * n, "multiplicities do not sum to dim"
    chis = {}
    chis[-2] = np.full(n, float(m))
    chis[-1] = (a - m) + FixC.astype(float) - m
    chis[0]  = FixX.astype(float) - FixC.astype(float) - (a - m)
    dims = {-2: d2, -1: d1, 0: d0}
    for v in (-2, -1, 0):
        assert abs(chis[v][G.e] - dims[v]) < 1e-9, f"chi_E(e) != dim E for E={v}"
    return dict(n=n, abelian=G.abelian, dim=n * n, k=G.k(), Z=len(G.centre),
                a=a, m=m, Csz=Csz, dims=dims, chis=chis, cent_sz=cent_sz)

def isotypic(G, ce):
    """(d_rho, m_rho) per eigenspace, EXACT, by characters.  Returns dict E -> list of (d,m)."""
    cl, chi, d, cls_of = G.chars()
    out = {}
    for v in (-2, -1, 0):
        if ce['dims'][v] == 0: out[v] = []; continue
        cf = np.array([ce['chis'][v][c[0]] for c in cl])          # class function
        mult, dd = G.decompose(cf)
        assert np.max(np.abs(mult - np.round(mult))) < 1e-6, "non-integer multiplicity"
        mult = np.round(mult).astype(int)
        assert mult.min() >= 0, "negative multiplicity"
        s = int(np.sum(dd * mult))
        assert s == ce['dims'][v], f"SELF-CHECK FAILED sum d*m = {s} != dim {ce['dims'][v]}"
        out[v] = [(int(dd[r]), int(mult[r])) for r in range(len(mult)) if mult[r] > 0]
    return out

def subset_sums(items, target):
    """can we pick k_rho in [0,m_rho] with sum d_rho k_rho = target?  return a witness or None"""
    reach = {0: []}
    for (d, m) in items:
        nr = {}
        for s, path in reach.items():
            for k in range(m + 1):
                t = s + d * k
                if t <= target and t not in nr: nr[t] = path + [k]
        reach = nr
    return reach.get(target)

def fixed_record_exists(iso, dims):
    """Is there a TRANSPORT-FIXED record?  R must lie in the commutant of every A_h, i.e.
       R|_E = sum_rho I_{d_rho} (x) R_{E,rho}.  Clause (iv) needs Tr(P_E R) = 0, i.e.
       sum_rho d_rho k_rho = dim(E)/2 where k_rho = #(+1) in the multiplicity space.
       Clause (iii) needs R|_E non-scalar on SOME E, i.e. some k_rho not in {0, m_rho}
       or two blocks disagreeing."""
    wit = {}
    for v, items in iso.items():
        D = dims[v]
        if D == 0: continue
        if D % 2: return False, None, f"eigenspace {v} has ODD dim {D}: clause (iv) impossible"
        w = subset_sums(items, D // 2)
        if w is None: return False, None, f"no balanced choice on eigenspace {v}"
        wit[v] = w
    # non-triviality: at least one eigenspace where R|_E is NOT +-I
    nontriv = False
    for v, w in wit.items():
        items = iso[v]
        if len(items) > 1:
            if any(0 < w[i] < items[i][1] for i in range(len(items))): nontriv = True
            # different blocks disagreeing (one full, one empty) also non-scalar
            full = [i for i in range(len(items)) if w[i] == items[i][1]]
            empt = [i for i in range(len(items)) if w[i] == 0]
            if full and empt: nontriv = True
        else:
            if 0 < w[0] < items[0][1]: nontriv = True
    return (True if nontriv else False), wit, ("" if nontriv else "balanced but only SCALAR: clause (iii) fails")

def phi(iso, dims):
    """fraction of the record parameter space that transport FIXES.
       all Hermitian directions in the commutant of H : sum_E dim(E)^2
       transport-fixed ones                          : sum_E sum_rho m_rho^2"""
    allw = sum(dims[v] ** 2 for v in dims)
    fix  = sum(sum(m * m for (_, m) in iso[v]) for v in iso)
    return fix, allw, fix / allw


# ---------------------------------------------------------------- carrier builders
def minimal_torus(G):
    """1 vertex, 2 edges, 1 face.  Returns H, list of transport permutations (one per h), dim."""
    n = G.n; D = n * n
    def num(g1, g2): return g1 * n + g2
    perms = []
    for h in range(n):
        p = np.zeros(D, dtype=np.int64)
        for g1 in range(n):
            c1 = G.conj(h, g1)
            for g2 in range(n):
                p[num(g1, g2)] = num(c1, G.conj(h, g2))
        perms.append(p)
    A = np.zeros((D, D))
    for p in perms:
        A[p, np.arange(D)] += 1.0 / n
    B = np.zeros((D, D))
    for g1 in range(n):
        for g2 in range(n):
            if G.mt[g1, g2] == G.mt[g2, g1]: B[num(g1, g2), num(g1, g2)] = 1.0
    H = -(A + B)
    return H, perms, D

def torus_1x2(G):
    """L1=1, L2=2 square-lattice torus: V=2, E=4, F=2, as in LANE_O35.  H = -(sum_v A_v + B).
       Transport = GLOBAL conjugation on all four edges (the diagonal gauge action)."""
    n = G.n; nE = 4; D = n ** nE
    # edges: h(0,0)=0, v(0,0)=1, h(0,1)=2, v(0,1)=3
    eid = {(0, 0, 0): 0, (0, 0, 1): 1, (0, 1, 0): 2, (0, 1, 1): 3}
    faces = []
    for c in (0, 1):
        faces.append([(eid[(0, c, 0)], +1), (eid[(0, (c + 1) % 2, 1)], +1),
                      (eid[(0, c, 0)], -1), (eid[(0, c, 1)], -1)])
    verts = []
    for c in (0, 1):
        verts.append(([eid[(0, c, 0)], eid[(0, c, 1)]],
                      [eid[(0, (c - 1) % 2, 0)], eid[(0, c, 1)]]))
    def cfg(i):
        out = []
        for _ in range(nE): out.append(i % n); i //= n
        return tuple(reversed(out))
    def num(t):
        i = 0
        for g in t: i = i * n + g
        return i
    cfgs = [cfg(i) for i in range(D)]
    B = np.zeros((D, D))
    for i, t in enumerate(cfgs):
        ok = True
        for f in faces:
            p = G.e
            for (k, s) in f: p = G.mt[p, t[k] if s > 0 else G.iv[t[k]]]
            if p != G.e: ok = False; break
        if ok: B[i, i] = 1.0
    Asum = np.zeros((D, D))
    for (outs, ins) in verts:
        for h in range(n):
            for i, t in enumerate(cfgs):
                tt = list(t)
                for k in outs: tt[k] = G.mt[h, tt[k]]
                for k in ins:  tt[k] = G.mt[tt[k], G.iv[h]]
                Asum[num(tuple(tt)), i] += 1.0 / n
    perms = []
    for h in range(n):
        p = np.zeros(D, dtype=np.int64)
        for i, t in enumerate(cfgs):
            p[i] = num(tuple(G.conj(h, x) for x in t))
        perms.append(p)
    H = -(Asum + B)
    return H, perms, D

# ---------------------------------------------------------------- eigenspaces / records
def eigblocks(H, tol=1e-8):
    w, V = np.linalg.eigh(H)
    out = []; i = 0
    while i < len(w):
        j = i
        while j + 1 < len(w) and abs(w[j + 1] - w[i]) < tol: j += 1
        out.append((float(w[i]), V[:, i:j + 1]))
        i = j + 1
    return out

def perm_apply(p, M):
    """A_h M A_h^dagger for a permutation A_h (A_h e_i = e_{p[i]})."""
    return M[np.ix_(p, p)]

def check_clauses(R, blocks, tol=1e-7):
    n = R.shape[0]
    ok_i = np.linalg.norm(R - R.conj().T) < tol and np.linalg.norm(R @ R - np.eye(n)) < tol
    ok_iii = False; ok_iv = True
    for _, Q in blocks:
        S = Q.conj().T @ R @ Q
        m = S.shape[0]
        if np.linalg.norm(S - (np.trace(S) / m) * np.eye(m)) > tol: ok_iii = True
        if abs(np.trace(S)) > tol * max(1.0, m): ok_iv = False
    return ok_i, ok_iii, ok_iv

def generic_record(blocks, rng, n):
    R = np.zeros((n, n))
    for _, Q in blocks:
        m = Q.shape[1]
        if m % 2: return None
        A = rng.normal(size=(m, m)); Qo, _ = np.linalg.qr(A)
        s = np.array([1.0] * (m // 2) + [-1.0] * (m // 2))
        R += Q @ (Qo * s) @ Qo.conj().T @ Q.conj().T
    return R

def gauge_record(blocks, perms, rng, n, tol=1e-7):
    """A record inside the commutant of EVERY A_h: average a random Hermitian over the
       transport group, split into minimal projections, pick a half-dimensional subset."""
    R = np.zeros((n, n))
    for _, Q in blocks:
        m = Q.shape[1]
        if m % 2: return None
        M = rng.normal(size=(n, n)); M = (M + M.T) / 2
        Mb = sum(perm_apply(p, M) for p in perms) / len(perms)
        S = Q.conj().T @ Mb @ Q
        S = (S + S.conj().T) / 2
        w, V = np.linalg.eigh(S)
        projs = []; i = 0
        while i < m:
            j = i
            while j + 1 < m and abs(w[j + 1] - w[i]) < 1e-7: j += 1
            projs.append(V[:, i:j + 1]); i = j + 1
        ranks = [p.shape[1] for p in projs]
        sel = subset_sums([(r, 1) for r in ranks], m // 2)
        if sel is None: return None
        P = np.zeros((m, m))
        for idx, take in enumerate(sel):
            if take: P += projs[idx] @ projs[idx].conj().T
        R += Q @ (2 * P - np.eye(m)) @ Q.conj().T
    return R

def moved(R, perms, tol=1e-7):
    best = 0.0
    for p in perms:
        c = perm_apply(p, R) - R
        best = max(best, float(np.linalg.norm(c)))
    return best > tol, best

# ---------------------------------------------------------------- run
