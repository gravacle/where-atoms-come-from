"""
LANE G — GROUP AXIS REFUTER — library, written from scratch (numpy only).
No corpus code reused. Every convention published here.

PUBLISHED CONVENTIONS
---------------------
CARRIER K1 (S1 sec.1, re-encoded here, not imported as code):
  vertices  v0..v4 indexed 0..4;  root r = v0
  edges (source, target), index 0..5:
      e1 = (0,1)  e2 = (1,2)  e3 = (2,0)      -- boundary of the filled face F
      e4 = (0,3)  e5 = (3,4)  e6 = (4,0)      -- the unfilled triangle
  faces: F attached along the closed edge path e1.e2.e3   (2-cell index 0)

  d1[v,e] = +1 if v is the target of e, -1 if v is the source of e, 0 otherwise
            (so d1 is 5x6, columns are edge boundaries)
  d2[e,f] = coefficient of e in the attaching 1-chain of f   (6x1)

FIBRE AND CONNECTION (generalised from S1 sec.3, which fixes rank one):
  fibre L_v = C^n at every vertex, same n at every vertex
  connection U_e in U(n) on each edge; transport along e:u->v is z |-> U_e z,
  reverse traversal transports by U_e^{-1} = U_e^dagger

LOOPS:
  gamma_F = e1,e2,e3 traversed forward, based at v0 -> vertex set {0,1,2}
  gamma_C = e4,e5,e6 traversed forward, based at v0 -> vertex set {0,3,4}
  Based holonomy at a vertex w on the loop: start at w, traverse the loop's
  cyclic order once, return to w.  W^(v0) = U_e3 U_e2 U_e1 (rightmost acts first).
  At rank one all based holonomies of one loop coincide; at rank > 1 they are
  conjugate.  Using the BASED holonomy at each vertex is the unique gauge-covariant
  extension of the corpus operator "multiply by W(gamma) at vertices on gamma"
  (S4 choice C11 / W-01).  The naive alternative (use v0's holonomy everywhere) is
  implemented too, as NAIVE, and is shown non-gauge-invariant in g2.

LOOP TRANSPORT OPERATOR on Gamma(L) = direct sum over v of C^n:
  (M_gamma s)(v) = W^(v)(gamma) s(v)   if v on gamma
                 = s(v)                otherwise

BRANCH OVERLAP (S3 sec.4.1 / sec.4.2, generalised):
  Z_k = < M_F^k s , M_C^k s >   with <a,b> = a^dagger b
      = sum_v  s(v)^dagger (A_v)^k (B_v)^k s(v),
        A_v = (W_F^(v))^dagger if v on gamma_F else I
        B_v =  W_C^(v)         if v on gamma_C else I

  At n=1 with W_F=e^{if}, W_C=e^{ic} this is exactly
  Z_k = (uv)^k p0 + u^k (p1+p2) + v^k (p3+p4),  u = conj(W_F), v = W_C.

SCHEDULES (S3 sec.3.5):
  schedule A (uniform):   k_n = 1,  lambda_A = log|Z_1|
  schedule B (canonical): k_n = n,  lambda_B = lim (1/N) sum_{n=1..N} log|Z_n|

SEEDS: every random draw uses numpy.random.default_rng(seed) with the seed printed.
"""

import numpy as np

# ----------------------------------------------------------------- carrier ---

EDGES_K1 = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]
FACES_K1 = [[(0, +1), (1, +1), (2, +1)]]          # face 0 = F, along e1.e2.e3
LOOP_F = [(0, +1), (1, +1), (2, +1)]              # edge index, direction
LOOP_C = [(3, +1), (4, +1), (5, +1)]
NV_K1 = 5


def incidence(nv, edges, faces):
    d1 = np.zeros((nv, len(edges)), dtype=int)
    for j, (s, t) in enumerate(edges):
        d1[s, j] -= 1
        d1[t, j] += 1
    d2 = np.zeros((len(edges), len(faces)), dtype=int)
    for f, path in enumerate(faces):
        for (e, sgn) in path:
            d2[e, f] += sgn
    return d1, d2


def betti(nv, edges, faces):
    d1, d2 = incidence(nv, edges, faces)
    r1 = np.linalg.matrix_rank(d1.astype(float)) if d1.size else 0
    r2 = np.linalg.matrix_rank(d2.astype(float)) if d2.size else 0
    b0 = nv - r1
    b1 = len(edges) - r1 - r2
    b2 = len(faces) - r2
    return dict(V=nv, E=len(edges), F=len(faces), rank_d1=int(r1), rank_d2=int(r2),
                b0=int(b0), b1=int(b1), b2=int(b2),
                chi=nv - len(edges) + len(faces))


def loop_vertices(edges, loop):
    """ordered vertex sequence of the loop, starting at its base."""
    vs = []
    for (e, sgn) in loop:
        s, t = edges[e]
        vs.append(s if sgn > 0 else t)
    return vs


def based_holonomies(U, edges, loop):
    """dict vertex -> based holonomy W^(v) of the loop, as an n x n unitary.
    U is a list of n x n unitaries indexed by edge."""
    vs = loop_vertices(edges, loop)
    L = len(loop)
    out = {}
    for start in range(L):
        n = U[0].shape[0]
        W = np.eye(n, dtype=complex)
        for j in range(L):
            e, sgn = loop[(start + j) % L]
            Ue = U[e] if sgn > 0 else U[e].conj().T
            W = Ue @ W                      # rightmost acts first
        out[vs[start]] = W
    return out


# ------------------------------------------------------- character extraction -

def _eig_unitary(M):
    """eigen-decomposition of a unitary M = P diag(a) P^dagger with P unitary."""
    n = M.shape[0]
    if np.allclose(M, np.eye(n), atol=1e-13):
        return None, None                 # signal: identity, basis is free
    a, P = np.linalg.eig(M)
    # re-orthonormalise (numpy eig on unitary gives a numerically-unitary P but
    # degenerate eigenvalues can spoil it); use QR on the eigenvector blocks.
    P, _ = np.linalg.qr(P)
    a = np.diag(P.conj().T @ M @ P)
    return a, P


def characters(U, edges, loopF, loopC, s, naive=False):
    """Return (zeta, coeff) : Z_k = sum_j coeff[j] * zeta[j]**k, exactly.

    s is a list of length nv of C^n vectors (the ready state, sum |s(v)|^2 = 1).
    Derivation:  s^d A^k B^k s  with A = P diag(al) P^d, B = Q diag(be) Q^d gives
        sum_{m,n} (s^d P)_m (P^d Q)_{mn} (Q^d s)_n * (al_m be_n)^k .
    """
    n = U[0].shape[0]
    hF = based_holonomies(U, edges, loopF)
    hC = based_holonomies(U, edges, loopC)
    if naive:
        b0F = loop_vertices(edges, loopF)[0]
        b0C = loop_vertices(edges, loopC)[0]
        hF = {v: hF[b0F] for v in hF}
        hC = {v: hC[b0C] for v in hC}
    I = np.eye(n, dtype=complex)
    zetas, coeffs = [], []
    for v in range(len(s)):
        A = hF[v].conj().T if v in hF else I
        B = hC[v] if v in hC else I
        al, P = _eig_unitary(A)
        be, Q = _eig_unitary(B)
        if P is None and Q is None:
            P = Q = I
            al = be = np.ones(n, dtype=complex)
        elif P is None:
            P = Q
            al = np.ones(n, dtype=complex)
        elif Q is None:
            Q = P
            be = np.ones(n, dtype=complex)
        sv = s[v]
        left = sv.conj() @ P            # (s^d P)_m
        right = Q.conj().T @ sv         # (Q^d s)_n
        PQ = P.conj().T @ Q
        for m in range(n):
            for nn in range(n):
                c = left[m] * PQ[m, nn] * right[nn]
                if abs(c) > 1e-14:
                    zetas.append(al[m] * be[nn])
                    coeffs.append(c)
    return np.array(zetas), np.array(coeffs)


def merge_characters(zeta, coeff, tol=1e-10):
    """merge numerically equal characters; drop zero coefficients."""
    zs, cs = [], []
    for z, c in zip(zeta, coeff):
        hit = None
        for i, z0 in enumerate(zs):
            if abs(z - z0) < tol:
                hit = i
                break
        if hit is None:
            zs.append(z)
            cs.append(c)
        else:
            cs[hit] += c
    keep = [i for i in range(len(zs)) if abs(cs[i]) > tol]
    return np.array([zs[i] for i in keep]), np.array([cs[i] for i in keep])


# ------------------------------------------------------------ direct transport -

def Z_direct(U, edges, loopF, loopC, s, k, naive=False):
    """Z_k computed by literal matrix action, no character shortcut."""
    n = U[0].shape[0]
    hF = based_holonomies(U, edges, loopF)
    hC = based_holonomies(U, edges, loopC)
    if naive:
        b0F = loop_vertices(edges, loopF)[0]
        b0C = loop_vertices(edges, loopC)[0]
        hF = {v: hF[b0F] for v in hF}
        hC = {v: hC[b0C] for v in hC}
    tot = 0.0 + 0.0j
    for v in range(len(s)):
        a = s[v].copy()
        b = s[v].copy()
        if v in hF:
            Wk = np.linalg.matrix_power(hF[v], k)
            a = Wk @ a
        if v in hC:
            Wk = np.linalg.matrix_power(hC[v], k)
            b = Wk @ b
        tot += a.conj() @ b
    return tot


def Z_from_chars(zeta, coeff, ks):
    ks = np.asarray(ks)
    return (coeff[None, :] * zeta[None, :] ** ks[:, None]).sum(axis=1)


def lambda_B(zeta, coeff, N=200000):
    ks = np.arange(1, N + 1)
    Z = Z_from_chars(zeta, coeff, ks)
    a = np.abs(Z)
    a = np.where(a < 1e-300, 1e-300, a)
    return float(np.mean(np.log(a)))


def lambda_A(zeta, coeff):
    return float(np.log(abs(Z_from_chars(zeta, coeff, [1])[0])))


# --------------------------------------------------------------- Mahler route -

def mahler_torus(exps, coeff, ngrid=2400):
    """(1/(2pi)^d) int log| sum_j c_j exp(i w_j . theta) | dtheta, d = 1 or 2,
    by uniform midpoint grid.  exps is an array (J,d) of integer exponent vectors."""
    exps = np.asarray(exps)
    d = exps.shape[1]
    if d == 1:
        th = (np.arange(ngrid) + 0.5) * 2 * np.pi / ngrid
        val = np.zeros(ngrid, dtype=complex)
        for c, w in zip(coeff, exps):
            val += c * np.exp(1j * w[0] * th)
        return float(np.mean(np.log(np.abs(val) + 1e-300)))
    th = (np.arange(ngrid) + 0.5) * 2 * np.pi / ngrid
    ph = (np.arange(ngrid) + 0.5) * 2 * np.pi / ngrid
    TH, PH = np.meshgrid(th, ph, indexing='ij')
    val = np.zeros_like(TH, dtype=complex)
    for c, w in zip(coeff, exps):
        val += c * np.exp(1j * (w[0] * TH + w[1] * PH))
    return float(np.mean(np.log(np.abs(val) + 1e-300)))


# ------------------------------------------------------------------- states ---

def state_rank1(p):
    """rank-one ready state with vertex weights p (list of 5 non-negative reals)."""
    return [np.array([np.sqrt(pi)], dtype=complex) for pi in p]


def normalise(s):
    tot = sum(float(np.vdot(x, x).real) for x in s)
    return [x / np.sqrt(tot) for x in s]


def class_weights(s, edges, loopF, loopC):
    """pushforward onto the four vertex classes (a,b) = (on F?, on C?)."""
    vF = set(loop_vertices(edges, loopF))
    vC = set(loop_vertices(edges, loopC))
    out = {(0, 0): 0.0, (1, 0): 0.0, (0, 1): 0.0, (1, 1): 0.0}
    for v, x in enumerate(s):
        a = 1 if v in vF else 0
        b = 1 if v in vC else 0
        out[(a, b)] += float(np.vdot(x, x).real)
    return out


def u1_conn(f, c, nv_edges=6):
    """rank-one U(1) connection with a1=f, a4=c and the rest zero: W_F=e^{if}, W_C=e^{ic}."""
    U = [np.array([[1.0 + 0j]]) for _ in range(nv_edges)]
    U[0] = np.array([[np.exp(1j * f)]])
    U[3] = np.array([[np.exp(1j * c)]])
    return U


def diag_conn(fs, cs, nv_edges=6):
    """abelian rank-n connection: W_F = diag(e^{i fs}), W_C = diag(e^{i cs}),
    carried entirely on e1 and e4 (all other edges identity)."""
    n = len(fs)
    U = [np.eye(n, dtype=complex) for _ in range(nv_edges)]
    U[0] = np.diag(np.exp(1j * np.asarray(fs)))
    U[3] = np.diag(np.exp(1j * np.asarray(cs)))
    return U


def su2(theta, axis):
    """exp(-i theta/2 * axis.sigma), axis a unit 3-vector."""
    ax = np.asarray(axis, dtype=float)
    ax = ax / np.linalg.norm(ax)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    n_s = ax[0] * sx + ax[1] * sy + ax[2] * sz
    return np.cos(theta / 2) * np.eye(2, dtype=complex) - 1j * np.sin(theta / 2) * n_s


def su2_conn(WF, WC, nv_edges=6):
    U = [np.eye(2, dtype=complex) for _ in range(nv_edges)]
    U[0] = WF
    U[3] = WC
    return U


def gauge_transform(U, edges, g):
    """a_e -> g_target U_e g_source^dagger."""
    return [g[t] @ U[j] @ g[s].conj().T for j, (s, t) in enumerate(edges)]


def gauge_state(s, g):
    return [g[v] @ s[v] for v in range(len(s))]


# ----------------------------------------- high-accuracy Mahler, exact routes -

def mahler1_jensen(coeff, exps=None):
    """m(P) for a one-variable Laurent polynomial with coefficients `coeff` at
    integer exponents `exps` (default 0..J-1), by Jensen:
        m = log|lead| + sum_{|root|>1} log|root|.
    Exact up to the root-finder."""
    coeff = np.asarray(coeff, dtype=complex)
    if exps is None:
        exps = np.arange(len(coeff))
    exps = np.asarray(exps, dtype=int)
    shift = exps.min()
    exps = exps - shift                      # z^shift has Mahler measure 0
    deg = exps.max()
    a = np.zeros(deg + 1, dtype=complex)
    for c, e in zip(coeff, exps):
        a[e] += c
    poly = a[::-1]                            # numpy roots wants descending
    nz = np.nonzero(np.abs(poly) > 1e-15)[0]
    poly = poly[nz[0]:]
    lead = poly[0]
    r = np.roots(poly) if len(poly) > 1 else np.array([])
    return float(np.log(abs(lead)) + np.sum(np.log(np.abs(r[np.abs(r) > 1.0]))))


def mahler_box(A_coeff, A_exps, B_coeff, B_exps, ngrid=2000001):
    """m( A(x) + y B(x) ) = (1/2pi) int log max(|A(x)|,|B(x)|) dx    (Jensen in y).
    Integrand is continuous whenever A and B have no common zero on |x|=1, so a
    plain midpoint rule converges fast -- no log singularity."""
    th = (np.arange(ngrid) + 0.5) * 2 * np.pi / ngrid
    A = np.zeros(ngrid, dtype=complex)
    B = np.zeros(ngrid, dtype=complex)
    for c, e in zip(A_coeff, A_exps):
        A += c * np.exp(1j * e * th)
    for c, e in zip(B_coeff, B_exps):
        B += c * np.exp(1j * e * th)
    m = np.maximum(np.abs(A), np.abs(B))
    return float(np.mean(np.log(m + 1e-300)))


def mahler_4class(p00, p10, p01, p11, ngrid=2000001):
    """lambda_B for the four vertex classes at rank one, unit charge.
    Characters 1, u, v, uv with exponent vectors (0,0),(-1,0),(0,1),(-1,1) in (f,c).
    P = p00 + p10 x + y(p01 + p11 x)   after clearing x^{-1} (Mahler-invariant)."""
    return mahler_box([p00, p10], [0, 1], [p01, p11], [0, 1], ngrid)
