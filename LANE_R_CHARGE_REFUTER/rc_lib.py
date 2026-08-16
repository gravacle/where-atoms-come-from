"""
LANE_R_CHARGE_REFUTER -- independent library, written from scratch for this lane.
No code, constant or convention is copied from LANE_S5_CHARGE_CODE or from any S4 lane.

PUBLISHED CONVENTIONS (this lane's own; stated here so every number below is re-derivable)

  Carrier K1 (S1_CARRIER_K1_V001.md sections 1,3):
    vertices v0..v4 -> indices 0..4
    edges e1..e6    -> (0,1) (1,2) (2,0) (0,3) (3,4) (4,0)   [source, target]
    face  F         -> attached along e1+e2+e3
    d1[v,e] = +1 if v = target(e), -1 if v = source(e), 0 else      (5 x 6)
    d2[e,F] = coefficient of e in the attaching 1-chain of F        (6 x 1)

  Connection:  U_e = exp(i a_e);  W_F = exp(i(a1+a2+a3));  W_C = exp(i(a4+a5+a6)).
  S2/S3/S4 characters:  u = conj(W_F) = exp(-i f),  v = W_C = exp(i c),  f,c real mod 2pi.
  I write the character exponent of a vertex as E_v = (m_v, n_v) in Z^2 and
      chi_{E}(f,c) = u^m v^n = exp( i ( -m f + n c ) ).
  Unit charge:  E_v = (a_v, b_v) with a_v = [v on gamma_F], b_v = [v on gamma_C].
  Per-vertex charge q_v in Z:  E_v = q_v (a_v, b_v).      (the ONLY reading that is
      gauge-invariant without an edge transport; see rc_1.py section 2)

  Transport product:  Z_k = sum_v p_v chi_{E_v}(f,c)^k = sum_v p_v exp(i k (-m_v f + n_v c)).
  Schedule k_n; Omega_N = prod_{n=1..N} Z_{k_n};  lambda = lim (1/N) log|Omega_N|.
      Schedule A: k_n = 1.  Schedule B (canonical clock): k_n = n.
      Schedule B[M]: k_n = M n  (M a positive integer) -- an admissible clock, used below.

  Relation lattice   L     = { (m,n) in Z^2 : chi_{(m,n)} = 1 }.
  Difference lattice Delta = < E_x - E_y : x,y in S >,  S = supp(p).
  Formation group    G     = < chi_x / chi_y : x,y in S >  <= U(1).

  Rational connections: f = 2 pi A / M, c = 2 pi Bq / M with integers A,Bq,M, so that
      L = { (m,n) : (-m A + n Bq) = 0 mod M }  is EXACT integer arithmetic.

  Mahler measure route for lambda_B at a generic connection (L = 0):
      reduce E to a basis of Delta -> integer coordinates c_v;
      the map theta |-> (<b1,theta>,<b2,theta>) : T^2 -> T^{rank} is a surjective
      homomorphism of compact groups, so it pushes Haar to Haar, hence
      lambda_B = m( sum_v p_v X^{c1_v} Y^{c2_v} )  exactly (rank 2),
      lambda_B = m( sum_v p_v X^{c_v} )            exactly (rank 1),
      lambda_B = 0                                          (rank 0).
      1-variable m: EXACT by Jensen (numpy.roots).
      2-variable m: exact-in-Y Jensen + midpoint quadrature in X, Nx published at call site.

  Floating point: float64 everywhere.  Seeds: numpy.random.default_rng, stated at call site.
"""

import numpy as np

# ---------------------------------------------------------------- carrier K1

K1_EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]
K1_FACES = [{0: 1, 1: 1, 2: 1}]          # F = e1+e2+e3
K1_GAMMA_F = [0, 1, 2]                    # edge indices of the face loop
K1_GAMMA_C = [3, 4, 5]                    # edge indices of the free cycle


def d1_matrix(nv, edges):
    d = np.zeros((nv, len(edges)), dtype=np.int64)
    for j, (s, t) in enumerate(edges):
        d[s, j] -= 1
        d[t, j] += 1
    return d


def d2_matrix(ne, faces):
    d = np.zeros((ne, len(faces)), dtype=np.int64)
    for j, ch in enumerate(faces):
        for e, coef in ch.items():
            d[e, j] += coef
    return d


def vertices_on(edges, loop):
    s = set()
    for e in loop:
        s.add(edges[e][0])
        s.add(edges[e][1])
    return s


def class_vectors(nv, edges, gF, gC):
    """unit-charge exponent vectors E_v = (a_v, b_v)."""
    A = vertices_on(edges, gF)
    B = vertices_on(edges, gC)
    return np.array([[1 if v in A else 0, 1 if v in B else 0] for v in range(nv)],
                    dtype=np.int64)


# ------------------------------------------------------- the transport product

def Zk(E, p, f, c, k):
    """Z_k = sum_v p_v exp(i k (-m_v f + n_v c)).  E is (V,2) integer, p is (V,)."""
    ph = k * (-E[:, 0] * f + E[:, 1] * c)
    return np.sum(p * np.exp(1j * ph))


def Zk_matrix(E, s, f, c, k):
    """Same object built as W-01's operator inner product, no closed form used:
       M_F = diag(W_F^{m_v}), M_C = diag(W_C^{n_v}) acting on Gamma(L) = C^V,
       Z_k = < M_F^k s, M_C^k s > = sum_v conj(W_F)^{k m_v} W_C^{k n_v} p_v.
       (E_v = (m_v,n_v); at unit charge m_v = a_v, n_v = b_v.)"""
    V = len(s)
    WF = np.exp(1j * f)
    WC = np.exp(1j * c)
    MF = np.diag(np.array([WF ** int(E[j, 0]) for j in range(V)]))
    MC = np.diag(np.array([WC ** int(E[j, 1]) for j in range(V)]))
    a = np.linalg.matrix_power(MF, k) @ s
    b = np.linalg.matrix_power(MC, k) @ s
    return np.vdot(a, b)


# --------------------------------------------------------------- 2D lattices

def _hnf2(rows):
    """Hermite-style reduction of a set of integer row vectors in Z^2.
       Returns a list of 0,1 or 2 basis rows (echelon, non-negative pivots)."""
    R = [list(map(int, r)) for r in rows if any(r)]
    if not R:
        return []
    # eliminate on first coordinate by integer gcd (Euclid on rows)
    while True:
        nz = [r for r in R if r[0] != 0]
        if len(nz) <= 1:
            break
        nz.sort(key=lambda r: abs(r[0]))
        piv = nz[0]
        newR = []
        for r in R:
            if r is piv:
                newR.append(r)
            elif r[0] != 0:
                q = r[0] // piv[0]
                newR.append([r[0] - q * piv[0], r[1] - q * piv[1]])
            else:
                newR.append(r)
        if all(x == y for x, y in zip(sorted(map(tuple, newR)), sorted(map(tuple, R)))):
            break
        R = [r for r in newR if any(r)]
    first = [r for r in R if r[0] != 0]
    rest = [r for r in R if r[0] == 0]
    basis = []
    if first:
        b = first[0]
        if b[0] < 0:
            b = [-b[0], -b[1]]
        basis.append(b)
    g = 0
    for r in rest:
        g = np.gcd(g, abs(r[1]))
    if g:
        basis.append([0, int(g)])
    return basis


def lattice_basis(vectors):
    return _hnf2(vectors)


def delta_lattice(E, p, tol=0.0):
    """Delta = < E_x - E_y : x,y in supp(p) >."""
    idx = [i for i in range(len(p)) if p[i] > tol]
    diffs = []
    for i in idx[1:]:
        diffs.append(E[i] - E[idx[0]])
    return lattice_basis(diffs)


def in_L(mn, A, Bq, M):
    """membership in L = {(m,n): -mA + nBq = 0 mod M} for f=2piA/M, c=2piBq/M."""
    return (-mn[0] * A + mn[1] * Bq) % M == 0


def delta_subset_L(basis, A, Bq, M):
    return all(in_L(b, A, Bq, M) for b in basis)


# ------------------------------------------------------------ Mahler measures

def mahler1(coeffs):
    """m(sum_j coeffs[j] x^j), exact by Jensen. coeffs a real/complex 1-D array,
       leading entry may be zero (trailing monomials are stripped: m is unchanged)."""
    c = np.array(coeffs, dtype=complex)
    nz = np.nonzero(np.abs(c) > 0)[0]
    c = c[nz[0]:nz[-1] + 1]
    if len(c) == 1:
        return float(np.log(abs(c[0])))
    r = np.roots(c[::-1])                    # numpy wants highest degree first
    lead = abs(c[-1])
    return float(np.log(lead) + np.sum(np.log(np.abs(r[np.abs(r) > 1.0]))))


def mahler2(terms, Nx):
    """m( sum_v p_v X^{i_v} Y^{j_v} ) by midpoint quadrature in X (Nx nodes,
       t_i = (i+0.5)/Nx, X = exp(2 pi i t_i)) and exact Jensen in Y."""
    js = [j for (_, _, j) in terms]
    jmin, jmax = min(js), max(js)
    deg = jmax - jmin
    tot = 0.0
    for i in range(Nx):
        t = (i + 0.5) / Nx
        X = np.exp(2j * np.pi * t)
        coef = np.zeros(deg + 1, dtype=complex)
        for (pv, iv, jv) in terms:
            coef[jv - jmin] += pv * X ** iv
        tot += mahler1(coef)
    return tot / Nx


def lambda_B_generic(E, p, Nx=16384):
    """lambda_B at a connection with L = 0 (generic), exact route via the Delta basis."""
    idx = [i for i in range(len(p)) if p[i] > 0]
    E0 = E[idx[0]]
    D = [E[i] - E0 for i in idx]
    basis = lattice_basis(D)
    r = len(basis)
    # aggregate equal exponent vectors
    agg = {}
    for i in idx:
        agg[tuple(E[i] - E0)] = agg.get(tuple(E[i] - E0), 0.0) + p[i]
    if r == 0:
        return 0.0
    if r == 1:
        b = np.array(basis[0], dtype=np.int64)
        # coordinate of d in the basis
        coords = {}
        for d, w in agg.items():
            dv = np.array(d, dtype=np.int64)
            k = dv[0] // b[0] if b[0] != 0 else dv[1] // b[1]
            assert np.all(dv == k * b), (dv, b)
            coords[int(k)] = coords.get(int(k), 0.0) + w
        kmin, kmax = min(coords), max(coords)
        coef = np.zeros(kmax - kmin + 1)
        for k, w in coords.items():
            coef[k - kmin] += w
        return mahler1(coef)
    Bm = np.array(basis, dtype=np.int64).T           # columns b1,b2 ; solve Bm @ x = d
    det = int(round(np.linalg.det(Bm)))
    terms = []
    for d, w in agg.items():
        x = np.linalg.solve(Bm.astype(float), np.array(d, dtype=float))
        xi = np.rint(x).astype(int)
        assert np.allclose(x, xi, atol=1e-9), (d, x)
        terms.append((w, int(xi[0]), int(xi[1])))
    return mahler2(terms, Nx)


def lambda_direct(E, p, f, c, N, schedule=lambda n: n):
    """(1/N) sum_{n=1..N} log|Z_{k_n}|, no closed form used."""
    tot = 0.0
    for n in range(1, N + 1):
        z = abs(Zk(E, p, f, c, schedule(n)))
        tot += np.log(z) if z > 0 else -np.inf
    return tot / N
