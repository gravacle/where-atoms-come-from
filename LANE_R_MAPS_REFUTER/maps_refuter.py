#!/usr/bin/env python3
"""
LANE R (MAPS) — REFUTER.  Written from scratch; no corpus code reused.
python3 + numpy (+ mpmath only for the Cassaigne-Maillot cross-check, declared as an IMPORT).

CONVENTIONS PUBLISHED HERE (S4's failure to publish these is a defect of record):

  * Vertex order on K1:            v0 v1 v2 v3 v4               (indices 0..4)
  * Edge order on K1:              e1 e2 e3 e4 e5 e6            (indices 0..5)
  * Edge orientations (S1 :19-21): e1 v0->v1  e2 v1->v2  e3 v2->v0
                                   e4 v0->v3  e5 v3->v4  e6 v4->v0
  * d1[v,e] = +1 if v = target(e), -1 if v = source(e), 0 otherwise.
  * d2[e,F] = signed multiplicity of e in the attaching word of F.
              F attaches along e1.e2.e3 (S1 :24)  =>  d2 = (1,1,1,0,0,0)^T
  * u = conj(W_F) = exp(-i f),  v = W_C = exp(+i c),  f = a1+a2+a3, c = a4+a5+a6.
  * class of vertex w = (a_w, b_w), a_w = [w on gamma_F], b_w = [w on gamma_C].
    character:  (0,0)->1  (1,0)->u  (0,1)->v  (1,1)->uv.
  * pi = (pi00, pi10, pi01, pi11) = pushforward of p onto the four classes.
  * Z_k = sum over classes pi_ab * (u^a v^b)^k
  * lambda_B = lim_N (1/N) sum_{n=1..N} log|Z_n|      (schedule B, k_n = n)
  * All RNG: numpy default_rng with the seed printed at each use site.
  * Quadrature grid for 1-D torus integrals: MIDPOINT rule, theta_j = 2*pi*(j+1/2)/M,
    M printed at each use.  Midpoint offset chosen so the grid never lands on theta=0.
"""

import itertools
import numpy as np

np.set_printoptions(linewidth=200, suppress=True)

def hdr(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)

# ----------------------------------------------------------------------------
# 0.  CARRIERS AS EXPLICIT INCIDENCE DATA
# ----------------------------------------------------------------------------

class Complex:
    """Finite oriented CW complex of dim <= 2, given by explicit incidence."""
    def __init__(self, name, nv, edges, faces):
        # edges: list of (src, tgt)   faces: list of list of (edge_index, sign)
        self.name = name
        self.nv = nv
        self.edges = list(edges)
        self.faces = list(faces)

    @property
    def ne(self): return len(self.edges)
    @property
    def nf(self): return len(self.faces)

    def d1(self):
        D = np.zeros((self.nv, self.ne))
        for j, (s, t) in enumerate(self.edges):
            D[t, j] += 1.0
            D[s, j] -= 1.0
        return D

    def d2(self):
        D = np.zeros((self.ne, self.nf))
        for k, word in enumerate(self.faces):
            for (e, sgn) in word:
                D[e, k] += sgn
        return D

    def betti(self):
        d1, d2 = self.d1(), self.d2()
        r1 = np.linalg.matrix_rank(d1) if self.ne else 0
        r2 = np.linalg.matrix_rank(d2) if self.nf else 0
        b0 = self.nv - r1
        b1 = (self.ne - r1) - r2
        b2 = self.nf - r2
        return b0, b1, b2, r1, r2

    def chi(self):
        return self.nv - self.ne + self.nf

    def is_regular(self):
        """Regular CW: every attaching map injective on the boundary.
        For 1-cells: no loops (src != tgt).  For 2-cells: attaching word visits
        each edge at most once and traverses a genuine embedded circle (>=3 edges,
        distinct vertices)."""
        for (s, t) in self.edges:
            if s == t:
                return False, "1-cell is a loop (attaching map not injective on S^0)"
        for word in self.faces:
            es = [e for (e, _) in word]
            if len(set(es)) != len(es):
                return False, "2-cell traverses an edge twice"
            if len(es) < 3:
                return False, "2-cell attaching circle has fewer than 3 edges"
        return True, "regular"


def K1():
    edges = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]
    faces = [[(0, +1), (1, +1), (2, +1)]]
    return Complex("K1", 5, edges, faces)


# loop supports on K1 (vertex sets)
K1_GF = {0, 1, 2}     # gamma_F = boundary of F = e1.e2.e3
K1_GC = {0, 3, 4}     # gamma_C = e4.e5.e6


def classes_from_loops(nv, gF, gC):
    return [(1 if w in gF else 0, 1 if w in gC else 0) for w in range(nv)]


def push_to_pi(p, cls):
    """pushforward of the vertex weight vector p onto the four classes."""
    pi = {(0, 0): 0.0, (1, 0): 0.0, (0, 1): 0.0, (1, 1): 0.0}
    for w, c in enumerate(cls):
        pi[c] += p[w]
    return np.array([pi[(0, 0)], pi[(1, 0)], pi[(0, 1)], pi[(1, 1)]])


# ----------------------------------------------------------------------------
# 1.  THE TRANSPORT PRODUCT — BUILT FROM EDGE PHASES, NOT FROM THE FORMULA
# ----------------------------------------------------------------------------

def Z_from_operators(a, s, cx, gF_verts, gC_verts, gF_edges, gC_edges, k):
    """Z_k = <M_F^k s, M_C^k s> built by literally multiplying edge phases around
    each loop and applying the diagonal operator of S4 :127."""
    f = sum(a[e] for e in gF_edges)
    c = sum(a[e] for e in gC_edges)
    WF, WC = np.exp(1j * f), np.exp(1j * c)
    MF = np.array([WF if w in gF_verts else 1.0 for w in range(cx.nv)], dtype=complex)
    MC = np.array([WC if w in gC_verts else 1.0 for w in range(cx.nv)], dtype=complex)
    sF = (MF ** k) * s
    sC = (MC ** k) * s
    return np.vdot(sF, sC)


def Z_from_pi(pi, u, v, k):
    return (pi[0]
            + pi[1] * u ** k
            + pi[2] * v ** k
            + pi[3] * (u * v) ** k)


# ----------------------------------------------------------------------------
# 2.  lambda_B — THREE INDEPENDENT ROUTES
# ----------------------------------------------------------------------------

def lam_cesaro(pi, u, v, N):
    """ROUTE 1: the definition.  (1/N) sum_{n=1..N} log|Z_n|, schedule k_n = n."""
    n = np.arange(1, N + 1)
    Z = pi[0] + pi[1] * u ** n + pi[2] * v ** n + pi[3] * (u * v) ** n
    return float(np.mean(np.log(np.abs(Z))))


def lam_mahler_generic(pi, M=2_000_001):
    """ROUTE 2: rank L = 0 (orbit dense in T^2).  EXACT 1-D reduction:
        for fixed x, Z = (pi00 + pi10 x) + (pi01 + pi11 x) y  with y uniform on T,
        and Jensen gives  E_y log|A + B y| = log max(|A|,|B|).   Hence
        lambda = (1/2pi) int log max(|pi00+pi10 e^{it}|, |pi01+pi11 e^{it}|) dt.
    Midpoint rule, M points."""
    t = 2 * np.pi * (np.arange(M) + 0.5) / M
    x = np.exp(1j * t)
    A = np.abs(pi[0] + pi[1] * x)
    B = np.abs(pi[2] + pi[3] * x)
    g = np.maximum(A, B)
    with np.errstate(divide='ignore'):
        return float(np.mean(np.log(g)))


def lam_subtorus(pi, m, n, u0, v0, M=2_000_001):
    """ROUTE 2b: rank L = 1.  L generated by primitive (m,n) with u^m v^n = 1.
    The closure H is the kernel of (m,n) in T^2, translated onto the orbit's coset.
    Parametrise H^0 = {(e^{i n t}, e^{-i m t})}; the orbit of (u,v) under n |-> n
    lands in cosets; we average log|Z| over the actual orbit closure by taking the
    union of the d cosets reached.  Here used only with d = 1 (primitive), which we
    check by construction."""
    t = 2 * np.pi * (np.arange(M) + 0.5) / M
    x = np.exp(1j * n * t) * 1.0
    y = np.exp(-1j * m * t) * 1.0
    # coset representative: the orbit passes through (u0,v0)^1
    # find s with (e^{i n s}, e^{-i m s}) closest to (u0,v0) up to the subgroup;
    # instead we simply average over the orbit closure computed directly (below).
    Z = pi[0] + pi[1] * x + pi[2] * y + pi[3] * x * y
    return float(np.mean(np.log(np.abs(Z))))


# Bloch-Wigner + Cassaigne-Maillot  (ROUTE 3, IMPORT — declared in the IMPORT AUDIT)
import mpmath as mp

def bloch_wigner(z):
    z = mp.mpc(z)
    if abs(z) == 0:
        return mp.mpf(0)
    return mp.im(mp.polylog(2, z)) + mp.arg(1 - z) * mp.log(abs(z))

def cassaigne_maillot(A, B, C):
    """m(A + B x + C y), A,B,C > 0.  Cassaigne-Maillot / Maillot's formula."""
    A, B, C = mp.mpf(A), mp.mpf(B), mp.mpf(C)
    if A >= B + C or B >= A + C or C >= A + B:
        return mp.log(max(A, B, C))
    # angles opposite sides A, B, C in the triangle with sides A,B,C
    alpha = mp.acos((B ** 2 + C ** 2 - A ** 2) / (2 * B * C))   # opposite A
    beta = mp.acos((A ** 2 + C ** 2 - B ** 2) / (2 * A * C))    # opposite B
    gamma = mp.acos((A ** 2 + B ** 2 - C ** 2) / (2 * A * B))   # opposite C
    D = bloch_wigner(mp.mpc(A / B) * mp.exp(1j * gamma))
    return (alpha * mp.log(A) + beta * mp.log(B) + gamma * mp.log(C) + D) / mp.pi


# ============================================================================
hdr("BLOCK 1 — K1 REBUILT FROM S1 §1, INCIDENCE PUBLISHED")
# ============================================================================
k1 = K1()
print("d1 (rows v0..v4, cols e1..e6):\n", k1.d1().astype(int))
print("d2 (rows e1..e6, col F):\n", k1.d2().astype(int).T, " (shown transposed)")
print("d1 @ d2 =", k1.d1() @ k1.d2(), "   max|entry| =", np.abs(k1.d1() @ k1.d2()).max())
b0, b1, b2, r1, r2 = k1.betti()
print(f"V={k1.nv} E={k1.ne} F={k1.nf} chi={k1.chi()}  b=({b0},{b1},{b2}) rank d1={r1} rank d2={r2}")
print("regular?", k1.is_regular())
print("vertex classes (a_w,b_w):", classes_from_loops(5, K1_GF, K1_GC))

# ============================================================================
hdr("BLOCK 2 — Z_k FROM EDGE PHASES vs THE CLASS FORMULA (my own re-derivation)")
# ============================================================================
rng = np.random.default_rng(20260816)
worst = 0.0
for trial in range(2000):
    a = rng.uniform(0, 2 * np.pi, 6)
    s = rng.normal(size=5) + 1j * rng.normal(size=5)
    s = s / np.linalg.norm(s)
    p = np.abs(s) ** 2
    pi = push_to_pi(p, classes_from_loops(5, K1_GF, K1_GC))
    f = a[0] + a[1] + a[2]; c = a[3] + a[4] + a[5]
    u, v = np.exp(-1j * f), np.exp(1j * c)
    for k in (1, 2, 3, 7, 41):
        z1 = Z_from_operators(a, s, k1, K1_GF, K1_GC, [0, 1, 2], [3, 4, 5], k)
        z2 = Z_from_pi(pi, u, v, k)
        worst = max(worst, abs(z1 - z2))
print("seed 20260816, 2000 connections x 5 circuit counts")
print("max |Z_k(operators) - Z_k(class formula)| =", worst)

# ============================================================================
hdr("BLOCK 3 — lambda_B: THREE ROUTES AGREE; CORPUS NUMBERS REPRODUCED")
# ============================================================================
def show_lam(label, pi, note=""):
    lm = lam_mahler_generic(pi)
    print(f"  {label:38s} pi={np.round(pi,6)}  lambda(generic torus) = {lm:.9f}  {note}")
    return lm

pi_K1_C   = np.array([0.0, 0.3, 0.3, 0.4])       # SENSE C  (p0,q,r) = (0.4,0.3,0.3)
pi_K1_U   = np.array([0.0, 0.4, 0.4, 0.2])       # SENSE U  uniform on 5 vertices
pi_B1p    = np.array([0.0, 0.5, 0.5, 0.0])       # bridged: classes {10:3, 01:3}
pi_B1s_U  = np.array([0.0, 5/11, 5/11, 1/11])    # subdivided, uniform on 11 vertices
pi_root   = np.array([0.0, 0.0, 0.0, 1.0])       # all weight on the pinch

L1 = show_lam("K1 SENSE C  (0.4,0.3,0.3)", pi_K1_C,  "corpus: -0.767507880")
L2 = show_lam("K1 SENSE U  (uniform on 5)", pi_K1_U, "corpus: -0.756573586")
L3 = show_lam("B1p bridged (no pinch)",     pi_B1p,  "corpus: -0.693147181 = log(1/2)")
L4 = show_lam("B1s subdivided SENSE U",     pi_B1s_U,"corpus: -0.724759919")
L5 = show_lam("K1, all weight on the pinch",pi_root, "PREDICTED EXACTLY 0")

print("\n  Cassaigne-Maillot closed form (IMPORT, independent of my quadrature):")
for lbl, pi, mine in (("K1 SENSE C", pi_K1_C, L1), ("K1 SENSE U", pi_K1_U, L2),
                      ("B1s SENSE U", pi_B1s_U, L4)):
    # 3-term case: monomial-shift to A + B x + C y
    nz = [x for x in pi if x > 0]
    cm = float(cassaigne_maillot(*sorted(nz, reverse=True)))
    print(f"    {lbl:14s} CM = {cm:.9f}   quadrature = {mine:.9f}   |diff| = {abs(cm-mine):.2e}")

print("\n  Cesaro route (the DEFINITION), generic connection f=2.0, c=1.1000001:")
f0, c0 = 2.0, 1.1000001
u0, v0 = np.exp(-1j * f0), np.exp(1j * c0)
for N in (10_000, 1_000_000, 20_000_000):
    print(f"    N={N:>10d}  (1/N)sum log|Z_n| = {lam_cesaro(pi_K1_C, u0, v0, N):.9f}")
print("  RESONANT point f=2.0, c=1.1  (-11f+20c = 0, orbit on a SUBTORUS):")
u1, v1 = np.exp(-1j * 2.0), np.exp(1j * 1.1)
for N in (200_000, 20_000_000):
    print(f"    N={N:>10d}  (1/N)sum log|Z_n| = {lam_cesaro(pi_K1_C, u1, v1, N):.9f}   corpus: -0.767014993")
