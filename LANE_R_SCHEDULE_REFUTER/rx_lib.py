"""
LANE R (schedule refuter) -- independent re-derivation, written from scratch.
No code reused from any other lane. python3 + numpy only.

CARRIER K1 -- PUBLISHED CONVENTIONS (this file is the publication)
-----------------------------------------------------------------
Vertices  v0..v4                                   (5 vertices)
Edges     e1 = v0->v1 , e2 = v1->v2 , e3 = v2->v0   (filled triangle, face F)
          e4 = v0->v3 , e5 = v3->v4 , e6 = v4->v0   (unfilled triangle, cycle C)
The two triangles meet only at v0 (the pinch / root).

d1 (edge->vertex boundary, columns = edges e1..e6, rows = vertices v0..v4):
        e1  e2  e3  e4  e5  e6
   v0 [ -1   0  +1  -1   0  +1 ]
   v1 [ +1  -1   0   0   0   0 ]
   v2 [  0  +1  -1   0   0   0 ]
   v3 [  0   0   0  +1  -1   0 ]
   v4 [  0   0   0   0  +1  -1 ]

d2 (face->edge boundary, single column = face F, rows = edges e1..e6):
   [ +1, +1, +1, 0, 0, 0 ]^T

d2 IS NEVER USED BY ANY NUMBER IN THIS FILE.  That is W-03's ruling ("topology
is inert" is an analytic fact of the transport, not a finding) and this lane
reproduces it by construction: no function below takes d2 as an argument.

CONNECTION AND TRANSPORT
------------------------
U(1) connection a = (a1..a6) on edges.  W_F = exp(i(a1+a2+a3)),  W_C = exp(i(a4+a5+a6)).
Following S3 sec 4.1 we parametrise by f = arg W_F, c = arg W_C, and set
   u = conj(W_F) = e^{-i f},   v = W_C = e^{+i c}.

Vertex class (a,b) = (touched by loop F, touched by loop C):
   v0 -> (1,1)   [pinch]      v1,v2 -> (1,0)      v3,v4 -> (0,1)
   class (0,0) = spectator, empty on K1.

Branch operators on C^5 are DIAGONAL (S3 sec 2: every vertex phase of the ready
state cancels at U(1) -- true at U(1), FALSE at SU(2) per W-03; we stay at U(1)):
   M_dF = diag( e^{+i f}, e^{+i f}, e^{+i f}, 1, 1 )      (charge 1)
   M_c  = diag( e^{+i c}, 1, 1, e^{+i c}, e^{+i c} )

Ready state s = sum_j sqrt(p_j) |j>,  p a probability vector on the 5 vertices.

One-cell overlap   Z_k = < M_dF^k s , M_c^k s >
                       = e^{ik(c-f)} p0 + e^{-ikf}(p1+p2) + e^{+ikc}(p3+p4)
(both forms implemented and cross-checked).

RECORD FUNCTIONAL / RATE
------------------------
Schedule (k_n).   Omega_N = prod_{n=1..N} Z_{k_n}.
lambda(schedule) = lim_N (1/N) sum_{n=1..N} log|Z_{k_n}|.
Schedule A: k_n = 1  -> lambda_A = log|Z_1|.
Schedule B: k_n = n  -> lambda_B = average of log|Z| over the closure of the
                        cyclic group <(u,v)> in T^2 (Weyl).
Alpha family: k_n = floor(n^alpha).   alpha=0 -> k_n=1 (schedule A).
"""

import numpy as np

# ---------------------------------------------------------------- incidence
D1 = np.array([
    [-1,  0, +1, -1,  0, +1],
    [+1, -1,  0,  0,  0,  0],
    [ 0, +1, -1,  0,  0,  0],
    [ 0,  0,  0, +1, -1,  0],
    [ 0,  0,  0,  0, +1, -1],
], dtype=int)

D2 = np.array([[+1], [+1], [+1], [0], [0], [0]], dtype=int)

LOOP_F_EDGES = [0, 1, 2]   # e1,e2,e3
LOOP_C_EDGES = [3, 4, 5]   # e4,e5,e6

# vertex -> (in F, in C)
VCLASS = [(1, 1), (1, 0), (1, 0), (0, 1), (0, 1)]


def check_complex():
    """d1 @ d2 == 0 and each declared loop is a cycle (d1 @ indicator == 0)."""
    out = {}
    out["d1_d2"] = int(np.abs(D1 @ D2).max())
    for name, ed in (("F", LOOP_F_EDGES), ("C", LOOP_C_EDGES)):
        x = np.zeros((6, 1), dtype=int)
        x[ed, 0] = 1
        out["cycle_" + name] = int(np.abs(D1 @ x).max())
    # the two triangles share exactly vertex v0
    sf = set(np.nonzero(np.abs(D1[:, LOOP_F_EDGES]).sum(axis=1))[0])
    sc = set(np.nonzero(np.abs(D1[:, LOOP_C_EDGES]).sum(axis=1))[0])
    out["shared_vertices"] = sorted(sf & sc)
    return out


def holonomies(a):
    """a = 6-vector of edge angles -> (W_F, W_C)."""
    a = np.asarray(a, dtype=float)
    return (np.exp(1j * a[LOOP_F_EDGES].sum()),
            np.exp(1j * a[LOOP_C_EDGES].sum()))


# ---------------------------------------------------------------- transport
def Z_matrix(k, f, c, p):
    """Z_k by literal diagonal matrix action on C^5 (no closed form)."""
    p = np.asarray(p, dtype=float)
    s = np.sqrt(p).astype(complex)
    dF = np.array([np.exp(1j * f)] * 3 + [1.0, 1.0], dtype=complex)
    dC = np.array([np.exp(1j * c), 1.0, 1.0, np.exp(1j * c), np.exp(1j * c)],
                  dtype=complex)
    bF = (dF ** k) * s
    bC = (dC ** k) * s
    return np.vdot(bF, bC)          # <bF, bC> = sum conj(bF)*bC


def Z_closed(k, f, c, p):
    """Z_k closed form (S3 sec 4.1)."""
    p = np.asarray(p, dtype=float)
    p0 = p[0]
    q = p[1] + p[2]
    r = p[3] + p[4]
    k = np.asarray(k)
    return (np.exp(1j * k * (c - f)) * p0
            + np.exp(-1j * k * f) * q
            + np.exp(1j * k * c) * r)


def logabsZ(k, f, c, p):
    z = np.abs(Z_closed(k, f, c, p))
    with np.errstate(divide="ignore"):
        return np.log(z)


# ---------------------------------------------------------------- schedules
def schedule_alpha(N, alpha):
    """k_n = floor(n^alpha), n = 1..N.  alpha = 0 gives k_n = 1 (schedule A)."""
    n = np.arange(1, N + 1, dtype=float)
    if alpha == 0.0:
        return np.ones(N, dtype=np.int64)
    return np.floor(n ** alpha).astype(np.int64)


def lam_of_schedule(ks, f, c, p):
    """(1/N) sum log|Z_{k_n}|  -- exact finite-stage value."""
    ks = np.asarray(ks, dtype=np.int64)
    return float(np.mean(logabsZ(ks.astype(float), f, c, p)))


def lam_alpha(N, alpha, f, c, p):
    return lam_of_schedule(schedule_alpha(N, alpha), f, c, p)


# ---------------------------------------------------------------- orbit closure
def orbit_order(f, c, maxm=100000, tol=1e-12):
    """Smallest m>=1 with (u^m, v^m) = (1,1) within tol; None if none <= maxm."""
    for m in range(1, maxm + 1):
        if abs(np.exp(-1j * f * m) - 1) < tol and abs(np.exp(1j * c * m) - 1) < tol:
            return m
    return None


def lambda_B_finite_orbit(m, f, c, p):
    """Exact lambda_B when the closure is the finite cyclic group of order m."""
    ks = np.arange(m, dtype=float)
    return float(np.mean(logabsZ(ks, f, c, p)))
