#!/usr/bin/env python3
"""
REFUTER LANE (charge axis) -- part 1.
Independent re-derivation of K1, the charged transport, and the exponent data.
Written from scratch; no code reused from any other lane.

PUBLISHED CONVENTIONS (this file is the publication).
  vertices v0..v4  ->  indices 0..4
  edges e1..e6     ->  (0,1) (1,2) (2,0) (0,3) (3,4) (4,0)      [S1 sec.1, source->target]
  face  F          ->  attached along e1+e2+e3                   [S1 sec.1]
  d1[v,e] = +1 if v = target(e), -1 if v = source(e), else 0     (5x6)
  d2[e,F] = coefficient of e in boundary of F                    (6x1)
  gamma_F = e1+e2+e3 (bounds F)      gamma_C = e4+e5+e6 (free cycle)
  a_v = [v on gamma_F]   b_v = [v on gamma_C]     (SET membership -- S4 CHOICE C11)
  u = conj(W_F) = e^{-if},  v = W_C = e^{ic}
  per-vertex charge q_v in Z; U(1) acts on the fibre L_v by z -> z^{q_v}
  exponent point E_v = q_v * (a_v, b_v) in Z^2
  Z_k = <M_F^k s, M_C^k s> = sum_v p_v exp( i k ( -q_v a_v f + q_v b_v c ) )
  Delta(S) = < E_x - E_y : x,y in S >   (difference lattice, a subgroup of Z^2)
  L        = { (m,n) in Z^2 : -m f + n c = 0 mod 2pi }   (relation lattice)
  G        = < chi_x/chi_y > <= U(1)  ;  G = image of Delta under (m,n) -> e^{i(-mf+nc)}
SEEDS: numpy.random.default_rng(70251101) -- section 3 only. Nothing else is random.
"""
import numpy as np
from fractions import Fraction
from math import gcd, log, pi, cos

np.set_printoptions(linewidth=150)

# ---------------------------------------------------------------- 1. the carrier
EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]
V, E = 5, 6
d1 = np.zeros((V, E), dtype=np.int64)
for j, (s, t) in enumerate(EDGES):
    d1[s, j] -= 1
    d1[t, j] += 1
d2 = np.zeros((E, 1), dtype=np.int64)
for j in (0, 1, 2):
    d2[j, 0] = 1

print("=" * 96)
print("R1.  K1 REBUILT FROM S1 SECTION 1.  d1, d2 PUBLISHED.")
print("=" * 96)
print("d1 (rows v0..v4, cols e1..e6) =")
print(d1)
print("d2^T (cols e1..e6) =", d2.T)
print("d1 @ d2 =", (d1 @ d2).T, "   max|entry| =", int(np.abs(d1 @ d2).max()))
b0 = V - np.linalg.matrix_rank(d1)
b2 = d2.shape[1] - np.linalg.matrix_rank(d2)
b1 = E - np.linalg.matrix_rank(d1) - np.linalg.matrix_rank(d2)
print(f"V={V} E={E} F=1  chi={V-E+1}  b0={b0} b1={b1} b2={b2}   (ranks over Q)")

gamma_F = np.array([1, 1, 1, 0, 0, 0], dtype=np.int64)
gamma_C = np.array([0, 0, 0, 1, 1, 1], dtype=np.int64)
print("gamma_F is a cycle:", not (d1 @ gamma_F).any(), " bounds:", True, "(= d2 column)")
print("gamma_C is a cycle:", not (d1 @ gamma_C).any(),
      " bounds:", bool(np.linalg.matrix_rank(np.hstack([d2, gamma_C.reshape(-1, 1)])) == 1))

a_v = np.zeros(V, dtype=np.int64)
b_v = np.zeros(V, dtype=np.int64)
for j in (0, 1, 2):
    a_v[EDGES[j][0]] = 1; a_v[EDGES[j][1]] = 1
for j in (3, 4, 5):
    b_v[EDGES[j][0]] = 1; b_v[EDGES[j][1]] = 1
print("a_v =", a_v, "  b_v =", b_v)
print("classes: v0=(1,1)  v1,v2=(1,0)  v3,v4=(0,1)")

# ------------------------------------------------- 2. charged transport, matrix vs closed form
def Z_matrix(f, c, q, p_amp, k):
    """Z_k from the OPERATOR definition: diagonal multiplication by W(gamma)^{q_v} at
    vertices on gamma.  s has amplitudes p_amp (complex, arbitrary phases)."""
    WF, WC = np.exp(1j * f), np.exp(1j * c)
    MF = np.diag([WF ** (q[v] * a_v[v]) for v in range(V)])
    MC = np.diag([WC ** (q[v] * b_v[v]) for v in range(V)])
    sF = np.linalg.matrix_power(MF, k) @ p_amp
    sC = np.linalg.matrix_power(MC, k) @ p_amp
    return np.vdot(sF, sC)

def Z_closed(f, c, q, p, k):
    """Z_k from the exponent points E_v = q_v (a_v,b_v)."""
    tot = 0j
    for v in range(V):
        m, n = q[v] * a_v[v], q[v] * b_v[v]
        tot += p[v] * np.exp(1j * k * (-m * f + n * c))
    return tot

rng = np.random.default_rng(70251101)
worst = 0.0
for _ in range(400):
    f, c = rng.uniform(0, 2 * pi, 2)
    q = rng.integers(-3, 4, V)
    amp = rng.uniform(0.1, 1.0, V)
    ph = rng.uniform(0, 2 * pi, V)
    s = amp * np.exp(1j * ph)
    p = amp ** 2
    k = int(rng.integers(1, 6))
    worst = max(worst, abs(Z_matrix(f, c, q, s, k) - Z_closed(f, c, q, p, k)))
print(f"\nmatrix operator vs closed form, 400 random (f,c,q,section,k), seed 70251101:"
      f"  max |diff| = {worst:.3e}")
print("  (vertex phases of s cancel: confirmed -- only p_v = |s_v|^2 enters)")

# ---------------------------------------------------------------- 3. exact lattice tools
def lattice_rank(vecs):
    """rank over Q of a set of integer 2-vectors (exact integer arithmetic)."""
    vs = [v for v in vecs if v != (0, 0)]
    if not vs:
        return 0
    for i in range(len(vs)):
        for j in range(i + 1, len(vs)):
            if vs[i][0] * vs[j][1] - vs[i][1] * vs[j][0] != 0:
                return 2
    return 1

def delta_rank(points):
    """rank of the DIFFERENCE lattice = affine dimension over Q of the point set."""
    pts = sorted(set(points))
    if len(pts) <= 1:
        return 0
    base = pts[0]
    return lattice_rank([(p[0] - base[0], p[1] - base[1]) for p in pts])

def E_points(q):
    return [(q[v] * a_v[v], q[v] * b_v[v]) for v in range(V)]

print("\n" + "=" * 96)
print("R2.  S4-1's OWN CONFIGURATION, REPRODUCED.  (unit charge, all 15 support subsets)")
print("=" * 96)
CORNERS = {(0, 0): (0, 0), (1, 0): (1, 0), (0, 1): (0, 1), (1, 1): (1, 1)}
keys = list(CORNERS)
tally = {0: 0, 1: 0, 2: 0}
from itertools import combinations
for r in range(1, 5):
    for sub in combinations(keys, r):
        tally[delta_rank([CORNERS[k] for k in sub])] += 1
print(f"  rank2={tally[2]}  rank1={tally[1]}  rank0={tally[0]}   (S4-1 predicts 5 / 6 / 4)")
print("  reproduced:", (tally[2], tally[1], tally[0]) == (5, 6, 4))
