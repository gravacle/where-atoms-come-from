#!/usr/bin/env python3
"""Exact/finite checks for the SDCP proof draft (standard library only)."""

from fractions import Fraction as Q
from itertools import product, combinations
from math import exp


def det_fraction(a):
    a = [row[:] for row in a]
    n = len(a)
    out = Q(1)
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col]), None)
        assert pivot is not None
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            out = -out
        p = a[col][col]
        out *= p
        for j in range(col, n):
            a[col][j] /= p
        for r in range(col + 1, n):
            c = a[r][col]
            if c:
                for j in range(col, n):
                    a[r][j] -= c * a[col][j]
    return out


def outer(u, v):
    return [[u[i] * v[j] for j in range(3)] for i in range(3)]


def madd(a, b):
    return [[a[i][j] + b[i][j] for j in range(3)] for i in range(3)]


def flatten_sym(m):
    return [m[0][0], m[1][1], m[2][2], m[0][1], m[0][2], m[1][2]]


# An orthonormal realization of v_a=P e_a in V=1^perp.
v = [
    [Q(1, 2), Q(1, 2), Q(1, 2)],
    [Q(1, 2), Q(-1, 2), Q(-1, 2)],
    [Q(-1, 2), Q(1, 2), Q(-1, 2)],
    [Q(-1, 2), Q(-1, 2), Q(1, 2)],
]

# Tight tetrahedral frame: sum_a v_a v_a^T = I_3.
f0 = [[Q(0) for _ in range(3)] for _ in range(3)]
for va in v:
    f0 = madd(f0, outer(va, va))
assert f0 == [[Q(int(i == j)) for j in range(3)] for i in range(3)]

# The six pair derivatives v_a v_b^T+v_b v_a^T span Sym^2(V).
edges = list(combinations(range(4), 2))
cols = []
for a, b in edges:
    cols.append(flatten_sym(madd(outer(v[a], v[b]), outer(v[b], v[a]))))
matrix = [[cols[c][r] for c in range(6)] for r in range(6)]
det6 = det_fraction(matrix)
assert det6 != 0

# Exact uniform-family expectations over the 16 four-bit outcomes.
states = list(product((-1, 1), repeat=4))
assert len(states) == 16
for a in range(4):
    assert sum(Q(s[a], 16) for s in states) == 0

# Pair sufficient statistics have identity covariance at J=0.
phis = [[s[a] * s[b] for a, b in edges] for s in states]
for e in range(6):
    for f in range(6):
        val = sum(Q(row[e] * row[f], 16) for row in phis)
        assert val == Q(int(e == f))

# F_theta(0)=I and dF/dJ_ab equals the symmetric tetrahedral dyad.
def x_of(s):
    return [sum(v[a][i] * s[a] for a in range(4)) for i in range(3)]

fisher = [[Q(0) for _ in range(3)] for _ in range(3)]
for s in states:
    x = x_of(s)
    for i in range(3):
        for j in range(3):
            fisher[i][j] += x[i] * x[j] / 16
assert fisher == f0

for e, (a, b) in enumerate(edges):
    deriv = [[Q(0) for _ in range(3)] for _ in range(3)]
    for s in states:
        x = x_of(s)
        phi = s[a] * s[b]
        for i in range(3):
            for j in range(3):
                deriv[i][j] += x[i] * x[j] * phi / 16
    target = madd(outer(v[a], v[b]), outer(v[b], v[a]))
    assert deriv == target

# Global-flip symmetry keeps every one-port marginal exactly one half.
jvals = [0.13, -0.21, 0.08, 0.17, -0.04, 0.11]
weights = []
for k, s in enumerate(states):
    weights.append(exp(sum(jvals[e] * phis[k][e] for e in range(6))))
z = sum(weights)
for a in range(4):
    p_plus = sum(w for w, s in zip(weights, states) if s[a] == 1) / z
    assert abs(p_plus - 0.5) < 1e-14

print("PASS SDCP exact checks")
print(f"tetra_pair_deformation_det={det6}")
print("pair_covariance_rank=6; spatial_fisher_tangent_rank=6")
print("global_flip_single_port_marginals=1/2")
