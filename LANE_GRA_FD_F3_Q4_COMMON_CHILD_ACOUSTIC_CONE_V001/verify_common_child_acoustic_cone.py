#!/usr/bin/env python3
"""Independent finite/algebraic checks for CCMAC V001."""

from itertools import product
from math import comb, sqrt
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
THEOREM = (ROOT / "THEOREM.md").read_text()
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def front(depth):
    return [m for m in product(range(depth + 1), repeat=4)
            if sum(m) == depth]


# FD01--FD04: finite incidence identity and sibling-root criterion.
for depth in (0, 1, 2, 4):
    parents = front(depth)
    children = front(depth + 1)
    child_index = {m: i for i, m in enumerate(children)}
    B = np.zeros((len(children), len(parents)), dtype=int)
    for j, m in enumerate(parents):
        for a in range(4):
            c = list(m)
            c[a] += 1
            B[child_index[tuple(c)], j] = 1
    sibling = np.zeros((len(parents), len(parents)), dtype=int)
    for i, left in enumerate(parents):
        for j, right in enumerate(parents):
            if i == j:
                continue
            diff = np.array(right) - np.array(left)
            sibling[i, j] = int(
                np.count_nonzero(diff == 1) == 1
                and np.count_nonzero(diff == -1) == 1
                and np.count_nonzero(diff) == 2
            )
    check(np.array_equal(B.T @ B, 4 * np.eye(len(parents), dtype=int) + sibling),
          f"finite BdagB identity N={depth}")
    check(len(parents) == comb(depth + 3, 3),
          f"stars and bars N={depth}")


# FD09--FD11: tetrahedral/root algebra.
one = np.ones(4)
P = np.eye(4) - np.outer(one, one) / 4
n = (2 / sqrt(3)) * P
gram = n.T @ n
check(np.allclose(np.diag(gram), 1), "tetrahedral unit norms")
check(np.allclose(gram - np.diag(np.diag(gram)),
                  (np.ones((4, 4)) - np.eye(4)) * (-1 / 3)),
      "tetrahedral cross Gram")
check(np.allclose(n.sum(axis=1), 0), "tetrahedral zero sum")

roots = []
for a in range(4):
    for b in range(a + 1, 4):
        roots.append(n[:, b] - n[:, a])
root_moment = sum(np.outer(alpha, alpha) for alpha in roots)
check(np.allclose(root_moment, (16 / 3) * P),
      "A3 root second moment on V")

# Six root dyads span Sym^2(V).
basis_v = np.linalg.qr(P[:, :3])[0][:, :3]
dyads = []
for alpha in roots:
    q = basis_v.T @ np.outer(alpha, alpha) @ basis_v
    dyads.append([q[0, 0], q[1, 1], q[2, 2],
                  sqrt(2) * q[0, 1], sqrt(2) * q[0, 2],
                  sqrt(2) * q[1, 2]])
check(np.linalg.matrix_rank(np.array(dyads)) == 6,
      "six root dyads span symmetric tangent")


# FD21--FD24: primitive covolume and tetrahedron volume.
primitive = np.column_stack([
    n[:, 3] - n[:, 0],
    n[:, 3] - n[:, 1],
    n[:, 3] - n[:, 2],
])
primitive_v = basis_v.T @ primitive
covolume = abs(np.linalg.det(primitive_v))
check(np.isclose(covolume, 16 / (3 * sqrt(3))),
      "A3 relational cell covolume")

tet_edges = basis_v.T @ np.column_stack([
    n[:, 1] - n[:, 0],
    n[:, 2] - n[:, 0],
    n[:, 3] - n[:, 0],
])
tet_volume = abs(np.linalg.det(tet_edges)) / 6
check(np.isclose(tet_volume, 8 / (9 * sqrt(3))),
      "unit relational tetrahedron volume")
check(np.isclose(tet_volume / covolume, 1 / 6),
      "bulk volume/front leading coefficient")


# FD08: isolated common-child block.
for delta, transfer in ((3.0, 0.2), (2.1, 0.7)):
    H = np.array([[0, 0, -transfer],
                  [0, 0, -transfer],
                  [-transfer, -transfer, delta]], dtype=float)
    exact = np.linalg.eigvalsh(H)
    expected = np.array([
        (delta - sqrt(delta * delta + 8 * transfer * transfer)) / 2,
        0.0,
        (delta + sqrt(delta * delta + 8 * transfer * transfer)) / 2,
    ])
    check(np.allclose(exact, expected),
          f"three-mode exact spectrum Delta={delta},t={transfer}")

# FD07: exact low spectral branch and stated operator-norm remainder.
K = B.T @ B
delta = 5.0
transfer = 0.1
evals, evecs = np.linalg.eigh(K)
f_evals = (delta - np.sqrt(delta**2 + 4 * transfer**2 * evals)) / 2
fK = (evecs * f_evals) @ evecs.T
approx = (-(transfer**2 / delta) * K
          + (transfer**4 / delta**3) * (K @ K))
remainder = np.linalg.norm(fK - approx, ord=2)
bound = (2 * abs(transfer)**6 / delta**5
         * np.linalg.norm(K, ord=2)**3)
check(4 * transfer**2 * np.linalg.norm(K, ord=2) / delta**2 <= 0.5,
      "Schur low-branch frozen smallness domain")
check(remainder <= bound * (1 + 1e-9),
      "Schur low-branch operator-norm remainder bound")


# FD12--FD18: bulk symbols and acoustic continuum scaling.
def adjacency_symbol(k, scale):
    return 2 * sum(np.cos(scale * np.dot(k, alpha)) for alpha in roots)


def incidence_symbol(k, scale):
    return abs(sum(np.exp(1j * scale * np.dot(k, n[:, a]))
                   for a in range(4))) ** 2


k = basis_v @ np.array([0.31, -0.17, 0.23])
for scale in (0.2, 0.1, 0.05):
    check(np.isclose(incidence_symbol(k, scale),
                     4 + adjacency_symbol(k, scale)),
          f"incidence/adjacency symbol identity a={scale}")

# Even scalar band has zero first derivative.
eps = 1e-6
direction = basis_v[:, 0]
gradient = (adjacency_symbol(eps * direction, 1.0)
            - adjacency_symbol(-eps * direction, 1.0)) / (2 * eps)
check(abs(gradient) < 1e-10, "scalar sibling band zero linear term")

# Lambda/a^2 approaches (16/3)|k|^2 with O(a^2) error.
def laplacian_symbol(k, scale):
    return 2 * sum(1 - np.cos(scale * np.dot(k, alpha))
                   for alpha in roots)


scales = np.array([0.2, 0.1, 0.05, 0.025])
target = (16 / 3) * np.dot(k, k)
errors = np.array([abs(laplacian_symbol(k, a) / a**2 - target)
                   for a in scales])
ratios = errors[:-1] / errors[1:]
check(np.all((ratios > 3.8) & (ratios < 4.2)),
      "collective cone error is second order in lattice scale")


# Minimal sector/common-cone screen.
P_A = np.ones((6, 6)) / 6
opposites = {0: 5, 5: 0, 1: 4, 4: 1, 2: 3, 3: 2}
O = np.zeros((6, 6))
for i, j in opposites.items():
    O[i, j] = 1
P_T = (np.eye(6) - O) / 2
P_E = (np.eye(6) + O) / 2 - P_A
check(np.allclose(P_A + P_E + P_T, np.eye(6)),
      "pair-sector projectors complete")
check(tuple(round(np.trace(p)) for p in (P_A, P_E, P_T)) == (1, 2, 3),
      "pair-sector dimensions 1+2+3")
check(all(np.allclose(p @ p, p) for p in (P_A, P_E, P_T)),
      "pair-sector projectors idempotent")
check(all(np.allclose(left @ right, 0)
          for i, left in enumerate((P_A, P_E, P_T))
          for j, right in enumerate((P_A, P_E, P_T)) if i != j),
      "pair-sector projectors orthogonal")

chi = (2.0, 3.0, 5.0)
kappa_common = tuple(7.0 * x for x in chi)
speeds = tuple(kap / ch for kap, ch in zip(kappa_common, chi))
check(np.allclose(speeds, speeds[0]), "common-cone ratio criterion")
speeds_split = tuple(kap / ch for kap, ch in zip((14.0, 18.0, 40.0), chi))
check(not np.allclose(speeds_split, speeds_split[0]),
      "S4 permits split sector cones")


# Claim ceilings must remain explicit.
for phrase in (
    "conditional on `DETUNED-Q4-CARRIER-LIFT`",
    "current F3 microscopic parent has already",
    "does not contradict",
    "actual physical volume still requires",
    "Finite `S4` symmetry alone does not force",
    "separately supplied action antecedent",
    "`Q4-PAIR-FIELD-LIFT`",
    "finite nonzero physical speed is an additional scale binding",
    "does not adopt a successor action",
):
    check(phrase in THEOREM, f"claim ceiling: {phrase}")

print(f"SUMMARY {checks}/{checks} PASS")
print("VERDICT EXACT_COMMON_CHILD_AND_A3_ALGEBRA_PASS__COLLECTIVE_PHASE_AND_GRAVITY_REMAIN_OPEN")
