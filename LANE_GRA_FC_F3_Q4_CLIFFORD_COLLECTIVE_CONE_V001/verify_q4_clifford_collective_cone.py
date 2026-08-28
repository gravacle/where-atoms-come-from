#!/usr/bin/env python3
"""Exact/numerical algebra checks for GRA-FC-Q4CCC-V001.

The script verifies only the finite relational and symbol algebra stated in
THEOREM.md.  It supplies no empirical or physical-instantiation weight.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
THEOREM = (HERE / "THEOREM.md").read_text(encoding="utf-8")

passed = 0


def check(condition: bool, label: str) -> None:
    global passed
    if not condition:
        raise AssertionError(label)
    passed += 1
    print(f"PASS {label}")


# Tetrahedral relational frame.
n = np.array(
    [
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ]
) / math.sqrt(3.0)

gram = n @ n.T
expected_gram = np.full((4, 4), -1.0 / 3.0)
np.fill_diagonal(expected_gram, 1.0)
check(np.allclose(gram, expected_gram), "tetrahedral Gram")
check(np.allclose(n.sum(axis=0), 0.0), "tetrahedral zero sum")
check(np.allclose(n.T @ n, (4.0 / 3.0) * np.eye(3)), "frame second moment")

roots = []
for a in range(4):
    for b in range(a + 1, 4):
        roots.append(n[b] - n[a])
roots = np.asarray(roots)

second = sum(np.outer(alpha, alpha) for alpha in roots)
check(np.allclose(second, (16.0 / 3.0) * np.eye(3)), "root second moment")

# Six dyads span Sym^2(R^3).
dyad_rows = []
for x, y, z in roots:
    dyad_rows.append([x * x, y * y, z * z, 2 * x * y, 2 * x * z, 2 * y * z])
dyad_rows = np.asarray(dyad_rows)
check(np.linalg.matrix_rank(dyad_rows, tol=1e-12) == 6, "six dyads full rank")

# Exact root-basis covolume formula.
basis = np.vstack([n[i] - n[3] for i in range(3)])
basis_gram = basis @ basis.T
covolume = math.sqrt(float(np.linalg.det(basis_gram)))
check(
    math.isclose(covolume, 16.0 / (3.0 * math.sqrt(3.0)), rel_tol=1e-12),
    "A3 relational covolume",
)

# Existing scalar even hopping has zero linear term and the predicted Hessian.
def scalar_even(k: np.ndarray, a_step: float = 0.2, t: float = 1.0) -> float:
    return 2.0 * t * sum(1.0 - math.cos(a_step * float(k @ alpha)) for alpha in roots)


h = 1e-5
gradient = np.array(
    [
        (scalar_even(h * np.eye(3)[i]) - scalar_even(-h * np.eye(3)[i])) / (2.0 * h)
        for i in range(3)
    ]
)
check(np.linalg.norm(gradient) < 1e-10, "even scalar stencil zero linear term")

scalar_hessian = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        ei = np.eye(3)[i]
        ej = np.eye(3)[j]
        scalar_hessian[i, j] = (
            scalar_even(h * (ei + ej))
            - scalar_even(h * (ei - ej))
            - scalar_even(h * (-ei + ej))
            + scalar_even(-h * (ei + ej))
        ) / (4.0 * h * h)
expected_hessian = 2.0 * (0.2**2) * second
check(np.allclose(scalar_hessian, expected_hessian, rtol=2e-5, atol=2e-7), "even scalar stencil quadratic Hessian")

# Clifford symbol.
sigma = np.array(
    [
        [[0, 1], [1, 0]],
        [[0, -1j], [1j, 0]],
        [[1, 0], [0, -1]],
    ],
    dtype=complex,
)


def gamma(v: np.ndarray) -> np.ndarray:
    return sum(float(v[i]) * sigma[i] for i in range(3))


for i in range(3):
    for j in range(3):
        anticom = gamma(np.eye(3)[i]) @ gamma(np.eye(3)[j]) + gamma(np.eye(3)[j]) @ gamma(np.eye(3)[i])
        check(np.allclose(anticom, 2.0 * (i == j) * np.eye(2)), f"Clifford pair {i}{j}")


def front_vector(k: np.ndarray, a_step: float) -> np.ndarray:
    return (3.0 / (16.0 * a_step)) * sum(
        alpha * math.sin(a_step * float(k @ alpha)) for alpha in roots
    )


k = np.array([0.31, -0.22, 0.17])
errors = []
for a_step in (0.2, 0.1, 0.05, 0.025):
    f = front_vector(k, a_step)
    errors.append(np.linalg.norm(f - k))
    H = gamma(f)
    check(np.allclose(H, H.conj().T), f"Hermitian Bloch symbol a={a_step}")
    check(np.allclose(H @ H, float(f @ f) * np.eye(2)), f"exact Clifford square a={a_step}")

ratios = [errors[i] / errors[i + 1] for i in range(len(errors) - 1)]
check(all(3.8 < ratio < 4.2 for ratio in ratios), "infrared error is second order in lattice scale")

# Pair-weight variations span the complete first-order co-metric tangent.
c0 = 3.0 / 16.0
B0 = sum(c0 * np.outer(alpha, alpha) for alpha in roots)
check(np.allclose(B0, np.eye(3)), "symmetric pair weights give identity principal map")

metric_variations = []
for alpha in roots:
    dB = np.outer(alpha, alpha)
    dQ = B0.T @ dB + dB.T @ B0
    metric_variations.append(
        [dQ[0, 0], dQ[1, 1], dQ[2, 2], dQ[0, 1], dQ[0, 2], dQ[1, 2]]
    )
check(np.linalg.matrix_rank(np.asarray(metric_variations), tol=1e-12) == 6, "pair weights span co-metric tangent")

required_phrases = (
    "bypasses rather than changes",
    "RGRL-B remains open",
    "does not derive a constitutive",
    "does not prove that nature instantiates",
)
for phrase in required_phrases:
    check(phrase in THEOREM, f"claim ceiling: {phrase}")

print(f"SUMMARY {passed}/{passed} PASS")
print("VERDICT EXACT_Q4_RELATIONAL_AND_CLIFFORD_SYMBOL_ALGEBRA_PASS__PHYSICAL_VERTEX_AND_GRAVITY_OPEN")
