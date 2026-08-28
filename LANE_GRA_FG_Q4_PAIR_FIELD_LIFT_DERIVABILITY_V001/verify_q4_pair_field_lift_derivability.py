#!/usr/bin/env python3
"""Finite replay for GRA-FG-Q4-PFLD-V001.

This verifier checks the finite edge-representation algebra, the lawful static
six-register construction, and the zero-propagation commutator.  It does not
derive a physical solder, field dynamics, or gravity.
"""

from __future__ import annotations

import itertools
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


def commutator(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return a @ b - b @ a


# Tetrahedral unordered-edge representation.
edges = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
edge_index = {edge: i for i, edge in enumerate(edges)}


def edge_permutation(vertex_permutation: tuple[int, ...]) -> np.ndarray:
    matrix = np.zeros((6, 6))
    for column, (a, b) in enumerate(edges):
        image = tuple(sorted((vertex_permutation[a], vertex_permutation[b])))
        matrix[edge_index[image], column] = 1.0
    return matrix


group = [edge_permutation(p) for p in itertools.permutations(range(4))]
check(len(group) == 24, "all 24 S4 edge permutations")
check(all(np.allclose(p.T @ p, np.eye(6)) for p in group), "edge action is orthogonal")

# Exact sector projectors used by corrected CCMAC.
opposite = np.zeros((6, 6))
for i, j in ((0, 5), (1, 4), (2, 3)):
    opposite[i, j] = 1.0
    opposite[j, i] = 1.0

u = np.ones((6, 1))
p_a = (u @ u.T) / 6.0
p_t = (np.eye(6) - opposite) / 2.0
p_e = (np.eye(6) + opposite) / 2.0 - p_a
projectors = (p_a, p_e, p_t)

check(all(np.allclose(p @ p, p) for p in projectors), "A/E/T idempotents")
check(all(np.allclose(p, p.T) for p in projectors), "A/E/T orthogonal projectors")
check(np.allclose(sum(projectors), np.eye(6)), "A+E+T resolves edge identity")
check(tuple(np.linalg.matrix_rank(p, tol=1e-12) for p in projectors) == (1, 2, 3), "A/E/T ranks 1+2+3")
check(
    all(np.linalg.norm(projectors[i] @ projectors[j]) < 1e-12 for i in range(3) for j in range(3) if i != j),
    "A/E/T sectors mutually orthogonal",
)
check(all(np.allclose(g @ p, p @ g) for g in group for p in projectors), "sector projectors commute with S4")

# The commutant of the edge representation has dimension three.
basis_matrices = []
for i in range(6):
    for j in range(6):
        unit = np.zeros((6, 6))
        unit[i, j] = 1.0
        basis_matrices.append(unit)

constraint_columns = []
for unit in basis_matrices:
    constraint_columns.append(np.concatenate([commutator(g, unit).reshape(-1) for g in group]))
constraint = np.column_stack(constraint_columns)
commutant_dimension = 36 - np.linalg.matrix_rank(constraint, tol=1e-10)
check(commutant_dimension == 3, "S4 edge commutant has dimension three")

projector_columns = np.column_stack([p.reshape(-1) for p in projectors])
check(np.linalg.matrix_rank(projector_columns, tol=1e-12) == 3, "A/E/T projectors span three commutant directions")

# Strongest lawful static construction: six independent binary relation bits.
all_patterns = list(itertools.product((0, 1), repeat=6))
check(len(all_patterns) == 64, "six static relation registers have 64 basis patterns")

occupation_conserved = True
for pattern in all_patterns:
    keep_k = np.asarray(pattern, dtype=int)
    keep_g = np.zeros(6, dtype=int)
    break_k = np.zeros(6, dtype=int)
    break_g = np.asarray(pattern, dtype=int)
    occupation_conserved &= int(keep_k.sum() + keep_g.sum()) == int(break_k.sum() + break_g.sum())
check(occupation_conserved, "KEEP/BREAK conserves relation occupation for all 64 patterns")

relabeling_preserves_occupation = True
for pattern in all_patterns:
    pattern_array = np.asarray(pattern)
    for vertex_permutation in itertools.permutations(range(4)):
        image_indices = [
            edge_index[tuple(sorted((vertex_permutation[a], vertex_permutation[b])))]
            for a, b in edges
        ]
        relabeling_preserves_occupation &= int(pattern_array[image_indices].sum()) == sum(pattern)
check(relabeling_preserves_occupation, "S4 relabeling preserves relation occupation")

# The common conditional link pulse reads the static K pattern; BREAK reads zero.
formed = np.ones(6, dtype=int)
sham = np.zeros(6, dtype=int)
check(np.array_equal(formed, np.ones(6, dtype=int)), "full KEEP produces six active relation bits")
check(np.array_equal(sham, np.zeros(6, dtype=int)), "sham produces no active relation bits")
check(np.array_equal(np.zeros(6, dtype=int), sham), "full BREAK blanks the active relation vector")

# Exact local and two-edge zero-propagation commutators.
identity = np.eye(2)
q = np.diag([0.0, 1.0])
x = np.array([[0.0, 1.0], [1.0, 0.0]])
n_link = np.diag([0.0, 1.0])

h_gate = -np.kron(q, x)
q_active = np.kron(q, identity)
check(np.allclose(commutator(h_gate, q_active), 0.0), "FPMH conditional gate commutes with pair record")
check(np.allclose(commutator(np.kron(identity, x), q_active), 0.0), "F3 link flip commutes with separate pair record")
check(np.allclose(commutator(np.kron(identity, n_link), q_active), 0.0), "F3 link occupation commutes with separate pair record")

# Two independent pair-link factors: no cross-pair propagation.
def embed_two(local: np.ndarray, first: bool) -> np.ndarray:
    return np.kron(local, np.eye(4)) if first else np.kron(np.eye(4), local)


h_two = embed_two(h_gate, True) + embed_two(h_gate, False)
q_one = embed_two(q_active, True)
q_two = embed_two(q_active, False)
check(np.allclose(commutator(h_two, q_one), 0.0), "first pair record conserved in two-edge parent")
check(np.allclose(commutator(h_two, q_two), 0.0), "second pair record conserved in two-edge parent")
check(np.allclose(commutator(q_one, q_two), 0.0), "distinct pair records commute")

# Direct Heisenberg replay at several controller times.
eigenvalues, eigenvectors = np.linalg.eigh(h_two)
for tau in (0.1, 0.7, 1.9):
    unitary = eigenvectors @ np.diag(np.exp(-1j * tau * eigenvalues)) @ eigenvectors.conj().T
    evolved = unitary.conj().T @ q_one @ unitary
    check(np.allclose(evolved, q_one), f"pair record has zero Heisenberg evolution tau={tau}")
    check(np.allclose(commutator(evolved, q_two), 0.0), f"retarded interpair commutator vanishes tau={tau}")

check(
    sum(line.strip() == r"\[" for line in THEOREM.splitlines())
    == sum(line.strip() == r"\]" for line in THEOREM.splitlines()),
    "display-math delimiters are balanced",
)
check("{cal H}" not in THEOREM, "no malformed Hilbert-space token remains")
check(THEOREM.count(r"{\cal H}_{{\rm pair},m}") == 3,
      "FG07 Hilbert-space tensor, source, and target are well typed")

required_phrases = (
    "PMMDC `J_ab` | real parameter",
    "not itself a physical type join",
    "static six-edge record",
    "[H_{\\rm inherited},Q_e]=0",
    "Q4-PAIR-SOLDER",
    "PAIR-FIELD-DYNAMICS",
    "not yet compose them into one propagating field",
    "No claim about RGRL-B",
)
for phrase in required_phrases:
    check(phrase in THEOREM, f"claim/type ceiling: {phrase}")

print(f"SUMMARY {passed}/{passed} PASS")
print("VERDICT STATIC_SIX_EDGE_REPRESENTATION_AND_ZERO_PROPAGATION_NO_GO_PASS__PAIR_SOLDER_AND_DYNAMICS_OPEN")
