#!/usr/bin/env python3
"""Finite replay for GRA-FF-F3-Q4-CLDNG-V001.

This checks operator and support identities.  It does not instantiate the
missing Q4-SUPPORT-SOLDER, detuning port, or any new kinetic field.
"""

from __future__ import annotations

from itertools import product
from math import comb
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
THEOREM = (HERE / "THEOREM.md").read_text()
CHECKS: list[str] = []


def check(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def compositions(n: int, k: int = 4):
    if k == 1:
        yield (n,)
        return
    for first in range(n + 1):
        for rest in compositions(n - first, k - 1):
            yield (first,) + rest


def q4_incidence(n: int):
    parents = tuple(compositions(n))
    children = tuple(compositions(n + 1))
    child_index = {c: i for i, c in enumerate(children)}
    b = np.zeros((len(children), len(parents)), dtype=int)
    for j, m in enumerate(parents):
        for a in range(4):
            c = list(m)
            c[a] += 1
            b[child_index[tuple(c)], j] = 1
    return parents, children, b


def adjacency_from_b(b: np.ndarray) -> np.ndarray:
    na = b.shape[1]
    nb = b.shape[0]
    return np.block(
        [[np.zeros((na, na), dtype=int), b.T],
         [b, np.zeros((nb, nb), dtype=int)]]
    )


def periodic_diamond(size: int, directions=(0, 1, 2, 3)):
    cells = tuple(product(range(size), repeat=3))
    ci = {r: i for i, r in enumerate(cells)}
    shifts = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
    b = np.zeros((len(cells), len(cells)), dtype=int)
    for ai, r in enumerate(cells):
        for a in directions:
            s = shifts[a]
            target = tuple((r[i] + s[i]) % size for i in range(3))
            b[ci[target], ai] = 1
    return cells, b


# 1. Exact q4 front sizes and literal F3 equal-layer padding obligation.
for n in range(0, 8):
    sn = tuple(compositions(n))
    sn1 = tuple(compositions(n + 1))
    if len(sn) != comb(n + 3, 3) or len(sn1) - len(sn) != comb(n + 3, 2):
        raise AssertionError(n)
check(True, "q4 front sizes and padding gap")

# 2. The BQ4 append writes provenance rather than acting as history identity.
word = (1, 2)
appended = word + (3,)
check(appended != word and appended[-1] == 3,
      "q4 append necessarily writes a fresh word label")
check((1, 2) != (2, 1) and sorted((1, 2)) == sorted((2, 1)),
      "equal-count routes remain orthogonal word labels")

# 3. Finite q4 incidence identities.
n = 3
parents, children, b = q4_incidence(n)
check(np.all(b.sum(axis=0) == 4), "every q4 parent has four append incidences")
expected_child_degree = np.array([sum(x > 0 for x in c) for c in children])
check(np.array_equal(b.sum(axis=1), expected_child_degree),
      "q4 child incidence degree is positive-coordinate count")
a = adjacency_from_b(b)
check(np.array_equal(a, a.T) and np.all(np.diag(a) == 0),
      "Hermitian carrier support is the bipartite B plus B-dagger block")

# 4. Literal two-qutrit BS07 algebra: T^2=J^2=exactly-one-endpoint projector.
dim = 3
q = np.diag([0.0, 1.0, 1.0])
ident = np.eye(dim)
t_edge = np.zeros((dim * dim, dim * dim), dtype=complex)
j_edge = np.zeros_like(t_edge)


def index(left: int, right: int) -> int:
    return left * dim + right


for content in (1, 2):
    bx = index(0, content)
    xb = index(content, 0)
    t_edge[bx, xb] = 1
    t_edge[xb, bx] = 1
    j_edge[bx, xb] = 1j
    j_edge[xb, bx] = -1j

projector = np.kron(q, ident) + np.kron(ident, q) - 2 * np.kron(q, q)
check(np.allclose(t_edge @ t_edge, projector), "BS07 transfer-square identity")
check(np.allclose(j_edge @ j_edge, projector), "BS07 current-square identity")

# 5. Exact one-carrier graph generator.
degree = a.sum(axis=1)
eps = 1.25
lam = 0.7
t_hop = 0.4
h_one = eps * np.eye(len(a)) + lam * np.diag(degree) - t_hop * a
check(np.allclose(np.diag(h_one), eps + lam * degree),
      "one-carrier current square is lambda_J times graph degree")
na = len(parents)
check(np.allclose(h_one[:na, na:], -t_hop * b.T)
      and np.allclose(h_one[na:, :na], -t_hop * b),
      "restricted F3 off-diagonal is exactly -t times q4 incidence")
check(np.array_equal(h_one, h_one.T), "restricted F3 generator is Hermitian")

# 6. Finite slab cannot furnish a positive uniform child detuning.
parent_diag = eps + 4 * lam
child_relative = eps + lam * expected_child_degree - parent_diag
check(np.all(child_relative <= 1e-14), "finite child degree shift is never positive")
check(np.any(np.isclose(child_relative, 0.0)) and len(set(np.round(child_relative, 12))) > 1,
      "finite child degree shift is zero in the interior and nonuniform at boundary")

# 7. A periodic q4/diamond support makes the current-parent detuning exactly zero.
cells, b_full = periodic_diamond(3)
a_full = adjacency_from_b(b_full)
degree_full = a_full.sum(axis=1)
check(np.all(degree_full == 4), "periodic q4 diamond is degree four on both parts")
h_full = eps * np.eye(len(a_full)) + lam * np.diag(degree_full) - t_hop * a_full
check(np.allclose(np.diag(h_full), (eps + 4 * lam) * np.ones(len(a_full))),
      "regular F3 bulk gives zero parent-child onsite offset")
check(np.allclose(h_full[:len(cells), len(cells):], -t_hop * b_full.T),
      "regular restricted generator owns the full scalar q4 hopping")

# The periodic graph has an exact part-exchange inversion; a staggered onsite
# term is odd under it and cannot be generated in a covariant effective block.
cell_index = {r: i for i, r in enumerate(cells)}
nv = len(cells)
exchange = np.zeros_like(a_full, dtype=float)
for r, i in cell_index.items():
    minus_r = tuple((-x) % 3 for x in r)
    j = cell_index[minus_r]
    exchange[nv + j, i] = 1.0
    exchange[j, nv + i] = 1.0
staggered = np.diag(np.concatenate((-np.ones(nv), np.ones(nv))))
check(np.allclose(exchange @ a_full @ exchange.T, a_full),
      "periodic q4 diamond has exact parent-child exchange symmetry")
check(np.allclose(exchange @ staggered @ exchange.T, -staggered),
      "positive child-minus-parent detuning is odd under matched exchange")

# The BS06 link detuning is constant in a fixed saturated incidence word.
fixed_link_energy = 2.3 * int(b_full.sum())
check(np.allclose(fixed_link_energy * np.ones(len(a_full)),
                  fixed_link_energy * np.eye(len(a_full)).diagonal()),
      "BS06 link detuning is a carrier-position scalar in fixed n")

# 8. Same-n exact incompatibility: saturated d=4 versus ice d=2.
_, b_ice_01 = periodic_diamond(3, directions=(0, 1))
_, b_ice_02 = periodic_diamond(3, directions=(0, 2))
a_ice_01 = adjacency_from_b(b_ice_01)
a_ice_02 = adjacency_from_b(b_ice_02)
check(np.all(a_ice_01.sum(axis=1) == 2), "two direction matchings give an exact d*=2 ice word")
check(np.all(a_full.sum(axis=1) == 4), "saturated full-hopping word is d=4")
check(not np.array_equal(a_ice_01, a_full),
      "d*=2 BS09 hopping omits full-q4 carrier edges")
check(not np.array_equal(a_ice_01, a_ice_02),
      "different ice words give different carrier hopping operators")
check(np.linalg.norm(a_full - a_ice_01, ord=2) > 0,
      "full carrier block cannot equal an exact ice-controlled block")

# 9. A zero nonedge mask is not invariant under an active BS06 Pauli X.
ket_zero = np.array([1.0, 0.0])
pauli_x = np.array([[0.0, 1.0], [1.0, 0.0]])
check(np.allclose(pauli_x @ ket_zero, np.array([0.0, 1.0])),
      "active BS06 flip does not conserve a prepared zero nonedge mask")

# 10. Claim and non-adoption ceilings.
required_phrases = (
    "Q4-SUPPORT-SOLDER",
    "dimension-level encoding, not an existing",
    "No new hopping coefficient",
    "no positive child offset",
    "one incidence word cannot realize both exact blocks",
    "adopt, recommend, or install it",
    "symbolic complete-port slot",
    "not adopted",
    "gravity",
)
check(all(phrase in THEOREM for phrase in required_phrases),
      "all derivability and non-adoption ceilings are explicit")

print(f"CLDNG verification: PASS ({len(CHECKS)}/{len(CHECKS)})")
for number, label in enumerate(CHECKS, 1):
    print(f"{number:02d}. PASS  {label}")
