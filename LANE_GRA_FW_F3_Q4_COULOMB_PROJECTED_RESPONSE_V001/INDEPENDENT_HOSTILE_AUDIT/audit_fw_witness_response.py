#!/usr/bin/env python3
"""Independent hostile replay of the frozen FW FV-WITNESS response.

This executable does not import or execute the FW builder or the FO builder.
It independently rebuilds the cyclic diamond quotient with frozenset states,
enumerates the 180-state ring component, constructs the two explicitly scoped
FV-WITNESS source pieces, and evaluates the exact k=0 response in the number
field Q(sqrt(2),sqrt(3)).
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
LANE = Path(__file__).resolve().parent.parent

CORE = {
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/DEPENDENCIES.sha256":
        "4b503ff10457cde42c2fe7acbcb60c87ba4befef2736a0674c2ae5b6667c92e5",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/README.md":
        "f9d22ebfe2641994c07a8ec9772b2edb4dd5891843afbcefe8d47c99bce63ea1",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/RESULT.md":
        "81d8f732d8395d757c5405c11c093156bf3a1c2dfae4670cda2db061a5c7e262",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/SELF_AUDIT.md":
        "66c020fd7f033ba3a51c6408a4f750cdf02bb04b681d5f792480c39e04123917",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/THEOREM.md":
        "db3e12d50fd1cb41cddc722a0445cdeaef6a52d49704fa6df1028dfd9abcba1b",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/VERIFICATION.txt":
        "1350c262df70da2f501c89b121e767160f2943e7466ecc5f6d92f1c4f75381d6",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/verify_projected_response.py":
        "87152814b07eeef30794626a003adc8b97f16eeb79252ee99a94d307e822aad8",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/MANIFEST.sha256":
        "17499bb43606bd110657218585bb9b8bebe73358cb49da6276ea8ae7d42d0c47",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/SEAL.sha256":
        "1f829571f82aa16eebac53492470b6a0a883b8401a198ee137a34f9a6b052c9b",
}

DEPENDENCIES = {
    "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001/THEOREM.md":
        "44fc28edc9820d2b4ea67cef9f83beef60e53bfe6407b582ffbeccfe42f756c5",
    "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001/INDEPENDENT_AUDIT.md":
        "84d8c02c3e560198f6a9fae04f5ee81bc72c98354e02f1f58d506a8c3171c453",
    "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001/verify_finite_tt_four_point.py":
        "fb44d45290c0530098c0e8f9593dff1c0f8149d42598f842374b61004a8ff6c2",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/THEOREM.md":
        "6fc221a31151340b91a946d33e442971c1373500e067c354b6c610e3964edb1c",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/RESULT.md":
        "b5d4c3de99aa4e100519c19a9b74de487b47c1a2d3671204e77740bd9094771a",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/verify_projected_source_rank.py":
        "0e93d84f9eb7cf7fdd62b5a14d5c6705c74841899dd1676bdd7e7a41eb971a00",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/MANIFEST.sha256":
        "651a66b9afd7545b04aa80e5f90952fda9327d011ecd19b973aa80ff51a739f3",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/SEAL.sha256":
        "8301da6bbc026d0e14d985592c5dabe3d91072957c4aaa4b1bebf1f45aadd894",
}


checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


# -------------------------------------------------------------------------
# Minimal exact field Q(sqrt(2), sqrt(3)).


@dataclass(frozen=True)
class K:
    """a + b sqrt(2) + c sqrt(3) + d sqrt(6), with rational coefficients."""

    a: F = F(0)
    b: F = F(0)
    c: F = F(0)
    d: F = F(0)

    @staticmethod
    def coerce(value) -> "K":
        if isinstance(value, K):
            return value
        return K(F(value))

    def __add__(self, other):
        other = K.coerce(other)
        return K(self.a + other.a, self.b + other.b,
                 self.c + other.c, self.d + other.d)

    __radd__ = __add__

    def __neg__(self):
        return K(-self.a, -self.b, -self.c, -self.d)

    def __sub__(self, other):
        return self + (-K.coerce(other))

    def __rsub__(self, other):
        return K.coerce(other) - self

    def __mul__(self, other):
        other = K.coerce(other)
        a, b, c, d = self.a, self.b, self.c, self.d
        e, f, g, h = other.a, other.b, other.c, other.d
        return K(
            a*e + 2*b*f + 3*c*g + 6*d*h,
            a*f + b*e + 3*c*h + 3*d*g,
            a*g + c*e + 2*b*h + 2*d*f,
            a*h + d*e + b*g + c*f,
        )

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, (int, F)):
            return K(self.a / other, self.b / other,
                     self.c / other, self.d / other)
        raise TypeError("only rational division is used in this audit")

    def is_zero(self) -> bool:
        return self == K()

    def value(self) -> float:
        return (float(self.a) + float(self.b)*2**0.5 +
                float(self.c)*3**0.5 + float(self.d)*6**0.5)


ONE = K(1)
ZERO = K()
S2 = K(b=F(1))
S3 = K(c=F(1))
S6 = K(d=F(1))
IS2 = S2 / 2
IS3 = S3 / 3
IS6 = S6 / 6
check(S2*S2 == K(2) and S3*S3 == K(3) and S6*S6 == K(6) and
      S2*S3 == S6 and S2*S6 == 2*S3 and S3*S6 == 3*S2,
      "exact Q(sqrt2,sqrt3) multiplication table is internally certified")


def kmat(rows):
    return [[K.coerce(value) for value in row] for row in rows]


def kzeros(n, m):
    return [[ZERO for _ in range(m)] for _ in range(n)]


def keye(n):
    return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]


def kadd(a, b):
    return [[x+y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]


def ksub(a, b):
    return [[x-y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]


def kscale(s, a):
    s = K.coerce(s)
    return [[s*x for x in row] for row in a]


def ktranspose(a):
    return [list(row) for row in zip(*a)]


def kmul(a, b):
    bt = ktranspose(b)
    return [[sum((x*y for x, y in zip(row, column)), ZERO)
             for column in bt] for row in a]


def kmatvec(a, v):
    return [sum((x*y for x, y in zip(row, v)), ZERO) for row in a]


def kdot(a, b):
    return sum((x*y for x, y in zip(a, b)), ZERO)


def kouter(a, b):
    return [[x*y for y in b] for x in a]


def ktrace(a):
    return sum((a[i][i] for i in range(len(a))), ZERO)


def kiszero_matrix(a):
    return all(x.is_zero() for row in a for x in row)


def kmatrix_equal(a, b):
    return kiszero_matrix(ksub(a, b))


def rational_rank(rows) -> int:
    matrix = [[F(value) for value in row] for row in rows]
    if not matrix:
        return 0
    row = 0
    for column in range(len(matrix[0])):
        pivot = next((i for i in range(row, len(matrix))
                      if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        value = matrix[row][column]
        matrix[row] = [x/value for x in matrix[row]]
        for i in range(len(matrix)):
            if i == row or not matrix[i][column]:
                continue
            factor = matrix[i][column]
            matrix[i] = [x-factor*y for x, y in zip(matrix[i], matrix[row])]
        row += 1
        if row == len(matrix):
            break
    return row


def rational_det(rows) -> F:
    matrix = [[F(value) for value in row] for row in rows]
    determinant = F(1)
    for column in range(len(matrix)):
        pivot = next((i for i in range(column, len(matrix))
                      if matrix[i][column]), None)
        if pivot is None:
            return F(0)
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            determinant *= -1
        value = matrix[column][column]
        determinant *= value
        for i in range(column+1, len(matrix)):
            factor = matrix[i][column] / value
            for j in range(column+1, len(matrix)):
                matrix[i][j] -= factor*matrix[column][j]
    return determinant


def polynomial_mul(a, b):
    result = [F(0)] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def determinant_polynomial_xI_minus(matrix):
    """Leibniz determinant, coefficients in ascending powers of x."""
    n = len(matrix)
    result = [F(0)] * (n+1)
    for permutation in permutations(range(n)):
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(n) for j in range(i+1, n))
        term = [F(1)]
        for i, j in enumerate(permutation):
            factor = [-F(matrix[i][j])]
            if i == j:
                factor.append(F(1))
            term = polynomial_mul(term, factor)
        sign = -1 if inversions % 2 else 1
        for degree, coefficient in enumerate(term):
            result[degree] += sign*coefficient
    return result


# -------------------------------------------------------------------------
# Custody first: repaired FV-PURE and frozen FW-WITNESS bytes.


for relative, expected in {**CORE, **DEPENDENCIES}.items():
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(),
          f"regular-file custody: {relative}")
    check(digest(path) == expected, f"frozen hash custody: {relative}")

fw_manifest = (LANE / "MANIFEST.sha256").read_text().splitlines()
fw_seal = (LANE / "SEAL.sha256").read_text().splitlines()
check(any(line.startswith(CORE[
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/THEOREM.md"])
    and line.endswith("  THEOREM.md") for line in fw_manifest),
    "FW manifest binds the frozen witness theorem")
check(fw_seal == [CORE[
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/MANIFEST.sha256"]
    + "  MANIFEST.sha256",
    CORE["LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/VERIFICATION.txt"]
    + "  VERIFICATION.txt"], "FW seal binds manifest and replay transcript")

fv_theorem = (ROOT / "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/THEOREM.md").read_text()
fv_audit = (ROOT / "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md").read_text()
check("S10 / FV-PURE" in fv_theorem and
      "FV-PURE` is stronger than FU `S1`--`S9`" in fv_theorem,
      "repaired FV dependency exposes S10/FV-PURE as an added premise")
check("PASS_AFTER_FV_PURE_PREMISE_AND_BYTE_HYGIENE_REPAIR" in fv_audit,
      "FV independent audit disposition belongs to the repaired core")


# -------------------------------------------------------------------------
# Independent frozenset construction of the FO quotient and component.


CELL_COUNT = 30
SHIFTS = (0, 1, 5, 19)
LINK_COUNT = 120
VERTEX_COUNT = 60

edges = []
labels = []
incidence = defaultdict(list)
for cell in range(CELL_COUNT):
    for label, shift in enumerate(SHIFTS):
        edge = len(edges)
        a = cell
        b = CELL_COUNT + (cell + shift) % CELL_COUNT
        edges.append((a, b))
        labels.append(label)
        incidence[a].append((b, edge))
        incidence[b].append((a, edge))

check(len(edges) == LINK_COUNT and len(incidence) == VERTEX_COUNT,
      "independent quotient has 120 links and 60 vertices")
check(all(len(incidence[v]) == 4 for v in range(VERTEX_COUNT)),
      "independent quotient is coordination four")


def canonical_cycle(order):
    order = tuple(order)
    images = []
    for direction in (order, tuple(reversed(order))):
        images.extend(direction[k:]+direction[:k] for k in range(len(order)))
    return min(images)


cycles = set()
for root in range(VERTEX_COUNT):
    stack = [(root, (root,), ())]
    while stack:
        vertex, visited_vertices, visited_edges = stack.pop()
        if len(visited_edges) == 6:
            if vertex == root:
                cycles.add(canonical_cycle(visited_edges))
            continue
        for neighbor, edge in incidence[vertex]:
            if edge in visited_edges:
                continue
            if neighbor == root:
                if len(visited_edges) == 5:
                    stack.append((neighbor, visited_vertices+(neighbor,),
                                  visited_edges+(edge,)))
            elif neighbor not in visited_vertices:
                stack.append((neighbor, visited_vertices+(neighbor,),
                              visited_edges+(edge,)))

cycles = tuple(sorted(cycles))
check(len(cycles) == 120, "independent simple-cycle census gives 120 hexagons")
missing_counts = Counter()
patterns = []
for cycle in cycles:
    counts = Counter(labels[e] for e in cycle)
    if not (len(counts) == 3 and sorted(counts.values()) == [2, 2, 2]):
        raise AssertionError("non-elementary cycle entered independent census")
    missing = next(iter(set(range(4))-set(counts)))
    missing_counts[missing] += 1
    first = frozenset(cycle[::2])
    second = frozenset(cycle[1::2])
    patterns.append((frozenset(cycle), first, second, missing))
check(all(len(Counter(labels[e] for e in cycle)) == 3 and
          sorted(Counter(labels[e] for e in cycle).values()) == [2, 2, 2]
          for cycle in cycles),
      "every independent cycle has one missing q4 label")
check(missing_counts == Counter({0: 30, 1: 30, 2: 30, 3: 30}),
      "all four missing-label ring orientations occur thirty times")


def state_key(state):
    return sum(1 << edge for edge in state)


def degree_tuple(state):
    return tuple(sum(edge in state for _, edge in incidence[v])
                 for v in range(VERTEX_COUNT))


def translate(state, amount=1):
    return frozenset(4*((edge//4+amount) % CELL_COUNT) + edge % 4
                     for edge in state)


base = frozenset(edge for edge, label in enumerate(labels)
                 if label in (0, 1))
seed_loop = frozenset((84, 11, 9, 114, 112, 39, 37, 87))
seed = base.symmetric_difference(seed_loop)
check(all(d == 2 for d in degree_tuple(base)) and
      all(d == 2 for d in degree_tuple(seed)),
      "independent frozen and winding-seed states obey exact ice")

queue = deque([seed])
component = {seed}
transition_counter = 0
while queue:
    state = queue.popleft()
    for cycle, first, second, _ in patterns:
        restriction = state & cycle
        if restriction != first and restriction != second:
            continue
        transition_counter += 1
        successor = state.symmetric_difference(cycle)
        if successor not in component:
            component.add(successor)
            queue.append(successor)

states = tuple(sorted(component, key=state_key))
index = {state: i for i, state in enumerate(states)}
check(len(states) == 180, "independent winding component has 180 states")
check(transition_counter == 840, "independent component has 420 undirected ring transitions")
check(all(all(d == 2 for d in degree_tuple(state)) for state in states),
      "every independently generated component state remains in ice")
check(all(translate(state) in component for state in states),
      "independent component is closed under cyclic translation")

full = frozenset(range(LINK_COUNT))
complement = {state.symmetric_difference(full) for state in component}
check(component.isdisjoint(complement) and len(complement) == 180,
      "selected winding component excludes its complement partner")

H = np.zeros((180, 180), dtype=np.int64)
ring_missing = {}
for row, state in enumerate(states):
    for cycle, first, second, missing in patterns:
        restriction = state & cycle
        if restriction != first and restriction != second:
            continue
        column = index[state.symmetric_difference(cycle)]
        H[row, column] -= 1
        ring_missing[(row, column)] = missing
check(np.array_equal(H, H.T) and np.count_nonzero(np.triu(H)) == 420,
      "independent H6 matrix is symmetric with 420 edges")
check(np.min(H) == -1 and np.max(H) == 0,
      "no duplicate ring path changes an H6 matrix element")
translation_permutation = np.array([index[translate(state)] for state in states])
check(np.array_equal(H, H[np.ix_(translation_permutation,
                                 translation_permutation)]),
      "independent H6 matrix commutes exactly with cyclic translation")

unassigned = set(states)
orbits = []
while unassigned:
    representative = min(unassigned, key=state_key)
    orbit = []
    current = representative
    while current not in orbit:
        orbit.append(current)
        current = translate(current)
    orbits.append(tuple(orbit))
    unassigned.difference_update(orbit)
check(len(orbits) == 6 and all(len(orbit) == 30 for orbit in orbits),
      "independent translation action has six free length-thirty orbits")

orbit_of = {state: i for i, orbit in enumerate(orbits) for state in orbit}
H0 = np.zeros((6, 6), dtype=np.int64)
for i, orbit in enumerate(orbits):
    row = index[orbit[0]]
    for column, value in enumerate(H[row]):
        if value:
            H0[i, orbit_of[states[column]]] += value

EXPECTED_H0 = np.array((
    (0, -1, -1, -1, -1, -2),
    (-1, 0, -1, -1, 0, -1),
    (-1, -1, 0, 0, -1, -1),
    (-1, -1, 0, 0, -1, -1),
    (-1, 0, -1, -1, 0, -1),
    (-2, -1, -1, -1, -1, 0),
), dtype=np.int64)
check(np.array_equal(H0, EXPECTED_H0),
      "independent orbit quotient gives the exact FW 6x6 H0")

expected_charpoly = polynomial_mul(
    polynomial_mul([F(0), F(0), F(1)], [F(4), F(-4), F(1)]),
    [F(-4), F(4), F(1)])
check(determinant_polynomial_xI_minus(H0.tolist()) == expected_charpoly,
      "exact H0 characteristic polynomial is x^2(x-2)^2(x^2+4x-4)")

H0K = kmat(H0.tolist())
I6 = keye(6)
g = [K(F(1, 2)), S2/4, S2/4, S2/4, S2/4, K(F(1, 2))]
E0 = -2-2*S2
check(kdot(g, g) == ONE and kmatvec(H0K, g) == [E0*x for x in g],
      "closed algebraic ground vector has E0=-2-2sqrt2")

numerical_eigenvalues = np.linalg.eigvalsh(H.astype(float))
check(abs(numerical_eigenvalues[0]-E0.value()) < 2e-12 and
      numerical_eigenvalues[1]-numerical_eigenvalues[0] > 0.36,
      "full independent 180-state spectrum has the same unique sector ground")


# -------------------------------------------------------------------------
# Independent FV family witness and restricted FV-WITNESS source.


family_rows = [
    [F(-1), F(1), F(0), F(0), F(0), F(0)],
    [F(-1), F(0), F(1), F(0), F(0), F(0)],
]
for sx, sy, sz in ((1, 1, 1), (1, -1, -1),
                   (-1, 1, -1), (-1, -1, 1)):
    family_rows.append([
        F(231, 8), F(231, 8), F(231, 8),
        F(-189*sx*sy, 8), F(-189*sx*sz, 8), F(-189*sy*sz, 8),
    ])
check(rational_rank(family_rows) == 6,
      "independent FV family matrix-element witness has off-shell rank six")
check(rational_det(family_rows) == F(-4678629417, 256),
      "independent FV witness determinant is -4678629417/256")

# D16 is sixteen times Q_pair/Ud.  R3 is three times Q_ring/J6.
SIGNS = ((1, 1, 1), (1, -1, -1),
         (-1, 1, -1), (-1, -1, 1))
PAIRS = tuple(combinations(range(4), 2))
ROOT_NUM = {}
for a, b in PAIRS:
    x, y, z = (SIGNS[b][i]-SIGNS[a][i] for i in range(3))
    ROOT_NUM[(a, b)] = (x*x, y*y, z*z, 2*x*y, 2*x*z, 2*y*z)

D16 = np.zeros((6, 180, 180), dtype=np.int64)
for row, state in enumerate(states):
    for vertex in range(VERTEX_COUNT):
        z_by_label = {labels[edge]: (1 if edge not in state else -1)
                      for _, edge in incidence[vertex]}
        for a, b in PAIRS:
            product = z_by_label[a]*z_by_label[b]
            D16[:, row, row] -= product*np.array(ROOT_NUM[(a, b)], dtype=np.int64)

R3 = np.zeros((6, 180, 180), dtype=np.int64)
for (row, column), missing in ring_missing.items():
    sx, sy, sz = SIGNS[missing]
    R3[:, row, column] = (-11, -11, -11,
                           9*sx*sy, 9*sx*sz, 9*sy*sz)

check(np.array_equal(D16, D16.transpose(0, 2, 1)),
      "independent direct source is exactly diagonal Hermitian")
check(np.array_equal(R3, R3.transpose(0, 2, 1)) and
      np.count_nonzero(R3[0]) == 840,
      "independent irreducible ring source is Hermitian on all 420 transitions")
check(all(np.array_equal(source,
                         source[np.ix_(translation_permutation,
                                       translation_permutation)])
          for source in tuple(D16)+tuple(R3)),
      "all twelve independently constructed source pieces commute with translation")
check(not np.any(D16[3:]) and
      np.array_equal(R3[0], R3[1]) and np.array_equal(R3[1], R3[2]),
      "direct and ring pieces separate exactly into A1+E and A1+T2")

# Exact normalization attack: FV11 times -8/63 gives FW04.
for missing, (sx, sy, sz) in enumerate(SIGNS):
    fv11 = [F(231, 8), F(231, 8), F(231, 8),
            F(-189*sx*sy, 8), F(-189*sx*sz, 8), F(-189*sy*sz, 8)]
    converted = [F(-8, 63)*x for x in fv11]
    fw04 = [F(-11, 3), F(-11, 3), F(-11, 3),
            F(3*sx*sy), F(3*sx*sz), F(3*sy*sz)]
    check(converted == fw04,
          f"FV11-to-J6 normalization is exact for missing label {missing}")

I180 = np.eye(180, dtype=np.int64)
check(np.array_equal(D16[0]+D16[1]+D16[2], 960*I180),
      "direct A1 sum is exactly 60 rho I after division by sixteen")
check(np.array_equal(R3[0]+R3[1]+R3[2], 33*H),
      "ring A1 sum is exactly 11 H after division by three")
check(np.array_equal(2*D16[0]-D16[1]-D16[2], 1536*I180) and
      not np.any(2*R3[0]-R3[1]-R3[2]),
      "component Ec combination is exactly 16sqrt6 rho times identity")
check(np.array_equal(H@R3[5], R3[5]@H),
      "component Tyz witness direction exactly commutes with H")


def gram_int(matrices):
    flattened = np.asarray(matrices, dtype=np.int64).reshape(len(matrices), -1)
    return flattened @ flattened.T


def rank_from_gram(matrices):
    return rational_rank(gram_int(matrices).tolist())


check(rank_from_gram([16*I180, *D16]) == 2,
      "direct component source has rank one modulo identity")
check(rank_from_gram(R3) == 4,
      "ring component source has exact A1+T2 rank four")
for rho in (1, 2):
    # 48 Q = 3 rho D16 + 16 R3.
    Q48 = 3*rho*D16 + 16*R3
    check(rank_from_gram([48*I180, *Q48]) == 6,
          f"full component FV-WITNESS span has rank five modulo identity at rho={rho}")


# Exact commutator Gram split before tensor-basis contraction.
CD = np.array([H@operator-operator@H for operator in D16])
CR = np.array([H@operator-operator@H for operator in R3])
GDD = gram_int(CD)
GRR = gram_int(CR)
GDR = CD.reshape(6, -1) @ CR.reshape(6, -1).T

B = kzeros(6, 6)
B[0][0] = B[0][1] = B[0][2] = IS3
B[1][0], B[1][1] = IS2, -IS2
B[2][0], B[2][1], B[2][2] = IS6, IS6, -2*IS6
B[3][3] = B[4][4] = B[5][5] = IS2


def transform_rational_gram(matrix, denominator):
    exact = kmat([[F(int(x), denominator) for x in row] for row in matrix])
    return kmul(kmul(B, exact), ktranspose(B))


direct_comm = transform_rational_gram(GDD, 16*16)
ring_comm = transform_rational_gram(GRR, 3*3)
cross_comm = transform_rational_gram(GDR+GDR.T, 16*3)
expected_direct_comm = kzeros(6, 6)
expected_direct_comm[1][1] = K(960)
expected_direct_comm[1][2] = expected_direct_comm[2][1] = -960*S3
expected_direct_comm[2][2] = K(2880)
expected_ring_comm = kzeros(6, 6)
expected_ring_comm[3][3] = expected_ring_comm[4][4] = K(25920)
check(kmatrix_equal(direct_comm, expected_direct_comm),
      "exact direct commutator Gram is the rank-one E block")
check(kmatrix_equal(ring_comm, expected_ring_comm),
      "exact ring commutator Gram contains only Txy and Txz")
check(kiszero_matrix(cross_comm),
      "direct/ring commutator cross Gram vanishes exactly")
check(expected_direct_comm[1][1]*expected_direct_comm[2][2] ==
      expected_direct_comm[1][2]*expected_direct_comm[2][1] and
      not expected_direct_comm[1][1].is_zero(),
      "exact E commutator block has rank one")
check(not expected_ring_comm[3][3].is_zero() and
      not expected_ring_comm[4][4].is_zero(),
      "exact ring commutator block has two orthogonal nonzero directions")
check(kiszero_matrix(cross_comm),
      "rank-one E plus rank-two T2 proves adH rank three for rho nonzero")


# -------------------------------------------------------------------------
# Exact 6x6 k=0 source blocks and projector response, no eigenvectors.


def reduce_source(source):
    reduced = np.zeros((6, 6, 6), dtype=np.int64)
    for component_index in range(6):
        for i, orbit in enumerate(orbits):
            row = index[orbit[0]]
            for j, target_orbit in enumerate(orbits):
                reduced[component_index, i, j] = sum(
                    int(source[component_index, row, index[state]])
                    for state in target_orbit)
    return reduced


D160 = reduce_source(D16)
R30 = reduce_source(R3)
check(np.array_equal(D160, D160.transpose(0, 2, 1)) and
      np.array_equal(R30, R30.transpose(0, 2, 1)),
      "exact orbit-sum reduction preserves source Hermiticity")
check(np.all(D160 % 4 == 0),
      "direct quarter-lattice and ring third-lattice give denominator-twelve k0 blocks")


def source_basis(rho):
    coordinates = []
    for component_index in range(6):
        coordinates.append(kmat([
            [F(rho*int(D160[component_index, i, j]), 16) +
             F(int(R30[component_index, i, j]), 3)
             for j in range(6)] for i in range(6)
        ]))
    A = kscale(IS3, kadd(kadd(coordinates[0], coordinates[1]), coordinates[2]))
    E1 = kscale(IS2, ksub(coordinates[0], coordinates[1]))
    E2 = kscale(IS6, ksub(kadd(coordinates[0], coordinates[1]),
                           kscale(2, coordinates[2])))
    return [A, E1, E2,
            kscale(IS2, coordinates[3]),
            kscale(IS2, coordinates[4]),
            kscale(IS2, coordinates[5])]


H02 = kmul(H0K, H0K)
P0 = kscale(F(1, 8), kmul(
    ksub(kadd(H02, kscale(4, H0K)), kscale(4, I6)),
    ksub(H0K, kscale(2, I6))))
P2 = kscale(F(1, 16), kmul(
    ksub(kadd(H02, kscale(4, H0K)), kscale(4, I6)), H0K))
check(kmatrix_equal(kmul(P0, P0), P0) and
      kmatrix_equal(kmul(P2, P2), P2) and
      kiszero_matrix(kmul(P0, P2)),
      "rational polynomial P0/P2 are exact orthogonal projectors")
check(ktrace(P0) == K(2) and ktrace(P2) == K(2),
      "both degenerate responding energy projectors have rank two")


def exact_residues(rho):
    q = source_basis(rho)
    vectors = [kmatvec(operator, g) for operator in q]
    residue_zero = [[kdot(vectors[a], kmatvec(P0, vectors[b]))
                     for b in range(6)] for a in range(6)]
    residue_two = [[kdot(vectors[a], kmatvec(P2, vectors[b]))
                    for b in range(6)] for a in range(6)]
    r_one = [ZERO, rho*IS2, -rho*S6/2, -3*IS2, -3*IS2, ZERO]
    r_two = [ZERO, ZERO, ZERO, 3*IS2, -3*IS2, ZERO]
    return q, vectors, residue_zero, residue_two, r_one, r_two


# Every residue entry is a polynomial of degree at most two in rho because the
# source is affine-linear.  Three exact interpolation points prove FW16 for
# arbitrary rho, rather than merely sampling the physical scale.
interpolation_packets = [exact_residues(value) for value in (0, 1, 2)]
for interpolation_rho, packet in zip((0, 1, 2), interpolation_packets):
    _, _, residue_zero, residue_two, r_one, r_two = packet
    check(kmatrix_equal(residue_zero, kouter(r_one, r_one)) and
          kmatrix_equal(residue_two, kouter(r_two, r_two)),
          f"exact residue polynomial matches at rho={interpolation_rho}")
check(all(kiszero_matrix(kadd(ksub(q2, kscale(2, q1)), q0))
          for q0, q1, q2 in zip(interpolation_packets[0][0],
                                interpolation_packets[1][0],
                                interpolation_packets[2][0])),
      "source blocks are affine in rho, so three-point residue interpolation is exact")

for rho in (1, 2):
    Q, source_vectors, residue0, residue2, r1, r2 = exact_residues(rho)
    check(kmatrix_equal(Q[0], kscale(IS3, kadd(kscale(60*rho, I6),
                                                kscale(11, H0K)))),
          f"exact k0 A1 identity/H relation holds at rho={rho}")
    Ec = kadd(kscale(S3/2, Q[1]), kscale(F(1, 2), Q[2]))
    Ea = kadd(kscale(F(-1, 2), Q[1]), kscale(S3/2, Q[2]))
    check(kmatrix_equal(Ec, kscale(16*S6*rho, I6)),
          f"exact k0 Ec identity relation holds at rho={rho}")
    check(kiszero_matrix(ksub(kmul(H0K, Q[5]), kmul(Q[5], H0K))),
          f"exact k0 Tyz commutation holds at rho={rho}")
    check(kmatvec(Q[5], g) == [(3*S2-6)*x for x in g],
          f"Tyz ground eigenvalue is 3sqrt2-6 at rho={rho}")

    expected_r1 = kouter(r1, r1)
    expected_r2 = kouter(r2, r2)
    check(kmatrix_equal(residue0, expected_r1),
          f"degeneracy-safe P0 residue is exact rank one at rho={rho}")
    check(kmatrix_equal(residue2, expected_r2),
          f"degeneracy-safe P2 residue is exact rank one at rho={rho}")
    check(kdot(r1, r1) == K(2*rho*rho+9) and
          kdot(r2, r2) == K(9) and kdot(r1, r2) == ZERO,
          f"exact residue norms and orthogonality hold at rho={rho}")

    for vector in source_vectors:
        centered = [x-y*kdot(g, vector) for x, y in zip(vector, g)]
        remainder = [x-y-z for x, y, z in zip(
            centered, kmatvec(P0, vector), kmatvec(P2, vector))]
        check(all(value.is_zero() for value in remainder),
              f"one centered source has no spectral support outside P0+P2 at rho={rho}")

    delta1 = 2+2*S2
    delta2 = 4+2*S2
    moment1 = kscale(-2, kadd(kscale(delta1, expected_r1),
                              kscale(delta2, expected_r2)))
    direct_m0 = kzeros(6, 6)
    direct_m1 = kzeros(6, 6)
    for a in range(6):
        ad = ksub(kmul(H0K, Q[a]), kmul(Q[a], H0K))
        for b in range(6):
            comm0 = ksub(kmul(Q[a], Q[b]), kmul(Q[b], Q[a]))
            comm1 = ksub(kmul(ad, Q[b]), kmul(Q[b], ad))
            direct_m0[a][b] = kdot(g, kmatvec(comm0, g))
            direct_m1[a][b] = kdot(g, kmatvec(comm1, g))
    check(kiszero_matrix(direct_m0),
          f"direct exact equal-time commutator M0 vanishes at rho={rho}")
    check(kmatrix_equal(direct_m1, moment1),
          f"direct exact M1 equals minus twice the gap-weighted residues at rho={rho}")
    check(kdot(r1, r2) == ZERO,
          f"M1 and retarded spectral image have common exact rank two at rho={rho}")

    # Static Kubo eigenvalues follow without numerical diagonalization because
    # the two residue vectors are orthogonal.
    kappa1 = K(2*rho*rho+9)*(S2-1)
    kappa2 = K(F(9, 2))*(2-S2)
    check(kappa1.value() > 0 and kappa2.value() > 0,
          f"both exact static Kubo eigenvalues are positive at rho={rho}")
    check((S2-1)*delta1 == K(2) and
          ((2-S2)/2)*delta2 == K(2),
          f"static Kubo gap normalization identities hold at rho={rho}")

    Tplus = kscale(IS2, kadd(Q[3], Q[4]))
    dark = ksub(kscale(3, Ea), kscale(S2*rho, Tplus))
    dark_g = kmatvec(dark, g)
    expectation = kdot(g, dark_g)
    centered_dark_g = [x-expectation*y for x, y in zip(dark_g, g)]
    check(all(value.is_zero() for value in centered_dark_g),
          f"nonconserved combination is exactly ground-dark at rho={rho}")
    dark_comm = ksub(kmul(H0K, dark), kmul(dark, H0K))
    check(not kiszero_matrix(dark_comm),
          f"ground-dark combination is not conserved at rho={rho}")


# Documentary attack: the exact result must remain a witness-only k=0 screen.
theorem = (LANE / "THEOREM.md").read_text()
result = (LANE / "RESULT.md").read_text()
self_audit = (LANE / "SELF_AUDIT.md").read_text()
joined = " ".join((theorem+result+self_audit).split())
for phrase in (
    "Q_{\\rm FV-WITNESS}",
    "Q_{\\rm diag}^{(2)}+Q_{\\rm diag}^{(4)}+Q_{\\rm diag}^{(6)}",
    "not an upper bound",
    "finite-sector, homogeneous-source result",
    "not a Ward identity",
    "massless tensor",
    "gravity",
    "numerical value of `G`",
):
    check(phrase in joined, f"frozen scope retains: {phrase}")

for forbidden in (
    "the complete FV source has ground-state retarded rank two",
    "rank two excludes a six-channel CTP packet",
    "Tyz conservation is a Ward identity",
    "the finite poles prove a graviton",
    "FW derives gravity",
):
    check(forbidden not in joined, f"forbidden FW promotion absent: {forbidden}")

check("FV-WITNESS" in theorem and "Q_diag^(2,4,6)" in theorem and
      "complete fixed-order H6 source response" in theorem,
      "successor calculation is the omitted generated-diagonal completion")

print(f"SUMMARY {checks}/{checks} independent hostile-audit checks passed")
print("HIERARCHY FV_family=6 component_mod_I=5 adH=3 ground_retarded=2 M1=2")
print("POLES exact gaps 2+2sqrt2 and 4+2sqrt2; projector residues rank1+rank1")
print("VERDICT PASS for FV-WITNESS only; complete Qdiag/CTP/Ward/helicity/gravity open")
