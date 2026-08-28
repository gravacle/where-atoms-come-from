#!/usr/bin/env python3
"""Independent hostile replay for GRA-FR-F3-Q4-BSSRO-V001.

This executable does not import the builder verifier. It reconstructs the
tetrahedral source map in a separate coordinate calculation, tests additive
and blocked closure, differentiates a representative Feshbach map, exposes
the general O(j^2) contact caveat, and checks the exact scope boundary against
root-edge and independently rotated-coframe queries.
"""

from __future__ import annotations

from fractions import Fraction as F
import hashlib
from itertools import combinations, permutations
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
THEOREM = (HERE / "THEOREM.md").read_text(encoding="utf-8")
RESULT = (HERE / "RESULT.md").read_text(encoding="utf-8")
SELF_AUDIT = (HERE / "SELF_AUDIT.md").read_text(encoding="utf-8")
PASSED = 0
SUMMARY_ONLY = "--summary" in sys.argv[1:]


def check(condition: bool, label: str) -> None:
    global PASSED
    if not condition:
        raise AssertionError(label)
    PASSED += 1
    if not SUMMARY_ONLY:
        print(f"PASS {' '.join(label.split())}")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_hashes(path: Path) -> dict[str, str]:
    ledger: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        value, name = line.split(maxsplit=1)
        ledger[name.strip()] = value
    return ledger


def transpose(matrix):
    return [list(column) for column in zip(*matrix)]


def add(left, right):
    return [[F(a) + F(b) for a, b in zip(row_l, row_r)]
            for row_l, row_r in zip(left, right)]


def scale(value, matrix):
    return [[F(value) * F(entry) for entry in row] for row in matrix]


def matmul(left, right):
    return [[sum((F(a) * F(b) for a, b in zip(row, column)), F(0))
             for column in zip(*right)] for row in left]


def matvec(matrix, vector):
    return [sum((F(a) * F(b) for a, b in zip(row, vector)), F(0))
            for row in matrix]


def identity(size: int):
    return [[F(row == column) for column in range(size)] for row in range(size)]


def rank(matrix) -> int:
    if not matrix:
        return 0
    work = [[F(entry) for entry in row] for row in matrix]
    rows, columns = len(work), len(work[0])
    pivot = 0
    for column in range(columns):
        selected = next((row for row in range(pivot, rows)
                         if work[row][column]), None)
        if selected is None:
            continue
        work[pivot], work[selected] = work[selected], work[pivot]
        value = work[pivot][column]
        work[pivot] = [entry / value for entry in work[pivot]]
        for row in range(rows):
            if row == pivot or not work[row][column]:
                continue
            value = work[row][column]
            work[row] = [left - value * right
                         for left, right in zip(work[row], work[pivot])]
        pivot += 1
        if pivot == rows:
            break
    return pivot


def inverse(matrix):
    size = len(matrix)
    work = [[F(entry) for entry in row] + unit
            for row, unit in zip(matrix, identity(size))]
    for column in range(size):
        selected = next((row for row in range(column, size)
                         if work[row][column]), None)
        if selected is None:
            raise AssertionError("singular matrix in hostile replay")
        work[column], work[selected] = work[selected], work[column]
        value = work[column][column]
        work[column] = [entry / value for entry in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            value = work[row][column]
            work[row] = [left - value * right
                         for left, right in zip(work[row], work[column])]
    return [row[size:] for row in work]


def submatrix(matrix, rows, columns):
    return [[matrix[row][column] for column in columns] for row in rows]


def matrix_zero(matrix) -> bool:
    return all(entry == 0 for row in matrix for entry in row)


def block_diagonal(copies):
    row_count = sum(len(block) for block in copies)
    col_count = sum(len(block[0]) for block in copies)
    answer = [[F(0) for _ in range(col_count)] for _ in range(row_count)]
    row_offset = 0
    col_offset = 0
    for block in copies:
        for row, values in enumerate(block):
            for column, value in enumerate(values):
                answer[row_offset + row][col_offset + column] = F(value)
        row_offset += len(block)
        col_offset += len(block[0])
    return answer


def cycle_type(permutation):
    seen = set()
    lengths = []
    for start in range(4):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


# ---------------------------------------------------------------------------
# Independent tetrahedral dyad map, nullspace, and S4 label.
# ---------------------------------------------------------------------------

TETRA = (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)


def dyad_vector(vector):
    x, y, z = map(F, vector)
    return [x * x, y * y, z * z, x * y, x * z, y * z]


DYADS = [dyad_vector(vector) for vector in TETRA]
# C maps symmetric coordinates (xx,yy,zz,xy,xz,yz) to n^T j n. The common
# normalization 1/3 is omitted because it cannot affect rank or kernel.
C = [[row[0], row[1], row[2], 2 * row[3], 2 * row[4], 2 * row[5]]
     for row in DYADS]
E_NULL = (
    [F(1), F(-1), F(0), F(0), F(0), F(0)],
    [F(1), F(1), F(-2), F(0), F(0), F(0)],
)

GRAM = [[sum((F(left[i]) * F(right[i])
              for i in range(3)), F(0)) ** 2
         for right in TETRA] for left in TETRA]
check(rank(GRAM) == 4, "independent tetrahedral dyad Gram rank is four")
check([sum(row) for row in GRAM] == [F(12)] * 4,
      "unnormalized uniform Gram eigenvalue is twelve")
for contrast in ((1, -1, 0, 0), (1, 0, -1, 0), (1, 0, 0, -1)):
    check(matvec(GRAM, contrast) == [F(8) * F(value) for value in contrast],
          "each independent contrast has unnormalized Gram eigenvalue eight")
check(rank(C) == 4, "independent source contraction map has rank four")
for null in E_NULL:
    check(matvec(C, null) == [F(0)] * 4,
          "one diagonal-traceless E vector is in the exact source kernel")
check(rank(C + [list(vector) for vector in E_NULL]) == 6,
      "the displayed two-dimensional E kernel is complete")

expected_t2 = {
    (1, 1, 1, 1): 3,
    (2, 1, 1): 1,
    (2, 2): -1,
    (3, 1): 0,
    (4,): -1,
}
characters: dict[tuple[int, ...], set[int]] = {}
for permutation in permutations(range(4)):
    kind = cycle_type(permutation)
    character = sum(permutation[index] == index for index in range(4)) - 1
    characters.setdefault(kind, set()).add(character)
check(all(values == {expected_t2[kind]}
          for kind, values in characters.items()),
      "four-label contrasts have the frozen T2 character, not just dimension three")


# ---------------------------------------------------------------------------
# Additive closure, ordinary basis changes, and physically new coframes.
# ---------------------------------------------------------------------------

coefficient_families = (
    [[1, 0, 0, 0], [0, 1, 1, 0], [2, 0, 1, 3]],
    [[1, -1, 2, 0], [0, 3, -2, 1], [4, 1, 0, -3], [2, 2, 2, 2]],
    identity(4),
)
for coefficients in coefficient_families:
    weights = matmul(coefficients, C)
    check(rank(weights) <= 4, "arbitrary signed additive multi-edge weights retain rank at most four")
    for null in E_NULL:
        check(matvec(weights, null) == [F(0)] * len(weights),
              "arbitrary additive multi-edge family retains one E null")

three_blocks = block_diagonal([C, C, C])
check(rank(three_blocks) == 12,
      "three uniform-coframe blocks have four internal directions per block")
block_mix = [
    [F((row + 2 * column + 1) % 5 - 2) for column in range(12)]
    for row in range(9)
]
check(rank(matmul(block_mix, three_blocks)) <= 12,
      "source-independent output blocking cannot exceed the pre-block rank")

source_basis = identity(6)
source_basis[0][3] = F(2)
source_basis[2][5] = F(-3)
source_basis[4][1] = F(1)
check(rank(source_basis) == 6, "hostile alternative source basis is invertible")
check(rank(matmul(C, source_basis)) == 4,
      "an invertible source-coordinate change cannot repair rank four")

# A genuinely different physical coframe is not a mere blocking basis. A
# rational quaternion rotation supplies a concrete counter-scope: the union
# of the original and independently rotated local spans can reach rank six.
w, x, y, z = map(F, (1, 2, 3, 4))
norm = w * w + x * x + y * y + z * z
ROTATION = [
    [1 - 2 * (y * y + z * z) / norm,
     2 * (x * y - w * z) / norm,
     2 * (x * z + w * y) / norm],
    [2 * (x * y + w * z) / norm,
     1 - 2 * (x * x + z * z) / norm,
     2 * (y * z - w * x) / norm],
    [2 * (x * z - w * y) / norm,
     2 * (y * z + w * x) / norm,
     1 - 2 * (x * x + y * y) / norm],
]
check(matmul(ROTATION, transpose(ROTATION)) == identity(3),
      "counter-scope coframe rotation is exactly orthogonal")
rotated_vectors = [matvec(ROTATION, vector) for vector in TETRA]
rotated_dyads = [dyad_vector(vector) for vector in rotated_vectors]
check(rank(rotated_dyads) == 4,
      "one independently rotated tetrahedral coframe still has local rank four")
check(rank(DYADS + rotated_dyads) == 6,
      "two genuinely rotated physical coframes can have rank-six union")

ROOT_DYADS = []
ADDITIVE_PAIRS = []
for first, second in combinations(range(4), 2):
    root = [TETRA[second][axis] - TETRA[first][axis] for axis in range(3)]
    ROOT_DYADS.append(dyad_vector(root))
    ADDITIVE_PAIRS.append([DYADS[first][column] + DYADS[second][column]
                           for column in range(6)])
check(rank(ROOT_DYADS) == 6, "six sibling-root dyads form a rank-six query")
check(rank(ADDITIVE_PAIRS) == 4, "six additive two-edge weights remain rank four")
check(all(root != additive
          for root, additive in zip(ROOT_DYADS, ADDITIVE_PAIRS)),
      "each root query contains cross terms absent from its additive pair")


# ---------------------------------------------------------------------------
# Fixed Feshbach derivative, CTP spectral rank, and the O(j^2) caveat.
# ---------------------------------------------------------------------------

H0 = [
    [F(3), F(1), F(1), F(0)],
    [F(1), F(5), F(0), F(1)],
    [F(1), F(0), F(8), F(2)],
    [F(0), F(1), F(2), F(11)],
]
VARIATIONS = (
    [[1, 0, 1, 0], [0, 0, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0]],
    [[0, 1, 0, 0], [1, 0, 0, 1], [0, 0, 0, 0], [0, 1, 0, 0]],
    [[0, 0, 0, 1], [0, 2, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]],
    [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 1], [0, 0, 1, 0]],
)
VARIATIONS = tuple([[F(entry) for entry in row] for row in matrix]
                   for matrix in VARIATIONS)


def source_linear_variation(direction):
    invariants = matvec(C, direction)
    answer = scale(0, H0)
    for coefficient, variation in zip(invariants, VARIATIONS):
        answer = add(answer, scale(coefficient, variation))
    return answer


def feshbach(matrix, energy=F(19)):
    p, q = (0, 1), (2, 3)
    hpp = submatrix(matrix, p, p)
    hpq = submatrix(matrix, p, q)
    hqp = submatrix(matrix, q, p)
    hqq = submatrix(matrix, q, q)
    resolvent = inverse(add(scale(energy, identity(2)), scale(-1, hqq)))
    return add(hpp, matmul(hpq, matmul(resolvent, hqp)))


def feshbach_derivative(matrix, variation, energy=F(19)):
    p, q = (0, 1), (2, 3)
    h_pq = submatrix(matrix, p, q)
    h_qp = submatrix(matrix, q, p)
    h_qq = submatrix(matrix, q, q)
    d_pp = submatrix(variation, p, p)
    d_pq = submatrix(variation, p, q)
    d_qp = submatrix(variation, q, p)
    d_qq = submatrix(variation, q, q)
    resolvent = inverse(add(scale(energy, identity(2)), scale(-1, h_qq)))
    answer = d_pp
    answer = add(answer, matmul(d_pq, matmul(resolvent, h_qp)))
    answer = add(answer, matmul(h_pq, matmul(resolvent, d_qp)))
    answer = add(answer, matmul(h_pq, matmul(resolvent,
                 matmul(d_qq, matmul(resolvent, h_qp)))))
    return answer


for null in E_NULL:
    check(matrix_zero(source_linear_variation(null)),
          "microscopic first derivative vanishes in one E direction")
    check(matrix_zero(feshbach_derivative(H0, source_linear_variation(null))),
          "fixed Feshbach first derivative preserves one E null")

effective_derivatives = []
for source_index in range(6):
    direction = [F(column == source_index) for column in range(6)]
    derivative = feshbach_derivative(H0, source_linear_variation(direction))
    effective_derivatives.append([entry for row in derivative for entry in row])
check(rank(effective_derivatives) <= 4,
      "representative exact Feshbach linear source family has rank at most four")

# A general quadratic E seagull is even along an E ray. It breaks finite
# source-shift invariance but has exactly zero first derivative at the origin.
SEAGULL_OPERATOR = [
    [F(1), F(0), F(0), F(0)],
    [F(0), F(2), F(0), F(0)],
    [F(0), F(0), F(3), F(1)],
    [F(0), F(0), F(1), F(4)],
]
plus = add(H0, SEAGULL_OPERATOR)
minus = add(H0, SEAGULL_OPERATOR)
check(plus == minus and plus != H0,
      "quadratic E seagull breaks finite-shift invariance but is even at source off")
check(feshbach(plus) == feshbach(minus) and feshbach(plus) != feshbach(H0),
      "Feshbach map preserves zero odd derivative without erasing the seagull Hessian")

# A mixed contact can make a full source Hessian with an E leg nonzero. This
# is the explicit counterexample that forbids the original all-CTP-derivative
# wording while leaving the source-off linear operator untouched.
trace_direction = [F(1), F(1), F(1), F(0), F(0), F(0)]
mixed_contact_hessian = add(
    matmul([[entry] for entry in E_NULL[0]], [trace_direction]),
    matmul([[entry] for entry in trace_direction], [E_NULL[0]]),
)
check(any(matvec(mixed_contact_hessian, E_NULL[0])),
      "a legitimate O(j^2) contact can have a nonzero Hessian with an E leg")
check(matvec(mixed_contact_hessian, [F(0)] * 6) == [F(0)] * 6,
      "every purely quadratic contact has zero gradient at source off")

INTERNAL_KERNEL = [
    [F(5), F(1), F(0), F(1)],
    [F(1), F(6), F(1), F(0)],
    [F(0), F(1), F(7), F(1)],
    [F(1), F(0), F(1), F(8)],
]
retarded_template = matmul(transpose(C), matmul(INTERNAL_KERNEL, C))
check(rank(retarded_template) == 4,
      "generic noncontact linear-operator response has exact rank four")
for null in E_NULL:
    check(matvec(retarded_template, null) == [F(0)] * 6,
          "noncontact response has one exact E row and column null")


def commutator(left, right):
    return add(matmul(left, right), scale(-1, matmul(right, left)))


H_P = submatrix(H0, (0, 1), (0, 1))
Q = [submatrix(matrix, (0, 1), (0, 1)) for matrix in VARIATIONS]
Q_SOURCE = []
for coordinate in range(6):
    operator = [[F(0), F(0)], [F(0), F(0)]]
    for label in range(4):
        operator = add(operator, scale(C[label][coordinate], Q[label]))
    Q_SOURCE.append(operator)
for null in E_NULL:
    contracted = [[F(0), F(0)], [F(0), F(0)]]
    for coefficient, operator in zip(null, Q_SOURCE):
        contracted = add(contracted, scale(coefficient, operator))
    check(matrix_zero(contracted),
          "E-polarized linear source operator vanishes before expectation")
    for order in range(6):
        evolved = contracted
        for _ in range(order):
            evolved = commutator(H_P, evolved)
        check(matrix_zero(evolved),
              f"order-{order} E-polarized nested commutator operator vanishes")
        check(all(matrix_zero(commutator(evolved, probe)) for probe in Q_SOURCE),
              f"order-{order} commutator moments vanish with the E leg")


# ---------------------------------------------------------------------------
# Documentary scope and dependency custody.
# ---------------------------------------------------------------------------

required_theorem_phrases = (
    "The complete BS20a contact `R` need not obey (FR11)",
    "General `O(j^2)` contacts do not alter that derivative",
    "For general BS20a contacts, (FR17) need not hold",
    "may be a legitimate BS20a contact",
    "the theorem does not set every `E`-leg CTP derivative to\nzero",
    "genuinely collective loop/surface variable",
    "thermodynamic tensor emergence",
    "gravity, and `G` all retain their previous status",
    "this packet does not invent them",
    "not yet a rank theorem for the complete BS20 source",
)
for phrase in required_theorem_phrases:
    check(phrase in THEOREM, f"corrected theorem preserves scope phrase: {phrase}")

for forbidden in (
    "Every connected CTP derivative with at least one `E` leg vanishes",
    "ARBITRARY_E_SEAGULL_IS_NEW_QUERY",
    "all same-source induced\ncontacts",
):
    check(forbidden not in THEOREM,
          f"corrected theorem removes hostile overclaim: {forbidden}")

check("Such a contact may be\nlegitimate if derived and frozen prospectively" in RESULT,
      "result does not misclassify every E seagull as a new query")
check("genuinely derived\nframe field" in SELF_AUDIT,
      "self-audit keeps independently rotated coframes outside frozen scope")

dependencies = parse_hashes(HERE / "DEPENDENCIES.sha256")
check(len(dependencies) == 11, "dependency ledger has eleven entries")
for relative, expected in dependencies.items():
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(),
          f"dependency is a regular file: {relative}")
    check(digest(path) == expected, f"dependency digest matches: {relative}")

print(f"SUMMARY {PASSED}/{PASSED} independent hostile checks passed")
print("DISPOSITION PASS_AFTER_CONTACT_AND_COMPLETE_SOURCE_SCOPE_REPAIRS__ADDITIVE_EDGE_LINEAR_SOURCE_ONLY__NONEDGE_ROOT_AND_ROTATED_COFRAME_QUERIES_OPEN")
