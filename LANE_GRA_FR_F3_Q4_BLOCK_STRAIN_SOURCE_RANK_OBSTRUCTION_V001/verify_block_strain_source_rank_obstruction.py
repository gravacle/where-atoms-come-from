#!/usr/bin/env python3
"""Exact checks for the frozen additive q4 block-strain source obstruction.

The calculation uses rational arithmetic only.  It tests the microscopic
tetrahedral edge-dyad map, its S4 type and nullspace, additive multi-edge
closure, source-quotient response factorization, a representative exact
Feshbach Schur complement, induced and independently supplied contact terms,
the six A3 root-dyad contrast, and operator-level commutator moments.
"""

from __future__ import annotations

from fractions import Fraction as F
import hashlib
from itertools import permutations, product
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent.parent
LANE = Path(__file__).resolve().parent
THEOREM_TEXT = (LANE / "THEOREM.md").read_text(encoding="utf-8")
SUMMARY_ONLY = "--summary" in sys.argv[1:]

DEPENDENCIES = {
    "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md":
        "00eba581b90fb9f0b25e3fad1362b055049824897433ba06e356ab9b1f6c76ec",
    "LANE_GRA_FC_F3_Q4_CLIFFORD_COLLECTIVE_CONE_V001/THEOREM.md":
        "28b6319e3187337da8ebef2212b030ff6e5b9f8168d9844ae172d94f3e0641a6",
    "LANE_GRA_FE_F3_Q4_DIAMOND_ICE_CARRIER_JOIN_V001/THEOREM.md":
        "4cc63e3e5853b4250a2a5b78256d41b83b195cf527901819a62a04ef53f8d932",
    "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md":
        "5e68e4a8c62ad89cff309781a5cb54e071092e82594e60c04e7992414bc18dbe",
    "LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/THEOREM.md":
        "05f4a619a6f80aa40c48570ab4035ab874426502a31468a08f435e66610bd769",
    "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/THEOREM.md":
        "cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98",
    "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md":
        "78f0687c9f597c96d235224dad45c204d12d7e6c973c270e3305a901efc75b25",
    "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001/THEOREM.md":
        "44fc28edc9820d2b4ea67cef9f83beef60e53bfe6407b582ffbeccfe42f756c5",
    "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001/INDEPENDENT_AUDIT.md":
        "84d8c02c3e560198f6a9fae04f5ee81bc72c98354e02f1f58d506a8c3171c453",
    "LANE_GRA_FQ_F3_Q4_COLLECTIVE_METRIC_ORIGIN_SCREEN_V001/THEOREM.md":
        "07445c035ed4c5167a5a20280c4db69a5101eeb71831cdeb126b29702d04b69d",
    "LANE_GRA_FQ_F3_Q4_COLLECTIVE_METRIC_ORIGIN_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "91aa35170432684a47278e46ee2b9d56658a43acc8acbbb84480d047cdbe6dcf",
}


checks = 0


def check(condition: bool, message: str) -> None:
    global checks
    if not condition:
        raise AssertionError(message)
    checks += 1
    if not SUMMARY_ONLY:
        print(f"PASS {' '.join(message.split())}")


def sha256(path: Path) -> str:
    if path.is_symlink() or not path.is_file():
        raise AssertionError(f"missing, non-file, or symlinked dependency: {path}")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def zeros(rows: int, columns: int):
    return [[F(0) for _ in range(columns)] for _ in range(rows)]


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def add(left, right):
    return [[a + b for a, b in zip(row_l, row_r)]
            for row_l, row_r in zip(left, right)]


def scale(value, matrix):
    return [[value * item for item in row] for row in matrix]


def multiply(left, right):
    right_t = transpose(right)
    return [[sum((a * b for a, b in zip(row, column)), F(0))
             for column in right_t] for row in left]


def identity(size: int):
    return [[F(int(row == column)) for column in range(size)]
            for row in range(size)]


def matrix_rank(matrix) -> int:
    work = [list(map(F, row)) for row in matrix]
    if not work:
        return 0
    rows, columns = len(work), len(work[0])
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, rows)
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][column]
        work[pivot_row] = [item / value for item in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [a - factor * b
                         for a, b in zip(work[row], work[pivot_row])]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def inverse(matrix):
    size = len(matrix)
    work = [list(map(F, row)) + ident
            for row, ident in zip(matrix, identity(size))]
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if work[row][column]), None)
        if pivot is None:
            raise AssertionError("singular exact matrix")
        work[column], work[pivot] = work[pivot], work[column]
        value = work[column][column]
        work[column] = [item / value for item in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [a - factor * b
                         for a, b in zip(work[row], work[column])]
    return [row[size:] for row in work]


def submatrix(matrix, row_indices, column_indices):
    return [[matrix[row][column] for column in column_indices]
            for row in row_indices]


def commutator(left, right):
    return add(multiply(left, right), scale(F(-1), multiply(right, left)))


def mat_equal(left, right):
    return left == right


def mat_zero(matrix):
    return all(item == 0 for row in matrix for item in row)


def vector_matrix(signs):
    return [[F(signs[row] * signs[column], 3) for column in range(3)]
            for row in range(3)]


def symvec(matrix):
    return (matrix[0][0], matrix[1][1], matrix[2][2],
            matrix[0][1], matrix[0][2], matrix[1][2])


def frobenius(left, right):
    return sum((left[row][column] * right[row][column]
                for row in range(3) for column in range(3)), F(0))


def cycle_type(permutation):
    seen = set()
    lengths = []
    for start in range(4):
        if start in seen:
            continue
        length = 0
        value = start
        while value not in seen:
            seen.add(value)
            length += 1
            value = permutation[value]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


for relative, expected in DEPENDENCIES.items():
    check(sha256(ROOT / relative) == expected,
          f"dependency custody matches {relative}")


# The four q4 append vectors.  Their normalization is exact because all
# matrices are outer products and therefore rational.
SIGNS = (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)
DYADS = tuple(vector_matrix(signs) for signs in SIGNS)
DYAD_VECTORS = tuple(symvec(matrix) for matrix in DYADS)

check(all(frobenius(dyad, dyad) == 1 for dyad in DYADS),
      "each normalized tetrahedral edge dyad has unit Frobenius norm")
check(all(frobenius(DYADS[a], DYADS[b]) == F(1, 9)
          for a in range(4) for b in range(4) if a != b),
      "distinct tetrahedral edge dyads have exact Gram overlap one ninth")
gram = [[frobenius(left, right) for right in DYADS] for left in DYADS]
check(matrix_rank(gram) == 4, "tetrahedral edge-dyad Gram matrix has rank four")
uniform = [F(1)] * 4
check([sum(row) for row in gram] == [F(4, 3)] * 4,
      "uniform edge-dyad direction is the A1 Gram eigenvector")
for contrast in ((1, -1, 0, 0), (1, 0, -1, 0), (1, 0, 0, -1)):
    image = [sum((row[i] * contrast[i] for i in range(4)), F(0))
             for row in gram]
    check(image == [F(8, 9) * value for value in contrast],
          "each independent contrast is a T2 Gram eigenvector")

# Character of the four-dyad permutation module minus its uniform A1 line.
expected_t2 = {
    (1, 1, 1, 1): 3,
    (2, 1, 1): 1,
    (2, 2): -1,
    (3, 1): 0,
    (4,): -1,
}
observed = {}
for permutation in permutations(range(4)):
    kind = cycle_type(permutation)
    character = sum(int(permutation[index] == index) for index in range(4)) - 1
    observed.setdefault(kind, set()).add(character)
check(all(values == {expected_t2[kind]} for kind, values in observed.items()),
      "contrast character is exactly the frozen S4 T2 character")

# A symmetric source is encoded as (xx,yy,zz,xy,xz,yz).  Frobenius
# contraction with a dyad doubles its off-diagonal coordinates.
CONTRACTION = [
    [F(1), F(1), F(1), F(2), F(2), F(2)],
    [F(1), F(1), F(1), F(-2), F(-2), F(2)],
    [F(1), F(1), F(1), F(-2), F(2), F(-2)],
    [F(1), F(1), F(1), F(2), F(-2), F(-2)],
]
E_NULL = (
    [F(1), F(-1), F(0), F(0), F(0), F(0)],
    [F(1), F(1), F(-2), F(0), F(0), F(0)],
)

check(matrix_rank(CONTRACTION) == 4,
      "microscopic strain-to-edge contraction map has rank four")
for vector in E_NULL:
    check([sum((row[i] * vector[i] for i in range(6)), F(0))
           for row in CONTRACTION] == [F(0)] * 4,
          "one diagonal-traceless E generator is an exact source null")
check(matrix_rank(CONTRACTION + [list(vector) for vector in E_NULL]) == 6,
      "the two displayed E nulls complete the six-dimensional source space")

sum_dyads = zeros(3, 3)
for dyad in DYADS:
    sum_dyads = add(sum_dyads, dyad)
check(sum_dyads == scale(F(4, 3), identity(3)),
      "four incident additive dyads sum to the isotropic A1 tensor")

# Exhaust every nonnegative label multiplicity vector of total degree <= 8.
additive_vectors = []
additive_null_checks = [True, True]
for counts in product(range(9), repeat=4):
    if sum(counts) > 8:
        continue
    vector = [sum((F(counts[a]) * DYAD_VECTORS[a][column]
                   for a in range(4)), F(0)) for column in range(6)]
    additive_vectors.append(vector)
    for null_index, null in enumerate(E_NULL):
        # The nulls are diagonal, so raw symmetric coordinates have the
        # ordinary diagonal contraction here.
        additive_null_checks[null_index] &= (
            sum((vector[i] * null[i] for i in range(3)), F(0)) == 0
        )
for null_index, condition in enumerate(additive_null_checks, start=1):
    check(condition,
          f"every additive weight through degree eight annihilates E null {null_index}")
check(matrix_rank(additive_vectors) == 4,
      "all additive FQ17a weights through eight occurrences still span rank four")

# Six blocked A3 sibling-root dyads are a different rank-six source.  They
# contain cross terms absent from D_a + D_b.
ROOT_DYADS = []
ADDITIVE_PAIRS = []
for a in range(4):
    for b in range(a + 1, 4):
        difference = [SIGNS[b][i] - SIGNS[a][i] for i in range(3)]
        root_matrix = [[F(difference[i] * difference[j], 3)
                        for j in range(3)] for i in range(3)]
        ROOT_DYADS.append(list(symvec(root_matrix)))
        ADDITIVE_PAIRS.append(list(symvec(add(DYADS[a], DYADS[b]))))
check(matrix_rank(ROOT_DYADS) == 6,
      "six A3 sibling-root dyads span the full symmetric tensor space")
check(matrix_rank(ADDITIVE_PAIRS) == 4,
      "six additive two-edge weights remain in the microscopic rank-four span")
check(any(root != additive for root, additive in zip(ROOT_DYADS, ADDITIVE_PAIRS)),
      "root dyads differ from additive edge dyads by nonzero cross terms")

# Any connected or CTP derivative generated from x_a=(Cj)_a factors through
# the four-dimensional quotient.  A generic nonsingular internal response K
# therefore pulls back to C^T K C with rank four and the two exact nulls.
KERNEL = [
    [F(3), F(1), F(0), F(1)],
    [F(1), F(4), F(1), F(0)],
    [F(0), F(1), F(5), F(1)],
    [F(1), F(0), F(1), F(6)],
]
check(matrix_rank(KERNEL) == 4, "internal four-invariant test kernel is nonsingular")
PULLED_RESPONSE = multiply(transpose(CONTRACTION),
                           multiply(KERNEL, CONTRACTION))
check(matrix_rank(PULLED_RESPONSE) == 4,
      "generic pulled-back two-point or CTP Hessian has rank four")
for null in E_NULL:
    column = [[value] for value in null]
    check(mat_zero(multiply(PULLED_RESPONSE, column)),
          "pulled-back response annihilates one E source direction")

# A representative exact Feshbach Schur complement is invariant under both
# E shifts because the complete microscopic matrix depends only on Cj.
BASE = [
    [F(2), F(1), F(1), F(0)],
    [F(1), F(3), F(0), F(1)],
    [F(1), F(0), F(5), F(1)],
    [F(0), F(1), F(1), F(7)],
]
VARIATIONS = (
    [[F(1), F(0), F(1), F(0)], [F(0), F(0), F(0), F(0)],
     [F(1), F(0), F(0), F(0)], [F(0), F(0), F(0), F(0)]],
    [[F(0), F(1), F(0), F(0)], [F(1), F(0), F(0), F(1)],
     [F(0), F(0), F(0), F(0)], [F(0), F(1), F(0), F(0)]],
    [[F(0), F(0), F(0), F(1)], [F(0), F(1), F(1), F(0)],
     [F(0), F(1), F(0), F(0)], [F(1), F(0), F(0), F(0)]],
    [[F(1), F(0), F(0), F(0)], [F(0), F(1), F(0), F(0)],
     [F(0), F(0), F(1), F(1)], [F(0), F(0), F(1), F(0)]],
)


def microscopic(source):
    invariants = [sum((CONTRACTION[a][i] * source[i]
                       for i in range(6)), F(0)) for a in range(4)]
    result = [row[:] for row in BASE]
    for coefficient, variation in zip(invariants, VARIATIONS):
        result = add(result, scale(coefficient, variation))
    return result


def feshbach(matrix, energy=F(13)):
    p = (0, 1)
    q = (2, 3)
    h_pp = submatrix(matrix, p, p)
    h_pq = submatrix(matrix, p, q)
    h_qp = submatrix(matrix, q, p)
    h_qq = submatrix(matrix, q, q)
    resolvent = inverse(add(scale(energy, identity(2)), scale(F(-1), h_qq)))
    return add(h_pp, multiply(h_pq, multiply(resolvent, h_qp)))


source = [F(2), F(-1), F(3), F(1), F(-2), F(1)]
reference_micro = microscopic(source)
reference_eff = feshbach(reference_micro)
for null in E_NULL:
    shifted = [value + shift for value, shift in zip(source, null)]
    check(mat_equal(microscopic(shifted), reference_micro),
          "complete microscopic source matrix is invariant under one E shift")
    check(mat_equal(feshbach(microscopic(shifted)), reference_eff),
          "exact Feshbach resolvent and fold preserve one E shift")

# Contacts induced from the same four invariants retain the nulls.  An
# independently supplied quadratic E contact can make a seagull Hessian, but
# its first derivative at j=0 is exactly zero and it has no spectral operator.
DERIVED_CONTACT_HESSIAN = PULLED_RESPONSE
for null in E_NULL:
    check(mat_zero(multiply(DERIVED_CONTACT_HESSIAN, [[x] for x in null])),
          "same-source induced contact retains one E null")

explicit_e_seagull = add(
    multiply([[x] for x in E_NULL[0]], [E_NULL[0]]),
    multiply([[x] for x in E_NULL[1]], [E_NULL[1]]),
)
check(matrix_rank(explicit_e_seagull) == 2,
      "an independently inserted quadratic E seagull can have Hessian rank two")
zero_source = [F(0)] * 6
explicit_contact_gradient = [
    sum((explicit_e_seagull[row][column] * zero_source[column]
         for column in range(6)), F(0)) for row in range(6)
]
check(explicit_contact_gradient == [F(0)] * 6,
      "an O(j^2) E seagull has no linear source-conjugate operator at j=0")

# Operator-valued Q_mu factors through C^T.  Consequently both E-contracted
# Q operators vanish before a state is selected, and every nested-commutator
# moment with an E leg vanishes.
O = (
    [[F(1), F(1), F(0)], [F(1), F(0), F(0)], [F(0), F(0), F(-1)]],
    [[F(0), F(1), F(1)], [F(1), F(2), F(0)], [F(1), F(0), F(1)]],
    [[F(2), F(0), F(1)], [F(0), F(-1), F(1)], [F(1), F(1), F(0)]],
    [[F(1), F(0), F(0)], [F(0), F(3), F(1)], [F(0), F(1), F(-2)]],
)
Q = []
for mu in range(6):
    operator = zeros(3, 3)
    for a in range(4):
        operator = add(operator, scale(CONTRACTION[a][mu], O[a]))
    Q.append(operator)

for null in E_NULL:
    contracted = zeros(3, 3)
    for coefficient, operator in zip(null, Q):
        contracted = add(contracted, scale(coefficient, operator))
    check(mat_zero(contracted), "E-contracted effective source operator is identically zero")

H_TEST = [[F(1), F(1), F(0)], [F(1), F(4), F(1)], [F(0), F(1), F(7)]]
for order in range(7):
    evolved = Q
    for _ in range(order):
        evolved = [commutator(H_TEST, operator) for operator in evolved]
    nested_nulls = []
    moment_nulls = []
    for null in E_NULL:
        left = zeros(3, 3)
        for coefficient, operator in zip(null, evolved):
            left = add(left, scale(coefficient, operator))
        nested_nulls.append(mat_zero(left))
        moment_nulls.extend(mat_zero(commutator(left, probe)) for probe in Q)
    check(all(nested_nulls),
          f"order-{order} nested source operators vanish on both E directions")
    check(all(moment_nulls),
          f"order-{order} commutator moments vanish with either E leg")

check(matrix_rank(CONTRACTION) < 6,
      "frozen additive source fails the six-off-shell-direction prerequisite")
check(matrix_rank(ROOT_DYADS) == 6 and matrix_rank(CONTRACTION) == 4,
      "a rank-six root-edge query is mathematically distinct from the frozen source")

check("The complete BS20a contact `R` need not obey (FR11)" in THEOREM_TEXT,
      "theorem separates general O(j^2) contacts from invariant contacts")
check("D_EH_{\\rm eff}(E,j)\\big|_{j=0}=0" in THEOREM_TEXT,
      "theorem promotes only the unconditional source-off linear Feshbach null")
check("For general BS20a contacts, (FR17) need not hold" in THEOREM_TEXT,
      "theorem does not overpromote finite-shift invariance to all contacts")
check("may be a legitimate BS20a contact" in THEOREM_TEXT,
      "theorem does not misclassify every E seagull as a new query")
check("does not set every `E`-leg CTP derivative to" in THEOREM_TEXT,
      "theorem preserves the general contact Hessian caveat")
check("Every connected CTP derivative with at least one `E` leg vanishes" not in THEOREM_TEXT,
      "original all-CTP-derivative overclaim is absent")
check("this packet does not invent them" in THEOREM_TEXT,
      "theorem does not silently assign non-edge linear weights")
check("not yet a rank theorem for the complete BS20 source" in THEOREM_TEXT,
      "theorem preserves the complete-source ceiling")

print(f"SUMMARY {checks}/{checks} exact checks passed")
print("DISPOSITION ADDITIVE_EDGE_SUPPORTED_Q4_SOURCE_RANK4_A1_PLUS_T2__E_NULL2__COMPLETE_NONEDGE_LINEAR_WEIGHT_CENSUS_OPEN")
