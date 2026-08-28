#!/usr/bin/env python3
"""Exact finite replay for GRA-FN-F3-Q4-ITFSB-V001.

The verifier uses only integer and rational arithmetic for the ice geometry,
uniform moments, Hessian, representation eigenvalues, complement boundary,
generic broken-background rank, and symmetrized-mixture identity.  It does
not test state preparation, a continuum pole, physical metric solder, or
gravity.
"""

from __future__ import annotations

import hashlib
import itertools
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
THEOREM_BYTES = (HERE / "THEOREM.md").read_bytes()
THEOREM = THEOREM_BYTES.decode("utf-8")
PASSED = 0


DEPENDENCIES = {
    "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/THEOREM.md":
        "cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98",
    "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/INDEPENDENT_AUDIT.md":
        "c52eab9d701d1c6e82f1d7ec395841f4d2810e96cccbc3e2504760b6742e81e4",
    "LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/THEOREM.md":
        "98e2b3bc7a1c998d7839dc1a6b435cc1c8ed6d5a622ba45f63571be9ef646452",
    "LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/INDEPENDENT_AUDIT.md":
        "327bf6a4476c4c6382757dc156a96c6032233d34c25c1f7935e2582acf6c607a",
}


def check(condition: bool, label: str) -> None:
    global PASSED
    if not condition:
        raise AssertionError(label)
    PASSED += 1
    print(f"PASS {label}")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dot(x: tuple[F, ...], y: tuple[F, ...]) -> F:
    return sum((a * b for a, b in zip(x, y)), F(0))


def outer(x: tuple[F, ...], y: tuple[F, ...]) -> tuple[tuple[F, ...], ...]:
    return tuple(tuple(a * b for b in y) for a in x)


def madd(*matrices: tuple[tuple[F, ...], ...]) -> tuple[tuple[F, ...], ...]:
    return tuple(
        tuple(sum((matrix[i][j] for matrix in matrices), F(0)) for j in range(len(matrices[0][0])))
        for i in range(len(matrices[0]))
    )


def mscale(c: F, matrix: tuple[tuple[F, ...], ...]) -> tuple[tuple[F, ...], ...]:
    return tuple(tuple(c * value for value in row) for row in matrix)


def mmul(
    left: tuple[tuple[F, ...], ...],
    right: tuple[tuple[F, ...], ...],
) -> tuple[tuple[F, ...], ...]:
    return tuple(
        tuple(
            sum((left[i][k] * right[k][j] for k in range(len(right))), F(0))
            for j in range(len(right[0]))
        )
        for i in range(len(left))
    )


def transpose(matrix: tuple[tuple[F, ...], ...]) -> tuple[tuple[F, ...], ...]:
    return tuple(tuple(matrix[j][i] for j in range(len(matrix))) for i in range(len(matrix[0])))


def diag(values: tuple[F, ...]) -> tuple[tuple[F, ...], ...]:
    return tuple(
        tuple(values[i] if i == j else F(0) for j in range(len(values)))
        for i in range(len(values))
    )


I3 = diag((F(1), F(1), F(1)))


def flatten_sym(matrix: tuple[tuple[F, ...], ...]) -> tuple[F, ...]:
    return (
        matrix[0][0], matrix[1][1], matrix[2][2],
        matrix[0][1], matrix[0][2], matrix[1][2],
    )


def rank(columns: list[tuple[F, ...]]) -> int:
    if not columns:
        return 0
    matrix = [list(row) for row in zip(*columns)]
    rows = len(matrix)
    cols = len(matrix[0])
    pivot_row = 0
    for col in range(cols):
        pivot = next((row for row in range(pivot_row, rows) if matrix[row][col] != 0), None)
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][col]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        for row in range(rows):
            if row == pivot_row or matrix[row][col] == 0:
                continue
            factor = matrix[row][col]
            matrix[row] = [a - factor * b for a, b in zip(matrix[row], matrix[pivot_row])]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


for relative, expected in DEPENDENCIES.items():
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(), f"dependency exists: {relative}")
    check(digest(path) == expected, f"dependency frozen: {relative}")

for forbidden, label in (
    (bytes((0x0D,)), "carriage return"),
    (bytes((0x08,)), "backspace"),
    (bytes((0x0C,)), "form feed"),
):
    check(forbidden not in THEOREM_BYTES, f"theorem contains no {label} control byte")

for required in (
    "SECOND_JET_NOT_LINEAR_TANGENT",
    "T2_TERM_IS_VECTOR_MEAN_DYAD",
    "COMPLEMENT_PRESERVING_FISHER_TANGENT_E_ONLY",
    "PHYSICAL_METRIC_SOLDER_REMAINS_OPEN",
    "not `O(3)`-equivariant",
    "external classical record",
    "not an independently derived helicity-two",
    "six-parameter Fisher metric",
    "six propagating degrees of",
    "conditional on `MAXWELL-IR`",
):
    check(required in THEOREM, f"promotion ceiling present: {required}")


# The six two-in/two-out sign states and the three doubled orthonormal axes.
ice4 = tuple(
    tuple(F(value) for value in state)
    for state in itertools.product((-1, 1), repeat=4)
    if sum(state) == 0
)
axes4 = (
    (F(1, 2), F(1, 2), F(-1, 2), F(-1, 2)),
    (F(1, 2), F(-1, 2), F(1, 2), F(-1, 2)),
    (F(1, 2), F(-1, 2), F(-1, 2), F(1, 2)),
)
check(len(ice4) == 6, "ice fiber has exactly six states")
check(all(dot(axis, axis) == 1 for axis in axes4), "three ice axes have unit norm")
check(all(dot(axes4[i], axes4[j]) == 0 for i in range(3) for j in range(i)),
      "three ice axes are mutually orthogonal")

coordinates = tuple(tuple(dot(axis, state) for axis in axes4) for state in ice4)
expected_coordinates = {
    (F(2), F(0), F(0)), (F(-2), F(0), F(0)),
    (F(0), F(2), F(0)), (F(0), F(-2), F(0)),
    (F(0), F(0), F(2)), (F(0), F(0), F(-2)),
}
check(set(coordinates) == expected_coordinates, "ice fiber is exactly three antipodal doubled axes")
check(all(sum(value != 0 for value in state) == 1 for state in coordinates),
      "every ice state occupies exactly one axis")


# Exact uniform second and fourth moments.
def moment(indices: tuple[int, ...]) -> F:
    return sum(
        (prod(state[index] for index in indices) for state in coordinates),
        F(0),
    ) / len(coordinates)


def prod(values) -> F:
    result = F(1)
    for value in values:
        result *= value
    return result


cov0 = tuple(tuple(moment((i, j)) for j in range(3)) for i in range(3))
check(cov0 == mscale(F(4, 3), I3), "uniform Fisher covariance is exactly 4/3 I")
check(all(moment((i, j, k)) == 0 for i in range(3) for j in range(3) for k in range(3)),
      "all uniform cubic moments vanish")
first_derivative = tuple(
    tuple(tuple(
        moment((i, j, k))
        - moment((i, j)) * moment((k,))
        - moment((i, k)) * moment((j,))
        - moment((j, k)) * moment((i,))
        + 2 * moment((i,)) * moment((j,)) * moment((k,))
        for k in range(3)) for j in range(3)) for i in range(3)
)
check(all(value == 0 for matrix in first_derivative for row in matrix for value in row),
      "all 27 first derivatives of the uniform Fisher covariance vanish")
check(all(
    moment((i, j, k, ell)) == (F(16, 3) if i == j == k == ell else F(0))
    for i in range(3) for j in range(3) for k in range(3) for ell in range(3)
), "uniform fourth moment is exact axial delta 16/3")


# The Hessian from exact cumulants agrees with FN12 component by component.
for a, b, c, d in itertools.product(range(3), repeat=4):
    cumulant = (
        moment((a, b, c, d))
        - cov0[a][b] * cov0[c][d]
        - cov0[a][c] * cov0[b][d]
        - cov0[a][d] * cov0[b][c]
    )
    expected = (
        (F(16, 3) if a == b == c == d else F(0))
        - F(16, 9) * (int(a == b and c == d) + int(a == c and b == d) + int(a == d and b == c))
    )
    if cumulant != expected:
        raise AssertionError(f"Hessian component mismatch {(a, b, c, d)}")
check(True, "all 81 exact Fisher Hessian components match FN12")


def hessian_operator(source: tuple[tuple[F, ...], ...]) -> tuple[tuple[F, ...], ...]:
    trace = sum(source[i][i] for i in range(3))
    diagonal = diag(tuple(source[i][i] for i in range(3)))
    return madd(
        mscale(F(16, 3), diagonal),
        mscale(F(-16, 9) * trace, I3),
        mscale(F(-32, 9), source),
    )


a1 = I3
e_basis = (
    diag((F(1), F(-1), F(0))),
    diag((F(1), F(1), F(-2))),
)
t2_basis = (
    ((F(0), F(1), F(0)), (F(1), F(0), F(0)), (F(0), F(0), F(0))),
    ((F(0), F(0), F(1)), (F(0), F(0), F(0)), (F(1), F(0), F(0))),
    ((F(0), F(0), F(0)), (F(0), F(0), F(1)), (F(0), F(1), F(0))),
)
check(hessian_operator(a1) == mscale(F(-32, 9), a1), "Hessian A1 eigenvalue is -32/9")
check(all(hessian_operator(source) == mscale(F(16, 9), source) for source in e_basis),
      "Hessian E eigenvalue is 16/9")
check(all(hessian_operator(source) == mscale(F(-32, 9), source) for source in t2_basis),
      "Hessian T2 eigenvalue is -32/9")
sym2_basis = [a1, *e_basis, *t2_basis]
check(rank([flatten_sym(hessian_operator(source)) for source in sym2_basis]) == 6,
      "polarized Hessian is an isomorphism on Sym2")

# Replay the full S4 action inherited from coordinate permutations in R4.
# In the ice-axis basis every action is an exact signed permutation matrix.
s4_actions = []
for permutation in itertools.permutations(range(4)):
    action = tuple(
        tuple(
            dot(axes4[i], tuple(axes4[j][permutation[k]] for k in range(4)))
            for j in range(3)
        )
        for i in range(3)
    )
    s4_actions.append(action)
check(len(set(s4_actions)) == 24, "ice-axis representation realizes all 24 S4 actions faithfully")
check(all(mmul(action, transpose(action)) == I3 for action in s4_actions),
      "all inherited S4 actions are exactly orthogonal")
check(all(
    hessian_operator(mmul(mmul(action, source), transpose(action)))
    == mmul(mmul(action, hessian_operator(source)), transpose(action))
    for action in s4_actions for source in sym2_basis
), "Hessian intertwines every inherited S4 action")
check(F(16, 9) != F(-32, 9), "E and T2 ell=2 eigenvalues differ, excluding O3 equivariance")

for vector in (
    (F(1), F(0), F(0)),
    (F(1), F(1), F(0)),
    (F(1), F(2), F(3)),
    (F(-2), F(3), F(5)),
):
    source = outer(vector, vector)
    quadratic = mscale(F(1, 2), hessian_operator(source))
    check(sum(quadratic[i][i] for i in range(3)) == F(-16, 9) * dot(vector, vector),
          f"nonzero source {vector} has nonzero second coefficient")


# The function-module character contains T2 exactly once.
class_sizes = (1, 6, 3, 8, 6)
six_state_character = (6, 2, 2, 0, 0)
t2_character = (3, 1, -1, 0, -1)
multiplicity = sum(size * left * right for size, left, right in zip(
    class_sizes, six_state_character, t2_character
)) // 24
check(multiplicity == 1, "six-state diagonal function module contains one equivariant T2 query")


# Exact arbitrary-state covariance, complement boundary, and mixture identity.
p_plus = (F(1, 4), F(7, 30), F(1, 5))
p_minus = (F(1, 12), F(1, 10), F(2, 15))
check(sum(p_plus + p_minus, F(0)) == 1, "test distribution is normalized")
w = tuple(p_plus[i] + p_minus[i] for i in range(3))
m = tuple(2 * (p_plus[i] - p_minus[i]) for i in range(3))


def covariance_from_probabilities(plus: tuple[F, ...], minus: tuple[F, ...]):
    weights = tuple(plus[i] + minus[i] for i in range(3))
    mean = tuple(2 * (plus[i] - minus[i]) for i in range(3))
    return madd(mscale(F(4), diag(weights)), mscale(F(-1), outer(mean, mean)))


covariance = covariance_from_probabilities(p_plus, p_minus)
formula = madd(mscale(F(4), diag(w)), mscale(F(-1), outer(m, m)))
check(covariance == formula, "arbitrary-state covariance is 4 diag(w) minus mmT")
complement_covariance = covariance_from_probabilities(p_minus, p_plus)
check(complement_covariance == covariance, "complemented state has the same conditional covariance")
mixture_plus = tuple((p_plus[i] + p_minus[i]) / 2 for i in range(3))
mixture_covariance = covariance_from_probabilities(mixture_plus, mixture_plus)
check(mixture_covariance == madd(covariance, outer(m, m)),
      "unlabelled complement mixture adds between-branch mmT")
check(mixture_covariance == mscale(F(4), diag(w)),
      "unlabelled complement mixture cancels every T2 off-diagonal")

# Complement-preserving normalized tangents span E only; scalar adds A1.
dw1 = flatten_sym(mscale(F(4), diag((F(1), F(0), F(-1)))))
dw2 = flatten_sym(mscale(F(4), diag((F(0), F(1), F(-1)))))
scalar = flatten_sym(I3)
check(rank([dw1, dw2]) == 2, "complement-preserving normalized Fisher tangent has rank two")
check(rank([dw1, dw2, scalar]) == 3, "independent scalar raises symmetric tangent only to rank three")


def dm_column(mean: tuple[F, ...], index: int) -> tuple[F, ...]:
    direction = tuple(F(int(i == index)) for i in range(3))
    response = mscale(F(-1), madd(outer(mean, direction), outer(direction, mean)))
    return flatten_sym(response)


generic_mean = (F(1, 3), F(1, 4), F(1, 5))
dm_columns = [dm_column(generic_mean, index) for index in range(3)]
check(rank(dm_columns + [dw1, dw2]) == 5, "generic normalized broken-background metric map has rank five")
check(rank(dm_columns + [dw1, dw2, scalar]) == 6,
      "generic broken background plus independent scalar has rank six")

for nongeneric_mean in (
    (F(0), F(1, 4), F(1, 5)),
    (F(1, 3), F(0), F(1, 5)),
    (F(1, 3), F(1, 4), F(0)),
):
    columns = [dm_column(nongeneric_mean, index) for index in range(3)]
    check(rank(columns + [dw1, dw2, scalar]) < 6,
          f"vanishing mean component blocks full rank: {nongeneric_mean}")

# Symbolic determinant of the off-diagonal dm block, represented as a tiny
# monomial dictionary, is exactly 2 m1 m2 m3 after the overall minus sign.
# The underlying unsigned matrix has determinant -2 m1 m2 m3.
permutations = tuple(itertools.permutations(range(3)))
unsigned_entries = (
    ((1, F(1)), (0, F(1)), None),
    ((2, F(1)), None, (0, F(1))),
    (None, (2, F(1)), (1, F(1))),
)


def parity(permutation: tuple[int, ...]) -> int:
    inversions = sum(permutation[i] > permutation[j] for i in range(3) for j in range(i + 1, 3))
    return -1 if inversions % 2 else 1


polynomial: dict[tuple[int, int, int], F] = {}
for permutation in permutations:
    entries = [unsigned_entries[row][permutation[row]] for row in range(3)]
    if any(entry is None for entry in entries):
        continue
    exponents = [0, 0, 0]
    coefficient = F(parity(permutation))
    for variable, factor in entries:
        exponents[variable] += 1
        coefficient *= factor
    key = tuple(exponents)
    polynomial[key] = polynomial.get(key, F(0)) + coefficient
check(polynomial == {(1, 1, 1): F(-2)}, "unsigned off-diagonal determinant is -2 m1 m2 m3")
check(-polynomial[(1, 1, 1)] == 2, "response-block determinant is 2 m1 m2 m3")


print(f"SUMMARY {PASSED}/{PASSED} PASS")
