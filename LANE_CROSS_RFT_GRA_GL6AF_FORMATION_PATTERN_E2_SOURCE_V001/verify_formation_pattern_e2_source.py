#!/usr/bin/env python3
"""Exact standard-library replay of the GL6AF formation-pattern theorem."""

from itertools import combinations, product


checks = 0


def require(condition, message):
    global checks
    if not condition:
        raise AssertionError(message)
    checks += 1


# Polynomials in the formed-link Bloch z, represented coefficient-wise.
def add(left, right):
    out = dict(left)
    for degree, value in right.items():
        out[degree] = out.get(degree, 0) + value
    return {degree: value for degree, value in out.items() if value}


def scale(poly, value):
    return {degree: value * coefficient for degree, coefficient in poly.items()
            if value * coefficient}


def multiply(left, right):
    out = {}
    for a, avalue in left.items():
        for b, bvalue in right.items():
            out[a + b] = out.get(a + b, 0) + avalue * bvalue
    return {degree: value for degree, value in out.items() if value}


PAIRS = tuple(combinations(range(4), 2))
E = (
    (1, 1),
    (-1, 0),
    (0, -1),
    (0, -1),
    (-1, 0),
    (1, 1),
)


def response_without_minus4hx(pattern):
    """Return D/(-4 h x) as a 6x6 z-polynomial matrix."""
    zvalue = ({1: 1} if formed else {0: 1} for formed in pattern)
    zvalue = tuple(zvalue)
    matrix = [[{} for _ in PAIRS] for _ in PAIRS]
    for source_index, source in enumerate(PAIRS):
        for read_index, read in enumerate(PAIRS):
            if source == read:
                matrix[read_index][source_index] = {
                    0: pattern[source[0]] + pattern[source[1]]}
                continue
            shared = set(source) & set(read)
            if len(shared) != 1:
                continue
            common = next(iter(shared))
            source_other = next(iter(set(source) - {common}))
            read_other = next(iter(set(read) - {common}))
            matrix[read_index][source_index] = scale(
                multiply(zvalue[source_other], zvalue[read_other]),
                pattern[common],
            )
    return matrix


def project_e(matrix):
    out = [[{} for _ in range(2)] for _ in range(2)]
    for read_basis in range(2):
        for source_basis in range(2):
            value = {}
            for read in range(6):
                for source in range(6):
                    value = add(value, scale(
                        matrix[read][source],
                        E[read][read_basis] * E[source][source_basis],
                    ))
            out[read_basis][source_basis] = value
    return out


def outer(vector):
    return [[vector[row] * vector[column] for column in range(2)]
            for row in range(2)]


zero = [[{} for _ in range(2)] for _ in range(2)]
rank_by_count = {0: set(), 1: set(), 2: set(), 3: set(), 4: set()}

for pattern in product((0, 1), repeat=4):
    matrix = response_without_minus4hx(pattern)
    projected = project_e(matrix)
    count = sum(pattern)
    if count <= 1:
        require(projected == zero, "zero/one formed E null")
        rank_by_count[count].add(0)
    elif count == 2:
        formed_pair = tuple(index for index, formed in enumerate(pattern)
                            if formed)
        row = E[PAIRS.index(formed_pair)]
        expected = [[scale({0: 1, 1: -1}, 4 * entry)
                     for entry in outer(row)[matrix_row]]
                    for matrix_row in range(2)]
        require(projected == expected, "two-formed rank-one block")
        determinant = add(
            multiply(projected[0][0], projected[1][1]),
            scale(multiply(projected[0][1], projected[1][0]), -1),
        )
        require(not determinant, "two-formed determinant zero")
        require(any(projected[row_index][column_index]
                    for row_index in range(2) for column_index in range(2)),
                "two-formed block nonzero polynomial")
        rank_by_count[count].add(1)
    else:
        determinant = add(
            multiply(projected[0][0], projected[1][1]),
            scale(multiply(projected[0][1], projected[1][0]), -1),
        )
        require(bool(determinant), "three/four formed determinant nonzero")
        rank_by_count[count].add(2)

require(rank_by_count == {0: {0}, 1: {0}, 2: {1}, 3: {2}, 4: {2}},
        "rank threshold census")

# Representative closed forms after removing the common -4 h x.
three = project_e(response_without_minus4hx((1, 1, 1, 0)))
three_scalar = {0: 3, 1: -2, 2: -1}
three_expected = [[scale(three_scalar, value) for value in (2, 1)],
                  [scale(three_scalar, value) for value in (1, 2)]]
require(three == three_expected, "three-record representative")

four = project_e(response_without_minus4hx((1, 1, 1, 1)))
four_scalar = {0: 4, 2: -4}
four_expected = [[scale(four_scalar, value) for value in (2, 1)],
                 [scale(four_scalar, value) for value in (1, 2)]]
require(four == four_expected, "four-record representative")

# The six two-record branches use exactly three noncollinear covector lines.
lines = set()
for pair in PAIRS:
    row = E[PAIRS.index(pair)]
    canonical = row
    if canonical[0] < 0 or (canonical[0] == 0 and canonical[1] < 0):
        canonical = tuple(-value for value in canonical)
    lines.add(canonical)
require(len(lines) == 3, "three distinct two-record lines")
require(any(a * d - b * c != 0 for (a, b) in lines for (c, d) in lines),
        "two-record lines span E2")

print(f"PASS GL6AF exact formation-pattern checks {checks}/{checks}")
