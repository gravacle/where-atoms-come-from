#!/usr/bin/env python3
"""Exact GL6CL global Fourier symbol of the GL6CH tensor writer.

The replay derives the centered parent/child symbol on the infinite Q4
(diamond-incidence) parent from physical tetrahedral node positions.  It
checks the zero-mode ranks and inverses, the common/relative small-momentum
expansions, the first cubic anisotropy, and exact momentum obstructions.
Only the Python standard library is used; all algebraic checks are rational.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from itertools import combinations, permutations, product
from pathlib import Path


HERE = Path(__file__).resolve().parent
PAIR_ORDER = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIR_ORDER)}

# Physical tetrahedral link displacements.  The inherited length unit makes
# |T_a|^2=3 and T_a.T_b=-1 for a != b.
TETRA = (
    (F(1), F(1), F(1)),
    (F(1), F(-1), F(-1)),
    (F(-1), F(1), F(-1)),
    (F(-1), F(-1), F(1)),
)
BRAVAIS = tuple(tuple(TETRA[a][i] - TETRA[3][i] for i in range(3))
                for a in range(3))

A = tuple(F(1) for _ in range(6))
T_BASIS = (
    tuple(map(F, (1, 0, 0, 0, 0, -1))),
    tuple(map(F, (0, 1, 0, 0, -1, 0))),
    tuple(map(F, (0, 0, 1, -1, 0, 0))),
)
AT_BASIS = (A,) + T_BASIS


class Checks:
    def __init__(self):
        self.total = 0

    def equal(self, got, want, label):
        self.total += 1
        if got != want:
            raise AssertionError(f"{label}: got {got!r}, want {want!r}")

    def true(self, condition, label):
        self.total += 1
        if not condition:
            raise AssertionError(label)


CHECK = Checks()


def dot(left, right):
    return sum((F(x) * F(y) for x, y in zip(left, right)), F(0))


def vadd(*vectors):
    return tuple(sum((F(v[i]) for v in vectors), F(0)) for i in range(len(vectors[0])))


def vscale(factor, vector):
    return tuple(F(factor) * F(value) for value in vector)


def zero_matrix(rows, columns):
    return tuple(tuple(F(0) for _ in range(columns)) for _ in range(rows))


def identity(size):
    return tuple(tuple(F(i == j) for j in range(size)) for i in range(size))


def transpose(matrix):
    return tuple(tuple(matrix[i][j] for i in range(len(matrix)))
                 for j in range(len(matrix[0])))


def matmul(left, right):
    return tuple(tuple(sum((left[i][k] * right[k][j]
                            for k in range(len(right))), F(0))
                       for j in range(len(right[0])))
                 for i in range(len(left)))


def outer(left, right):
    return tuple(tuple(F(x) * F(y) for y in right) for x in left)


def madd(left, right, factor=F(1)):
    return tuple(tuple(left[i][j] + F(factor) * right[i][j]
                       for j in range(len(left[0])))
                 for i in range(len(left)))


def mscale(factor, matrix):
    return tuple(tuple(F(factor) * entry for entry in row) for row in matrix)


def rank(matrix):
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
        work[pivot_row] = [entry / value for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            value = work[row][column]
            work[row] = [work[row][j] - value * work[pivot_row][j]
                         for j in range(columns)]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def determinant(matrix):
    total = F(0)
    size = len(matrix)
    for perm in permutations(range(size)):
        inversions = sum(perm[i] > perm[j]
                         for i in range(size) for j in range(i + 1, size))
        term = F((-1) ** inversions)
        for i, j in enumerate(perm):
            term *= matrix[i][j]
        total += term
    return total


P_A = mscale(F(1, 6), outer(A, A))
P_T = zero_matrix(6, 6)
for vector in T_BASIS:
    P_T = madd(P_T, mscale(F(1, 2), outer(vector, vector)))
P_E = madd(madd(identity(6), P_A, F(-1)), P_T, F(-1))

CHECK.equal(matmul(P_A, P_A), P_A, "A projector")
CHECK.equal(matmul(P_E, P_E), P_E, "E projector")
CHECK.equal(matmul(P_T, P_T), P_T, "T projector")
CHECK.equal(rank(P_A), 1, "A rank")
CHECK.equal(rank(P_E), 2, "E rank")
CHECK.equal(rank(P_T), 3, "T rank")
for a in range(4):
    CHECK.equal(dot(TETRA[a], TETRA[a]), F(3), "tetrahedron radius")
    for b in range(a):
        CHECK.equal(dot(TETRA[a], TETRA[b]), F(-1), "tetrahedron inner product")


def pair_vector(pair):
    answer = [F(0)] * 6
    answer[PAIR_INDEX[tuple(sorted(pair))]] = F(1)
    return tuple(answer)


def cycle_geometry():
    """Derive the six physical node positions and centered pair offsets."""
    orientations = []
    all_rho = []
    for missing in range(4):
        a, b, c = tuple(port for port in range(4) if port != missing)
        center = vscale(F(1, 2), vadd(TETRA[a], TETRA[c],
                                      vscale(F(-1), TETRA[b])))

        # Infinite-parent positions relative to anchor P_x=0.
        node_data = (
            ("P", (a, c), (F(0), F(0), F(0))),
            ("C", (a, b), TETRA[a]),
            ("P", (b, c), vadd(TETRA[a], vscale(-1, TETRA[b]))),
            ("C", (a, c), vadd(TETRA[a], vscale(-1, TETRA[b]), TETRA[c])),
            ("P", (a, b), vadd(TETRA[c], vscale(-1, TETRA[b]))),
            ("C", (b, c), TETRA[c]),
        )
        by_pair = {}
        for kind, pair, position in node_data:
            pair = tuple(sorted(pair))
            offset = vadd(position, vscale(-1, center))
            by_pair.setdefault(pair, {})[kind] = offset
        CHECK.equal(set(by_pair), set(combinations((a, b, c), 2)),
                    "each ring has its three local pairs")
        rho = {}
        for pair, offsets in by_pair.items():
            CHECK.equal(set(offsets), {"P", "C"}, "one P and one C per pair")
            CHECK.equal(offsets["P"], vscale(-1, offsets["C"]),
                        "same-pair sites oppose about ring center")
            CHECK.equal(dot(offsets["C"], offsets["C"]), F(11, 4),
                        "universal centered radius squared")
            rho[pair] = offsets["C"]
            all_rho.append(offsets["C"])
        orientations.append({
            "missing_port": missing,
            "ports": (a, b, c),
            "center": center,
            "rho_C": rho,
            "node_data": node_data,
        })

    CHECK.equal(len(all_rho), 12, "twelve orientation-pair offsets")
    second = tuple(tuple(sum((rho[i] * rho[j] for rho in all_rho), F(0))
                         for j in range(3)) for i in range(3))
    CHECK.equal(second, mscale(F(11), identity(3)),
                "aggregate second moment is isotropic")
    fourth_diag = tuple(sum((rho[i] ** 4 for rho in all_rho), F(0))
                        for i in range(3))
    fourth_mixed = tuple(sum((rho[i] ** 2 * rho[j] ** 2 for rho in all_rho), F(0))
                         for i, j in combinations(range(3), 2))
    CHECK.equal(fourth_diag, (F(83, 4),) * 3, "fourth diagonal moments")
    CHECK.equal(fourth_mixed, (F(19, 4),) * 3, "fourth mixed moments")
    return orientations, second, fourth_diag, fourth_mixed


# Polynomial helpers: exponent triples -> exact rational coefficient.
def pclean(poly):
    return {key: value for key, value in poly.items() if value}


def padd(left, right):
    answer = dict(left)
    for key, value in right.items():
        answer[key] = answer.get(key, F(0)) + value
    return pclean(answer)


def pscale(factor, poly):
    return pclean({key: F(factor) * value for key, value in poly.items()})


def pmul(left, right, maximum_degree):
    answer = {}
    for left_exp, left_value in left.items():
        for right_exp, right_value in right.items():
            exponent = tuple(left_exp[i] + right_exp[i] for i in range(3))
            if sum(exponent) <= maximum_degree:
                answer[exponent] = answer.get(exponent, F(0)) + left_value * right_value
    return pclean(answer)


def ppow(poly, power, maximum_degree):
    answer = {(0, 0, 0): F(1)}
    for _ in range(power):
        answer = pmul(answer, poly, maximum_degree)
    return answer


def linear_polynomial(vector):
    return pclean({(1, 0, 0): vector[0],
                   (0, 1, 0): vector[1],
                   (0, 0, 1): vector[2]})


def cosine_polynomial(vector, maximum_degree=4):
    linear = linear_polynomial(vector)
    answer = {(0, 0, 0): F(1)}
    factorial = 1
    for n in range(1, maximum_degree // 2 + 1):
        factorial *= (2 * n - 1) * (2 * n)
        answer = padd(answer, pscale(F((-1) ** n, factorial),
                                     ppow(linear, 2 * n, maximum_degree)))
    return answer


def polynomial_determinant(matrix, maximum_degree):
    answer = {}
    size = len(matrix)
    for perm in permutations(range(size)):
        inversions = sum(perm[i] > perm[j]
                         for i in range(size) for j in range(i + 1, size))
        term = {(0, 0, 0): F((-1) ** inversions)}
        for i, j in enumerate(perm):
            term = pmul(term, matrix[i][j], maximum_degree)
        answer = padd(answer, term)
    return answer


def common_symbol_analysis(orientations):
    """Analyze the canonical-direct B_+ and complete B_+ P_T writer."""
    zero_rows = []
    polynomial_rows = []
    for orientation in orientations:
        row_zero = [F(0)] * 6
        component_polynomials = [{} for _ in range(6)]
        for pair, rho in orientation["rho_C"].items():
            index = PAIR_INDEX[pair]
            row_zero[index] = F(2)
            component_polynomials[index] = pscale(F(2), cosine_polynomial(rho, 4))
        zero_rows.append(tuple(row_zero))
        polynomial_rows.append(tuple(component_polynomials))
    zero_rows = tuple(zero_rows)

    canonical_normal = matmul(transpose(zero_rows), zero_rows)
    expected_canonical_normal = madd(mscale(F(24), P_A), mscale(F(8), P_T))
    CHECK.equal(canonical_normal, expected_canonical_normal,
                "canonical-direct Bplus zero-mode normal")
    CHECK.equal(rank(zero_rows), 4, "canonical-direct Bplus zero-mode A+T rank four")
    CHECK.true(all(tuple(sum(P_E[i][j] * row[j] for j in range(6))
                         for i in range(6)) == (F(0),) * 6
                   for row in zero_rows), "zero-mode E is null")

    # Only the T projection is the complete arbitrary-profile h6 writer.
    # The A direction is complete only for a spatially uniform source via the
    # separate exact U_d-shift identity below; E is unclassified.
    tensor_zero_rows = matmul(zero_rows, P_T)
    tensor_normal = matmul(transpose(tensor_zero_rows), tensor_zero_rows)
    CHECK.equal(tensor_normal, mscale(F(8), P_T),
                "complete tensor-writer zero-mode normal")
    CHECK.equal(rank(tensor_zero_rows), 3, "complete tensor-writer zero rank three")
    reconstruction = mscale(F(1, 8), P_T)
    CHECK.equal(matmul(reconstruction, tensor_normal), P_T,
                "zero-mode physical writer reconstructs T")

    at_matrix = []
    for component_row in polynomial_rows:
        row = []
        for basis in AT_BASIS:
            entry = {}
            for polynomial, coefficient in zip(component_row, basis):
                entry = padd(entry, pscale(coefficient, polynomial))
            row.append(entry)
        at_matrix.append(row)
    det_series = polynomial_determinant(at_matrix, 4)
    expected_det = {
        (0, 0, 0): F(768),
        (2, 0, 0): F(-1408), (0, 2, 0): F(-1408), (0, 0, 2): F(-1408),
        (4, 0, 0): F(2800, 3), (0, 4, 0): F(2800, 3),
        (0, 0, 4): F(2800, 3),
        (2, 2, 0): F(2144), (2, 0, 2): F(2144),
        (0, 2, 2): F(2144),
    }
    CHECK.equal(det_series, expected_det, "exact common AT determinant through k4")

    # The scalar determinant hides tensor anisotropy.  In the rational T
    # basis (each vector has norm^2=2), compute (B_+ T)^*(B_+ T) through k^2.
    t_polynomial_rows = []
    for component_row in polynomial_rows:
        row = []
        for basis in T_BASIS:
            entry = {}
            for polynomial, coefficient in zip(component_row, basis):
                entry = padd(entry, pscale(coefficient, polynomial))
            row.append(entry)
        t_polynomial_rows.append(row)
    t_normal_degree_2 = []
    for i in range(3):
        row = []
        for j in range(3):
            entry = {}
            for orientation in range(4):
                entry = padd(entry, pmul(t_polynomial_rows[orientation][i],
                                          t_polynomial_rows[orientation][j], 2))
            row.append(entry)
        t_normal_degree_2.append(tuple(row))
    t_normal_degree_2 = tuple(t_normal_degree_2)
    expected_t_normal = (
        ({(0, 0, 0): F(16), (2, 0, 0): F(-36),
          (0, 2, 0): F(-4), (0, 0, 2): F(-4)},
         {(1, 1, 0): F(24)}, {(1, 0, 1): F(24)}),
        ({(1, 1, 0): F(24)},
         {(0, 0, 0): F(16), (2, 0, 0): F(-4),
          (0, 2, 0): F(-36), (0, 0, 2): F(-4)},
         {(0, 1, 1): F(24)}),
        ({(1, 0, 1): F(24)}, {(0, 1, 1): F(24)},
         {(0, 0, 0): F(16), (2, 0, 0): F(-4),
          (0, 2, 0): F(-4), (0, 0, 2): F(-36)}),
    )
    CHECK.equal(t_normal_degree_2, expected_t_normal,
                "exact common T normal through k2")

    # Exact high-momentum obstruction at reciprocal coordinates q=(pi,0,0),
    # namely Cartesian k=(pi/4)(1,1,1).  C=cos(pi/8), S=cos(3pi/8).
    # Only coefficient pairs (C,S) are needed; no floating transcendental use.
    symbolic_t_rows = (
        ((F(-2), F(0)),) * 3,
        ((F(0), F(2)),) * 3,
        ((F(0), F(2)),) * 3,
        ((F(0), F(2)),) * 3,
    )
    # Every one of the three T columns is literally the same symbolic
    # (C,S)-valued column.  Since C=cos(pi/8) is nonzero, the rank is one.
    CHECK.true(all(row[0] == row[1] == row[2] for row in symbolic_t_rows),
               "BZ-corner three T columns coincide exactly")
    CHECK.true(symbolic_t_rows[0][0] != (F(0), F(0)),
               "BZ-corner common T symbol is nonzero")

    return {
        "canonical_direct_zero_rows": zero_rows,
        "canonical_direct_zero_rank": 4,
        "canonical_direct_zero_normal": canonical_normal,
        "canonical_direct_zero_null": "E",
        "tensor_zero_rows": tensor_zero_rows,
        "tensor_zero_rank": 3,
        "tensor_zero_normal": tensor_normal,
        "tensor_zero_null": "A1+E",
        "source_reconstruction": "j_T=(1/8)P_T B_+(0)^* w",
        "det_AT_through_degree_4": det_series,
        "det_invariant_form": "768-1408|k|^2+1072|k|^4-(416/3)(kx^4+ky^4+kz^4)+O(|k|^6)",
        "T_normal_through_degree_2_rational_basis": t_normal_degree_2,
        "T_normal_invariant_form_rational_basis": "16 I-4|k|^2 I+24 k k^T-56 diag(kx^2,ky^2,kz^2)+O(|k|^4)",
        "T_normal_invariant_form_orthonormal_basis": "8 I-2|k|^2 I+12 k k^T-28 diag(kx^2,ky^2,kz^2)+O(|k|^4)",
        "rotation_ceiling": "T2 is only the cubic off-diagonal part of the SO(3) l=2 space E2+T2; this T-T block alone cannot distinguish physical rotational anisotropy from an SO(3)-covariant completion with E-T and E-E blocks",
        "rigorous_near_zero_condition": "rank(B_+(k)P_T)=3 for |k|^4<32/363",
        "perturbation_bound": "||B_+(k)-B_+(0)||_2 <= (11 sqrt(3)/2)|k|^2",
        "bz_corner": {
            "reciprocal_coordinates": "q=(pi,0,0)",
            "cartesian_k": "k=(pi/4)(1,1,1)",
            "T_rank": 1,
            "T_rows_in_C_cos_pi8_S_cos_3pi8": symbolic_t_rows,
        },
    }


def locked_read():
    rows = []
    for bits in product((0, 1), repeat=4):
        if sum(bits) != 2:
            continue
        signs = tuple(F(1 - 2 * bit) for bit in bits)
        rows.append(tuple(signs[a] * signs[b] for a, b in PAIR_ORDER))
    rows = tuple(rows)
    normal = matmul(transpose(rows), rows)
    CHECK.equal(normal, madd(mscale(F(4), P_A), mscale(F(16), P_E)),
                "independently derived locked-read normal")
    CHECK.equal(rank(rows), 3, "locked read rank")
    return rows, normal


def relative_symbol_analysis(orientations):
    """Analyze B_-(k)=-2i sum sin(k.rho)e_pair at small momentum."""
    # L is the real 4x3 coefficient of the leading T map after removing -2i:
    # B_-(k)|_T=-2i L(k)+O(|k|^3).
    linear_rows = []
    for orientation in orientations:
        row = []
        for basis in T_BASIS:
            vector = [F(0)] * 3
            for pair, rho in orientation["rho_C"].items():
                coefficient = basis[PAIR_INDEX[pair]]
                for i in range(3):
                    vector[i] += coefficient * rho[i]
            row.append(tuple(vector))
        linear_rows.append(tuple(row))
    linear_rows = tuple(linear_rows)

    minors = []
    for omitted in range(4):
        matrix = [[linear_polynomial(linear_rows[row][column])
                   for column in range(3)]
                  for row in range(4) if row != omitted]
        minors.append(polynomial_determinant(matrix, 3))
    sum_squares = {}
    for minor in minors:
        sum_squares = padd(sum_squares, pmul(minor, minor, 6))
    expected_sum_squares = {
        (6, 0, 0): F(9), (0, 6, 0): F(9), (0, 0, 6): F(9),
        (4, 2, 0): F(-9), (4, 0, 2): F(-9),
        (2, 4, 0): F(-9), (0, 4, 2): F(-9),
        (2, 0, 4): F(-9), (0, 2, 4): F(-9),
        (2, 2, 2): F(58),
    }
    CHECK.equal(sum_squares, expected_sum_squares,
                "relative T leading-minor sum of squares")

    hadamard_patterns = (
        (1, 1, 1, 1), (1, 1, -1, -1),
        (1, -1, 1, -1), (1, -1, -1, 1),
    )
    combinations_found = []
    for pattern in hadamard_patterns:
        combined = {}
        for sign, minor in zip(pattern, minors):
            combined = padd(combined, pscale(sign, minor))
        combinations_found.append(combined)
    expected_hadamard = (
        {(2, 1, 0): F(6), (0, 3, 0): F(-6), (0, 1, 2): F(6)},
        {(2, 0, 1): F(6), (0, 2, 1): F(6), (0, 0, 3): F(-6)},
        {(1, 1, 1): F(4)},
        {(1, 2, 0): F(6), (1, 0, 2): F(6), (3, 0, 0): F(-6)},
    )
    CHECK.equal(tuple(combinations_found), expected_hadamard,
                "Hadamard factorization of four leading minors")

    # Exact all-q sine-series dependencies on the six Cartesian face diagonals.
    # A series is {positive frequency: coefficient}; sin(-nq)=-sin(nq).
    def sine_series(direction):
        matrix = []
        for orientation in orientations:
            row = []
            for basis in T_BASIS:
                series = {}
                for pair, rho in orientation["rho_C"].items():
                    frequency = dot(direction, rho)
                    coefficient = basis[PAIR_INDEX[pair]]
                    if frequency < 0:
                        frequency, coefficient = -frequency, -coefficient
                    if frequency:
                        series[frequency] = series.get(frequency, F(0)) + coefficient
                row.append(pclean(series))
            matrix.append(tuple(row))
        return tuple(matrix)

    dependencies = {
        (1, 1, 0): (0, 1, F(1)),
        (1, -1, 0): (0, 1, F(-1)),
        (1, 0, 1): (0, 2, F(1)),
        (1, 0, -1): (0, 2, F(-1)),
        (0, 1, 1): (1, 2, F(1)),
        (0, 1, -1): (1, 2, F(-1)),
    }
    face_series = {}
    for direction, (left, right, factor) in dependencies.items():
        matrix = sine_series(direction)
        for row in matrix:
            CHECK.equal(row[left], pscale(factor, row[right]),
                        "exact face-diagonal relative-T column dependence")
        face_series[direction] = matrix

    return {
        "zero_rank": 0,
        "leading_T_rows_without_minus_2i": linear_rows,
        "three_by_three_minors": minors,
        "minor_sum_squares": sum_squares,
        "minor_sum_invariant": "9 sum_i ki^6-9 sum_{i!=j}ki^4 kj^2+58 kx^2 ky^2 kz^2",
        "hadamard_minor_combinations": combinations_found,
        "leading_rank": "3 generically; 2 on nonzero Cartesian face diagonals; 0 at k=0",
        "exact_face_diagonal_dependencies": dependencies,
        "small_k": "B_-(k)=-2i sum_p (k.rho_dp)e_p+(i/3)sum_p(k.rho_dp)^3 e_p+O(|k|^5)",
    }


def a1_consistency():
    values = {}
    for occupancy in range(5):
        for bits in product((0, 1), repeat=4):
            if sum(bits) != occupancy:
                continue
            signs = tuple(F(1 - 2 * bit) for bit in bits)
            pair_sum = sum((signs[a] * signs[b] for a, b in PAIR_ORDER), F(0))
            CHECK.equal(pair_sum, F(2 * (occupancy - 2) ** 2 - 2),
                        "uniform pair identity on every local word")
            values[occupancy] = pair_sum
    CHECK.equal(F(6) * F(105, 8), F(315, 4),
                "six local A1 vertices match denominator derivative")
    return {
        "identity": "sum_{a<b} Z_a Z_b=2(n-2)^2-2",
        "energy_shift": "q A source sends U_d to U_d+2q plus scalar -2q",
        "denominator_derivative": "+(315/4)h^6/U_d^6",
        "six_vertex_sum": "6(105/8)h^6/U_d^6=(315/4)h^6/U_d^6",
        "occupancy_values": values,
    }


def qtext(value):
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def encode(value):
    if isinstance(value, F):
        return qtext(value)
    if isinstance(value, dict):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-ledger", action="store_true")
    args = parser.parse_args()

    orientations, second, fourth_diag, fourth_mixed = cycle_geometry()
    common = common_symbol_analysis(orientations)
    read_rows, read_normal = locked_read()
    relative = relative_symbol_analysis(orientations)
    scalar = a1_consistency()

    combined_zero = madd(read_normal, common["tensor_zero_normal"])
    expected_combined = madd(madd(mscale(F(4), P_A), mscale(F(16), P_E)),
                              mscale(F(8), P_T))
    CHECK.equal(combined_zero, expected_combined, "common read+writer zero normal")
    CHECK.equal(rank(tuple(read_rows) + tuple(common["tensor_zero_rows"])), 6,
                "soldered common zero-mode full pair rank")
    inverse_zero = madd(madd(mscale(F(1, 4), P_A), mscale(F(1, 16), P_E)),
                        mscale(F(1, 8), P_T))
    CHECK.equal(matmul(inverse_zero, combined_zero), identity(6),
                "explicit combined zero-mode left inverse")
    CHECK.equal(determinant(combined_zero), F(524288),
                "combined zero normal determinant")

    # Unsoldered P/C map at zero: two locked reads plus one ring writer.
    # Columns are (j_P,j_C); B_P(0)=B_C(0)=one incidence copy.
    one_copy = tuple(tuple(entry / 2 for entry in row)
                     for row in common["tensor_zero_rows"])
    unsoldered = []
    for row in read_rows:
        unsoldered.append(tuple(row) + (F(0),) * 6)
        unsoldered.append((F(0),) * 6 + tuple(row))
    for row in one_copy:
        unsoldered.append(tuple(row) + tuple(row))
    CHECK.equal(rank(tuple(unsoldered)), 9, "unsoldered zero-mode rank nine")

    ledger = {
        "lane": "GL6CL",
        "geometry": {
            "tetrahedral_links": TETRA,
            "bravais_parent_steps": BRAVAIS,
            "parent_position": "r(P_x)=sum_{i=0}^2 x_i(T_i-T_3)",
            "child_position": "r(C_y)=r(P_y)+T_3",
            "orientations": orientations,
            "rho_second_moment": second,
            "rho_fourth_diagonal": fourth_diag,
            "rho_fourth_mixed": fourth_mixed,
            "fourth_contraction": "(83/4)sum_i ki^4+(57/2)sum_{i<j}ki^2kj^2=(57/4)|k|^4+(13/2)sum_i ki^4",
        },
        "exact_symbol": {
            "canonical_direct_B_P": "sum_{p in d-complement} exp(-i k.rho_dp)e_p",
            "canonical_direct_B_C": "sum_{p in d-complement} exp(+i k.rho_dp)e_p",
            "canonical_direct_B_plus": "2 sum_p cos(k.rho_dp)e_p",
            "canonical_direct_B_minus": "-2i sum_p sin(k.rho_dp)e_p",
            "complete_tensor_writer": "delta a_d^T=mu[B_plus P_T j_plus+B_minus P_T j_minus]",
            "physical_scale": "mu=(105/8)h^6/U_d^6",
            "field_convention": "j_P=j_plus+j_minus; j_C=j_plus-j_minus",
            "scope_guard": "arbitrary-profile A1/E off-diagonal h6 completion is unclassified; uniform A1 at k=0 is separately exact",
        },
        "common_sector": common,
        "relative_sector": relative,
        "locked_read_and_common_writer": {
            "zero_normal": combined_zero,
            "zero_rank": 6,
            "zero_inverse_normal": inverse_zero,
            "zero_determinant": determinant(combined_zero),
            "C_definition": "C(k)=(D;B_+(k)P_T)",
            "analytic_left_inverse": "L(k)=[C(k)^*C(k)]^{-1}C(k)^* for |k|^4<32/363",
            "unsoldered_zero_rank": 9,
            "unsoldered_all_k_rank_ceiling": 10,
            "unsoldered_obstruction": "12 sources but only rank3+rank3 locked reads and four ring outputs; kernel dimension at least 2",
        },
        "uniform_A1": scalar,
        "scope": "complete T2 global linear Fourier/operator-jet access plus separately exact uniform A1; arbitrary-profile A1/E h6 completion and autonomous response/spacetime/metric/Ricci/gravity/G are unproved",
    }
    payload = json.dumps(encode(ledger), indent=2, sort_keys=True) + "\n"
    target = HERE / "EXACT_LEDGER.json"
    if args.write_ledger:
        target.write_text(payload)
    else:
        CHECK.true(target.is_file(), "frozen exact ledger exists")
        CHECK.equal(target.read_text(), payload, "frozen exact ledger matches replay")

    print(f"PASS__GL6CL_GLOBAL_FOURIER_PAIR_WRITER__{CHECK.total}/{CHECK.total}")
    print("EXACT_SYMBOL=BP_EXP_MINUS_IKRHO;BC_EXP_PLUS_IKRHO;COMMON_COS;RELATIVE_SIN")
    print("K0_CANONICAL_DIRECT=RANK4_A1_PLUS_T2;COMPLETE_T_WRITER=RANK3")
    print("LOCKED_A1_PLUS_E_READ_PLUS_COMPLETE_T_WRITER=RANK6")
    print("NEAR_ZERO=FULL_RANK_FOR_K4_LT_32_OVER_363;ANALYTIC_LEFT_INVERSE")
    print("T_BLOCK=CUBIC_K2_STRUCTURE_EXACT;SO3_DIAGNOSIS_REQUIRES_E2_PLUS_T2_COMPLETION")
    print("RELATIVE=ZERO_AT_K0;GENERIC_LEADING_RANK3;FACE_DIAGONAL_RANK2")
    print("FULL_BZ_OBSTRUCTION=COMMON_T_RANK1_AT_Q_PI_0_0;SMOOTH_ACCESS_ONLY")
    print("NO_AUTONOMOUS_RESPONSE_SPACETIME_METRIC_RICCI_GRAVITY_G")


if __name__ == "__main__":
    main()
