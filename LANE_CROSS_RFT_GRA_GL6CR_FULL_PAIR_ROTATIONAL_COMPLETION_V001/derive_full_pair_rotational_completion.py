#!/usr/bin/env python3
"""Exact full-six-pair cubic-to-rotational completion calculation.

The calculation classifies reciprocal, inversion-even, S4-covariant
quadratic symbols on the six tetrahedral pair coordinates.  It then uses the
inherited EW/GJ pair-to-symmetric-tensor solder to identify the exact SO(3)
subspace.  Nothing in this replay assumes a stationary phase, a continuum,
or an Einstein endpoint.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path


HERE = Path(__file__).resolve().parent
PAIRS = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIRS)}
PAIR_SYM = tuple((i, j) for i in range(6) for j in range(i, 6))
MONOMIALS = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
TETRA = (
    (F(1), F(1), F(1)),
    (F(1), F(-1), F(-1)),
    (F(-1), F(1), F(-1)),
    (F(-1), F(-1), F(1)),
)


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


def zero(rows, columns):
    return tuple(tuple(F(0) for _ in range(columns)) for _ in range(rows))


def identity(size):
    return tuple(tuple(F(i == j) for j in range(size)) for i in range(size))


def transpose(matrix):
    return tuple(tuple(matrix[i][j] for i in range(len(matrix)))
                 for j in range(len(matrix[0])))


def matmul(left, right):
    return tuple(tuple(sum((left[i][q] * right[q][j]
                            for q in range(len(right))), F(0))
                       for j in range(len(right[0])))
                 for i in range(len(left)))


def matvec(matrix, vector):
    return tuple(sum((matrix[i][j] * vector[j]
                      for j in range(len(vector))), F(0))
                 for i in range(len(matrix)))


def mscale(factor, matrix):
    return tuple(tuple(F(factor) * value for value in row) for row in matrix)


def madd(left, right):
    return tuple(tuple(left[i][j] + right[i][j]
                       for j in range(len(left[0])))
                 for i in range(len(left)))


def outer(left, right):
    return tuple(tuple(F(x) * F(y) for y in right) for x in left)


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


def solve_columns(columns, target):
    """Solve a full-column-rank rectangular system columns*x=target."""
    rows = len(target)
    count = len(columns)
    work = [[F(columns[column][row]) for column in range(count)]
            + [F(target[row])] for row in range(rows)]
    pivot_rows = []
    pivot_row = 0
    for column in range(count):
        pivot = next((row for row in range(pivot_row, rows)
                      if work[row][column]), None)
        if pivot is None:
            raise AssertionError("basis columns are not independent")
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][column]
        work[pivot_row] = [entry / value for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            value = work[row][column]
            work[row] = [work[row][j] - value * work[pivot_row][j]
                         for j in range(count + 1)]
        pivot_rows.append(pivot_row)
        pivot_row += 1
    for row in range(rows):
        if all(work[row][column] == 0 for column in range(count)):
            if work[row][-1] != 0:
                raise AssertionError("target is outside column span")
    return tuple(work[pivot_rows[column]][-1] for column in range(count))


def nullspace(matrix):
    """Exact right nullspace, with matrix supplied as rows."""
    work = [list(map(F, row)) for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    pivots = []
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
        pivots.append(column)
        pivot_row += 1
        if pivot_row == rows:
            break
    free = [column for column in range(columns) if column not in pivots]
    answer = []
    for free_column in free:
        vector = [F(0)] * columns
        vector[free_column] = F(1)
        for row, pivot in reversed(list(enumerate(pivots))):
            vector[pivot] = -sum((work[row][column] * vector[column]
                                  for column in free), F(0))
        answer.append(tuple(vector))
    return tuple(answer)


def permutation_matrix(permutation, size):
    matrix = [[F(0)] * size for _ in range(size)]
    for source, target in enumerate(permutation):
        matrix[target][source] = F(1)
    return tuple(tuple(row) for row in matrix)


def group_data():
    tetra_columns = transpose(TETRA[:3])
    # Signed tetrahedral matrices allow a simple exact inverse: G^{-1}=G^T.
    gram = matmul(transpose(tetra_columns), tetra_columns)
    CHECK.equal(gram, tuple(tuple(F(3 if i == j else -1)
                                  for j in range(3)) for i in range(3)),
                "tetrahedron Gram")

    # Generic exact square inverse.
    def inverse(matrix):
        size = len(matrix)
        work = [list(matrix[i]) + list(identity(size)[i])
                for i in range(size)]
        for column in range(size):
            pivot = next(row for row in range(column, size)
                         if work[row][column])
            work[column], work[pivot] = work[pivot], work[column]
            value = work[column][column]
            work[column] = [entry / value for entry in work[column]]
            for row in range(size):
                if row == column:
                    continue
                value = work[row][column]
                work[row] = [work[row][j] - value * work[column][j]
                             for j in range(2 * size)]
        return tuple(tuple(row[size:]) for row in work)

    tetra_inverse = inverse(tetra_columns)
    group = []
    for permutation in permutations(range(4)):
        target_columns = transpose(tuple(TETRA[permutation[a]]
                                         for a in range(3)))
        rotation = matmul(target_columns, tetra_inverse)
        CHECK.equal(matmul(transpose(rotation), rotation), identity(3),
                    "port action is orthogonal")
        pair_permutation = []
        for pair in PAIRS:
            image = tuple(sorted((permutation[pair[0]], permutation[pair[1]])))
            pair_permutation.append(PAIR_INDEX[image])
        pair_action = permutation_matrix(tuple(pair_permutation), 6)
        group.append((permutation, rotation, pair_action))
    CHECK.equal(len(group), 24, "S4 group order")
    return tuple(group)


def elementary_symmetric(size, i, j):
    matrix = [[F(0)] * size for _ in range(size)]
    matrix[i][j] = F(1)
    matrix[j][i] = F(1)
    if i == j:
        matrix[i][j] = F(1)
    return tuple(tuple(row) for row in matrix)


def momentum_matrix(monomial):
    i, j = monomial
    matrix = [[F(0)] * 3 for _ in range(3)]
    if i == j:
        matrix[i][i] = F(1)
    else:
        matrix[i][j] = matrix[j][i] = F(1, 2)
    return tuple(tuple(row) for row in matrix)


def momentum_coefficients(matrix):
    return tuple(matrix[i][j] if i == j else 2 * matrix[i][j]
                 for i, j in MONOMIALS)


def momentum_from_coefficients(coefficients):
    matrix = [[F(0)] * 3 for _ in range(3)]
    for value, (i, j) in zip(coefficients, MONOMIALS):
        if i == j:
            matrix[i][i] = F(value)
        else:
            matrix[i][j] = matrix[j][i] = F(value) / 2
    return tuple(tuple(row) for row in matrix)


def vector_from_factor(pair_matrix, momentum_quadratic):
    momentum = momentum_coefficients(momentum_quadratic)
    return tuple(pair_matrix[i][j] * momentum[m]
                 for i, j in PAIR_SYM for m in range(6))


def vadd(left, right):
    return tuple(left[i] + right[i] for i in range(len(left)))


def vscale(factor, vector):
    return tuple(F(factor) * value for value in vector)


def reynolds(seed_pair, seed_momentum, group):
    answer = (F(0),) * (len(PAIR_SYM) * len(MONOMIALS))
    for _, rotation, pair_action in group:
        # K(k)=sum_g P_g^T F(R_g k) P_g is fixed by the covariance action.
        pair_image = matmul(matmul(transpose(pair_action), seed_pair),
                            pair_action)
        momentum_image = matmul(matmul(transpose(rotation), seed_momentum),
                                rotation)
        answer = vadd(answer, vector_from_factor(pair_image, momentum_image))
    return answer


def invariant_basis(group):
    basis = []
    labels = []
    for i, j in PAIR_SYM:
        pair_seed = elementary_symmetric(6, i, j)
        for monomial in MONOMIALS:
            candidate = reynolds(pair_seed, momentum_matrix(monomial), group)
            if not any(candidate):
                continue
            old_rank = rank(basis)
            if rank(basis + [candidate]) > old_rank:
                basis.append(candidate)
                labels.append(f"Reynolds[K_{i}{j}*k_{monomial[0]}k_{monomial[1]}]")
    CHECK.equal(len(basis), 9, "full six-pair S4 quadratic invariant dimension")
    return tuple(basis), tuple(labels)


def covariance_sides(vector, rotation, pair_action):
    """Return coefficient vectors for K(Rk) and P K(k) P^T."""
    left = (F(0),) * len(vector)
    right = (F(0),) * len(vector)
    for pair_slot, (a, b) in enumerate(PAIR_SYM):
        coefficients = vector[pair_slot * 6:(pair_slot + 1) * 6]
        momentum = momentum_from_coefficients(coefficients)
        pair_matrix = elementary_symmetric(6, a, b)
        left = vadd(left, vector_from_factor(
            pair_matrix,
            matmul(matmul(transpose(rotation), momentum), rotation)))
        right = vadd(right, vector_from_factor(
            matmul(matmul(pair_action, pair_matrix), transpose(pair_action)),
            momentum))
    return left, right


def symmetric_tensor_solder():
    vertices = tuple(tuple(entry / 2 for entry in vector) for vector in TETRA)
    tensors = []
    for a, b in PAIRS:
        tensors.append(tuple(tuple(vertices[a][i] * vertices[b][j]
                                   + vertices[b][i] * vertices[a][j]
                                   for j in range(3)) for i in range(3)))
    coordinate_matrix = tuple(tuple(
        tensors[column][i][j]
        for column in range(6))
        for i, j in ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)))
    CHECK.equal(rank(coordinate_matrix), 6, "pair-to-Sym2 solder rank six")
    return tuple(tensors), coordinate_matrix


def verify_solder_equivariance(group, tensors):
    for permutation, rotation, _ in group:
        for source, (a, b) in enumerate(PAIRS):
            image_pair = tuple(sorted((permutation[a], permutation[b])))
            image = PAIR_INDEX[image_pair]
            transformed = matmul(matmul(rotation, tensors[source]),
                                  transpose(rotation))
            CHECK.equal(transformed, tensors[image],
                        "pair-memory solder is S4-equivariant")


def trace(tensor):
    return sum((tensor[i][i] for i in range(3)), F(0))


def frobenius(left, right):
    return sum((left[i][j] * right[i][j]
                for i in range(3) for j in range(3)), F(0))


def sym_product(left, right):
    product = matmul(left, right)
    reverse = matmul(right, left)
    return mscale(F(1, 2), madd(product, reverse))


def rotational_bases(tensors):
    """Four self-adjoint SO(3) quadratic bilinears on Sym2(R3)."""
    bases = []
    labels = (
        "r2_Frobenius",
        "r2_trace_trace",
        "hk_dot_gk",
        "trace_h_kgk_plus_trace_g_khk",
    )
    for kind in range(4):
        vector = [F(0)] * (len(PAIR_SYM) * len(MONOMIALS))
        for pair_slot, (a, b) in enumerate(PAIR_SYM):
            h, g = tensors[a], tensors[b]
            if kind == 0:
                momentum = mscale(frobenius(h, g), identity(3))
            elif kind == 1:
                momentum = mscale(trace(h) * trace(g), identity(3))
            elif kind == 2:
                momentum = sym_product(h, g)
            else:
                momentum = madd(mscale(trace(h), g), mscale(trace(g), h))
            coefficients = momentum_coefficients(momentum)
            for monomial_slot, value in enumerate(coefficients):
                vector[pair_slot * 6 + monomial_slot] = value
        bases.append(tuple(vector))
    CHECK.equal(rank(bases), 4, "SO3 self-adjoint quadratic basis dimension")
    return tuple(bases), labels


def constant_dimensions(group, tensors):
    pair_sym_bases = []
    for i, j in PAIR_SYM:
        seed = elementary_symmetric(6, i, j)
        average = zero(6, 6)
        for _, _, pair_action in group:
            average = madd(average, matmul(matmul(transpose(pair_action), seed),
                                           pair_action))
        vector = tuple(average[a][b] for a, b in PAIR_SYM)
        if rank(pair_sym_bases + [vector]) > rank(pair_sym_bases):
            pair_sym_bases.append(vector)
    CHECK.equal(len(pair_sym_bases), 3, "S4 constant pair-kernel dimension")

    rotation_constants = []
    for kind in range(2):
        matrix = [[F(0)] * 6 for _ in range(6)]
        for a in range(6):
            for b in range(6):
                matrix[a][b] = (frobenius(tensors[a], tensors[b]) if kind == 0
                                else trace(tensors[a]) * trace(tensors[b]))
        rotation_constants.append(tuple(matrix[i][j] for i, j in PAIR_SYM))
    CHECK.equal(rank(rotation_constants), 2,
                "SO3 constant symmetric-tensor kernel dimension")
    return len(pair_sym_bases), len(rotation_constants)


def representation_character_dimension(group):
    total = F(0)
    for permutation, _, _ in group:
        fixed = sum(permutation[i] == i for i in range(4))
        square = tuple(permutation[permutation[i]] for i in range(4))
        fixed_square = sum(square[i] == i for i in range(4))
        chi_v = F(fixed - 1)
        chi_v2 = F(fixed_square - 1)
        chi_sym_v = (chi_v * chi_v + chi_v2) / 2

        fixed_pairs = sum(tuple(sorted((permutation[a], permutation[b]))) == (a, b)
                          for a, b in PAIRS)
        fixed_pairs_square = sum(
            tuple(sorted((square[a], square[b]))) == (a, b) for a, b in PAIRS)
        chi_pair = F(fixed_pairs)
        chi_pair2 = F(fixed_pairs_square)
        chi_sym_pair = (chi_pair * chi_pair + chi_pair2) / 2
        total += chi_sym_v * chi_sym_pair
    answer = total / len(group)
    CHECK.equal(answer, F(9), "character count agrees with Reynolds rank")
    return answer


def reconstruct_matrix_at(vector, momentum):
    monomial_values = tuple(F(momentum[i]) * F(momentum[j])
                            for i, j in MONOMIALS)
    matrix = [[F(0)] * 6 for _ in range(6)]
    for pair_slot, (a, b) in enumerate(PAIR_SYM):
        value = sum((vector[pair_slot * 6 + m] * monomial_values[m]
                     for m in range(6)), F(0))
        matrix[a][b] = matrix[b][a] = value
    return tuple(tuple(row) for row in matrix)


def t2_block_coefficients(vector):
    """Return A,B,C for raw-T block A r2 I+B D+C O, and verify the fit."""
    t_basis = (
        (F(1), F(0), F(0), F(0), F(0), F(-1)),
        (F(0), F(1), F(0), F(0), F(-1), F(0)),
        (F(0), F(0), F(1), F(-1), F(0), F(0)),
    )
    block = [[[F(0)] * 6 for _ in range(3)] for _ in range(3)]
    for monomial_slot in range(6):
        pair_matrix = [[F(0)] * 6 for _ in range(6)]
        for pair_slot, (a, b) in enumerate(PAIR_SYM):
            pair_matrix[a][b] = pair_matrix[b][a] = vector[pair_slot * 6 + monomial_slot]
        for i in range(3):
            for j in range(3):
                block[i][j][monomial_slot] = sum((
                    t_basis[i][a] * pair_matrix[a][b] * t_basis[j][b]
                    for a in range(6) for b in range(6)), F(0))
    a_coefficient = block[0][0][1]  # ky^2 coefficient away from D_00.
    b_coefficient = block[0][0][0] - a_coefficient
    c_coefficient = block[0][1][3]  # kx ky coefficient in O_01.
    for i in range(3):
        for j in range(3):
            expected = [F(0)] * 6
            if i == j:
                for m in range(3):
                    expected[m] += a_coefficient
                expected[i] += b_coefficient
            else:
                pair = tuple(sorted((i, j)))
                expected[MONOMIALS.index(pair)] += c_coefficient
            CHECK.equal(tuple(block[i][j]), tuple(expected),
                        "cubic T2 block has A r2 I+B D+C O form")
    return a_coefficient, b_coefficient, c_coefficient


def direct_cubic_gauge_constraints(invariant, solder):
    """Exact coefficient constraints for K(k)(k odot xi)=0.

    The six-by-six solder maps pair coordinates to symmetric tensors.  For
    each Cartesian xi basis vector, solve the three coefficients of the pair
    coordinate that represent k odot xi.  Multiplication by each quadratic
    invariant basis symbol then gives a homogeneous cubic polynomial.  Its
    ten monomial coefficients, for all six outputs and three xi directions,
    are the complete direct gauge-null constraints.
    """
    solder_columns = transpose(solder)
    gauge_pair_coefficients = []
    for xi_slot in range(3):
        by_momentum = []
        for momentum_slot in range(3):
            momentum = identity(3)[momentum_slot]
            direction = identity(3)[xi_slot]
            gauge_tensor = tuple(tuple(
                momentum[i] * direction[j] + direction[i] * momentum[j]
                for j in range(3)) for i in range(3))
            gauge_coordinates = tuple(
                gauge_tensor[i][j]
                for i, j in ((0, 0), (1, 1), (2, 2),
                             (0, 1), (0, 2), (1, 2)))
            by_momentum.append(solve_columns(solder_columns,
                                             gauge_coordinates))
        gauge_pair_coefficients.append(tuple(by_momentum))

    cubic_monomials = tuple(
        (a, b, 3 - a - b)
        for a in range(4) for b in range(4 - a)
    )
    CHECK.equal(len(cubic_monomials), 10,
                "ten homogeneous cubic momentum monomials")
    cubic_index = {monomial: slot
                   for slot, monomial in enumerate(cubic_monomials)}

    rows = []
    for xi_slot in range(3):
        for output in range(6):
            coefficient_rows = [[F(0)] * len(cubic_monomials)
                                for _ in invariant]
            for basis_slot, vector in enumerate(invariant):
                for pair_slot, (left, right) in enumerate(PAIR_SYM):
                    for momentum_slot, (p, q) in enumerate(MONOMIALS):
                        coefficient = vector[pair_slot * 6 + momentum_slot]
                        if not coefficient:
                            continue
                        for input_pair in range(6):
                            if output == left and input_pair == right:
                                matrix_coefficient = coefficient
                            elif output == right and input_pair == left:
                                matrix_coefficient = coefficient
                            else:
                                continue
                            for linear_slot in range(3):
                                linear_coefficient = (
                                    gauge_pair_coefficients[xi_slot]
                                                           [linear_slot]
                                                           [input_pair]
                                )
                                if not linear_coefficient:
                                    continue
                                exponent = [0, 0, 0]
                                exponent[p] += 1
                                exponent[q] += 1
                                exponent[linear_slot] += 1
                                coefficient_rows[basis_slot][
                                    cubic_index[tuple(exponent)]
                                ] += matrix_coefficient * linear_coefficient
            for monomial_slot in range(len(cubic_monomials)):
                row = tuple(coefficient_rows[basis_slot][monomial_slot]
                            for basis_slot in range(len(invariant)))
                if any(row):
                    rows.append(row)

    constraint_rank = rank(rows)
    null = nullspace(rows)
    CHECK.equal(constraint_rank, 8,
                "direct cubic gauge-null constraint rank eight")
    CHECK.equal(len(null), 1,
                "direct cubic gauge-null subspace is one-dimensional")
    return tuple(rows), null[0], cubic_monomials


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-ledger", action="store_true")
    args = parser.parse_args()

    group = group_data()
    invariant, invariant_labels = invariant_basis(group)
    for vector in invariant:
        for _, rotation, pair_action in group:
            left, right = covariance_sides(vector, rotation, pair_action)
            CHECK.equal(left, right, "Reynolds basis satisfies exact S4 covariance")
    character_dimension = representation_character_dimension(group)
    tensors, solder = symmetric_tensor_solder()
    verify_solder_equivariance(group, tensors)
    rotational, rotational_labels = rotational_bases(tensors)
    constant_cubic, constant_rotational = constant_dimensions(group, tensors)

    coordinates = tuple(solve_columns(invariant, vector)
                        for vector in rotational)
    coordinate_matrix = transpose(coordinates)  # 9 x 4
    CHECK.equal(rank(coordinate_matrix), 4,
                "rotational subspace rank four inside cubic space")
    residual_covectors = nullspace(transpose(coordinate_matrix))
    CHECK.equal(len(residual_covectors), 5,
                "five independent cubic-to-rotational matching conditions")
    for residual in residual_covectors:
        CHECK.equal(tuple(sum((residual[i] * coordinate_matrix[i][j]
                               for i in range(9)), F(0)) for j in range(4)),
                    (F(0),) * 4, "rotational residual annihilates all four bases")

    t2_coefficients = tuple(t2_block_coefficients(vector) for vector in invariant)
    t2_mismatch = tuple(coefficients[1] + coefficients[2]
                        for coefficients in t2_coefficients)
    CHECK.true(any(t2_mismatch), "generic cubic T2 block has a nonzero mismatch")
    for coordinate in coordinates:
        CHECK.equal(sum((t2_mismatch[i] * coordinate[i]
                         for i in range(9)), F(0)), F(0),
                    "every full rotational basis satisfies the T2 B+C condition")

    # Einstein/Fierz-Pauli static symbol in the four rotational bilinears.
    # 1/2[r2 h:g-r2 trh trg-2(hk).(gk)
    #     +trh(kgk)+trg(khk)].
    einstein_coefficients = (F(1, 2), F(-1, 2), F(-1), F(1, 2))
    einstein_vector = tuple(sum((einstein_coefficients[j] * rotational[j][i]
                                  for j in range(4)), F(0))
                            for i in range(len(rotational[0])))
    einstein_matrix = reconstruct_matrix_at(einstein_vector, (2, 3, 5))
    CHECK.equal(rank(einstein_matrix), 3,
                "static Einstein reference has rank-three gauge quotient")
    gauge_constraints = (
        (F(2), F(0), F(1), F(0)),
        (F(0), F(0), F(1), F(2)),
        (F(0), F(2), F(0), F(2)),
    )
    CHECK.equal(rank(gauge_constraints), 3,
                "rotational family has three independent gauge-null constraints")
    CHECK.equal(matvec(gauge_constraints, (F(1), F(-1), F(-2), F(1))),
                (F(0),) * 3, "gauge-null ray is the Einstein coefficient ray")
    momentum = (F(2), F(3), F(5))
    solder_columns = transpose(solder)
    for direction in identity(3):
        gauge_tensor = tuple(tuple(momentum[i] * direction[j]
                                   + direction[i] * momentum[j]
                                   for j in range(3)) for i in range(3))
        gauge_coordinates = tuple(gauge_tensor[i][j]
                                  for i, j in ((0, 0), (1, 1), (2, 2),
                                               (0, 1), (0, 2), (1, 2)))
        pair_coordinates = solve_columns(solder_columns, gauge_coordinates)
        CHECK.equal(matvec(einstein_matrix, pair_coordinates), (F(0),) * 6,
                    "Einstein reference kills each longitudinal gauge mode")
    einstein_invariant_coordinates = solve_columns(invariant, einstein_vector)
    for residual in residual_covectors:
        CHECK.equal(sum((residual[i] * einstein_invariant_coordinates[i]
                         for i in range(9)), F(0)), F(0),
                    "Einstein reference obeys all five rotational conditions")

    # The shortest algebraic route does not need to impose SO(3) first.
    # Apply the longitudinal gauge null directly to the full nine-dimensional
    # cubic response family.  Exact cubic-coefficient elimination has rank
    # eight, and its sole null vector is the Einstein/Fierz--Pauli ray.
    direct_gauge_rows, direct_gauge_ray, cubic_monomials = (
        direct_cubic_gauge_constraints(invariant, solder)
    )
    for row in direct_gauge_rows:
        CHECK.equal(sum((row[i] * einstein_invariant_coordinates[i]
                         for i in range(9)), F(0)), F(0),
                    "Einstein reference obeys every direct cubic gauge constraint")
    ratio = next(direct_gauge_ray[i] / einstein_invariant_coordinates[i]
                 for i in range(9) if einstein_invariant_coordinates[i])
    CHECK.equal(direct_gauge_ray,
                tuple(ratio * value for value in einstein_invariant_coordinates),
                "direct cubic gauge-null ray equals Einstein ray")

    # The full SO3 family restricts to a two-dimensional T2-T2 block.  This
    # recovers the dimensional content of GL6CO while retaining E2 mixing.
    t_basis = (
        (F(1), F(0), F(0), F(0), F(0), F(-1)),
        (F(0), F(1), F(0), F(0), F(-1), F(0)),
        (F(0), F(0), F(1), F(-1), F(0), F(0)),
    )
    restricted = []
    for vector in rotational:
        samples = []
        for momentum in ((1, 0, 0), (0, 1, 0), (0, 0, 1),
                         (1, 1, 0), (1, 0, 1), (0, 1, 1)):
            matrix = reconstruct_matrix_at(vector, momentum)
            samples.extend(sum((t_basis[i][a] * matrix[a][b] * t_basis[j][b]
                                for a in range(6) for b in range(6)), F(0))
                               for i in range(3) for j in range(i, 3))
        restricted.append(tuple(samples))
    CHECK.equal(rank(restricted), 2,
                "full rotational family has two-dimensional T2 restriction")

    ledger = {
        "lane": "GL6CR",
        "scope": "full six-pair algebraic symmetry classification only; no phase, continuum, Ricci, gravity, or G",
        "pair_order": PAIRS,
        "metric_solder": {
            "definition": "D_C(e_ab)=v_a odot v_b with v_a=T_a/2",
            "coordinate_order": ("xx", "yy", "zz", "xy", "xz", "yz"),
            "matrix": solder,
            "rank": 6,
        },
        "constant_symbols": {
            "S4_dimension": constant_cubic,
            "SO3_dimension": constant_rotational,
            "consequence": "one independent constant splitting distinguishes E2 and T2 traceless response; rotational completion removes it",
        },
        "quadratic_symbols": {
            "S4_dimension_character": character_dimension,
            "S4_dimension_reynolds": len(invariant),
            "SO3_self_adjoint_dimension": len(rotational),
            "matching_codimension": len(residual_covectors),
            "invariant_basis_labels": invariant_labels,
            "rotational_basis_labels": rotational_labels,
            "rotational_coordinates_in_invariant_basis": coordinates,
            "five_residual_covectors_on_invariant_coordinates": residual_covectors,
            "T2_block_ABC_by_invariant_basis": t2_coefficients,
            "T2_extension_mismatch_B_plus_C": t2_mismatch,
            "T2_guard": "B+C=0 is only one projection of the five full-six-channel conditions",
            "necessary_and_sufficient_test": "a full S4 quadratic symbol is SO3-extendible iff all five displayed residual covectors vanish",
        },
        "general_SO3_bilinear": (
            "a |k|^2 h:g + b |k|^2 tr(h)tr(g) + "
            "c (h k).(g k) + d[tr(h)(k g k)+tr(g)(k h k)]"
        ),
        "T2_projection": {
            "dimension": rank(restricted),
            "consequence": "a T2-only test sees only two of the four rotational coefficients and cannot certify the full tensor operator",
        },
        "einstein_reference": {
            "coordinates_in_SO3_basis": einstein_coefficients,
            "coordinates_in_invariant_basis": einstein_invariant_coordinates,
            "generic_momentum": (2, 3, 5),
            "rank": rank(einstein_matrix),
            "gauge_constraints_on_general_SO3_coefficients_abcd": gauge_constraints,
            "unique_gauge_null_ray": (1, -1, -2, 1),
            "guard": "held-out algebraic 1PI reference only; the microscopic object currently derived is a connected response and must be completed and lawfully inverted",
        },
        "direct_cubic_gauge_shortcut": {
            "condition": "K^(2)(k) x(k,xi)=0 for every k and xi, where D_C x=k odot xi",
            "homogeneous_cubic_monomials": cubic_monomials,
            "constraint_rank_on_S4_nine_space": rank(direct_gauge_rows),
            "null_dimension": 1,
            "unique_null_ray_in_invariant_basis": direct_gauge_ray,
            "equals_Einstein_ray": True,
            "consequence": "on the full cubic S4 family, the direct longitudinal Ward/gauge null implies SO3 rotational completion and the Einstein/Fierz-Pauli ray; orientation averaging is not an independent algebraic prerequisite if this null is derived physically",
            "guard": "the gauge null is a target condition here, not yet derived for the complete same-state F3 1PI kernel",
        },
    }

    def encode(value):
        if isinstance(value, F):
            return (str(value.numerator) if value.denominator == 1
                    else f"{value.numerator}/{value.denominator}")
        if isinstance(value, dict):
            return {str(key): encode(item) for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            return [encode(item) for item in value]
        return value

    payload = json.dumps(encode(ledger), indent=2, sort_keys=True) + "\n"
    target = HERE / "EXACT_LEDGER.json"
    if args.write_ledger:
        target.write_text(payload)
    else:
        CHECK.true(target.is_file(), "frozen ledger exists")
        CHECK.equal(target.read_text(), payload, "frozen ledger matches exact replay")

    print(f"PASS__GL6CR_FULL_PAIR_ROTATIONAL_COMPLETION__{CHECK.total}/{CHECK.total}")
    print("CONSTANT_DIMENSIONS=S4_3;SO3_2;ONE_SPLITTING_CONDITION")
    print("QUADRATIC_DIMENSIONS=S4_9;SO3_4;FIVE_MATCHING_CONDITIONS")
    print("T2_RESTRICTION_DIMENSION=2;FULL_E2_T2_COMPLETION_MANDATORY")
    print("EINSTEIN_REFERENCE=ONE_RAY_INSIDE_SO3_FAMILY;RESPONSE_TO_1PI_OPEN")
    print("DIRECT_CUBIC_GAUGE_NULL=RANK8;UNIQUE_EINSTEIN_RAY;SO3_FOLLOWS")


if __name__ == "__main__":
    main()
