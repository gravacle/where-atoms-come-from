#!/usr/bin/env python3
"""Exact GL6CO invariant cycle-response and T-block matching replay.

Classifies every real reciprocal translation/inversion/S4-invariant
four-orientation cycle symbol through quadratic coordinate momentum, pulls it
back through the complete GL6CL T2 writer, and derives the exact condition
for its T-T block to admit an SO(3)-covariant symmetric-tensor completion.
The GL6BV h2 contact is reconstructed as a separately typed same-state block.
All checks use exact rational arithmetic and the Python standard library.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from itertools import combinations, permutations, product
from pathlib import Path


HERE = Path(__file__).resolve().parent
TETRA = (
    (F(1), F(1), F(1)),
    (F(1), F(-1), F(-1)),
    (F(-1), F(1), F(-1)),
    (F(-1), F(-1), F(1)),
)
PAIR_ORDER = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIR_ORDER)}
T_PAIR_BASIS = (
    tuple(map(F, (1, 0, 0, 0, 0, -1))),
    tuple(map(F, (0, 1, 0, 0, -1, 0))),
    tuple(map(F, (0, 0, 1, -1, 0, 0))),
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


def dot(left, right):
    return sum((F(x) * F(y) for x, y in zip(left, right)), F(0))


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
    return tuple(dot(row, vector) for row in matrix)


def madd(left, right, factor=F(1)):
    return tuple(tuple(left[i][j] + F(factor) * right[i][j]
                       for j in range(len(left[0])))
                 for i in range(len(left)))


def mscale(factor, matrix):
    return tuple(tuple(F(factor) * value for value in row) for row in matrix)


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


def inverse(matrix):
    size = len(matrix)
    work = [list(map(F, row)) + list(identity(size)[i])
            for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next(row for row in range(column, size) if work[row][column])
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


# Polynomial helpers.  A polynomial is {exponent triple: rational value}.
def pclean(poly):
    return {key: value for key, value in poly.items() if value}


def padd(left, right):
    answer = dict(left)
    for key, value in right.items():
        answer[key] = answer.get(key, F(0)) + value
    return pclean(answer)


def pscale(factor, poly):
    return pclean({key: F(factor) * value for key, value in poly.items()})


def pmul(left, right, maximum_degree=2):
    answer = {}
    for left_exp, left_value in left.items():
        for right_exp, right_value in right.items():
            exponent = tuple(left_exp[i] + right_exp[i] for i in range(3))
            if sum(exponent) <= maximum_degree:
                answer[exponent] = answer.get(exponent, F(0)) + left_value * right_value
    return pclean(answer)


ONE = {(0, 0, 0): F(1)}
ZERO = {}
KVAR = (
    {(1, 0, 0): F(1)},
    {(0, 1, 0): F(1)},
    {(0, 0, 1): F(1)},
)
R2 = padd(padd(pmul(KVAR[0], KVAR[0]), pmul(KVAR[1], KVAR[1])),
           pmul(KVAR[2], KVAR[2]))
D_POLY = tuple(tuple(pmul(KVAR[i], KVAR[i]) if i == j else ZERO
                     for j in range(3)) for i in range(3))
O_POLY = tuple(tuple(ZERO if i == j else pmul(KVAR[i], KVAR[j])
                     for j in range(3)) for i in range(3))
Q_POLY = (pmul(KVAR[1], KVAR[2]),
          pmul(KVAR[2], KVAR[0]),
          pmul(KVAR[0], KVAR[1]))


def pmat_zero(rows, columns):
    return tuple(tuple(ZERO for _ in range(columns)) for _ in range(rows))


def pmat_add(left, right, factor=F(1)):
    return tuple(tuple(padd(left[i][j], pscale(factor, right[i][j]))
                       for j in range(len(left[0])))
                 for i in range(len(left)))


def pmat_scale(factor, matrix):
    return tuple(tuple(pscale(factor, value) for value in row) for row in matrix)


def pmat_transpose(matrix):
    return tuple(tuple(matrix[i][j] for i in range(len(matrix)))
                 for j in range(len(matrix[0])))


def pmat_mul(left, right, maximum_degree=2):
    return tuple(tuple(
        sum_polynomials(pmul(left[i][q], right[q][j], maximum_degree)
                        for q in range(len(right)))
        for j in range(len(right[0]))) for i in range(len(left)))


def sum_polynomials(polynomials):
    answer = {}
    for polynomial in polynomials:
        answer = padd(answer, polynomial)
    return answer


def constant_pmatrix(matrix):
    return tuple(tuple(pscale(value, ONE) for value in row) for row in matrix)


def poly_outer(left, right):
    return tuple(tuple(pmul(x, y) for y in right) for x in left)


def transformed_variables(matrix):
    return tuple(sum_polynomials(pscale(matrix[i][j], KVAR[j])
                                 for j in range(3)) for i in range(3))


def polynomial_matrix_at_transformed_k(kind, orthogonal):
    kp = transformed_variables(orthogonal)
    r2 = sum_polynomials(pmul(x, x) for x in kp)
    diagonal = tuple(tuple(pmul(kp[i], kp[i]) if i == j else ZERO
                           for j in range(3)) for i in range(3))
    off = tuple(tuple(ZERO if i == j else pmul(kp[i], kp[j])
                      for j in range(3)) for i in range(3))
    q = (pmul(kp[1], kp[2]), pmul(kp[2], kp[0]), pmul(kp[0], kp[1]))
    if kind == "AA":
        out = [list(row) for row in pmat_zero(4, 4)]
        out[0][0] = r2
        return tuple(tuple(row) for row in out)
    if kind == "AV":
        out = [list(row) for row in pmat_zero(4, 4)]
        for i in range(3):
            out[0][i + 1] = q[i]
            out[i + 1][0] = q[i]
        return tuple(tuple(row) for row in out)
    if kind == "VI":
        out = [list(row) for row in pmat_zero(4, 4)]
        for i in range(3):
            out[i + 1][i + 1] = r2
        return tuple(tuple(row) for row in out)
    if kind == "VD":
        out = [list(row) for row in pmat_zero(4, 4)]
        for i in range(3):
            for j in range(3):
                out[i + 1][j + 1] = diagonal[i][j]
        return tuple(tuple(row) for row in out)
    if kind == "VO":
        out = [list(row) for row in pmat_zero(4, 4)]
        for i in range(3):
            for j in range(3):
                out[i + 1][j + 1] = off[i][j]
        return tuple(tuple(row) for row in out)
    raise ValueError(kind)


def basis_symbol(kind):
    return polynomial_matrix_at_transformed_k(kind, identity(3))


def conjugate_polymatrix(matrix, orthogonal4):
    left = constant_pmatrix(orthogonal4)
    right = constant_pmatrix(transpose(orthogonal4))
    return pmat_mul(pmat_mul(left, matrix), right)


def permutation_matrix(permutation):
    size = len(permutation)
    matrix = [[F(0)] * size for _ in range(size)]
    for source, target in enumerate(permutation):
        matrix[target][source] = F(1)
    return tuple(tuple(row) for row in matrix)


def group_classification():
    tetra_columns = transpose(TETRA[:3])
    tetra_inverse = inverse(tetra_columns)
    u = (F(1, 2),) * 4
    q_cycle = tuple(tuple(TETRA[d][i] / 2 for i in range(3))
                    for d in range(4))
    solder = tuple(tuple([u[row]] + list(q_cycle[row])) for row in range(4))
    CHECK.equal(matmul(transpose(solder), solder), identity(4),
                "cycle A+V solder is orthogonal")

    group = []
    character_sum = F(0)
    constant_invariant_sum = F(0)
    basis_kinds = ("AA", "AV", "VI", "VD", "VO")
    for perm in permutations(range(4)):
        target_columns = transpose(tuple(TETRA[perm[a]] for a in range(3)))
        orthogonal = matmul(target_columns, tetra_inverse)
        CHECK.equal(matmul(transpose(orthogonal), orthogonal), identity(3),
                    "tetrahedral port action is orthogonal")
        CHECK.true(all(sum(value != 0 for value in row) == 1 for row in orthogonal),
                   "tetrahedral action is signed-permutation in Cartesian chart")
        pmat = permutation_matrix(perm)
        transformed_cycle = matmul(matmul(transpose(solder), pmat), solder)
        expected_cycle = tuple(
            tuple(F(1) if i == j == 0 else
                  (orthogonal[i - 1][j - 1] if i and j else F(0))
                  for j in range(4)) for i in range(4))
        CHECK.equal(transformed_cycle, expected_cycle,
                    "cycle permutation decomposes as A plus Cartesian V")

        # Exact polynomial covariance of all five proposed quadratic bases.
        for kind in basis_kinds:
            left = polynomial_matrix_at_transformed_k(kind, orthogonal)
            right = conjugate_polymatrix(basis_symbol(kind), expected_cycle)
            CHECK.equal(left, right, f"quadratic invariant basis covariance {kind}")

        # Character count of Hom_S4(Sym2(V),Sym2(A+V)).
        fixed = sum(perm[i] == i for i in range(4))
        perm2 = tuple(perm[perm[i]] for i in range(4))
        fixed2 = sum(perm2[i] == i for i in range(4))
        chi_r, chi_r2 = F(fixed), F(fixed2)
        chi_v, chi_v2 = chi_r - 1, chi_r2 - 1
        chi_sym_r = (chi_r * chi_r + chi_r2) / 2
        chi_sym_v = (chi_v * chi_v + chi_v2) / 2
        character_sum += chi_sym_r * chi_sym_v
        constant_invariant_sum += chi_sym_r
        group.append((perm, orthogonal, pmat))
    quadratic_dimension = character_sum / 24
    constant_dimension = constant_invariant_sum / 24
    CHECK.equal(quadratic_dimension, F(5), "five quadratic invariant coefficients")
    CHECK.equal(constant_dimension, F(2), "two constant symmetric invariants")
    return {
        "group_order": len(group),
        "cycle_solder_A_plus_V": solder,
        "quadratic_invariant_dimension": quadratic_dimension,
        "constant_invariant_dimension": constant_dimension,
        "quadratic_basis": (
            "A-A: alpha |k|^2",
            "A-V: eta (ky kz,kz kx,kx ky)",
            "V-V: b |k|^2 I",
            "V-V: c diag(kx^2,ky^2,kz^2)",
            "V-V: d (k k^T-diag(kx^2,ky^2,kz^2))",
        ),
    }, solder


def vadd(*vectors):
    return tuple(sum((F(v[i]) for v in vectors), F(0)) for i in range(len(vectors[0])))


def vscale(factor, vector):
    return tuple(F(factor) * value for value in vector)


def cycle_rhos():
    orientations = []
    for missing in range(4):
        a, b, c = tuple(port for port in range(4) if port != missing)
        orientations.append({
            tuple(sorted((a, b))): vscale(F(1, 2), vadd(TETRA[a], TETRA[b], vscale(-1, TETRA[c]))),
            tuple(sorted((a, c))): vscale(F(1, 2), vadd(TETRA[a], vscale(-1, TETRA[b]), TETRA[c])),
            tuple(sorted((b, c))): vscale(F(1, 2), vadd(vscale(-1, TETRA[a]), TETRA[b], TETRA[c])),
        })
    return tuple(orientations)


def cosine_degree_two(rho):
    linear = sum_polynomials(pscale(rho[i], KVAR[i]) for i in range(3))
    return padd(ONE, pscale(F(-1, 2), pmul(linear, linear)))


def writer_matrix_rational_t():
    rows = []
    for orientation in cycle_rhos():
        components = [ZERO for _ in range(6)]
        for pair, rho in orientation.items():
            components[PAIR_INDEX[pair]] = pscale(F(2), cosine_degree_two(rho))
        row = []
        for basis in T_PAIR_BASIS:
            row.append(sum_polynomials(pscale(coefficient, polynomial)
                                       for coefficient, polynomial in zip(basis, components)))
        rows.append(tuple(row))
    return tuple(rows)


def av_to_cycle(symbol_av, solder):
    return pmat_mul(pmat_mul(constant_pmatrix(solder), symbol_av),
                    constant_pmatrix(transpose(solder)))


def pullback(writer, cycle_symbol):
    answer = pmat_mul(pmat_mul(pmat_transpose(writer), cycle_symbol), writer)
    # Rational T-pair basis vectors have norm^2=2.  Divide the coordinate
    # Gram by two to express it in the orthonormal basis t_i/sqrt(2).
    return pmat_scale(F(1, 2), answer)


def expected_matrix(kind):
    out = [[ZERO for _ in range(3)] for _ in range(3)]
    if kind == "KAPPA":
        for i in range(3):
            out[i][i] = padd(pscale(F(8), ONE),
                             padd(pscale(F(-2), R2),
                                  pscale(F(-28), D_POLY[i][i])))
        kk = pmat_add(D_POLY, O_POLY)
        out = [list(row) for row in tuple(tuple(row) for row in out)]
        for i in range(3):
            for j in range(3):
                out[i][j] = padd(out[i][j], pscale(F(12), kk[i][j]))
        return tuple(tuple(row) for row in out)
    if kind == "B":
        for i in range(3):
            out[i][i] = pscale(F(8), R2)
    elif kind == "C":
        for i in range(3):
            out[i][i] = pscale(F(8), D_POLY[i][i])
    elif kind == "D":
        out = [list(row) for row in pmat_scale(F(8), O_POLY)]
    return tuple(tuple(row) for row in out)


def composition(solder):
    writer = writer_matrix_rational_t()
    pv_av = [[ZERO for _ in range(4)] for _ in range(4)]
    for i in range(1, 4):
        pv_av[i][i] = ONE
    pv_av = tuple(tuple(row) for row in pv_av)

    symbols = {
        "KAPPA": pv_av,
        "ALPHA": basis_symbol("AA"),
        "ETA": basis_symbol("AV"),
        "B": basis_symbol("VI"),
        "C": basis_symbol("VD"),
        "D": basis_symbol("VO"),
    }
    pullbacks = {}
    for name, symbol in symbols.items():
        pulled = pullback(writer, av_to_cycle(symbol, solder))
        pullbacks[name] = pulled
        if name in ("ALPHA", "ETA"):
            CHECK.equal(pulled, pmat_zero(3, 3),
                        f"{name} cycle block first enters T pullback above k2")
        else:
            CHECK.equal(pulled, expected_matrix(name),
                        f"exact composed coefficient for {name}")

    # Coefficient vector in basis |k|^2 I, D, O.
    # H_T=8 kappa I + A |k|^2 I+B D+C O+O(k4).
    formula = {
        "A": "-2 kappa+8 b",
        "B": "-16 kappa+8 c",
        "C": "12 kappa+8 d",
        "cubic_mismatch_B_plus_C": "-4 kappa+8(c+d)",
    }

    # The T-T restriction of a general SO(3)-covariant symmetric-tensor
    # quadratic has basis r2 I and (r2 I-D)+O.  In (A,B,C) coordinates these
    # are (1,0,0) and (1,-1,1), a rank-two plane with B+C=0.
    so3_columns = ((F(1), F(1)), (F(0), F(-1)), (F(0), F(1)))
    CHECK.equal(rank(so3_columns), 2, "SO3-extendible T-block plane rank two")
    left_null = (F(0), F(1), F(1))
    CHECK.equal(tuple(sum(left_null[row] * so3_columns[row][column]
                          for row in range(3)) for column in range(2)),
                (F(0), F(0)), "SO3 T-block plane has left null B+C")
    # Direct substitution proves the matching condition and the stronger
    # Fierz-Pauli/Einstein-reference algebraic shape.
    CHECK.equal(F(-4) + F(8) * (F(1, 4) + F(1, 4)), F(0),
                "positive witness satisfies one SO3 relation")
    witness_a = F(-2) + F(8) * F(1, 4)
    witness_b = F(-16) + F(8) * F(1, 4)
    witness_c = F(12) + F(8) * F(1, 4)
    CHECK.equal((witness_a, witness_b, witness_c), (F(0), F(-14), F(14)),
                "positive witness gives algebraic D-minus-O reference ray")
    return {
        "writer_unscaled_rational_T_basis": writer,
        "basis_pullbacks_orthonormal_T": pullbacks,
        "composed_coefficients": formula,
        "SO3_TT_general_form": "u |k|^2 I+v[(|k|^2 I-D)+O]",
        "SO3_extension_iff": "B+C=0 iff c+d=kappa/2",
        "condition_count": "one independent linear condition; not automatic and not impossible",
        "bare_constant_cycle_response": "b=c=d=0 gives B+C=-4 kappa, so fails for nonzero kappa",
        "positive_matching_witness": "kappa=1, alpha=1, eta=0, b=c=d=1/4; K_VV=I+(|k|^2 I+k k^T)/4 is positive near zero",
        "reference_FP_TT_shape": "2D-k k^T=D-O",
        "reference_FP_additional_condition": "A=0 iff b=kappa/4, in addition to c+d=kappa/2",
        "reference_warning": "algebraic T-block proportionality is not a Ricci proof; full E2-T2/E2-E2 completion, contact, and response-to-1PI inversion remain required",
    }


def contact_block():
    tau = []
    for exceptional in range(4):
        tau.append(tuple(F(-1 if exceptional in pair else 1)
                         for pair in PAIR_ORDER))
    # Coordinates of Q_a=tau_a tau_a^T/6 in orthonormal t_i/sqrt(2).
    q_matrices = []
    for vector in tau:
        coordinates = tuple(dot(basis, vector) for basis in T_PAIR_BASIS)
        q_matrices.append(tuple(tuple(coordinates[i] * coordinates[j] / 12
                                      for j in range(3)) for i in range(3)))
    CHECK.true(all(sum(q[i][i] for i in range(3)) == 1 for q in q_matrices),
               "each contact Q_a is rank-one normalized")

    # theta_a=k.T_a; form sum theta_a^2 Q_a.
    weighted = pmat_zero(3, 3)
    i2 = {}
    for a in range(4):
        theta = sum_polynomials(pscale(TETRA[a][i], KVAR[i]) for i in range(3))
        theta2 = pmul(theta, theta)
        i2 = padd(i2, theta2)
        contribution = tuple(tuple(pscale(q_matrices[a][i][j], theta2)
                                   for j in range(3)) for i in range(3))
        weighted = pmat_add(weighted, contribution)
    CHECK.equal(i2, pscale(F(4), R2), "tetrahedral character norm I2=4|k|2")
    expected_weighted = pmat_add(
        tuple(tuple(pscale(F(4, 3), R2) if i == j else ZERO
                    for j in range(3)) for i in range(3)),
        pmat_scale(F(8, 3), O_POLY))
    CHECK.equal(weighted, expected_weighted,
                "exact contact sum theta_a^2 Q_a")

    # zeta=4p-4/3.  Coefficients in A|k|2 I+B D+C O.
    contact = {
        "state_parameter": "p=Pr[sigma_v(a)=sigma_w(a)]",
        "zeta": "4p-4/3",
        "gamma": "(8/3)(4p-1)",
        "constant_common_T": "8+gamma=(16/3)(1+2p)",
        "quadratic_A": "-(4/3)zeta-4/9=(4/3)(1-4p)",
        "quadratic_B": "0",
        "quadratic_C": "(4/3)zeta-8/9=(8/3)(2p-1)",
        "SO3_mismatch_B_plus_C": "(8/3)(2p-1)",
        "contact_alone_SO3_condition": "p=1/2",
        "overall_scale": "g_ct=h^2/(4U_d^3) in the connected-functional sign convention",
        "same_state_guard": "p must be evaluated in the same stationary state as the H6 cycle response before blocks may be added",
    }
    CHECK.equal(F(4, 3) * F(199, 96) - F(8, 9), F(15, 8),
                "Q4 orbit witness contact mismatch is 15/8")
    contact["Q4_orbit_witness_only"] = (
        "p=109/128 gives zeta=199/96 and mismatch 15/8; this orbit mixture "
        "is not proved stationary under the H6 Hamiltonian")
    return contact, q_matrices, weighted


def full_tensor_reference_kill_test():
    """Reconstruct one full-solder linearized-Einstein reference matrix.

    This is a held-out algebraic diagnostic only.  It proves that matching a
    T-T subblock cannot replace the missing A/E and mixed blocks.
    """
    v = tuple(tuple(entry / 2 for entry in vector) for vector in TETRA)

    def solder(pair_coordinates):
        return tuple(tuple(sum((
            F(pair_coordinates[index]) *
            (v[a][i] * v[b][j] + v[b][i] * v[a][j])
            for index, (a, b) in enumerate(PAIR_ORDER)), F(0))
            for j in range(3)) for i in range(3))

    pair_basis = (
        (1, 1, 1, 1, 1, 1),
        (1, -1, 0, 0, -1, 1),
        (1, 1, -2, -2, 1, 1),
        (1, 0, 0, 0, 0, -1),
        (0, 1, 0, 0, -1, 0),
        (0, 0, 1, -1, 0, 0),
    )
    tensors = tuple(solder(coordinate) for coordinate in pair_basis)
    CHECK.equal(tensors[0], mscale(F(-1), identity(3)),
                "inherited solder maps A to minus identity")
    expected_tensors = (
        ((0, 0, 0), (0, 0, -1), (0, -1, 0)),
        ((0, 0, -1), (0, 0, 0), (-1, 0, 0)),
        ((0, -1, 0), (-1, 0, 0), (0, 0, 0)),
    )
    CHECK.equal(tensors[3:], tuple(tuple(tuple(map(F, row)) for row in matrix)
                                    for matrix in expected_tensors),
                "inherited solder maps T to minus offdiagonal tensors")

    momentum = (F(2), F(3), F(5))
    momentum_squared = dot(momentum, momentum)

    def einstein_symbol(tensor):
        trace = sum((tensor[i][i] for i in range(3)), F(0))
        kh = tuple(sum((momentum[q] * tensor[q][j] for q in range(3)), F(0))
                   for j in range(3))
        khk = dot(momentum, kh)
        return tuple(tuple(F(1, 2) * (
            momentum_squared * tensor[i][j] + momentum[i] * momentum[j] * trace
            - momentum[i] * kh[j] - momentum[j] * kh[i]
            - F(i == j) * (momentum_squared * trace - khk))
            for j in range(3)) for i in range(3))

    images = tuple(einstein_symbol(tensor) for tensor in tensors)
    bilinear = tuple(tuple(sum((tensors[i][a][b] * images[j][a][b]
                                for a in range(3) for b in range(3)), F(0))
                           for j in range(6)) for i in range(6))
    expected = tuple(tuple(map(F, row)) for row in (
        (-38, 5, 37, 15, 10, 6),
        (5, 100, 20, -30, 20, 0),
        (37, 20, 4, -30, -20, 24),
        (15, -30, -30, 4, -6, -10),
        (10, 20, -20, -6, 9, -15),
        (6, 0, 24, -10, -15, 25),
    ))
    CHECK.equal(bilinear, expected, "full six-direction Einstein reference matrix")
    CHECK.equal(rank(bilinear), 3, "full Einstein reference has rank-three gauge quotient")
    d_minus_o = tuple(tuple(
        momentum[i] ** 2 if i == j else -momentum[i] * momentum[j]
        for j in range(3)) for i in range(3))
    CHECK.equal(tuple(tuple(row[j] for j in range(3, 6))
                      for row in bilinear[3:6]), d_minus_o,
                "T-T reference block is D minus O")
    CHECK.true(any(bilinear[i][j] for i in range(3) for j in range(3, 6)),
               "full reference has indispensable A/E to T cross blocks")
    return {
        "pair_basis": ("A", "E1", "E2", "T1", "T2", "T3"),
        "orthonormal_T_solder": (
            "if x are coordinates on T_i/sqrt(2), then "
            "(h_yz,h_zx,h_xy)=(-1/sqrt(2))(x1,x2,x3); the common congruence "
            "factor does not change matching relations"
        ),
        "test_momentum": momentum,
        "bilinear_matrix": bilinear,
        "rank": 3,
        "T_T_block": d_minus_o,
        "disposition": "held-out full-solder kill test only; cross blocks and rank-three quotient are mandatory for any later Ricci claim",
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

    classification, solder = group_classification()
    composed = composition(solder)
    contact, contact_q, contact_weighted = contact_block()
    full_reference = full_tensor_reference_kill_test()

    ledger = {
        "lane": "GL6CO",
        "classification": classification,
        "general_centered_cycle_symbol": {
            "A_plus_V_basis": "S^T K_cycle(k) S",
            "normalization": "K_cycle=K_bare with K_bare_cc'=2 Re <0|T_c R T_c'|0>; it contains no writer factor",
            "relation_to_GL6CM": "GL6CM CM10 equals lambda_T^2 K_bare; the actual common pair-source-to-cycle-amplitude derivative here is mu B_T",
            "units_and_scale": "kappa,alpha,eta,b,c,d have inverse-energy units and scale as O(1/J) in the isolated spectral regime",
            "formula": "[[alpha |k|^2, eta q(k)^T],[eta q(k), kappa I+b|k|^2 I+cD+dO]]+O(|k|^4)",
            "q": "(ky kz,kz kx,kx ky)",
            "D": "diag(kx^2,ky^2,kz^2)",
            "O": "k k^T-D",
            "stationary_common_null": "K_cycle(0)u=0 removes the constant A1 coefficient",
            "positivity": "kappa>=0 and, for a nonnegative analytic A branch, alpha>=0; higher-order Schur constraints not classified",
            "constant_tensor_guard": "nonzero kappa is local susceptibility/background-potential data; no zero-derivative, background-stationarity, masslessness, or Einstein-Hilbert condition is proved",
        },
        "writer_pullback": composed,
        "common_source_normalization": {
            "GL6CL_coordinates": "j_P=j_plus+j_minus; j_C=j_plus-j_minus, so j_plus=(j_P+j_C)/2",
            "orthonormal_common_coordinate": "jhat_plus=(j_P+j_C)/sqrt(2)=sqrt(2)j_plus",
            "normalized_writer": "delta a=(mu/sqrt(2))B_T jhat_plus",
            "normalized_cycle_hessian": "Hhat=(mu^2/2)B_T^*K_cycle B_T",
        },
        "h2_contact_separate_block": {
            **contact,
            "Q_a_orthonormal_T_coordinates": contact_q,
            "sum_theta_squared_Q_a": contact_weighted,
        },
        "conditional_same_state_total": {
            "mu": "(105/8)h^6/U_d^6",
            "SO3_extension_condition": "(mu^2/2)[-4 kappa+8(c+d)]+[h^2/(4U_d^3)](8/3)(2p-1)=0",
            "reference_FP_additional_condition": "(mu^2/2)[-2 kappa+8b]+[h^2/(4U_d^3)](4/3)(1-4p)=0",
            "guard": "displayed sum is lawful only after one stationary state and one complete source-first functional own both blocks; no cancellation or phase is asserted",
        },
        "strong_lock_power_counting": {
            "J": "(63/8)h^6/U_d^5",
            "spectral_coefficients": "kappa,b,c,d=O(1/J)",
            "cycle_pullback": "mu^2/J=O(h^6/U_d^7)",
            "contact": "h^2/(4U_d^3)=O(h^2/U_d^3)",
            "consequence": "isolated strong-lock cycle repetition cannot generically cancel a nonzero leading contact mismatch order by order; a match requires finite-ratio/collective stationary response, an independently vanishing leading relation, or another same-order block",
            "guard": "power counting only; no phase transition is claimed",
        },
        "full_tensor_reference_kill_test": full_reference,
        "ceiling": "exact symmetry/matching theorem only; coefficients kappa,b,c,d and same-state p are not derived, GL6CN diagonal h6 source is not included, E2 completion and 1PI inversion are absent, and no phase/spacetime/metric/Ricci/Einstein/gravity/G claim is made",
    }
    payload = json.dumps(encode(ledger), indent=2, sort_keys=True) + "\n"
    target = HERE / "EXACT_LEDGER.json"
    if args.write_ledger:
        target.write_text(payload)
    else:
        CHECK.true(target.is_file(), "frozen exact ledger exists")
        CHECK.equal(target.read_text(), payload, "frozen exact ledger matches replay")

    print(f"PASS__GL6CO_CYCLE_RESPONSE_TENSOR_MATCHING__{CHECK.total}/{CHECK.total}")
    print("CYCLE_SYMBOL=2_CONSTANT_INVARIANTS;5_QUADRATIC_INVARIANTS;STATIONARY_A1_NULL")
    print("T_PULLBACK=8KAPPA_I+A_R2_I+B_D+C_O;A=-2KAPPA+8B;B=-16KAPPA+8C;C=12KAPPA+8D")
    print("SO3_TENSOR_EXTENSION=ONE_CONDITION_C_PLUS_D_EQ_KAPPA_OVER_2")
    print("REFERENCE_FP_SHAPE=SECOND_CONDITION_B_EQ_KAPPA_OVER_4;NOT_RICCI_PROOF")
    print("H2_CONTACT=SEPARATE_SAME_STATE_BLOCK;NO_CANCELLATION_OR_PHASE_ASSUMED")
    print("NO_SPACETIME_METRIC_RICCI_EINSTEIN_GRAVITY_G")


if __name__ == "__main__":
    main()
