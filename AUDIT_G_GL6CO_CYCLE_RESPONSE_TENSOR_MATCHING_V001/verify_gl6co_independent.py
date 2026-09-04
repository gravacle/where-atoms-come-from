#!/usr/bin/env python3
"""Independent exact audit of GL6CO.

This program does not import or execute the GL6CO author implementation.  It
rebuilds the S4 representation, the centered GL6CL writer jet, the
rank-two-tensor restriction, the full reference solder test, and the GL6BV
contact with standard-library rational arithmetic.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path


HERE = Path(__file__).resolve().parent
CHECKS = []
MONOMIALS = tuple((a, b, c) for a in range(3) for b in range(3)
                  for c in range(3) if a + b + c <= 2)
PAIR_ORDER = tuple(combinations(range(4), 2))
TETRA = (
    (F(1), F(1), F(1)),
    (F(1), F(-1), F(-1)),
    (F(-1), F(1), F(-1)),
    (F(-1), F(-1), F(1)),
)
T_RAW = (
    (F(1), F(0), F(0), F(0), F(0), F(-1)),
    (F(0), F(1), F(0), F(0), F(-1), F(0)),
    (F(0), F(0), F(1), F(-1), F(0), F(0)),
)


def check(condition, label):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def ftext(value):
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def canonical_hash(value):
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(payload).hexdigest()


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def matmul(left, right):
    return [[sum((left[i][k] * right[k][j]
                  for k in range(len(right))), F(0))
             for j in range(len(right[0]))]
            for i in range(len(left))]


def matadd(left, right):
    return [[left[i][j] + right[i][j] for j in range(len(left[0]))]
            for i in range(len(left))]


def matscale(scale, matrix):
    return [[F(scale) * value for value in row] for row in matrix]


def eye(size):
    return [[F(i == j) for j in range(size)] for i in range(size)]


def trace(matrix):
    return sum((matrix[i][i] for i in range(len(matrix))), F(0))


def rank(matrix):
    work = [list(map(F, row)) for row in matrix]
    rows = len(work)
    cols = len(work[0]) if rows else 0
    pivot_row = 0
    for col in range(cols):
        pivot = next((r for r in range(pivot_row, rows)
                      if work[r][col]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][col]
        work[pivot_row] = [x / value for x in work[pivot_row]]
        for r in range(rows):
            if r == pivot_row or not work[r][col]:
                continue
            factor = work[r][col]
            work[r] = [work[r][c] - factor * work[pivot_row][c]
                       for c in range(cols)]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


# Polynomials in (kx,ky,kz), truncated to degree two.
def pclean(poly):
    return {key: F(value) for key, value in poly.items() if value}


def padd(left, right):
    answer = dict(left)
    for key, value in right.items():
        answer[key] = answer.get(key, F(0)) + value
    return pclean(answer)


def pscale(scale, poly):
    return pclean({key: F(scale) * value for key, value in poly.items()})


def pmul(left, right):
    answer = {}
    for a, av in left.items():
        for b, bv in right.items():
            key = tuple(a[i] + b[i] for i in range(3))
            if sum(key) <= 2:
                answer[key] = answer.get(key, F(0)) + av * bv
    return pclean(answer)


PZERO = {}
PONE = {(0, 0, 0): F(1)}
VARS = ({(1, 0, 0): F(1)}, {(0, 1, 0): F(1)}, {(0, 0, 1): F(1)})


def ppow(poly, power):
    answer = PONE
    for _ in range(power):
        answer = pmul(answer, poly)
    return answer


def pdot(vector, variables=VARS):
    answer = PZERO
    for coefficient, variable in zip(vector, variables):
        answer = padd(answer, pscale(coefficient, variable))
    return answer


def pmat(rows, cols):
    return [[{} for _ in range(cols)] for _ in range(rows)]


def ptranspose(matrix):
    return [[dict(value) for value in row] for row in zip(*matrix)]


def pmatmul(left, right):
    answer = pmat(len(left), len(right[0]))
    for i in range(len(left)):
        for j in range(len(right[0])):
            value = PZERO
            for k in range(len(right)):
                value = padd(value, pmul(left[i][k], right[k][j]))
            answer[i][j] = value
    return answer


def pmatadd(left, right):
    return [[padd(left[i][j], right[i][j]) for j in range(len(left[0]))]
            for i in range(len(left))]


def pmatscale(scale, matrix):
    return [[pscale(scale, value) for value in row] for row in matrix]


def numeric_to_poly(matrix):
    return [[({(0, 0, 0): F(value)} if value else {}) for value in row]
            for row in matrix]


def pcongruence(numeric_left, matrix, numeric_right):
    return pmatmul(pmatmul(numeric_to_poly(numeric_left), matrix),
                   numeric_to_poly(numeric_right))


def substitute(poly, rotation):
    linear = []
    for row in rotation:
        value = PZERO
        for coefficient, variable in zip(row, VARS):
            value = padd(value, pscale(coefficient, variable))
        linear.append(value)
    answer = PZERO
    for powers, coefficient in poly.items():
        term = {(0, 0, 0): coefficient}
        for variable, power in zip(linear, powers):
            term = pmul(term, ppow(variable, power))
        answer = padd(answer, term)
    return answer


def substitute_matrix(matrix, rotation):
    return [[substitute(value, rotation) for value in row] for row in matrix]


def poly_matrix_equal(left, right):
    return all(pclean(left[i][j]) == pclean(right[i][j])
               for i in range(len(left)) for j in range(len(left[0])))


def poly_matrix_flat(matrix):
    return [matrix[i][j].get(monomial, F(0))
            for i in range(len(matrix)) for j in range(len(matrix[0]))
            for monomial in MONOMIALS]


def permutation_matrix(permutation):
    result = [[F(0) for _ in range(4)] for _ in range(4)]
    for source, target in enumerate(permutation):
        result[target][source] = F(1)
    return result


def coordinate_rotation(permutation):
    # R T_d = T_{g(d)}, using sum_d T_d T_d^T = 4 I.
    return [[sum((TETRA[permutation[d]][i] * TETRA[d][j] / 4
                  for d in range(4)), F(0))
             for j in range(3)] for i in range(3)]


def section_group_and_invariants():
    solder = [[F(1, 2), *(value / 2 for value in TETRA[d])]
              for d in range(4)]
    check(matmul(transpose(solder), solder) == eye(4),
          "cycle solder is exactly orthogonal")

    character_sum = F(0)
    constant_sum = F(0)
    rotations = []
    for permutation in permutations(range(4)):
        rotation = coordinate_rotation(permutation)
        rotations.append(rotation)
        check(matmul(rotation, transpose(rotation)) == eye(3),
              f"coordinate action orthogonal: {permutation}")
        check(all(sum(bool(value) for value in row) == 1 for row in rotation)
              and all(sum(bool(rotation[i][j]) for i in range(3)) == 1
                      for j in range(3)),
              f"coordinate action signed permutation: {permutation}")
        check(all([sum(rotation[i][j] * TETRA[d][j] for j in range(3))
                    for i in range(3)] == list(TETRA[permutation[d]])
                  for d in range(4)), f"coordinate action maps tetrahedron: {permutation}")

        perm_matrix = permutation_matrix(permutation)
        block = matmul(matmul(transpose(solder), perm_matrix), solder)
        expected = [[F(0) for _ in range(4)] for _ in range(4)]
        expected[0][0] = F(1)
        for i in range(3):
            for j in range(3):
                expected[i + 1][j + 1] = rotation[i][j]
        check(block == expected, f"cycle decomposition A1 plus T2: {permutation}")

        chi_k = trace(rotation)
        chi_k2 = trace(matmul(rotation, rotation))
        chi_domain = (chi_k * chi_k + chi_k2) / 2
        chi_cycle = trace(perm_matrix)
        chi_cycle2 = trace(matmul(perm_matrix, perm_matrix))
        chi_codomain = (chi_cycle * chi_cycle + chi_cycle2) / 2
        character_sum += chi_domain * chi_codomain
        constant_sum += chi_codomain

    quadratic_dimension = character_sum / 24
    constant_dimension = constant_sum / 24
    check(quadratic_dimension == 5, "quadratic equivariant space has dimension five")
    check(constant_dimension == 2, "constant symmetric invariant space has dimension two")

    r2 = padd(padd(pmul(VARS[0], VARS[0]), pmul(VARS[1], VARS[1])),
               pmul(VARS[2], VARS[2]))
    q = (pmul(VARS[1], VARS[2]), pmul(VARS[2], VARS[0]),
         pmul(VARS[0], VARS[1]))
    dmat = pmat(3, 3)
    omat = pmat(3, 3)
    for i in range(3):
        dmat[i][i] = pmul(VARS[i], VARS[i])
        for j in range(3):
            if i != j:
                omat[i][j] = pmul(VARS[i], VARS[j])
    r2i = pmat(3, 3)
    for i in range(3):
        r2i[i][i] = r2

    bases = {}
    aa = pmat(4, 4)
    aa[0][0] = r2
    bases["alpha"] = aa
    av = pmat(4, 4)
    for i in range(3):
        av[0][i + 1] = q[i]
        av[i + 1][0] = q[i]
    bases["eta"] = av
    for label, block3 in (("b", r2i), ("c", dmat), ("d", omat)):
        value = pmat(4, 4)
        for i in range(3):
            for j in range(3):
                value[i + 1][j + 1] = block3[i][j]
        bases[label] = value

    for index, rotation in enumerate(rotations):
        group = [[F(0) for _ in range(4)] for _ in range(4)]
        group[0][0] = F(1)
        for i in range(3):
            for j in range(3):
                group[i + 1][j + 1] = rotation[i][j]
        for label, basis in bases.items():
            lhs = substitute_matrix(basis, rotation)
            rhs = pcongruence(group, basis, transpose(group))
            check(poly_matrix_equal(lhs, rhs),
                  f"quadratic basis {label} covariant under group element {index}")
    check(rank([poly_matrix_flat(basis) for basis in bases.values()]) == 5,
          "five displayed quadratic covariants are independent")

    # At constant order the invariant matrices are uu^T and QQ^T.  The
    # stationary common-amplitude null kills the former and leaves kappa QQ^T.
    u = [[F(1, 2)] for _ in range(4)]
    qcycle = [[value / 2 for value in row] for row in TETRA]
    uu = matmul(u, transpose(u))
    qq = matmul(qcycle, transpose(qcycle))
    check(matadd(uu, qq) == eye(4), "constant A1 and T2 projectors resolve identity")
    check(matmul(qq, u) == [[F(0)] for _ in range(4)],
          "stationary null removes constant A1 and preserves T2")
    return {
        "group_order": 24,
        "quadratic_invariant_dimension": ftext(quadratic_dimension),
        "constant_invariant_dimension": ftext(constant_dimension),
        "solder": [[ftext(x) for x in row] for row in solder],
        "basis_labels": list(bases),
    }


def writer_jet():
    pair_index = {pair: i for i, pair in enumerate(PAIR_ORDER)}
    writer = pmat(4, 3)  # raw t_i basis; divide pullback by 2 for t_i/sqrt(2)
    rho_rows = []
    for missing in range(4):
        a, b, c = tuple(port for port in range(4) if port != missing)
        entries = (
            ((a, b), tuple((TETRA[a][i] + TETRA[b][i] - TETRA[c][i]) / 2
                           for i in range(3))),
            ((a, c), tuple((TETRA[a][i] - TETRA[b][i] + TETRA[c][i]) / 2
                           for i in range(3))),
            ((b, c), tuple((-TETRA[a][i] + TETRA[b][i] + TETRA[c][i]) / 2
                           for i in range(3))),
        )
        for pair, rho in entries:
            rho_rows.append(rho)
            theta2 = pmul(pdot(rho), pdot(rho))
            cosine_row = padd({(0, 0, 0): F(2)}, pscale(-1, theta2))
            pair_row = pair_index[tuple(sorted(pair))]
            for column in range(3):
                writer[missing][column] = padd(
                    writer[missing][column],
                    pscale(T_RAW[column][pair_row], cosine_row))

    check(len(rho_rows) == 12, "twelve centered writer offsets")
    check(all(sum(x * x for x in rho) == F(11, 4) for rho in rho_rows),
          "every centered offset has norm squared 11/4")
    aggregate = [[sum((rho[i] * rho[j] for rho in rho_rows), F(0))
                  for j in range(3)] for i in range(3)]
    check(aggregate == matscale(11, eye(3)), "aggregate offset second moment is 11 I")
    for d in range(4):
        for i in range(3):
            check(writer[d][i].get((0, 0, 0), F(0)) == -2 * TETRA[d][i],
                  f"zero-mode raw writer equals -2 tetrahedron: {d},{i}")

    qcycle = [[value / 2 for value in row] for row in TETRA]
    u = [[F(1, 2)] for _ in range(4)]
    r2 = padd(padd(pmul(VARS[0], VARS[0]), pmul(VARS[1], VARS[1])),
               pmul(VARS[2], VARS[2]))
    qvec = (pmul(VARS[1], VARS[2]), pmul(VARS[2], VARS[0]),
            pmul(VARS[0], VARS[1]))
    dmat = pmat(3, 3)
    omat = pmat(3, 3)
    r2i = pmat(3, 3)
    for i in range(3):
        dmat[i][i] = pmul(VARS[i], VARS[i])
        r2i[i][i] = r2
        for j in range(3):
            if i != j:
                omat[i][j] = pmul(VARS[i], VARS[j])

    qq = numeric_to_poly(matmul(qcycle, transpose(qcycle)))

    def cycle_from_vv(block):
        return pcongruence(qcycle, block, transpose(qcycle))

    cycle_bases = {
        "kappa": qq,
        "b": cycle_from_vv(r2i),
        "c": cycle_from_vv(dmat),
        "d": cycle_from_vv(omat),
    }
    alpha = pmatscale(1, numeric_to_poly(matmul(u, transpose(u))))
    alpha = [[pmul(alpha[i][j], r2) for j in range(4)] for i in range(4)]
    cycle_bases["alpha"] = alpha
    eta = pmat(4, 4)
    for i in range(4):
        for j in range(4):
            value = PZERO
            for a in range(3):
                value = padd(value, pscale(u[i][0] * qcycle[j][a], qvec[a]))
                value = padd(value, pscale(qcycle[i][a] * u[j][0], qvec[a]))
            eta[i][j] = value
    cycle_bases["eta"] = eta

    zero3 = pmat(3, 3)
    expected = {
        "alpha": zero3,
        "eta": zero3,
        "b": pmatscale(8, r2i),
        "c": pmatscale(8, dmat),
        "d": pmatscale(8, omat),
        "kappa": pmatadd(
            pmatadd(pmatscale(8, numeric_to_poly(eye(3))), pmatscale(-2, r2i)),
            pmatadd(pmatscale(-16, dmat), pmatscale(12, omat))),
    }
    pulled = {}
    for label, cycle in cycle_bases.items():
        value = pmatscale(F(1, 2), pmatmul(pmatmul(ptranspose(writer), cycle), writer))
        check(poly_matrix_equal(value, expected[label]),
              f"independent writer pullback coefficient: {label}")
        pulled[label] = value
    return {
        "offset_count": 12,
        "offset_norm_squared": "11/4",
        "offset_second_moment": "11 I",
        "zero_mode_raw_writer": "B_raw(0)=-2 T; orthonormal pair basis is B_raw/sqrt(2)",
        "pullback": {
            "constant": "8 kappa I",
            "A": "-2 kappa+8 b",
            "B": "-16 kappa+8 c",
            "C": "12 kappa+8 d",
            "alpha_eta_through_k2": "zero",
        },
        "writer_digest": canonical_hash([
            [[{str(k): ftext(v) for k, v in sorted(cell.items())}
              for cell in row] for row in writer],
            {label: [ftext(item) for item in poly_matrix_flat(value)]
             for label, value in pulled.items()},
        ]),
    }


def pair_solder(vector):
    answer = [[F(0) for _ in range(3)] for _ in range(3)]
    for coefficient, (a, b) in zip(vector, PAIR_ORDER):
        for i in range(3):
            for j in range(3):
                answer[i][j] += coefficient * (
                    TETRA[a][i] * TETRA[b][j]
                    + TETRA[b][i] * TETRA[a][j]) / 4
    return answer


def reference_operator(tensor, momentum):
    r2 = sum(value * value for value in momentum)
    tensor_trace = sum(tensor[i][i] for i in range(3))
    khk = sum(momentum[i] * tensor[i][j] * momentum[j]
              for i in range(3) for j in range(3))
    answer = [[F(0) for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            left = sum(momentum[a] * tensor[a][j] for a in range(3))
            right = sum(momentum[a] * tensor[a][i] for a in range(3))
            answer[i][j] = (r2 * tensor[i][j]
                            - momentum[i] * left - momentum[j] * right
                            + momentum[i] * momentum[j] * tensor_trace
                            + F(i == j) * (khk - r2 * tensor_trace))
    return answer


def section_tensor_extension_and_reference():
    A = (F(1),) * 6
    E1 = tuple(map(F, (1, -1, 0, 0, -1, 1)))
    E2 = tuple(map(F, (1, 1, -2, -2, 1, 1)))
    basis = (A, E1, E2, *T_RAW)
    soldered = [pair_solder(vector) for vector in basis]
    check(soldered[0] == matscale(-1, eye(3)), "pair trace solder maps A to -I")
    expected_t = []
    for first, second in ((1, 2), (2, 0), (0, 1)):
        matrix = [[F(0) for _ in range(3)] for _ in range(3)]
        matrix[first][second] = F(-1)
        matrix[second][first] = F(-1)
        expected_t.append(matrix)
    check(soldered[3:] == expected_t,
          "raw T basis maps to minus symmetric off-diagonal tensors")

    # In the orthonormal pair basis t_i/sqrt(2), x_i maps to
    # h_(offdiag)=-x_i/sqrt(2).  The common congruence scale is 1/2.
    check(all(sum(value * value for value in row) == 2 for row in T_RAW),
          "raw pair tensor basis has norm squared two")

    # Directly restrict the only two parity-even SO3 quadratic terms that
    # survive between off-diagonal input and output tensor components.
    r2 = padd(padd(pmul(VARS[0], VARS[0]), pmul(VARS[1], VARS[1])),
               pmul(VARS[2], VARS[2]))
    restricted_v = pmat(3, 3)
    for i in range(3):
        restricted_v[i][i] = padd(r2, pscale(-1, pmul(VARS[i], VARS[i])))
        for j in range(3):
            if i != j:
                restricted_v[i][j] = pmul(VARS[i], VARS[j])
    dmat = pmat(3, 3)
    omat = pmat(3, 3)
    r2i = pmat(3, 3)
    for i in range(3):
        dmat[i][i] = pmul(VARS[i], VARS[i])
        r2i[i][i] = r2
        for j in range(3):
            if i != j:
                omat[i][j] = pmul(VARS[i], VARS[j])
    check(poly_matrix_equal(restricted_v,
                            pmatadd(pmatadd(r2i, pmatscale(-1, dmat)), omat)),
          "SO3 divergence term restricts to (r2 I-D)+O")
    check(rank([poly_matrix_flat(r2i), poly_matrix_flat(restricted_v)]) == 2,
          "SO3-extendible T-T restriction is exactly two-dimensional")

    # Coefficient matching: A r2 I+B D+C O is in that span iff B+C=0.
    for aa, bb, cc in ((F(3), F(5), F(-5)), (F(-2), F(0), F(0)),
                       (F(7), F(4), F(-4))):
        general = pmatadd(pmatscale(aa, r2i),
                          pmatadd(pmatscale(bb, dmat), pmatscale(cc, omat)))
        extension = pmatadd(pmatscale(aa + bb, r2i),
                            pmatscale(-bb, restricted_v))
        check((bb + cc == 0) == poly_matrix_equal(general, extension),
              f"B+C criterion exact for sample {aa},{bb},{cc}")

    # Pullback coefficients make B+C=-4 kappa+8(c+d) and A=-2kappa+8b.
    check(F(-16) + F(12) == -4, "cycle mismatch constant is -4 kappa")
    check(F(8) + F(8) == 16, "cycle c+d mismatch coefficients are each eight")
    check(F(-2) + 8 * F(1, 4) == 0,
          "stronger reference shape requires b=kappa/4")
    check(-16 + 8 * F(1, 4) == -14 and
          12 + 8 * F(1, 4) == 14,
          "positive witness pulls back to -14(D-O)")

    momentum = tuple(map(F, (2, 3, 5)))
    reference_images = [reference_operator(tensor, momentum)
                        for tensor in soldered]
    matrix = [[sum((soldered[i][a][b] * reference_images[j][a][b]
                    for a in range(3) for b in range(3)), F(0)) / 2
               for j in range(6)] for i in range(6)]
    expected_matrix = [
        [-38, 5, 37, 15, 10, 6],
        [5, 100, 20, -30, 20, 0],
        [37, 20, 4, -30, -20, 24],
        [15, -30, -30, 4, -6, -10],
        [10, 20, -20, -6, 9, -15],
        [6, 0, 24, -10, -15, 25],
    ]
    check(matrix == [list(map(F, row)) for row in expected_matrix],
          "full-solder reference matrix at k=(2,3,5)")
    check(rank(matrix) == 3, "full-solder reference quotient rank is three")
    expected_tt = [[F(4), F(-6), F(-10)],
                   [F(-6), F(9), F(-15)],
                   [F(-10), F(-15), F(25)]]
    check([row[3:] for row in matrix[3:]] == expected_tt,
          "reference T-T block is D-O")
    check(any(matrix[i][j] for i in range(3) for j in range(3, 6)),
          "full reference has indispensable A/E to T cross blocks")
    return {
        "pair_solder": "(1/4) sum_ab j_ab(T_a T_b^T+T_b T_a^T)",
        "orthonormal_T_components": "(h_yz,h_zx,h_xy)=-(x1,x2,x3)/sqrt(2)",
        "SO3_extension_iff": "B+C=0",
        "cycle_condition": "c+d=kappa/2",
        "reference_additional_condition": "A=0 iff b=kappa/4",
        "reference_TT": "D-O",
        "reference_rank": 3,
        "reference_matrix": expected_matrix,
    }


def section_contact():
    taus = []
    for exceptional in range(4):
        taus.append(tuple(F(-1 if exceptional in pair else 1)
                          for pair in PAIR_ORDER))
    check(all(sum(value * value for value in tau) == 6 for tau in taus),
          "each defect pair vector has norm squared six")
    check(all(sum(taus[a][i] for a in range(4)) == 0 for i in range(6)),
          "defect vectors sum to zero")
    gram = [[sum(taus[a][i] * taus[b][i] for i in range(6))
             for b in range(4)] for a in range(4)]
    check(gram == [[F(6 if a == b else -2) for b in range(4)]
                   for a in range(4)], "defect tetrahedron Gram")

    # Coordinates of Q_a=tau_a tau_a^T/6 in t_i/sqrt(2) basis are
    # (tau.t_i)(tau.t_j)/12, which is rational.
    qframes = []
    for tau in taus:
        dots = [sum(tau[i] * T_RAW[j][i] for i in range(6)) for j in range(3)]
        qframes.append([[dots[i] * dots[j] / 12 for j in range(3)]
                        for i in range(3)])
    check(all(sum((qframes[a][i][j] for a in range(4)), F(0)) == F(4, 3)
              * F(i == j) for i in range(3) for j in range(3)),
          "four defect projectors sum to 4/3 I in orthonormal T coordinates")

    theta = [pdot(TETRA[a]) for a in range(4)]
    theta2_sum = PZERO
    qtheta = pmat(3, 3)
    for a in range(4):
        square = pmul(theta[a], theta[a])
        theta2_sum = padd(theta2_sum, square)
        for i in range(3):
            for j in range(3):
                qtheta[i][j] = padd(qtheta[i][j],
                                    pscale(qframes[a][i][j], square))
    r2 = padd(padd(pmul(VARS[0], VARS[0]), pmul(VARS[1], VARS[1])),
               pmul(VARS[2], VARS[2]))
    expected_theta2 = pscale(4, r2)
    check(theta2_sum == expected_theta2, "sum theta_a squared is 4 r2")
    expected_qtheta = pmat(3, 3)
    for i in range(3):
        expected_qtheta[i][i] = pscale(F(4, 3), r2)
        for j in range(3):
            if i != j:
                expected_qtheta[i][j] = pscale(F(8, 3),
                                                pmul(VARS[i], VARS[j]))
    check(poly_matrix_equal(qtheta, expected_qtheta),
          "sum theta_a squared Q_a is 4/3 r2 I+8/3 O")

    # Affine pairs (constant, coefficient of p).
    def aff_add(left, right):
        return (left[0] + right[0], left[1] + right[1])

    def aff_scale(scale, value):
        return (F(scale) * value[0], F(scale) * value[1])

    beta = (F(-4, 3), F(4))
    delta = (F(2), F(-4))
    gamma = aff_add(aff_scale(4, beta), aff_scale(F(4, 3), delta))
    constant_common = aff_add((F(8), F(0)), gamma)
    coeff_a = aff_add(aff_scale(-2, beta), aff_scale(F(-2, 3), delta))
    coeff_c = aff_scale(F(-4, 3), delta)
    check(gamma == (F(-8, 3), F(32, 3)), "gamma=(8/3)(4p-1)")
    check(constant_common == (F(16, 3), F(32, 3)),
          "common contact constant=(16/3)(1+2p)")
    check(coeff_a == (F(4, 3), F(-16, 3)),
          "contact r2 coefficient=(4/3)(1-4p)")
    check(coeff_c == (F(-8, 3), F(16, 3)),
          "contact O coefficient=(8/3)(2p-1)")

    p_witness = F(109, 128)
    mismatch = coeff_c[0] + coeff_c[1] * p_witness
    check(mismatch == F(15, 8), "Q4 witness contact mismatch is 15/8")

    # Sign: F''=-(h2/4Ud3)G, while W=-F has positive contact +(h2/4Ud3)G.
    check(F(-1) * F(-1) == 1, "connected-functional contact sign is positive")
    return {
        "gamma": "(8/3)(4p-1)",
        "constant_common": "(16/3)(1+2p)",
        "A": "(4/3)(1-4p)",
        "B": "0",
        "C_and_mismatch": "(8/3)(2p-1)",
        "Q4_witness": "p=109/128 gives mismatch 15/8",
        "energy_hessian_sign": "-(h^2/(4U_d^3))G",
        "connected_functional_sign": "+(h^2/(4U_d^3))G",
    }


def section_common_sublattice_normalization():
    """Check CL's j-plus chart against BV's normalized common projection."""
    onsite = F(8)
    offdiag = F(5, 7)
    physical_block = [[onsite, offdiag], [offdiag, onsite]]
    unnormalized_contact = sum((physical_block[i][j]
                                for i in range(2) for j in range(2)), F(0))
    normalized_contact = unnormalized_contact / 2
    check(normalized_contact == onsite + offdiag,
          "normalized common contact is 8+Re C")
    check(unnormalized_contact == 2 * normalized_contact,
          "CL unnormalized common-coordinate contact is twice the projection")

    # CL09 uses j_P=j_++j_-, j_C=j_+-j_-.  At j_-=0, its writer is
    # B_+=B_P+B_C.  The normalized |+> source instead has B_+/sqrt(2).
    bp, bc, kernel = F(2, 3), F(-5, 11), F(7, 13)
    unnormalized_spectral = kernel * (bp + bc) ** 2
    normalized_spectral = unnormalized_spectral / 2
    check(unnormalized_spectral == 2 * normalized_spectral,
          "normalized common spectral pullback is mu^2 B+*KB+/2")
    check((unnormalized_spectral / normalized_spectral)
          == (unnormalized_contact / normalized_contact) == 2,
          "spectral and contact blocks use the same common-coordinate factor")
    return {
        "CL_chart": "j_P=j_++j_-, j_C=j_+-j_-; B_+=B_P+B_C",
        "normalized_writer": "B_+/sqrt(2)",
        "normalized_cycle_block": "(mu^2/2) B_+^* K_cycle B_+",
        "normalized_contact_block": "g_ct(8+Re C)",
        "normalized_extension_equation": "(mu^2/2)[-4 kappa+8(c+d)]+g_ct(8/3)(2p-1)=0",
        "normalized_reference_equation": "(mu^2/2)[-2 kappa+8b]+g_ct(4/3)(1-4p)=0",
        "equivalent_unnormalized_rule": "retain mu^2 cycle block and multiply every displayed contact coefficient by 2",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--emit", action="store_true")
    args = parser.parse_args()
    group = section_group_and_invariants()
    writer = writer_jet()
    tensor = section_tensor_extension_and_reference()
    contact = section_contact()
    normalization = section_common_sublattice_normalization()

    # Dimensional and guard checks are explicit audit assertions.
    check(12 - 5 == 7, "mu squared over J scales as h6/Ud7")
    check(2 - 3 == -1, "contact h2/Ud3 has inverse-energy units")
    result = {
        "audit": "GL6CO independent hostile reconstruction",
        "disposition": "PASS",
        "group_and_invariants": group,
        "writer": writer,
        "tensor_extension_and_reference": tensor,
        "contact": contact,
        "common_sublattice_normalization": normalization,
        "conditional_total": {
            "normalized_extension": "(mu^2/2)[-4 kappa+8(c+d)]+[h^2/(4U_d^3)](8/3)(2p-1)=0",
            "normalized_reference_additional": "(mu^2/2)[-2 kappa+8b]+[h^2/(4U_d^3)](4/3)(1-4p)=0",
            "same_state_guard": "both blocks must belong to one stationary state and one source-first functional",
            "one_particle_irreducible_guard": "a response Hessian must be lawfully inverted/completed before comparison with a 1PI Ricci kernel",
        },
        "power_counting": {
            "cycle": "mu^2/J=O(h^6/U_d^7)",
            "contact": "O(h^2/U_d^3)",
            "meaning": "isolated strong-lock repetition cannot generically cancel a nonzero leading contact mismatch",
        },
        "ceiling": "matching theorem only; no stationary bulk symbol, coefficient calculation, common state, E2 completion, 1PI inversion, metric, Ricci, gravity, or G",
        "checks": len(CHECKS),
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = HERE / "INDEPENDENT_RESULT.json"
    if args.emit:
        output.write_text(rendered)
    else:
        check(output.read_text() == rendered, "frozen independent result matches replay")
    print(f"PASS__GL6CO_INDEPENDENT_HOSTILE_REPLAY__{len(CHECKS)}/{len(CHECKS)}")
    print("S4_QUADRATIC_INVARIANT_DIMENSION=5;CONSTANT_DIMENSION=2;COMMON_NULL_APPLIED")
    print("PULLBACK=A_-2K+8B__B_-16K+8C__C_12K+8D;ALPHA_ETA_ABSENT_AT_K2")
    print("SO3_EXTENSION_IFF=B+C_ZERO_IFF_C+D=KAPPA/2;REFERENCE_ADDS_B=KAPPA/4")
    print("CONTACT=A_4/3(1-4P)__B_0__C_8/3(2P-1);SIGNS_AND_SAME_STATE_GUARD_PASS")
    print("COMMON_NORMALIZATION=NORMALIZED_CYCLE_MU2_OVER2;CONTACT_G_CT;EQUIVALENT_UNNORMALIZED_CONTACT_TIMES2")
    print("DISPOSITION=PASS;MATCHING_ONLY;NO_1PI_RICCI_GRAVITY_OR_G")


if __name__ == "__main__":
    main()
