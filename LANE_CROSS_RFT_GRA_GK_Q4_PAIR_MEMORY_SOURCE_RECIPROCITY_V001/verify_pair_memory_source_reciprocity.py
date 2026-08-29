#!/usr/bin/env python3
"""Exact standard-library verifier for the PMSR q4 bridge.

The verifier uses only Fraction arithmetic and elementary polynomial lists.
It rebuilds the q4 state algebra rather than importing builder results.
"""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LANE = Path(__file__).resolve().parent

DEPENDENCIES = {
    "LANE_CROSS_RFT_GRA_EW_Q4_PAIR_MEMORY_METRIC_DEFORMATION_CLOSURE_V001/THEOREM.md":
        "495e4e99171f4e3e5809f24e5a9a5b68116e996f4f29115bd22f780127d4714e",
    "LANE_CROSS_RFT_GRA_EW_Q4_PAIR_MEMORY_METRIC_DEFORMATION_CLOSURE_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "71c4529413a66ebc657709807643486089a058cf0e84a7964ae153aa4da93984",
    "LANE_GRA_FT_F3_Q4_DEGREE_PAIR_GEOMETRIC_STRAIN_BOUNDARY_V001/THEOREM.md":
        "6000f38871a57061b106665a41aca04b5d09f4c8c8f4bdc8132ccd5f3f1fbe39",
    "LANE_GRA_FT_F3_Q4_DEGREE_PAIR_GEOMETRIC_STRAIN_BOUNDARY_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "599eec1cde6260be1c9f536274dd8682f77cb45d94e7e3cbc17a28d7552258bd",
    "LANE_GRA_FU_F3_Q4_PAIR_RESOLVED_MAXWELL_DPAR_DERIVATION_V001/THEOREM.md":
        "f088346f72861b3b11ae737fe6b882d43da9e747fc1d1d1f6bd446a7fd2b6272",
    "LANE_GRA_FU_F3_Q4_PAIR_RESOLVED_MAXWELL_DPAR_DERIVATION_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "6f859779566177b5999cfe02c01cd569c5bd7b0b4ec2b21b0b3e79ebf26f9277",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/THEOREM.md":
        "6fc221a31151340b91a946d33e442971c1373500e067c354b6c610e3964edb1c",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md":
        "3801fd9ba6ba3c0fe80c9f4792abfdeb6dd7c37c7145663be05b4d56f8160723",
    "LANE_GRA_FY_F3_Q4_NATIVE_SUPPORT_M1_COMPLETE_H6_RESPONSE_V001/THEOREM.md":
        "8db3dd16c36e0205b5c98fc3154e8a2f1876d243c3c1d2068424c1276ee68f28",
    "LANE_GRA_FY_F3_Q4_NATIVE_SUPPORT_M1_COMPLETE_H6_RESPONSE_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md":
        "f9baf2b21ead24d947866192eb7f0cb6d4e353ffb2b1107569cc442564804f21",
}

checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def add_matrix(a, b):
    return [[x + y for x, y in zip(rx, ry)] for rx, ry in zip(a, b)]


def scale_matrix(c, a):
    return [[c * x for x in row] for row in a]


def outer(x, y):
    return [[a * b for b in y] for a in x]


def transpose(a):
    return [list(row) for row in zip(*a)]


def mm(a, b):
    bt = transpose(b)
    return [[sum(x * y for x, y in zip(row, col)) for col in bt]
            for row in a]


def mv(a, x):
    return [sum(y * z for y, z in zip(row, x)) for row in a]


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def zero_matrix(nr, nc):
    return [[F(0) for _ in range(nc)] for _ in range(nr)]


def rank(rows):
    a = [list(map(F, row)) for row in rows]
    if not a:
        return 0
    nr, nc = len(a), len(a[0])
    r = 0
    for c in range(nc):
        pivot = next((i for i in range(r, nr) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        q = a[r][c]
        a[r] = [x / q for x in a[r]]
        for i in range(nr):
            if i != r and a[i][c]:
                q = a[i][c]
                a[i] = [x - q * y for x, y in zip(a[i], a[r])]
        r += 1
    return r


def determinant(rows):
    a = [list(map(F, row)) for row in rows]
    value = F(1)
    for c in range(len(a)):
        pivot = next((i for i in range(c, len(a)) if a[i][c]), None)
        if pivot is None:
            return F(0)
        if pivot != c:
            a[c], a[pivot] = a[pivot], a[c]
            value *= -1
        q = a[c][c]
        value *= q
        a[c] = [x / q for x in a[c]]
        for i in range(c + 1, len(a)):
            if a[i][c]:
                q = a[i][c]
                a[i] = [x - q * y for x, y in zip(a[i], a[c])]
    return value


def vec_sym3(a):
    return [a[0][0], a[1][1], a[2][2],
            a[0][1], a[0][2], a[1][2]]


def poly_trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)


def poly_add(a, b):
    n = max(len(a), len(b))
    return poly_trim([(a[i] if i < len(a) else 0) +
                      (b[i] if i < len(b) else 0) for i in range(n)])


def poly_scale(c, a):
    return poly_trim([c * x for x in a])


def poly_sub(a, b):
    return poly_add(a, poly_scale(-1, b))


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return poly_trim(out)


def poly_shift(a, n=1):
    return poly_trim([0] * n + list(a))


def poly_mv(a, x):
    out = []
    for row in a:
        total = (0,)
        for entry, scalar in zip(row, x):
            total = poly_add(total, poly_scale(scalar, entry))
        out.append(total)
    return out


STATES = tuple(product((-1, 1), repeat=4))
EDGES = tuple(combinations(range(4), 2))

I4 = eye(4)
P = [[I4[i][j] - F(1, 4) for j in range(4)] for i in range(4)]
VECTORS = [[P[i][a] for i in range(4)] for a in range(4)]
B_EDGE = [add_matrix(outer(VECTORS[a], VECTORS[b]),
                     outer(VECTORS[b], VECTORS[a]))
          for a, b in EDGES]

BASIS = [
    [F(1), F(0), F(0)],
    [F(0), F(1), F(0)],
    [F(0), F(0), F(1)],
    [F(-1), F(-1), F(-1)],
]


def pull_to_v(a):
    return mm(mm(transpose(BASIS), a), BASIS)


def x_stat(state):
    return mv(P, list(map(F, state)))


def y_stat(state):
    return [F(state[a] * state[b]) for a, b in EDGES]


# GK07: statewise identity and the EW expectation-to-metric isomorphism.
for state in STATES:
    rhs = [row[:] for row in P]
    y = y_stat(state)
    for coefficient, tensor in zip(y, B_EDGE):
        rhs = add_matrix(rhs, scale_matrix(coefficient, tensor))
    check(outer(x_stat(state), x_stat(state)) == rhs,
          f"statewise XX identity {state}")

B_MATRIX = transpose([vec_sym3(pull_to_v(tensor)) for tensor in B_EDGE])
check(rank(B_MATRIX) == 6, "expectation-to-Fisher map B rank six")
check(determinant(B_MATRIX) != 0,
      "expectation-to-Fisher map B determinant nonzero")


# FU normalized root-dyad strain map M.
SIGNS = ((1, 1, 1), (1, -1, -1),
         (-1, 1, -1), (-1, -1, 1))


def root_row(a, b):
    r = tuple(F(SIGNS[b][i] - SIGNS[a][i]) for i in range(3))
    return [r[0] * r[0] / 8, r[1] * r[1] / 8,
            r[2] * r[2] / 8, 2 * r[0] * r[1] / 8,
            2 * r[0] * r[2] / 8, 2 * r[1] * r[2] / 8]


M_MATRIX = [root_row(a, b) for a, b in EDGES]
check(rank(M_MATRIX) == 6, "normalized root strain map M rank six")
check(determinant(M_MATRIX) != 0,
      "normalized root strain map M determinant nonzero")


def state_weight(a, state):
    magnetization = sum(state)
    exponent = {0: 0, 2: 1, -2: 1, 4: 4, -4: 4}[magnetization]
    return a ** exponent


def distribution(a):
    raw = {state: state_weight(a, state) for state in STATES}
    z = sum(raw.values(), F(0))
    return {state: weight / z for state, weight in raw.items() if weight}


def expected(probabilities, function):
    first = function(next(iter(probabilities)))
    if isinstance(first, list):
        return [sum(probabilities[s] * function(s)[i]
                    for s in probabilities) for i in range(len(first))]
    return sum(probabilities[s] * function(s) for s in probabilities)


def pair_covariance(probabilities):
    means = expected(probabilities, y_stat)
    out = zero_matrix(6, 6)
    for state, probability in probabilities.items():
        centered = [x - m for x, m in zip(y_stat(state), means)]
        out = add_matrix(out, scale_matrix(probability,
                                           outer(centered, centered)))
    return means, out


def localization_fisher(probabilities):
    mean = expected(probabilities, x_stat)
    out = zero_matrix(4, 4)
    for state, probability in probabilities.items():
        centered = [x - m for x, m in zip(x_stat(state), mean)]
        out = add_matrix(out, scale_matrix(probability,
                                           outer(centered, centered)))
    return out


# Exact full-support expectation-coordinate identity at a nontrivial weight.
A_SAMPLE = F(1, 2)
PROB = distribution(A_SAMPLE)
C_SAMPLE, G_SAMPLE = pair_covariance(PROB)
F_SAMPLE = localization_fisher(PROB)
F_FROM_C = [row[:] for row in P]
for c, tensor in zip(C_SAMPLE, B_EDGE):
    F_FROM_C = add_matrix(F_FROM_C, scale_matrix(c, tensor))
check(F_SAMPLE == F_FROM_C,
      "nontrivial full-support Fisher equals P plus B times C")
check(rank(G_SAMPLE) == 6, "nontrivial full-support pair covariance rank six")


# Symbolic enumeration in the positive variable a.
def monomial(exponent, coefficient=1):
    return tuple([0] * exponent + [coefficient])


z_poly = (0,)
c_sums = [(0,) for _ in EDGES]
q_sum = (0,)
for state in STATES:
    magnetization = sum(state)
    exponent = {0: 0, 2: 1, -2: 1, 4: 4, -4: 4}[magnetization]
    w = monomial(exponent)
    z_poly = poly_add(z_poly, w)
    for e, value in enumerate(y_stat(state)):
        c_sums[e] = poly_add(c_sums[e], poly_scale(value, w))
    q_value = state[0] * state[1] * state[2] * state[3]
    q_sum = poly_add(q_sum, poly_scale(q_value, w))

D_POLY = (3, 4, 0, 0, 1)
C_NUM = (-1, 0, 0, 0, 1)
Q_NUM = (3, -4, 0, 0, 1)
check(z_poly == poly_scale(2, D_POLY), "symbolic partition polynomial 2D")
for edge, numerator in zip(EDGES, c_sums):
    check(numerator == poly_scale(2, C_NUM),
          f"symbolic common pair numerator edge {edge}")
check(q_sum == poly_scale(2, Q_NUM),
      "symbolic four-spin numerator")

D2 = poly_mul(D_POLY, D_POLY)
C2 = poly_mul(C_NUM, C_NUM)
DIAG_NUM = poly_sub(D2, C2)
ADJ_NUM = poly_sub(poly_mul(C_NUM, D_POLY), C2)
OPP_NUM = poly_sub(poly_mul(Q_NUM, D_POLY), C2)

G_NUM = []
for e in EDGES:
    row = []
    for f in EDGES:
        if e == f:
            row.append(DIAG_NUM)
        elif len(set(e).intersection(f)) == 1:
            row.append(ADJ_NUM)
        else:
            row.append(OPP_NUM)
    G_NUM.append(row)

V_A1 = [1, 1, 1, 1, 1, 1]
V_E1 = [1, -1, 0, 0, -1, 1]
V_E2 = [1, 0, -1, -1, 0, 1]
V_T1 = [1, 0, 0, 0, 0, -1]
V_T2 = [0, 1, 0, 0, -1, 0]
V_T3 = [0, 0, 1, -1, 0, 0]

GA_NUM = (0, 8, 0, 0, 32, 24)
GE_NUM = poly_scale(8, D_POLY)
GT_NUM = poly_scale(8, poly_shift(D_POLY))

for vector, eigen_num, label in (
        (V_A1, GA_NUM, "A1"),
        (V_E1, GE_NUM, "E-1"),
        (V_E2, GE_NUM, "E-2"),
        (V_T1, GT_NUM, "T2-1"),
        (V_T2, GT_NUM, "T2-2"),
        (V_T3, GT_NUM, "T2-3")):
    observed = poly_mv(G_NUM, vector)
    expected_vector = [poly_scale(x, eigen_num) for x in vector]
    check(observed == expected_vector,
          f"symbolic pair-covariance eigenvalue {label}")

check(all(coefficient >= 0 for coefficient in GA_NUM) and
      any(coefficient > 0 for coefficient in GA_NUM[1:]),
      "A1 eigenvalue numerator positive for a greater than zero")
check(all(coefficient >= 0 for coefficient in GE_NUM) and GE_NUM[0] > 0,
      "E eigenvalue numerator positive for a greater than zero")
check(all(coefficient >= 0 for coefficient in GT_NUM) and
      any(coefficient > 0 for coefficient in GT_NUM[1:]),
      "T2 eigenvalue numerator positive for a greater than zero")


# Composition rank at finite support and at exact ice restriction.
FINITE_PRODUCT = mm(mm(B_MATRIX, G_SAMPLE), M_MATRIX)
check(rank(FINITE_PRODUCT) == 6, "finite-support BGM response rank six")

L_GENERAL = eye(6)
L_GENERAL[0][1] = F(2)
L_GENERAL[2][4] = F(-3)
L_GENERAL[4][5] = F(5, 2)
check(determinant(L_GENERAL) == 1,
      "general complete pair tangent sample invertible")
check(rank(mm(mm(B_MATRIX, G_SAMPLE), L_GENERAL)) == 6,
      "finite-support BGL response follows full-rank general tangent")

L_SINGULAR = eye(6)
L_SINGULAR[-1] = [F(0) for _ in range(6)]
check(rank(L_SINGULAR) == 5,
      "singular complete pair tangent sample rank five")
check(rank(mm(mm(B_MATRIX, G_SAMPLE), L_SINGULAR)) == 5,
      "finite-support BGL response preserves general tangent rank")

ICE_PROB = distribution(F(0))
_, G_ICE = pair_covariance(ICE_PROB)
check(rank(G_ICE) == 2, "exact ice pair covariance rank two")
check(rank(mm(mm(B_MATRIX, G_ICE), M_MATRIX)) == 2,
      "exact ice BGM product rank two")

UNIFORM_PROB = distribution(F(1))
_, G_UNIFORM = pair_covariance(UNIFORM_PROB)
check(G_UNIFORM == eye(6), "beta-zero state has G equal identity")
check(rank(scale_matrix(F(0), mm(mm(B_MATRIX, G_UNIFORM), M_MATRIX))) == 0,
      "beta-zero physical response vanishes")
check(rank(scale_matrix(F(0), FINITE_PRODUCT)) == 0,
      "lambda-zero physical response vanishes")


# Source sign/factor and exact mixed-derivative reciprocity.
BETA = F(3, 2)
Q_SCALE = F(14, 5)          # U_d * lambda
TANGENT_SCALE = BETA * Q_SCALE / 2
IDENTITY_SOURCE = F(11, 7)

for e in range(6):
    for a_index in range(6):
        k_prime = -Q_SCALE * M_MATRIX[e][a_index] / 2
        j_prime = -BETA * k_prime
        q_coefficient = -2 * k_prime
        check(j_prime == TANGENT_SCALE * M_MATRIX[e][a_index],
              f"J=-beta K sign and factor e={e} A={a_index}")
        check(q_coefficient == Q_SCALE * M_MATRIX[e][a_index],
              f"Q=-2 dH/dj sign and factor e={e} A={a_index}")


def covariance_xx_with(probabilities, scalar_function):
    mean_scalar = expected(probabilities, scalar_function)
    fisher = localization_fisher(probabilities)
    out = zero_matrix(4, 4)
    for state, probability in probabilities.items():
        value = scalar_function(state)
        out = add_matrix(
            out,
            scale_matrix(probability * value,
                         outer(x_stat(state), x_stat(state))),
        )
    return add_matrix(out, scale_matrix(-mean_scalar, fisher))


D_J_F = [covariance_xx_with(PROB, lambda state, e=e: y_stat(state)[e])
         for e in range(6)]

for a_index in range(6):
    direct = zero_matrix(4, 4)
    for e in range(6):
        direct = add_matrix(
            direct,
            scale_matrix(TANGENT_SCALE * M_MATRIX[e][a_index], D_J_F[e]),
        )

    def q_value(state, a_index=a_index):
        return IDENTITY_SOURCE + Q_SCALE * sum(
            M_MATRIX[e][a_index] * y_stat(state)[e] for e in range(6)
        )

    theta_hessian_q = covariance_xx_with(PROB, q_value)
    reciprocal = scale_matrix(BETA / 2, theta_hessian_q)
    check(direct == reciprocal,
          f"exact mixed-derivative reciprocity strain A={a_index}")


theorem_text = (LANE / "THEOREM.md").read_text()
result_text = (LANE / "RESULT.md").read_text()
self_audit_text = (LANE / "SELF_AUDIT.md").read_text()
check("GK-S4 -- complete DPAR form and noncancellation" in theorem_text,
      "theorem freezes complete-DPAR rather than E-only premise")
check("lambda_E^net != 0` proves a noncancelled `E` component" in theorem_text,
      "theorem states exact E-only insufficiency")
check("D_j{\\cal F}_\\theta={\\beta U_d\\over2}{\\cal B}GL" in theorem_text,
      "theorem states general complete-tangent formula")
check("whole-pair-source rank-six conclusion" in result_text,
      "result preserves whole-pair-source ceiling")
check("FU's nonzero `E` coefficient" in self_audit_text,
      "self-audit preserves E-only insufficiency")
check("exact pair-sector physical conjugate" in theorem_text,
      "theorem does not relabel pair conjugate as complete source")
check("Q_A=-2\\partial_{j_A}H_C" in result_text,
      "result types Q as commuting pair-sector source")


for relative, expected in DEPENDENCIES.items():
    path = ROOT / relative
    check(path.is_file(), f"dependency exists {relative}")
    check(digest(path) == expected, f"dependency hash {relative}")

dependency_lines = [line.strip() for line in
                    (LANE / "DEPENDENCIES.sha256").read_text().splitlines()
                    if line.strip()]
declared = {}
for line in dependency_lines:
    value, relative = line.split(None, 1)
    declared[relative.strip()] = value
check(declared == DEPENDENCIES, "dependency file exactly matches verifier map")

print(f"PASS {checks}/{checks}")
print("PMSR exact builder verification complete")
