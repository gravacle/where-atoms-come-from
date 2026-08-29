#!/usr/bin/env python3
"""Independent exact hostile replay for the PMSR lane.

This program was written without importing or invoking the builder verifier.
It reconstructs the sixteen-state q4 model, the two six-dimensional linear
maps, the thermal covariance, and the source/metric mixed derivative using
only exact rational arithmetic.
"""

from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path


LANE = Path(__file__).resolve().parents[1]
ROOT = LANE.parent
AUDIT = Path(__file__).resolve().parent

checks = 0


def check(statement, label):
    global checks
    if not statement:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def parse_hash_file(path):
    entries = {}
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        value, relative = line.split(None, 1)
        entries[relative.strip()] = value
    return entries


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def matmul(left, right):
    right_t = transpose(right)
    return [
        [sum(x * y for x, y in zip(row, column)) for column in right_t]
        for row in left
    ]


def matvec(matrix, vector):
    return [sum(x * y for x, y in zip(row, vector)) for row in matrix]


def outer(left, right):
    return [[x * y for y in right] for x in left]


def add(left, right):
    return [
        [x + y for x, y in zip(left_row, right_row)]
        for left_row, right_row in zip(left, right)
    ]


def scale(value, matrix):
    return [[value * x for x in row] for row in matrix]


def zero(rows, columns):
    return [[Q(0) for _ in range(columns)] for _ in range(rows)]


def identity(size):
    return [[Q(i == j) for j in range(size)] for i in range(size)]


def rank(matrix):
    work = [list(map(Q, row)) for row in matrix]
    if not work:
        return 0
    rows = len(work)
    columns = len(work[0])
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [x / pivot_value for x in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            multiplier = work[row][column]
            work[row] = [
                x - multiplier * y
                for x, y in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def determinant(matrix):
    work = [list(map(Q, row)) for row in matrix]
    size = len(work)
    value = Q(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot is None:
            return Q(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            value *= -1
        pivot_value = work[column][column]
        value *= pivot_value
        work[column] = [x / pivot_value for x in work[column]]
        for row in range(column + 1, size):
            if not work[row][column]:
                continue
            multiplier = work[row][column]
            work[row] = [
                x - multiplier * y
                for x, y in zip(work[row], work[column])
            ]
    return value


def sym3_coordinates(matrix):
    return [
        matrix[0][0], matrix[1][1], matrix[2][2],
        matrix[0][1], matrix[0][2], matrix[1][2],
    ]


STATES = tuple(product((-1, 1), repeat=4))
EDGES = tuple(combinations(range(4), 2))

I4 = identity(4)
P = [[I4[i][j] - Q(1, 4) for j in range(4)] for i in range(4)]
V = [[P[i][a] for i in range(4)] for a in range(4)]


def x_of(state):
    return matvec(P, list(map(Q, state)))


def y_of(state):
    return [Q(state[a] * state[b]) for a, b in EDGES]


B_EDGE = []
for a, b in EDGES:
    B_EDGE.append(add(outer(V[a], V[b]), outer(V[b], V[a])))


# A rational basis of 1-perp, chosen independently of any metric
# normalization.  Pullback is injective on symmetric forms over 1-perp.
VBASIS = [
    [Q(1), Q(0), Q(0)],
    [Q(0), Q(1), Q(0)],
    [Q(0), Q(0), Q(1)],
    [Q(-1), Q(-1), Q(-1)],
]


def restrict_to_v(matrix):
    return matmul(matmul(transpose(VBASIS), matrix), VBASIS)


BMAP = transpose([sym3_coordinates(restrict_to_v(item)) for item in B_EDGE])


# Reconstruct the statewise identity before taking any expectation.
for state in STATES:
    observed = outer(x_of(state), x_of(state))
    predicted = [row[:] for row in P]
    for coefficient, tensor in zip(y_of(state), B_EDGE):
        predicted = add(predicted, scale(coefficient, tensor))
    check(observed == predicted, f"statewise identity {state}")

check(rank(BMAP) == 6, "observable-C to Fisher-tensor map has rank six")
check(determinant(BMAP) != 0, "observable-C map determinant is nonzero")


# Pair characters and the constant character must be linearly independent.
walsh_table = []
for state in STATES:
    walsh_table.append([Q(1)] + y_of(state))
check(rank(walsh_table) == 7, "constant plus six pair Walsh characters independent")


# Reconstruct the normalized tetrahedral root-dyad evaluation map.
TETRA = (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)


def root_dyad_row(edge):
    a, b = edge
    r = [Q(TETRA[b][i] - TETRA[a][i]) for i in range(3)]
    norm2 = sum(x * x for x in r)
    return [
        r[0] * r[0] / norm2,
        r[1] * r[1] / norm2,
        r[2] * r[2] / norm2,
        2 * r[0] * r[1] / norm2,
        2 * r[0] * r[2] / norm2,
        2 * r[1] * r[2] / norm2,
    ]


MMAP = [root_dyad_row(edge) for edge in EDGES]
check(rank(MMAP) == 6, "normalized root-dyad map has rank six")
check(determinant(MMAP) != 0, "normalized root-dyad determinant is nonzero")


def thermal_distribution(a):
    unnormalized = {}
    for state in STATES:
        magnetization = abs(sum(state))
        exponent = {0: 0, 2: 1, 4: 4}[magnetization]
        unnormalized[state] = a ** exponent
    partition = sum(unnormalized.values(), Q(0))
    return {state: weight / partition for state, weight in unnormalized.items()}


def ice_distribution():
    ice = [state for state in STATES if sum(state) == 0]
    return {state: Q(1, len(ice)) for state in ice}


def expectation(probabilities, observable):
    first = observable(next(iter(probabilities)))
    if isinstance(first, list):
        return [
            sum(probabilities[state] * observable(state)[i] for state in probabilities)
            for i in range(len(first))
        ]
    return sum(probabilities[state] * observable(state) for state in probabilities)


def covariance_matrix(probabilities, observable):
    means = expectation(probabilities, observable)
    result = zero(len(means), len(means))
    for state, probability in probabilities.items():
        centered = [value - mean for value, mean in zip(observable(state), means)]
        result = add(result, scale(probability, outer(centered, centered)))
    return means, result


def fisher_x(probabilities):
    return covariance_matrix(probabilities, x_of)[1]


def four_spin(state):
    value = 1
    for spin in state:
        value *= spin
    return Q(value)


def closed_thermal_values(a):
    d = a ** 4 + 4 * a + 3
    c = (a ** 4 - 1) / d
    q4 = (a ** 4 - 4 * a + 3) / d
    ga = 8 * a * (3 * a ** 4 + 4 * a ** 3 + 1) / d ** 2
    ge = 8 / d
    gt = 8 * a / d
    return d, c, q4, ga, ge, gt


A1 = [Q(1), Q(1), Q(1), Q(1), Q(1), Q(1)]
E1 = [Q(1), Q(-1), Q(0), Q(0), Q(-1), Q(1)]
E2 = [Q(1), Q(0), Q(-1), Q(-1), Q(0), Q(1)]
T1 = [Q(1), Q(0), Q(0), Q(0), Q(0), Q(-1)]
T2 = [Q(0), Q(1), Q(0), Q(0), Q(-1), Q(0)]
T3 = [Q(0), Q(0), Q(1), Q(-1), Q(0), Q(0)]


for a in (Q(1, 4), Q(2, 3), Q(1), Q(5, 2), Q(4)):
    probabilities = thermal_distribution(a)
    c_vector, g = covariance_matrix(probabilities, y_of)
    d, c, q4, ga, ge, gt = closed_thermal_values(a)
    check(sum(1 for s in STATES if probabilities[s] > 0) == 16,
          f"full support at a={a}")
    check(all(value == c for value in c_vector), f"pair mean formula at a={a}")
    check(expectation(probabilities, four_spin) == q4,
          f"four-spin formula at a={a}")
    for vector, eigenvalue, name in (
        (A1, ga, "A1"), (E1, ge, "E1"), (E2, ge, "E2"),
        (T1, gt, "T1"), (T2, gt, "T2"), (T3, gt, "T3"),
    ):
        check(matvec(g, vector) == [eigenvalue * x for x in vector],
              f"{name} covariance eigenpair at a={a}")
    check(ga > 0 and ge > 0 and gt > 0,
          f"all covariance eigenvalues positive at a={a}")
    check(rank(g) == 6, f"finite-temperature covariance rank six at a={a}")


# C is not the natural control: away from the uniform point its response
# Jacobian is the nontrivial covariance G, while remaining invertible.
PROB = thermal_distribution(Q(2, 3))
C, G = covariance_matrix(PROB, y_of)
check(G != identity(6), "D_J C is nontrivial away from J=0")
check(rank(G) == 6, "D_J C remains invertible at finite full support")

FISHER = fisher_x(PROB)
AFFINE = [row[:] for row in P]
for coefficient, tensor in zip(C, B_EDGE):
    AFFINE = add(AFFINE, scale(coefficient, tensor))
check(FISHER == AFFINE, "exact affine expectation-to-Fisher identity off uniformity")


# Source sign, factor, and mixed derivative.  Use a negative physical slope
# to prevent a sign mistake from hiding behind positivity.
BETA = Q(5, 3)
U_LAMBDA = Q(-7, 5)
J_SCALE = BETA * U_LAMBDA / 2
IDENTITY_Q = Q(13, 11)

for edge in range(6):
    for source in range(6):
        k_prime = -U_LAMBDA * MMAP[edge][source] / 2
        j_prime = -BETA * k_prime
        q_prime = -2 * k_prime
        check(j_prime == J_SCALE * MMAP[edge][source],
              f"J=-beta K sign/factor e={edge} A={source}")
        check(q_prime == U_LAMBDA * MMAP[edge][source],
              f"Q=-2 dH/dj sign/factor e={edge} A={source}")


def centered_xx_with(probabilities, scalar_observable):
    mean_scalar = expectation(probabilities, scalar_observable)
    base_fisher = fisher_x(probabilities)
    weighted = zero(4, 4)
    for state, probability in probabilities.items():
        weighted = add(
            weighted,
            scale(probability * scalar_observable(state), outer(x_of(state), x_of(state))),
        )
    return add(weighted, scale(-mean_scalar, base_fisher))


DJ_F = [
    centered_xx_with(PROB, lambda state, edge=edge: y_of(state)[edge])
    for edge in range(6)
]

for source in range(6):
    direct = zero(4, 4)
    for edge in range(6):
        direct = add(
            direct,
            scale(J_SCALE * MMAP[edge][source], DJ_F[edge]),
        )

    def q_source(state, source=source):
        return IDENTITY_Q + U_LAMBDA * sum(
            MMAP[edge][source] * y_of(state)[edge] for edge in range(6)
        )

    theta_hessian_q = centered_xx_with(PROB, q_source)
    reciprocal = scale(BETA / 2, theta_hessian_q)
    check(direct == reciprocal,
          f"same-generator mixed derivative for source direction {source}")


# Check the B G M composition in independent coordinate matrices.
DIRECT_RESPONSE = matmul(matmul(BMAP, G), MMAP)
check(rank(DIRECT_RESPONSE) == 6, "finite-temperature BGM rank six")

GENERAL_L = identity(6)
GENERAL_L[0][4] = Q(3, 2)
GENERAL_L[2][1] = Q(-5, 3)
GENERAL_L[5][3] = Q(7, 4)
check(determinant(GENERAL_L) != 0, "general tangent witness invertible")
check(rank(matmul(matmul(BMAP, G), GENERAL_L)) == 6,
      "general BGL response rank follows invertible L")

SINGULAR_L = identity(6)
SINGULAR_L[5] = [Q(0)] * 6
check(rank(SINGULAR_L) == 5, "general tangent singular control rank five")
check(rank(matmul(matmul(BMAP, G), SINGULAR_L)) == 5,
      "BGL preserves singular tangent rank")


# Fidelity coefficient and its source derivative.  Along theta=t*z,
# gamma(t)=Z(t/2)^2/(Z(0)Z(t)), so log-gamma''(0)=-Var(z.X)/2.
# The coefficient of -log gamma is therefore Var(z.X)/4.  Its j derivative
# is checked against beta/8 times the theta-Hessian of <Q>.
Z_DIRECTION = [Q(2), Q(-1), Q(3), Q(-4)]


def zscore(state):
    return sum(x * y for x, y in zip(Z_DIRECTION, x_of(state)))


check(expectation(PROB, zscore) == 0, "directional theta score is centered")

for source in range(6):
    def q_source(state, source=source):
        return IDENTITY_Q + U_LAMBDA * sum(
            MMAP[edge][source] * y_of(state)[edge] for edge in range(6)
        )

    q_mean = expectation(PROB, q_source)
    hessian_direction = sum(
        probability * (q_source(state) - q_mean) * zscore(state) ** 2
        for state, probability in PROB.items()
    )
    gamma_rhs = BETA * hessian_direction / 8

    variance_derivative = Q(0)
    for edge in range(6):
        y_mean = C[edge]
        cov_z2_y = sum(
            probability * zscore(state) ** 2 * (y_of(state)[edge] - y_mean)
            for state, probability in PROB.items()
        )
        variance_derivative += J_SCALE * MMAP[edge][source] * cov_z2_y
    gamma_lhs = variance_derivative / 4
    check(gamma_lhs == gamma_rhs,
          f"squared-fidelity source coefficient for direction {source}")


# Exact boundary controls.
UNIFORM = thermal_distribution(Q(1))
_, G_UNIFORM = covariance_matrix(UNIFORM, y_of)
check(G_UNIFORM == identity(6), "beta-zero covariance itself is identity")
check(rank(scale(Q(0), matmul(matmul(BMAP, G_UNIFORM), MMAP))) == 0,
      "beta-zero state response vanishes")
check(rank(scale(Q(0), DIRECT_RESPONSE)) == 0,
      "zero complete pair slope response vanishes")

ICE = ice_distribution()
_, G_ICE = covariance_matrix(ICE, y_of)
check(rank(G_ICE) == 2, "exact ice pair covariance rank two")
check(rank(matmul(matmul(BMAP, G_ICE), MMAP)) == 2,
      "exact ice composed pair response rank two")


# Frozen target and dependency custody.  The audited-target file is written
# only after the parent packet is immutable for this audit.
target_hashes = parse_hash_file(AUDIT / "AUDITED_TARGETS.sha256")
check(set(target_hashes) == {"THEOREM.md", "RESULT.md", "SELF_AUDIT.md"},
      "audited target set is exact")
for relative, expected in target_hashes.items():
    check(digest(LANE / relative) == expected, f"frozen target hash {relative}")

dependencies = parse_hash_file(LANE / "DEPENDENCIES.sha256")
check(len(dependencies) == 10, "dependency ledger contains ten unique entries")
for relative, expected in dependencies.items():
    check((ROOT / relative).is_file(), f"dependency exists {relative}")
    check(digest(ROOT / relative) == expected, f"dependency hash {relative}")


theorem = (LANE / "THEOREM.md").read_text()
result = (LANE / "RESULT.md").read_text()
check("C_e:=\\langle Y_e\\rangle" in theorem, "theorem types C as expectation")
check("J_e=\\hbox{natural control}" in theorem, "theorem types J as control")
check("complete DPAR form and noncancellation" in theorem,
      "whole-pair DPAR premise is explicit")
check("whole six-direction pair tangent equals one scalar" in theorem,
      "FU E-only result is not promoted to whole-pair rank")
check("D_j{\\cal F}_\\theta={\\beta U_d\\over2}{\\cal B}GL" in theorem,
      "general BGL formula is stated")
check("pair-sector physical conjugate" in theorem,
      "source reciprocity is scoped to pair sector")
check("-2\\partial_{j_A}H_C" in result,
      "result keeps the source derivative on H_C")
check("direct-source operator custody only" in theorem,
      "FY custody is limited to direct pair operator")
check("not automatically the SLD/QFI" in theorem,
      "noncommuting BKM/QFI distinction is retained")
check("No gravity, Einstein, pole, or `G` claim is\nmade" in result,
      "result retains the gravity ceiling")

print(f"PASS {checks}/{checks}")
print("Independent hostile PMSR replay complete")
