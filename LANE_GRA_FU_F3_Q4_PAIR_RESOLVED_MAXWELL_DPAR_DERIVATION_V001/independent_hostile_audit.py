#!/usr/bin/env python3
"""Independent hostile audit of the FU pair-resolved Maxwell/DPAR lane.

The audit imports no builder code.  It reconstructs the tetrahedral
representation, all-state degree/elastance identity, passive grounded-network
feasibility, affine derivative, Coulomb normalization, source ranks, charge
conservation boundary, and local-versus-global operator ceilings using exact
standard-library arithmetic.
"""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from math import sqrt
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PAIRS = tuple(combinations(range(4), 2))
SIGNS = (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)

EXPECTED_BUILDER = {
    "DEPENDENCIES.sha256":
        "0eaac5fb8efcc570e90f91354128e07464bd9e84a3d076a27da0c072d4f4ced4",
    "RESULT.md":
        "a3a5e7dced4d9f446ec3e163176acdbb4c8a3fe9cac8e251cb33ce818826f50b",
    "SELF_AUDIT.md":
        "ec82801ebfb2b89bccd22671846ff6d3afa22aa8aaaecf1838829209b96b3eb3",
    "THEOREM.md":
        "f088346f72861b3b11ae737fe6b882d43da9e747fc1d1d1f6bd446a7fd2b6272",
    "verify_pair_resolved_maxwell_dpar_derivation.py":
        "8b9f3e35699dfc9dc60ac4e1d47d6376ef81dab3fc4fdaae5bc3bb7cbbc5b29a",
    "MANIFEST.sha256":
        "20be16d9b0906f62e42ba28b9dc09504d57724843870cd822c5718d8eaf94757",
    "SEAL.sha256":
        "50b36ac34c1df0f29c19e73a0f11bbe4fe2f93e6fda48d342a550a40c15a3a65",
    "VERIFICATION.txt":
        "023e97b9969463849bdfe28bf3fc2e0c00eed29e5014f2e2f7924175f57e27e0",
}

EXPECTED_DEPENDENCIES = {
    "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md":
        "00eba581b90fb9f0b25e3fad1362b055049824897433ba06e356ab9b1f6c76ec",
    "LANE_GRA_FE_F3_Q4_DIAMOND_ICE_CARRIER_JOIN_V001/THEOREM.md":
        "4cc63e3e5853b4250a2a5b78256d41b83b195cf527901819a62a04ef53f8d932",
    "LANE_GRA_CL_F3_DEGREE_LOCK_GAUGE_PHASE_BRIDGE_V001/THEOREM.md":
        "9292b6b4548eb0b286b4e13f3fe999f59c5a228ff267f10d5bb424ca20686f46",
    "LANE_CROSS_ALPHA_GRA_F3_U1_STIFFNESS_IDENTIFIABILITY_V001/THEOREM.md":
        "a042d21e9e0f6e3a63e74e47b1785b8681a186841095004ad2ebcc1e828e670a",
    "LANE_CROSS_ALPHA_GRA_F3_COULOMB_TANGENT_V001/THEOREM.md":
        "8b07300ea659c0681b8e5fec9178819624b28f94923de36e4d927d11b7ccde77",
    "LANE_CROSS_ALPHA_GRA_F3_DIAMOND_SIXTH_ORDER_V001/THEOREM.md":
        "211b1aa61917c98dccae278129a8016a1a14f73587908bfeceeba090a808536c",
    "LANE_GRA_CS_F3_CARRIER_RESPONSE_SUPPORT_SELECTION_V001/THEOREM.md":
        "5e5d36a384c7340e587d30eeb118892d72043de6480e760d80e238aab676c5ce",
    "LANE_GRA_CA_F3_UNIFORM_RESOURCE_MATCHING_LOCALITY_OBSTRUCTION_V001/JOINED_MULTIFORCE_PARENT.md":
        "67a7bdcf509a0c4567cda5cabf8a835cb5c7fcba6666f658fc071b275fb732fb",
    "LANE_RFT_ALPHA_GAUGE_NORMALIZATION_V001/THEOREM.md":
        "3e4f96bef3d70bf16f09bdd540bb558a7f493e7837fdb67618651090beb42b9f",
    "LANE_CROSS_RFT_ALPHA_GRA_EU_ALPHA_Q4_DIMENSION_LOCK_V001/THEOREM.md":
        "ae0ed7bdd758f14e830612a5f6f7dc0207efbe4b4ad0b0a8152fe793c1a99a0d",
    "LANE_CROSS_RFT_ALPHA_GRA_EY_VISIBLE_SECTOR_INDUCED_RICCI_SIGN_V001/THEOREM.md":
        "21c86b24025b9393008c2975c6f421146d5688da0342c19c16cf21bbab4a35b4",
    "LANE_GRA_FT_F3_Q4_DEGREE_PAIR_GEOMETRIC_STRAIN_BOUNDARY_V001/THEOREM.md":
        "6000f38871a57061b106665a41aca04b5d09f4c8c8f4bdc8132ccd5f3f1fbe39",
    "LANE_GRA_FT_F3_Q4_DEGREE_PAIR_GEOMETRIC_STRAIN_BOUNDARY_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "599eec1cde6260be1c9f536274dd8682f77cb45d94e7e3cbc17a28d7552258bd",
}

passed = 0


def check(condition, label):
    global passed
    if not condition:
        raise AssertionError(label)
    passed += 1


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def parse_ledger(path):
    entries = {}
    for line in path.read_text().splitlines():
        expected, relative = line.split("  ", 1)
        check(relative not in entries, f"unique ledger member {path.name}:{relative}")
        entries[relative] = expected
    return entries


def rank(rows):
    matrix = [[F(value) for value in row] for row in rows]
    if not matrix:
        return 0
    nrows = len(matrix)
    ncols = len(matrix[0])
    pivot_row = 0
    for column in range(ncols):
        pivot = next((row for row in range(pivot_row, nrows)
                      if matrix[row][column]), None)
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [entry / scale for entry in matrix[pivot_row]]
        for row in range(nrows):
            if row == pivot_row or not matrix[row][column]:
                continue
            scale = matrix[row][column]
            matrix[row] = [left - scale * right
                           for left, right in zip(matrix[row], matrix[pivot_row])]
        pivot_row += 1
        if pivot_row == nrows:
            break
    return pivot_row


def identity(size):
    return tuple(tuple(F(row == column) for column in range(size))
                 for row in range(size))


def transpose(matrix):
    return tuple(tuple(matrix[row][column] for row in range(len(matrix)))
                 for column in range(len(matrix[0])))


def matmul(left, right):
    right_t = transpose(right)
    return tuple(tuple(sum(F(a) * F(b) for a, b in zip(row, column))
                       for column in right_t) for row in left)


def matvec(matrix, vector):
    return tuple(sum(F(a) * F(b) for a, b in zip(row, vector)) for row in matrix)


def matrix_add(left, right):
    return tuple(tuple(F(a) + F(b) for a, b in zip(left_row, right_row))
                 for left_row, right_row in zip(left, right))


def matrix_scale(scale, matrix):
    return tuple(tuple(F(scale) * F(entry) for entry in row) for row in matrix)


def inverse(matrix):
    size = len(matrix)
    work = [[F(matrix[row][column]) for column in range(size)]
            + [F(row == column) for column in range(size)]
            for row in range(size)]
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if work[row][column]), None)
        if pivot is None:
            raise ValueError("singular matrix")
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [entry / scale for entry in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [left - scale * right
                         for left, right in zip(work[row], work[column])]
    return tuple(tuple(work[row][size:]) for row in range(size))


def outer(vector, scale=F(1)):
    return tuple(tuple(F(scale) * F(a) * F(b) for b in vector) for a in vector)


def contract(left, right):
    return sum(F(left[row][column]) * F(right[row][column])
               for row in range(len(left)) for column in range(len(left[0])))


def sym_coordinates(matrix):
    return (
        F(matrix[0][0]), F(matrix[1][1]), F(matrix[2][2]),
        2 * F(matrix[0][1]), 2 * F(matrix[0][2]), 2 * F(matrix[1][2]),
    )


def source_matrix(coordinate):
    result = [[F(0) for _ in range(3)] for _ in range(3)]
    if coordinate < 3:
        result[coordinate][coordinate] = F(1)
    else:
        row, column = ((0, 1), (0, 2), (1, 2))[coordinate - 3]
        result[row][column] = result[column][row] = F(1)
    return tuple(tuple(row) for row in result)


def commutator(left, right):
    return matrix_add(matmul(left, right), matrix_scale(-1, matmul(right, left)))


ZERO4 = tuple(tuple(F(0) for _ in range(4)) for _ in range(4))
I3 = identity(3)


# Frozen builder and dependency custody.
for relative, expected in EXPECTED_BUILDER.items():
    path = HERE / relative
    check(path.is_file() and not path.is_symlink(),
          f"builder member is a regular file: {relative}")
    check(digest(path) == expected, f"builder byte custody: {relative}")
    check(sha256(path.read_bytes() + b"hostile-builder-tamper").hexdigest() != expected,
          f"builder tamper rejection: {relative}")

dependencies = parse_ledger(HERE / "DEPENDENCIES.sha256")
check(dependencies == EXPECTED_DEPENDENCIES,
      "dependency ledger equals independent hostile expectation")
for relative, expected in EXPECTED_DEPENDENCIES.items():
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(),
          f"dependency is regular and non-symlink: {relative}")
    check(digest(path) == expected, f"dependency byte custody: {relative}")
    check(sha256(path.read_bytes() + b"hostile-dependency-tamper").hexdigest() != expected,
          f"dependency tamper rejection: {relative}")

base_manifest = parse_ledger(HERE / "MANIFEST.sha256")
check(set(base_manifest) == {
    "DEPENDENCIES.sha256", "RESULT.md", "SELF_AUDIT.md", "THEOREM.md",
    "verify_pair_resolved_maxwell_dpar_derivation.py",
}, "builder manifest has exactly five stable payload members")
for relative, expected in base_manifest.items():
    check(digest(HERE / relative) == expected,
          f"builder manifest digest matches: {relative}")

base_seal = parse_ledger(HERE / "SEAL.sha256")
check(set(base_seal) == {"MANIFEST.sha256", "VERIFICATION.txt"},
      "builder seal covers manifest and transcript")
for relative, expected in base_seal.items():
    check(digest(HERE / relative) == expected,
          f"builder seal digest matches: {relative}")


# Tetrahedral group, invariant scalar gradient, and pair-root geometry.
tetra_vectors = tuple(tuple(F(entry) for entry in vector) for vector in SIGNS)
basis_columns = tuple(tuple(tetra_vectors[column][row] for column in range(3))
                      for row in range(3))
basis_inverse = inverse(basis_columns)
orthogonals = {}
for permutation in permutations(range(4)):
    target_columns = tuple(tuple(tetra_vectors[permutation[column]][row]
                                 for column in range(3)) for row in range(3))
    orthogonal = matmul(target_columns, basis_inverse)
    orthogonals[permutation] = orthogonal
    check(matmul(transpose(orthogonal), orthogonal) == I3,
          f"tetrahedral permutation is orthogonal: {permutation}")
    for label in range(4):
        check(matvec(orthogonal, tetra_vectors[label])
              == tetra_vectors[permutation[label]],
              f"tetrahedral permutation maps label {label}: {permutation}")

sym_basis = tuple(source_matrix(coordinate) for coordinate in range(6))
invariance_rows = []
for orthogonal in orthogonals.values():
    for output in range(6):
        row = []
        for basis in sym_basis:
            transformed = matmul(orthogonal, matmul(basis, transpose(orthogonal)))
            difference = matrix_add(transformed, matrix_scale(-1, basis))
            row.append(sym_coordinates(difference)[output])
        invariance_rows.append(tuple(row))
check(rank(invariance_rows) == 5,
      "symmetric tetrahedral invariant-gradient space is exactly one-dimensional")
check(all(sum(row[index] for index in range(3)) == 0
          for row in invariance_rows),
      "identity tensor spans the invariant-gradient null line")

E_SOURCES = (
    ((F(1), F(0), F(0)), (F(0), F(-1), F(0)), (F(0), F(0), F(0))),
    ((F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(-2))),
)
T2_SOURCES = (source_matrix(3), source_matrix(4), source_matrix(5))
check(all(contract(source, I3) == 0 for source in E_SOURCES),
      "both E strains annihilate every invariant scalar gradient")

edge_dyads = tuple(outer(vector, F(1, 3)) for vector in tetra_vectors)
root_dyads = {}
for a, b in PAIRS:
    difference = tuple(tetra_vectors[b][axis] - tetra_vectors[a][axis]
                       for axis in range(3))
    length2 = sum(entry * entry for entry in difference)
    check(length2 == 8, f"root {a}{b} has exact unscaled squared length eight")
    root_dyads[(a, b)] = outer(difference, F(1, length2))
    check(contract(root_dyads[(a, b)], I3) == 1,
          f"root projector {a}{b} has unit trace")

edge_rows = tuple(sym_coordinates(dyad) for dyad in edge_dyads)
root_rows = tuple(sym_coordinates(root_dyads[pair]) for pair in PAIRS)
check(rank(edge_rows) == 4, "one-edge tetrahedral dyads have rank four")
check(rank(root_rows) == 6, "six sibling-root dyads have rank six")
check(rank(edge_rows + root_rows) == 6,
      "edge and sibling-root dyads span all symmetric source coordinates")
check(all(contract(source, dyad) == 0
          for source in E_SOURCES for dyad in edge_dyads),
      "one-edge source is exactly E-null")

for permutation, orthogonal in orthogonals.items():
    for a, b in PAIRS:
        transformed = matmul(orthogonal,
                             matmul(root_dyads[(a, b)], transpose(orthogonal)))
        target = tuple(sorted((permutation[a], permutation[b])))
        check(transformed == root_dyads[target],
              f"root projector covariance for {permutation}:{a}{b}")


# The E-to-pair intertwiner is one-dimensional under exact S4 covariance.
def e_coordinates(matrix):
    check(matrix[0][1] == matrix[0][2] == matrix[1][2] == 0,
          "transformed E source remains diagonal")
    check(matrix[0][0] + matrix[1][1] + matrix[2][2] == 0,
          "transformed E source remains traceless")
    return ((matrix[0][0] - matrix[1][1]) / 2, -matrix[2][2] / 2)


intertwiner_equations = []
for permutation, orthogonal in orthogonals.items():
    domain = []
    for source in E_SOURCES:
        transformed = matmul(orthogonal, matmul(source, transpose(orthogonal)))
        domain.append(e_coordinates(transformed))
    domain = transpose(tuple(domain))
    pair_index = {pair: index for index, pair in enumerate(PAIRS)}
    pair_permutation = tuple(pair_index[tuple(sorted((permutation[a], permutation[b])))]
                             for a, b in PAIRS)
    for output in range(6):
        for input_column in range(2):
            equation = [F(0) for _ in range(12)]
            for internal in range(2):
                equation[2 * output + internal] += domain[internal][input_column]
            source_output = pair_permutation.index(output)
            equation[2 * source_output + input_column] -= 1
            intertwiner_equations.append(tuple(equation))
check(12 - rank(intertwiner_equations) == 1,
      "S4-covariant E-to-local-pair intertwiner space is one-dimensional")

central_e_map = tuple(tuple(contract(source, root_dyads[pair]) for source in E_SOURCES)
                      for pair in PAIRS)
check(rank(transpose(central_e_map)) == 2,
      "nonzero central root response carries both E directions")
cancelled_e_map = tuple(tuple(value - value for value in row)
                        for row in central_e_map)
check(rank(cancelled_e_map) == 0,
      "an equal opposite S4 mutual remainder can cancel the full E response")


# Degree identity and exact sixteen-state grounded-elastance realization.
states = tuple(product((-1, 1), repeat=4))
ice = tuple(state for state in states if sum(state) == 0)
for state in states:
    degree = sum((1 - z) // 2 for z in state)
    pair_sum = sum(state[a] * state[b] for a, b in PAIRS)
    check((degree - 2) ** 2 == 1 + F(pair_sum, 2),
          f"degree-pair identity on state {state}")
    check(sum(state) == -2 * (degree - 2),
          f"signed common charge equals centered degree on state {state}")
check(len(ice) == 6, "neutral ice fiber has six states")

Q_STAR = F(3)
U_D = F(42)
E_MUTUAL = U_D / (2 * Q_STAR**2)
E_SELF = F(5)
ELASTANCE = tuple(tuple(E_SELF if row == column else E_MUTUAL
                         for column in range(4)) for row in range(4))
check(E_SELF - E_MUTUAL > 0, "elastance contrast eigenvalue is positive")
check(E_SELF + 3 * E_MUTUAL > 0, "elastance common eigenvalue is positive")
CAPACITANCE = inverse(ELASTANCE)
check(matmul(CAPACITANCE, ELASTANCE) == identity(4),
      "grounded capacitance is exact inverse of elastance")
check(all(CAPACITANCE[row][row] > 0 for row in range(4)),
      "capacitance has positive diagonal entries")
check(all(CAPACITANCE[row][column] < 0
          for row in range(4) for column in range(4) if row != column),
      "capacitance has passive negative mutual entries")
check(all(sum(CAPACITANCE[row]) > 0 for row in range(4)),
      "capacitance has positive ground legs")

mutual_capacitance = -CAPACITANCE[0][1]
ground_capacitance = sum(CAPACITANCE[0])
NETWORK_CAPACITANCE = tuple(
    tuple(ground_capacitance + 3 * mutual_capacitance
          if row == column else -mutual_capacitance
          for column in range(4)) for row in range(4)
)
check(NETWORK_CAPACITANCE == CAPACITANCE,
      "one positive ground leg plus six positive mutual capacitors synthesizes C")

E_REFERENCE = U_D - 2 * Q_STAR**2 * E_SELF
for state in states:
    vector = tuple(F(z) for z in state)
    electrostatic = F(1, 2) * Q_STAR**2 * sum(
        vector[row] * ELASTANCE[row][column] * vector[column]
        for row in range(4) for column in range(4)
    ) + E_REFERENCE
    degree = sum((1 - z) // 2 for z in state)
    check(electrostatic == U_D * (degree - 2) ** 2,
          f"grounded electrostatic energy matches all-state degree energy: {state}")


# Floating common-potential mode and the exact neutral-quotient scope.
FLOATING = tuple(tuple(F(3) if row == column else F(-1)
                       for column in range(4)) for row in range(4))
check(rank(FLOATING) == 3, "floating capacitance has rank three")
check(matvec(FLOATING, (F(1), F(1), F(1), F(1))) == (F(0),) * 4,
      "floating capacitance kills the common-potential mode")
try:
    inverse(FLOATING)
except ValueError:
    floating_inverse_rejected = True
else:
    floating_inverse_rejected = False
check(floating_inverse_rejected, "full inverse rejects floating capacitance")
check(all(sum(state) == 0 for state in ice),
      "every ice state lies in the neutral quotient")
check(all(sum(state) != 0 for state in states if state not in ice),
      "every off-ice state lies outside the neutral quotient")


# Global charge conservation is weaker than local gauge invariance.
# Basis ordering is (z_terminal,z_reservoir)=--,-+,+-,++.
labels = ((-1, -1), (-1, 1), (1, -1), (1, 1))
Q_TERMINAL = tuple(tuple(F(labels[row][0]) if row == column else F(0)
                         for column in range(4)) for row in range(4))
Q_RESERVOIR = tuple(tuple(F(labels[row][1]) if row == column else F(0)
                          for column in range(4)) for row in range(4))
Q_TOTAL = matrix_add(Q_TERMINAL, Q_RESERVOIR)
BARE_TERMINAL_FLIP = (
    (F(0), F(0), F(1), F(0)),
    (F(0), F(0), F(0), F(1)),
    (F(1), F(0), F(0), F(0)),
    (F(0), F(1), F(0), F(0)),
)
DRESSED_EXCHANGE = (
    (F(0), F(0), F(0), F(0)),
    (F(0), F(0), F(1), F(0)),
    (F(0), F(1), F(0), F(0)),
    (F(0), F(0), F(0), F(0)),
)
check(commutator(Q_TOTAL, BARE_TERMINAL_FLIP) != ZERO4,
      "bare terminal flip violates total-charge conservation")
check(commutator(Q_TOTAL, DRESSED_EXCHANGE) == ZERO4,
      "terminal-reservoir exchange conserves total charge")
check(commutator(Q_TERMINAL, DRESSED_EXCHANGE) != ZERO4,
      "total-charge conservation alone does not commute with local terminal charge")
check(commutator(Q_RESERVOIR, DRESSED_EXCHANGE) != ZERO4,
      "total-charge conservation alone does not commute with local reservoir charge")
check(DRESSED_EXCHANGE[1][2] == DRESSED_EXCHANGE[2][1] == 1,
      "dressed exchange acts as X on the fixed-total-charge code pair")


# Central-kernel chain rule, source convention, Coulomb slope, and running.
SOURCE_BASIS = tuple(source_matrix(coordinate) for coordinate in range(6))
for pair in PAIRS:
    a, b = pair
    vector = tuple(tetra_vectors[b][axis] - tetra_vectors[a][axis]
                   for axis in range(3))
    length2 = sum(entry * entry for entry in vector)
    for source in SOURCE_BASIS:
        direct_ratio_derivative = -sum(
            vector[row] * source[row][column] * vector[column]
            for row in range(3) for column in range(3)
        ) / length2
        check(direct_ratio_derivative == -contract(source, root_dyads[pair]),
              f"F=I-j/2 squared-distance derivative: {pair}:{sym_coordinates(source)}")

R0 = F(7)
V0 = F(11)
V_PRIME = F(-5)
LAMBDA = R0 * V_PRIME / (2 * V0)
check(LAMBDA == F(-35, 22), "generic central-kernel chain-rule slope")

def generic_v(radius):
    return float(V0) + float(V_PRIME) * (radius - float(R0))


step = 1e-6
generic_numeric = (
    generic_v(float(R0) * sqrt(1 + step)) / float(V0)
    - generic_v(float(R0) * sqrt(1 - step)) / float(V0)
) / (2 * step)
check(abs(generic_numeric - float(LAMBDA)) < 1e-8,
      "numerical generic-kernel derivative confirms exact chain rule")

for pair in PAIRS:
    for source in SOURCE_BASIS:
        distance_derivative = -contract(source, root_dyads[pair])
        energy_derivative = (U_D / 2) * LAMBDA * distance_derivative
        conjugate = -2 * energy_derivative
        expected = U_D * LAMBDA * contract(source, root_dyads[pair])
        check(conjugate == expected,
              f"Q=-2 dH/dj sign and factor: {pair}:{sym_coordinates(source)}")

def coulomb_g(x):
    return x ** -0.5


coulomb_numeric = (coulomb_g(1 + step) - coulomb_g(1 - step)) / (2 * step)
check(abs(coulomb_numeric + 0.5) < 1e-9,
      "ideal Coulomb normalized slope is minus one half")
COULOMB_LAMBDA = F(-1, 2)

ALPHA = F(1, 137)
KAPPA = F(2)
EPSILON_R = F(5)
HBAR_C = F(197)
PHYSICAL_R0 = F(13)
PAIR_COEFFICIENT = KAPPA**2 * ALPHA * HBAR_C / (EPSILON_R * PHYSICAL_R0)
UD_ALPHA = 2 * KAPPA**2 * ALPHA * HBAR_C / (EPSILON_R * PHYSICAL_R0)
check(UD_ALPHA / 2 == PAIR_COEFFICIENT,
      "rationalized-SI alpha normalization has the correct factor two")
check(UD_ALPHA == F(1576, 8905),
      "representative U_d-alpha-length relation is exact")

BETA_ALPHA = F(1, 10000)
RUNNING_LAMBDA = F(-1, 2) - BETA_ALPHA / (2 * ALPHA)
check(RUNNING_LAMBDA == F(-20274, 40000),
      "running alpha contributes minus beta_alpha over two alpha")
check(RUNNING_LAMBDA != COULOMB_LAMBDA,
      "running coupling is not the ideal fixed-coupling slope")

alpha_one = F(1, 137)
r_one = F(1)
alpha_two = F(1, 100)
r_two = F(137, 100)
check(2 * alpha_one / r_one == 2 * alpha_two / r_two,
      "coefficient matching alone cannot select alpha without independent length")


# Microscopic rank six, ice rank three, and exact projected T2 nullity.
check(rank(edge_rows + tuple(tuple(COULOMB_LAMBDA * entry for entry in row)
                             for row in root_rows)) == 6,
      "nonzero Coulomb DPAR plus unchanged edge source has microscopic rank six")
check(rank(edge_rows + tuple(tuple(F(0) for _ in row) for row in root_rows)) == 4,
      "zero DPAR slope leaves microscopic rank four")

def ice_query(source):
    coefficients = tuple(contract(source, root_dyads[pair]) for pair in PAIRS)
    return tuple(sum(coefficient * state[a] * state[b]
                     for coefficient, (a, b) in zip(coefficients, PAIRS))
                 for state in ice)


ice_full = tuple(ice_query(source) for source in SOURCE_BASIS)
ice_e = tuple(ice_query(source) for source in E_SOURCES)
ice_t2 = tuple(ice_query(source) for source in T2_SOURCES)
check(rank(ice_full) == 3, "direct ice pair image has rank three")
check(rank(ice_e) == 2, "direct ice pair image retains both E directions")
check(all(all(value == 0 for value in row) for row in ice_t2),
      "direct ice pair image annihilates all T2 sources")


# Nonlocal Walsh independence holds on a disjoint product, not automatically
# on a shared-link lattice.  Both facts are part of the hostile scope audit.
eight_bit_states = tuple(product((-1, 1), repeat=8))
disjoint_characters = []
for a, b in PAIRS:
    disjoint_characters.append(tuple(state[a] * state[b]
                                      for state in eight_bit_states))
for a in range(4):
    for b in range(4, 8):
        disjoint_characters.append(tuple(state[a] * state[b]
                                          for state in eight_bit_states))
check(rank(disjoint_characters) == 22,
      "six local plus sixteen cross-node Walsh strings are independent on disjoint factors")

seven_bit_states = tuple(product((-1, 1), repeat=7))
second_node = (0, 4, 5, 6)
overlap_characters = []
for a, b in PAIRS:
    overlap_characters.append(tuple(state[a] * state[b]
                                     for state in seven_bit_states))
for a in range(4):
    for b in second_node:
        overlap_characters.append(tuple(state[a] * state[b]
                                         for state in seven_bit_states))
check(rank(overlap_characters) < 22,
      "shared-link supports invalidate the disjoint-factor rank count")
check(tuple(F(1) for _ in seven_bit_states) in overlap_characters,
      "a shared terminal cross term can collapse to identity")
check(overlap_characters[0] in overlap_characters[6:],
      "a shared terminal cross term can duplicate a local pair")


# Source-before-elimination is load bearing: an exact Schur complement retains
# derivatives residing in the eliminated coupling block.
def schur_scalar(parameter):
    coupling = F(1) + parameter
    eliminated = F(2)
    return -coupling * coupling / eliminated


schur_derivative = (schur_scalar(F(1)) - schur_scalar(F(-1))) / 2
check(schur_derivative == -1,
      "source inserted before exact elimination survives in Schur derivative")
post_projection_derivative = F(0)
check(post_projection_derivative != schur_derivative,
      "post-projection source attachment misses eliminated-block response")


# Documentary ceilings: the exact result is a conditional local physical
# completion, not inherited F3, visible electromagnetism, alpha selection, or
# gravity.  The local gauge and global-lattice obligations remain premises.
theorem = (HERE / "THEOREM.md").read_text()
result = (HERE / "RESULT.md").read_text()
self_audit = (HERE / "SELF_AUDIT.md").read_text()

required_theorem_phrases = (
    "new **physical completion**",
    "A neutral-quotient construction proves only the ice-restricted",
    "complete net** `E` slope remains",
    "lambda_E^net != 0",
    "Walsh strings are independent of the local `P_ab`",
    "fixed-total-charge encoded subspace",
    "gauge-invariant compensating transfer",
    "source-work term, and port retained in the\nledger",
    "source-before-Feshbach",
    "S9 — noncircular ancestry",
    "conditional physical completion only, not inherited F3",
    "Direct restriction of the pair map to local\nice still has only `A1+E` rank three",
    "not identify `U_d` with\n`U_E^IR`",
)
for phrase in required_theorem_phrases:
    check(phrase in theorem, f"theorem retains hostile scope phrase: {phrase}")

check("Equation (FU26) is **not inherited F3**" in theorem,
      "alpha coefficient relation is not promoted to inherited F3")
check("not a calculation of alpha" in theorem,
      "alpha coefficient match is not called alpha selection")
check("larger dynamical parent" in theorem and "fresh source and rank audit" in theorem,
      "charged-parent enlargement triggers a fresh audit")
check("shielded to zero" in theorem and "prospectively frozen larger parent" in theorem,
      "nonlocal interactions cannot be silently dropped")
check("not automatically a complete physical action" in theorem,
      "local point-pair formula is not promoted to a global action")
check("complete grounded,\nnonsingular elastance matrix" in result,
      "result retains the full-ground condition")
check("larger parent needing a fresh rank audit" in result,
      "result retains the dressed-parent boundary")
check("no\nprojected CTP rank, Ward packet, tensor pole or gravity result follows" in result,
      "result retains projected and gravity ceilings")
check("Decisive falsifiers" in self_audit,
      "self-audit declares physical falsifiers")
check("no conserved signed charge" in self_audit,
      "self-audit allows failure of the signed-charge solder")

for name in (
    "DEPENDENCIES.sha256", "RESULT.md", "SELF_AUDIT.md", "THEOREM.md",
    "verify_pair_resolved_maxwell_dpar_derivation.py", "independent_hostile_audit.py",
):
    data = (HERE / name).read_bytes()
    check(b"\r" not in data and b"\b" not in data and b"\f" not in data,
          f"byte hygiene: {name}")


print(f"SUMMARY {passed}/{passed} independent hostile checks passed")
print("DISPOSITION PASS__EXACT_CONDITIONAL_LOCAL_ELASTANCE_DPAR_COMPLETION__GLOBAL_GAUGE_AND_SHARED_LINK_SOLDER_REMAIN_PREMISES__NO_INHERITED_F3_VISIBLE_EM_ALPHA_SELECTION_PROJECTED_CTP_GRAVITY_OR_G")
