#!/usr/bin/env python3
"""Exact and numerical verifier for the FU pair-resolved Maxwell/DPAR lane.

The exact checks use only Python's standard library and Fraction arithmetic.
Floating-point checks are restricted to differentiating the displayed smooth
central kernels and are backed by exact symbolic formulas.
"""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from math import sqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LANE = Path(__file__).resolve().parent

DEPENDENCIES = {
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

MANIFEST_FILES = {
    "DEPENDENCIES.sha256",
    "RESULT.md",
    "SELF_AUDIT.md",
    "THEOREM.md",
    "verify_pair_resolved_maxwell_dpar_derivation.py",
}

AUDIT_FILES = {
    "AUDIT_MANIFEST.sha256",
    "AUDIT_SEAL.sha256",
    "INDEPENDENT_HOSTILE_AUDIT.md",
    "INDEPENDENT_HOSTILE_VERIFICATION.txt",
    "independent_hostile_audit.py",
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


def rank(rows):
    matrix = [[F(value) for value in row] for row in rows]
    if not matrix:
        return 0
    nr, nc = len(matrix), len(matrix[0])
    pivot_row = 0
    for column in range(nc):
        pivot = next((row for row in range(pivot_row, nr)
                      if matrix[row][column]), None)
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        for row in range(nr):
            if row != pivot_row and matrix[row][column]:
                scale = matrix[row][column]
                matrix[row] = [left - scale * right
                               for left, right in
                               zip(matrix[row], matrix[pivot_row])]
        pivot_row += 1
        if pivot_row == nr:
            break
    return pivot_row


def identity(size):
    return tuple(tuple(F(row == column) for column in range(size))
                 for row in range(size))


def transpose(matrix):
    return tuple(tuple(matrix[row][column] for row in range(len(matrix)))
                 for column in range(len(matrix[0])))


def matmul(left, right):
    return tuple(tuple(sum(left[row][k] * right[k][column]
                           for k in range(len(right)))
                       for column in range(len(right[0])))
                 for row in range(len(left)))


def matscale(scale, matrix):
    return tuple(tuple(F(scale) * value for value in row) for row in matrix)


def matadd(left, right):
    return tuple(tuple(a + b for a, b in zip(lrow, rrow))
                 for lrow, rrow in zip(left, right))


def inverse(matrix):
    size = len(matrix)
    work = [list(map(F, matrix[row])) + list(identity(size)[row])
            for row in range(size)]
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if work[row][column]), None)
        if pivot is None:
            raise ValueError("singular matrix")
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [value / scale for value in work[column]]
        for row in range(size):
            if row != column and work[row][column]:
                scale = work[row][column]
                work[row] = [left - scale * right
                             for left, right in zip(work[row], work[column])]
    return tuple(tuple(work[row][size:]) for row in range(size))


def matvec(matrix, vector):
    return tuple(sum(value * component for value, component in zip(row, vector))
                 for row in matrix)


def dot(left, right):
    return sum(F(a) * F(b) for a, b in zip(left, right))


def dyad_coordinates(vector):
    x, y, z = map(F, vector)
    return (x*x, y*y, z*z, 2*x*y, 2*x*z, 2*y*z)


def symmetric_matrix(coordinates):
    xx, yy, zz, xy, xz, yz = map(F, coordinates)
    return ((xx, xy, xz), (xy, yy, yz), (xz, yz, zz))


def symmetric_coordinates(matrix):
    return (matrix[0][0], matrix[1][1], matrix[2][2],
            matrix[0][1], matrix[0][2], matrix[1][2])


# Dependency bytes and declared semantic custody.
for relative, expected in DEPENDENCIES.items():
    path = ROOT / relative
    check(path.is_file() and digest(path) == expected,
          f"dependency custody {relative}")
    check(sha256(path.read_bytes() + b"tamper").hexdigest() != expected,
          f"dependency tamper rejection {relative}")

texts = {relative: (ROOT / relative).read_text()
         for relative in DEPENDENCIES}
bs = texts["LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md"]
fe = texts["LANE_GRA_FE_F3_Q4_DIAMOND_ICE_CARRIER_JOIN_V001/THEOREM.md"]
cl = texts["LANE_GRA_CL_F3_DEGREE_LOCK_GAUGE_PHASE_BRIDGE_V001/THEOREM.md"]
u1 = texts["LANE_CROSS_ALPHA_GRA_F3_U1_STIFFNESS_IDENTIFIABILITY_V001/THEOREM.md"]
ct = texts["LANE_CROSS_ALPHA_GRA_F3_COULOMB_TANGENT_V001/THEOREM.md"]
d6 = texts["LANE_CROSS_ALPHA_GRA_F3_DIAMOND_SIXTH_ORDER_V001/THEOREM.md"]
cs = texts["LANE_GRA_CS_F3_CARRIER_RESPONSE_SUPPORT_SELECTION_V001/THEOREM.md"]
ca = texts["LANE_GRA_CA_F3_UNIFORM_RESOURCE_MATCHING_LOCALITY_OBSTRUCTION_V001/JOINED_MULTIFORCE_PARENT.md"]
agn = texts["LANE_RFT_ALPHA_GAUGE_NORMALIZATION_V001/THEOREM.md"]
eu = texts["LANE_CROSS_RFT_ALPHA_GRA_EU_ALPHA_Q4_DIMENSION_LOCK_V001/THEOREM.md"]
ey = texts["LANE_CROSS_RFT_ALPHA_GRA_EY_VISIBLE_SECTOR_INDUCED_RICCI_SIGN_V001/THEOREM.md"]
ft = texts["LANE_GRA_FT_F3_Q4_DEGREE_PAIR_GEOMETRIC_STRAIN_BOUNDARY_V001/THEOREM.md"]
ft_audit = texts["LANE_GRA_FT_F3_Q4_DEGREE_PAIR_GEOMETRIC_STRAIN_BOUNDARY_V001/INDEPENDENT_HOSTILE_AUDIT.md"]

check("No distance or dimension occurs" in bs,
      "BS keeps the degree parent metric free")
check("declared completion slot" in bs and "symbolic slot alone does not prove" in bs,
      "BS port slot is not a frozen physical completion")
check("not automatically a material lattice" in fe
      and "bind `a_*` to an absolute physical length" in fe,
      "FE withholds physical coexistence and absolute length")
check("a discrete divergence constraint" in cl and "staggered background charge" in cl,
      "CL identifies U_d as a discrete Gauss-charge penalty")
check("exact **Gauss-charge penalty**" in u1
      and "not thereby a Maxwell electric stiffness" in u1,
      "U1SI separates U_d from transverse Maxwell stiffness")
check("does not identify the emergent" in ct
      and "U(1) with visible electromagnetism" in ct,
      "Coulomb-tangent lane withholds visible-EM identity")
check("not evidence that F3 selects any spatial support" in d6
      and "visible electromagnetism" in d6,
      "sixth-order lane retains supplied-support and visible-EM ceilings")
check("(J^\psi_{uv})^2" in cs and "occupied--blank boundary cost" in cs,
      "existing current square is edge-local boundary cost")
check("post-emergence matching" in ca and "not the pregeometric action" in ca,
      "visible-U1 action remains a post-emergence match")
check("alpha(\mu)={e(\mu)^2\over4\pi}" in agn,
      "alpha normalization dependency owns the rationalized invariant")
check("Circular-Maxwell control" in eu and "Different-front control" in eu,
      "EU forbids circular and different-front Maxwell imports")
check("physical record-to-metric" in ey
      and "Alpha itself is not an extra" in ey,
      "EY retains solder and no-double-counting ceilings")
check("g'(1)=" in ft and "lambda" in ft
      and "exact rank-six statement in this lane is microscopic" in ft,
      "FT owns only conditional microscopic DPAR rank six")
check("neither inherited nor adopted" in ft_audit and "rank three" in ft_audit,
      "FT hostile audit retains unadopted and projected-rank ceilings")


# Tetrahedral geometry and exact S4 action.
SIGNS = ((1, 1, 1), (1, -1, -1),
         (-1, 1, -1), (-1, -1, 1))
PAIRS = tuple(combinations(range(4), 2))
EDGE = tuple(dyad_coordinates(vector) for vector in SIGNS)
ROOT = {}
for a, b in PAIRS:
    vector = tuple(SIGNS[b][axis] - SIGNS[a][axis] for axis in range(3))
    length2 = dot(vector, vector)
    ROOT[(a, b)] = tuple(value / length2 for value in dyad_coordinates(vector))
    check(length2 == 8, f"unnormalized root {a+1}{b+1} has length squared eight")

check(rank(EDGE) == 4, "tetrahedral one-edge dyads have rank four")
check(rank(tuple(ROOT[pair] for pair in PAIRS)) == 6,
      "normalized sibling-root dyads have rank six")

# Orthogonal matrices realizing all 24 label permutations.
source_columns = tuple(tuple(F(SIGNS[column][row]) for column in range(3))
                       for row in range(3))
source_inverse = inverse(source_columns)
orthogonals = {}
for permutation in permutations(range(4)):
    target_columns = tuple(
        tuple(F(SIGNS[permutation[column]][row]) for column in range(3))
        for row in range(3)
    )
    orthogonal = matmul(target_columns, source_inverse)
    orthogonals[permutation] = orthogonal
    check(matmul(transpose(orthogonal), orthogonal) == identity(3),
          f"label permutation {permutation} is orthogonal")
    for label, vector in enumerate(SIGNS):
        check(matvec(orthogonal, tuple(map(F, vector)))
              == tuple(map(F, SIGNS[permutation[label]])),
              f"permutation {permutation} maps tetrahedral vector {label+1}")
    for a, b in PAIRS:
        vector = tuple(F(SIGNS[b][axis] - SIGNS[a][axis]) for axis in range(3))
        mapped = matvec(orthogonal, vector)
        target = tuple(sorted((permutation[a], permutation[b])))
        check(dyad_coordinates(mapped)
              == tuple(8 * value for value in ROOT[target]),
              f"permutation {permutation} covariantly maps root {a+1}{b+1}")


# The invariant subspace of symmetric gradient tensors is exactly span{I}.
SYMMETRIC_BASIS = tuple(symmetric_matrix(tuple(F(i == coordinate)
                                                   for i in range(6)))
                        for coordinate in range(6))
invariance_equations = []
for orthogonal in orthogonals.values():
    for output_coordinate in range(6):
        row = []
        for basis in SYMMETRIC_BASIS:
            transformed = matmul(transpose(orthogonal),
                                 matmul(basis, orthogonal))
            difference = matadd(basis, matscale(-1, transformed))
            row.append(symmetric_coordinates(difference)[output_coordinate])
        invariance_equations.append(tuple(row))

check(rank(invariance_equations) == 5,
      "tetrahedral invariant symmetric-gradient subspace is one-dimensional")
A1 = (F(1), F(1), F(1), F(0), F(0), F(0))
E_BASIS = ((F(1), F(-1), F(0), F(0), F(0), F(0)),
           (F(1), F(1), F(-2), F(0), F(0), F(0)))
T2_BASIS = ((F(0), F(0), F(0), F(1), F(0), F(0)),
            (F(0), F(0), F(0), F(0), F(1), F(0)),
            (F(0), F(0), F(0), F(0), F(0), F(1)))
check(all(dot(A1, equation) == 0 for equation in invariance_equations),
      "identity tensor spans the invariant gradient line")
check(all(dot(A1, source) == 0 for source in E_BASIS),
      "both diagonal-traceless E sources annihilate a lumped gradient")


# Degree algebra and lumped/common-mode no-go.
ALL_STATES = tuple(product((-1, 1), repeat=4))
ICE = tuple(state for state in ALL_STATES if sum(state) == 0)
for state in ALL_STATES:
    degree = sum((1 - z) // 2 for z in state)
    pair_sum = sum(state[a] * state[b] for a, b in PAIRS)
    check(F((degree - 2) ** 2) == F(1) + F(pair_sum, 2),
          f"degree-pair identity for state {state}")
    check(sum(state) ** 2 == 4 * (degree - 2) ** 2,
          f"total signed mode is the lumped degree charge for state {state}")
check(len(ICE) == 6 and all(sum(state) == 0 for state in ICE),
      "the six ice states have zero lumped total charge")
check(all(dot(source, A1) == 0 for source in E_BASIS),
      "an S4 common-mode coefficient has exact E null")


# Pair-resolved elastance matches the source-off degree operator exactly.
Q = F(3)
Q2 = Q * Q
U_D = F(10)
E_MUTUAL = U_D / (2 * Q2)
E_SELF = F(2)
E_REFERENCE = U_D - 2 * Q2 * E_SELF
ELASTANCE = tuple(tuple(E_SELF if row == column else E_MUTUAL
                         for column in range(4)) for row in range(4))

check(E_SELF - E_MUTUAL > 0 and E_SELF + 3 * E_MUTUAL > 0,
      "example S4 elastance is positive on contrast and common modes")
CAPACITANCE = inverse(ELASTANCE)
check(matmul(ELASTANCE, CAPACITANCE) == identity(4),
      "grounded example capacitance and elastance are exact inverses")
check(Q2 * E_MUTUAL == U_D / 2,
      "off-diagonal elastance matches the U_d over two pair coefficient")

# A floating four-terminal capacitance has only a neutral-quotient inverse.
# It cannot define the off-ice energies used by the full operator theorem.
FLOATING_CAPACITANCE = tuple(tuple(F(3) if row == column else F(-1)
                                   for column in range(4))
                             for row in range(4))
check(rank(FLOATING_CAPACITANCE) == 3,
      "floating S4 capacitance has one common-potential zero mode")
try:
    inverse(FLOATING_CAPACITANCE)
except ValueError:
    floating_rejected = True
else:
    floating_rejected = False
check(floating_rejected,
      "singular floating capacitance is rejected as a full-domain inverse")
check(all(sum(state) == 0 for state in ICE)
      and any(sum(state) != 0 for state in ALL_STATES if state not in ICE),
      "neutral quotient covers ice but not the off-ice defect sectors")

for state in ALL_STATES:
    z = tuple(map(F, state))
    quadratic = F(1, 2) * Q2 * dot(z, matvec(ELASTANCE, z)) + E_REFERENCE
    degree = sum((1 - value) // 2 for value in state)
    check(quadratic == U_D * (degree - 2) ** 2,
          f"pair-resolved elastance reproduces degree energy for {state}")

# q_* Z solder makes bare X charged.  A compensating transfer preserves the
# total charge on the fixed-total-charge encoded pair.
CHARGE_TOTAL = tuple(tuple(F(((-1, -1), (-1, 1), (1, -1), (1, 1))[row][0]
                              + ((-1, -1), (-1, 1), (1, -1), (1, 1))[row][1])
                            if row == column else F(0)
                            for column in range(4)) for row in range(4))
BARE_X = ((F(0), F(0), F(1), F(0)),
          (F(0), F(0), F(0), F(1)),
          (F(1), F(0), F(0), F(0)),
          (F(0), F(1), F(0), F(0)))
DRESSED_X = ((F(0), F(0), F(0), F(0)),
             (F(0), F(0), F(1), F(0)),
             (F(0), F(1), F(0), F(0)),
             (F(0), F(0), F(0), F(0)))
bare_commutator = matadd(matmul(CHARGE_TOTAL, BARE_X),
                         matscale(-1, matmul(BARE_X, CHARGE_TOTAL)))
dressed_commutator = matadd(matmul(CHARGE_TOTAL, DRESSED_X),
                            matscale(-1, matmul(DRESSED_X, CHARGE_TOTAL)))
check(bare_commutator != tuple(tuple(F(0) for _ in range(4)) for _ in range(4)),
      "bare terminal X fails total-charge conservation")
check(dressed_commutator == tuple(tuple(F(0) for _ in range(4)) for _ in range(4)),
      "compensating terminal-reservoir transfer conserves total charge")
check(DRESSED_X[1][2] == DRESSED_X[2][1] == 1,
      "dressed transfer acts as Pauli X on the fixed-total-charge code pair")

# Exact derivative of inverse: D(C^-1)=-C^-1(DC)C^-1.
TEST_LAMBDA = F(3, 7)
test_source = E_BASIS[0]
d_elastance = [[F(0) for _ in range(4)] for _ in range(4)]
for a, b in PAIRS:
    derivative = -(U_D / (2 * Q2)) * TEST_LAMBDA * dot(test_source, ROOT[(a, b)])
    d_elastance[a][b] = derivative
    d_elastance[b][a] = derivative
d_elastance = tuple(tuple(row) for row in d_elastance)
d_capacitance = matscale(-1, matmul(CAPACITANCE,
                                   matmul(d_elastance, CAPACITANCE)))
recovered_d_elastance = matscale(-1, matmul(ELASTANCE,
                                            matmul(d_capacitance, ELASTANCE)))
check(recovered_d_elastance == d_elastance,
      "capacitance derivative recovers elastance derivative exactly")


# Affine derivative, central-kernel chain rule, and source normalization.
SOURCE_COORDINATES = tuple(tuple(F(index == coordinate) for index in range(6))
                           for coordinate in range(6))
for pair in PAIRS:
    a, b = pair
    root_vector = tuple(F(SIGNS[b][axis] - SIGNS[a][axis]) for axis in range(3))
    for source in SOURCE_COORDINATES:
        matrix = symmetric_matrix(source)
        direct = -dot(root_vector, matvec(matrix, root_vector)) / dot(root_vector, root_vector)
        check(direct == -dot(source, ROOT[pair]),
              f"F=I-j/2 affine ratio derivative for pair {a+1}{b+1}, source {source}")

R0 = F(5)
V0 = F(17)
V_PRIME = F(3)
GENERIC_LAMBDA = R0 * V_PRIME / (2 * V0)
check(GENERIC_LAMBDA == F(15, 34),
      "generic central-kernel slope formula is exact")

def generic_g(x):
    return (2.0 + 3.0 * float(R0) * sqrt(x)) / float(V0)

h = 1e-6
numeric_generic = (generic_g(1.0 + h) - generic_g(1.0 - h)) / (2 * h)
check(abs(numeric_generic - float(GENERIC_LAMBDA)) < 1e-9,
      "numerical central-kernel derivative matches exact lambda")

for pair in PAIRS:
    for source in SOURCE_COORDINATES:
        ratio_derivative = -dot(source, ROOT[pair])
        energy_derivative = (U_D / 2) * GENERIC_LAMBDA * ratio_derivative
        conjugate = -2 * energy_derivative
        expected = U_D * GENERIC_LAMBDA * dot(source, ROOT[pair])
        check(conjugate == expected,
              f"DPAR Q=-2 dH/dj sign and factor for pair {pair}, source {source}")


# Ideal Coulomb and rationalized alpha normalization.
COULOMB_LAMBDA = F(-1, 2)
check(COULOMB_LAMBDA == F(-1, 2), "ideal Coulomb symbolic lambda is minus one half")

def coulomb_g(x):
    return x ** -0.5

numeric_coulomb = (coulomb_g(1.0 + h) - coulomb_g(1.0 - h)) / (2 * h)
check(abs(numeric_coulomb + 0.5) < 1e-9,
      "numerical Coulomb derivative matches exact minus one half")

ALPHA = F(1, 137)
KAPPA = F(2)
EPSILON_R = F(3)
HBAR_C = F(197)
PHYSICAL_R0 = F(5)
PAIR_COEFFICIENT = KAPPA**2 * ALPHA * HBAR_C / (EPSILON_R * PHYSICAL_R0)
UD_COULOMB = 2 * KAPPA**2 * ALPHA * HBAR_C / (EPSILON_R * PHYSICAL_R0)
check(UD_COULOMB / 2 == PAIR_COEFFICIENT,
      "rationalized SI alpha normalization reproduces the pair coefficient")
check(UD_COULOMB == F(1576, 2055),
      "representative U_d-alpha-length normalization is exact rational arithmetic")

BETA_ALPHA = F(1, 10000)
RUNNING_LAMBDA = F(-1, 2) - BETA_ALPHA / (2 * ALPHA)
check(RUNNING_LAMBDA != COULOMB_LAMBDA,
      "running coupling changes the ideal fixed-coupling Coulomb slope")


# Source-rank composition and the projected ceiling.
ROOT_ROWS = tuple(ROOT[pair] for pair in PAIRS)
check(rank(EDGE + ROOT_ROWS) == 6,
      "unchanged edge source plus nonzero DPAR roots has microscopic rank six")
check(rank(EDGE + tuple(tuple(COULOMB_LAMBDA * value for value in row)
                        for row in ROOT_ROWS)) == 6,
      "Coulomb lambda minus one half retains exact microscopic rank six")
check(rank(EDGE + tuple(tuple(F(0) for _ in row) for row in ROOT_ROWS)) == 4,
      "zero pair slope leaves the inherited microscopic rank four")

def pair_query_values(source):
    coefficients = tuple(dot(source, ROOT[pair]) for pair in PAIRS)
    return tuple(sum(coefficient * state[a] * state[b]
                     for coefficient, (a, b) in zip(coefficients, PAIRS))
                 for state in ICE)

ice_all = tuple(pair_query_values(source) for source in SOURCE_COORDINATES)
ice_e = tuple(pair_query_values(source) for source in E_BASIS)
ice_t2 = tuple(pair_query_values(source) for source in T2_BASIS)
check(rank(ice_all) == 3, "direct ice pair map remains exact A1+E rank three")
check(rank(ice_e) == 2, "both centered ice E directions remain present")
check(all(all(value == 0 for value in row) for row in ice_t2),
      "direct ice pair map kills all T2 directions")


# Common-mode terms cannot cancel E; an equal and opposite mutual remainder can.
check(all(dot(source, A1) == 0 for source in E_BASIS),
      "uniform common-mode source cannot cancel a pair E response")
central_e_map = tuple(tuple(COULOMB_LAMBDA * dot(source, ROOT[pair])
                            for pair in PAIRS) for source in E_BASIS)
cancelled_e_map = tuple(tuple(value - value for value in row)
                        for row in central_e_map)
check(rank(central_e_map) == 2, "nonzero central slope has exact E rank two")
check(rank(cancelled_e_map) == 0,
      "an equal opposite local mutual remainder can physically cancel E")


# Cross-node Walsh operators are independent of the six local pair strings.
EIGHT_BIT_STATES = tuple(product((-1, 1), repeat=8))
characters = []
for a, b in PAIRS:
    characters.append(tuple(state[a] * state[b] for state in EIGHT_BIT_STATES))
for a in range(4):
    for b in range(4, 8):
        characters.append(tuple(state[a] * state[b] for state in EIGHT_BIT_STATES))
check(rank(characters) == 22,
      "six local and sixteen cross-node pair Walsh strings are independent")


# Documentary ceilings and byte hygiene.
theorem = (LANE / "THEOREM.md").read_text()
result = (LANE / "RESULT.md").read_text()
self_audit = (LANE / "SELF_AUDIT.md").read_text()
check("A1`-only first variation and exact\n`E` nullity two" in theorem,
      "theorem states the lumped A1-only E-null result")
check("new **physical completion**" in theorem
      or "new\n**physical completion**" in theorem,
      "theorem does not call the capacitance completion inherited F3")
check("lambda_E^net != 0" in theorem and "real physical cancellation" in theorem,
      "theorem makes complete-source noncancellation load bearing")
check("Cross-node Walsh operators" not in theorem
      and "Walsh strings are independent" in theorem,
      "theorem retains nonlocal operator independence without verifier leakage")
check("not identify `U_d` with\n`U_E^IR`" in theorem,
      "theorem keeps Gauss penalty distinct from infrared stiffness")
check("Direct restriction of the pair map to local\nice still has only `A1+E` rank three" in theorem,
      "theorem retains the projected rank-three ceiling")
check("S9 — noncircular ancestry" in theorem
      and "conditional physical completion only, not inherited F3" in theorem,
      "theorem makes the circularity ceiling decisive")
check("not an\ninherited F3 alpha calculation" in result,
      "result withholds inherited alpha promotion")
check("microscopic source has exact rank\nsix before Feshbach" in result,
      "result limits rank closure to the microscopic source")
check("neutral quotient may instead be used **only after restriction" in theorem
      and "cannot discharge (FU29)" in theorem,
      "theorem scopes a neutral quotient to the ice restriction")
check("bare inherited `X_a` is a gauge-invariant physical-charge flip" in theorem
      and "bare inherited `X_a`" in theorem,
      "theorem withholds bare-X charged-dynamics promotion")
check("[Q_{\\rm tot},\\widetilde X_a]=0" in theorem
      and "fixed-total-charge encoded subspace" in theorem,
      "theorem requires a charge-conserving dressed flip and exact code")
check("Was `U_d` relabeled as Maxwell electric stiffness?" in self_audit
      and "No." in self_audit,
      "self-audit explicitly tests the U1SI reconciliation")
check("Decisive falsifiers" in self_audit and "cancel the entire `E` slope" in self_audit,
      "self-audit includes physical failure conditions")
check("derives gravity" not in theorem.lower()
      and "calculates alpha" not in theorem.lower(),
      "theorem contains no gravity or alpha overpromotion phrase")

for name in MANIFEST_FILES:
    data = (LANE / name).read_bytes()
    check(b"\r" not in data and b"\b" not in data and b"\f" not in data,
          f"byte hygiene {name}")

verification_text = (
    "GRA_FU_PAIR_RESOLVED_MAXWELL_DPAR_DERIVATION: PASS\n"
    "Lumped tetrahedral charging: A1-only first variation, E null two\n"
    "Pair-resolved central elastance: DPAR derived conditionally\n"
    "Ideal fixed-coupling Coulomb: lambda=-1/2 and U_d-alpha-r0 normalization checked\n"
    "FT composition: exact microscopic rank six only under complete noncircular solder and noncancellation\n"
    "IR Maxwell stiffness, visible-U1 identity, projected CTP rank, gravity, and G: not claimed\n"
)
check((LANE / "VERIFICATION.txt").read_text() == verification_text,
      "verification transcript is exact")

manifest = {}
for line in (LANE / "MANIFEST.sha256").read_text().splitlines():
    expected, name = line.split("  ", 1)
    manifest[name] = expected
check(set(manifest) == MANIFEST_FILES, "manifest member set is exact")
check(all(digest(LANE / name) == expected for name, expected in manifest.items()),
      "manifest hashes match")
check(all(not (LANE / name).is_symlink() for name in manifest),
      "manifest members are not symlinks")

seal = {}
for line in (LANE / "SEAL.sha256").read_text().splitlines():
    expected, name = line.split("  ", 1)
    seal[name] = expected
check(set(seal) == {"MANIFEST.sha256", "VERIFICATION.txt"},
      "seal covers manifest and verification transcript")
check(all(digest(LANE / name) == expected for name, expected in seal.items()),
      "seal hashes match")
check(all(not (LANE / name).is_symlink() for name in seal),
      "seal members are not symlinks")

lane_files = {path.name for path in LANE.iterdir() if path.is_file()}
expected_files = (MANIFEST_FILES | AUDIT_FILES
                  | {"MANIFEST.sha256", "SEAL.sha256", "VERIFICATION.txt"})
check(lane_files == expected_files, "lane file set is exact")

print("GRA_FU_PAIR_RESOLVED_MAXWELL_DPAR_DERIVATION: PASS")
print(f"Checks: {checks}/{checks}")
print("Lumped tetrahedral charging: A1-only first variation, E null two")
print("Pair-resolved central elastance: DPAR derived conditionally")
print("Ideal fixed-coupling Coulomb: lambda=-1/2 and normalization checked")
print("FT composition: microscopic rank six under complete solder and noncancellation")
print("IR stiffness, visible U1, projected CTP rank, gravity, and G: not claimed")
