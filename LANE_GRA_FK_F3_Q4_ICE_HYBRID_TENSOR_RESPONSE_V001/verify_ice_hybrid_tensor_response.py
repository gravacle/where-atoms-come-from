#!/usr/bin/env python3
"""Finite exact replay for GRA-FK-F3-Q4-IHTR-V001.

The verifier checks the local ice representation, the abstract hybrid
Sym^2(V) representation rank, the symmetric-ice Fisher-query ceiling, the
inherited ring-flip parity algebra, and the compressed linked-ring response.
It does not test a physical metric solder, a thermodynamic phase, or gravity.
"""

from __future__ import annotations

import hashlib
import itertools
from fractions import Fraction
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
THEOREM_BYTES = (HERE / "THEOREM.md").read_bytes()
THEOREM = THEOREM_BYTES.decode("utf-8")
PASSED = 0


DEPENDENCIES = {
    "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md":
        "5e68e4a8c62ad89cff309781a5cb54e071092e82594e60c04e7992414bc18dbe",
    "LANE_GRA_FE_F3_Q4_DIAMOND_ICE_CARRIER_JOIN_V001/THEOREM.md":
        "4cc63e3e5853b4250a2a5b78256d41b83b195cf527901819a62a04ef53f8d932",
    "LANE_GRA_FH_F3_Q4_FINITE_PROGRAMMED_SUPPORT_SOLDER_V001/THEOREM.md":
        "2b88febc569efa0de0238e8000d018bf3f798a8ebed2e4ff1327f053d6bd9284",
    "LANE_GRA_FH_F3_Q4_FINITE_PROGRAMMED_SUPPORT_SOLDER_V001/INDEPENDENT_REAUDIT.md":
        "5c275748d54743ef44098f74c4c5698aead0845d51e6c2dcf32a1bef63f0c7bf",
    "LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/THEOREM.md":
        "05f4a619a6f80aa40c48570ab4035ab874426502a31468a08f435e66610bd769",
    "LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/INDEPENDENT_AUDIT.md":
        "44690ad431c85af7a4947a431a4c57ad4cd8b19a346e3d26720b341f77256f90",
    "LANE_CROSS_RFT_GRA_EW_Q4_PAIR_MEMORY_METRIC_DEFORMATION_CLOSURE_V001/THEOREM.md":
        "495e4e99171f4e3e5809f24e5a9a5b68116e996f4f29115bd22f780127d4714e",
}


def check(condition: bool, label: str) -> None:
    global PASSED
    if not condition:
        raise AssertionError(label)
    PASSED += 1
    print(f"PASS {label}")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def orthonormal_columns(matrix: np.ndarray, tol: float = 1e-11) -> np.ndarray:
    u, singular, _ = np.linalg.svd(matrix, full_matrices=False)
    return u[:, singular > tol]


def symmetric_vector(matrix: np.ndarray) -> np.ndarray:
    return np.array([
        matrix[0, 0], matrix[1, 1], matrix[2, 2],
        np.sqrt(2) * matrix[0, 1],
        np.sqrt(2) * matrix[0, 2],
        np.sqrt(2) * matrix[1, 2],
    ])


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


# Local ice states and observables.
ice = np.array([
    state for state in itertools.product((-1.0, 1.0), repeat=4)
    if sum(state) == 0
])
edges = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
edge_index = {edge: index for index, edge in enumerate(edges)}
one = ice.copy()
pairs = np.array([[state[a] * state[b] for a, b in edges] for state in ice])
constant = np.ones((len(ice), 1))

check(ice.shape == (6, 4), "local d*=2 ice fiber has six states")
check(np.allclose(one.sum(axis=1), 0), "one-link ice sum vanishes")
check(np.linalg.matrix_rank(one) == 3, "one-link span has rank three")
check(np.allclose(pairs.sum(axis=1), -2), "uniform pair sum is fixed at -2")

opposites = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))
check(
    all(np.array_equal(pairs[:, edge_index[a]], pairs[:, edge_index[b]]) for a, b in opposites),
    "opposite pair observables are identical in ice",
)
check(np.linalg.matrix_rank(pairs) == 3, "pair span has rank three A1+E")
centered_pairs = pairs - pairs.mean(axis=0, keepdims=True)
check(np.allclose(pairs.mean(axis=0), -1 / 3), "every pair mean is -1/3 at the uniform ice point")
check(np.linalg.matrix_rank(centered_pairs) == 2, "normalized pair tangent has rank two E only")
check(
    np.linalg.matrix_rank(np.column_stack((constant, one, centered_pairs))) == 6,
    "constant plus one-link and centered-pair spans every diagonal ice function",
)


# Characters of the full, one-link, pair, and centered-pair modules.
representatives = (
    (0, 1, 2, 3),
    (1, 0, 2, 3),
    (1, 0, 3, 2),
    (1, 2, 0, 3),
    (1, 2, 3, 0),
)


def state_action(permutation: tuple[int, ...]) -> np.ndarray:
    result = np.zeros((6, 6))
    lookup = {tuple(state): index for index, state in enumerate(ice)}
    for column, state in enumerate(ice):
        moved = [0.0] * 4
        for old in range(4):
            moved[permutation[old]] = state[old]
        result[lookup[tuple(moved)], column] = 1.0
    return result


def character(functions: np.ndarray) -> tuple[int, ...]:
    basis = orthonormal_columns(functions)
    values = []
    for permutation in representatives:
        trace = np.trace(basis.T @ state_action(permutation) @ basis)
        values.append(int(round(float(trace))))
    return tuple(values)


check(character(np.eye(6)) == (6, 2, 2, 0, 0), "six-state module character")
check(character(one) == (3, 1, -1, 0, -1), "one-link module is T2")
check(character(pairs) == (3, 1, 3, 0, 1), "pair module is A1+E")
check(character(centered_pairs) == (2, 0, 2, -1, 0), "centered pair module is E")


# The five nonconstant statistics are a saturated normalized family.
matchings = np.column_stack((
    pairs[:, edge_index[(0, 1)]],
    pairs[:, edge_index[(0, 2)]],
    pairs[:, edge_index[(0, 3)]],
))
statistics = np.column_stack((one[:, :3], matchings[:, 0] - matchings[:, 2], matchings[:, 1] - matchings[:, 2]))
statistics -= statistics.mean(axis=0, keepdims=True)
covariance = statistics.T @ statistics / 6
check(np.linalg.matrix_rank(statistics) == 5, "ice exponential-family statistic rank is five")
check(np.all(np.linalg.eigvalsh(covariance) > 1e-10), "uniform ice Fisher covariance is positive definite")
check(np.linalg.matrix_rank(np.column_stack((constant, statistics))) == 6, "ice family is locally saturated")


# Hybrid equivariant map to Sym^2(V).
p4 = np.eye(4) - np.ones((4, 4)) / 4
eigvals, eigvecs = np.linalg.eigh(p4)
v_basis = eigvecs[:, eigvals > 0.5]
tetra = [p4[:, a] for a in range(4)]

one_coefficients = (
    np.array([1.0, 0.0, 0.0, -1.0]),
    np.array([0.0, 1.0, 0.0, -1.0]),
    np.array([0.0, 0.0, 1.0, -1.0]),
)
one_tensors = []
for coefficient in one_coefficients:
    tensor4 = p4 @ np.diag(coefficient) @ p4
    one_tensors.append(v_basis.T @ tensor4 @ v_basis)

e_vectors = []
for matching_weights in ((1.0, -1.0, 0.0), (1.0, 1.0, -2.0)):
    vector = np.zeros(6)
    for weight, (first, second) in zip(matching_weights, opposites):
        vector[edge_index[first]] = weight
        vector[edge_index[second]] = weight
    e_vectors.append(vector)

pair_tensors = []
for coefficients in e_vectors:
    tensor4 = np.zeros((4, 4))
    for coefficient, (a, b) in zip(coefficients, edges):
        tensor4 += coefficient * (
            np.outer(tetra[a], tetra[b]) + np.outer(tetra[b], tetra[a])
        )
    pair_tensors.append(v_basis.T @ tensor4 @ v_basis)

check(np.linalg.matrix_rank(np.column_stack([symmetric_vector(t) for t in one_tensors])) == 3,
      "one-link T2 map has rank three")
check(all(abs(np.trace(t)) < 1e-10 for t in one_tensors), "one-link T2 tensors are trace-free")
check(np.linalg.matrix_rank(np.column_stack([symmetric_vector(t) for t in pair_tensors])) == 2,
      "ice pair E map has rank two")
check(all(abs(np.trace(t)) < 1e-10 for t in pair_tensors), "ice pair E tensors are trace-free")

hybrid_columns = [symmetric_vector(np.eye(3))]
hybrid_columns.extend(symmetric_vector(t) for t in pair_tensors)
hybrid_columns.extend(symmetric_vector(t) for t in one_tensors)
check(np.linalg.matrix_rank(np.column_stack(hybrid_columns)) == 6,
      "abstract A1 plus pair E plus one-link T2 candidate spans Sym2(V)")

# Exact physical ceiling for the most direct symmetric ice Fisher query.
fisher0_4 = one.T @ one / len(ice)
check(np.allclose(fisher0_4, (4 / 3) * p4),
      "uniform ice Fisher covariance is exactly (4/3) P")

theta_fisher_derivatives = []
for coefficient in one_coefficients:
    score = one @ coefficient
    derivative4 = sum(
        np.outer(state, state) * weight for state, weight in zip(ice, score)
    ) / len(ice)
    theta_fisher_derivatives.append(v_basis.T @ derivative4 @ v_basis)
check(all(np.allclose(derivative, 0) for derivative in theta_fisher_derivatives),
      "symmetric ice Fisher query has zero first-order one-link T2 derivative")
check(
    all(
        sum(int(state[a] * state[b] * state[c]) for state in ice) == 0
        for a in range(4) for b in range(4) for c in range(4)
    ),
    "all cubic ice moments vanish exactly by global flip",
)

pair_fisher_derivatives = []
for coefficients, target in zip(e_vectors, pair_tensors):
    score = pairs @ coefficients
    derivative4 = sum(
        np.outer(state, state) * weight for state, weight in zip(ice, score)
    ) / len(ice)
    derivative_v = v_basis.T @ derivative4 @ v_basis
    pair_fisher_derivatives.append(derivative_v)
    if not np.allclose(derivative_v, (8 / 3) * target):
        raise AssertionError("pair E Fisher derivative normalization")
check(True, "symmetric ice pair-E Fisher derivative is exactly (8/3) M(y)")
check(
    np.linalg.matrix_rank(
        np.column_stack([symmetric_vector(t) for t in pair_fisher_derivatives])
    ) == 2,
    "symmetric ice first-order Fisher-metric tangent has E rank two only",
)

# A uniform pair coupling is only a normalization constant inside ice.
uniform_pair_score = pairs @ np.ones(6)
check(np.allclose(uniform_pair_score, -2), "uniform pair coupling is constant in ice")
check(np.allclose(uniform_pair_score - uniform_pair_score.mean(), 0),
      "uniform pair coupling supplies no normalized A1 tangent")


# Local ring moves span every surviving nonconstant representation sector.
one_differences = []
pair_differences = []
for state in ice:
    for p, q in itertools.combinations(range(4), 2):
        if state[p] == state[q]:
            continue
        moved = state.copy()
        moved[p] *= -1
        moved[q] *= -1
        before_pair = np.array([state[a] * state[b] for a, b in edges])
        after_pair = np.array([moved[a] * moved[b] for a, b in edges])
        one_differences.append(moved - state)
        pair_differences.append(after_pair - before_pair)

check(np.linalg.matrix_rank(np.array(one_differences)) == 3,
      "allowed local ring moves span the one-link T2 directions")
check(np.linalg.matrix_rank(np.array(pair_differences)) == 2,
      "allowed local ring moves span the pair E directions")
diamond_ring_pairs = set()
for label_triple in itertools.combinations(range(4), 3):
    diamond_ring_pairs.update(itertools.combinations(label_triple, 2))
check(
    diamond_ring_pairs == set(itertools.combinations(range(4), 2)),
    "diamond three-label hexagons expose all six local ring-edge pairs",
)


# Independently re-sum the inherited symmetric sixth-order coefficient.
sigma = tuple(-1 if index % 2 == 0 else 1 for index in range(6))


def virtual_gap(selected: tuple[int, ...]) -> int:
    chosen = set(selected)
    total = 0
    for vertex in range(6):
        charge = sum(sigma[edge] for edge in ((vertex - 1) % 6, vertex) if edge in chosen)
        total += charge * charge
    return total


j6_reduced = Fraction(0, 1)
for order in itertools.permutations(range(6)):
    denominator = 1
    for length in range(1, 6):
        gap = virtual_gap(order[:length])
        if gap <= 0:
            raise AssertionError("proper ring prefix returned to ice")
        denominator *= gap
    j6_reduced += Fraction(1, denominator)
check(j6_reduced == Fraction(63, 8), "inherited symmetric J6 coefficient is 63/8")


# Exact parity algebra on one six-link flippable ring.
dimension = 64
alternating_0 = sum(1 << edge for edge in range(0, 6, 2))
alternating_1 = ((1 << 6) - 1) ^ alternating_0
b_ring = np.zeros((dimension, dimension), dtype=complex)
b_ring[alternating_0, alternating_1] = 1
b_ring[alternating_1, alternating_0] = 1


def walsh(mask: int) -> np.ndarray:
    diagonal = [(-1) ** bin(state & mask).count("1") for state in range(dimension)]
    return np.diag(diagonal).astype(complex)


for mask in range(64):
    parity = bin(mask).count("1") % 2
    lhs = b_ring @ walsh(mask)
    rhs = ((-1) ** parity) * walsh(mask) @ b_ring
    if not np.allclose(lhs, rhs):
        raise AssertionError(f"ring parity failed for mask {mask}")
check(True, "ring flip obeys exact Walsh intersection parity for all 64 subsets")

h = 0.31
u_degree = 1.7
j6 = 63 * h**6 / (8 * u_degree**5)
h_ring = -j6 * b_ring
odd = walsh(1 << 0)
even = walsh((1 << 0) | (1 << 1))
check(np.allclose(h_ring @ odd - odd @ h_ring, 2 * j6 * odd @ b_ring),
      "odd Walsh commutator is 2 J6 W B")
check(np.allclose(h_ring @ even - even @ h_ring, 0),
      "even Walsh observable commutes with the ring term")

second_odd = walsh(1 << 2)
first_commutator = h_ring @ odd - odd @ h_ring
double_commutator = first_commutator @ second_odd - second_odd @ first_commutator
check(np.allclose(double_commutator, -4 * j6 * odd @ second_odd @ b_ring),
      "two odd observables have exact -4 J6 A D B nested commutator")
check(np.linalg.norm(double_commutator) > 1e-12,
      "common-ring cross-observable operator response is nonzero")


# Exact compressed two-state response.
indices = (alternating_0, alternating_1)
h2 = h_ring[np.ix_(indices, indices)]
a2 = odd[np.ix_(indices, indices)]
d2 = second_odd[np.ix_(indices, indices)]
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
check(np.allclose(h2, -j6 * sigma_x), "compressed ring Hamiltonian is -J6 sigma_x")
check(np.allclose(abs(a2), abs(sigma_z)) and np.allclose(abs(d2), abs(sigma_z)),
      "odd one-link/crossing-pair observables compress to signed sigma_z")


def spectral_response(hamiltonian: np.ndarray, a_op: np.ndarray, d_op: np.ndarray, z: complex) -> complex:
    energies, vectors = np.linalg.eigh(hamiltonian)
    ground = vectors[:, 0]
    e0 = energies[0]
    total = 0.0 + 0.0j
    for index in range(1, len(energies)):
        excited = vectors[:, index]
        gap = energies[index] - e0
        total += np.vdot(ground, a_op @ excited) * np.vdot(excited, d_op @ ground) / (z - gap)
        total -= np.vdot(ground, d_op @ excited) * np.vdot(excited, a_op @ ground) / (z + gap)
    return total


kappa = 0.19
z = 1j * kappa
numeric_response = spectral_response(h2, a2, d2, z)
o_a = float(np.real(a2[0, 0]))
o_d = float(np.real(d2[0, 0]))
analytic_response = o_a * o_d * (1 / (z - 2 * j6) - 1 / (z + 2 * j6))
check(np.allclose(numeric_response, analytic_response, atol=1e-12),
      "compressed cross-vertex response is o_A o_D R_(2J6)")
check(abs(numeric_response) > 1e-12, "compressed common-ring kernel is nonzero")
check(np.isclose((1 / (z - 2 * j6) - 1 / (z + 2 * j6)).real,
                 -4 * j6 / (kappa**2 + 4 * j6**2)),
      "imaginary-axis ring response has the stated finite-gap value")


# Claim and type ceilings.
required_phrases = (
    "six-mode solder does not survive ice",
    "independently owned scalar",
    "No new microscopic interaction",
    "representation-isomorphism candidate",
    "first-order `T2` no-go",
    "physical metric rank is not",
    "global degree-two ice sector is empty",
    "supplied compatible",
    "explicit flippability premise",
    "full many-ring",
    "massless",
    "spin-one ice photon",
    "universal stress coupling",
    "gravity",
)
for phrase in required_phrases:
    check(phrase in THEOREM, f"claim ceiling present: {phrase}")

check(
    sum(line.strip() == r"\[" for line in THEOREM.splitlines())
    == sum(line.strip() == r"\]" for line in THEOREM.splitlines()),
    "display-math delimiters are balanced",
)

print(f"SUMMARY {PASSED}/{PASSED} PASS")
print("VERDICT ICE_MODULE_AND_ABSTRACT_HYBRID_REPRESENTATION_PASS__SYMMETRIC_FISHER_PAIR_E_PASS_BUT_T2_FIRST_ORDER_ZERO__RING_RESPONSE_PASS__PHYSICAL_METRIC_MASSLESS_TENSOR_AND_GRAVITY_OPEN")
