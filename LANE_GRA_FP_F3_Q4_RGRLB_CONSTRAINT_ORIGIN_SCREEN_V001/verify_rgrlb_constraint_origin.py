#!/usr/bin/env python3
"""Exact finite replay for GRA-FP-F3-Q4-RCOS-V001.

This verifier checks the inherited q4 ice constraint algebra, local S4
representation content, pair-relation boundary, finite incidence ranks, and
the continuum symbol/degree-of-freedom mismatch.  It does not test a
thermodynamic limit, a tensor pole, a microscopic metric solder, or gravity.
"""

from __future__ import annotations

import hashlib
import itertools
import math
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
THEOREM_BYTES = (HERE / "THEOREM.md").read_bytes()
THEOREM = THEOREM_BYTES.decode("utf-8")
PASSED = 0


DEPENDENCIES = {
    "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md":
        "00eba581b90fb9f0b25e3fad1362b055049824897433ba06e356ab9b1f6c76ec",
    "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md":
        "5e68e4a8c62ad89cff309781a5cb54e071092e82594e60c04e7992414bc18dbe",
    "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "a91caa20d16b0a1194333f9b51d96546a4ea24d55e23bf1f04c7d249641af8db",
    "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/THEOREM.md":
        "cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98",
    "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/INDEPENDENT_AUDIT.md":
        "c52eab9d701d1c6e82f1d7ec395841f4d2810e96cccbc3e2504760b6742e81e4",
    "LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/THEOREM.md":
        "98e2b3bc7a1c998d7839dc1a6b435cc1c8ed6d5a622ba45f63571be9ef646452",
    "LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/INDEPENDENT_AUDIT.md":
        "327bf6a4476c4c6382757dc156a96c6032233d34c25c1f7935e2582acf6c607a",
    "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md":
        "78f0687c9f597c96d235224dad45c204d12d7e6c973c270e3305a901efc75b25",
    "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/INDEPENDENT_AUDIT.md":
        "53893c7198241f0f8f6aa766f3595fb75b83d208581833c32656b28d7c7f02b9",
    "FINAL_GRAVITY_REAL_WORLD_THEOREM_V001.md":
        "1caabded24b861932b319ed715556a5d4123b2cff5ea3004676e12c4c76de155",
    "FINAL_GRAVITY_SAME_WORLD_COMPOSITION_GATE_MATRIX_V001.md":
        "8bde0b0a8bff8561b929d8a32412d6e3a51d810c89f16f163a1f35488a3772e6",
}


def check(condition: bool, label: str) -> None:
    global PASSED
    if not condition:
        raise AssertionError(label)
    PASSED += 1
    print(f"PASS {label}")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank(matrix: list[list[int | Fraction]]) -> int:
    """Exact rational row rank."""
    if not matrix:
        return 0
    a = [[Fraction(value) for value in row] for row in matrix]
    rows = len(a)
    cols = len(a[0])
    pivot_row = 0
    for col in range(cols):
        pivot = next((r for r in range(pivot_row, rows) if a[r][col]), None)
        if pivot is None:
            continue
        a[pivot_row], a[pivot] = a[pivot], a[pivot_row]
        scale = a[pivot_row][col]
        a[pivot_row] = [value / scale for value in a[pivot_row]]
        for row in range(rows):
            if row != pivot_row and a[row][col]:
                factor = a[row][col]
                a[row] = [left - factor * right for left, right in zip(a[row], a[pivot_row])]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def matmul(left: list[list[int | Fraction]], right: list[list[int | Fraction]]) -> list[list[Fraction]]:
    return [
        [sum(Fraction(a) * Fraction(b) for a, b in zip(row, col)) for col in zip(*right)]
        for row in left
    ]


def matvec(matrix: list[list[int]], vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in matrix)


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

required_text = (
    "PROJECTED_ICE_HAS_ONE_SCALAR_U1_GAUSS_SPECIES",
    "PAIR_RELATIONS_ARE_ALGEBRAIC_NOT_NEW_FIRST_CLASS_GENERATORS",
    "S4_MODULE_MATCH_IS_NOT_CONTINUUM_SPIN2_EQUIVALENCE",
    "RGRLB_NOT_MICROSCOPICALLY_DERIVED_FROM_CURRENT_Q4_ICE_BRANCH",
    "FINITE_SCREEN_NOT_THERMODYNAMIC_NO_GO",
)
for token in required_text:
    check(token in THEOREM, f"theorem binds disposition token: {token}")


# Six local ice states and their one-link/pair function ranks.
ice = [state for state in itertools.product((-1, 1), repeat=4) if sum(state) == 0]
edges = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
edge_index = {edge: index for index, edge in enumerate(edges)}
one = [list(state) for state in ice]
pairs = [[state[a] * state[b] for a, b in edges] for state in ice]
centered_pairs = [[Fraction(value) + Fraction(1, 3) for value in row] for row in pairs]

check(len(ice) == 6, "q4 d*=2 local ice fiber has six states")
check(all(sum(state) == 0 for state in ice), "local scalar ice constraint holds")
check(rank(one) == 3, "one-link nonconstant module has rank three")
check(rank(pairs) == 3, "pair module has rank three including its constant")
check(rank(centered_pairs) == 2, "centered pair tangent has rank two")
check(rank([[1] + one_row + pair_row for one_row, pair_row in zip(one, centered_pairs)]) == 6,
      "constant plus one-link plus centered-pair functions saturate six-state diagonal algebra")

opposites = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))
check(
    all(
        all(row[edge_index[first]] == row[edge_index[second]] for row in pairs)
        for first, second in opposites
    ),
    "three opposite-pair identities hold on every ice state",
)
check(all(sum(row) == -2 for row in pairs), "uniform pair sum is fixed at minus two")

# Four affine relations on the six formal pair coordinates.
relations: list[list[int]] = []
for first, second in opposites:
    row = [0] * 7
    row[1 + edge_index[first]] = 1
    row[1 + edge_index[second]] = -1
    relations.append(row)
relations.append([2] + [1] * 6)
augmented_pairs = [[1] + row for row in pairs]
relation_values = matmul(augmented_pairs, [list(column) for column in zip(*relations)])
check(all(all(value == 0 for value in row) for row in relation_values),
      "all four affine pair relations vanish identically after ice projection")
check(rank([row[1:] for row in relations]) == 4,
      "pair-relation coefficient gradients have formal rank four")
check(rank(augmented_pairs) == 3,
      "six formal pairs restricted by four affine relations leave A1 plus E rank three")


# S4 characters in class order 1, (12), (12)(34), (123), (1234).
representatives = (
    (0, 1, 2, 3),
    (1, 0, 2, 3),
    (1, 0, 3, 2),
    (1, 2, 0, 3),
    (1, 2, 3, 0),
)


def fixed_points(permutation: tuple[int, ...]) -> int:
    return sum(index == moved for index, moved in enumerate(permutation))


def permute_state(state: tuple[int, ...], permutation: tuple[int, ...]) -> tuple[int, ...]:
    moved = [0] * 4
    for old, value in enumerate(state):
        moved[permutation[old]] = value
    return tuple(moved)


def compose(first: tuple[int, ...], second: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(first[second[index]] for index in range(4))


matchings = (
    frozenset((frozenset((0, 1)), frozenset((2, 3)))),
    frozenset((frozenset((0, 2)), frozenset((1, 3)))),
    frozenset((frozenset((0, 3)), frozenset((1, 2)))),
)


def permute_matching(matching: frozenset[frozenset[int]], permutation: tuple[int, ...]) -> frozenset[frozenset[int]]:
    return frozenset(frozenset(permutation[index] for index in pair) for pair in matching)


six_state_character = tuple(
    sum(permute_state(state, permutation) == state for state in ice)
    for permutation in representatives
)
vector_character = tuple(fixed_points(permutation) - 1 for permutation in representatives)
pair_character = tuple(
    sum(permute_matching(matching, permutation) == matching for matching in matchings)
    for permutation in representatives
)
formal_pair_character = tuple(
    sum(
        frozenset((permutation[a], permutation[b])) == frozenset((a, b))
        for a, b in edges
    )
    for permutation in representatives
)
sym2_character = tuple(
    (chi * chi + (fixed_points(compose(permutation, permutation)) - 1)) // 2
    for chi, permutation in zip(vector_character, representatives)
)

check(six_state_character == (6, 2, 2, 0, 0), "six-state permutation character is A1+E+T2")
check(vector_character == (3, 1, -1, 0, -1), "one-link standard character is T2")
check(pair_character == (3, 1, 3, 0, 1), "three matching character is A1+E")
check(tuple(value - 1 for value in pair_character) == (2, 0, 2, -1, 0),
      "centered pair character is E")
check(formal_pair_character == (6, 2, 2, 0, 0),
      "six formal pair coordinates carry A1+E+T2")
check(sym2_character == six_state_character, "Sym2(V) and ice functions coincide only as S4 modules")
check(tuple(1 + value for value in vector_character) == (4, 2, 0, 1, 0),
      "formal scalar-plus-vector constraint packet has A1+T2 character")
relation_character = tuple(
    1 + formal - image
    for formal, image in zip(formal_pair_character, pair_character)
)
check(relation_character == tuple(1 + value for value in vector_character),
      "four affine pair relations themselves carry A1+T2")


# Connected periodic coordination-four bipartite graph incidence ranks.
def periodic_diamond_incidence(length: int) -> tuple[list[list[int]], int, int]:
    coordinates = list(itertools.product(range(length), repeat=3))
    vertices = [(side,) + coordinate for side in (0, 1) for coordinate in coordinates]
    vertex_index = {vertex: index for index, vertex in enumerate(vertices)}
    offsets = ((0, 0, 0), (-1, 0, 0), (0, -1, 0), (0, 0, -1))
    graph_edges = []
    for coordinate in coordinates:
        for offset in offsets:
            target = tuple((coordinate[i] + offset[i]) % length for i in range(3))
            graph_edges.append(((0,) + coordinate, (1,) + target))
    incidence = [[0 for _ in graph_edges] for _ in vertices]
    for column, (tail, head) in enumerate(graph_edges):
        incidence[vertex_index[tail]][column] = 1
        incidence[vertex_index[head]][column] = -1
    return incidence, len(vertices), len(graph_edges)


for length, expected_rank in ((2, 15), (3, 53), (4, 127)):
    incidence, vertex_count, edge_count = periodic_diamond_incidence(length)
    check(vertex_count == 2 * length**3 and edge_count == 4 * length**3,
          f"L={length} periodic q4 graph has expected V/E counts")
    check(rank(incidence) == expected_rank == vertex_count - 1,
          f"L={length} connected incidence has exactly one global dependency")
    check(all(sum(row[column] for row in incidence) == 0 for column in range(edge_count)),
          f"L={length} oriented incidence columns sum to zero")


# A directed six-cycle models the exact local closed-ring Gauss algebra.
cycle_incidence = [[0 for _ in range(6)] for _ in range(6)]
for edge in range(6):
    cycle_incidence[edge][edge] = 1
    cycle_incidence[(edge + 1) % 6][edge] = -1
plus_circulation = (1, 1, 1, 1, 1, 1)
minus_circulation = tuple(-value for value in plus_circulation)
single_flip = (-1, 1, 1, 1, 1, 1)

check(matvec(cycle_incidence, plus_circulation) == (0,) * 6,
      "positive closed circulation satisfies every local Gauss sum")
check(matvec(cycle_incidence, minus_circulation) == (0,) * 6,
      "ring-reversed circulation satisfies every local Gauss sum")
check(matvec(cycle_incidence, single_flip) != (0,) * 6,
      "one microscopic link flip creates a nonzero Gauss-defect pair")
check(rank(cycle_incidence) == 5, "six-cycle incidence has one global Gauss dependency")
ring_swap = [[0, 1], [1, 0]]
zero_generator = [[0, 0], [0, 0]]
check(matmul(zero_generator, ring_swap) == matmul(ring_swap, zero_generator),
      "projected ring swap commutes with all restricted Gauss generators")


# Exact continuum principal-symbol discriminator at one nonzero momentum.
kx, ky, kz = 1, 2, 3
k2 = kx * kx + ky * ky + kz * kz
maxwell_gauss = [[kx, ky, kz]]
tensor_vector_constraint = [
    [kx, 0, 0, ky, kz, 0],
    [0, ky, 0, kx, 0, kz],
    [0, 0, kz, 0, kx, ky],
]
tensor_scalar_constraint = [[
    kx * kx - k2,
    ky * ky - k2,
    kz * kz - k2,
    2 * kx * ky,
    2 * kx * kz,
    2 * ky * kz,
]]
tensor_gauge_map = [
    [2 * kx, 0, 0],
    [0, 2 * ky, 0],
    [0, 0, 2 * kz],
    [ky, kx, 0],
    [kz, 0, kx],
    [0, kz, ky],
]
tt_conditions = tensor_vector_constraint + [[1, 1, 1, 0, 0, 0]]
combined_phase_constraints = [row + [0] * 6 for row in tensor_scalar_constraint]
combined_phase_constraints += [[0] * 6 + row for row in tensor_vector_constraint]

check(rank(maxwell_gauss) == 1, "Maxwell Gauss symbol has one scalar constraint at nonzero k")
check(3 - rank(maxwell_gauss) == 2, "Maxwell transverse configuration space has dimension two")
check(rank(tensor_vector_constraint) == 3, "rank-two momentum-constraint symbol has vector rank three")
check(rank(tensor_scalar_constraint) == 1, "rank-two Hamiltonian/scalar symbol has independent rank one")
check(rank(combined_phase_constraints) == 4, "rank-two spatial phase-space constraint packet has rank four")
check(rank(tensor_gauge_map) == 3, "spatial symmetric-tensor gauge map has three independent vector parameters")
check(all(value == 0 for row in matmul(tensor_scalar_constraint, tensor_gauge_map) for value in row),
      "linear scalar curvature row annihilates spatial gauge directions")
check(rank(tt_conditions) == 4 and 6 - rank(tt_conditions) == 2,
      "transverse-tracefree symmetric tensor space has dimension two")
check(6 - 2 * rank(maxwell_gauss) == 4,
      "one Maxwell first-class constraint leaves four phase dimensions")
check(12 - 2 * rank(combined_phase_constraints) == 4,
      "four rank-two first-class constraints leave four phase dimensions")

# One scalar Fourier amplitude cannot become four independent lapse/shift amplitudes.
scalar_fourier_jet = [[1], [kx], [ky], [kz]]
independent_scalar_vector_packet = [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
]
check(rank(scalar_fourier_jet) == 1, "one scalar Fourier gauge amplitude and its gradient have rank one")
check(rank(independent_scalar_vector_packet) == 4,
      "independent scalar-plus-vector gauge amplitudes have rank four")

# Continuum helicity characters distinguish the S4-isomorphic T2 copies.
theta = math.pi / 2
photon_character = 2 * math.cos(theta)
tt_character = 2 * math.cos(2 * theta)
check(abs(photon_character) < 1e-12, "helicity-one character vanishes at pi/2")
check(abs(tt_character + 2) < 1e-12, "helicity-two character equals minus two at pi/2")
check(abs(photon_character - tt_character) > 1,
      "continuum rotation characters separate vector and TT channels")


print(f"SUMMARY {PASSED}/{PASSED} exact checks passed")
print("DISPOSITION PROJECTED_ICE_ONE_SCALAR_U1_ONLY__NO_INHERITED_RGRLB_RANK2_CONSTRAINT_PACKET")
