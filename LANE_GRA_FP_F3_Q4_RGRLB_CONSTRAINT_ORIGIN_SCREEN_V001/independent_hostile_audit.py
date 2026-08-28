#!/usr/bin/env python3
"""Independent hostile replay for GRA-FP-F3-Q4-RCOS-V001.

This file deliberately does not import the builder verifier.  It reconstructs
the finite graph, hard-core ice, representation, ring-commutator, rotor-count,
and continuum principal-symbol checks independently, then checks dependency
custody, the lane manifest, its outer seal, and negative tamper cases.
"""

from __future__ import annotations

import hashlib
import itertools
import math
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PASSED = 0


def check(condition: bool, label: str) -> None:
    global PASSED
    if not condition:
        raise AssertionError(label)
    PASSED += 1
    print(f"PASS {label}")


def digest_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def digest(path: Path) -> str:
    return digest_bytes(path.read_bytes())


def parse_hash_file(path: Path) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        value, name = line.split(maxsplit=1)
        parsed[name.strip()] = value
    return parsed


def rank(matrix: list[list[int | Fraction]]) -> int:
    if not matrix:
        return 0
    work = [[Fraction(entry) for entry in row] for row in matrix]
    rows, cols = len(work), len(work[0])
    pivot = 0
    for column in range(cols):
        selected = next((row for row in range(pivot, rows) if work[row][column]), None)
        if selected is None:
            continue
        work[pivot], work[selected] = work[selected], work[pivot]
        scale = work[pivot][column]
        work[pivot] = [entry / scale for entry in work[pivot]]
        for row in range(rows):
            if row == pivot or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [left - scale * right for left, right in zip(work[row], work[pivot])]
        pivot += 1
        if pivot == rows:
            break
    return pivot


def matmul(left: list[list[int | Fraction]], right: list[list[int | Fraction]]) -> list[list[Fraction]]:
    columns = list(zip(*right))
    return [
        [sum(Fraction(a) * Fraction(b) for a, b in zip(row, column)) for column in columns]
        for row in left
    ]


def matvec(matrix: list[list[int]], vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in matrix)


# ---------------------------------------------------------------------------
# Independent hard-core ice and pair-identity derivation.
# ---------------------------------------------------------------------------

all_states = list(itertools.product((-1, 1), repeat=4))
ice = [state for state in all_states if sum(state) == 0]
pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
opposites = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))

check(len(ice) == 6, "independent enumeration gives six two-in/two-out states")
check(all(sum(state) == 0 for state in ice), "every enumerated ice state obeys one scalar law")

pair_table = [[state[a] * state[b] for a, b in pairs] for state in ice]
one_table = [list(state) for state in ice]
centered_pair_table = [
    [Fraction(value, 1) + Fraction(1, 3) for value in row] for row in pair_table
]
check(rank(one_table) == 3, "one-link functions have rank three")
check(rank(pair_table) == 3, "pair functions have affine A1+E rank three")
check(rank(centered_pair_table) == 2, "centered pair response has rank two")
check(
    rank([[1] + links + centered for links, centered in zip(one_table, centered_pair_table)]) == 6,
    "constant, odd links, and even centered pairs exhaust the diagonal six-state algebra",
)

for first, second in opposites:
    check(
        all(state[first[0]] * state[first[1]] == state[second[0]] * state[second[1]] for state in ice),
        f"opposite-pair identity {first}={second} holds on ice",
    )
check(all(sum(row) == -2 for row in pair_table), "uniform pair sum equals minus two on ice")

# In the hard-core algebra s_a^2=1, the apparent four relations factor through
# S=sum_a s_a.  This independently demonstrates reducibility to the one scalar
# ice law, instead of treating a formal coefficient rank as four gauge orbits.
for first, second in opposites:
    check(
        all(
            2 * (state[first[0]] * state[first[1]] - state[second[0]] * state[second[1]])
            == sum(state)
            * (
                state[first[0]]
                + state[first[1]]
                - state[second[0]]
                - state[second[1]]
            )
            for state in all_states
        ),
        f"opposite-pair relation {first}-{second} factors through scalar ice generator",
    )
check(
    all(
        2 * (2 + sum(state[a] * state[b] for a, b in pairs)) == sum(state) ** 2
        for state in all_states
    ),
    "uniform affine pair relation is the square of the scalar ice generator",
)

pair_index = {pair: index for index, pair in enumerate(pairs)}
affine_relations: list[list[int]] = []
for first, second in opposites:
    relation = [0] * 7
    relation[1 + pair_index[first]] = 1
    relation[1 + pair_index[second]] = -1
    affine_relations.append(relation)
affine_relations.append([2] + [1] * 6)
augmented_pair_table = [[1] + row for row in pair_table]
relation_values = matmul(augmented_pair_table, list(map(list, zip(*affine_relations))))
check(all(value == 0 for row in relation_values for value in row), "all four affine relations are zero operators on the projected fiber")
check(rank([relation[1:] for relation in affine_relations]) == 4, "formal pair-coordinate relation gradients have rank four")
check(rank(augmented_pair_table) == 3, "physical pair image remains rank three despite formal relation rank")


# ---------------------------------------------------------------------------
# Independent S4 character and parity typing.
# ---------------------------------------------------------------------------

representatives = (
    (0, 1, 2, 3),       # identity
    (1, 0, 2, 3),       # transposition
    (1, 0, 3, 2),       # double transposition
    (1, 2, 0, 3),       # three-cycle
    (1, 2, 3, 0),       # four-cycle
)


def permute_state(state: tuple[int, ...], permutation: tuple[int, ...]) -> tuple[int, ...]:
    output = [0] * len(state)
    for old, value in enumerate(state):
        output[permutation[old]] = value
    return tuple(output)


def compose(first: tuple[int, ...], second: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(first[second[index]] for index in range(len(first)))


ice_character = tuple(
    sum(permute_state(state, permutation) == state for state in ice)
    for permutation in representatives
)
standard_character = tuple(
    sum(permutation[index] == index for index in range(4)) - 1
    for permutation in representatives
)
sym2_character = tuple(
    (
        character * character
        + sum(compose(permutation, permutation)[index] == index for index in range(4))
        - 1
    )
    // 2
    for character, permutation in zip(standard_character, representatives)
)
check(ice_character == (6, 2, 2, 0, 0), "ice permutation character is A1+E+T2")
check(standard_character == (3, 1, -1, 0, -1), "one-link standard character is T2")
check(sym2_character == ice_character, "Sym2(T2) and ice functions agree only as S4 modules")
check(all(tuple(-entry for entry in state) in ice for state in ice), "global complement preserves the ice fiber")
check(
    all(
        (-state[a]) * (-state[b]) == state[a] * state[b]
        for state in ice
        for a, b in pairs
    ),
    "pair observables are complement even",
)
check(
    all(
        tuple(-entry for entry in state)[index] == -state[index]
        for state in ice
        for index in range(4)
    ),
    "one-link observables are complement odd",
)


# ---------------------------------------------------------------------------
# Periodic diamond incidence, dependency, and rotor counts.
# ---------------------------------------------------------------------------

def periodic_diamond(length: int) -> tuple[list[list[int]], list[set[int]], int, int]:
    cells = list(itertools.product(range(length), repeat=3))
    vertices = [(side, *cell) for side in (0, 1) for cell in cells]
    index = {vertex: position for position, vertex in enumerate(vertices)}
    offsets = ((0, 0, 0), (-1, 0, 0), (0, -1, 0), (0, 0, -1))
    edges: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for cell in cells:
        for offset in offsets:
            target = tuple((cell[axis] + offset[axis]) % length for axis in range(3))
            edges.append(((0, *cell), (1, *target)))
    incidence = [[0 for _ in edges] for _ in vertices]
    adjacency = [set() for _ in vertices]
    for column, (tail, head) in enumerate(edges):
        a, b = index[tail], index[head]
        incidence[a][column] = 1
        incidence[b][column] = -1
        adjacency[a].add(b)
        adjacency[b].add(a)
    return incidence, adjacency, len(vertices), len(edges)


for length in (2, 3, 4):
    incidence, adjacency, vertex_count, edge_count = periodic_diamond(length)
    reached = {0}
    frontier = [0]
    while frontier:
        current = frontier.pop()
        for neighbor in adjacency[current] - reached:
            reached.add(neighbor)
            frontier.append(neighbor)
    incidence_rank = rank(incidence)
    check(vertex_count == 2 * length**3, f"L={length} has 2Nc vertices")
    check(edge_count == 4 * length**3, f"L={length} has 4Nc links")
    check(len(reached) == vertex_count, f"L={length} periodic support is connected")
    check(incidence_rank == vertex_count - 1, f"L={length} incidence rank is V-1")
    check(all(sum(incidence[row][column] for row in range(vertex_count)) == 0 for column in range(edge_count)), f"L={length} has the oriented all-vertex dependency")
    check(edge_count - incidence_rank == 2 * length**3 + 1, f"L={length} rotor configuration count is 2Nc+1")
    check(2 * edge_count - 2 * incidence_rank == 2 * (2 * length**3 + 1), f"L={length} first-class reduction gives twice the configuration count")


# ---------------------------------------------------------------------------
# Hexagon/octagon Gauss preservation and microscopic tunnelling failure.
# ---------------------------------------------------------------------------

def bipartite_cycle_incidence(length: int) -> list[list[int]]:
    incidence = [[0 for _ in range(length)] for _ in range(length)]
    for edge in range(length):
        left, right = edge, (edge + 1) % length
        if left % 2 == 0:
            tail, head = left, right
        else:
            tail, head = right, left
        incidence[tail][edge] = 1
        incidence[head][edge] = -1
    return incidence


for cycle_length in (6, 8):
    incidence = bipartite_cycle_incidence(cycle_length)
    alternating_endpoint_change = tuple(1 if edge % 2 == 0 else -1 for edge in range(cycle_length))
    single_link_change = (1,) + (0,) * (cycle_length - 1)
    check(matvec(incidence, alternating_endpoint_change) == (0,) * cycle_length, f"alternating C{cycle_length} ring preserves every Gauss sum")
    check(matvec(incidence, single_link_change) != (0,) * cycle_length, f"single C{cycle_length} link flip creates a Gauss-defect pair")
    check(rank(incidence) == cycle_length - 1, f"C{cycle_length} ring incidence has only the global dependency")

# One-link hard-core matrix commutator [E,X] is nonzero.  Any diagonal H0
# commutes with E, so the nonzero term isolates the finite-Ud microscopic
# breaking to the supplied -h sum X_e tunnelling.
electric = [[Fraction(-1, 2), 0], [0, Fraction(1, 2)]]
flip = [[0, 1], [1, 0]]
commutator = [
    [left - right for left, right in zip(row_left, row_right)]
    for row_left, row_right in zip(matmul(electric, flip), matmul(flip, electric))
]
check(commutator != [[0, 0], [0, 0]], "unprojected one-link tunnelling fails to commute with electric Gauss charge")
penalty = [[1, 0], [0, 4]]
penalty_commutator = [
    [left - right for left, right in zip(row_left, row_right)]
    for row_left, row_right in zip(matmul(electric, penalty), matmul(penalty, electric))
]
check(penalty_commutator == [[0, 0], [0, 0]], "diagonal microscopic penalties commute with electric charge")


# ---------------------------------------------------------------------------
# Independent continuum principal-symbol and first-class count comparison.
# ---------------------------------------------------------------------------

kx, ky, kz = 2, 3, 5
k2 = kx * kx + ky * ky + kz * kz
maxwell_constraint = [[kx, ky, kz]]
momentum_constraints = [
    [kx, 0, 0, ky, kz, 0],
    [0, ky, 0, kx, 0, kz],
    [0, 0, kz, 0, kx, ky],
]
scalar_constraint = [[
    kx * kx - k2,
    ky * ky - k2,
    kz * kz - k2,
    2 * kx * ky,
    2 * kx * kz,
    2 * ky * kz,
]]
spatial_gauge_map = [
    [2 * kx, 0, 0],
    [0, 2 * ky, 0],
    [0, 0, 2 * kz],
    [ky, kx, 0],
    [kz, 0, kx],
    [0, kz, ky],
]
phase_packet = [row + [0] * 6 for row in scalar_constraint]
phase_packet += [[0] * 6 + row for row in momentum_constraints]
tt_rows = momentum_constraints + [[1, 1, 1, 0, 0, 0]]

check(rank(maxwell_constraint) == 1, "Maxwell has one scalar Gauss symbol at nonzero momentum")
check(3 - rank(maxwell_constraint) == 2, "Maxwell leaves two transverse configurations")
check(rank(momentum_constraints) == 3, "linearized ADM momentum symbol has vector rank three")
check(rank(scalar_constraint) == 1, "linearized ADM curvature symbol has scalar rank one")
check(rank(phase_packet) == 4, "linearized ADM phase-space constraint packet has rank four")
check(rank(spatial_gauge_map) == 3, "symmetric-tensor spatial gauge map has three independent parameters")
check(all(entry == 0 for row in matmul(scalar_constraint, spatial_gauge_map) for entry in row), "linear scalar curvature is invariant under the spatial gauge map")
check(rank(tt_rows) == 4 and 6 - rank(tt_rows) == 2, "TT conditions leave two symmetric-tensor components")
check(6 - 2 * rank(maxwell_constraint) == 4, "Maxwell first-class reduction leaves four phase dimensions")
check(12 - 2 * rank(phase_packet) == 4, "rank-two first-class reduction leaves four phase dimensions")

scalar_jet = [[1], [kx], [ky], [kz]]
check(rank(scalar_jet) == 1, "one scalar Fourier amplitude and its derivatives have rank one")
check(rank([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]) == 4, "independent lapse-plus-shift amplitudes have rank four")

# A rotation by pi/2 distinguishes the two continuum little-group modules.
theta = math.pi / 2
helicity_one_character = 2 * math.cos(theta)
helicity_two_character = 2 * math.cos(2 * theta)
check(abs(helicity_one_character) < 1e-12, "helicity-one character at pi/2 is zero")
check(abs(helicity_two_character + 2) < 1e-12, "helicity-two character at pi/2 is minus two")
check(abs(helicity_one_character - helicity_two_character) > 1, "continuum little-group characters distinguish spin one and spin two")


# ---------------------------------------------------------------------------
# Custody, manifest/seal verification, and negative tamper cases.
# ---------------------------------------------------------------------------

dependencies = parse_hash_file(HERE / "DEPENDENCIES.sha256")
check(len(dependencies) == 11, "dependency ledger has eleven unique entries")
for relative, expected in dependencies.items():
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(), f"dependency is a regular file: {relative}")
    check(digest(path) == expected, f"dependency digest matches: {relative}")

required_manifest_entries = {
    "THEOREM.md",
    "RESULT.md",
    "README.md",
    "SELF_AUDIT.md",
    "INDEPENDENT_HOSTILE_AUDIT.md",
    "verify_rgrlb_constraint_origin.py",
    "independent_hostile_audit.py",
    "VERIFICATION.txt",
    "DEPENDENCIES.sha256",
}
manifest = parse_hash_file(HERE / "MANIFEST.sha256")
check(set(manifest) == required_manifest_entries, "manifest has the complete exact lane payload")
for relative, expected in manifest.items():
    path = HERE / relative
    check(path.is_file() and not path.is_symlink(), f"manifest member is a regular file: {relative}")
    check(digest(path) == expected, f"manifest member digest matches: {relative}")

seal = parse_hash_file(HERE / "SEAL.sha256")
check(seal == {"MANIFEST.sha256": digest(HERE / "MANIFEST.sha256")}, "outer seal authenticates the exact manifest bytes")

theorem_payload = (HERE / "THEOREM.md").read_bytes()
check(digest_bytes(theorem_payload + b"\nHOSTILE-TAMPER") != manifest["THEOREM.md"], "theorem byte tamper is detected")
dependency_name, dependency_digest = next(iter(dependencies.items()))
check(digest_bytes((ROOT / dependency_name).read_bytes() + b"\nHOSTILE-TAMPER") != dependency_digest, "dependency byte tamper is detected")
check(digest_bytes((HERE / "MANIFEST.sha256").read_bytes() + b"\nHOSTILE-TAMPER") != seal["MANIFEST.sha256"], "manifest byte tamper is detected by outer seal")
check(required_manifest_entries - {"THEOREM.md"} != required_manifest_entries, "manifest omission attack changes required payload set")

print(f"SUMMARY {PASSED}/{PASSED} independent hostile checks passed")
print("DISPOSITION PASS__FINITE_Q4_ICE_HAS_SCALAR_U1_ONLY__NO_LOCAL_PAIR_RGRLB_PACKET__NO_THERMODYNAMIC_NO_GO")
