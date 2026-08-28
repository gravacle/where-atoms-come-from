#!/usr/bin/env python3
"""Exact checks for the inherited F3/q4 even-channel kernel boundary."""

from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import combinations, permutations, product
from math import prod
from pathlib import Path


checks = 0


def check(condition, message):
    global checks
    if not condition:
        raise AssertionError(message)
    checks += 1


# ---------------------------------------------------------------------------
# Exact alternating-octagon coefficient at symmetric detuning.

cycle_length = 8
cycle_edges = tuple((vertex, (vertex + 1) % cycle_length)
                    for vertex in range(cycle_length))
cycle_initial = tuple(1 if edge % 2 == 0 else 0
                      for edge in range(cycle_length))


def cycle_gap_coefficient(toggled):
    edge_delta = [0] * cycle_length
    for edge in toggled:
        edge_delta[edge] = 1 if cycle_initial[edge] == 0 else -1
    vertex_delta = [0] * cycle_length
    for edge, (first, second) in enumerate(cycle_edges):
        vertex_delta[first] += edge_delta[edge]
        vertex_delta[second] += edge_delta[edge]
    return sum(value * value for value in vertex_delta)


octagon_classes = Counter()
for order in permutations(range(cycle_length)):
    toggled = set()
    denominators = []
    for edge in order[:-1]:
        toggled.add(edge)
        denominators.append(cycle_gap_coefficient(toggled))
    octagon_classes[tuple(sorted(denominators))] += 1

expected_octagon_classes = Counter({
    (2, 2, 2, 2, 2, 2, 2): 512,
    (2, 2, 2, 2, 2, 2, 4): 1280,
    (2, 2, 2, 2, 2, 4, 4): 2816,
    (2, 2, 2, 2, 4, 4, 4): 4672,
    (2, 2, 2, 2, 4, 4, 6): 1152,
    (2, 2, 2, 4, 4, 4, 4): 5632,
    (2, 2, 2, 4, 4, 4, 6): 3456,
    (2, 2, 2, 4, 4, 6, 6): 2304,
    (2, 2, 4, 4, 4, 4, 4): 4096,
    (2, 2, 4, 4, 4, 4, 6): 4608,
    (2, 2, 4, 4, 4, 6, 6): 5184,
    (2, 2, 4, 4, 6, 6, 6): 3456,
    (2, 2, 4, 4, 6, 6, 8): 1152,
})
check(sum(octagon_classes.values()) == 40320,
      "all 8! alternating-octagon paths")
check(octagon_classes == expected_octagon_classes,
      "thirteen exact octagon denominator classes")
octagon_coefficient = sum(
    Fraction(multiplicity, prod(denominators))
    for denominators, multiplicity in octagon_classes.items()
)
check(octagon_coefficient == Fraction(429, 16),
      "J8=429/16 in h^8/Ud^7 units")


@lru_cache(None)
def octagon_subset_sum(mask):
    if mask == (1 << cycle_length) - 1:
        return Fraction(1)
    total = Fraction(0)
    for edge in range(cycle_length):
        if (mask >> edge) & 1:
            continue
        successor = mask | (1 << edge)
        factor = Fraction(1)
        if successor != (1 << cycle_length) - 1:
            toggled = {position for position in range(cycle_length)
                       if (successor >> position) & 1}
            factor /= cycle_gap_coefficient(toggled)
        total += factor * octagon_subset_sum(successor)
    return total


check(octagon_subset_sum(0) == octagon_coefficient,
      "independent subset recursion reproduces octagon coefficient")


def cycle_component_count(mask):
    """Connected components of a nonempty proper edge subset of C8."""
    return sum(
        1
        for edge in range(cycle_length)
        if ((mask >> edge) & 1) and
        not ((mask >> ((edge - 1) % cycle_length)) & 1)
    )


proper_masks = range(1, (1 << cycle_length) - 1)
check(all(cycle_gap_coefficient({edge for edge in range(cycle_length)
                                 if (mask >> edge) & 1}) ==
          2 * cycle_component_count(mask)
          for mask in proper_masks),
      "independent boundary-component formula gives every proper-subset gap")
check(cycle_gap_coefficient(set(range(cycle_length))) == 0,
      "full alternating octagon returns exactly to the ice manifold")
check((-1) ** cycle_length * (-1) ** (cycle_length - 1) == -1,
      "eight flips and seven negative resolvents give -J8")


# ---------------------------------------------------------------------------
# Exact order-eight parity/operator classification.


def integer_partitions(total, maximum=None):
    if total == 0:
        yield ()
        return
    maximum = total if maximum is None else min(maximum, total)
    for first in range(maximum, 0, -1):
        for tail in integer_partitions(total - first, first):
            yield (first,) + tail


partitions_eight = tuple(integer_partitions(8))
diagonal_partitions = {
    partition for partition in partitions_eight
    if all(value % 2 == 0 for value in partition)
}
check(diagonal_partitions == {
    (8,), (6, 2), (4, 4), (4, 2, 2), (2, 2, 2, 2)
}, "complete diagonal multiplicity partitions at order eight")

odd_count_six = {
    partition for partition in partitions_eight
    if sum(value % 2 for value in partition) == 6
}
odd_count_eight = {
    partition for partition in partitions_eight
    if sum(value % 2 for value in partition) == 8
}
check(odd_count_six == {(3, 1, 1, 1, 1, 1),
                        (2, 1, 1, 1, 1, 1, 1)},
      "hexagon dressing is one repeated edge")
check(odd_count_eight == {(1, 1, 1, 1, 1, 1, 1, 1)},
      "new order-eight transition is one eight-edge cycle")


# ---------------------------------------------------------------------------
# A literal coordination-four, girth-six hostile support: the incidence graph
# of PG(2,3).  Its two degree-two ice states differ by one alternating hexagon.


def normalized_projective_vectors(field_order):
    representatives = []
    seen = set()
    for vector in product(range(field_order), repeat=3):
        if vector == (0, 0, 0):
            continue
        first = next(value for value in vector if value)
        inverse = pow(first, -1, field_order)
        normalized = tuple((value * inverse) % field_order for value in vector)
        if normalized not in seen:
            seen.add(normalized)
            representatives.append(normalized)
    return tuple(sorted(representatives))


projective_points = normalized_projective_vectors(3)
point_count = len(projective_points)
host_edges = tuple(
    (point_index, point_count + line_index)
    for point_index, point in enumerate(projective_points)
    for line_index, normal in enumerate(projective_points)
    if sum(x * y for x, y in zip(point, normal)) % 3 == 0
)
host_vertices = tuple(range(2 * point_count))
host_incidence = {vertex: [] for vertex in host_vertices}
for edge_index, (first, second) in enumerate(host_edges):
    host_incidence[first].append(edge_index)
    host_incidence[second].append(edge_index)

check(len(projective_points) == 13, "PG(2,3) has thirteen points and lines")
check(len(host_edges) == 52, "PG(2,3) incidence graph has fifty-two edges")
check(all(len(host_incidence[vertex]) == 4 for vertex in host_vertices),
      "host is coordination four")


def graph_girth(edges):
    vertices = sorted({vertex for edge in edges for vertex in edge})
    neighbors = {vertex: [] for vertex in vertices}
    for first, second in edges:
        neighbors[first].append(second)
        neighbors[second].append(first)
    best = len(vertices) + 1
    for root in vertices:
        distance = {root: 0}
        parent = {root: None}
        queue = [root]
        for vertex in queue:
            for neighbor in neighbors[vertex]:
                if neighbor not in distance:
                    distance[neighbor] = distance[vertex] + 1
                    parent[neighbor] = vertex
                    queue.append(neighbor)
                elif parent[vertex] != neighbor:
                    best = min(best,
                               distance[vertex] + distance[neighbor] + 1)
    return best


check(graph_girth(host_edges) == 6, "host girth is six")

ice_state = frozenset({
    0, 1, 4, 6, 8, 10, 12, 13, 16, 18, 20, 21, 25,
    27, 29, 31, 32, 34, 38, 39, 41, 42, 45, 47, 50, 51,
})
flippable_hexagon = frozenset({0, 3, 16, 17, 28, 29})
switched_state = ice_state ^ flippable_hexagon


def degrees(state):
    result = Counter()
    for edge_index in state:
        for vertex in host_edges[edge_index]:
            result[vertex] += 1
    return result


check(all(degrees(ice_state)[vertex] == 2 for vertex in host_vertices),
      "source state obeys degree-two ice")
check(all(degrees(switched_state)[vertex] == 2 for vertex in host_vertices),
      "hexagon-switched state obeys degree-two ice")
hexagon_vertices = Counter(
    vertex for edge_index in flippable_hexagon
    for vertex in host_edges[edge_index]
)
check(len(hexagon_vertices) == 6 and
      set(hexagon_vertices.values()) == {2},
      "declared six edges form a simple cycle")
check(len(flippable_hexagon & ice_state) == 3,
      "declared hexagon has three occupied edges")
for vertex in hexagon_vertices:
    incident_cycle_edges = [edge for edge in host_incidence[vertex]
                            if edge in flippable_hexagon]
    check(len(incident_cycle_edges) == 2 and
          sum(edge in ice_state for edge in incident_cycle_edges) == 1,
          "declared hexagon is alternating at every vertex")


# Canonical colored incidence keys are precomputed for speed.  A key consists
# of k occupation bits followed by one bit for each pairwise edge incidence.


@lru_cache(None)
def index_pairs(size):
    return tuple((first, second)
                 for first in range(size)
                 for second in range(first + 1, size))


@lru_cache(None)
def canonical_lookup(size):
    pair_list = index_pairs(size)
    permutation_list = tuple(permutations(range(size)))
    lookup = {}
    bit_count = size + len(pair_list)
    for bits in product((0, 1), repeat=bit_count):
        occupations = bits[:size]
        adjacency = [[0] * size for _ in range(size)]
        for offset, (first, second) in enumerate(pair_list):
            adjacency[first][second] = bits[size + offset]
            adjacency[second][first] = bits[size + offset]
        images = []
        for permutation in permutation_list:
            image = tuple(occupations[permutation[position]]
                          for position in range(size))
            image += tuple(adjacency[permutation[first]][permutation[second]]
                           for first, second in pair_list)
            images.append(image)
        lookup[bits] = min(images)
    return lookup


edge_endpoints = tuple(frozenset(edge) for edge in host_edges)


def cluster_key(selected, state):
    size = len(selected)
    bits = tuple(int(edge in state) for edge in selected)
    bits += tuple(int(bool(edge_endpoints[selected[first]] &
                           edge_endpoints[selected[second]]))
                  for first, second in index_pairs(size))
    return canonical_lookup(size)[bits]


def cluster_census(state, size):
    counts = Counter()
    representatives = {}
    for selected in combinations(range(len(host_edges)), size):
        key = cluster_key(selected, state)
        counts[key] += 1
        representatives.setdefault(key, selected)
    return counts, representatives


censuses = {}
representatives = {}
for size in (2, 3, 4):
    source_counts, source_representatives = cluster_census(ice_state, size)
    target_counts, _ = cluster_census(switched_state, size)
    check(source_counts == target_counts,
          f"colored {size}-edge census is invariant under hostile hex flip")
    censuses[size] = source_counts
    representatives[size] = source_representatives

check(len(censuses[2]) == 6, "six colored two-edge types")
check(len(censuses[3]) == 18, "eighteen colored three-edge types")
check(len(censuses[4]) == 55, "fifty-five colored four-edge types")


def virtual_gap(selected, state, parity_mask):
    degree_delta = Counter()
    for local_index, edge_index in enumerate(selected):
        if not ((parity_mask >> local_index) & 1):
            continue
        delta = -1 if edge_index in state else 1
        for vertex in host_edges[edge_index]:
            degree_delta[vertex] += delta
    return sum(value * value for value in degree_delta.values())


def irreducible_word_weight(selected, state, multiplicities):
    size = len(selected)

    @lru_cache(None)
    def recurse(counts, parity_mask):
        used = sum(counts)
        if used == 8:
            return Fraction(1)
        total = Fraction(0)
        for local_index in range(size):
            if counts[local_index] >= multiplicities[local_index]:
                continue
            successor_counts = list(counts)
            successor_counts[local_index] += 1
            successor_mask = parity_mask ^ (1 << local_index)
            if used + 1 < 8:
                # With at most four distinct edges on a girth-six host, the
                # only possible return to the ice manifold is empty parity.
                if successor_mask == 0:
                    continue
                gap = virtual_gap(selected, state, successor_mask)
                total += (Fraction(1, gap) *
                          recurse(tuple(successor_counts), successor_mask))
            else:
                total += recurse(tuple(successor_counts), successor_mask)
        return total

    return recurse((0,) * size, 0)


def diagonal_partition_weight(selected, state):
    size = len(selected)
    if size == 2:
        multiplicity_assignments = {(6, 2), (2, 6), (4, 4)}
    elif size == 3:
        multiplicity_assignments = set(permutations((4, 2, 2)))
    elif size == 4:
        multiplicity_assignments = {(2, 2, 2, 2)}
    else:
        raise AssertionError("unsupported diagonal cluster size")
    return sum((irreducible_word_weight(selected, state, assignment)
                for assignment in multiplicity_assignments), Fraction(0))


def irreducible_diagonal_eight(state):
    total = Fraction(0)
    for size in (2, 3, 4):
        counts, reps = cluster_census(state, size)
        for key, multiplicity in counts.items():
            total += (multiplicity *
                      diagonal_partition_weight(reps[key], state))
    return total


source_diagonal_eight = irreducible_diagonal_eight(ice_state)
target_diagonal_eight = irreducible_diagonal_eight(switched_state)
check(source_diagonal_eight == target_diagonal_eight,
      "exact hostile irreducible order-eight diagonal sum is scalar")
check(source_diagonal_eight == Fraction(2526594309109, 13608000),
      "hostile exact rational diagonal checksum")


# ---------------------------------------------------------------------------
# Literal alternating octagon on periodic diamond support.


linear_size = 4
cells = tuple(product(range(linear_size), repeat=3))
cell_index = {cell: index for index, cell in enumerate(cells)}
cell_count = len(cells)
shifts = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
diamond_edges = []
direction_edges = [[] for _ in shifts]
for cell in cells:
    for direction, shift in enumerate(shifts):
        target = tuple((cell[axis] + shift[axis]) % linear_size
                       for axis in range(3))
        direction_edges[direction].append(len(diamond_edges))
        diamond_edges.append((cell_index[cell],
                              cell_count + cell_index[target]))

diamond_ice = frozenset(direction_edges[0] + direction_edges[1])
diamond_octagon = (0, 15, 12, 11, 8, 7, 4, 3)
diamond_vertex_counts = Counter(
    vertex for edge_index in diamond_octagon
    for vertex in diamond_edges[edge_index]
)
check(len(diamond_vertex_counts) == 8 and
      set(diamond_vertex_counts.values()) == {2},
      "periodic diamond witness contains a simple octagon")
octagon_occupations = tuple(int(edge in diamond_ice)
                            for edge in diamond_octagon)
check(all(octagon_occupations[index] !=
          octagon_occupations[(index + 1) % 8]
          for index in range(8)),
      "periodic diamond octagon is alternating")
diamond_switched = diamond_ice ^ frozenset(diamond_octagon)
diamond_degrees = Counter()
for edge_index in diamond_switched:
    for vertex in diamond_edges[edge_index]:
        diamond_degrees[vertex] += 1
check(all(diamond_degrees[vertex] == 2
          for vertex in range(2 * cell_count)),
      "octagon flip preserves periodic diamond ice")


# ---------------------------------------------------------------------------
# Pole-boundary algebra: a finite first-order kernel correction cannot make
# the geometric Bethe-Salpeter denominator that only appears after resummation.


for bubble, kernel in ((Fraction(2), Fraction(1, 10)),
                       (Fraction(5, 3), Fraction(-1, 7)),
                       (Fraction(11, 4), Fraction(3, 20))):
    first_order = bubble + bubble * kernel * bubble
    resummed = bubble / (1 - kernel * bubble)
    check(first_order == bubble * (1 + kernel * bubble),
          "first-order connected correction is polynomial in kernel")
    check(resummed - first_order ==
          bubble * (kernel * bubble) ** 2 / (1 - kernel * bubble),
          "new Bethe-Salpeter denominator requires higher ladders")


root = Path(__file__).resolve().parent
theorem = (root / "THEOREM.md").read_text()
result = (root / "RESULT.md").read_text()
self_audit = (root / "SELF_AUDIT.md").read_text()
independent_audit = (root / "INDEPENDENT_AUDIT.md").read_text()
joined = theorem + result + self_audit + independent_audit
theorem_flat = " ".join(theorem.split())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


workspace = root.parent
dependency_hashes = {
    "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md":
        "00eba581b90fb9f0b25e3fad1362b055049824897433ba06e356ab9b1f6c76ec",
    "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md":
        "5e68e4a8c62ad89cff309781a5cb54e071092e82594e60c04e7992414bc18dbe",
    "LANE_CROSS_ALPHA_GRA_F3_DIAMOND_SIXTH_ORDER_V001/THEOREM.md":
        "211b1aa61917c98dccae278129a8016a1a14f73587908bfeceeba090a808536c",
    "LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/THEOREM.md":
        "05f4a619a6f80aa40c48570ab4035ab874426502a31468a08f435e66610bd769",
    "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/THEOREM.md":
        "cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98",
    "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/INDEPENDENT_AUDIT.md":
        "c52eab9d701d1c6e82f1d7ec395841f4d2810e96cccbc3e2504760b6742e81e4",
    "LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/THEOREM.md":
        "98e2b3bc7a1c998d7839dc1a6b435cc1c8ed6d5a622ba45f63571be9ef646452",
}
for relative_path, expected_hash in dependency_hashes.items():
    path = workspace / relative_path
    check(path.is_file() and sha256(path) == expected_hash,
          f"dependency custody pinned: {relative_path}")

for required in (
    "J_8={429h^8\\over16U_d^7}",
    "V_8=0",
    "dressed length-six transitions",
    "new alternating length-eight ring transitions",
    "not the first inherited non-Gaussian interaction",
    "g_{TT}^{\\rm match}",
    "strictly at the single-insertion expression",
    "does not establish a new isolated tensor pole",
    "fully dressed massless 1PI four-point function",
    "two-particle irreducible",
    "G^{(4),c}_{TT}",
    "full Bethe--Salpeter resummation",
    "No new microscopic interaction",
):
    check(required in theorem_flat, f"theorem retains load-bearing claim: {required}")

for forbidden in (
    "order eight proves a graviton",
    "the octagon is the graviton",
    "F3 derives the diamond lattice",
    "g_TT is known numerically",
    "Gaussian Maxwell is the exact microscopic Hamiltonian",
    "gravity is closed by J8",
    "finite-order Hamiltonian can never create a pole",
    "fully dressed 1PI kernel is analytic",
):
    check(forbidden not in joined,
          f"forbidden promotion absent: {forbidden}")

for required in (
    "PASS_WITH_SCOPED_CONTINUUM_CORRECTION",
    "429/16",
    "V_8=0",
    "8!=40320",
    "strict single-insertion",
    "not a nonperturbative no-pole theorem",
    "connected TT four-point",
    "two-particle-irreducible",
    "Ward",
    "residue",
):
    check(required in independent_audit,
          f"independent audit retains required finding: {required}")


print(f"Inherited TT-kernel boundary verification: {checks} passed, 0 failed")
print(f"J8 symmetric coefficient: {octagon_coefficient} h^8/Ud^7")
print(f"PG(2,3) D8 checksum: {source_diagonal_eight}")
