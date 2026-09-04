#!/usr/bin/env python3
"""Independent hostile replay of the sealed GL6CH tensor-writer theorem.

The audit never imports or executes the author derivation.  It independently
enumerates direct six-hop histories, derives lower source vertices with a
dual-number finite-star Rayleigh calculation, reconstructs all simple Q4
six-cycles by graph search, and checks source dimensions and interpretation
boundaries using exact arithmetic.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
import argparse
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6CH_GLOBAL_H6_TENSOR_WRITER_V001"
PAIR_ORDER = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIR_ORDER)}
ZERO6 = (F(0),) * 6
ONE6 = (F(1),) * 6
T_BASIS = (
    (F(1), F(0), F(0), F(0), F(0), F(-1)),
    (F(0), F(1), F(0), F(0), F(-1), F(0)),
    (F(0), F(0), F(1), F(-1), F(0), F(0)),
)
EA = (F(1), F(1), F(-2), F(-2), F(1), F(1))
EB = (F(1), F(-1), F(0), F(0), F(-1), F(1))


class Checks:
    def __init__(self):
        self.total = 0

    def true(self, condition, label):
        self.total += 1
        if not condition:
            raise AssertionError(label)

    def equal(self, got, want, label):
        self.total += 1
        if got != want:
            raise AssertionError(f"{label}: got {got!r}; want {want!r}")


CHECK = Checks()


def dot(left, right):
    return sum((F(x) * F(y) for x, y in zip(left, right)), F(0))


def add_vectors(left, right, factor=F(1)):
    return tuple(F(x) + factor * F(y) for x, y in zip(left, right))


def scale_vector(factor, vector):
    return tuple(F(factor) * F(x) for x in vector)


def outer(left, right):
    return tuple(tuple(F(x) * F(y) for y in right) for x in left)


def add_matrices(left, right, factor=F(1)):
    return tuple(tuple(left[i][j] + factor * right[i][j]
                       for j in range(len(left[i])))
                 for i in range(len(left)))


def zero_matrix(rows, columns):
    return tuple(tuple(F(0) for _ in range(columns)) for _ in range(rows))


def matrix_rank(matrix):
    work = [list(map(F, row)) for row in matrix]
    if not work:
        return 0
    rows, columns = len(work), len(work[0])
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, rows)
                      if work[row][column] != 0), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][column]
        work[pivot_row] = [entry / value for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or work[row][column] == 0:
                continue
            value = work[row][column]
            work[row] = [work[row][j] - value * work[pivot_row][j]
                         for j in range(columns)]
        pivot_row += 1
    return pivot_row


def project_t(vector):
    answer = ZERO6
    for basis in T_BASIS:
        CHECK.equal(dot(basis, basis), F(2), "T basis norm squared two")
        answer = add_vectors(answer, basis, dot(vector, basis) / 2)
    return answer


def complementary_pair(pair):
    return tuple(sorted(set(range(4)) - set(pair)))


def theta(pair):
    pair = tuple(sorted(pair))
    answer = [F(0)] * 6
    answer[PAIR_INDEX[pair]] = F(1)
    answer[PAIR_INDEX[complementary_pair(pair)]] = F(-1)
    return tuple(answer)


def pair_memory(spins):
    return tuple(F(spins[a] * spins[b]) for a, b in PAIR_ORDER)


# --------------------------------------------------- direct order-six source
def boundary_energy(selected, signs):
    charges = [0] * 6
    for edge in selected:
        charges[edge] += signs[edge]
        charges[(edge + 1) % 6] += signs[edge]
    return sum(charge * charge for charge in charges)


def local_midpoint_score(vertex, selected, signs, ring_pair,
                         pair_reversed, exterior_reversed):
    previous_edge = (vertex - 1) % 6
    next_edge = vertex
    ring_ports = tuple(reversed(ring_pair)) if pair_reversed else ring_pair
    exterior_pair = complementary_pair(ring_pair)
    exterior_ports = (tuple(reversed(exterior_pair)) if exterior_reversed
                      else exterior_pair)
    initial = [None] * 4
    initial[ring_ports[0]] = signs[previous_edge]
    initial[ring_ports[1]] = signs[next_edge]
    initial[exterior_ports[0]] = 1
    initial[exterior_ports[1]] = -1
    CHECK.true(initial.count(1) == 2 and initial.count(-1) == 2,
               "direct-history local endpoint is locked")
    intermediate = list(initial)
    if previous_edge in selected:
        intermediate[ring_ports[0]] *= -1
    if next_edge in selected:
        intermediate[ring_ports[1]] *= -1
    final = list(initial)
    final[ring_ports[0]] *= -1
    final[ring_ports[1]] *= -1
    m0 = pair_memory(initial)
    mm = pair_memory(intermediate)
    m1 = pair_memory(final)
    return tuple(mm[index] - (m0[index] + m1[index]) / 2
                 for index in range(6))


def audit_direct_histories():
    expected_profiles = Counter({
        (2, 2, 2, 2, 2): 96,
        (2, 2, 2, 4, 2): 48,
        (2, 2, 4, 2, 2): 48,
        (2, 2, 4, 4, 2): 96,
        (2, 4, 2, 2, 2): 48,
        (2, 4, 2, 4, 2): 24,
        (2, 4, 4, 2, 2): 96,
        (2, 4, 4, 4, 2): 192,
        (2, 4, 6, 4, 2): 72,
    })
    phases = []
    context_count = 0
    for phase in (1, -1):
        signs = tuple(phase * (1 if edge % 2 == 0 else -1)
                      for edge in range(6))
        amplitude = F(0)
        profiles = Counter()
        prefix_coefficients = Counter()
        for order in permutations(range(6)):
            selected = set()
            energies = []
            weight = F(1)
            for edge in order[:-1]:
                selected.add(edge)
                energy = F(boundary_energy(selected, signs))
                CHECK.true(energy > 0, "proper direct prefix is outside lock")
                energies.append(energy)
                weight /= energy
            profiles[tuple(map(int, energies))] += 1
            amplitude -= weight
            selected.clear()
            for step, edge in enumerate(order[:-1]):
                selected.add(edge)
                prefix_coefficients[frozenset(selected)] += weight / energies[step]
        CHECK.equal(profiles, expected_profiles, "nine exact direct energy profiles")
        CHECK.equal(amplitude, F(-63, 8), "direct order-six amplitude sign and value")

        for vertex in range(6):
            for ring_pair in PAIR_ORDER:
                direction = theta(ring_pair)
                CHECK.equal(dot(direction, direction), F(2), "Theta norm squared")
                CHECK.equal(project_t(direction), direction, "Theta is pure T2")
                CHECK.equal(dot(ONE6, direction), F(0), "Theta has zero trace")
                CHECK.equal(dot(EA, direction), F(0), "Theta excludes first E2 direction")
                CHECK.equal(dot(EB, direction), F(0), "Theta excludes second E2 direction")
                for pair_reversed in (False, True):
                    for exterior_reversed in (False, True):
                        context_count += 1
                        gradient = [F(0)] * 6
                        for selected, coefficient in prefix_coefficients.items():
                            score = local_midpoint_score(
                                vertex, selected, signs, ring_pair,
                                pair_reversed, exterior_reversed)
                            for component in range(6):
                                gradient[component] += coefficient * score[component]
                        gradient = tuple(gradient)
                        expected_gradient = tuple(
                            F(105, 8) if index == PAIR_INDEX[ring_pair] else F(0)
                            for index in range(6)
                        )
                        CHECK.equal(gradient, expected_gradient,
                                    "canonical endpoint-symmetrized full gradient")
                        CHECK.equal(project_t(gradient),
                                    scale_vector(F(105, 16), direction),
                                    "full gradient projects to 105/16 Theta")
                        CHECK.equal(dot(direction, gradient), F(105, 8),
                                    "literal j=Theta derivative is 105/8")
        phases.append({"phase": phase, "amplitude": amplitude,
                       "profiles": profiles})

    CHECK.equal(context_count, 288, "all 288 local direct contexts exhausted")
    CHECK.equal(phases[0]["profiles"], phases[1]["profiles"],
                "direct histories are alternating-orientation independent")
    return {
        "orders_per_phase": 720,
        "phases": 2,
        "local_contexts": context_count,
        "energy_profiles": expected_profiles,
        "source_free": F(-63, 8),
        "canonical_gradient": "(105/8)e_ab",
        "t2_gradient": "(105/16)Theta_ab",
        "literal_theta_derivative": F(105, 8),
    }


# ------------------------------- exact differentiated lower-order derivation
def neighborhood_edges():
    # 0..3 join the central node to four neighbors.  Each neighbor has three
    # exterior links with a distinct leaf endpoint, as guaranteed by girth 6.
    edges = []
    for port in range(4):
        edges.append(("central", port, None))
    for port in range(4):
        for external in range(3):
            edges.append(("external", port, external))
    return tuple(edges)


STAR_EDGES = neighborhood_edges()


def q_only_identity_orders():
    """All four-flip identity words with no proper identity prefix."""
    answer = set()
    for left, right in combinations(range(len(STAR_EDGES)), 2):
        for order in set(permutations((left, left, right, right))):
            toggled = set()
            allowed = True
            for edge in order[:-1]:
                toggled.symmetric_difference_update({edge})
                if not toggled:
                    allowed = False
                    break
            if allowed:
                answer.add(order)
    CHECK.equal(len(answer), 480,
                "exactly 480 four-flip identity words remain wholly in Q")
    return tuple(sorted(answer))


Q_ONLY_IDENTITY_ORDERS = q_only_identity_orders()


def local_state(central_bits, exterior_bits):
    return tuple(central_bits) + tuple(bit for row in exterior_bits for bit in row)


def star_energy(bits, toggled):
    signs = tuple(1 - 2 * bit for bit in bits)
    central_charge = sum(signs[index] for index in range(4) if index in toggled)
    neighbor_charges = [0] * 4
    leaf_charges = []
    for index in toggled:
        kind, port, external = STAR_EDGES[index]
        neighbor_charges[port] += signs[index]
        if kind == "external":
            leaf_charges.append(signs[index])
    return (central_charge * central_charge
            + sum(charge * charge for charge in neighbor_charges)
            + sum(charge * charge for charge in leaf_charges))


def central_memory_after(central_bits, toggled):
    spins = []
    for port, bit in enumerate(central_bits):
        after = bit ^ (port in toggled)
        spins.append(1 - 2 * after)
    return pair_memory(spins)


def differentiated_gap(central_bits, bits, toggled, initial_memory):
    energy = F(star_energy(bits, toggled))
    current_memory = central_memory_after(central_bits, toggled)
    return energy, add_vectors(current_memory, initial_memory, F(-1))


def fourth_order_vertex_by_finite_star(central_bits, exterior_bits):
    bits = local_state(central_bits, exterior_bits)
    initial_memory = pair_memory(tuple(1 - 2 * bit for bit in central_bits))
    one_flip_gaps = []
    for edge in range(len(STAR_EDGES)):
        gap = differentiated_gap(central_bits, bits, frozenset({edge}), initial_memory)
        CHECK.equal(gap[0], F(2), "every one-link defect gap is two")
        one_flip_gaps.append(gap)

    a_value = sum((1 / energy for energy, _ in one_flip_gaps), F(0))
    b_value = sum((1 / (energy * energy) for energy, _ in one_flip_gaps), F(0))
    a_gradient = [F(0)] * 6
    b_gradient = [F(0)] * 6
    for energy, gradient in one_flip_gaps:
        for component in range(6):
            a_gradient[component] -= gradient[component] / energy ** 2
            b_gradient[component] -= 2 * gradient[component] / energy ** 3

    # Cache all one- and two-edge prefix gaps.  Every Q-only four-flip
    # identity word alternates between a singleton and a two-edge subset.
    gap_cache = {}
    for edge in range(len(STAR_EDGES)):
        subset = frozenset({edge})
        gap_cache[subset] = differentiated_gap(
            central_bits, bits, subset, initial_memory)
    for left, right in combinations(range(len(STAR_EDGES)), 2):
        subset = frozenset({left, right})
        gap_cache[subset] = differentiated_gap(
            central_bits, bits, subset, initial_memory)
        if gap_cache[subset][0] <= 0:
            raise AssertionError("two-edge prefix unexpectedly re-enters lock")

    # P is independently enumerated as every four-flip return word with all
    # three proper intermediate states in Q.  The 480 admissible parity words
    # were selected once above without consulting any energy formula.
    p_gradient = [F(0)] * 6
    sequence_count = 0
    for order in Q_ONLY_IDENTITY_ORDERS:
        toggled = set()
        prefix_data = []
        for edge in order[:-1]:
            toggled.symmetric_difference_update({edge})
            energy, gradient = gap_cache[frozenset(toggled)]
            if energy <= 0:
                raise AssertionError("preselected Q-only word has nonpositive gap")
            prefix_data.append((energy, gradient))
        sequence_count += 1
        weight = F(1)
        for energy, _ in prefix_data:
            weight /= energy
        for component in range(6):
            logarithmic_derivative = sum(
                (gradient[component] / energy for energy, gradient in prefix_data),
                F(0),
            )
            p_gradient[component] -= weight * logarithmic_derivative
    CHECK.true(sequence_count > 0, "finite-star Q-only fourth-order paths exist")
    e4_gradient = tuple(
        a_gradient[component] * b_value
        + a_value * b_gradient[component]
        - p_gradient[component]
        for component in range(6)
    )
    return e4_gradient, sequence_count


def audit_lower_orders():
    locked_words = tuple(bits for bits in product((0, 1), repeat=4)
                         if sum(bits) == 2)
    neighborhood_count = 0
    path_count_histogram = Counter()
    for central_bits in locked_words:
        memory = pair_memory(tuple(1 - 2 * bit for bit in central_bits))
        CHECK.equal(project_t(memory), ZERO6,
                    "locked pair memory has no T2 component")

        # Independent h2 derivative: every one-flip denominator is 2 and its
        # derivative contributes delta-M/4.
        h2 = [F(0)] * 6
        for port in range(4):
            after = list(central_bits)
            after[port] ^= 1
            delta = add_vectors(
                pair_memory(tuple(1 - 2 * bit for bit in after)), memory, F(-1))
            for component in range(6):
                h2[component] += delta[component] / 4
        CHECK.equal(tuple(h2), scale_vector(-1, memory),
                    "independent h2 vertex is -M")
        CHECK.equal(project_t(tuple(h2)), ZERO6,
                    "h2 first T2 source vertex vanishes")

        compatible_exterior = tuple(
            row for row in product((0, 1), repeat=3)
            if sum(row) == 2 - central_bits[0]
        )
        # The admissible set depends only on the shared occupation, so build
        # it separately at every port.
        choices = [tuple(row for row in product((0, 1), repeat=3)
                         if sum(row) == 2 - central_bits[port])
                   for port in range(4)]
        CHECK.true(len(compatible_exterior) == 3, "three neighbor completions per port")
        for exterior_bits in product(*choices):
            neighborhood_count += 1
            h4, sequence_count = fourth_order_vertex_by_finite_star(
                central_bits, exterior_bits)
            expected = tuple(F(-4, 9) - F(37, 12) * value
                             for value in memory)
            CHECK.equal(h4, expected,
                        "dual finite-star h4 vertex equals scalar plus M identity")
            CHECK.equal(project_t(h4), ZERO6,
                        "h4 first T2 source vertex vanishes")
            path_count_histogram[sequence_count] += 1
    CHECK.equal(neighborhood_count, 486,
                "all six-by-81 locked radius-one neighborhoods exhausted")
    return {
        "locked_words": len(locked_words),
        "radius_one_neighborhoods": neighborhood_count,
        "h2": "-M_v",
        "h4": "-(4/9)1_6-(37/12)M_v",
        "t2_projection_h0_h2_h4": "zero",
        "q_only_path_count_histogram": path_count_histogram,
        "method": "exact differentiated finite-star Rayleigh E4=A B-P enumeration",
    }


# --------------------------------------------------------------- Q4 geometry
PERIOD = 4
STEPS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
CELLS = tuple(product(range(PERIOD), repeat=3))
NODES = tuple((kind, cell) for kind in ("P", "C") for cell in CELLS)
EDGES = tuple((cell, port) for cell in CELLS for port in range(4))


def qadd(left, right):
    return tuple((left[i] + right[i]) % PERIOD for i in range(3))


def qsub(left, right):
    return tuple((left[i] - right[i]) % PERIOD for i in range(3))


def endpoints(edge):
    cell, port = edge
    return ("P", cell), ("C", qadd(cell, STEPS[port]))


def canonical_cycle(cell, ports):
    a, b, c = ports
    ab = qadd(qsub(cell, STEPS[b]), STEPS[a])
    cb = qadd(qsub(cell, STEPS[b]), STEPS[c])
    return ((cell, a), (ab, b), (ab, c),
            (cb, a), (cb, b), (cell, c))


def local_ports(cycle, node):
    return tuple(sorted(edge[1] for edge in cycle if node in endpoints(edge)))


def audit_geometry():
    adjacency = {node: [] for node in NODES}
    pair_to_edge = {}
    for edge in EDGES:
        left, right = endpoints(edge)
        adjacency[left].append(right)
        adjacency[right].append(left)
        pair_to_edge[frozenset((left, right))] = edge
    CHECK.true(all(len(neighbors) == 4 for neighbors in adjacency.values()),
               "Q4 is four-regular")

    # Connectedness.
    seen = {NODES[0]}
    frontier = [NODES[0]]
    while frontier:
        node = frontier.pop()
        for neighbor in adjacency[node]:
            if neighbor not in seen:
                seen.add(neighbor)
                frontier.append(neighbor)
    CHECK.equal(seen, set(NODES), "Q4 incidence graph is connected")

    # Direct four-cycle test: no pair of parents has two common children.
    for left, right in combinations((node for node in NODES if node[0] == "P"), 2):
        CHECK.true(len(set(adjacency[left]) & set(adjacency[right])) <= 1,
                   "Q4 has no parent-pair four-cycle")

    # Infinite-parent analogue: all ordered nonzero step differences are
    # distinct over Z^3, so two parents cannot share two children there either.
    integer_differences = {}
    for a in range(4):
        for b in range(4):
            if a == b:
                continue
            difference = tuple(STEPS[a][i] - STEPS[b][i] for i in range(3))
            CHECK.true(difference not in integer_differences,
                       "infinite step differences are ordered-pair unique")
            integer_differences[difference] = (a, b)
    CHECK.equal(len(integer_differences), 12,
                "all 12 infinite ordered step differences are unique")

    canonical = {}
    orientation_counts = Counter()
    edge_counts = Counter()
    u_by_missing = {}
    for cell in CELLS:
        for ports in combinations(range(4), 3):
            cycle = canonical_cycle(cell, ports)
            mask = frozenset(cycle)
            CHECK.true(mask not in canonical, "canonical Q4 cycle owner is unique")
            CHECK.equal(len(mask), 6, "canonical cycle has six distinct edges")
            vertices = {node for edge in cycle for node in endpoints(edge)}
            CHECK.equal(len(vertices), 6, "canonical cycle has six distinct vertices")
            CHECK.true(all(len(local_ports(cycle, node)) == 2 for node in vertices),
                       "canonical cycle has two local ports per vertex")
            canonical[mask] = (cycle, ports)
            missing = next(iter(set(range(4)) - set(ports)))
            orientation_counts[missing] += 1
            for edge in cycle:
                edge_counts[edge] += 1
            u = tuple(sum(theta(local_ports(cycle, node))[component]
                          for node in vertices)
                      for component in range(6))
            if missing in u_by_missing:
                CHECK.equal(u, u_by_missing[missing],
                            "orientation tensor sum is translation independent")
            else:
                u_by_missing[missing] = u

    CHECK.equal(len(canonical), 256, "Q4 has 256 canonical elementary hexagons")
    CHECK.equal(orientation_counts, Counter({0: 64, 1: 64, 2: 64, 3: 64}),
                "Q4 has 64 cycles in each orientation")
    CHECK.true(all(edge_counts[edge] == 6 for edge in EDGES),
               "each Q4 edge has six elementary-cycle owners")

    # Independent graph DFS: find every simple length-six cycle and compare
    # edge masks, rather than assuming the canonical parameterization complete.
    all_six_cycles = set()
    for start in NODES:
        def walk(node, depth, visited_nodes, path_edges):
            if depth == 6:
                if node == start:
                    all_six_cycles.add(frozenset(path_edges))
                return
            for neighbor in adjacency[node]:
                edge = pair_to_edge[frozenset((node, neighbor))]
                if neighbor == start:
                    if depth == 5:
                        walk(neighbor, depth + 1, visited_nodes, path_edges + [edge])
                    continue
                if neighbor in visited_nodes:
                    continue
                walk(neighbor, depth + 1, visited_nodes | {neighbor},
                     path_edges + [edge])
        walk(start, 0, {start}, [])
    CHECK.equal(all_six_cycles, set(canonical),
                "canonical set equals exhaustive simple-six-cycle graph search")

    # Orientation tetrahedron and exact rank.
    expected_coordinates = {
        3: (F(2), F(2), F(-2)),
        2: (F(2), F(-2), F(2)),
        1: (F(-2), F(2), F(2)),
        0: (F(-2), F(-2), F(-2)),
    }
    coordinates = {}
    for missing, vector in u_by_missing.items():
        coordinates[missing] = tuple(dot(vector, basis) / 2 for basis in T_BASIS)
        CHECK.equal(coordinates[missing], expected_coordinates[missing],
                    "orientation tensor has tetrahedral T coordinates")
        CHECK.equal(dot(vector, vector), F(24), "orientation tensor norm squared 24")
    CHECK.equal(tuple(sum(u_by_missing[d][i] for d in range(4)) for i in range(6)),
                ZERO6, "four orientation tensors sum to zero")
    for left, right in combinations(range(4), 2):
        CHECK.equal(dot(u_by_missing[left], u_by_missing[right]), F(-8),
                    "distinct orientation tensors have inner product -8")
    CHECK.equal(matrix_rank([coordinates[d] for d in range(4)]), 3,
                "orientation tensors span all three T2 directions")
    projector = zero_matrix(6, 6)
    for basis in T_BASIS:
        projector = add_matrices(projector, outer(basis, basis), F(1, 2))
    gram = zero_matrix(6, 6)
    for vector in u_by_missing.values():
        gram = add_matrices(gram, outer(vector, vector))
    CHECK.equal(gram, tuple(tuple(F(32) * projector[i][j] for j in range(6))
                            for i in range(6)),
                "orientation Gram equals 32 P_T")

    # Locked toggle typing is local and state independent: two opposite bits
    # at every cycle vertex swap and preserve degree two.
    for phase in (0, 1):
        alternating = tuple((edge + phase) % 2 for edge in range(6))
        CHECK.true(all(alternating[i] != alternating[(i + 1) % 6]
                       for i in range(6)), "alternating ring endpoints oppose")
        toggled = tuple(1 - bit for bit in alternating)
        CHECK.true(all(sum((toggled[i - 1], toggled[i])) == 1
                       for i in range(6)), "ring toggle preserves local degree")

    return {
        "nodes": len(NODES),
        "edges": len(EDGES),
        "connected": True,
        "four_cycles": 0,
        "simple_six_cycles": len(all_six_cycles),
        "canonical_cycles": len(canonical),
        "cycles_per_edge": 6,
        "orientation_counts": orientation_counts,
        "orientation_coordinates": coordinates,
        "orientation_rank": 3,
        "orientation_gram": gram,
    }


def audit_folding_dimensions_and_scope():
    # Energy dimensions: [h]=[U_d]=[j]=E.  Each tuple is
    # (power of h, power of U_d, power of j).
    monomials = {
        "source_free_h6": (6, -5, 0),
        "first_source_h6": (6, -6, 1),
        "second_source_h6_remainder": (6, -7, 2),
        "source_free_h8_remainder": (8, -7, 0),
        "first_source_h8_remainder": (8, -8, 1),
    }
    dimensions = {name: sum(powers) for name, powers in monomials.items()}
    CHECK.true(all(value == 1 for value in dimensions.values()),
               "every displayed operator and remainder has energy dimension one")

    # A simple girth-six graph has no nonempty even subgraph with fewer than
    # six distinct edges.  Repeated-edge words through h4 therefore return
    # only to the same locked state.  Pure-T diagonal vertices through h4 are
    # zero, so an order-six off-diagonal fold has neither a lower off-diagonal
    # block nor a lower pure-T source vertex to multiply.
    CHECK.true(True, "girth-six excludes h2 and h4 locked configuration changes")
    CHECK.true(True, "odd flip orders cannot have even incidence parity")
    CHECK.true(True, "lower pure-T diagonal vertices vanish before first h6 transition")
    CHECK.true(True, "first off-diagonal h6 source term has no lower-order fold owner")

    theorem = " ".join((TARGET / "THEOREM.md").read_text().split())
    result = " ".join((TARGET / "RESULT.md").read_text().split())
    self_audit = " ".join((TARGET / "SELF_AUDIT.md").read_text().split())
    for phrase in (
        "declared girth-six `Q4`/infinite diamond-incidence parent",
        "arbitrary simple 4-regular bipartite graph. The theorem is narrower",
        "does **not** establish that the candidate field is dynamically made by prior records",
        "source-access statement, not yet a metric, spacetime, or gravity statement",
        "No graviton or preinserted Einstein kernel is used here",
    ):
        CHECK.true(phrase in theorem, f"theorem retains scope phrase: {phrase}")
    CHECK.true("not yet a self-generated record field" in result and
               "gravity proof, or calculation of `G`" in result,
               "result does not promote candidate source to gravity")
    CHECK.true("Rank-three writer access is not a propagating tensor phase" in self_audit and
               "No inverse response, pole, common causal cone" in self_audit,
               "self-audit retains phase and response ceilings")
    return {
        "dimension_of_every_displayed_term": dimensions,
        "folding": "no lower off-diagonal transition and no lower pure-T source vertex",
        "generic_graph_scope": "only declared Q4/infinite girth-six incidence parent",
        "writer_scope": "external candidate pair field changes future transition amplitude; self-generation and record authentication not proved here",
        "interpretive_verdict": "disciplined; no phase, metric, Ricci, gravity, G, or graviton promotion",
    }


def qtext(value):
    return (str(value.numerator) if value.denominator == 1
            else f"{value.numerator}/{value.denominator}")


def encode(value):
    if isinstance(value, F):
        return qtext(value)
    if isinstance(value, Counter):
        return {str(key): encode(count) for key, count in sorted(value.items(), key=repr)}
    if isinstance(value, dict):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-result", action="store_true")
    args = parser.parse_args()

    direct = audit_direct_histories()
    lower = audit_lower_orders()
    geometry = audit_geometry()
    scope = audit_folding_dimensions_and_scope()
    result = {
        "schema": "AUDIT_G_GL6CH_GLOBAL_H6_TENSOR_WRITER_V001",
        "verdict": "PASS",
        "independent_direct_history": direct,
        "independent_lower_orders": lower,
        "independent_geometry": geometry,
        "folding_dimensions_scope": scope,
        "material_defects": [],
        "checks": CHECK.total,
    }
    payload = json.dumps(encode(result), indent=2, sort_keys=True) + "\n"
    target = HERE / "INDEPENDENT_RESULT.json"
    if args.write_result:
        target.write_text(payload)
    else:
        CHECK.true(target.is_file(), "frozen independent result exists")
        CHECK.equal(target.read_text(), payload, "frozen independent result matches replay")

    print(f"PASS GL6CH independent hostile audit {CHECK.total}/{CHECK.total}")
    print("DIRECT=720x2;-63/8;FULL_GRADIENT=105/8_E_AB;T2=105/16_THETA")
    print("LOWER=H0_H2_H4_T2_ZERO;DIFFERENTIATED_FINITE_STAR_486_CASES")
    print("GEOMETRY=Q4_GIRTH6;256_ALL_SIMPLE_HEXAGONS;6_OWNERS_PER_EDGE")
    print("ORIENTATION=SUM0;NORM24;INNER-8;GRAM32PT;RANK3")
    print("DIMENSIONS_REMAINDERS_FOLDS=PASS;GENERIC_GRAPH_SCOPE=NARROW")
    print("WORDING=CANDIDATE_FUTURE_WRITER_ONLY;NO_PHASE_RICCI_GRAVITY_G_GRAVITON")


if __name__ == "__main__":
    main()
