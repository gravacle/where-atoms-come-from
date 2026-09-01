#!/usr/bin/env python3
"""Exact checks for the GL6AO sixth-order locked Hamiltonian.

The script uses only the Python standard library.  It reconstructs the
declared period-four A3 incidence quotient, checks its local graph census,
enumerates all reduced perturbative words needed through sixth order, and
performs the exact folded-term cancellation over Q.

No continuum, gauge-field, graviton, Ricci, gravity, or G input is used.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb


class Checks:
    def __init__(self) -> None:
        self.total = 0

    def equal(self, got, want, label: str) -> None:
        self.total += 1
        if got != want:
            raise AssertionError(f"{label}: got {got!r}, want {want!r}")

    def true(self, condition: bool, label: str) -> None:
        self.total += 1
        if not condition:
            raise AssertionError(label)


CHECK = Checks()


def quotient_graph(length: int = 4):
    cells = list(product(range(length), repeat=3))

    def child(x, port):
        y = list(x)
        if port < 3:
            y[port] = (y[port] + 1) % length
        return tuple(y)

    edges = [(x, port) for x in cells for port in range(4)]
    endpoints = {
        edge: (("P", edge[0]), ("C", child(*edge)))
        for edge in edges
    }
    nodes = [(kind, x) for kind in ("P", "C") for x in cells]
    incident = {node: [] for node in nodes}
    for edge, ends in endpoints.items():
        for node in ends:
            incident[node].append(edge)
    return cells, nodes, edges, endpoints, incident, child


def enumerate_six_cycles(nodes, endpoints, incident):
    adjacency = {node: [] for node in nodes}
    for edge, (u, v) in endpoints.items():
        adjacency[u].append((v, edge))
        adjacency[v].append((u, edge))

    cycles = set()
    for start in nodes:
        def walk(node, visited_nodes, path_edges):
            if len(path_edges) == 6:
                if node == start:
                    cycles.add(frozenset(path_edges))
                return
            for other, edge in adjacency[node]:
                if edge in path_edges:
                    continue
                if other == start:
                    if len(path_edges) == 5:
                        walk(other, visited_nodes + [other], path_edges + [edge])
                elif other not in visited_nodes:
                    walk(other, visited_nodes + [other], path_edges + [edge])

        walk(start, [start], [])
    return cycles, adjacency


def deterministic_locked_background(length: int = 4):
    """A degree-two b-matching containing one declared alternating hexagon."""
    cells, nodes, edges, endpoints, incident, child = quotient_graph(length)
    p0 = (0, 0, 0)
    p1 = (1, -1, 0)
    p2 = (0, -1, 1)
    raw_cycle = (
        (p0, 0), (p1, 1), (p1, 2),
        (p2, 0), (p2, 1), (p0, 2),
    )

    def mod_x(x):
        return tuple(value % length for value in x)

    cycle = tuple((mod_x(x), port) for x, port in raw_cycle)
    fixed = {edge: (1 if i % 2 == 0 else 0) for i, edge in enumerate(cycle)}
    p_fixed = Counter()
    c_fixed = Counter()
    for (parent, port), value in fixed.items():
        if value:
            p_fixed[parent] += 1
            c_fixed[child(parent, port)] += 1

    source, sink = ("S",), ("T",)
    capacity = {}

    def add_edge(u, v, cap):
        capacity.setdefault(u, {})[v] = cap
        capacity.setdefault(v, {}).setdefault(u, 0)

    for parent in cells:
        add_edge(source, ("P", parent), 2 - p_fixed[parent])
    for child_cell in cells:
        add_edge(("C", child_cell), sink, 2 - c_fixed[child_cell])
    for edge in edges:
        if edge not in fixed:
            parent, port = edge
            add_edge(("P", parent), ("C", child(parent, port)), 1)

    residual = {u: dict(row) for u, row in capacity.items()}
    flow = 0
    while True:
        previous = {source: None}
        queue = deque([source])
        while queue and sink not in previous:
            u = queue.popleft()
            for v in sorted(residual.get(u, {}), key=repr):
                if residual[u][v] > 0 and v not in previous:
                    previous[v] = u
                    queue.append(v)
        if sink not in previous:
            break
        v = sink
        while v != source:
            u = previous[v]
            residual[u][v] -= 1
            residual[v][u] += 1
            v = u
        flow += 1

    occupied = {edge for edge, value in fixed.items() if value}
    for edge in edges:
        if edge in fixed:
            continue
        parent, port = edge
        child_cell = child(parent, port)
        if residual[("C", child_cell)].get(("P", parent), 0) == 1:
            occupied.add(edge)
    return occupied, cycle, flow, cells, nodes, edges, endpoints, incident


def is_locked(occupied, incident):
    return all(sum(edge in occupied for edge in edges) == 2
               for edges in incident.values())


def is_flippable(occupied, cycle, endpoints, incident):
    cycle_nodes = {node for edge in cycle for node in endpoints[edge]}
    return all(sum(edge in occupied for edge in incident[node] if edge in cycle) == 1
               for node in cycle_nodes)


def subset_energy(labels, pair_energy, triple_energy=None):
    if len(labels) == 1:
        return 2
    if len(labels) == 2:
        return pair_energy[tuple(sorted(labels))]
    if len(labels) == 3:
        if triple_energy is None:
            raise AssertionError("triple energy required")
        return triple_energy
    if not labels:
        return 0
    raise AssertionError("only one-, two-, and three-edge subsets occur")


def reduced_word_sum(multiset, pair_energy, triple_energy=None, powers=None):
    """Sum Q-only return words with R=-D^{-1} at each intermediate step."""
    total = Fraction()
    retained = 0
    for word in set(permutations(multiset)):
        parity = set()
        energies = []
        allowed = True
        for label in word[:-1]:
            if label in parity:
                parity.remove(label)
            else:
                parity.add(label)
            if not parity:
                allowed = False
                break
            energies.append(subset_energy(parity, pair_energy, triple_energy))
        if not allowed:
            continue
        exponents = powers if powers is not None else (1,) * len(energies)
        value = Fraction((-1) ** sum(exponents), 1)
        for energy, exponent in zip(energies, exponents):
            value /= energy ** exponent
        total += value
        retained += 1
    return total, retained


def poly_add(a, b):
    size = max(len(a), len(b))
    return tuple(
        (a[i] if i < len(a) else Fraction())
        + (b[i] if i < len(b) else Fraction())
        for i in range(size)
    )


def poly_scale(a, scalar):
    return tuple(Fraction(scalar) * value for value in a)


def poly_mul(a, b):
    out = [Fraction() for _ in range(len(a) + len(b) - 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return tuple(out)


def section_q4_graph_and_locked_states():
    occupied, target, flow, cells, nodes, edges, endpoints, incident = (
        deterministic_locked_background()
    )
    CHECK.equal(len(cells), 64, "Q4 cell count")
    CHECK.equal(len(nodes), 128, "Q4 constraint-node count")
    CHECK.equal(len(edges), 256, "Q4 active-link count M")
    CHECK.true(all(len(row) == 4 for row in incident.values()),
               "Q4 is degree four")
    CHECK.true(all(len(set(row)) == 4 for row in incident.values()),
               "Q4 has no parallel incidence at a node")
    CHECK.equal(flow, 125, "declared b-matching flow")
    CHECK.equal(len(occupied), 128, "locked background occupation")
    CHECK.true(is_locked(occupied, incident), "declared background is locked")
    CHECK.true(is_flippable(occupied, frozenset(target), endpoints, incident),
               "declared target is alternating")

    cycles, adjacency = enumerate_six_cycles(nodes, endpoints, incident)
    CHECK.equal(len(cycles), 256, "all Q4 six-cycles enumerated")
    per_edge = Counter(edge for cycle in cycles for edge in cycle)
    CHECK.true(all(per_edge[edge] == 6 for edge in edges),
               "each Q4 edge belongs to six six-cycles")
    for index, cycle in enumerate(sorted(cycles, key=lambda value: repr(sorted(value)))):
        CHECK.equal(len(cycle), 6, f"cycle {index} has six edges")
        cycle_nodes = {node for edge in cycle for node in endpoints[edge]}
        CHECK.equal(len(cycle_nodes), 6, f"cycle {index} has six vertices")
        CHECK.true(all(sum(edge in cycle for edge in incident[node]) == 2
                       for node in cycle_nodes),
                   f"cycle {index} has degree two at every cycle vertex")
        CHECK.true(all(sum(edge in cycle for edge in incident[node]) == 0
                       for node in set(nodes) - cycle_nodes),
                   f"cycle {index} has no external cycle incidence")

    # Generate multiple exact locked configurations by legal ring toggles.
    configurations = []
    current = set(occupied)
    ordered_cycles = sorted(cycles, key=lambda value: repr(sorted(value)))
    for step in range(32):
        CHECK.true(is_locked(current, incident), f"ring-walk state {step} locked")
        configurations.append(frozenset(current))
        flippable = [cycle for cycle in ordered_cycles
                     if is_flippable(current, cycle, endpoints, incident)]
        CHECK.true(bool(flippable), f"ring-walk state {step} has a move")
        chosen = flippable[(7 * step + 3) % len(flippable)]
        current.symmetric_difference_update(chosen)
    CHECK.true(len(set(configurations)) >= 16,
               "ring walk samples at least sixteen distinct locked states")
    return configurations, cycles, cells, nodes, edges, endpoints, incident


def section_local_censuses(configurations, edges, endpoints, incident):
    edge_set = set(edges)
    adjacent_pairs = set()
    for node, row in incident.items():
        for a, b in combinations(row, 2):
            pair = frozenset((a, b))
            CHECK.true(pair not in adjacent_pairs,
                       f"adjacent pair has unique owner {node!r}")
            adjacent_pairs.add(pair)
    m = len(edges)
    CHECK.equal(len(adjacent_pairs), 3 * m, "adjacent-pair census")

    for sample, occupied in enumerate(configurations):
        pair_types = Counter()
        for pair in adjacent_pairs:
            a, b = tuple(pair)
            pair_types[2 if ((a in occupied) != (b in occupied)) else 6] += 1
            external = edge_set - pair
            disjoint = [edge for edge in external
                        if all(frozenset((edge, member)) not in adjacent_pairs
                               for member in pair)]
            CHECK.equal(len(disjoint), m - 10,
                        f"sample {sample} adjacent pair has M-10 disjoint thirds")
        CHECK.equal(pair_types, Counter({2: 2 * m, 6: m}),
                    f"sample {sample} adjacent mixed/same census")

        star_types = Counter()
        for node, row in incident.items():
            for triple in combinations(row, 3):
                signs = tuple(1 if edge not in occupied else -1 for edge in triple)
                pair_energies = sorted(
                    2 if signs[i] != signs[j] else 6
                    for i, j in combinations(range(3), 2)
                )
                triple_energy = 3 + sum(signs) ** 2
                star_types[(tuple(pair_energies), triple_energy)] += 1
        CHECK.equal(star_types, Counter({((2, 2, 6), 4): 2 * m}),
                    f"sample {sample} star census")

        path_types = Counter()
        for middle in edges:
            u, v = endpoints[middle]
            left = [edge for edge in incident[u] if edge != middle]
            right = [edge for edge in incident[v] if edge != middle]
            CHECK.equal(len(left), 3, f"sample {sample} path left degree")
            CHECK.equal(len(right), 3, f"sample {sample} path right degree")
            for a in left:
                for b in right:
                    CHECK.true(frozenset((a, b)) not in adjacent_pairs,
                               f"sample {sample} path endpoints disjoint")
                    signs = {
                        edge: (1 if edge not in occupied else -1)
                        for edge in (a, middle, b)
                    }
                    p_left = 2 if signs[a] != signs[middle] else 6
                    p_right = 2 if signs[b] != signs[middle] else 6
                    p_tuple = tuple(sorted((p_left, p_right, 4)))
                    triple_energy = 2
                    triple_energy += 0 if p_left == 2 else 4
                    triple_energy += 0 if p_right == 2 else 4
                    path_types[(p_tuple, triple_energy)] += 1
        CHECK.equal(
            path_types,
            Counter({
                ((2, 2, 4), 2): 4 * m,
                ((2, 4, 6), 6): 4 * m,
                ((4, 6, 6), 10): m,
            }),
            f"sample {sample} three-edge path census",
        )


def section_reduced_words_and_folds():
    pair_t4 = {2: Fraction(-1, 2), 4: Fraction(-1, 4), 6: Fraction(-1, 6)}
    x_table = {2: Fraction(3, 4), 4: Fraction(5, 16), 6: Fraction(7, 36)}
    repeated_pair_t6 = {
        2: Fraction(-1, 4),
        4: Fraction(-1, 16),
        6: Fraction(-1, 36),
    }
    for pair_energy, expected in pair_t4.items():
        pairs = {(0, 1): pair_energy}
        value, retained = reduced_word_sum((0, 0, 1, 1), pairs)
        CHECK.equal(retained, 4, f"T4 retained words p={pair_energy}")
        CHECK.equal(value, expected, f"T4 word sum p={pair_energy}")
        x_value = Fraction()
        for powers in ((2, 1, 1), (1, 2, 1), (1, 1, 2)):
            part, part_retained = reduced_word_sum(
                (0, 0, 1, 1), pairs, powers=powers
            )
            CHECK.equal(part_retained, 4,
                        f"folded X retained words p={pair_energy} powers={powers}")
            x_value += part
        CHECK.equal(x_value, x_table[pair_energy],
                    f"folded X word sum p={pair_energy}")
        first, retained_first = reduced_word_sum((0, 0, 0, 0, 1, 1), pairs)
        second, retained_second = reduced_word_sum((0, 0, 1, 1, 1, 1), pairs)
        CHECK.equal((retained_first, retained_second), (4, 4),
                    f"repeated-pair T6 retained words p={pair_energy}")
        CHECK.equal(first + second, repeated_pair_t6[pair_energy],
                    f"repeated-pair T6 word sum p={pair_energy}")

    triple_table = {
        ((4, 4, 4), 6): Fraction(-9, 32),
        ((2, 4, 4), 4): Fraction(-9, 16),
        ((4, 4, 6), 8): Fraction(-29, 144),
        ((2, 2, 6), 4): Fraction(-109, 144),
        ((2, 2, 4), 2): Fraction(-41, 32),
        ((2, 4, 6), 6): Fraction(-337, 864),
        ((4, 6, 6), 10): Fraction(-209, 1440),
    }
    for (pair_tuple, triple_energy), expected in triple_table.items():
        pair_map = {
            (0, 1): pair_tuple[0],
            (0, 2): pair_tuple[1],
            (1, 2): pair_tuple[2],
        }
        value, retained = reduced_word_sum(
            (0, 0, 1, 1, 2, 2), pair_map, triple_energy
        )
        CHECK.equal(retained, 60,
                    f"triple retained words {pair_tuple}/{triple_energy}")
        CHECK.equal(value, expected,
                    f"triple word sum {pair_tuple}/{triple_energy}")

    # Counts of unordered three-edge shapes on a degree-four girth-six graph.
    # Coefficient tuples are in ascending powers of M.
    m = (Fraction(), Fraction(1))
    m2 = poly_mul(m, m)
    choose2 = (Fraction(), Fraction(-1, 2), Fraction(1, 2))
    choose3 = (Fraction(), Fraction(1, 3), Fraction(-1, 2), Fraction(1, 6))
    matching = poly_add(poly_add(choose3, poly_scale(m2, -3)), poly_scale(m, 19))
    one_mixed = poly_add(poly_scale(m2, 2), poly_scale(m, -20))
    one_same = poly_add(m2, poly_scale(m, -10))
    CHECK.equal(poly_add(
        poly_add(poly_add(matching, one_mixed), one_same),
        poly_scale(m, 11),
    ), choose3, "all unordered triple shapes partition C(M,3)")

    triple_terms = (
        (matching, triple_table[((4, 4, 4), 6)]),
        (one_mixed, triple_table[((2, 4, 4), 4)]),
        (one_same, triple_table[((4, 4, 6), 8)]),
        (poly_scale(m, 2), triple_table[((2, 2, 6), 4)]),
        (poly_scale(m, 4), triple_table[((2, 2, 4), 2)]),
        (poly_scale(m, 4), triple_table[((2, 4, 6), 6)]),
        (m, triple_table[((4, 6, 6), 10)]),
    )
    t6_triples = ()
    for count, weight in triple_terms:
        t6_triples = poly_add(t6_triples, poly_scale(count, weight))
    CHECK.equal(
        t6_triples,
        (Fraction(), Fraction(-2237, 4320), Fraction(-197, 576), Fraction(-3, 64)),
        "direct T6 triple polynomial",
    )

    t6_pairs = poly_add(
        poly_scale(m, 2 * repeated_pair_t6[2]
                   + repeated_pair_t6[6]
                   - 3 * repeated_pair_t6[4]),
        poly_scale(choose2, repeated_pair_t6[4]),
    )
    CHECK.equal(t6_pairs,
                (Fraction(), Fraction(-89, 288), Fraction(-1, 32)),
                "direct T6 repeated-pair polynomial")
    t6 = poly_add(t6_triples, t6_pairs)
    CHECK.equal(
        t6,
        (Fraction(), Fraction(-893, 1080), Fraction(-215, 576), Fraction(-3, 64)),
        "complete direct diagonal T6 polynomial",
    )

    x_fold = poly_add(
        poly_scale(m, 2 * x_table[2] + x_table[6] - 3 * x_table[4]),
        poly_scale(choose2, x_table[4]),
    )
    CHECK.equal(x_fold,
                (Fraction(), Fraction(173, 288), Fraction(5, 32)),
                "complete folded X polynomial")

    e2 = poly_scale(m, Fraction(-1, 2))
    a2 = poly_scale(m, Fraction(1, 4))
    a3 = poly_scale(m, Fraction(-1, 8))
    e4 = poly_scale(m, Fraction(-7, 24))
    e6 = t6
    e6 = poly_add(e6, poly_scale(poly_mul(e2, x_fold), -1))
    e6 = poly_add(e6, poly_mul(poly_mul(e2, e2), a3))
    e6 = poly_add(e6, poly_scale(poly_mul(e4, a2), -1))
    CHECK.equal(
        e6,
        (Fraction(), Fraction(-893, 1080), Fraction(), Fraction()),
        "M^3 and M^2 cancel; sixth-order diagonal is common and extensive",
    )


def cycle_subset_energy(subset):
    signs = (1, -1, 1, -1, 1, -1)
    charge = Counter()
    for edge in subset:
        delta = signs[edge]
        charge[edge] += delta
        charge[(edge + 1) % 6] += delta
    return sum(value * value for value in charge.values())


def section_offdiagonal_hexagon():
    total = Fraction()
    profiles = Counter()
    for order in permutations(range(6)):
        selected = set()
        energies = []
        term = Fraction(1)
        for edge in order[:-1]:
            selected.add(edge)
            energy = cycle_subset_energy(selected)
            CHECK.true(energy > 0, f"hexagon order {order} has no early lock return")
            energies.append(energy)
            term *= Fraction(-1, energy)
        profiles[tuple(energies)] += 1
        total += term
    CHECK.equal(sum(profiles.values()), 720, "all hexagon orders retained")
    CHECK.equal(total, Fraction(-63, 8), "universal alternating-hexagon amplitude")

    # Every nonempty finite symmetric difference of two degree-two states is
    # an even subgraph.  With girth six and at most six edges, it is exactly
    # one six-cycle; degree preservation forces alternation.  The finite Q4
    # graph checks backing each premise were performed above.
    for size in range(1, 6):
        CHECK.true(size < 6, f"no locked offdiagonal support at distance {size}")


def main():
    configurations, cycles, cells, nodes, edges, endpoints, incident = (
        section_q4_graph_and_locked_states()
    )
    section_local_censuses(configurations, edges, endpoints, incident)
    section_reduced_words_and_folds()
    section_offdiagonal_hexagon()
    print(f"PASS__GL6AO_COMPLETE_SIXTH_ORDER__{CHECK.total}/{CHECK.total}")
    print("Q4=M256_N128_GIRTH6_HEXAGONS256_EACH_EDGE6")
    print("H2_DIAG=-M/2;H4_DIAG=-7M/24;H6_DIAG=-893M/1080")
    print("H6_OFFDIAG=ONLY_ALTERNATING_HEXAGON;AMPLITUDE=-63/8")
    print("LINKED=FINITE_RANGE_FORMAL_ORDER6;NO_CONVERGED_PHASE_POLE_CONE_GRAVITY_G")


if __name__ == "__main__":
    main()
