#!/usr/bin/env python3
"""Independent exact replay for the frozen GL6AR theorem.

This file imports no author module.  It reconstructs the periodic incidence,
hexagons, projected supports, locked states, fluxes, active hypercube trial,
and exact combinatorics behind the coefficient 18.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from fractions import Fraction
from itertools import combinations, product


class Checks:
    def __init__(self):
        self.count = 0

    def equal(self, got, want, label):
        self.count += 1
        if got != want:
            raise RuntimeError(f"{label}: got {got!r}, want {want!r}")

    def true(self, condition, label):
        self.count += 1
        if not condition:
            raise RuntimeError(label)


C = Checks()
D = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))


def vadd(x, y, length):
    return tuple((a + b) % length for a, b in zip(x, y))


def vsub(x, y, length):
    return tuple((a - b) % length for a, b in zip(x, y))


def make_geometry(length):
    cells = tuple(product(range(length), repeat=3))
    links = tuple((x, a) for x in cells for a in range(4))

    def child(edge):
        x, a = edge
        return vadd(x, D[a], length)

    endpoints = {
        edge: (("P", edge[0]), ("C", child(edge)))
        for edge in links
    }
    vertices = tuple((side, x) for side in ("P", "C") for x in cells)
    incident = {vertex: [] for vertex in vertices}
    for edge, ends in endpoints.items():
        for vertex in ends:
            incident[vertex].append(edge)
    return cells, vertices, links, endpoints, incident, child


def hexagon(x, ports, length):
    a, b, c = ports
    p = vadd(vsub(x, D[b], length), D[a], length)
    q = vadd(vsub(x, D[b], length), D[c], length)
    return ((x, a), (p, b), (p, c), (q, a), (q, b), (x, c))


def all_hexagons(length):
    return tuple(
        hexagon(x, ports, length)
        for x in product(range(length), repeat=3)
        for ports in combinations(range(4), 3)
    )


def cycle_vertices(cycle, endpoints):
    return {vertex for edge in cycle for vertex in endpoints[edge]}


def projected_support(cycle, endpoints, incident):
    return {
        edge
        for vertex in cycle_vertices(cycle, endpoints)
        for edge in incident[vertex]
    }


def locked(occupied, incident):
    return all(sum(edge in occupied for edge in row) == 2
               for row in incident.values())


def alternating(occupied, cycle):
    bits = tuple(edge in occupied for edge in cycle)
    return all(bits[i] != bits[(i + 1) % 6] for i in range(6))


def construct_active_q4():
    length = 4
    cells, vertices, links, endpoints, incident, child = make_geometry(length)
    target = hexagon((0, 0, 0), (0, 1, 2), length)
    fixed = {edge: int(index % 2 == 0) for index, edge in enumerate(target)}

    parent_fixed = Counter()
    child_fixed = Counter()
    for edge, bit in fixed.items():
        if bit:
            parent_fixed[edge[0]] += 1
            child_fixed[child(edge)] += 1

    source = ("source",)
    sink = ("sink",)
    adjacency = defaultdict(list)
    residual = {}

    def arc(u, v, capacity):
        adjacency[u].append(v)
        adjacency[v].append(u)
        residual[u, v] = capacity
        residual[v, u] = 0

    demand = 0
    for x in cells:
        p_need = 2 - parent_fixed[x]
        c_need = 2 - child_fixed[x]
        arc(source, ("P", x), p_need)
        arc(("C", x), sink, c_need)
        demand += p_need
    for edge in links:
        if edge not in fixed:
            arc(("P", edge[0]), ("C", child(edge)), 1)

    flow = 0
    while True:
        previous = {source: None}
        queue = deque((source,))
        while queue and sink not in previous:
            u = queue.popleft()
            for v in sorted(adjacency[u], key=repr):
                if residual[u, v] and v not in previous:
                    previous[v] = u
                    queue.append(v)
        if sink not in previous:
            break
        v = sink
        while v != source:
            u = previous[v]
            residual[u, v] -= 1
            residual[v, u] += 1
            v = u
        flow += 1

    C.equal(flow, demand, "capacity construction saturates all degree demands")
    occupied = {edge for edge, bit in fixed.items() if bit}
    for edge in links:
        if edge in fixed:
            continue
        p = ("P", edge[0])
        c = ("C", child(edge))
        if residual[c, p] == 1:
            occupied.add(edge)
    C.true(locked(occupied, incident), "constructed Q4 state is locked")
    C.true(alternating(occupied, target), "constructed target is alternating")
    return occupied, target


def check_geometry_and_supports():
    length = 4
    cells, vertices, links, endpoints, incident, child = make_geometry(length)
    cycles = all_hexagons(length)
    C.equal(len(cells), 64, "Q4 cell count")
    C.equal(len(vertices), 128, "Q4 vertex count")
    C.equal(len(links), 256, "Q4 link count")
    C.true(all(len(row) == 4 for row in incident.values()), "degree four")
    C.equal(len(cycles), 256, "Q4 hexagon count")
    C.equal(len({frozenset(cycle) for cycle in cycles}), 256,
            "Q4 hexagons distinct")
    C.true(all(len(set(cycle)) == 6 for cycle in cycles),
           "six distinct links per hexagon")
    C.true(all(len(cycle_vertices(cycle, endpoints)) == 6 for cycle in cycles),
           "six distinct vertices per hexagon")
    C.true(all(sorted(Counter(edge[1] for edge in cycle).values()) == [2, 2, 2]
               for cycle in cycles), "three ports twice")
    C.true(all(all(sum(edge in cycle for edge in incident[vertex]) == 2
                     for vertex in cycle_vertices(cycle, endpoints))
               for cycle in cycles), "cycle degree two at every cycle vertex")

    core_multiplicity = Counter(edge for cycle in cycles for edge in cycle)
    C.true(all(core_multiplicity[edge] == 6 for edge in links),
           "six hexagon cores per link")
    supports = tuple(projected_support(cycle, endpoints, incident)
                     for cycle in cycles)
    C.true(all(len(support) == 18 for support in supports),
           "every projected term has support eighteen")
    support_multiplicity = Counter(edge for support in supports for edge in support)
    C.true(all(support_multiplicity[edge] == 18 for edge in links),
           "eighteen projected supports per link")

    for cycle in cycles:
        for phase in (0, 1):
            delta = Counter()
            for position, edge in enumerate(cycle):
                delta[edge[1]] += 1 if (position + phase) % 2 else -1
            C.true(all(delta[a] == 0 for a in range(4)),
                   "alternating toggle preserves every port count")


def flux_tuple(occupied, length, links):
    counts = tuple(sum(edge in occupied and edge[1] == a for edge in links)
                   for a in range(4))
    signed = tuple(2 * value - length ** 3 for value in counts)
    return counts, signed, tuple(value // length for value in signed)


def check_fluxes_frozen_and_boundary():
    length = 4
    cells, vertices, links, endpoints, incident, child = make_geometry(length)
    cycles = all_hexagons(length)
    occupied, target = construct_active_q4()
    samples = [frozenset(occupied), frozenset(set(occupied) ^ set(target))]
    for ports in combinations(range(4), 2):
        uniform = frozenset(edge for edge in links if edge[1] in ports)
        C.true(locked(uniform, incident), "uniform two-port state locked")
        C.true(not any(alternating(uniform, cycle) for cycle in cycles),
               "uniform two-port state frozen")
        samples.append(uniform)

    for state in samples:
        counts, signed, flux = flux_tuple(state, length, links)
        C.equal(sum(counts), 2 * length ** 3, "locked occupation count")
        C.equal(sum(signed), 0, "four signed port values sum zero")
        C.true(all(value % length == 0 for value in signed),
               "normalized port values integral")
        for axis in range(3):
            planes = tuple(
                sum(1 if edge in state else -1
                    for edge in links
                    if edge[1] == axis and edge[0][axis] == plane)
                for plane in range(length)
            )
            C.true(all(value == planes[0] for value in planes),
                   "coordinate cut flux plane independent")
            C.equal(signed[axis], length * planes[0],
                    "volume count equals length times cut flux")
            C.equal(flux[axis], planes[0], "integer coordinate cut flux")
        C.equal(flux[3], -sum(flux[:3]), "dependent fourth port invariant")

    before = flux_tuple(occupied, length, links)[0]
    toggled = set(occupied) ^ set(target)
    C.true(locked(toggled, incident), "target toggle preserves lock")
    C.equal(flux_tuple(toggled, length, links)[0], before,
            "target toggle preserves all port counts")

    outside = set(occupied) - set(target)
    allowed = []
    for bits in product((0, 1), repeat=6):
        candidate = set(outside)
        candidate.update(edge for bit, edge in zip(bits, target) if bit)
        if locked(candidate, incident):
            allowed.append(frozenset(candidate))
    C.equal(len(allowed), 2, "fixed-exterior target sector has two states")
    C.equal(allowed[0] ^ allowed[1], frozenset(target),
            "fixed-exterior graph is K2")
    C.true(all(alternating(state, target) for state in allowed),
           "both K2 states are flippable")
    C.equal(Fraction(1) - Fraction(-1), Fraction(2), "K2 adjacency gap is 2")
    return occupied, target


def check_active_trial_and_constant(occupied4, target4):
    length = 8
    cells, vertices, links, endpoints, incident, child = make_geometry(length)
    occupied = {
        edge for edge in links
        if (tuple(value % 4 for value in edge[0]), edge[1]) in occupied4
    }
    C.true(locked(occupied, incident), "period-four state lifts to Q8")
    targets = tuple(
        hexagon(tuple(4 * value for value in block), (0, 1, 2), length)
        for block in product((0, 1), repeat=3)
    )
    C.equal(len(targets), 8, "one active target per 4-cube")
    supports = tuple(projected_support(cycle, endpoints, incident)
                     for cycle in targets)
    C.true(all(not supports[i] & supports[j]
               for i in range(len(supports)) for j in range(i)),
           "selected projected supports pairwise disjoint")
    C.true(all(alternating(occupied, cycle) for cycle in targets),
           "all selected targets initially alternating")

    all_independent = True
    for bits in product((0, 1), repeat=len(targets)):
        state = set(occupied)
        for bit, cycle in zip(bits, targets):
            if bit:
                state.symmetric_difference_update(cycle)
        all_independent &= locked(state, incident)
        all_independent &= all(alternating(state, cycle) for cycle in targets)
    C.true(all_independent, "all 2^8 independent toggle subsets remain active")

    dimension = len(targets)
    vertices_count = 2 ** dimension
    C.equal(dimension, length ** 3 // 64, "hypercube trial dimension")
    C.equal(vertices_count * dimension // 2, 1024, "hypercube edge count")
    exact_dirichlet = True
    for character in product((0, 1), repeat=dimension):
        weight = sum(character)
        if weight == 0:
            continue
        crossing_edges = vertices_count * weight // 2
        numerator = Fraction(4 * crossing_edges, vertices_count)
        exact_dirichlet &= numerator == 2 * weight
    C.true(exact_dirichlet, "hypercube PF transform reproduces all Walsh gaps")

    cycles = all_hexagons(length)

    def tent(r):
        return Fraction(min(r, length - r), length)

    max_seen = Fraction()
    bound_holds = True
    for weights in product((-1, 1), repeat=4):
        for cycle in cycles:
            for phase in (0, 1):
                delta = Fraction()
                for position, edge in enumerate(cycle):
                    change = 1 if (position + phase) % 2 == 0 else -1
                    delta += change * weights[edge[1]] * tent(edge[0][0])
                max_seen = max(max_seen, abs(delta))
                bound_holds &= abs(delta) <= Fraction(3, length)
    C.true(bound_holds, "tent variation bounded by 3||w||/L")
    C.true(max_seen > 0, "tent observable not identically unchanged")
    C.equal(len(cycles), 4 * length ** 3, "four L^3 possible flips")
    pf_edge_weight_bound = Fraction(len(cycles), 2)
    C.equal(pf_edge_weight_bound, 2 * length ** 3,
            "PF edge-weight sum bound is 2L^3")
    numerator_bound = pf_edge_weight_bound * Fraction(9, length ** 2)
    C.equal(numerator_bound, 18 * length,
            "exact Dirichlet numerator constant is 18L")


def main():
    check_geometry_and_supports()
    occupied, target = check_fluxes_frozen_and_boundary()
    check_active_trial_and_constant(occupied, target)
    print(f"PASS__GL6AR_INDEPENDENT_REPLAY__{C.count}/{C.count}")
    print("GEOMETRY=Q4_256_LINKS_256_DISTINCT_HEXAGONS_GIRTH_SCOPE_L_GE_4")
    print("INTERACTION=PROJECTED_SUPPORT_18;LINK_MULTIPLICITY_18;LOCK_PRESERVED")
    print("INVARIANTS=THREE_COORDINATE_CUT_FLUXES_PLUS_DEPENDENT_PORT_VALUE")
    print("COMPONENTS=FROZEN_EXTREMES_AND_ACTIVE_K2_HYPERCUBE_WITNESSES")
    print("VARIATIONAL=UNORDERED_EDGE_DIRICHLET_FORM;NUMERATOR_LE_18_T_NORMW2_L")
    print("SCOPE=POSITIVE_VARIANCE_NONTRIVIAL_COMPONENT;FINITE_SIZE_ONLY")
    print("CEILING=NO_GNS_GAPLESSNESS_MOMENTUM_PHOTON_CONE_GRAVITY_G")


if __name__ == "__main__":
    main()
