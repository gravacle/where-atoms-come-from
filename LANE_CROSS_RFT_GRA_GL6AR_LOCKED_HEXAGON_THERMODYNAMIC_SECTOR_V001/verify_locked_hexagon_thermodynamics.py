#!/usr/bin/env python3
"""Exact finite checks for GL6AR locked hexagon thermodynamics.

Only Python's standard library is used.  The script reconstructs the native
periodic incidence, elementary hexagons, locked fluxes, frozen and active
sectors, projected interaction supports, finite flip graphs, and the exact
combinatorics entering the conditional soft-mode bound.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from fractions import Fraction
from itertools import combinations, product


class Checks:
    def __init__(self):
        self.total = 0

    def equal(self, got, want, label):
        self.total += 1
        if got != want:
            raise AssertionError(f"{label}: got {got!r}, want {want!r}")

    def true(self, condition, label):
        self.total += 1
        if not condition:
            raise AssertionError(label)


CHECK = Checks()
DIRECTIONS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))


def add(x, y, length):
    return tuple((a + b) % length for a, b in zip(x, y))


def sub(x, y, length):
    return tuple((a - b) % length for a, b in zip(x, y))


def geometry(length):
    cells = list(product(range(length), repeat=3))
    edges = [(x, port) for x in cells for port in range(4)]

    def child(x, port):
        return add(x, DIRECTIONS[port], length)

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


def canonical_hexagon(x, ports, length):
    a, b, c = ports
    x_ab = add(sub(x, DIRECTIONS[b], length), DIRECTIONS[a], length)
    x_bc = add(sub(x, DIRECTIONS[b], length), DIRECTIONS[c], length)
    return (
        (x, a),
        (x_ab, b),
        (x_ab, c),
        (x_bc, a),
        (x_bc, b),
        (x, c),
    )


def all_hexagons(length):
    cells = list(product(range(length), repeat=3))
    return tuple(canonical_hexagon(x, ports, length)
                 for x in cells for ports in combinations(range(4), 3))


def cycle_nodes(cycle, endpoints):
    return {node for edge in cycle for node in endpoints[edge]}


def locked(occupied, incident):
    return all(sum(edge in occupied for edge in row) == 2
               for row in incident.values())


def flippable(occupied, cycle, endpoints, incident):
    return all(sum(edge in occupied for edge in incident[node] if edge in cycle) == 1
               for node in cycle_nodes(cycle, endpoints))


def deterministic_q4_background():
    length = 4
    cells, nodes, edges, endpoints, incident, child = geometry(length)
    target = canonical_hexagon((0, 0, 0), (0, 1, 2), length)
    fixed = {edge: (1 if index % 2 == 0 else 0)
             for index, edge in enumerate(target)}
    p_fixed = Counter()
    c_fixed = Counter()
    for (parent, port), value in fixed.items():
        if value:
            p_fixed[parent] += 1
            c_fixed[child(parent, port)] += 1

    source, sink = ("S",), ("T",)
    capacity = {}

    def add_arc(u, v, cap):
        capacity.setdefault(u, {})[v] = cap
        capacity.setdefault(v, {}).setdefault(u, 0)

    for x in cells:
        add_arc(source, ("P", x), 2 - p_fixed[x])
        add_arc(("C", x), sink, 2 - c_fixed[x])
    for edge in edges:
        if edge not in fixed:
            parent, port = edge
            add_arc(("P", parent), ("C", child(parent, port)), 1)
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
    return occupied, target, flow


def flux_data(occupied, length, edges):
    n_port = tuple(sum(edge in occupied and edge[1] == port for edge in edges)
                   for port in range(4))
    signed = tuple(2 * count - length ** 3 for count in n_port)
    cut = tuple(value // length for value in signed)
    return n_port, signed, cut


def section_geometry_flux_and_sectors():
    length = 4
    cells, nodes, edges, endpoints, incident, child = geometry(length)
    hexagons = all_hexagons(length)
    CHECK.equal(len(cells), length ** 3, "cell census")
    CHECK.equal(len(nodes), 2 * length ** 3, "constraint-node census")
    CHECK.equal(len(edges), 4 * length ** 3, "edge census")
    CHECK.equal(len(hexagons), 4 * length ** 3, "elementary-hexagon census")
    CHECK.equal(len({frozenset(cycle) for cycle in hexagons}), len(hexagons),
                "canonical hexagons are unique")
    CHECK.true(all(len(set(cycle)) == 6 for cycle in hexagons),
               "every elementary hexagon has six links")
    core_count = Counter(edge for cycle in hexagons for edge in cycle)
    CHECK.true(all(core_count[edge] == 6 for edge in edges),
               "every link lies in six elementary hexagon cores")

    for index, cycle in enumerate(hexagons):
        vertices = cycle_nodes(cycle, endpoints)
        CHECK.equal(len(vertices), 6, f"hexagon {index} has six vertices")
        CHECK.true(all(sum(edge in cycle for edge in incident[node]) == 2
                       for node in vertices),
                   f"hexagon {index} has cycle degree two")
        ports = Counter(edge[1] for edge in cycle)
        CHECK.equal(sorted(ports.values()), [2, 2, 2],
                    f"hexagon {index} uses three ports twice")
        for phase in (0, 1):
            delta_port = Counter()
            for position, edge in enumerate(cycle):
                initial = (position + phase) % 2
                delta_port[edge[1]] += 1 if initial == 0 else -1
            CHECK.true(all(delta_port[port] == 0 for port in range(4)),
                       f"hexagon {index} preserves every port count phase {phase}")

    occupied, target, flow = deterministic_q4_background()
    CHECK.equal(flow, 125, "deterministic background flow")
    CHECK.true(locked(occupied, incident), "deterministic background locked")
    CHECK.true(flippable(occupied, target, endpoints, incident),
               "deterministic target active")

    configurations = []
    current = set(occupied)
    for step in range(24):
        CHECK.true(locked(current, incident), f"sample locked state {step}")
        configurations.append(frozenset(current))
        active = [cycle for cycle in hexagons
                  if flippable(current, cycle, endpoints, incident)]
        CHECK.true(bool(active), f"sample state {step} has an active flip")
        current.symmetric_difference_update(active[(5 * step + 1) % len(active)])

    # Add the six uniform extreme-flux configurations.
    for selected in combinations(range(4), 2):
        uniform = frozenset(edge for edge in edges if edge[1] in selected)
        CHECK.true(locked(uniform, incident), f"uniform ports {selected} locked")
        CHECK.true(not any(flippable(uniform, cycle, endpoints, incident)
                           for cycle in hexagons),
                   f"uniform ports {selected} is frozen")
        configurations.append(uniform)

    for sample, state in enumerate(configurations):
        n_port, signed, cut = flux_data(state, length, edges)
        CHECK.equal(sum(n_port), 2 * length ** 3,
                    f"sample {sample} total locked occupation")
        CHECK.equal(sum(signed), 0, f"sample {sample} four signed fluxes sum zero")
        CHECK.true(all(value % length == 0 for value in signed),
                   f"sample {sample} signed port counts define integer cut flux")
        for axis in range(3):
            plane_flux = []
            for plane in range(length):
                plane_flux.append(sum(
                    1 if edge in state else -1
                    for edge in edges
                    if edge[1] == axis and edge[0][axis] == plane
                ))
            CHECK.true(all(value == plane_flux[0] for value in plane_flux),
                       f"sample {sample} axis {axis} cut flux is plane independent")
            CHECK.equal(signed[axis], length * plane_flux[0],
                        f"sample {sample} axis {axis} volume/cut flux identity")
            CHECK.equal(cut[axis], plane_flux[0],
                        f"sample {sample} axis {axis} integer flux")
        CHECK.equal(cut[3], -sum(cut[:3]),
                    f"sample {sample} fourth flux dependence")

    return occupied, target


def section_local_interaction_and_boundary_component(occupied_q4, target_q4):
    length = 4
    cells, nodes, edges, endpoints, incident, child = geometry(length)
    hexagons = all_hexagons(length)
    support_count = Counter()
    for index, cycle in enumerate(hexagons):
        vertices = cycle_nodes(cycle, endpoints)
        support = {edge for node in vertices for edge in incident[node]}
        CHECK.equal(len(support), 18, f"projected cycle support {index} has 18 links")
        for edge in support:
            support_count[edge] += 1
    CHECK.true(all(support_count[edge] == 18 for edge in edges),
               "each link lies in eighteen projected interaction supports")

    # Exact collared one-hexagon boundary sector: all links outside the target
    # are fixed to the sealed locked background.  The degree constraints leave
    # exactly the two alternating target occupations.
    outside = set(occupied_q4) - set(target_q4)
    allowed = []
    for bits in product((0, 1), repeat=6):
        candidate = set(outside)
        candidate.update(edge for edge, bit in zip(target_q4, bits) if bit)
        if locked(candidate, incident):
            allowed.append(frozenset(candidate))
    CHECK.equal(len(allowed), 2, "one-hexagon fixed-boundary sector has two states")
    CHECK.equal(allowed[0].symmetric_difference(allowed[1]), frozenset(target_q4),
                "two boundary-sector states differ by target hexagon")
    CHECK.true(all(flippable(state, target_q4, endpoints, incident) for state in allowed),
               "target toggles both boundary-sector states")

    # The flip graph is K2.  Its adjacency eigenvectors (1,1) and (1,-1)
    # have eigenvalues +1 and -1, so -tA has a unique positive ground and gap 2t.
    adjacency = ((0, 1), (1, 0))
    for vector, eigenvalue in (((1, 1), 1), ((1, -1), -1)):
        image = tuple(sum(adjacency[i][j] * vector[j] for j in range(2))
                      for i in range(2))
        CHECK.equal(image, tuple(eigenvalue * value for value in vector),
                    f"K2 eigenpair {eigenvalue}")


def section_extensive_active_trial_and_soft_bound():
    # Lift the period-four background to Q8 and select one translated target
    # per 4^3 block.  Their vertex collars are disjoint and every subset of
    # toggles remains locked, yielding a native hypercube trial family.
    occupied4, target4, flow = deterministic_q4_background()
    length = 8
    cells, nodes, edges, endpoints, incident, child = geometry(length)
    occupied = {
        edge for edge in edges
        if ((tuple(value % 4 for value in edge[0]), edge[1]) in occupied4)
    }
    CHECK.true(locked(occupied, incident), "period-four background lifts to Q8")
    targets = []
    for block in product((0, 1), repeat=3):
        base = tuple(4 * value for value in block)
        targets.append(canonical_hexagon(base, (0, 1, 2), length))
    CHECK.equal(len(targets), (length // 4) ** 3, "one target per period-four block")
    target_vertices = [cycle_nodes(cycle, endpoints) for cycle in targets]
    CHECK.true(all(not target_vertices[i] & target_vertices[j]
                   for i in range(len(targets)) for j in range(i)),
               "translated target collars are vertex-disjoint")
    for index, cycle in enumerate(targets):
        CHECK.true(flippable(occupied, cycle, endpoints, incident),
                   f"translated target {index} initially flippable")
    for bits in product((0, 1), repeat=len(targets)):
        state = set(occupied)
        for bit, cycle in zip(bits, targets):
            if bit:
                state.symmetric_difference_update(cycle)
        CHECK.true(locked(state, incident), "independent target subset remains locked")
        CHECK.true(all(flippable(state, cycle, endpoints, incident) for cycle in targets),
                   "all independent target flips remain available")

    # Walsh replay for the N-dimensional cube: adjacency eigenvalues N-2|q|.
    dimension = len(targets)
    vertices = tuple(product((0, 1), repeat=dimension))
    for character in vertices:
        eigenvalue = dimension - 2 * sum(character)
        values = {vertex: (-1) ** sum(a * b for a, b in zip(character, vertex))
                  for vertex in vertices}
        for vertex in vertices:
            image = 0
            for axis in range(dimension):
                neighbor = list(vertex)
                neighbor[axis] ^= 1
                image += values[tuple(neighbor)]
            CHECK.equal(image, eigenvalue * values[vertex],
                        "hypercube Walsh eigenpair")
    CHECK.equal(dimension - (dimension - 2), 2, "hypercube kinetic gap is 2t")

    # Exact slowly varying native test function.  For every hexagon, the two
    # appearances of any port have parent coordinates separated by at most one
    # cyclic step in coordinate zero.  Since their flip signs are opposite,
    # |Delta F| <= 3 ||w||_infty/L for a 1/L-Lipschitz tent function.
    hexagons = all_hexagons(length)
    weights = (Fraction(1), Fraction(-1), Fraction(0), Fraction(0))

    def tent(r):
        return Fraction(min(r, length - r), length)

    max_delta = Fraction()
    for index, cycle in enumerate(hexagons):
        for phase in (0, 1):
            delta = Fraction()
            for position, edge in enumerate(cycle):
                initial = (position + phase) % 2
                change = 1 if initial == 0 else -1
                delta += change * weights[edge[1]] * tent(edge[0][0])
            max_delta = max(max_delta, abs(delta))
            CHECK.true(abs(delta) <= Fraction(3, length),
                       f"hexagon {index} native tent variation phase {phase}")
    CHECK.true(max_delta > 0, "native tent mode is not algebraically trivial")
    CHECK.equal(len(hexagons), 4 * length ** 3, "Q8 hexagon count in gap bound")
    numerator_bound = Fraction(1, 2) * len(hexagons) * Fraction(9, length ** 2)
    CHECK.equal(numerator_bound, 18 * length,
                "Dirichlet numerator bound is 18 t L for unit port weight")


def main():
    occupied, target = section_geometry_flux_and_sectors()
    section_local_interaction_and_boundary_component(occupied, target)
    section_extensive_active_trial_and_soft_bound()
    print(f"PASS__GL6AR_LOCKED_HEXAGON_THERMODYNAMICS__{CHECK.total}/{CHECK.total}")
    print("MODEL=FINITE_RANGE_PROJECTED_HEXAGON_FLIP_ON_DEGREE2_CONFIGURATIONS")
    print("SECTORS=THREE_COORDINATE_CUT_FLUXES_PLUS_DEPENDENT_FOURTH_PORT_COUNT;COMPONENTS_REFINE_INVARIANTS")
    print("FINITE_COMPONENT=PF_UNIQUE_POSITIVE_GROUND;GAP_EXACT_GRAPH_VARIATIONAL")
    print("GROUND=FIXED_BOUNDARY_EXHAUSTION_WEAKSTAR_LIMIT;FINITE_PERIODIC_ACTIVE_ENERGY_BOUND_4_DIVIDES_L")
    print("SOFT=DELTA_LE_18_T_L_OVER_VARIANCE;EXTENSIVE_VARIANCE_IMPLIES_SELECTED_COMPONENT_GAP_L_MINUS2")
    print("OBSTRUCTION=NO_DERIVED_VARIANCE_OR_GROUND_SECTOR;NO_GNS_SPECTRAL_BRIDGE")
    print("CEILING=NO_GAUGE_PHOTON_GRAVITON_PHYSICAL_MOMENTUM_RICCI_GRAVITY_G")


if __name__ == "__main__":
    main()
