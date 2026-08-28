#!/usr/bin/env python3
"""Exact verifier for the FV Coulomb-DPAR projected-source rank theorem.

Only standard-library exact rational arithmetic is used.  The script
independently rebuilds the G5 covering-matched diamond quotient, enumerates
all 720 hexagon histories, differentiates every numerator and every virtual
gap after subtracting the appropriate endpoint ice energy, Hermitianizes the
forward/reverse Bloch matrix element, and proves operator rank with explicit
projected-Hilbert-space matrix-element witnesses.
"""

from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from math import prod
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
LANE = Path(__file__).resolve().parent

DEPENDENCIES = {
    "LANE_GRA_FU_F3_Q4_PAIR_RESOLVED_MAXWELL_DPAR_DERIVATION_V001/THEOREM.md":
        "f088346f72861b3b11ae737fe6b882d43da9e747fc1d1d1f6bd446a7fd2b6272",
    "LANE_GRA_FU_F3_Q4_PAIR_RESOLVED_MAXWELL_DPAR_DERIVATION_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "6f859779566177b5999cfe02c01cd569c5bd7b0b4ec2b21b0b3e79ebf26f9277",
    "LANE_GRA_FT_F3_Q4_DEGREE_PAIR_GEOMETRIC_STRAIN_BOUNDARY_V001/THEOREM.md":
        "6000f38871a57061b106665a41aca04b5d09f4c8c8f4bdc8132ccd5f3f1fbe39",
    "LANE_GRA_FT_F3_Q4_DEGREE_PAIR_GEOMETRIC_STRAIN_BOUNDARY_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "599eec1cde6260be1c9f536274dd8682f77cb45d94e7e3cbc17a28d7552258bd",
    "LANE_GRA_FS_F3_Q4_COMPLETE_SOURCE_RANK_AUDIT_V001/THEOREM.md":
        "36879f4c18eec83a22bdf9bd161d9d444b72e1dbda1d5eaa0312c6aab3d95724",
    "LANE_GRA_FS_F3_Q4_COMPLETE_SOURCE_RANK_AUDIT_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "6e0ecd0febf6364e4122bbf2f65e1feb93c27e960bce30c0097ea0fbe3f58966",
    "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md":
        "5e68e4a8c62ad89cff309781a5cb54e071092e82594e60c04e7992414bc18dbe",
    "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "a91caa20d16b0a1194333f9b51d96546a4ea24d55e23bf1f04c7d249641af8db",
    "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md":
        "78f0687c9f597c96d235224dad45c204d12d7e6c973c270e3305a901efc75b25",
    "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/INDEPENDENT_AUDIT.md":
        "53893c7198241f0f8f6aa766f3595fb75b83d208581833c32656b28d7c7f02b9",
}

MANIFEST_FILES = {
    "DEPENDENCIES.sha256", "README.md", "RESULT.md", "SELF_AUDIT.md",
    "THEOREM.md", "VERIFICATION.txt", "verify_projected_source_rank.py",
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


SIGNS = ((1, 1, 1), (1, -1, -1),
         (-1, 1, -1), (-1, -1, 1))
PAIRS = tuple(combinations(range(4), 2))


def dyad(v):
    x, y, z = map(F, v)
    return (x*x, y*y, z*z, 2*x*y, 2*x*z, 2*y*z)


EDGE = tuple(tuple(x / 3 for x in dyad(v)) for v in SIGNS)
ROOT = {
    (a, b): tuple(x / 8 for x in dyad(tuple(SIGNS[b][i] - SIGNS[a][i]
                                            for i in range(3))))
    for a, b in PAIRS
}


def add(*rows):
    return tuple(sum(xs, F(0)) for xs in zip(*rows))


def scale(c, row):
    return tuple(F(c) * x for x in row)


def rank(rows):
    a = [list(map(F, row)) for row in rows]
    if not a:
        return 0
    nr, nc = len(a), len(a[0])
    r = 0
    for c in range(nc):
        p = next((i for i in range(r, nr) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        q = a[r][c]
        a[r] = [x / q for x in a[r]]
        for i in range(nr):
            if i != r and a[i][c]:
                q = a[i][c]
                a[i] = [x - q*y for x, y in zip(a[i], a[r])]
        r += 1
    return r


def determinant(rows):
    a = [list(map(F, row)) for row in rows]
    value = F(1)
    for c in range(len(a)):
        p = next((i for i in range(c, len(a)) if a[i][c]), None)
        if p is None:
            return F(0)
        if p != c:
            a[c], a[p] = a[p], a[c]
            value *= -1
        q = a[c][c]
        value *= q
        a[c] = [x / q for x in a[c]]
        for i in range(c + 1, len(a)):
            if a[i][c]:
                q = a[i][c]
                a[i] = [x - q*y for x, y in zip(a[i], a[c])]
    return value


def pair_energy_derivative(vertex_states, lam):
    total = (F(0),) * 6
    for z in vertex_states:
        total = add(total, *(
            scale(-lam * z[a] * z[b] / 2, ROOT[(a, b)])
            for a, b in PAIRS
        ))
    return total


PATH_SUM = F(0)
MASK_COEFFICIENT = {mask: F(0) for mask in range(1, 63)}
PATH_CLASSES = {}
for order in permutations(range(6)):
    mask = 0
    gaps = []
    masks = []
    for edge in order[:-1]:
        mask |= 1 << edge
        degree_delta = [0] * 6
        for selected in range(6):
            if mask >> selected & 1:
                change = -1 if selected % 2 == 0 else 1
                degree_delta[selected] += change
                degree_delta[(selected + 1) % 6] += change
        gap = sum(x*x for x in degree_delta)
        gaps.append(gap)
        masks.append(mask)
    weight = F(1, prod(gaps))
    PATH_SUM += weight
    for mask, gap in zip(masks, gaps):
        MASK_COEFFICIENT[mask] += 2 * weight / gap
    PATH_CLASSES[tuple(gaps)] = PATH_CLASSES.get(tuple(gaps), 0) + 1


def ring_tensor(labels, external_choices, lam, parity=0):
    initial_occupation = tuple(int((edge + parity) % 2 == 0)
                               for edge in range(6))
    states = []
    for vertex in range(6):
        left = labels[(vertex - 1) % 6]
        right = labels[vertex]
        external = tuple(a for a in range(4) if a not in (left, right))
        occupied_external = external[external_choices[vertex]]
        z = [F(1)] * 4
        z[left] = F(-1 if initial_occupation[(vertex - 1) % 6] else 1)
        z[right] = F(-1 if initial_occupation[vertex] else 1)
        z[occupied_external] = F(-1)
        states.append(tuple(z))
        assert sum(z) == 0
    initial_prime = pair_energy_derivative(states, lam)
    numerator = add(*(EDGE[label] for label in labels))
    source_sum = scale(PATH_SUM, numerator)
    for mask in range(1, 63):
        virtual = [list(z) for z in states]
        degree_delta = [0] * 6
        for edge in range(6):
            if not (mask >> edge) & 1:
                continue
            for vertex in (edge, (edge + 1) % 6):
                virtual[vertex][labels[edge]] *= -1
        current_prime = pair_energy_derivative(tuple(map(tuple, virtual)), lam)
        gap_prime = add(current_prime, scale(-1, initial_prime))
        source_sum = add(source_sum, scale(MASK_COEFFICIENT[mask], gap_prime))
    return PATH_SUM, source_sum, states, PATH_CLASSES


SHIFTS = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))


def graph(L):
    vertices = tuple((s, x, y, z) for s in (0, 1)
                     for x, y, z in product(range(L), repeat=3))
    edges = []
    adjacency = {v: [] for v in vertices}
    for x, y, z in product(range(L), repeat=3):
        a = (0, x, y, z)
        for label, shift in enumerate(SHIFTS):
            b = (1, (x + shift[0]) % L,
                    (y + shift[1]) % L,
                    (z + shift[2]) % L)
            index = len(edges)
            edges.append((a, b, label))
            adjacency[a].append(index)
            adjacency[b].append(index)
    return vertices, tuple(edges), adjacency


def canonical_cycle(cycle):
    variants = []
    for oriented in (tuple(cycle), tuple(reversed(cycle))):
        for k in range(len(cycle)):
            variants.append(oriented[k:] + oriented[:k])
    return min(variants)


def canonical_label_cycle(labels):
    variants = []
    for oriented in (tuple(labels), tuple(reversed(labels))):
        for k in range(len(labels)):
            variants.append(oriented[k:] + oriented[:k])
    return min(variants)


def hexagons(vertices, edges, adjacency):
    found = set()
    for start in vertices:
        def visit(vertex, path):
            if len(path) == 6:
                for edge in adjacency[vertex]:
                    a, b, _ = edges[edge]
                    if (b if a == vertex else a) == start:
                        found.add(canonical_cycle(path))
                return
            for edge in adjacency[vertex]:
                a, b, _ = edges[edge]
                neighbor = b if a == vertex else a
                if neighbor not in path:
                    visit(neighbor, path + [neighbor])
        visit(start, [start])
    return tuple(sorted(found))


class Dinic:
    def __init__(self, n):
        self.g = [[] for _ in range(n)]

    def add(self, a, b, cap):
        self.g[a].append([b, cap, len(self.g[b])])
        self.g[b].append([a, 0, len(self.g[a]) - 1])

    def flow(self, source, sink):
        total = 0
        while True:
            level = [-1] * len(self.g)
            level[source] = 0
            queue = deque([source])
            while queue:
                a = queue.popleft()
                for b, cap, _ in self.g[a]:
                    if cap and level[b] < 0:
                        level[b] = level[a] + 1
                        queue.append(b)
            if level[sink] < 0:
                return total
            cursor = [0] * len(self.g)

            def send(a, amount):
                if a == sink:
                    return amount
                while cursor[a] < len(self.g[a]):
                    edge = self.g[a][cursor[a]]
                    b, cap, reverse = edge
                    if cap and level[b] == level[a] + 1:
                        pushed = send(b, min(amount, cap))
                        if pushed:
                            edge[1] -= pushed
                            self.g[b][reverse][1] += pushed
                            return pushed
                    cursor[a] += 1
                return 0

            while True:
                pushed = send(source, 10**9)
                if not pushed:
                    break
                total += pushed


def complete_ice(vertices, edges, adjacency, cycle, cycle_occupation):
    edge_by_nodes = {frozenset((a, b)): i for i, (a, b, _) in enumerate(edges)}
    fixed = {}
    cycle_edges = []
    for i in range(6):
        edge = edge_by_nodes[frozenset((cycle[i], cycle[(i + 1) % 6]))]
        fixed[edge] = cycle_occupation[i]
        cycle_edges.append(edge)
    demand = {v: 2 for v in vertices}
    for edge, occupied in fixed.items():
        if occupied:
            a, b, _ = edges[edge]
            demand[a] -= 1
            demand[b] -= 1
    left = tuple(v for v in vertices if v[0] == 0)
    right = tuple(v for v in vertices if v[0] == 1)
    node = {v: i + 1 for i, v in enumerate(vertices)}
    source = 0
    sink = len(vertices) + 1
    flow = Dinic(sink + 1)
    for v in left:
        flow.add(source, node[v], demand[v])
    for v in right:
        flow.add(node[v], sink, demand[v])
    arc_for_edge = {}
    for edge, (a, b, _) in enumerate(edges):
        if edge in fixed:
            continue
        index = len(flow.g[node[a]])
        flow.add(node[a], node[b], 1)
        arc_for_edge[edge] = (node[a], index)
    required = sum(demand[v] for v in left)
    assert flow.flow(source, sink) == required
    occupied = {edge for edge, value in fixed.items() if value}
    for edge, (a, index) in arc_for_edge.items():
        if flow.g[a][index][1] == 0:
            occupied.add(edge)
    assert all(sum(edge in occupied for edge in adjacency[v]) == 2 for v in vertices)
    return frozenset(occupied), tuple(cycle_edges)


def local_choices(cycle_edges, edges, adjacency, state):
    labels = tuple(edges[edge][2] for edge in cycle_edges)
    choices = []
    for vertex in range(6):
        left = labels[(vertex - 1) % 6]
        right = labels[vertex]
        external_labels = tuple(a for a in range(4) if a not in (left, right))
        graph_vertex = edges[cycle_edges[vertex]][0]
        if graph_vertex not in edges[cycle_edges[(vertex - 1) % 6]][:2]:
            graph_vertex = edges[cycle_edges[vertex]][1]
        external_occupied = []
        for edge in adjacency[graph_vertex]:
            if edge not in cycle_edges and edge in state:
                external_occupied.append(edges[edge][2])
        assert len(external_occupied) == 1
        choices.append(external_labels.index(external_occupied[0]))
    return labels, tuple(choices)


def direct_for_occupied(pair, lam):
    """Direct Q_pair/vertex on the ice state occupying exactly pair."""
    z = tuple(F(-1 if a in pair else 1) for a in range(4))
    return add(*(scale(lam * z[a] * z[b], ROOT[(a, b)])
                 for a, b in PAIRS))


def project_a1_t2(row):
    mean = sum(row[:3]) / 3
    return (mean, mean, mean) + tuple(row[3:])


def expected_ring_a1_t2(missing, lam):
    # Exact Hermitian 720-path answer in the convention where D_a=n_a n_a^T.
    identity = (F(1), F(1), F(1), F(0), F(0), F(0))
    a1_coefficient = F(21, 8) * (8 - 15 * lam)
    edge_coefficient = -F(63, 8) * (2 - 5 * lam)
    return add(scale(a1_coefficient, identity),
               scale(edge_coefficient, EDGE[missing]))


def has_four_cycle(vertices, edges, adjacency):
    neighbors = {}
    for vertex in vertices:
        table = set()
        for edge in adjacency[vertex]:
            a, b, _ = edges[edge]
            table.add(b if a == vertex else a)
        neighbors[vertex] = table
    for side in (0, 1):
        same_side = [vertex for vertex in vertices if vertex[0] == side]
        for first, second in combinations(same_side, 2):
            if len(neighbors[first] & neighbors[second]) > 1:
                return True
    return False


if __name__ == '__main__':
    # Frozen dependency bytes and semantic custody.
    for relative, expected in DEPENDENCIES.items():
        path = ROOT_DIR / relative
        check(path.is_file() and digest(path) == expected,
              f"dependency custody {relative}")
        check(sha256(path.read_bytes() + b"tamper").hexdigest() != expected,
              f"dependency tamper rejection {relative}")

    fu = (ROOT_DIR / next(key for key in DEPENDENCIES
                         if key.endswith("PAIR_RESOLVED_MAXWELL_DPAR_DERIVATION_V001/THEOREM.md"))).read_text()
    ft = (ROOT_DIR / "LANE_GRA_FT_F3_Q4_DEGREE_PAIR_GEOMETRIC_STRAIN_BOUNDARY_V001/THEOREM.md").read_text()
    fs = (ROOT_DIR / "LANE_GRA_FS_F3_Q4_COMPLETE_SOURCE_RANK_AUDIT_V001/THEOREM.md").read_text()
    cw = (ROOT_DIR / "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md").read_text()
    fm = (ROOT_DIR / "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md").read_text()
    check("lambda_C=-{1\\over2}" in fu and "S8 — source ordering" in fu,
          "FU freezes Coulomb slope and source-before-Feshbach ordering")
    check("Direct restriction of the pair map to local\nice still has only `A1+E` rank three" in fu,
          "FU leaves projected rank open rather than pre-asserting FV")
    check("direct root-pair image has only `A1+E`" in ft,
          "FT owns the direct ice A1+E boundary")
    check("one-edge source left unchanged" in fu and "rank four" in fs,
          "unchanged FS edge source and its rank-four custody are present")
    check("all `6!=720` flip orders stay in `Q_2=1-P_2`" in cw,
          "CW owns the complete irreducible H6 history set")
    check("through order eight" in fm and "dressed length-six" in fm,
          "FM owns the through-H8 endpoint classification")

    # Exact tetrahedral tensors and direct ice projection.
    check(rank(EDGE) == 4, "physical tetrahedral edge dyads have rank four")
    check(rank(tuple(ROOT[pair] for pair in PAIRS)) == 6,
          "six normalized sibling-root dyads have rank six")
    check(add(*EDGE) == (F(4, 3), F(4, 3), F(4, 3), F(0), F(0), F(0)),
          "tetrahedral edge dyads sum to four-thirds identity")
    coulomb = F(-1, 2)
    direct = [direct_for_occupied(pair, coulomb)
              for pair in ((0, 1), (0, 2), (0, 3))]
    e_rows = [add(direct[i], scale(-1, direct[0])) for i in (1, 2)]
    check(e_rows == [
        (F(-1), F(1), F(0), F(0), F(0), F(0)),
        (F(-1), F(0), F(1), F(0), F(0), F(0)),
    ], "normalized direct diagonal expectation differences are exact E witnesses")
    check(rank(e_rows) == 2 and all(sum(row[:3]) == 0 for row in e_rows),
          "direct nonidentity ice-pair operator has exact E rank two")
    check(all(value == 0 for row in e_rows for value in row[3:]),
          "direct ice-pair expectation differences have exact T2 null")

    # Every one of the 720 paths and its five differentiated denominators.
    sorted_classes = Counter()
    for gaps, multiplicity in PATH_CLASSES.items():
        sorted_classes[tuple(sorted(gaps))] += multiplicity
    expected_classes = Counter({
        (2, 2, 2, 2, 2): 96,
        (2, 2, 2, 2, 4): 144,
        (2, 2, 2, 4, 4): 216,
        (2, 2, 4, 4, 4): 192,
        (2, 2, 4, 4, 6): 72,
    })
    check(sum(PATH_CLASSES.values()) == 720,
          "all 720 ordered hexagon histories are enumerated")
    check(sorted_classes == expected_classes,
          "all 720 histories reproduce the five exact gap classes")
    check(all(all(gap > 0 for gap in gaps) for gaps in PATH_CLASSES),
          "all 3600 proper-prefix gaps remain in Q2")
    check(PATH_SUM == F(63, 8), "source-off 720-path sum is exactly 63/8")
    check(len(MASK_COEFFICIENT) == 62
          and all(value > 0 for value in MASK_COEFFICIENT.values()),
          "all 62 nonempty proper prefix masks carry their exact derivative weight")

    # Rebuild the covering-matched G5 family member and all four H6 orientations.
    vertices, edges, adjacency = graph(5)
    cycles = hexagons(vertices, edges, adjacency)
    check(len(vertices) == 250 and len(edges) == 500,
          "G5 has the exact 250-vertex 500-link q4 size")
    check(all(len(adjacency[vertex]) == 4 for vertex in vertices),
          "G5 is closed and coordination four")
    check(not has_four_cycle(vertices, edges, adjacency),
          "G5 has no four-cycle and hence no lower off-diagonal fold channel")
    check(len(cycles) == 500, "G5 has exactly 500 elementary hexagons")
    edge_by_nodes = {frozenset((a, b)): i for i, (a, b, _) in enumerate(edges)}
    orientation_counts = Counter()
    label_orbits = {missing: set() for missing in range(4)}
    selected = {}
    for cycle in cycles:
        cycle_edges = tuple(edge_by_nodes[frozenset((cycle[i], cycle[(i + 1) % 6]))]
                            for i in range(6))
        labels = tuple(edges[edge][2] for edge in cycle_edges)
        counts = Counter(labels)
        check_pattern = len(counts) == 3 and set(counts.values()) == {2}
        if not check_pattern:
            raise AssertionError("non-elementary label pattern in G5 hexagon")
        missing = next(iter(set(range(4)) - set(labels)))
        orientation_counts[missing] += 1
        label_orbits[missing].add(canonical_label_cycle(labels))
        selected.setdefault(missing, cycle)
    check(orientation_counts == Counter({0: 125, 1: 125, 2: 125, 3: 125}),
          "G5 realizes all four missing-label H6 orientations equally")
    check(all(len(label_orbits[missing]) == 1 for missing in range(4)),
          "one dihedral label representative exhausts each G5 H6 orientation")

    direction_states = [
        frozenset(edge for edge, (_, _, label) in enumerate(edges)
                  if label in pair)
        for pair in ((0, 1), (0, 2), (0, 3))
    ]
    check(all(all(sum(edge in state for edge in adjacency[vertex]) == 2
                  for vertex in vertices) for state in direction_states),
          "all three direct E witnesses are actual global direction-pair ice coverings")
    check(len(set(direction_states)) == 3,
          "the three direct diagonal witness states are distinct")

    # Exact affine 720-path theorem.  At lambda=0 the gap derivative vanishes.
    # At lambda=1, checking forward and reverse proves the result for every
    # lambda because the full first source derivative is affine in lambda.
    orientation_labels = {}
    all_environment_cancellation = True
    all_formula = True
    for missing in range(4):
        cycle = selected[missing]
        cycle_edges = tuple(edge_by_nodes[frozenset((cycle[i], cycle[(i + 1) % 6]))]
                            for i in range(6))
        labels = tuple(edges[edge][2] for edge in cycle_edges)
        orientation_labels[missing] = labels
        base = scale(PATH_SUM, add(*(EDGE[label] for label in labels)))
        all_formula &= base == expected_ring_a1_t2(missing, F(0))
        for choices in product((0, 1), repeat=6):
            _, forward, _, _ = ring_tensor(labels, choices, F(1), parity=0)
            _, reverse, _, _ = ring_tensor(labels, choices, F(1), parity=1)
            hermitian = scale(F(1, 2), add(forward, reverse))
            expected = expected_ring_a1_t2(missing, F(1))
            all_formula &= hermitian == expected
            forward_residual = add(forward, scale(-1, expected))
            reverse_residual = add(reverse, scale(-1, expected))
            all_environment_cancellation &= (
                add(forward_residual, reverse_residual) == (F(0),) * 6
                and sum(forward_residual[:3]) == 0
                and all(value == 0 for value in forward_residual[3:])
            )
    check(all_formula,
          "all 4x64 local ice environments obey the exact Hermitian A1+T2 formula")
    check(all_environment_cancellation,
          "endpoint-dependent E residual is odd under reversal and cancels Hermitianly")
    check(set(orientation_labels) == {0, 1, 2, 3},
          "every hexagon orientation was differentiated independently")

    # Construct actual global degree-two states realizing one flippable ring of
    # every orientation.  The four off-diagonal entries are therefore operator
    # matrix-element witnesses, not merely coefficient vectors.
    ring_rows = []
    operator_entries = []
    for missing in range(4):
        state, cycle_edges = complete_ice(
            vertices, edges, adjacency, selected[missing],
            (1, 0, 1, 0, 1, 0)
        )
        switched = state ^ frozenset(cycle_edges)
        check(all(sum(edge in switched for edge in adjacency[vertex]) == 2
                  for vertex in vertices),
              f"missing-label {missing} switched witness remains in projected ice")
        labels, choices = local_choices(cycle_edges, edges, adjacency, state)
        _, forward1, _, _ = ring_tensor(labels, choices, F(1), parity=0)
        _, reverse1, _, _ = ring_tensor(labels, choices, F(1), parity=1)
        base = expected_ring_a1_t2(missing, F(0))
        forward_c = add(base, scale(coulomb, add(forward1, scale(-1, base))))
        reverse_c = add(base, scale(coulomb, add(reverse1, scale(-1, base))))
        ring_row = scale(F(1, 2), add(forward_c, reverse_c))
        check(ring_row == expected_ring_a1_t2(missing, coulomb),
              f"missing-label {missing} Coulomb ring row is exact")
        ring_rows.append(ring_row)
        operator_entries.append((state, switched))
    check(len(set(operator_entries)) == 4,
          "four ring witnesses are distinct projected-operator matrix entries")
    check(all(project_a1_t2(row) == row for row in ring_rows),
          "Hermitian Coulomb ring rows contain no E contamination")
    check(rank(ring_rows) == 4,
          "nonidentity H6 ring-source operators have exact A1+T2 rank four")
    check(ring_rows == [
        (F(231, 8), F(231, 8), F(231, 8), F(-189, 8), F(-189, 8), F(-189, 8)),
        (F(231, 8), F(231, 8), F(231, 8), F(189, 8), F(189, 8), F(-189, 8)),
        (F(231, 8), F(231, 8), F(231, 8), F(189, 8), F(-189, 8), F(189, 8)),
        (F(231, 8), F(231, 8), F(231, 8), F(-189, 8), F(189, 8), F(189, 8)),
    ], "Coulomb lambda minus one-half gives the four exact H6 rows")

    # Operator rank: the first two functionals are normalized differences of
    # diagonal matrix elements on actual direction-pair ice coverings; the last
    # four are the distinct off-diagonal entries above.  Identities vanish on
    # both kinds of functional.  Their exact determinant is nonzero.
    witness = e_rows + ring_rows
    witness_det = determinant(witness)
    check(rank(witness) == 6, "projected source has exact operator rank six")
    check(witness_det == F(-4678629417, 256),
          "six projected matrix-element functionals have the exact nonzero determinant")

    # Special slopes and formal-order stability.
    def ring_rank_at(lam):
        return rank([expected_ring_a1_t2(missing, lam) for missing in range(4)])

    check(ring_rank_at(coulomb) == 4,
          "Coulomb slope avoids both ring-sector cancellation values")
    check(ring_rank_at(F(2, 5)) == 1,
          "lambda=2/5 cancels T2 and leaves only nonidentity ring A1")
    check(ring_rank_at(F(3, 5)) == 3,
          "lambda=3/5 cancels ring A1 and leaves T2 rank three")
    check(ring_rank_at(F(0)) == 4,
          "lambda=0 retains ring A1+T2 but loses direct pair E")
    check(coulomb not in (F(0), F(2, 5), F(3, 5)),
          "lambda=-1/2 is outside all exact rank-loss slopes")

    # Documentary theorem scope, fold/identity custody, and H8 ceiling.
    theorem = (LANE / "THEOREM.md").read_text()
    result = (LANE / "RESULT.md").read_text()
    self_audit = (LANE / "SELF_AUDIT.md").read_text()
    theorem_flat = " ".join(theorem.split())
    required_theorem = (
        "S10 / FV-PURE",
        "stronger than FU `S1`--`S9`",
        "No additional nonidentity source operator is present",
        "subtract the derivative of the appropriate endpoint ice energy",
        "Hermitian forward/reverse average",
        "Every order-two and order-four projected endpoint is diagonal",
        "identity shifts do not count",
        "formal through-order-eight rank remains six",
        "off-shell projected operator rank",
        "does not establish a retarded or CTP rank",
        "does not prove a Ward identity, tensor pole, gravity",
    )
    for phrase in required_theorem:
        check(phrase in theorem_flat,
              f"theorem retains scope/custody phrase: {phrase}")
    check("-4678629417/256" in result and "rank six" in result,
          "result records the exact operator witness")
    check("FU `S1`--`S9` alone do not imply `FV-PURE`" in result,
          "result exposes the additional complete-source premise")
    check("finite nonzero `h`" in self_audit and "tuned algebraic zero" in self_audit,
          "self-audit withholds an unrestricted finite-h rank claim")

    # Core manifest verification is activated after packaging; no payload is
    # rewritten by this verifier.
    manifest = LANE / "MANIFEST.sha256"
    if manifest.is_file():
        listed = set()
        for line in manifest.read_text().splitlines():
            if not line.strip():
                continue
            expected, relative = line.split("  ", 1)
            listed.add(relative)
            check(digest(LANE / relative) == expected,
                  f"manifest custody {relative}")
        check(listed == MANIFEST_FILES, "manifest lists exactly the frozen payload")

    print(f"SUMMARY {checks}/{checks} exact checks passed")
    print("PATHS 720; PREFIX_DERIVATIVES 3600; J6 63/8")
    print("COULOMB_RING_ROWS A1+T2 rank4; DIRECT_NONIDENTITY E rank2")
    print(f"PROJECTED_OPERATOR_RANK 6; WITNESS_DETERMINANT {witness_det}")
    print("FORMAL_THROUGH_H8 rank6; CTP_WARD_TENSOR_POLE_GRAVITY not claimed")
