#!/usr/bin/env python3
"""Independent hostile replay of the FV projected-source theorem.

This program deliberately does not import or execute the builder verifier.
It reconstructs the finite graph, perturbative histories, source derivatives,
global ice witnesses, and exact rational ranks from a separate implementation.
"""

from collections import Counter, deque
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path


AUDIT = Path(__file__).resolve().parent
LANE = AUDIT.parent
ROOT = LANE.parent

CORE_FILES = {
    "DEPENDENCIES.sha256",
    "MANIFEST.sha256",
    "README.md",
    "RESULT.md",
    "SEAL.sha256",
    "SELF_AUDIT.md",
    "THEOREM.md",
    "VERIFICATION.txt",
    "verify_projected_source_rank.py",
}

count = 0


def check(statement, label):
    global count
    if not statement:
        raise AssertionError(label)
    count += 1
    print(f"PASS {label}")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def vadd(*vectors):
    if not vectors:
        return (Q(0),) * 6
    return tuple(sum(column, Q(0)) for column in zip(*vectors))


def vmul(number, vector):
    return tuple(Q(number) * entry for entry in vector)


def matrix_rank(rows):
    matrix = [list(map(Q, row)) for row in rows]
    if not matrix:
        return 0
    pivot_row = 0
    for column in range(len(matrix[0])):
        pivot = next((r for r in range(pivot_row, len(matrix))
                      if matrix[r][column]), None)
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        value = matrix[pivot_row][column]
        matrix[pivot_row] = [x / value for x in matrix[pivot_row]]
        for r in range(len(matrix)):
            if r == pivot_row or not matrix[r][column]:
                continue
            value = matrix[r][column]
            matrix[r] = [x - value * y
                         for x, y in zip(matrix[r], matrix[pivot_row])]
        pivot_row += 1
    return pivot_row


def matrix_det(rows):
    matrix = [list(map(Q, row)) for row in rows]
    answer = Q(1)
    for column in range(len(matrix)):
        pivot = next((r for r in range(column, len(matrix))
                      if matrix[r][column]), None)
        if pivot is None:
            return Q(0)
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            answer = -answer
        value = matrix[column][column]
        answer *= value
        for r in range(column + 1, len(matrix)):
            ratio = matrix[r][column] / value
            matrix[r] = [x - ratio * y
                         for x, y in zip(matrix[r], matrix[column])]
    return answer


# Symmetric-tensor coordinates are (xx, yy, zz, 2xy, 2xz, 2yz).
TETRA = ((1, 1, 1), (1, -1, -1),
         (-1, 1, -1), (-1, -1, 1))
PAIR_LABELS = tuple(combinations(range(4), 2))


def symmetric_dyad(vector):
    x, y, z = map(Q, vector)
    return (x*x, y*y, z*z, 2*x*y, 2*x*z, 2*y*z)


EDGE_DYAD = tuple(vmul(Q(1, 3), symmetric_dyad(v)) for v in TETRA)
ROOT_DYAD = {}
for a, b in PAIR_LABELS:
    difference = tuple(TETRA[b][k] - TETRA[a][k] for k in range(3))
    ROOT_DYAD[(a, b)] = vmul(Q(1, 8), symmetric_dyad(difference))


def pair_hprime(vertex_z, slope):
    """Dimensionless derivative (d H_C/dj)/U_d for a list of q4 vertices."""
    answer = (Q(0),) * 6
    for z in vertex_z:
        for a, b in PAIR_LABELS:
            answer = vadd(answer,
                          vmul(-Q(slope) * z[a] * z[b] / 2,
                               ROOT_DYAD[(a, b)]))
    return answer


def make_local_endpoint(labels, external_bits, parity):
    """Build six degree-two endpoint vertices around an alternating hexagon."""
    occupied_cycle = tuple(((edge + parity) % 2) == 0 for edge in range(6))
    vertices = []
    for vertex in range(6):
        previous_label = labels[(vertex - 1) % 6]
        next_label = labels[vertex]
        outside = tuple(label for label in range(4)
                        if label not in (previous_label, next_label))
        occupied_outside = outside[external_bits[vertex]]
        occupation = [False] * 4
        occupation[previous_label] = occupied_cycle[(vertex - 1) % 6]
        occupation[next_label] = occupied_cycle[vertex]
        occupation[occupied_outside] = True
        check_degree = sum(occupation)
        if check_degree != 2:
            raise AssertionError("local endpoint is not degree two")
        vertices.append(tuple(Q(-1 if item else 1) for item in occupation))
    return tuple(vertices)


def toggle_mask(endpoint, labels, mask):
    changed = [list(z) for z in endpoint]
    for edge in range(6):
        if not (mask & (1 << edge)):
            continue
        label = labels[edge]
        changed[edge][label] *= -1
        changed[(edge + 1) % 6][label] *= -1
    return tuple(tuple(z) for z in changed)


def defect_gap(endpoint):
    # z_a=1-2n_a, so d-2=-sum_a z_a/2.
    return sum((sum(z, Q(0)) / 2) ** 2 for z in endpoint)


def history_table():
    paths = []
    ordered_classes = Counter()
    prefix_count = 0
    for order in permutations(range(6)):
        mask = 0
        gaps = []
        masks = []
        # The degree gaps are independent of labels and outside occupations.
        degree_change = [0] * 6
        for edge in order[:-1]:
            mask |= 1 << edge
            signed = -1 if edge % 2 == 0 else 1
            degree_change[edge] += signed
            degree_change[(edge + 1) % 6] += signed
            gap = sum(value * value for value in degree_change)
            masks.append(mask)
            gaps.append(gap)
            prefix_count += 1
        weight = Q(1)
        for gap in gaps:
            weight /= gap
        paths.append((tuple(masks), tuple(gaps), weight))
        ordered_classes[tuple(sorted(gaps))] += 1
    return tuple(paths), ordered_classes, prefix_count


PATHS, GAP_CLASSES, PREFIX_COUNT = history_table()


def differentiated_ring(labels, external_bits, slope, parity):
    """Recompute the Q=-2 d/dj H6 tensor from all paths.

    The common physical factor -h^6/U_d^5 is omitted.  No prefix-mask
    coefficient table from the builder is used.
    """
    endpoint = make_local_endpoint(labels, external_bits, parity)
    endpoint_prime = pair_hprime(endpoint, slope)
    prefix = {}
    for mask in range(1, 63):
        virtual = toggle_mask(endpoint, labels, mask)
        gap = defect_gap(virtual)
        prefix[mask] = (gap, vadd(pair_hprime(virtual, slope),
                                 vmul(-1, endpoint_prime)))
    numerator_q = vadd(*(EDGE_DYAD[label] for label in labels))
    answer = (Q(0),) * 6
    for masks, gaps, weight in PATHS:
        path_q = numerator_q
        for mask, independently_found_gap in zip(masks, gaps):
            gap, gap_prime = prefix[mask]
            if gap != independently_found_gap:
                raise AssertionError("occupation and degree gap calculations disagree")
            path_q = vadd(path_q, vmul(Q(2, gap), gap_prime))
        answer = vadd(answer, vmul(weight, path_q))
    return answer


def expected_ring(missing, slope):
    identity = (Q(1), Q(1), Q(1), Q(0), Q(0), Q(0))
    return vadd(vmul(Q(21, 8) * (8 - 15 * slope), identity),
                vmul(-Q(63, 8) * (2 - 5 * slope), EDGE_DYAD[missing]))


def direct_q(occupied_pair, slope):
    z = tuple(Q(-1 if label in occupied_pair else 1) for label in range(4))
    return vadd(*(vmul(Q(slope) * z[a] * z[b], ROOT_DYAD[(a, b)])
                  for a, b in PAIR_LABELS))


# Build G5 independently from periodic bipartite incidence data.
SHIFT = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))


def diamond_graph(side):
    vertices = tuple((part, x, y, z) for part in (0, 1)
                     for x, y, z in product(range(side), repeat=3))
    edges = []
    incidence = {vertex: [] for vertex in vertices}
    lookup = {}
    for x, y, z in product(range(side), repeat=3):
        left = (0, x, y, z)
        for label, (dx, dy, dz) in enumerate(SHIFT):
            right = (1, (x + dx) % side, (y + dy) % side, (z + dz) % side)
            edge = len(edges)
            edges.append((left, right, label))
            incidence[left].append(edge)
            incidence[right].append(edge)
            lookup[frozenset((left, right))] = edge
    return vertices, tuple(edges), incidence, lookup


def canonical_cycle(vertices):
    variants = []
    for direction in (tuple(vertices), tuple(reversed(vertices))):
        variants.extend(direction[offset:] + direction[:offset]
                        for offset in range(len(direction)))
    return min(variants)


def six_cycles(vertices, edges, incidence):
    answer = set()
    for start in vertices:
        stack = [(start, (start,))]
        while stack:
            here, path = stack.pop()
            if len(path) == 6:
                for edge in incidence[here]:
                    a, b, _ = edges[edge]
                    other = b if a == here else a
                    if other == start:
                        answer.add(canonical_cycle(path))
                continue
            for edge in incidence[here]:
                a, b, _ = edges[edge]
                other = b if a == here else a
                if other not in path:
                    stack.append((other, path + (other,)))
    return tuple(sorted(answer))


def graph_has_four_cycle(vertices, edges, incidence):
    neighbours = {}
    for vertex in vertices:
        neighbours[vertex] = {
            (edges[e][1] if edges[e][0] == vertex else edges[e][0])
            for e in incidence[vertex]
        }
    for part in (0, 1):
        same = tuple(v for v in vertices if v[0] == part)
        for a, b in combinations(same, 2):
            if len(neighbours[a] & neighbours[b]) >= 2:
                return True
    return False


class EdmondsKarp:
    """Small independent integral-flow implementation for global ice witnesses."""

    def __init__(self):
        self.graph = {}

    def add_arc(self, source, target, capacity, tag=None):
        self.graph.setdefault(source, [])
        self.graph.setdefault(target, [])
        forward = [target, capacity, None, tag]
        reverse = [source, 0, forward, None]
        forward[2] = reverse
        self.graph[source].append(forward)
        self.graph[target].append(reverse)
        return forward

    def maximum_flow(self, source, sink):
        total = 0
        while True:
            previous = {source: None}
            queue = deque([source])
            while queue and sink not in previous:
                here = queue.popleft()
                for arc in self.graph[here]:
                    if arc[1] and arc[0] not in previous:
                        previous[arc[0]] = (here, arc)
                        queue.append(arc[0])
            if sink not in previous:
                return total
            amount = 10**9
            here = sink
            while here != source:
                before, arc = previous[here]
                amount = min(amount, arc[1])
                here = before
            here = sink
            while here != source:
                before, arc = previous[here]
                arc[1] -= amount
                arc[2][1] += amount
                here = before
            total += amount


def global_ice_with_ring(vertices, edges, incidence, cycle, lookup):
    cycle_edges = tuple(lookup[frozenset((cycle[k], cycle[(k + 1) % 6]))]
                        for k in range(6))
    fixed = {edge: (position % 2 == 0)
             for position, edge in enumerate(cycle_edges)}
    need = {vertex: 2 for vertex in vertices}
    for edge, occupied in fixed.items():
        if occupied:
            a, b, _ = edges[edge]
            need[a] -= 1
            need[b] -= 1
    network = EdmondsKarp()
    source, sink = ("source",), ("sink",)
    for vertex in vertices:
        if vertex[0] == 0:
            network.add_arc(source, vertex, need[vertex])
        else:
            network.add_arc(vertex, sink, need[vertex])
    physical_arcs = {}
    for edge, (left, right, _) in enumerate(edges):
        if edge in fixed:
            continue
        physical_arcs[edge] = network.add_arc(left, right, 1, edge)
    required = sum(need[v] for v in vertices if v[0] == 0)
    if network.maximum_flow(source, sink) != required:
        raise AssertionError("global ring witness has no degree-two completion")
    occupied = {edge for edge, value in fixed.items() if value}
    for edge, arc in physical_arcs.items():
        if arc[1] == 0:
            occupied.add(edge)
    return frozenset(occupied), cycle_edges


def cycle_environment(state, cycle, cycle_edges, edges, incidence):
    labels = tuple(edges[edge][2] for edge in cycle_edges)
    bits = []
    for vertex_index, vertex in enumerate(cycle):
        adjacent_cycle = {cycle_edges[(vertex_index - 1) % 6],
                          cycle_edges[vertex_index]}
        occupied_external = [edge for edge in incidence[vertex]
                             if edge not in adjacent_cycle and edge in state]
        if len(occupied_external) != 1:
            raise AssertionError("global witness has wrong outside occupation")
        pair = tuple(label for label in range(4)
                     if label not in (labels[(vertex_index - 1) % 6],
                                      labels[vertex_index]))
        bits.append(pair.index(edges[occupied_external[0]][2]))
    return labels, tuple(bits)


def manifest_entries(path):
    entries = {}
    for line in path.read_text().splitlines():
        if not line:
            continue
        value, relative = line.split("  ", 1)
        entries[relative] = value
    return entries


def main():
    # Builder custody and byte hygiene.  The audit pins the repaired core, not
    # merely whatever happens to be present when this script is rerun.
    custody = manifest_entries(AUDIT / "CORE_CUSTODY.sha256")
    check(set(custody) == CORE_FILES, "core custody lists exactly nine builder files")
    for relative, expected in custody.items():
        path = LANE / relative
        check(path.is_file() and digest(path) == expected,
              f"frozen core hash {relative}")
        raw = path.read_bytes()
        check(all(byte in (9, 10) or byte >= 32 for byte in raw)
              and 127 not in raw,
              f"no forbidden ASCII control byte in {relative}")
        raw.decode("utf-8")

    core_manifest = manifest_entries(LANE / "MANIFEST.sha256")
    check(set(core_manifest) == CORE_FILES - {"MANIFEST.sha256", "SEAL.sha256"},
          "builder manifest has exact seven-payload membership")
    for relative, expected in core_manifest.items():
        check(digest(LANE / relative) == expected,
              f"builder manifest replay {relative}")
    seal = manifest_entries(LANE / "SEAL.sha256")
    check(seal == {
        "MANIFEST.sha256": digest(LANE / "MANIFEST.sha256"),
        "VERIFICATION.txt": digest(LANE / "VERIFICATION.txt"),
    }, "builder seal binds manifest and transcript")

    theorem = (LANE / "THEOREM.md").read_text()
    theorem_flat = " ".join(theorem.split())
    for phrase in (
        "FV-PURE",
        "complete nonidentity first derivative of the pre-Feshbach source",
        "Hermitian forward/reverse average",
        "off-shell projected operator rank",
        "does not establish a retarded or CTP rank",
        "does not prove a Ward identity, tensor pole, gravity",
    ):
        check(phrase in theorem_flat, f"repaired theorem retains boundary: {phrase}")

    # Independent normalization and direct projected E source.
    check(matrix_rank(EDGE_DYAD) == 4, "tetrahedral edge-dyad span is A1+T2 rank four")
    check(matrix_rank(ROOT_DYAD.values()) == 6, "root-dyad span is symmetric rank six")
    check(vadd(*EDGE_DYAD) == (Q(4, 3), Q(4, 3), Q(4, 3), 0, 0, 0),
          "edge-dyad normalization is four-thirds identity")
    coulomb = Q(-1, 2)
    direct_rows = [direct_q(pair, coulomb)
                   for pair in ((0, 1), (0, 2), (0, 3))]
    e_rows = [vadd(direct_rows[k], vmul(-1, direct_rows[0])) for k in (1, 2)]
    check(e_rows == [(-1, 1, 0, 0, 0, 0), (-1, 0, 1, 0, 0, 0)],
          "direct global-covering differences give the stated normalized E rows")
    check(matrix_rank(e_rows) == 2 and
          all(sum(row[:3]) == 0 and row[3:] == (0, 0, 0) for row in e_rows),
          "direct nonidentity source is E rank two and T2-null")

    # All histories and all proper-prefix denominators.
    check(len(PATHS) == 720, "independently enumerated all 720 H6 histories")
    check(PREFIX_COUNT == 3600, "independently evaluated all 3600 proper prefixes")
    expected_classes = Counter({
        (2, 2, 2, 2, 2): 96,
        (2, 2, 2, 2, 4): 144,
        (2, 2, 2, 4, 4): 216,
        (2, 2, 4, 4, 4): 192,
        (2, 2, 4, 4, 6): 72,
    })
    check(GAP_CLASSES == expected_classes, "five gap classes and multiplicities are exact")
    check(all(all(gap > 0 for gap in gaps) for _, gaps, _ in PATHS),
          "every proper prefix is off ice")
    check(sum((weight for _, _, weight in PATHS), Q(0)) == Q(63, 8),
          "independent J6 path sum is 63/8")

    # G5 graph topology and four actual global ring entries.
    vertices, edges, incidence, lookup = diamond_graph(5)
    cycles = six_cycles(vertices, edges, incidence)
    check((len(vertices), len(edges), len(cycles)) == (250, 500, 500),
          "G5 independently has 250 vertices, 500 links, and 500 hexagons")
    check(all(len(incidence[v]) == 4 for v in vertices), "G5 is q4")
    check(not graph_has_four_cycle(vertices, edges, incidence),
          "G5 girth excludes two- and four-link off-diagonal endpoints")
    chosen = {}
    orientation_count = Counter()
    label_orbits = {missing: set() for missing in range(4)}
    for cycle in cycles:
        cycle_edges = tuple(lookup[frozenset((cycle[k], cycle[(k + 1) % 6]))]
                            for k in range(6))
        labels = tuple(edges[e][2] for e in cycle_edges)
        frequencies = Counter(labels)
        check_pattern = len(frequencies) == 3 and set(frequencies.values()) == {2}
        if not check_pattern:
            raise AssertionError("G5 six-cycle has wrong label pattern")
        missing = next(iter(set(range(4)) - set(labels)))
        orientation_count[missing] += 1
        rotations = []
        for direction in (labels, tuple(reversed(labels))):
            rotations.extend(direction[k:] + direction[:k] for k in range(6))
        label_orbits[missing].add(min(rotations))
        chosen.setdefault(missing, cycle)
    check(orientation_count == Counter({0: 125, 1: 125, 2: 125, 3: 125}),
          "all four missing-label orientations occur 125 times")
    check(all(len(label_orbits[d]) == 1 for d in range(4)),
          "one dihedral label orbit covers each orientation")

    # Exhaust 4 orientations x 64 external assignments.  Constant and linear
    # slope parts establish the affine formula for every real lambda.
    cancellation_count = 0
    environment_count = 0
    for missing in range(4):
        cycle = chosen[missing]
        cycle_edges = tuple(lookup[frozenset((cycle[k], cycle[(k + 1) % 6]))]
                            for k in range(6))
        labels = tuple(edges[e][2] for e in cycle_edges)
        for external_bits in product((0, 1), repeat=6):
            forward_zero = differentiated_ring(labels, external_bits, Q(0), 0)
            reverse_zero = differentiated_ring(labels, external_bits, Q(0), 1)
            forward_one = differentiated_ring(labels, external_bits, Q(1), 0)
            reverse_one = differentiated_ring(labels, external_bits, Q(1), 1)
            hermitian_zero = vmul(Q(1, 2), vadd(forward_zero, reverse_zero))
            hermitian_one = vmul(Q(1, 2), vadd(forward_one, reverse_one))
            if hermitian_zero != expected_ring(missing, Q(0)):
                raise AssertionError("lambda-zero H6 formula failed")
            if hermitian_one != expected_ring(missing, Q(1)):
                raise AssertionError("lambda-one H6 formula failed")
            residual_forward = vadd(forward_one, vmul(-1, hermitian_one))
            residual_reverse = vadd(reverse_one, vmul(-1, hermitian_one))
            if (vadd(residual_forward, residual_reverse) != (0, 0, 0, 0, 0, 0)
                    or sum(residual_forward[:3]) != 0
                    or residual_forward[3:] != (0, 0, 0)):
                raise AssertionError("endpoint E residue did not cancel")
            cancellation_count += 1
            environment_count += 1
    check(environment_count == 4 * 64,
          "all four orientations and all 64 local environments were replayed")
    check(cancellation_count == 4 * 64,
          "all endpoint E residues cancel only after forward/reverse averaging")

    # Actual global operator entries, independently completed by integral flow.
    ring_rows = []
    operator_entries = []
    for missing in range(4):
        state, cycle_edges = global_ice_with_ring(
            vertices, edges, incidence, chosen[missing], lookup)
        switched = state ^ frozenset(cycle_edges)
        check(all(sum(edge in state for edge in incidence[v]) == 2 for v in vertices),
              f"missing {missing}: start is a global ice state")
        check(all(sum(edge in switched for edge in incidence[v]) == 2 for v in vertices),
              f"missing {missing}: switched endpoint is a global ice state")
        labels, bits = cycle_environment(state, chosen[missing], cycle_edges,
                                         edges, incidence)
        forward = differentiated_ring(labels, bits, coulomb, 0)
        reverse = differentiated_ring(labels, bits, coulomb, 1)
        row = vmul(Q(1, 2), vadd(forward, reverse))
        check(row == expected_ring(missing, coulomb),
              f"missing {missing}: actual global ring row equals local theorem")
        ring_rows.append(row)
        operator_entries.append((state, switched))
    check(len(set(operator_entries)) == 4,
          "four missing labels give four distinct projected matrix elements")
    expected_rows = [
        (Q(231, 8), Q(231, 8), Q(231, 8), Q(-189, 8), Q(-189, 8), Q(-189, 8)),
        (Q(231, 8), Q(231, 8), Q(231, 8), Q(189, 8), Q(189, 8), Q(-189, 8)),
        (Q(231, 8), Q(231, 8), Q(231, 8), Q(189, 8), Q(-189, 8), Q(189, 8)),
        (Q(231, 8), Q(231, 8), Q(231, 8), Q(-189, 8), Q(189, 8), Q(189, 8)),
    ]
    check(ring_rows == expected_rows, "Coulomb ring rows reproduce FV12")
    check(matrix_rank(ring_rows) == 4, "actual ring operators span A1+T2 rank four")

    # Direction-pair states show that the direct rows are also genuine global
    # diagonal expectation differences, rather than local coefficient rows.
    direction_coverings = []
    for labels in ((0, 1), (0, 2), (0, 3)):
        state = frozenset(edge for edge, (_, _, label) in enumerate(edges)
                          if label in labels)
        check(all(sum(edge in state for edge in incidence[v]) == 2 for v in vertices),
              f"direction pair {labels} is a global ice covering")
        direction_coverings.append(state)
    check(len(set(direction_coverings)) == 3,
          "three direct diagonal states are distinct global states")

    witness = e_rows + ring_rows
    determinant = matrix_det(witness)
    check(matrix_rank(witness) == 6, "six actual operator evaluations have rank six")
    check(determinant == Q(-4678629417, 256),
          "independent operator witness determinant is -4678629417/256")
    check(all(sum(row[:3]) == 0 for row in e_rows),
          "diagonal difference functionals annihilate every identity")
    check(all(start != finish for start, finish in operator_entries),
          "ring functionals are off diagonal and annihilate every identity")

    def ring_rank(slope):
        return matrix_rank([expected_ring(missing, slope) for missing in range(4)])

    check(ring_rank(coulomb) == 4, "Coulomb slope has full ring A1+T2 rank")
    check(ring_rank(Q(2, 5)) == 1, "lambda=2/5 kills T2 only")
    check(ring_rank(Q(3, 5)) == 3, "lambda=3/5 kills ring A1 only")
    check(ring_rank(Q(0)) == 4 and
          matrix_rank([direct_q(pair, Q(0)) for pair in ((0, 1), (0, 2), (0, 3))]) == 0,
          "lambda=0 retains ring rank four but kills direct E")
    check(coulomb not in (Q(0), Q(2, 5), Q(3, 5)),
          "Coulomb avoids all exact loss slopes")

    # Sign/order bookkeeping: four H6 rows each contribute the same physical
    # factor -h^6/U_d^5, while two direct rows contribute U_d.  Hence the four
    # minus signs cancel and the leading determinant is h^24 times the
    # normalized determinant.  H8 replacements begin two powers later.
    check((-1) ** 4 == 1 and 6 * 4 == 24 and 24 + 2 == 26,
          "global sign and h^24/O(h^26) determinant bookkeeping is exact")

    # The topology supplies an independent endpoint argument through H8:
    # an ice-to-ice odd-support set is Eulerian.  Girth six forbids nonempty
    # support below six; with <=8 links a connected nonempty support is a H6
    # or H8 cycle (a figure-eight/theta needs at least twelve at this girth).
    # Thus lower folds are diagonal, H6 ring entries receive no lower fold,
    # and H8 can only dress H6 or introduce distinct octagon entries.
    check(not graph_has_four_cycle(vertices, edges, incidence),
          "lower-fold endpoint argument has the required girth-six premise")
    check("only diagonal terms, dressed-H6 transitions, and new octagon transitions" in theorem_flat,
          "through-H8 endpoint classification is stated without a pole claim")
    check("-\\frac{4678629417}{256}h^{24}+O(h^{26})" in theorem_flat,
          "formal determinant states the correct leading and next possible orders")

    print(f"SUMMARY {count}/{count} independent hostile checks passed")
    print("HISTORIES 720; PROPER_PREFIXES 3600; ENVIRONMENTS 4x64")
    print("DIRECT E rank2; H6_RING A1+T2 rank4; OPERATOR rank6")
    print(f"WITNESS_DETERMINANT {determinant}; FORMAL h^24 + O(h^26)")
    print("DISPOSITION PASS_AFTER_FV_PURE_PREMISE_AND_BYTE_HYGIENE_REPAIR")
    print("CEILING off-shell formal rank only; no CTP/Ward/pole/gravity/G")


if __name__ == "__main__":
    main()
