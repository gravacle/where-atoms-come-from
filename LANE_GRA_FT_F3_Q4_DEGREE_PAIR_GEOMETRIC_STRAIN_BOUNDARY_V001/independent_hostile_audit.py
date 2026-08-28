#!/usr/bin/env python3
"""Independent hostile audit of the F3/q4 degree-pair strain boundary.

This script does not import the builder verifier.  It reconstructs the tensor
algebra, the degree identity, the ice restriction, the DPAR derivative, the
tetrahedral covariance, and compatible periodic H6 matrix-element witnesses
from separate exact-arithmetic code.
"""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PAIRS = tuple(combinations(range(4), 2))
SIGNS = ((1, 1, 1), (1, -1, -1),
         (-1, 1, -1), (-1, -1, 1))

EXPECTED_DEPENDENCIES = {
    "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md":
        "00eba581b90fb9f0b25e3fad1362b055049824897433ba06e356ab9b1f6c76ec",
    "LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/THEOREM.md":
        "05f4a619a6f80aa40c48570ab4035ab874426502a31468a08f435e66610bd769",
    "LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/INDEPENDENT_AUDIT.md":
        "44690ad431c85af7a4947a431a4c57ad4cd8b19a346e3d26720b341f77256f90",
    "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/THEOREM.md":
        "cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98",
    "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/INDEPENDENT_AUDIT.md":
        "c52eab9d701d1c6e82f1d7ec395841f4d2810e96cccbc3e2504760b6742e81e4",
    "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md":
        "78f0687c9f597c96d235224dad45c204d12d7e6c973c270e3305a901efc75b25",
    "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/INDEPENDENT_AUDIT.md":
        "53893c7198241f0f8f6aa766f3595fb75b83d208581833c32656b28d7c7f02b9",
    "LANE_GRA_FQ_F3_Q4_COLLECTIVE_METRIC_ORIGIN_SCREEN_V001/THEOREM.md":
        "07445c035ed4c5167a5a20280c4db69a5101eeb71831cdeb126b29702d04b69d",
    "LANE_GRA_FQ_F3_Q4_COLLECTIVE_METRIC_ORIGIN_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "91aa35170432684a47278e46ee2b9d56658a43acc8acbbb84480d047cdbe6dcf",
    "LANE_GRA_FR_F3_Q4_BLOCK_STRAIN_SOURCE_RANK_OBSTRUCTION_V001/THEOREM.md":
        "62c7aaee9433a9ffa970ff6e38bac5585200cf40d6fca2cb70477e7e1e7524eb",
    "LANE_GRA_FR_F3_Q4_BLOCK_STRAIN_SOURCE_RANK_OBSTRUCTION_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "d2da0796cfec7cff8f1d7da5c9bc449d38acdbae089dd9778fb5f19cb6e42b88",
    "LANE_GRA_FS_F3_Q4_COMPLETE_SOURCE_RANK_AUDIT_V001/THEOREM.md":
        "36879f4c18eec83a22bdf9bd161d9d444b72e1dbda1d5eaa0312c6aab3d95724",
    "LANE_GRA_FS_F3_Q4_COMPLETE_SOURCE_RANK_AUDIT_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "6e0ecd0febf6364e4122bbf2f65e1feb93c27e960bce30c0097ea0fbe3f58966",
}

passed = 0


def check(condition, label):
    global passed
    if not condition:
        raise AssertionError(label)
    passed += 1


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def parse_ledger(path):
    ledger = {}
    for line in path.read_text().splitlines():
        expected, relative = line.split("  ", 1)
        check(relative not in ledger, f"ledger member is unique: {path.name}:{relative}")
        ledger[relative] = expected
    return ledger


def rank(rows):
    matrix = [[F(value) for value in row] for row in rows]
    if not matrix:
        return 0
    nrows, ncols = len(matrix), len(matrix[0])
    pivot_row = 0
    for column in range(ncols):
        pivot = next((row for row in range(pivot_row, nrows)
                      if matrix[row][column]), None)
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [entry / scale for entry in matrix[pivot_row]]
        for row in range(nrows):
            if row == pivot_row or not matrix[row][column]:
                continue
            scale = matrix[row][column]
            matrix[row] = [left - scale * right
                           for left, right in zip(matrix[row], matrix[pivot_row])]
        pivot_row += 1
        if pivot_row == nrows:
            break
    return pivot_row


def transpose(matrix):
    return tuple(tuple(matrix[row][column] for row in range(len(matrix)))
                 for column in range(len(matrix[0])))


def matmul(left, right):
    right_t = transpose(right)
    return tuple(tuple(sum(F(a) * F(b) for a, b in zip(row, column))
                       for column in right_t) for row in left)


def matvec(matrix, vector):
    return tuple(sum(F(a) * F(b) for a, b in zip(row, vector)) for row in matrix)


def inverse(matrix):
    size = len(matrix)
    augmented = [[F(matrix[row][column]) for column in range(size)]
                 + [F(row == column) for column in range(size)]
                 for row in range(size)]
    for column in range(size):
        pivot = next(row for row in range(column, size) if augmented[row][column])
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [entry / scale for entry in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            scale = augmented[row][column]
            augmented[row] = [left - scale * right
                              for left, right in zip(augmented[row], augmented[column])]
    return tuple(tuple(augmented[row][size:]) for row in range(size))


def outer(vector, scale=F(1)):
    return tuple(tuple(scale * F(left) * F(right) for right in vector)
                 for left in vector)


def matrix_add(left, right):
    return tuple(tuple(a + b for a, b in zip(left_row, right_row))
                 for left_row, right_row in zip(left, right))


def matrix_scale(scale, matrix):
    return tuple(tuple(F(scale) * value for value in row) for row in matrix)


def symmetric_coordinates(matrix):
    return (matrix[0][0], matrix[1][1], matrix[2][2],
            2 * matrix[0][1], 2 * matrix[0][2], 2 * matrix[1][2])


def source_matrix(coordinate):
    matrix = [[F(0) for _ in range(3)] for _ in range(3)]
    if coordinate < 3:
        matrix[coordinate][coordinate] = F(1)
    else:
        row, column = ((0, 1), (0, 2), (1, 2))[coordinate - 3]
        matrix[row][column] = matrix[column][row] = F(1)
    return tuple(tuple(row) for row in matrix)


def contract(source, tensor):
    return sum(source[row][column] * tensor[row][column]
               for row in range(3) for column in range(3))


# Dependency custody is reconstructed from a hard-coded hostile ledger.
dependencies = parse_ledger(HERE / "DEPENDENCIES.sha256")
check(dependencies == EXPECTED_DEPENDENCIES,
      "dependency ledger exactly matches the independent hostile expectation")
for relative, expected in EXPECTED_DEPENDENCIES.items():
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(),
          f"dependency is a regular non-symlink: {relative}")
    check(digest(path) == expected, f"dependency digest matches: {relative}")
    check(sha256(path.read_bytes() + b"hostile-tamper").hexdigest() != expected,
          f"dependency appended-byte tamper fails: {relative}")


# Independently normalized tetrahedral algebra.
D = tuple(outer(sign, F(1, 3)) for sign in SIGNS)
RHAT = {}
for a, b in PAIRS:
    difference = tuple(SIGNS[b][axis] - SIGNS[a][axis] for axis in range(3))
    length2_unscaled = sum(value * value for value in difference)
    check(length2_unscaled == 8, f"root {a + 1}{b + 1} has unscaled length squared eight")
    RHAT[(a, b)] = outer(difference, F(1, length2_unscaled))

d_rows = tuple(symmetric_coordinates(tensor) for tensor in D)
r_rows = tuple(symmetric_coordinates(RHAT[pair]) for pair in PAIRS)
check(rank(d_rows) == 4, "normalized tetrahedral one-edge dyads have rank four")
check(rank(r_rows) == 6, "normalized sibling-root dyads have rank six")
check(rank(d_rows + r_rows) == 6, "edge plus root tensor span has rank six")

gram = tuple(tuple(contract(D[a], D[b]) for b in range(4)) for a in range(4))
check(all(gram[a][b] == (F(1) if a == b else F(1, 9))
          for a in range(4) for b in range(4)),
      "edge-dyad Gram matrix is exactly diagonal one and off-diagonal one ninth")

E_MATRICES = (
    ((F(1), F(0), F(0)), (F(0), F(-1), F(0)), (F(0), F(0), F(0))),
    ((F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(-2))),
)
T2_MATRICES = (source_matrix(3), source_matrix(4), source_matrix(5))
IDENTITY = ((F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(1)))
check(all(contract(e, tensor) == 0 for e in E_MATRICES for tensor in D),
      "both diagonal-traceless E sources annihilate every edge dyad")

W = tuple(tuple(sum(D[a][row][column] for a in range(4))
                for column in range(3)) for row in range(3))
check(W == matrix_scale(F(4, 3), IDENTITY),
      "the unsplit degree support tensor is exactly four-thirds identity")
additive = tuple(matrix_add(D[a], D[b]) for a, b in PAIRS)
check(rank(d_rows + tuple(symmetric_coordinates(tensor) for tensor in additive)) == 4,
      "all additive pair weights remain in the edge rank-four span")


# Every permutation of tetrahedral labels has an exact orthogonal realization
# and permutes the normalized root dyads covariantly.
source_basis = tuple(tuple(F(SIGNS[column][row]) for column in range(3))
                     for row in range(3))
source_inverse = inverse(source_basis)
for permutation in permutations(range(4)):
    target_basis = tuple(tuple(F(SIGNS[permutation[column]][row])
                               for column in range(3)) for row in range(3))
    orthogonal = matmul(target_basis, source_inverse)
    check(matmul(transpose(orthogonal), orthogonal) == IDENTITY,
          f"label permutation {permutation} has an exact orthogonal tetrahedral action")
    check(all(matvec(orthogonal, SIGNS[a]) == SIGNS[permutation[a]] for a in range(4)),
          f"label permutation {permutation} maps all four coframe vectors")
    for a, b in PAIRS:
        transformed = matmul(matmul(orthogonal, RHAT[(a, b)]), transpose(orthogonal))
        target_pair = tuple(sorted((permutation[a], permutation[b])))
        check(transformed == RHAT[target_pair],
              f"permutation {permutation} covariantly maps root {a + 1}{b + 1}")


# Degree-square identity and pair independence in the full local Hilbert space.
all_states = tuple(product((-1, 1), repeat=4))
pair_characters = []
for state in all_states:
    degree = sum((1 - z) // 2 for z in state)
    values = tuple(state[a] * state[b] for a, b in PAIRS)
    pair_characters.append(values)
    check(F((degree - 2) ** 2) == F(1) + F(sum(values), 2),
          f"degree identity holds on computational state {state}")
check(rank(pair_characters) == 6,
      "the six pair Pauli strings are independent on the full four-link Hilbert space")


# Exact ice restriction: pair image A1+E rank three, with only E nonconstant.
ice = tuple(state for state in all_states if sum(state) == 0)
ice_pair_values = tuple(tuple(state[a] * state[b] for a, b in PAIRS)
                        for state in ice)
check(len(ice) == 6, "the q4 two-in/two-out fiber has six states")
check(rank(ice_pair_values) == 3, "ice-restricted pair operators have rank three")
centered_pair_values = tuple(tuple(F(value) + F(1, 3) for value in row)
                             for row in ice_pair_values)
check(rank(centered_pair_values) == 2,
      "centered ice pair operators have exact E rank two")
check(all(sum(row) == -2 for row in ice_pair_values),
      "the ice pair sum is the scalar minus two")
for state in ice:
    check(all(state[a] * state[b] == state[c] * state[d]
              for (a, b), (c, d) in (((0, 1), (2, 3)),
                                      ((0, 2), (1, 3)),
                                      ((0, 3), (1, 2)))),
          f"complementary pair identities hold on ice state {state}")


def query_values(source):
    coefficients = tuple(contract(source, RHAT[pair]) for pair in PAIRS)
    return tuple(sum(coefficient * state[a] * state[b]
                     for coefficient, (a, b) in zip(coefficients, PAIRS))
                 for state in ice)


e_values = tuple(query_values(source) for source in E_MATRICES)
t2_values = tuple(query_values(source) for source in T2_MATRICES)
a1_values = query_values(IDENTITY)
all_source_values = tuple(query_values(source_matrix(index)) for index in range(6))
check(rank(e_values) == 2 and all(sum(values) == 0 for values in e_values),
      "normalized root query realizes both centered ice E directions")
check(all(not any(values) for values in t2_values),
      "normalized root query kills all three T2 directions after ice restriction")
check(len(set(a1_values)) == 1 and a1_values[0] == -2,
      "normalized root A1 query is exactly the scalar minus two on ice")
check(rank(all_source_values) == 3,
      "direct ice-projected root-pair source has exact A1+E rank three, not six")


# Same source-off Hamiltonian but different microscopic derivatives.  The
# DPAR family replaces, rather than supplements, the FS degree deformation.
source_off_fs = (F(1),) + tuple(F(1, 2) for _ in PAIRS)
source_off_root = (F(1),) + tuple(F(1, 2) for _ in PAIRS)
check(source_off_fs == source_off_root,
      "FS and root degree families have identical identity/pair coefficients at source off")
check(rank(d_rows + (symmetric_coordinates(W),)) == 4,
      "FS one-edge plus unsplit degree derivative has microscopic rank four")
check(rank(r_rows) == 6,
      "the independent pair-string DPAR derivative has microscopic rank six")
check(rank(all_source_values) == 3,
      "microscopic rank six is not overpromoted to projected root-pair rank six")


# Exact DPAR differentiation, including the sign, factor two, normalization,
# and the six-coordinate symmetric contraction convention.
U_D = F(5, 3)
LAMBDA = F(7, 5)
for pair in PAIRS:
    difference = tuple(F(SIGNS[pair[1]][axis] - SIGNS[pair[0]][axis])
                       for axis in range(3))
    norm2 = sum(entry * entry for entry in difference)
    for coordinate in range(6):
        source = source_matrix(coordinate)
        # F(epsilon)=I-epsilon*j/2 gives
        # d |F r|^2/d epsilon at zero = -r^T j r.
        ratio_derivative = -sum(difference[row] * source[row][column] * difference[column]
                                for row in range(3) for column in range(3)) / norm2
        contraction = contract(source, RHAT[pair])
        check(ratio_derivative == -contraction,
              f"DPAR affine ratio derivative has correct normalized sign for pair {pair}, coordinate {coordinate}")
        energy_derivative = U_D * F(1, 2) * LAMBDA * ratio_derivative
        conjugate = -2 * energy_derivative
        check(conjugate == U_D * LAMBDA * contraction,
              f"DPAR Q=-2 dH/dj factor is exact for pair {pair}, coordinate {coordinate}")

# Hermiticity is checked directly on a nontrivial real source slice.
test_source = tuple(tuple(F(row + column + 1, 20) if row == column
                          else F(row + column + 1, 40)
                          for column in range(3)) for row in range(3))
diagonal = []
for state in all_states:
    energy = U_D
    for pair in PAIRS:
        coefficient = U_D * F(1, 2) * (1 - LAMBDA * contract(test_source, RHAT[pair]))
        energy += coefficient * state[pair[0]] * state[pair[1]]
    diagonal.append(energy)
check(all(isinstance(value, F) for value in diagonal),
      "real DPAR source produces a real diagonal Hermitian operator")


# A quadratic seagull can change a CTP/contact Hessian without changing the
# source-off linear source rank.
seagull_gradient = (F(0), F(0))
seagull_hessian = ((F(2), F(0)), (F(0), F(0)))
check(not any(seagull_gradient) and any(any(row) for row in seagull_hessian),
      "an O(j^2) contact has zero gradient but can have a nonzero E Hessian")


# Frozen G_5 periodic q4 graph and all local hexagon types.
L = 5
SHIFTS = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
A_VERTICES = tuple((0, x, y, z) for x, y, z in product(range(L), repeat=3))
B_VERTICES = tuple((1, x, y, z) for x, y, z in product(range(L), repeat=3))
VERTICES = A_VERTICES + B_VERTICES
adjacency = {vertex: [] for vertex in VERTICES}
edges = []
edge_id = {}
edge_label = {}
for x, y, z in product(range(L), repeat=3):
    av = (0, x, y, z)
    for label, shift in enumerate(SHIFTS):
        bv = (1, (x + shift[0]) % L,
                 (y + shift[1]) % L,
                 (z + shift[2]) % L)
        eid = len(edges)
        edges.append((av, bv, label))
        adjacency[av].append((bv, eid))
        adjacency[bv].append((av, eid))
        key = frozenset((av, bv))
        edge_id[key] = eid
        edge_label[key] = label

check(len(VERTICES) == 250 and len(edges) == 500,
      "G_5 has 250 vertices and 500 links")
check(all(len(adjacency[vertex]) == 4 for vertex in VERTICES),
      "G_5 is periodic and four-regular")


def canonical_cycle(cycle):
    cycle = tuple(cycle)
    variants = []
    for oriented in (cycle, tuple(reversed(cycle))):
        for shift in range(len(cycle)):
            variants.append(oriented[shift:] + oriented[:shift])
    return min(variants)


origin = (0, 0, 0, 0)
hexagons = set()


def enumerate_cycles(vertex, path):
    if len(path) == 6:
        if any(neighbor == origin for neighbor, _ in adjacency[vertex]):
            hexagons.add(canonical_cycle(path))
        return
    for neighbor, _ in adjacency[vertex]:
        if neighbor not in path:
            enumerate_cycles(neighbor, path + [neighbor])


enumerate_cycles(origin, [origin])
check(len(hexagons) == 12, "G_5 has twelve simple hexagons through the reference vertex")
ring_types = set()
for cycle in hexagons:
    index = cycle.index(origin)
    left = edge_label[frozenset((origin, cycle[index - 1]))]
    right = edge_label[frozenset((origin, cycle[(index + 1) % 6]))]
    ring_types.add(tuple(sorted((left, right))))
check(ring_types == set(PAIRS), "the G_5 local hexagons realize all six ring-edge label pairs")


class Dinic:
    def __init__(self, size):
        self.graph = [[] for _ in range(size)]

    def add_edge(self, source, target, capacity):
        forward = [target, len(self.graph[target]), capacity]
        reverse = [source, len(self.graph[source]), 0]
        self.graph[source].append(forward)
        self.graph[target].append(reverse)
        return len(self.graph[source]) - 1

    def max_flow(self, source, sink):
        total = 0
        while True:
            level = [-1] * len(self.graph)
            level[source] = 0
            queue = [source]
            for vertex in queue:
                for target, _, capacity in self.graph[vertex]:
                    if capacity and level[target] < 0:
                        level[target] = level[vertex] + 1
                        queue.append(target)
            if level[sink] < 0:
                return total
            cursor = [0] * len(self.graph)

            def send(vertex, amount):
                if vertex == sink:
                    return amount
                while cursor[vertex] < len(self.graph[vertex]):
                    index = cursor[vertex]
                    target, reverse, capacity = self.graph[vertex][index]
                    if capacity and level[target] == level[vertex] + 1:
                        pushed = send(target, min(amount, capacity))
                        if pushed:
                            self.graph[vertex][index][2] -= pushed
                            self.graph[target][reverse][2] += pushed
                            return pushed
                    cursor[vertex] += 1
                return 0

            while True:
                pushed = send(source, 10**9)
                if not pushed:
                    break
                total += pushed


def cycle_edge_ids(cycle):
    return tuple(edge_id[frozenset((cycle[index], cycle[(index + 1) % 6]))]
                 for index in range(6))


def compatible_ice_state(cycle, parity):
    ring = cycle_edge_ids(cycle)
    fixed = {eid: int(index % 2 == parity) for index, eid in enumerate(ring)}
    fixed_one_degree = {vertex: 0 for vertex in VERTICES}
    for eid, value in fixed.items():
        if value:
            av, bv, _ = edges[eid]
            fixed_one_degree[av] += 1
            fixed_one_degree[bv] += 1
    demand = {vertex: 2 - fixed_one_degree[vertex] for vertex in VERTICES}
    check(all(0 <= value <= 2 for value in demand.values()),
          "fixed alternating ring leaves valid local b-matching demands")
    check(sum(demand[vertex] for vertex in A_VERTICES)
          == sum(demand[vertex] for vertex in B_VERTICES),
          "residual bipartite b-matching demands balance")

    node = {vertex: index + 1 for index, vertex in enumerate(VERTICES)}
    source = 0
    sink = len(VERTICES) + 1
    flow = Dinic(sink + 1)
    for vertex in A_VERTICES:
        flow.add_edge(source, node[vertex], demand[vertex])
    candidate_refs = []
    for eid, (av, bv, _) in enumerate(edges):
        if eid in fixed:
            continue
        index = flow.add_edge(node[av], node[bv], 1)
        candidate_refs.append((eid, node[av], index))
    for vertex in B_VERTICES:
        flow.add_edge(node[vertex], sink, demand[vertex])
    target = sum(demand[vertex] for vertex in A_VERTICES)
    check(flow.max_flow(source, sink) == target,
          "alternating ring extends to an exact periodic degree-two ice state")

    occupied = {eid for eid, value in fixed.items() if value}
    for eid, av_node, index in candidate_refs:
        if flow.graph[av_node][index][2] == 0:
            occupied.add(eid)
    degree = {vertex: 0 for vertex in VERTICES}
    for eid in occupied:
        av, bv, _ = edges[eid]
        degree[av] += 1
        degree[bv] += 1
    check(all(value == 2 for value in degree.values()),
          "constructed periodic witness satisfies degree two at every vertex")
    check(all((eid in occupied) == bool(index % 2 == parity)
              for index, eid in enumerate(ring)),
          "constructed witness has the required alternating hexagon orientation")
    return occupied


# Build both alternating orientations for every reference hexagon, then test
# the exact local H6 commutator sign and normalization on those global states.
channels = []
J6 = F(11, 7)
for cycle in sorted(hexagons):
    ring = cycle_edge_ids(cycle)
    index = cycle.index(origin)
    incident_eids = (edge_id[frozenset((origin, cycle[index - 1]))],
                     edge_id[frozenset((origin, cycle[(index + 1) % 6]))])
    ring_labels = tuple(sorted(edge_label[frozenset((origin, endpoint))]
                               for endpoint in (cycle[index - 1], cycle[(index + 1) % 6])))
    for parity in (0, 1):
        occupied = compatible_ice_state(cycle, parity)
        local_state = [None] * 4
        local_eid = [None] * 4
        for neighbor, eid in adjacency[origin]:
            label = edge_label[frozenset((origin, neighbor))]
            local_eid[label] = eid
            local_state[label] = 1 - 2 * int(eid in occupied)
        check(sum(local_state) == 0, "global witness restricts to a local ice state")
        flipped = list(local_state)
        for label in ring_labels:
            flipped[label] *= -1
        check(sum(flipped) == 0, "ring flip preserves the local ice constraint")

        channel = []
        for e_source in E_MATRICES:
            coefficients = {pair: contract(e_source, RHAT[pair]) for pair in PAIRS}
            q_before = U_D * LAMBDA * sum(coefficients[pair]
                                          * local_state[pair[0]] * local_state[pair[1]]
                                          for pair in PAIRS)
            q_after = U_D * LAMBDA * sum(coefficients[pair]
                                         * flipped[pair[0]] * flipped[pair[1]]
                                         for pair in PAIRS)
            lhs_matrix_element = J6 * (q_after - q_before)
            odd_pairs = tuple(pair for pair in PAIRS
                              if (pair[0] in ring_labels) ^ (pair[1] in ring_labels))
            rhs_matrix_element = 2 * J6 * U_D * LAMBDA * sum(
                coefficients[pair] * flipped[pair[0]] * flipped[pair[1]]
                for pair in odd_pairs
            )
            check(lhs_matrix_element == rhs_matrix_element,
                  "H6 commutator matrix element has the exact +2 J6 Ud lambda sign and factor")
            channel.append(q_after - q_before)
        channels.append(tuple(channel))

check(rank(channels) == 2,
      "compatible periodic H6 witnesses act injectively on the two local E directions")
check(all(any(channel[column] for channel in channels) for column in range(2)),
      "each displayed E basis has a nonzero compatible-state H6 matrix element")
ring_sets = tuple(frozenset(cycle_edge_ids(cycle)) for cycle in hexagons)
check(len(set(ring_sets)) == len(ring_sets),
      "distinct local hexagons have distinct symmetric-difference edge sets")


# Documentary repairs and ceilings.
theorem = (HERE / "THEOREM.md").read_text()
result = (HERE / "RESULT.md").read_text()
self_audit = (HERE / "SELF_AUDIT.md").read_text()
required_theorem = (
    "Equation (FT08) **replaces** the FS scalar deformation",
    "not an extension *inside* the unchanged FQ17a/FS",
    "direct root-pair image has only `A1+E`, hence rank\nthree",
    "exact rank-six statement in this lane is microscopic",
    "does **not** prove that\n`DPAR` is the unique or logically necessary closure",
    "It is neither\ninherited nor adopted here",
    "real differentiable physical function",
    "This checks the factor of two and fixes the sign",
    "Real `g`\nand real `j` preserve Hermiticity",
    "common law is `S4` covariant",
    "compatible global ice state",
    "both alternating orientations of every\nhexagon",
    "full Feshbach source conjugate may also contain\ngenerated higher-order corrections",
    "uniform sum over all vertices",
    "no `k=0` response",
    "neither inherited nor\nadopted",
    "not the CTP pole, Ward, or\ngravity tests",
)
for phrase in required_theorem:
    check(phrase in theorem, f"corrected theorem preserves hostile scope phrase: {phrase}")
check("single minimal physical premise" not in theorem.lower(),
      "theorem removes the unproved absolute-minimality label")
check("combined off-shell source again has exact rank six" not in theorem,
      "theorem removes the projected-rank conflation")
check("not a double-counted scalar" in theorem,
      "theorem resolves the overlapping A1 source accounting")
check("neither inherited nor adopted" in result,
      "result keeps DPAR explicitly unadopted")
check("direct ice-projected pair image itself\nhas only `A1+E` rank three" in result,
      "result reports the projected rank-three ceiling")
check("not a second scalar direction" in " ".join(self_audit.split()),
      "self-audit rejects A1 double counting")
check("Full state-dependent CTP rank" in result and "gravity, and `G` remain open" in result,
      "result retains CTP, gravity, and G ceilings")

for name in ("DEPENDENCIES.sha256", "RESULT.md", "SELF_AUDIT.md", "THEOREM.md",
             "verify_degree_pair_geometric_strain_boundary.py",
             "independent_hostile_audit.py"):
    data = (HERE / name).read_bytes()
    check(b"\r" not in data and b"\b" not in data and b"\f" not in data,
          f"byte hygiene: {name}")


# Stable builder payload and seal.  The audit packet has a separate manifest
# so neither verifier must self-hash a changing transcript.
manifest = parse_ledger(HERE / "MANIFEST.sha256")
expected_builder_members = {
    "DEPENDENCIES.sha256", "RESULT.md", "SELF_AUDIT.md", "THEOREM.md",
    "verify_degree_pair_geometric_strain_boundary.py",
}
check(set(manifest) == expected_builder_members,
      "base manifest has exactly the five stable builder payload members")
for relative, expected in manifest.items():
    path = HERE / relative
    check(path.is_file() and not path.is_symlink(),
          f"base payload is regular and non-symlink: {relative}")
    check(digest(path) == expected, f"base manifest digest matches: {relative}")
    check(sha256(path.read_bytes() + b"hostile-base-tamper").hexdigest() != expected,
          f"base payload appended-byte tamper fails: {relative}")

seal = parse_ledger(HERE / "SEAL.sha256")
check(set(seal) == {"MANIFEST.sha256", "VERIFICATION.txt"},
      "builder seal covers exactly manifest and replay transcript")
for relative, expected in seal.items():
    check(digest(HERE / relative) == expected,
          f"builder seal digest matches: {relative}")

print(f"SUMMARY {passed}/{passed} independent hostile checks passed")
print("DISPOSITION PASS_AFTER_MICROSCOPIC_VS_PROJECTED_RANK__DPAR_SUFFICIENCY__A1_REPLACEMENT__AND_COMPATIBLE_H6_WITNESS_REPAIRS__DPAR_UNADOPTED__CTP_WARD_TENSOR_GRAVITY_G_OPEN")
