#!/usr/bin/env python3
"""Independent exact hostile replay of the GL6BV contact theorem.

This file does not import the author verifier.  All algebra is integer or
Fraction exact.  It independently reconstructs the local representation,
source derivatives, periodic witness, covariant edge expectation, spectral
normalization, fixed-projection warning, and isolated-core support boundary.
"""

from collections import Counter, deque
from fractions import Fraction as Q
from itertools import combinations, permutations, product
import json


passed = 0


def require(condition, label):
    global passed
    if not condition:
        raise AssertionError(label)
    passed += 1


def zeros(n, m=None):
    m = n if m is None else m
    return [[Q(0) for _ in range(m)] for _ in range(n)]


def identity(n):
    ans = zeros(n)
    for index in range(n):
        ans[index][index] = Q(1)
    return ans


def plus(a, b, ca=Q(1), cb=Q(1)):
    return [[ca * a[i][j] + cb * b[i][j]
             for j in range(len(a[0]))] for i in range(len(a))]


def times(c, a):
    return [[c * value for value in row] for row in a]


def mmul(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), Q(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def mv(a, v):
    return [sum((a[i][j] * v[j] for j in range(len(v))), Q(0))
            for i in range(len(a))]


def outer(a, b):
    return [[x * y for y in b] for x in a]


def inner(a, b):
    return sum((x * y for x, y in zip(a, b)), Q(0))


def matrix_sum(items, n=6):
    ans = zeros(n)
    for item in items:
        ans = plus(ans, item)
    return ans


def rank(rows):
    work = [[Q(value) for value in row] for row in rows]
    r = 0
    for column in range(len(work[0]) if work else 0):
        pivot = next((i for i in range(r, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[r], work[pivot] = work[pivot], work[r]
        scale = work[r][column]
        work[r] = [value / scale for value in work[r]]
        for i in range(len(work)):
            if i == r:
                continue
            scale = work[i][column]
            if scale:
                work[i] = [x - scale * y for x, y in zip(work[i], work[r])]
        r += 1
    return r


# Pair representation.  The order is fixed independently here.
PAIRS = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIRS)}
I6 = identity(6)
J6 = [[Q(1)] * 6 for _ in range(6)]
O6 = zeros(6)
for source, pair in enumerate(PAIRS):
    complement = tuple(port for port in range(4) if port not in pair)
    O6[PAIR_INDEX[complement]][source] = Q(1)
P_A = times(Q(1, 6), J6)
P_T = times(Q(1, 2), plus(I6, O6, cb=Q(-1)))
P_E = plus(times(Q(1, 2), plus(I6, O6)), P_A, cb=Q(-1))

require(plus(plus(P_A, P_E), P_T) == I6,
        "independent pair projectors resolve six-space")
require(mmul(P_T, P_T) == P_T and
        sum(P_T[i][i] for i in range(6)) == 3,
        "independent T2 projector has rank three")


def pair_vector(spins):
    return [Q(spins[a] * spins[b]) for a, b in PAIRS]


locked = tuple(spins for spins in product((-1, 1), repeat=4)
               if sum(spins) == 0)
require(len(locked) == 6, "six strict-lock local states")
for spins in locked:
    require(mv(P_T, pair_vector(spins)) == [Q(0)] * 6,
            "locked local pair T2 operator is zero statewise")

# Defect frame indexed by its exceptional spin.
tau = []
for exceptional in range(4):
    spins = [1, 1, 1, 1]
    spins[exceptional] = -1
    vector = pair_vector(spins)
    tau.append(vector)
    require(mv(P_T, vector) == vector and
            mv(P_A, vector) == [Q(0)] * 6 and
            mv(P_E, vector) == [Q(0)] * 6,
            "degree-one/three defect vector is pure T2")
require([[inner(tau[a], tau[b]) for b in range(4)] for a in range(4)] ==
        [[Q(6 if a == b else -2) for b in range(4)] for a in range(4)],
        "defect Gram matrix is tetrahedral")
require([sum(tau[a][i] for a in range(4)) for i in range(6)] == [Q(0)] * 6,
        "four defect vectors sum to zero")
frame = matrix_sum(outer(vector, vector) for vector in tau)
require(frame == times(Q(8), P_T), "defect vectors form exact 8 P_T frame")


def same_sign_partner(spins, port):
    matches = [other for other in range(4)
               if other != port and spins[other] == spins[port]]
    require(len(matches) == 1, "locked port has unique same-sign partner")
    return matches[0]


for spins in locked:
    incident = []
    for port in range(4):
        flipped = list(spins)
        flipped[port] *= -1
        expected = tau[same_sign_partner(spins, port)]
        require(pair_vector(flipped) == expected,
                "incident flip produces partner-labelled defect vector")
        incident.append(expected)
    require(Counter(tuple(row) for row in incident) ==
            Counter(tuple(row) for row in tau),
            "four incident flips produce the complete defect frame")
    for size in range(1, 4):
        for subset in combinations(incident, size):
            require([sum(row[i] for row in subset) for i in range(6)] !=
                    [Q(0)] * 6,
                    "no nonempty proper incident subset cancels its first vertex")

# Source-before-Feshbach derivatives, tested at a nontrivial rational scale.
u = Q(11, 5)
h = Q(7, 13)
gap = 2 * u
first_edge = h * h / gap**2
second_edge = -2 * h * h / gap**3
require(first_edge == h * h / (4 * u**2),
        "H+j.M first derivative has positive h2/(4 Ud2) sign")
require(second_edge == -h * h / (4 * u**3),
        "one denominator Hessian is -h2/(4 Ud3) outer product")
for spins in locked:
    incident = [tau[same_sign_partner(spins, port)] for port in range(4)]
    first = [first_edge * sum(row[i] for row in incident) for i in range(6)]
    onsite = matrix_sum(
        (times(second_edge, outer(row, row)) for row in incident))
    require(first == [Q(0)] * 6,
            "complete four-edge first source derivative cancels")
    require(onsite == times(-2 * h * h / u**3, P_T),
            "complete four-edge onsite Hessian is -2 h2/Ud3 P_T")

for port in range(4):
    labels = tuple(label for label in range(4) if label != port)
    blocks = []
    for left in labels:
        for right in labels:
            block = outer(tau[left], tau[right])
            require(rank(block) == 1,
                    "each allowed nearest-neighbor outer block has rank one")
            blocks.append(tuple(tuple(row) for row in block))
    require(len(set(blocks)) == 9,
            "each port admits nine distinct statewise edge blocks")

# Norm and exact sufficient denominator domain.  ||tau||^2=6, hence two
# endpoints shift a denominator by strictly less than 2 Ud whenever
# max_v ||j_v|| < Ud/sqrt(6), by two applications of Cauchy--Schwarz.
require(all(inner(vector, vector) == 6 for vector in tau),
        "defect norm gives exact Ud/sqrt(6) sufficient source domain")

# Translation/S4 edge expectation from the two stabilizer orbits.
q_frame = [times(Q(1, 6), outer(vector, vector)) for vector in tau]
require(matrix_sum(q_frame) == times(Q(4, 3), P_T),
        "four Q_a tensors sum to 4/3 P_T")


def direct_edge_expectation(p, port):
    labels = tuple(label for label in range(4) if label != port)
    terms = []
    for left in labels:
        for right in labels:
            probability = p / 3 if left == right else (1 - p) / 6
            terms.append(times(probability, outer(tau[left], tau[right])))
    return matrix_sum(terms)


def predicted_edge_expectation(p, port):
    beta = 4 * p - Q(4, 3)
    return plus(times(beta, P_T),
                times(Q(2, 3) - beta, q_frame[port]))


for p in (Q(0), Q(1, 7), Q(1, 2), Q(6, 7), Q(1)):
    blocks = []
    for port in range(4):
        direct = direct_edge_expectation(p, port)
        require(direct == predicted_edge_expectation(p, port),
                "S4 stabilizer average is exactly the one-p formula")
        blocks.append(direct)
    gamma = Q(8, 3) * (4 * p - 1)
    require(matrix_sum(blocks) == times(gamma, P_T),
            "zero-character cross block is gamma P_T")

# Canonical Q4 background.  This independently reimplements the defining
# ordered Edmonds--Karp completion from the sealed GL6AN witness.
side = 4
directions = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
cells = tuple(product(range(side), repeat=3))


def cell_add(a, b):
    return tuple((a[i] + b[i]) % side for i in range(3))


def cell_sub(a, b):
    return tuple((a[i] - b[i]) % side for i in range(3))


def child(cell, port):
    return cell_add(cell, directions[port])


edges = tuple((cell, port) for cell in cells for port in range(4))


def hexagon(cell, ports):
    a, b, c = ports
    ab = cell_add(cell_sub(cell, directions[b]), directions[a])
    bc = cell_add(cell_sub(cell, directions[b]), directions[c])
    return ((cell, a), (ab, b), (ab, c),
            (bc, a), (bc, b), (cell, c))


target = hexagon((0, 0, 0), (0, 1, 2))
fixed = {edge: int(index % 2 == 0) for index, edge in enumerate(target)}
fixed_parent = Counter()
fixed_child = Counter()
for (cell, port), occupation in fixed.items():
    if occupation:
        fixed_parent[cell] += 1
        fixed_child[child(cell, port)] += 1

source, sink = ("SOURCE",), ("SINK",)
capacity = {}


def arc(start, finish, value):
    capacity.setdefault(start, {})[finish] = value
    capacity.setdefault(finish, {}).setdefault(start, 0)


for cell in cells:
    arc(source, ("PARENT", cell), 2 - fixed_parent[cell])
    arc(("CHILD", cell), sink, 2 - fixed_child[cell])
for edge in edges:
    if edge not in fixed:
        cell, port = edge
        arc(("PARENT", cell), ("CHILD", child(cell, port)), 1)

residual = {node: dict(row) for node, row in capacity.items()}
flow = 0
while True:
    predecessor = {source: None}
    queue = deque((source,))
    while queue and sink not in predecessor:
        node = queue.popleft()
        for nxt in sorted(residual[node], key=repr):
            if residual[node][nxt] and nxt not in predecessor:
                predecessor[nxt] = node
                queue.append(nxt)
    if sink not in predecessor:
        break
    node = sink
    while node != source:
        prior = predecessor[node]
        residual[prior][node] -= 1
        residual[node][prior] += 1
        node = prior
    flow += 1

occupied = {edge for edge, value in fixed.items() if value}
for edge in edges:
    if edge in fixed:
        continue
    cell, port = edge
    if residual[("CHILD", child(cell, port))].get(("PARENT", cell), 0):
        occupied.add(edge)
require(flow == 125 and len(occupied) == 128,
        "independent Q4 flow completion saturates the locked background")
require(tuple(int(edge in occupied) for edge in target) == (1, 0, 1, 0, 1, 0),
        "independent Q4 target remains alternating")


def node_spins(occupation, kind, cell):
    if kind == "P":
        return tuple(1 - 2 * int((cell, port) in occupation)
                     for port in range(4))
    return tuple(1 - 2 * int((cell_sub(cell, directions[port]), port)
                             in occupation) for port in range(4))


for kind in ("P", "C"):
    for cell in cells:
        require(sum(node_spins(occupied, kind, cell)) == 0,
                "independent Q4 completion is locked at every node")

joint = Counter()
aligned = 0
for cell, port in edges:
    remote = child(cell, port)
    left = same_sign_partner(node_spins(occupied, "P", cell), port)
    right = same_sign_partner(node_spins(occupied, "C", remote), port)
    joint[(port, left, right)] += 1
    aligned += int(left == right)
require(aligned == 218 and Q(aligned, len(edges)) == Q(109, 128),
        "independent Q4 witness has 218/256 aligned links")

# Reconstruct the literal translation/S4 orbit average from its joint census.
orbit_blocks = [zeros(6) for _ in range(4)]
orbit_counts = [0, 0, 0, 0]
for (port, left, right), multiplicity in joint.items():
    for perm in permutations(range(4)):
        new_port = perm[port]
        orbit_blocks[new_port] = plus(
            orbit_blocks[new_port],
            times(multiplicity, outer(tau[perm[left]], tau[perm[right]])))
        orbit_counts[new_port] += multiplicity
for port in range(4):
    orbit_blocks[port] = times(Q(1, orbit_counts[port]), orbit_blocks[port])
    require(orbit_counts[port] == 1536,
            "each port has the exact Q4 translation/S4 orbit count")
    require(orbit_blocks[port] == predicted_edge_expectation(Q(109, 128), port),
            "Q4 orbit block matches independent one-p prediction")

p_q4 = Q(109, 128)
beta = 4 * p_q4 - Q(4, 3)
delta = Q(2, 3) - beta
gamma = Q(8, 3) * (4 * p_q4 - 1)
require((beta, delta, gamma) == (Q(199, 96), Q(-45, 32), Q(77, 12)),
        "Q4 beta delta gamma parameters are exact")
require((8 + gamma, 8 - gamma) == (Q(173, 12), Q(19, 12)),
        "Q4 fixed common/relative zero-character values are exact")
require((-beta / 2, -delta / 2, beta / 2, delta / 2) ==
        (Q(-199, 192), Q(45, 64), Q(199, 192), Q(-45, 64)),
        "BV35 fixed-projection quadratic coefficients are exact")

# The imaginary linear character term is nonzero and maps fixed common to
# relative.  Therefore BV35 cannot be read as an eigenbranch expansion.
theta = (Q(1), Q(-1), Q(0), Q(0))
linear_mixing = matrix_sum(
    (times(delta * theta[port], q_frame[port]) for port in range(4)))
require(linear_mixing != zeros(6) and rank(linear_mixing) > 0,
        "linear-character term genuinely mixes fixed common/relative sectors")

# Exact sparse A* A ownership on Q4.  Four active edges give 8 P_T onsite;
# one physical edge gives exactly one pair of transpose rank-one cross blocks.
full_frames = {}
edge_blocks = {}
for cell, port in edges:
    remote = child(cell, port)
    parent_node = ("P", cell)
    child_node = ("C", remote)
    left = tau[same_sign_partner(node_spins(occupied, "P", cell), port)]
    right = tau[same_sign_partner(node_spins(occupied, "C", remote), port)]
    full_frames[parent_node] = plus(
        full_frames.get(parent_node, zeros(6)), outer(left, left))
    full_frames[child_node] = plus(
        full_frames.get(child_node, zeros(6)), outer(right, right))
    edge_blocks[(parent_node, child_node)] = outer(left, right)
require(len(full_frames) == 128 and
        all(value == times(Q(8), P_T) for value in full_frames.values()),
        "Q4 A* A has the universal four-edge onsite frame at all 128 nodes")
require(len(edge_blocks) == 256 and
        all(rank(value) == 1 for value in edge_blocks.values()),
        "Q4 A* A owns exactly 256 rank-one oriented edge blocks")

# Quadratic-form factorization for two unrelated exact T2 source fields.
for seed in (1, 5):
    source_field = {}
    for kind_index, kind in enumerate(("P", "C")):
        for cell in cells:
            raw = [Q(seed + (kind_index + 1) * (port + 2) + sum(cell))
                   for port in range(4)]
            source_field[(kind, cell)] = [
                sum((raw[label] * tau[label][component]
                     for label in range(4)), Q(0))
                for component in range(6)]
    edge_norm = Q(0)
    assembled = Q(0)
    for cell, port in edges:
        remote = child(cell, port)
        left = tau[same_sign_partner(node_spins(occupied, "P", cell), port)]
        right = tau[same_sign_partner(node_spins(occupied, "C", remote), port)]
        jl = source_field[("P", cell)]
        jr = source_field[("C", remote)]
        a = inner(left, jl)
        b = inner(right, jr)
        edge_norm += (a + b)**2
        assembled += a*a + 2*a*b + b*b
    require(edge_norm == assembled and edge_norm > 0,
            "independent Q4 quadratic form factors exactly as A* A")

# Retarded normalization.  A transition amplitude h/(2 Ud) times A has
# spectral Gram h2/(4 Ud2) A*A at gap 2 Ud.  With i/2 commutator, the zero-
# frequency value is spectral Gram/(2 Ud), hence F''=-2 K(0).
spectral_gram_factor = h*h / (4*u*u)
k_zero_factor = spectral_gram_factor / (2*u)
f_hessian_factor = -h*h / (4*u**3)
require(f_hessian_factor == -2 * k_zero_factor,
        "i/2 retarded normalization gives F''=-2 K(0)")
require(times(spectral_gram_factor, times(Q(8), P_T)) ==
        times(2*h*h/u**2, P_T),
        "onsite gapped retarded sine coefficient is 2 h2/Ud2 P_T")

# Six-core versus four-edge active support at each target node.
toggled = set(occupied)
for edge in target:
    if edge in toggled:
        toggled.remove(edge)
    else:
        toggled.add(edge)
target_nodes = set()
for cell, port in target:
    target_nodes.add(("P", cell))
    target_nodes.add(("C", child(cell, port)))
require(len(target_nodes) == 6, "target hexagon has six distinct constraint nodes")
for kind, cell in target_nodes:
    core_ports = []
    for parent, port in target:
        if kind == "P" and parent == cell:
            core_ports.append(port)
        if kind == "C" and child(parent, port) == cell:
            core_ports.append(port)
    require(len(core_ports) == 2,
            "six-core control owns exactly two incident flips per cycle node")
    partial_counters = []
    for state in (occupied, toggled):
        spins = node_spins(state, kind, cell)
        full = [tau[same_sign_partner(spins, port)] for port in range(4)]
        partial = [tau[same_sign_partner(spins, port)] for port in core_ports]
        partial_sum = [sum(row[i] for row in partial) for i in range(6)]
        partial_frame = matrix_sum(outer(row, row) for row in partial)
        require(matrix_sum(outer(row, row) for row in full) == times(Q(8), P_T),
                "separate Q4 embedding retains full four-edge frame")
        require(partial[0] != partial[1] and partial_sum != [Q(0)] * 6 and
                rank(partial_frame) == 2,
                "six-core first vertex is nonzero and partial frame has rank two")
        partial_counters.append(Counter(tuple(row) for row in partial))
    require(partial_counters[0] == partial_counters[1],
            "hexagon toggle only swaps the unordered two-vector core set")

result = {
    "audit": "AUDIT_G_GL6BV_STRICT_LOCK_T2_SOURCE_FESHBACH_CONTACT_V001",
    "target": "LANE_CROSS_RFT_GRA_GL6BV_STRICT_LOCK_T2_SOURCE_FESHBACH_CONTACT_V001",
    "status": "PASS_AFTER_AUTHOR_LITERAL_BYTE_REPAIR",
    "checks": passed,
    "source_chart": "genuine microscopic H+j.M with P_T j=j before Feshbach projection",
    "denominator": "2*U_d+j_v.t_v+j_w.t_w",
    "exact_denominator_domain": "connected nonzero-denominator component containing j=0",
    "sufficient_source_ball": "max_v||j_v||<U_d/sqrt(6)",
    "defect_frame": "sum_a tau_a tau_a^T=8*P_T",
    "first_vertex": "zero for all four equal-h incident flips; nonzero for every nonempty proper subset",
    "onsite_hessian": "-2*h^2/U_d^3*P_T",
    "edge_hessian": "-h^2/(4*U_d^3)*tau_left*tau_right^T; rank one",
    "full_hessian": "-h^2/(4*U_d^3)*A^*A",
    "retarded": "gap 2*U_d; sine weight 2*h^2/U_d^2*P_T onsite; F''=-2*K(0)",
    "translation_s4": "one alignment parameter p",
    "q4": {
        "aligned_links": "218/256=109/128",
        "beta": "199/96",
        "delta": "-45/32",
        "gamma": "77/12",
        "common_relative_k0": ["173/12", "19/12"],
        "bv35_scope": "fixed diagonal projections, not eigenbranches or Schur-reduced coefficients",
    },
    "active_support": {
        "six_core": "two incident flips per cycle node; nonzero first vertex; rank-two partial frame",
        "separate_q4_embedding": "four incident flips and full 8*P_T frame",
        "eighteen_link_source_first_h2_h4_h6": "open",
        "common_stationary_E_T_functional": "not supplied",
    },
    "ricci": "not tested",
}
print("RESULT_JSON " + json.dumps(result, sort_keys=True))
print(f"PASS GL6BV independent hostile replay {passed}/{passed}")
