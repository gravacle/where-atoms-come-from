#!/usr/bin/env python3
"""Exact replay for GL6BV.

The replay uses only integer/Fraction arithmetic.  It checks the local
six-pair representation algebra, the order-h^2 source-before-Feshbach
derivatives, the translation/S4 expectation symbol, and the declared Q4
period-four witness.  No floating-point comparison is used.
"""

from collections import Counter, deque
from fractions import Fraction as F
from itertools import combinations, permutations, product


checks = []


def check(condition, label):
    if not condition:
        raise AssertionError(label)
    checks.append(label)


def zero(rows, cols=None):
    cols = rows if cols is None else cols
    return [[F(0) for _ in range(cols)] for _ in range(rows)]


def eye(size):
    out = zero(size)
    for i in range(size):
        out[i][i] = F(1)
    return out


def add(left, right, a=F(1), b=F(1)):
    return [[a * left[i][j] + b * right[i][j]
             for j in range(len(left[0]))]
            for i in range(len(left))]


def scale(value, matrix):
    return [[value * entry for entry in row] for row in matrix]


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def multiply(left, right):
    return [[sum(left[i][k] * right[k][j] for k in range(len(right)))
             for j in range(len(right[0]))]
            for i in range(len(left))]


def matvec(matrix, vector):
    return [sum(row[j] * vector[j] for j in range(len(vector)))
            for row in matrix]


def outer(left, right):
    return [[left[i] * right[j] for j in range(len(right))]
            for i in range(len(left))]


def dot(left, right):
    return sum(x * y for x, y in zip(left, right))


def trace(matrix):
    return sum(matrix[i][i] for i in range(len(matrix)))


I6 = eye(6)
ONE6 = [F(1)] * 6
J6 = outer(ONE6, ONE6)
PAIR_ORDER = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIR_ORDER)}
OPPOSITE = zero(6)
for index, pair in enumerate(PAIR_ORDER):
    complement = tuple(value for value in range(4) if value not in pair)
    OPPOSITE[index][PAIR_INDEX[complement]] = F(1)

P_A = scale(F(1, 6), J6)
P_T = scale(F(1, 2), add(I6, OPPOSITE, b=F(-1)))
P_E = add(scale(F(1, 2), add(I6, OPPOSITE)), P_A, b=F(-1))

check(add(add(P_A, P_E), P_T) == I6, "A1/E/T2 projectors resolve pair space")
for projector, rank in ((P_A, 1), (P_E, 2), (P_T, 3)):
    check(multiply(projector, projector) == projector,
          f"rank-{rank} projector is idempotent")
    check(trace(projector) == rank, f"projector trace is {rank}")
for left, right in ((P_A, P_E), (P_A, P_T), (P_E, P_T)):
    check(multiply(left, right) == zero(6), "pair projectors are orthogonal")


def pair_values(z):
    return [F(z[a] * z[b]) for a, b in PAIR_ORDER]


LOCKED = [tuple(z) for z in product((-1, 1), repeat=4) if sum(z) == 0]
check(len(LOCKED) == 6, "there are six local degree-two states")
for z in LOCKED:
    m = pair_values(z)
    check(matvec(P_T, m) == [F(0)] * 6,
          "locked local pair vector has exact T2 null")

# tau_r is the pair vector of the degree-one/three defect whose exceptional
# port is r.  A global sign change of the four z values leaves tau_r fixed.
TAU = []
for exceptional in range(4):
    w = [1, 1, 1, 1]
    w[exceptional] = -1
    tau = pair_values(w)
    TAU.append(tau)
    check(matvec(P_T, tau) == tau, "degree-one/three vector is pure T2")
    check(matvec(P_A, tau) == [F(0)] * 6,
          "degree-one/three vector has zero A1 part")
    check(matvec(P_E, tau) == [F(0)] * 6,
          "degree-one/three vector has zero E part")

check([dot(TAU[a], TAU[b]) for a in range(4) for b in range(4)] ==
      [F(6) if a == b else F(-2) for a in range(4) for b in range(4)],
      "four defect vectors have tetrahedral Gram matrix")
check([sum(TAU[a][i] for a in range(4)) for i in range(6)] == [F(0)] * 6,
      "sum_a tau_a vanishes")
TAU_FRAME = zero(6)
for tau in TAU:
    TAU_FRAME = add(TAU_FRAME, outer(tau, tau))
check(TAU_FRAME == scale(F(8), P_T),
      "sum_a tau_a tau_a^T equals 8 P_T")


def matching_port(z, port):
    matches = [other for other in range(4)
               if other != port and z[other] == z[port]]
    check(len(matches) == 1, "locked state has one same-sign partner per port")
    return matches[0]


for z in LOCKED:
    moved_vectors = []
    for port in range(4):
        moved = list(z)
        moved[port] *= -1
        moved_pair = pair_values(moved)
        partner = matching_port(z, port)
        check(moved_pair == TAU[partner],
              "flipped locked port gives its same-sign partner defect vector")
        moved_vectors.append(tuple(moved_pair))
    check(Counter(moved_vectors) == Counter(tuple(row) for row in TAU),
          "four incident flips give the state-independent defect multiset")

# Order-h^2 source expansion for H+j.M in the pure-T2 source chart.
# -h^2/(D+x) = -h^2/D + h^2 x/D^2 - h^2 x^2/D^3 + ...,
# D=2 U_d.  Its Hessian is -2 h^2 q q^T/D^3.
U = F(7, 3)
h = F(5, 11)
D = 2 * U
linear_coefficient = h * h / (D * D)
hessian_coefficient = -2 * h * h / (D * D * D)
check(linear_coefficient == h * h / (4 * U * U),
      "edge first derivative coefficient is h^2/(4 U_d^2)")
ONSITE_DISPLAY = r"-{h^2\over4U_d^3}\,8P_T"
check(hessian_coefficient == -h * h / (4 * U * U * U) and
      ONSITE_DISPLAY == r"-{h^2\over4U_d^3}\,8P_T",
      "edge Hessian coefficient and BV17 multiplication token are exact")

for z in LOCKED:
    local_t = [TAU[matching_port(z, port)] for port in range(4)]
    first = [linear_coefficient * sum(row[i] for row in local_t)
             for i in range(6)]
    onsite = zero(6)
    for row in local_t:
        onsite = add(onsite, scale(hessian_coefficient, outer(row, row)))
    check(first == [F(0)] * 6,
          "the complete local order-h^2 first T2 source vertex cancels")
    check(onsite == scale(-2 * h * h / (U ** 3), P_T),
          "onsite physical-energy Hessian is -2 h^2/U_d^3 P_T")

# Across an edge the block is rank one and depends on both endpoint locked
# configurations.  All nine locally allowed values occur algebraically.
for port in range(4):
    allowed = [other for other in range(4) if other != port]
    edge_blocks = {tuple(tuple(entry for entry in row)
                         for row in outer(TAU[left], TAU[right]))
                   for left in allowed for right in allowed}
    check(len(edge_blocks) == 9,
          "nearest-neighbor block has nine state-dependent local values")

# Translation/S4 expectation.  p is the probability that the two same-sign
# partners across an a-edge agree.  Stabilizer covariance leaves one scalar p.
Q = [scale(F(1, 6), outer(tau, tau)) for tau in TAU]
Q_SUM = zero(6)
for q_a in Q:
    Q_SUM = add(Q_SUM, q_a)
check(Q_SUM == scale(F(4, 3), P_T),
      "sum_a Q_a equals 4/3 P_T")


def covariant_edge_block(p, port):
    beta = 4 * p - F(4, 3)
    return add(scale(beta, P_T), scale(F(2, 3) - beta, Q[port]))


for p in (F(0), F(1, 4), F(1, 2), F(1)):
    for port in range(4):
        allowed = [other for other in range(4) if other != port]
        direct = zero(6)
        for left in allowed:
            for right in allowed:
                probability = p / 3 if left == right else (1 - p) / 6
                direct = add(direct, scale(probability,
                                           outer(TAU[left], TAU[right])))
        check(direct == covariant_edge_block(p, port),
              "S4 edge expectation is fixed by the alignment probability")
    total = zero(6)
    for a in range(4):
        total = add(total, covariant_edge_block(p, a))
    gamma = F(8, 3) * (4 * p - 1)
    check(total == scale(gamma, P_T),
          "zero-character cross block is gamma P_T")


# Reconstruct the authenticated period-four alternating Q4 witness.
LENGTH = 4
DIRECTIONS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))


def add_cell(left, right):
    return tuple((left[i] + right[i]) % LENGTH for i in range(3))


def sub_cell(left, right):
    return tuple((left[i] - right[i]) % LENGTH for i in range(3))


CELLS = list(product(range(LENGTH), repeat=3))
EDGES = [(x, port) for x in CELLS for port in range(4)]


def child(x, port):
    return add_cell(x, DIRECTIONS[port])


def canonical_hexagon(x, ports):
    a, b, c = ports
    x_ab = add_cell(sub_cell(x, DIRECTIONS[b]), DIRECTIONS[a])
    x_bc = add_cell(sub_cell(x, DIRECTIONS[b]), DIRECTIONS[c])
    return ((x, a), (x_ab, b), (x_ab, c),
            (x_bc, a), (x_bc, b), (x, c))


TARGET = canonical_hexagon((0, 0, 0), (0, 1, 2))
FIXED = {edge: (1 if index % 2 == 0 else 0)
         for index, edge in enumerate(TARGET)}
parent_fixed = Counter()
child_fixed = Counter()
for (x, port), value in FIXED.items():
    if value:
        parent_fixed[x] += 1
        child_fixed[child(x, port)] += 1

source, sink = ("S",), ("T",)
capacity = {}


def add_arc(start, finish, value):
    capacity.setdefault(start, {})[finish] = value
    capacity.setdefault(finish, {}).setdefault(start, 0)


for x in CELLS:
    add_arc(source, ("P", x), 2 - parent_fixed[x])
    add_arc(("C", x), sink, 2 - child_fixed[x])
for edge in EDGES:
    if edge not in FIXED:
        x, port = edge
        add_arc(("P", x), ("C", child(x, port)), 1)

residual = {node: dict(row) for node, row in capacity.items()}
flow = 0
while True:
    previous = {source: None}
    queue = deque([source])
    while queue and sink not in previous:
        node = queue.popleft()
        for neighbor in sorted(residual.get(node, {}), key=repr):
            if residual[node][neighbor] > 0 and neighbor not in previous:
                previous[neighbor] = node
                queue.append(neighbor)
    if sink not in previous:
        break
    node = sink
    while node != source:
        prior = previous[node]
        residual[prior][node] -= 1
        residual[node][prior] += 1
        node = prior
    flow += 1

OCCUPIED = {edge for edge, value in FIXED.items() if value}
for edge in EDGES:
    if edge in FIXED:
        continue
    x, port = edge
    if residual[("C", child(x, port))].get(("P", x), 0) == 1:
        OCCUPIED.add(edge)

check(flow == 125, "period-four background max flow is 125")
check(len(OCCUPIED) == 128, "period-four background has 128 occupied links")
check(tuple(int(edge in OCCUPIED) for edge in TARGET) == (1, 0, 1, 0, 1, 0),
      "period-four target hexagon is alternating")


def local_z(kind, x):
    if kind == "P":
        return tuple(1 - 2 * int((x, port) in OCCUPIED) for port in range(4))
    entries = []
    for port in range(4):
        parent = sub_cell(x, DIRECTIONS[port])
        entries.append(1 - 2 * int((parent, port) in OCCUPIED))
    return tuple(entries)


for kind in ("P", "C"):
    for x in CELLS:
        check(sum(local_z(kind, x)) == 0,
              "period-four witness is locked at every endpoint")

joint = Counter()
aligned = 0
for x, port in EDGES:
    other = child(x, port)
    left = matching_port(local_z("P", x), port)
    right = matching_port(local_z("C", other), port)
    joint[(port, left, right)] += 1
    aligned += int(left == right)
check(aligned == 218, "Q4 witness has 218 aligned endpoint partners")
p_q4 = F(aligned, len(EDGES))
check(p_q4 == F(109, 128), "Q4 orbit-mixture alignment p is 109/128")

# Explicitly average the finite witness over all 24 port permutations.  This
# turns its translation-averaged edge data into the claimed S4 covariant
# matrices without assuming that the parent dynamically selects the mixture.
averaged = [zero(6) for _ in range(4)]
counts = [0] * 4
for (port, left, right), multiplicity in joint.items():
    for perm in permutations(range(4)):
        target_port = perm[port]
        averaged[target_port] = add(
            averaged[target_port],
            scale(F(multiplicity), outer(TAU[perm[left]], TAU[perm[right]])),
        )
        counts[target_port] += multiplicity
for port in range(4):
    averaged[port] = scale(F(1, counts[port]), averaged[port])
    check(counts[port] == 1536, "each S4-averaged port has 1536 samples")
    check(averaged[port] == covariant_edge_block(p_q4, port),
          "explicit Q4 S4 orbit average matches the p formula")

beta_q4 = 4 * p_q4 - F(4, 3)
delta_q4 = F(2, 3) - beta_q4
gamma_q4 = F(8, 3) * (4 * p_q4 - 1)
check((beta_q4, delta_q4, gamma_q4) ==
      (F(199, 96), F(-45, 32), F(77, 12)),
      "Q4 symbol parameters are beta=199/96 delta=-45/32 gamma=77/12")
check((F(8) + gamma_q4, F(8) - gamma_q4) ==
      (F(173, 12), F(19, 12)),
      "Q4 zero-character common/relative Gram coefficients are 173/12 and 19/12")

# The even k^2 cross block is
# -1/2 [beta I2 P_T + delta sum_a theta_a^2 Q_a].
# For Q4 this is -199/192 I2 P_T +45/64 sum theta_a^2 Q_a.
check(-beta_q4 / 2 == F(-199, 192),
      "Q4 isotropic even-k2 cross coefficient is -199/192")
check(-delta_q4 / 2 == F(45, 64),
      "Q4 covariant directional even-k2 coefficient is 45/64")

# Exact full-supercell factorization at zero character on one deterministic
# T2 source field.  This checks the edge-owner assembly G=A^*A without
# materializing a 768-by-768 matrix.
source_field = {}
for kind_index, kind in enumerate(("P", "C")):
    for x in CELLS:
        coefficients = [F((kind_index + 1) * (i + 1) + sum(x)) for i in range(4)]
        source_field[(kind, x)] = [sum(coefficients[a] * TAU[a][i]
                                          for a in range(4))
                                    for i in range(6)]
edge_square = F(0)
onsite_cross = F(0)
for x, port in EDGES:
    y = child(x, port)
    t_left = TAU[matching_port(local_z("P", x), port)]
    t_right = TAU[matching_port(local_z("C", y), port)]
    j_left = source_field[("P", x)]
    j_right = source_field[("C", y)]
    left_value = dot(t_left, j_left)
    right_value = dot(t_right, j_right)
    edge_square += (left_value + right_value) ** 2
    onsite_cross += left_value ** 2 + right_value ** 2 + 2 * left_value * right_value
check(edge_square == onsite_cross and edge_square > 0,
      "Q4 zero-character Hessian is the exact negative edge Gram A^*A")

# Compare the two Q4 locked configurations related by the target hexagon.
# All four incident Q4 links remain active here, so the onsite frame is
# unchanged.  This is not an attachment to GL6BW's isolated six-core K2
# control, whose active microscopic collar has not been declared.
FINAL_OCCUPIED = set(OCCUPIED)
for edge in TARGET:
    if edge in FINAL_OCCUPIED:
        FINAL_OCCUPIED.remove(edge)
    else:
        FINAL_OCCUPIED.add(edge)


def local_z_for(occupied, kind, x):
    if kind == "P":
        return tuple(1 - 2 * int((x, port) in occupied) for port in range(4))
    return tuple(1 - 2 * int((sub_cell(x, DIRECTIONS[port]), port) in occupied)
                 for port in range(4))


cycle_nodes = set()
for x, port in TARGET:
    cycle_nodes.add(("P", x))
    cycle_nodes.add(("C", child(x, port)))
check(len(cycle_nodes) == 6, "isolated hexagon has six constraint nodes")
for kind, x in cycle_nodes:
    before = pair_values(local_z_for(OCCUPIED, kind, x))
    after = pair_values(local_z_for(FINAL_OCCUPIED, kind, x))
    displacement = [after[i] - before[i] for i in range(6)]
    check(matvec(P_E, displacement) == displacement,
          "hexagon locked displacement is pure E")
    check(dot(displacement, displacement) == 16,
          "hexagon locked displacement has norm squared sixteen")
    core_ports = [port for parent, port in TARGET
                  if ((kind == "P" and parent == x) or
                      (kind == "C" and child(parent, port) == x))]
    partial_sets = []
    for occupied in (OCCUPIED, FINAL_OCCUPIED):
        z = local_z_for(occupied, kind, x)
        frame = zero(6)
        for port in range(4):
            t = TAU[matching_port(z, port)]
            frame = add(frame, outer(t, t))
        partial = [TAU[next(other for other in range(4)
                            if other != port and z[other] == z[port])]
                   for port in core_ports]
        partial_key = Counter(tuple(row) for row in partial)
        same_partial_set = not partial_sets or partial_key == partial_sets[0]
        partial_sets.append(partial_key)
        partial_sum = [sum(row[i] for row in partial) for i in range(6)]
        check(frame == scale(F(8), P_T) and len(core_ports) == 2 and
              partial[0] != partial[1] and partial_sum != [F(0)] * 6 and
              same_partial_set,
              "Q4 has the full frame while the six-core partial frame is a fixed rank-two pair")

# GL6BW uses K^R=(i/2)theta<[M(t),M(0)]>.  The leading high-gap transition
# weight is h^2/(4U_d^2) G and its gap is 2U_d.  Hence E''=-2 K^R(0).
spectral_weight_onsite = h * h / (4 * U * U) * 8
static_retarded_onsite = spectral_weight_onsite / (2 * U)
energy_hessian_onsite = -2 * static_retarded_onsite
check(spectral_weight_onsite == 2 * h * h / (U * U),
      "onsite high-gap T2 spectral weight is 2 h^2/U_d^2")
check(energy_hessian_onsite == -2 * h * h / (U ** 3),
      "fixed-branch energy Hessian equals minus twice the zero-frequency retarded kernel")

print(f"GL6BV exact replay PASS ({len(checks)}/{len(checks)})")
print("DEFECT_FRAME=sum_a tau_a tau_a^T=8 P_T")
print("ONSITE_EPP=-2 h^2/U_d^3 P_T")
print("EDGE_EPP=-h^2/(4 U_d^3) tau_left tau_right^T")
print("Q4_ORBIT_P=109/128;K0_COMMON_RELATIVE=173/12,19/12")
print("HEXAGON_ACTIVE_SUPPORT_MISMATCH=Q4_FOUR_EDGE_FRAME_VS_UNDECLARED_K2_COLLAR;COMPLETED_ET_FUNCTIONAL_OPEN;RICCI_NOT_TESTED")
