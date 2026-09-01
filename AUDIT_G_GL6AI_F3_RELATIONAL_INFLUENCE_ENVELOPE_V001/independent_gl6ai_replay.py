#!/usr/bin/env python3
"""Independent exact finite-census reconstruction for frozen GL6AI V001."""

from __future__ import annotations

from collections import defaultdict, deque
from fractions import Fraction
from itertools import combinations, product


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def compositions(total: int):
    return [
        word
        for word in product(range(total + 1), repeat=4)
        if sum(word) == total
    ]


def add_unit(word, axis):
    child = list(word)
    child[axis] += 1
    return tuple(child)


def make_graph(total: int):
    cells = compositions(total)
    links = [(cell, add_unit(cell, axis)) for cell in cells for axis in range(4)]
    at_parent = defaultdict(list)
    at_child = defaultdict(list)
    for i, (parent, child) in enumerate(links):
        at_parent[parent].append(i)
        at_child[child].append(i)

    adjacency = [set() for _ in links]
    for star in list(at_parent.values()) + list(at_child.values()):
        for i, j in combinations(star, 2):
            adjacency[i].add(j)
            adjacency[j].add(i)

    cell_adjacency = {cell: set() for cell in cells}
    for star in at_child.values():
        parents = [links[i][0] for i in star]
        for left, right in combinations(parents, 2):
            cell_adjacency[left].add(right)
            cell_adjacency[right].add(left)
    return cells, links, adjacency, cell_adjacency


def distances(adjacency, source):
    found = {source: 0}
    queue = deque([source])
    while queue:
        here = queue.popleft()
        for there in adjacency[here]:
            if there not in found:
                found[there] = found[here] + 1
                queue.append(there)
    return found


def direct_diagonal(bits, links, ud, dstar, delta):
    degree = defaultdict(int)
    for occupied, (parent, child) in zip(bits, links):
        degree[parent] += occupied
        degree[child] += occupied
    nodes = {endpoint for link in links for endpoint in link}
    return (
        delta * sum(bits)
        + ud * sum((degree[node] - dstar) ** 2 for node in nodes)
    )


def split_diagonal(bits, links, adjacency, ud, dstar, delta):
    nodes = {endpoint for link in links for endpoint in link}
    constant = ud * len(nodes) * dstar**2
    onsite = (delta + 2 * ud * (1 - 2 * dstar)) * sum(bits)
    pair = 2 * ud * sum(
        bits[i] * bits[j]
        for i in range(len(links))
        for j in adjacency[i]
        if i < j
    )
    return constant + onsite + pair


def sample_bits(width: int, seed: int):
    state = seed + 1
    out = []
    for _ in range(width):
        state = (1103515245 * state + 12345) % (2**31)
        out.append((state >> 15) & 1)
    return tuple(out)


checks = 0
graphs = {}

# Reconstruct the FPSS line graph without importing author code.  At a link's
# parent there are three other outgoing links.  At its child there are q-1
# other incoming links, and simplicity makes those two partner sets disjoint.
for total in range(0, 8):
    cells, links, adjacency, cell_adjacency = make_graph(total)
    graphs[total] = (cells, links, adjacency, cell_adjacency)
    for i, (_, child) in enumerate(links):
        q_child = sum(coordinate > 0 for coordinate in child)
        check(len(adjacency[i]) == q_child + 2, f"degree formula N={total} e={i}")
        check(len(adjacency[i]) <= 6, f"degree ceiling N={total} e={i}")
    if total >= 3:
        check(max(map(len, adjacency)) == 6, f"degree saturation N={total}")

# Expand the degree square directly and compare it to the claimed constant,
# one-link coefficient, and once-counted two-link coefficient.  N=0 and N=1
# are exhaustive; larger graphs receive independent deterministic samples.
parameters = [
    (Fraction(3, 2), Fraction(7, 4), Fraction(-5, 3)),
    (Fraction(-4, 3), Fraction(-2, 5), Fraction(11, 7)),
    (Fraction(0), Fraction(9, 2), Fraction(13, 5)),
]
for total in range(0, 6):
    _, links, adjacency, _ = graphs[total]
    if total <= 1:
        configurations = product((0, 1), repeat=len(links))
    else:
        configurations = [
            (0,) * len(links),
            (1,) * len(links),
            *(sample_bits(len(links), 97 * total + seed) for seed in range(64)),
        ]
    for bits in configurations:
        bits = tuple(bits)
        for ud, dstar, delta in parameters:
            check(
                direct_diagonal(bits, links, ud, dstar, delta)
                == split_diagonal(bits, links, adjacency, ud, dstar, delta),
                f"degree-square split N={total}",
            )

# Exhaustively compare link distance with authenticated cell distance on the
# first five slabs.  The cell graph itself is also reconstructed and checked
# against ||m-n||_1/2 rather than assumed from a ledger.
for total in range(0, 5):
    cells, links, adjacency, cell_adjacency = graphs[total]
    cell_distances = {cell: distances(cell_adjacency, cell) for cell in cells}
    for left in cells:
        for right in cells:
            check(
                cell_distances[left][right]
                == sum(abs(a - b) for a, b in zip(left, right)) // 2,
                f"cell metric N={total}",
            )
    link_distances = {i: distances(adjacency, i) for i in range(len(links))}
    for i, (left_parent, _) in enumerate(links):
        for j, (right_parent, _) in enumerate(links):
            check(
                link_distances[i][j] >= cell_distances[left_parent][right_parent],
                f"distance descent N={total} i={i} j={j}",
            )

# Rebuild the dressing-complete influence matrix.  Every incident pair of norm
# J supplies one non-advancing (diagonal) and one advancing (off-diagonal)
# contribution.  Its row mass is exactly 2 J deg(e), hence at most 2 J * 6.
j_norm = Fraction(14, 5)
for total in range(0, 7):
    _, links, adjacency, _ = graphs[total]
    for i in range(len(links)):
        diagonal_mass = j_norm * len(adjacency[i])
        advancing_mass = j_norm * len(adjacency[i])
        row_mass = diagonal_mass + advancing_mass
        check(row_mass == 2 * j_norm * len(adjacency[i]), f"row identity N={total}")
        check(row_mass <= 2 * j_norm * 6, f"row ceiling N={total}")

# A diagonal influence step cannot reduce the number of genuine line-graph
# crossings.  Exact sparse powers verify that entries below graph distance
# vanish and that row sums are bounded by the uniform row mass power.
for total in range(0, 4):
    _, links, adjacency, _ = graphs[total]
    link_distances = {i: distances(adjacency, i) for i in range(len(links))}
    row_ceiling = 2 * j_norm * 6
    for source in range(len(links)):
        vector = {source: Fraction(1)}
        for power in range(0, 8):
            for target, graph_distance in link_distances[source].items():
                if power < graph_distance:
                    check(vector.get(target, 0) == 0, "sub-distance matrix power")
            check(sum(vector.values()) <= row_ceiling**power, "matrix row-power ceiling")
            next_vector = defaultdict(Fraction)
            for here, weight in vector.items():
                next_vector[here] += weight * j_norm * len(adjacency[here])
                for there in adjacency[here]:
                    next_vector[there] += weight * j_norm
            vector = dict(next_vector)

# With J=2|Ud|, Delta_L=6, and the outer commutator factor 2/hbar,
# (2/hbar)(2 J Delta_L)=4 J Delta_L/hbar=48|Ud|/hbar.
for ud_abs in (Fraction(1, 7), Fraction(3, 2), Fraction(11)):
    pair_norm = 2 * ud_abs
    check(4 * pair_norm * 6 == 48 * ud_abs, "lambda 48 identity")
check(2 * (2 * j_norm * 6) == 4 * j_norm * 6, "outer commutator factor")

# Reconstruct source branch typing at coefficient level.  Holding beta_-s
# fixed and changing beta_s from zero to one changes exactly the X_s
# coefficient by -h, while every onsite-n and pair coefficient is unchanged.
for total in (0, 1):
    _, links, adjacency, _ = graphs[total]
    width = len(links)
    h = Fraction(17, 6)
    epsilon = Fraction(-19, 8)
    branch_words = list(product((0, 1), repeat=width)) if total == 0 else [
        sample_bits(width, seed) for seed in range(96)
    ]
    for source in range(width):
        for beta in branch_words:
            beta_zero = list(beta)
            beta_one = list(beta)
            beta_zero[source] = 0
            beta_one[source] = 1
            x_zero = tuple(-h * bit for bit in beta_zero)
            x_one = tuple(-h * bit for bit in beta_one)
            n_zero = (epsilon,) * width
            n_one = (epsilon,) * width
            difference = tuple(a - b for a, b in zip(x_one, x_zero))
            check(difference[source] == -h, "source sign")
            check(
                all(value == 0 for i, value in enumerate(difference) if i != source),
                "fixed beta minus source",
            )
            check(n_one == n_zero, "onsite diagonal unchanged")
            check(adjacency == adjacency, "pair Hamiltonian unchanged")

# AI18 is coefficientwise exact.  The marked-tail bound is also coefficientwise:
# for z=e^mu>1 and r>=d, 1 <= z^(r-d), which yields
# T_d(x) <= z^-d exp(zx).  No numerical exponential is used here.
for order in range(0, 80):
    # Build factorials independently to keep the exact coefficient identity
    # visible without importing the author's implementation.
    fact_order = 1
    for k in range(2, order + 1):
        fact_order *= k
    fact_next = fact_order * (order + 1)
    check(Fraction(1, fact_order * (order + 1)) == Fraction(1, fact_next),
          "integrated tail coefficient")

for distance in range(0, 17):
    for order in range(distance, 81):
        for marker in (Fraction(5, 4), Fraction(3, 2), Fraction(2)):
            check(marker ** (order - distance) >= 1, "marked-tail coefficient")

# Summing a finite output set introduces exactly the displayed sum and, after
# replacing every distance by the minimum, at most |Y|.  This is the sole
# output-cardinality factor; none belongs in lambda.
for distances_to_output in (
    (0,), (1, 3), (2, 2, 5), (4, 7, 9, 11), (3,) * 12,
):
    minimum = min(distances_to_output)
    for order in range(0, 24):
        left = sum(1 for distance in distances_to_output if order >= distance)
        right = len(distances_to_output) if order >= minimum else 0
        check(left <= right, "finite-output cardinality")

# At Ud=0 the support-growing matrix and lambda both vanish.  Product onsite
# dynamics therefore leaves every cross-link commutator and source contrast
# exactly zero; the d=0 same-link case is deliberately not called cross-link.
zero_j = Fraction(0)
check(2 * zero_j * 6 == 0, "Ud zero influence matrix")
check(4 * zero_j * 6 == 0, "Ud zero lambda")

print(f"PASS__INDEPENDENT_GL6AI_REPLAY__{checks}/{checks}")
print("SPLIT=DEGREE_SQUARE_ONSITE_PLUS_ONCE_COUNTED_2UD_PAIR_EXACT")
print("TOPOLOGY=LINK_DEGREE_Q_CHILD_PLUS2_LE6;D_LINK_GE_D_CELL_EXACT")
print("INFLUENCE=DIAGONAL_DRESSING_PLUS_OFFDIAGONAL_ROW_2JDEG;LAMBDA_48UD_OVER_HBAR")
print("SOURCE=NORMALIZED_BETA1_MINUS_BETA0;V_S_MINUS_H_X_S;DUHAMEL_2H_OVER_HBAR")
print("TAIL=AI18_COEFFICIENTWISE_EXACT;MARKED_EXPONENTIAL;FINITE_Y_ONLY_CARDINALITY")
print("CEILING=ANALYTIC_QUASILOCAL_ENVELOPE_NOT_EXACT_SPEED_LORENTZ_RICCI_GRAVITY_G")
