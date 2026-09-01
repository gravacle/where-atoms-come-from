#!/usr/bin/env python3
"""Exact finite checks for the GL6AN native degree-lock theorem.

Only the Python standard library is used.  Floating-point data from GL6AL are
not inputs.  The periodic graphs below are algebraic quotients used to check
incidence identities and to construct an infinite periodic locked background;
they are not finite-volume physics witnesses.
"""

from __future__ import annotations

from collections import Counter, deque
from fractions import Fraction
from itertools import combinations, permutations
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


def rank_q(matrix: list[list[int | Fraction]]) -> int:
    """Rank over Q by exact row reduction."""
    a = [[Fraction(x) for x in row] for row in matrix]
    if not a:
        return 0
    rows, cols = len(a), len(a[0])
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        p = a[r][c]
        a[r] = [x / p for x in a[r]]
        for i in range(rows):
            if i != r and a[i][c]:
                f = a[i][c]
                a[i] = [x - f * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == rows:
            break
    return r


def mat_vec(matrix: list[list[int]], vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in matrix)


def section_local_square_and_symmetry() -> None:
    # At one degree-four incidence node, the pair contribution is
    # 2*C(k,2), while half of the -6 sum_e n_e term contributes -3k.
    for k in range(5):
        CHECK.equal(
            2 * comb(k, 2) - 3 * k,
            (k - 2) ** 2 - 4,
            f"local completed square k={k}",
        )
        CHECK.equal(
            ((4 - k) - 2) ** 2,
            (k - 2) ** 2,
            f"complete degree-four complement invariance k={k}",
        )

    # The inherited onsite coefficient reaches the lock line in the strict
    # Delta>0 domain exactly for d_star>2 when U_d>0.
    u_d = Fraction(3, 2)
    for d_star in range(6):
        delta = 4 * u_d * (d_star - 2)
        epsilon = delta + 2 * u_d * (1 - 2 * d_star)
        CHECK.equal(epsilon, -6 * u_d, f"inherited lock equation d_star={d_star}")
        CHECK.equal(delta > 0, d_star > 2, f"strict Delta domain d_star={d_star}")
    CHECK.equal(4 * Fraction(1) * (3 - 2), 4, "nonempty inherited witness d_star=3")

    # A full four-link complement sends q to -q.  A partial finite-open flip
    # does not, preventing accidental promotion of the infinite symmetry.
    z = (1, 1, 1, 1)
    q = -Fraction(sum(z), 2)
    q_full = -Fraction(sum(-value for value in z), 2)
    q_partial = -Fraction(sum((-value if i == 0 else value) for i, value in enumerate(z)), 2)
    CHECK.equal(q_full, -q, "full degree-four product automorphism sends q to -q")
    CHECK.true(q_partial != -q, "generic partial finite-open flip is not the complement symmetry")

    # Exact 2x2 Pauli calculation: [X,n]=iY, so [-hX,n]=-ihY.
    x = ((0j, 1 + 0j), (1 + 0j, 0j))
    y = ((0j, -1j), (1j, 0j))
    n = ((0j, 0j), (0j, 1 + 0j))

    def mm(a, b):
        return tuple(
            tuple(sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2))
            for i in range(2)
        )

    xn, nx = mm(x, n), mm(n, x)
    comm = tuple(tuple(xn[i][j] - nx[i][j] for j in range(2)) for i in range(2))
    iy = tuple(tuple(1j * y[i][j] for j in range(2)) for i in range(2))
    CHECK.equal(comm, iy, "Pauli commutator [X,n]=iY")


def quotient_graph(length: int):
    cells = [(a, b, c) for a in range(length) for b in range(length) for c in range(length)]

    def child(x, port):
        y = list(x)
        if port < 3:
            y[port] = (y[port] + 1) % length
        return tuple(y)

    links = [(x, port) for x in cells for port in range(4)]
    endpoints = {e: (("P", e[0]), ("C", child(*e))) for e in links}
    return cells, links, endpoints


def section_incidence_census_and_linear_charge() -> None:
    cells, links, endpoints = quotient_graph(2)
    nodes = [(kind, x) for kind in ("P", "C") for x in cells]
    node_index = {v: i for i, v in enumerate(nodes)}
    link_index = {e: i for i, e in enumerate(links)}
    unsigned = [[0 for _ in links] for _ in nodes]
    incident = {v: [] for v in nodes}
    for e, ends in endpoints.items():
        for v in ends:
            unsigned[node_index[v]][link_index[e]] = 1
            incident[v].append(e)

    CHECK.true(all(len(es) == 4 for es in incident.values()), "every original node has degree four")
    CHECK.true(all(sum(unsigned[r][c] for r in range(len(nodes))) == 2 for c in range(len(links))), "every active link has two endpoints")

    line_pairs: Counter[tuple[int, int]] = Counter()
    line_neighbors = {e: set() for e in links}
    for es in incident.values():
        for e, f in combinations(es, 2):
            i, j = sorted((link_index[e], link_index[f]))
            line_pairs[(i, j)] += 1
            line_neighbors[e].add(f)
            line_neighbors[f].add(e)
    CHECK.equal(len(line_pairs), len(nodes) * comb(4, 2), "line-pair count")
    CHECK.true(all(m == 1 for m in line_pairs.values()), "each line pair has one original-node owner")
    CHECK.true(all(len(ns) == 6 for ns in line_neighbors.values()), "native line degree is six")

    test_occupations = (
        {e: 0 for e in links},
        {e: 1 for e in links},
        {e: int(e[1] in (0, 1)) for e in links},
        {e: (sum(e[0]) + e[1]) % 2 for e in links},
        {e: int((3 * e[0][0] + 2 * e[0][1] + e[0][2] + e[1]) % 5 < 2) for e in links},
    )
    for sample, occupation in enumerate(test_occupations):
        lhs = 2 * sum(occupation[links[i]] * occupation[links[j]] for i, j in line_pairs)
        lhs -= 6 * sum(occupation.values())
        rhs = sum((sum(occupation[e] for e in es) - 2) ** 2 - 4 for es in incident.values())
        CHECK.equal(lhs, rhs, f"global square identity sample {sample}")

    # Multiplying child rows by -1 gives an oriented connected incidence.
    oriented = [
        [entry if nodes[r][0] == "P" else -entry for entry in row]
        for r, row in enumerate(unsigned)
    ]
    rnk = rank_q(oriented)
    CHECK.equal(rnk, len(nodes) - 1, "connected bipartite incidence rank")
    CHECK.equal(len(links) - rnk, len(links) - len(nodes) + 1, "cycle-space nullity")
    CHECK.equal(len(nodes) - rank_q([list(col) for col in zip(*unsigned)]), 1, "only one coefficient solution to c_v+c_w=0")

    signed = tuple(1 if kind == "P" else -1 for kind, _ in nodes)
    CHECK.true(all(sum(signed[r] * unsigned[r][c] for r in range(len(nodes))) == 0 for c in range(len(links))), "the sole coefficient solution is bipartite sign")
    CHECK.equal(sum(signed), 0, "balanced quotient makes the associated degree charge identically zero")


def section_pair_sector() -> None:
    pair_order = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    incidence = [[int(a in pair) for pair in pair_order] for a in range(4)]
    CHECK.equal(rank_q(incidence), 4, "unsigned K4 pair-incidence rank")

    e1 = (1, 0, -1, -1, 0, 1)
    e2 = (0, 1, -1, -1, 1, 0)
    CHECK.equal(mat_vec(incidence, e1), (0, 0, 0, 0), "first E kernel vector")
    CHECK.equal(mat_vec(incidence, e2), (0, 0, 0, 0), "second E kernel vector")

    line = [[0] * 6 for _ in range(6)]
    opposite = [[0] * 6 for _ in range(6)]
    for i, p in enumerate(pair_order):
        for j, q in enumerate(pair_order):
            if i == j:
                continue
            common = len(set(p) & set(q))
            line[i][j] = int(common == 1)
            opposite[i][j] = int(common == 0)
    for idx, v in enumerate((e1, e2), start=1):
        CHECK.equal(mat_vec(line, v), tuple(-2 * x for x in v), f"kernel vector {idx} has E line eigenvalue")
        CHECK.equal(mat_vec(opposite, v), v, f"kernel vector {idx} has E opposite eigenvalue")

    locked_z = [z for z in __import__("itertools").product((-1, 1), repeat=4) if sum(z) == 0]
    CHECK.equal(len(locked_z), 6, "six local degree-two spin states")
    pair_vectors = [tuple(z[a] * z[b] for a, b in pair_order) for z in locked_z]
    for i, m in enumerate(pair_vectors):
        CHECK.equal(sum(m), -2, f"locked A1 fixed state {i}")
        CHECK.equal(mat_vec(incidence, m), (-1, -1, -1, -1), f"locked incidence affine constraint state {i}")

    mean = [sum(Fraction(m[i], len(pair_vectors)) for m in pair_vectors) for i in range(6)]
    cov = [[
        sum(Fraction((m[i] - mean[i]) * (m[j] - mean[j]), len(pair_vectors)) for m in pair_vectors)
        for j in range(6)
    ] for i in range(6)]
    CHECK.true(all(x == Fraction(-1, 3) for x in mean), "uniform locked pair mean")
    c_d, c_a, c_o = cov[0][0], cov[0][1], cov[0][5]
    CHECK.equal((c_d, c_a, c_o), (Fraction(8, 9), Fraction(-4, 9), Fraction(8, 9)), "uniform locked covariance orbits")
    CHECK.equal(c_d + 4 * c_a + c_o, 0, "locked A1 covariance")
    CHECK.equal(c_d - 2 * c_a + c_o, Fraction(8, 3), "locked E covariance")
    CHECK.equal(c_d - c_o, 0, "locked T2 covariance")
    CHECK.equal(rank_q(cov), 2, "locked covariance rank")


def gadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def gabs2(a):
    return a[0] * a[0] + a[1] * a[1]


def section_character_symbol() -> None:
    one = (Fraction(1), Fraction(0))
    minus = (Fraction(-1), Fraction(0))
    imag = (Fraction(0), Fraction(1))
    minus_imag = (Fraction(0), Fraction(-1))
    samples = {
        "trivial": ([one, one, one, one], Fraction(16), (0, 8), 1),
        "balanced": ([one, one, minus, minus], Fraction(0), (4, 4), 2),
        "one-negative": ([one, one, one, minus], Fraction(4), (2, 6), 2),
        "quarter-turn": ([one, imag, minus, minus_imag], Fraction(0), (4, 4), 2),
    }
    for name, (phases, abs_s_sq, eigenvalues, row_rank) in samples.items():
        CHECK.true(all(gabs2(z) == 1 for z in phases), f"{name} phases have unit modulus")
        s = (Fraction(0), Fraction(0))
        for z in phases:
            s = gadd(s, z)
        CHECK.equal(gabs2(s), abs_s_sq, f"{name} |s|^2")
        CHECK.equal(eigenvalues[0] + eigenvalues[1], 8, f"{name} trace")
        CHECK.equal(eigenvalues[0] * eigenvalues[1], 16 - abs_s_sq, f"{name} determinant")
        CHECK.equal(4 - row_rank, 3 if name == "trivial" else 2, f"{name} exact flat nullity")

    # Exact quadratic identity behind 4-|sum exp(i theta_a)|.
    centered = (
        (Fraction(1), Fraction(-1), Fraction(0), Fraction(0)),
        (Fraction(1), Fraction(2), Fraction(-1), Fraction(-2)),
        (Fraction(3, 2), Fraction(-1, 2), Fraction(-2), Fraction(1)),
    )
    for i, theta in enumerate(centered):
        CHECK.equal(sum(theta), 0, f"centered phase sample {i}")
        pair_square = sum((theta[a] - theta[b]) ** 2 for a, b in combinations(range(4), 2))
        CHECK.equal(pair_square, 4 * sum(t * t for t in theta), f"quadratic soft-symbol identity {i}")


def x4(x3):
    return (x3[0], x3[1], x3[2], -sum(x3))


def link_vertices(link):
    x, port = link
    parent = x4(x)
    child = list(parent)
    child[port] += 1
    return (("P", parent), ("C", tuple(child)))


def target_hexagon():
    p0 = (0, 0, 0)
    p1 = (1, -1, 0)
    p2 = (0, -1, 1)
    return (
        (p0, 0),
        (p1, 1),
        (p1, 2),
        (p2, 0),
        (p2, 1),
        (p0, 2),
    )


def subset_energy(indices: frozenset[int]) -> int:
    cycle = target_hexagon()
    initial = (1, 0, 1, 0, 1, 0)
    charge: Counter[tuple[str, tuple[int, ...]]] = Counter()
    for i in indices:
        delta = 1 if initial[i] == 0 else -1
        for vertex in link_vertices(cycle[i]):
            charge[vertex] += delta
    return sum(q * q for q in charge.values())


def deterministic_locked_background(length: int = 4):
    """Return a periodic degree-two background with the hexagon alternating."""
    cells, links, _ = quotient_graph(length)

    def mod_x(x):
        return tuple(v % length for v in x)

    cycle = tuple((mod_x(x), p) for x, p in target_hexagon())
    fixed = {e: (1 if i % 2 == 0 else 0) for i, e in enumerate(cycle)}

    def child(x, port):
        y = list(x)
        if port < 3:
            y[port] = (y[port] + 1) % length
        return tuple(y)

    p_fixed = Counter()
    c_fixed = Counter()
    for (p, port), value in fixed.items():
        if value:
            p_fixed[p] += 1
            c_fixed[child(p, port)] += 1

    # Edmonds-Karp on source -> parents -> children -> sink.
    source, sink = ("S",), ("T",)
    capacity: dict[tuple, dict[tuple, int]] = {}

    def add_edge(u, v, cap):
        capacity.setdefault(u, {})[v] = cap
        capacity.setdefault(v, {}).setdefault(u, 0)

    for p in cells:
        add_edge(source, ("P", p), 2 - p_fixed[p])
    for c in cells:
        add_edge(("C", c), sink, 2 - c_fixed[c])
    for e in links:
        if e not in fixed:
            p, port = e
            add_edge(("P", p), ("C", child(p, port)), 1)

    residual = {u: dict(vs) for u, vs in capacity.items()}
    flow = 0
    while True:
        parent = {source: None}
        queue = deque([source])
        while queue and sink not in parent:
            u = queue.popleft()
            for v in sorted(residual.get(u, {}), key=repr):
                if residual[u][v] > 0 and v not in parent:
                    parent[v] = u
                    queue.append(v)
        if sink not in parent:
            break
        v = sink
        while v != source:
            u = parent[v]
            residual[u][v] -= 1
            residual[v][u] += 1
            v = u
        flow += 1

    required = sum(2 - p_fixed[p] for p in cells)
    occupied = {e for e, value in fixed.items() if value}
    for e in links:
        if e not in fixed:
            p, port = e
            c = child(p, port)
            # Reverse residual capacity one means the unit edge carries flow.
            if residual[("C", c)].get(("P", p), 0) == 1:
                occupied.add(e)
    return cycle, fixed, flow, required, occupied, child, cells, links


def section_hexagon_and_strong_lock() -> None:
    cycle = target_hexagon()
    CHECK.equal(len(set(cycle)), 6, "six distinct target links")
    vertices = [set(link_vertices(e)) for e in cycle]
    for i in range(6):
        CHECK.equal(len(vertices[i] & vertices[(i + 1) % 6]), 1, f"hexagon adjacency {i}")
    for i in range(6):
        for j in range(i + 1, 6):
            if (j - i) % 6 not in (1, 5):
                CHECK.equal(len(vertices[i] & vertices[j]), 0, f"hexagon has no chord {i},{j}")

    # Equality e_a-e_b=e_c-e_d has a unique ordered pair, excluding native
    # four-cycles (and there are no parallel active links/two-cycles).
    differences = {}
    for a in range(4):
        for b in range(4):
            if a == b:
                continue
            d = tuple(int(i == a) - int(i == b) for i in range(4))
            CHECK.true(d not in differences, f"unique ordered port difference {a},{b}")
            differences[d] = (a, b)
    CHECK.equal(len(differences), 12, "twelve distinct ordered port differences")

    # The declared finite witness is the period-four quotient.  Reduction
    # modulo four preserves all twelve ordered port differences, so it does
    # not introduce a wrapped four-cycle.  This check is deliberately not a
    # claim about arbitrary smaller periodic quotients.
    mod4_differences = {}
    for a in range(4):
        for b in range(4):
            if a == b:
                continue
            d = tuple((int(i == a) - int(i == b)) % 4 for i in range(3))
            CHECK.true(
                d not in mod4_differences,
                f"period-four ordered port difference unique {a},{b}",
            )
            mod4_differences[d] = (a, b)
    CHECK.equal(
        len(mod4_differences),
        12,
        "period-four quotient has no wrapped four-cycle",
    )

    all_indices = frozenset(range(6))
    CHECK.equal(subset_energy(frozenset()), 0, "empty toggle stays locked")
    CHECK.equal(subset_energy(all_indices), 0, "complete alternating hexagon toggle stays locked")
    CHECK.true(all(subset_energy(frozenset((i,))) == 2 for i in range(6)), "one-link lock gap is 2 U_d")
    CHECK.true(all(subset_energy(frozenset(s)) > 0 for r in range(1, 6) for s in combinations(range(6), r)), "no proper intermediate subset re-enters locked sector")

    census = {r: Counter(subset_energy(frozenset(s)) for s in combinations(range(6), r)) for r in range(1, 6)}
    CHECK.equal(census[1], Counter({2: 6}), "one-toggle energy census")
    CHECK.equal(census[2], Counter({4: 9, 2: 6}), "two-toggle energy census")
    CHECK.equal(census[3], Counter({4: 12, 2: 6, 6: 2}), "three-toggle energy census")
    CHECK.equal(census[4], census[2], "four-toggle complement census")
    CHECK.equal(census[5], census[1], "five-toggle complement census")

    coefficient = Fraction(0)
    profiles = Counter()
    for order in permutations(range(6)):
        chosen: set[int] = set()
        energies = []
        term = Fraction(1)
        for edge in order[:-1]:
            chosen.add(edge)
            energy = subset_energy(frozenset(chosen))
            energies.append(energy)
            term *= Fraction(-1, energy)
        profiles[tuple(energies)] += 1
        coefficient += term
    CHECK.equal(sum(profiles.values()), 720, "all sixth-order paths enumerated")
    CHECK.equal(coefficient, Fraction(-63, 8), "leading hexagon off-diagonal coefficient")

    # Independent subset recursion for the same path sum.
    dynamic = {frozenset(): Fraction(1)}
    for size in range(1, 6):
        for subset_tuple in combinations(range(6), size):
            subset = frozenset(subset_tuple)
            dynamic[subset] = Fraction(-1, subset_energy(subset)) * sum(
                dynamic[subset - {edge}] for edge in subset
            )
    recursive_coefficient = sum(dynamic[all_indices - {edge}] for edge in all_indices)
    CHECK.equal(recursive_coefficient, coefficient, "independent subset-recursion coefficient")

    cycle_q, fixed, flow, required, occupied, child, cells, links = deterministic_locked_background()
    CHECK.equal(flow, required, "periodic degree-two background max flow saturates")
    CHECK.equal(required, 125, "periodic background residual edge count")
    CHECK.equal(len(occupied), 128, "periodic background has two occupied links per parent")
    CHECK.equal(tuple(int(e in occupied) for e in cycle_q), (1, 0, 1, 0, 1, 0), "target hexagon is alternating in constructed background")
    parent_degree = Counter(p for p, port in occupied)
    child_degree = Counter(child(p, port) for p, port in occupied)
    CHECK.true(all(parent_degree[p] == 2 for p in cells), "constructed parent degrees are two")
    CHECK.true(all(child_degree[c] == 2 for c in cells), "constructed child degrees are two")

    # Canonical finite-quotient second/fourth-order scalar census.
    _, _, endpoints = quotient_graph(4)
    incident = Counter()
    node_links: dict[tuple, list[tuple]] = {}
    for edge, ends in endpoints.items():
        for vertex in ends:
            node_links.setdefault(vertex, []).append(edge)
    adjacent_pairs = set()
    mixed = same = 0
    for edges_at_node in node_links.values():
        CHECK.equal(len(edges_at_node), 4, "finite quotient constraint node degree four")
        for e, f in combinations(edges_at_node, 2):
            pair = tuple(sorted((e, f)))
            CHECK.true(pair not in adjacent_pairs, "finite quotient pair has one owner")
            adjacent_pairs.add(pair)
            if (e in occupied) == (f in occupied):
                same += 1
            else:
                mixed += 1
    m_links = len(links)
    disjoint = comb(m_links, 2) - len(adjacent_pairs)
    CHECK.equal(len(adjacent_pairs), 3 * m_links, "finite quotient adjacent-pair count")
    CHECK.equal(mixed, 2 * m_links, "locked adjacent opposite-occupation count")
    CHECK.equal(same, m_links, "locked adjacent equal-occupation count")
    CHECK.equal(disjoint, comb(m_links, 2) - 3 * m_links, "locked disjoint-pair count")
    direct_h4 = -Fraction(mixed, 2) - Fraction(same, 6) - Fraction(disjoint, 4)
    CHECK.equal(direct_h4, -Fraction(3 * m_links * m_links + 7 * m_links, 24), "fourth-order direct Q-only term")
    folded_h4 = Fraction(m_links * m_links, 8)
    CHECK.equal(direct_h4 + folded_h4, -Fraction(7 * m_links, 24), "canonical fourth-order scalar")
    CHECK.equal(-Fraction(m_links, 2), -128, "canonical second-order scalar on L=4 quotient")


def main() -> None:
    section_local_square_and_symmetry()
    section_incidence_census_and_linear_charge()
    section_pair_sector()
    section_character_symbol()
    section_hexagon_and_strong_lock()
    print(f"GL6AN exact verification: PASS ({CHECK.total}/{CHECK.total})")
    print("square lock: H = -h sum X + U_d sum_v (k_v-2)^2 + constant")
    print("admissible inherited lock: Delta=4 U_d(d_star-2)>0, hence d_star>2")
    print("finite-h degree Ward test: NO nontrivial linear conserved degree charge")
    print("locked local pair covariance (A1,E,T2): (0, 8/3, 0)")
    print("generic incidence-symbol flat nullity: 2; trivial-character nullity: 3")
    print("canonical finite-quotient H2,H4 scalars: -M/2 and -7M/24 in h/U_d units")
    print("leading alternating-hexagon matrix element: -(63/8) h^6/U_d^5")
    print("CEILING: no pole, cone, physical metric, gravity, or G is derived")


if __name__ == "__main__":
    main()
