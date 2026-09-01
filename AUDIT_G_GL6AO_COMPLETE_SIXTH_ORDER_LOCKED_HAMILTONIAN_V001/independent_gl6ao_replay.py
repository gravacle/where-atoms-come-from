#!/usr/bin/env python3
"""Independent exact GL6AO hostile replay; imports no author verifier."""

from __future__ import annotations

import hashlib
import itertools
from collections import Counter, defaultdict, deque
from fractions import Fraction
from functools import lru_cache
from math import comb
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def mm(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def chain(*matrices):
    out = matrices[0]
    for matrix in matrices[1:]:
        out = mm(out, matrix)
    return out


def madd(*matrices):
    return [[sum(matrix[i][j] for matrix in matrices)
             for j in range(len(matrices[0][0]))]
            for i in range(len(matrices[0]))]


def mscale(matrix, scalar):
    return [[scalar * value for value in row] for row in matrix]


def scalar_block(matrix):
    return matrix[0][0]


def eigenseries(energies, perturbation, order=6):
    """Intermediate-normalized exact eigenpair series for a simple P state."""
    size = len(energies)
    psi = [[Fraction(0) for _ in range(size)] for _ in range(order + 1)]
    psi[0][0] = Fraction(1)
    coeff = [Fraction(0) for _ in range(order + 1)]
    for n in range(1, order + 1):
        v_previous = [sum(perturbation[i][j] * psi[n - 1][j]
                          for j in range(size)) for i in range(size)]
        coeff[n] = v_previous[0]
        for i in range(1, size):
            folded = sum(coeff[k] * psi[n - k][i] for k in range(1, n + 1))
            psi[n][i] = (folded - v_previous[i]) / energies[i]
    return coeff


def section_kato_formula() -> None:
    grades = (0, 1, 0, 1, 0, 1)
    energies = tuple(Fraction(value) for value in (0, 2, 3, 5, 7, 11))
    size = len(energies)
    p = [[Fraction(int(i == j == 0)) for j in range(size)] for i in range(size)]
    r = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    for i in range(1, size):
        r[i][i] = -Fraction(1, energies[i])

    for seed in (1, 3, 7):
        w = [[Fraction(0) for _ in range(size)] for _ in range(size)]
        for i in range(size):
            for j in range(i + 1, size):
                if grades[i] != grades[j]:
                    value = Fraction(((i + 2) * (j + 3) + seed) % 7 + 1)
                    if (i + j + seed) % 2:
                        value = -value
                    w[i][j] = w[j][i] = value

        b = scalar_block(chain(p, w, r, w, p))
        a2 = scalar_block(chain(p, w, mm(r, r), w, p))
        a3 = scalar_block(chain(p, w, chain(r, r, r), w, p))
        t4 = scalar_block(chain(p, w, r, w, r, w, r, w, p))
        d = t4 - a2 * b
        t6 = scalar_block(chain(p, w, r, w, r, w, r, w, r, w, r, w, p))
        x4 = sum((
            scalar_block(chain(p, w, mm(r, r), w, r, w, r, w, p)),
            scalar_block(chain(p, w, r, w, mm(r, r), w, r, w, p)),
            scalar_block(chain(p, w, r, w, r, w, mm(r, r), w, p)),
        ), Fraction(0))
        k6 = t6 - b * x4 + b * b * a3 - d * a2
        exact = eigenseries(energies, w)
        check(exact[1] == exact[3] == exact[5] == 0,
              "grading kills odd eigenseries blocks")
        check(exact[2] == b, "independent eigenseries reproduces K2")
        check(exact[4] == d, "independent eigenseries reproduces K4 fold")
        check(exact[6] == k6, "independent eigenseries reproduces K6 formula")


def quotient_graph(length=4):
    cells = tuple(itertools.product(range(length), repeat=3))

    def child(cell, port):
        out = list(cell)
        if port < 3:
            out[port] = (out[port] + 1) % length
        return tuple(out)

    edges = tuple((cell, port) for cell in cells for port in range(4))
    endpoints = {edge: (("P", edge[0]), ("C", child(*edge))) for edge in edges}
    adjacency = defaultdict(list)
    incident = defaultdict(list)
    for edge, (left, right) in endpoints.items():
        adjacency[left].append((right, edge))
        adjacency[right].append((left, edge))
        incident[left].append(edge)
        incident[right].append(edge)
    return cells, edges, endpoints, adjacency, incident


def enumerate_cycles(adjacency):
    cycles = set()
    for start in adjacency:
        def visit(node, visited, path):
            if len(path) == 6:
                if node == start:
                    cycles.add(frozenset(path))
                return
            for other, edge in adjacency[node]:
                if edge in path:
                    continue
                if other == start:
                    if len(path) == 5:
                        visit(other, visited, path + (edge,))
                elif other not in visited:
                    visit(other, visited | {other}, path + (edge,))
        visit(start, {start}, tuple())
    return cycles


def section_graph() -> tuple:
    cells, edges, endpoints, adjacency, incident = quotient_graph()
    nodes = set(adjacency)
    check(len(cells) == 64, "period-four cell count")
    check(len(nodes) == 128, "constraint-node count")
    check(len(edges) == 256, "active-link count")
    check(all(len(incident[node]) == 4 for node in nodes), "degree four")
    check(len(set(endpoints.values())) == len(edges), "simple incidence edges")

    reached = {next(iter(nodes))}
    queue = deque(reached)
    while queue:
        node = queue.popleft()
        for other, _ in adjacency[node]:
            if other not in reached:
                reached.add(other)
                queue.append(other)
    check(reached == nodes, "Q4 is connected")

    girth = 10**9
    for start in nodes:
        distance = {start: 0}
        parent = {start: None}
        queue = deque([start])
        while queue:
            node = queue.popleft()
            for other, _ in adjacency[node]:
                if other not in distance:
                    distance[other] = distance[node] + 1
                    parent[other] = node
                    queue.append(other)
                elif parent[node] != other:
                    girth = min(girth, distance[node] + distance[other] + 1)
    check(girth == 6, "Q4 girth is exactly six")

    cycles = enumerate_cycles(adjacency)
    check(len(cycles) == 256, "independent Q4 six-cycle count")
    core_incidence = Counter(edge for cycle in cycles for edge in cycle)
    check(all(core_incidence[edge] == 6 for edge in edges), "six cycle cores per link")

    support_incidence = Counter()
    for cycle in cycles:
        vertices = {vertex for edge in cycle for vertex in endpoints[edge]}
        check(len(vertices) == 6, "six vertices per simple hexagon")
        support = {edge for vertex in vertices for edge in incident[vertex]}
        check(len(support) == 18, "local projected cycle operator has 18-link support")
        for edge in support:
            support_incidence[edge] += 1
    check(all(support_incidence[edge] == 18 for edge in edges),
          "each link occurs in exactly 18 projected-cycle supports")

    # A length-six quotient cycle has three parent-to-parent steps.  No integer
    # coordinate displacement can reach nonzero 4Z in only three root steps.
    for steps in itertools.product(
            tuple(tuple(int(i == a) - int(i == b) for i in range(3))
                  for a in range(4) for b in range(4) if a != b), repeat=3):
        displacement = tuple(sum(step[i] for step in steps) for i in range(3))
        if all(value % 4 == 0 for value in displacement):
            check(displacement == (0, 0, 0), "period-four six-cycle closure lifts exactly")

    return edges, endpoints, incident, cycles


def section_shape_and_occupation_census(edges, incident):
    m = len(edges)
    adjacent = {frozenset((a, b)) for row in incident.values()
                for a, b in itertools.combinations(row, 2)}
    check(len(adjacent) == 3 * m, "adjacent edge-pair count")

    shapes = Counter()
    for a, b, c in itertools.combinations(edges, 3):
        count = sum(frozenset(pair) in adjacent for pair in ((a, b), (a, c), (b, c)))
        shapes[count] += 1
    expected = Counter({
        0: comb(m, 3) - 3 * m * m + 19 * m,
        1: 3 * m * (m - 10),
        2: 9 * m,
        3: 2 * m,
    })
    check(shapes == expected, "all unordered three-edge shapes independently enumerated")

    # Exhaust the six possible local degree-two patterns.  These identities
    # prove the global occupation counts for every locked configuration.
    for occupied in itertools.combinations(range(4), 2):
        occupied = set(occupied)
        pairs = Counter(
            "opposite" if ((a in occupied) != (b in occupied)) else "equal"
            for a, b in itertools.combinations(range(4), 2)
        )
        check(pairs == Counter({"opposite": 4, "equal": 2}),
              "local degree-two adjacent pair split")
        for middle in range(4):
            continuations = Counter(
                "opposite" if ((middle in occupied) != (other in occupied)) else "equal"
                for other in range(4) if other != middle
            )
            check(continuations == Counter({"opposite": 2, "equal": 1}),
                  "two opposite and one equal continuation per endpoint")
        stars = Counter()
        for triple in itertools.combinations(range(4), 3):
            signs = [1 if port not in occupied else -1 for port in triple]
            pair_energies = tuple(sorted(2 if signs[i] != signs[j] else 6
                                         for i, j in itertools.combinations(range(3), 2)))
            triple_energy = 3 + sum(signs) ** 2
            stars[(pair_energies, triple_energy)] += 1
        check(stars == Counter({((2, 2, 6), 4): 4}), "local star class")


def word_sum(multiplicities, pair_energy, triple_energy=None, powers=None):
    total_length = sum(multiplicities)
    if powers is None:
        powers = (1,) * (total_length - 1)

    @lru_cache(None)
    def recurse(counts, parity_tuple, step):
        if sum(counts) == 1:
            return Fraction(1)
        parity = set(parity_tuple)
        answer = Fraction(0)
        for edge, count in enumerate(counts):
            if count == 0:
                continue
            next_counts = list(counts)
            next_counts[edge] -= 1
            next_parity = set(parity)
            if edge in next_parity:
                next_parity.remove(edge)
            else:
                next_parity.add(edge)
            if not next_parity:
                continue
            if len(next_parity) == 1:
                energy = 2
            elif len(next_parity) == 2:
                energy = pair_energy[tuple(sorted(next_parity))]
            elif len(next_parity) == 3:
                energy = triple_energy
            else:
                raise AssertionError("unexpected parity support")
            exponent = powers[step]
            factor = Fraction((-1) ** exponent, energy ** exponent)
            answer += factor * recurse(tuple(next_counts), tuple(sorted(next_parity)), step + 1)
        return answer

    return recurse(tuple(multiplicities), tuple(), 0)


def padd(a, b):
    length = max(len(a), len(b))
    return tuple((a[i] if i < len(a) else Fraction(0))
                 + (b[i] if i < len(b) else Fraction(0)) for i in range(length))


def pscale(a, scalar):
    return tuple(Fraction(scalar) * value for value in a)


def pmul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        for j, right in enumerate(b):
            out[i + j] += left * right
    return tuple(out)


def section_words_and_cancellation():
    repeated = {}
    folded = {}
    fourth = {}
    for energy in (2, 4, 6):
        pair = {(0, 1): energy}
        fourth[energy] = word_sum((2, 2), pair)
        repeated[energy] = word_sum((4, 2), pair) + word_sum((2, 4), pair)
        folded[energy] = sum((word_sum((2, 2), pair, powers=powers)
                              for powers in ((2, 1, 1), (1, 2, 1), (1, 1, 2))),
                             Fraction(0))
    check(fourth == {2: Fraction(-1, 2), 4: Fraction(-1, 4), 6: Fraction(-1, 6)},
          "independent fourth-order pair weights")
    check(repeated == {2: Fraction(-1, 4), 4: Fraction(-1, 16), 6: Fraction(-1, 36)},
          "independent repeated-pair sixth-order weights")
    check(folded == {2: Fraction(3, 4), 4: Fraction(5, 16), 6: Fraction(7, 36)},
          "independent squared-resolvent fold weights")

    classes = (
        ((4, 4, 4), 6, Fraction(-9, 32)),
        ((2, 4, 4), 4, Fraction(-9, 16)),
        ((4, 4, 6), 8, Fraction(-29, 144)),
        ((2, 2, 6), 4, Fraction(-109, 144)),
        ((2, 2, 4), 2, Fraction(-41, 32)),
        ((2, 4, 6), 6, Fraction(-337, 864)),
        ((4, 6, 6), 10, Fraction(-209, 1440)),
    )
    triple_weights = []
    for energies, triple_energy, expected in classes:
        pair = {(0, 1): energies[0], (0, 2): energies[1], (1, 2): energies[2]}
        value = word_sum((2, 2, 2), pair, triple_energy)
        check(value == expected, "independent three-link sixth-order word weight")
        triple_weights.append(value)

    m = (Fraction(0), Fraction(1))
    m2 = pmul(m, m)
    choose2 = (Fraction(0), Fraction(-1, 2), Fraction(1, 2))
    choose3 = (Fraction(0), Fraction(1, 3), Fraction(-1, 2), Fraction(1, 6))
    counts = (
        padd(padd(choose3, pscale(m2, -3)), pscale(m, 19)),
        padd(pscale(m2, 2), pscale(m, -20)),
        padd(m2, pscale(m, -10)),
        pscale(m, 2), pscale(m, 4), pscale(m, 4), m,
    )
    t6_triples = tuple()
    for count, weight in zip(counts, triple_weights):
        t6_triples = padd(t6_triples, pscale(count, weight))
    check(t6_triples == (Fraction(0), Fraction(-2237, 4320),
                          Fraction(-197, 576), Fraction(-3, 64)),
          "independent direct triple polynomial")

    t6_pairs = padd(
        pscale(m, 2 * repeated[2] + repeated[6] - 3 * repeated[4]),
        pscale(choose2, repeated[4]),
    )
    check(t6_pairs == (Fraction(0), Fraction(-89, 288), Fraction(-1, 32)),
          "independent repeated-pair polynomial")
    t6 = padd(t6_triples, t6_pairs)
    check(t6 == (Fraction(0), Fraction(-893, 1080), Fraction(-215, 576),
                 Fraction(-3, 64)), "complete direct T6 diagonal")

    x4 = padd(
        pscale(m, 2 * folded[2] + folded[6] - 3 * folded[4]),
        pscale(choose2, folded[4]),
    )
    check(x4 == (Fraction(0), Fraction(173, 288), Fraction(5, 32)),
          "complete X4 fold polynomial")
    b = pscale(m, Fraction(-1, 2))
    d = pscale(m, Fraction(-7, 24))
    a2 = pscale(m, Fraction(1, 4))
    a3 = pscale(m, Fraction(-1, 8))
    k6 = padd(padd(padd(t6, pscale(pmul(b, x4), -1)), pmul(pmul(b, b), a3)),
                pscale(pmul(d, a2), -1))
    check(k6 == (Fraction(0), Fraction(-893, 1080), Fraction(0), Fraction(0)),
          "M cubed and squared cancel exactly")


def hex_energy(mask):
    signs = (1, -1, 1, -1, 1, -1)
    charges = [0] * 6
    for edge in range(6):
        if mask & (1 << edge):
            charges[edge] += signs[edge]
            charges[(edge + 1) % 6] += signs[edge]
    return sum(value * value for value in charges)


def section_hexagon():
    dynamic = {0: Fraction(1)}
    for size in range(1, 6):
        for mask in range(1, 64):
            if bin(mask).count("1") != size:
                continue
            energy = hex_energy(mask)
            check(energy > 0, "no proper alternating-hexagon subset returns to lock")
            dynamic[mask] = Fraction(-1, energy) * sum(
                (dynamic[mask ^ (1 << edge)] for edge in range(6) if mask & (1 << edge)),
                Fraction(0),
            )
    coefficient = sum((dynamic[63 ^ (1 << edge)] for edge in range(6)), Fraction(0))
    check(coefficient == Fraction(-63, 8), "independent subset recursion gives hexagon coefficient")


def section_custody_and_scope():
    rows = []
    for raw in (HERE / "AUDITED_TARGETS.sha256").read_text().splitlines():
        if not raw.strip():
            continue
        expected, relative = raw.split(maxsplit=1)
        rows.append(relative)
        target = ROOT / relative
        check(target.is_file(), f"audited target exists: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"audited target hash: {relative}")
    check(len(rows) == 11 and len(set(rows)) == 11, "exact unique author target count")

    theorem = (ROOT / "LANE_CROSS_RFT_GRA_GL6AO_COMPLETE_SIXTH_ORDER_LOCKED_HAMILTONIAN_V001/THEOREM.md").read_text()
    for marker in (
        "K_6 = T_6 - b X_4 + b^2 A_3 - d A_2",
        "-(893/1080) M",
        "<s'|K_6|s> = -63/8",
        "formal sixth-order interaction",
        "does **not** prove",
        "gravity, or `G`",
    ):
        check(marker in theorem, f"author theorem scope marker: {marker}")


def main():
    section_kato_formula()
    edges, endpoints, incident, cycles = section_graph()
    section_shape_and_occupation_census(edges, incident)
    section_words_and_cancellation()
    section_hexagon()
    section_custody_and_scope()
    print(f"PASS__INDEPENDENT_GL6AO_HOSTILE_REPLAY__{checks}/{checks}")
    print("KATO=K6_T6_MINUS_bX4_PLUS_b2A3_MINUS_dA2")
    print("Q4=M256_GIRTH6_HEX256_CORE6_SUPPORT18")
    print("DIAGONAL=COMMON_MINUS893_OVER1080_M_AFTER_M3_M2_CANCELLATION")
    print("OFFDIAGONAL=ONLY_ALTERNATING_HEXAGON_MINUS63_OVER8")
    print("LINKED=FORMAL_FINITE_RANGE_NOT_ALL_ORDERS_OR_PHASE")
    print("CEILING=NO_POLE_MOMENTUM_CONE_RICCI_GRAVITY_G")


if __name__ == "__main__":
    main()
