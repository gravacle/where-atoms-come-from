#!/usr/bin/env python3
"""Independent finite checks for GRA-FE-F3-Q4-DICJ-V001.

The verifier checks the new graph join and rederives the inherited symmetric
hexagon coefficient.  It does not simulate a thermodynamic U(1) phase or turn
the conditional Q4-CARRIER/EDGE-LIFT into a physical fact.
"""

from __future__ import annotations

from collections import deque
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb
from pathlib import Path


HERE = Path(__file__).resolve().parent
THEOREM = (HERE / "THEOREM.md").read_text()
CHECKS: list[str] = []


def check(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def add(v: tuple[int, ...], w: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(a + b for a, b in zip(v, w))


def sub(v: tuple[int, ...], w: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(a - b for a, b in zip(v, w))


E = tuple(tuple(1 if i == a else 0 for i in range(4)) for a in range(4))


def compositions(n: int, k: int = 4):
    if k == 1:
        yield (n,)
        return
    for first in range(n + 1):
        for rest in compositions(n - first, k - 1):
            yield (first,) + rest


def infinite_neighbors(v: tuple[int, ...], n: int) -> tuple[tuple[int, ...], ...]:
    total = sum(v)
    if total == n:
        return tuple(add(v, e) for e in E)
    if total == n + 1:
        return tuple(sub(v, e) for e in E)
    raise ValueError("vertex is outside the two-front slab")


def finite_neighbors(v: tuple[int, ...], n: int) -> tuple[tuple[int, ...], ...]:
    return tuple(w for w in infinite_neighbors(v, n) if min(w) >= 0)


def ball(neighbor, start: tuple[int, ...], radius: int):
    distance = {start: 0}
    queue = deque([start])
    while queue:
        v = queue.popleft()
        if distance[v] == radius:
            continue
        for w in neighbor(v):
            if w not in distance:
                distance[w] = distance[v] + 1
                queue.append(w)
    return distance


def direction(delta: tuple[int, ...], sign: int) -> int:
    target = sign
    hits = [i for i, x in enumerate(delta) if x == target]
    if len(hits) != 1 or any(x not in (0, target) for x in delta):
        raise AssertionError(f"not a signed append edge: {delta}")
    return hits[0]


def six_cycles_through_origin():
    start = (0, 0, 0, 0)
    n = 0
    found = set()

    def dfs(path):
        v = path[-1]
        if len(path) == 6:
            if start in infinite_neighbors(v, n):
                forward = tuple(path[1:])
                reverse = tuple(reversed(path[1:]))
                found.add(min(forward, reverse))
            return
        for w in infinite_neighbors(v, n):
            if w == start or w in path:
                continue
            dfs(path + [w])

    dfs([start])
    return [(start,) + cycle for cycle in found]


def periodic_neighbors(vertex, size: int):
    part, r = vertex
    shifts = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
    if part == 0:
        return tuple(
            (1, tuple((r[i] + shift[i]) % size for i in range(3)))
            for shift in shifts
        )
    return tuple(
        (0, tuple((r[i] - shift[i]) % size for i in range(3)))
        for shift in shifts
    )


def graph_girth(vertices, neighbors) -> int:
    best = 10**9
    for root in vertices:
        dist = {root: 0}
        parent = {root: None}
        queue = deque([root])
        while queue:
            v = queue.popleft()
            if 2 * dist[v] + 1 >= best:
                continue
            for w in neighbors(v):
                if w not in dist:
                    dist[w] = dist[v] + 1
                    parent[w] = v
                    queue.append(w)
                elif parent[v] != w:
                    best = min(best, dist[v] + dist[w] + 1)
    return best


def symmetric_hexagon_coefficient() -> Fraction:
    # Cycle edge i joins cycle vertices i and i+1.  Occupations alternate,
    # with even edges occupied, so sigma=1-2n is -1,+1,-1,+1,-1,+1.
    sigma = tuple(-1 if i % 2 == 0 else 1 for i in range(6))

    def gap(selected) -> int:
        selected = set(selected)
        vertex_charge = []
        for vertex in range(6):
            incident = ((vertex - 1) % 6, vertex)
            vertex_charge.append(sum(sigma[e] for e in incident if e in selected))
        return sum(q * q for q in vertex_charge)

    total = Fraction(0, 1)
    denominator_classes = {}
    for order in permutations(range(6)):
        gaps = tuple(gap(order[:r]) for r in range(1, 6))
        check_nonzero = all(g > 0 for g in gaps)
        if not check_nonzero:
            raise AssertionError(f"proper alternating prefix returned to ice: {order}")
        product_gap = 1
        for g in gaps:
            product_gap *= g
        total += Fraction(1, product_gap)
        denominator_classes[tuple(sorted(gaps))] = denominator_classes.get(tuple(sorted(gaps)), 0) + 1

    check(sum(denominator_classes.values()) == 720, "all 720 hexagon orders retained")
    check(len(denominator_classes) == 5, "five symmetric denominator classes")
    return total


# 1. Tetrahedral frame and FCC lattice.
tetra = (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)
check(tuple(sum(v[i] for v in tetra) for i in range(3)) == (0, 0, 0),
      "tetrahedral bond vectors sum to zero")
gram = tuple(tuple(sum(a * b for a, b in zip(v, w)) for w in tetra) for v in tetra)
check(all(gram[i][j] == (3 if i == j else -1) for i in range(4) for j in range(4)),
      "tetrahedral Gram is diag 1 and offdiag -1/3")

fcc_generators = ((0, 1, 1), (1, 0, 1), (1, 1, 0))
generated = {
    tuple(sum(coeff[j] * fcc_generators[j][i] for j in range(3)) for i in range(3))
    for coeff in product(range(-3, 4), repeat=3)
}
check(all(sum(v) % 2 == 0 for v in generated), "A3 roots land in the FCC parity lattice")
for v in product(range(-2, 3), repeat=3):
    if sum(v) % 2 == 0:
        coefficients = (
            (v[1] + v[2] - v[0]) // 2,
            (v[0] + v[2] - v[1]) // 2,
            (v[0] + v[1] - v[2]) // 2,
        )
        rebuilt = tuple(
            sum(coefficients[j] * fcc_generators[j][i] for j in range(3))
            for i in range(3)
        )
        if rebuilt != v:
            raise AssertionError((v, coefficients, rebuilt))
check(True, "FCC parity lattice is generated by the three A3 roots")

# 2. Finite slab incidence and boundary-degree formula.
n = 8
s_n = tuple(compositions(n))
s_np1 = tuple(compositions(n + 1))
check(len(s_n) == comb(n + 3, 3) and len(s_np1) == comb(n + 4, 3),
      "finite q4 front cardinalities")
check(all(len(finite_neighbors(m, n)) == 4 for m in s_n),
      "every finite-slab parent has four children")
check(all(len(finite_neighbors(c, n)) == sum(x > 0 for x in c) for c in s_np1),
      "child degree equals positive-coordinate count")
check(all(len(finite_neighbors(c, n)) == 4 for c in s_np1 if min(c) > 0),
      "deep-interior children have degree four")

children = {m: set(finite_neighbors(m, n)) for m in s_n}
check(all(len(children[a] & children[b]) <= 1 for a, b in combinations(s_n, 2)),
      "distinct parents share at most one child, excluding four-cycles")

# 3. Exact six-cycle witness and classification.
witness = (
    (0, 0, 0, 0),
    (1, 0, 0, 0),
    (1, -1, 0, 0),
    (1, -1, 1, 0),
    (0, -1, 1, 0),
    (0, 0, 1, 0),
)
check(all(witness[(i + 1) % 6] in infinite_neighbors(witness[i], 0) for i in range(6)),
      "explicit simple six-cycle")
cycles = six_cycles_through_origin()
check(len(cycles) == 12, "twelve undirected diamond hexagons through one vertex")
for cycle in cycles:
    plus = []
    minus = []
    closed = cycle + (cycle[0],)
    for i in range(6):
        delta = sub(closed[i + 1], closed[i])
        if i % 2 == 0:
            plus.append(direction(delta, +1))
        else:
            minus.append(direction(delta, -1))
    if sorted(plus) != sorted(minus) or len(set(plus)) != 3:
        raise AssertionError((cycle, plus, minus))
check(True, "every simple six-cycle uses three distinct append labels")

# 4. Rooted local exhaustion.
for radius in range(1, 6):
    base = (radius, radius, radius, radius)
    layer = 4 * radius
    infinite_ball = ball(lambda v: infinite_neighbors(v, layer), base, radius)
    finite_ball = ball(lambda v: finite_neighbors(v, layer), base, radius)
    if infinite_ball != finite_ball:
        raise AssertionError(f"local exhaustion failed at radius {radius}")
check(True, "nonnegative q4 slabs locally exhaust infinite diamond")

# 5. A safe periodic quotient is degree four and retains girth six.
size = 4
periodic_vertices = tuple(
    (part, r) for part in (0, 1) for r in product(range(size), repeat=3)
)
check(len(periodic_vertices) == 2 * size**3, "periodic quotient has 2L^3 vertices")
check(all(len(set(periodic_neighbors(v, size))) == 4 for v in periodic_vertices),
      "periodic quotient is simple and degree four")
check(graph_girth(periodic_vertices, lambda v: periodic_neighbors(v, size)) == 6,
      "L=4 periodic quotient retains girth six")

# 6. The d*=2 arrow rule is exactly two-in/two-out on either bipartition.
for occupied in combinations(range(4), 2):
    occupied = set(occupied)
    a_out = sum(i in occupied for i in range(4))
    a_in = 4 - a_out
    b_in = a_out
    b_out = 4 - b_in
    if (a_in, a_out, b_in, b_out) != (2, 2, 2, 2):
        raise AssertionError(occupied)
check(True, "d*=2 occupation is exactly the diamond two-in/two-out ice rule")

# 7. Independent inherited coefficient check.
j6 = symmetric_hexagon_coefficient()
check(j6 == Fraction(63, 8), "symmetric alternating-hexagon coefficient is 63/8")

# 8. Claim-ceiling guards.
check(
    sum(line.strip() == r"\[" for line in THEOREM.splitlines())
    == sum(line.strip() == r"\]" for line in THEOREM.splitlines()),
    "display-math delimiters are balanced",
)
required_phrases = (
    "Q4-CARRIER/EDGE-LIFT",
    "construction target, not a theorem",
    "support/phase join, not gravity closure",
    "Physical coexistence",
    "Support stability",
    "Thermodynamic and all-orders control",
    "Visible electromagnetism",
    "Tensor gravity",
)
check(all(phrase in THEOREM for phrase in required_phrases),
      "all physical and scientific ceilings are explicit")

print(f"Q4DICJ verification: PASS ({len(CHECKS)}/{len(CHECKS)})")
for index, label in enumerate(CHECKS, 1):
    print(f"{index:02d}. PASS  {label}")
