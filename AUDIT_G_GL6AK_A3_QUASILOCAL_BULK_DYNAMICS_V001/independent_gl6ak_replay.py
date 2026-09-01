#!/usr/bin/env python3
"""Independent hostile replay for the mutable GL6AK pre-freeze theorem.

This reconstruction deliberately imports neither the author ledger nor the
author verifier.  It checks the finite A3 incidence, collar embedding, shell
census, influence-matrix distance filtration, integrated factorial tail,
Folner property, and the multiplicity-free six-pair projectors.
"""

from __future__ import annotations

import hashlib
import itertools
import math
from collections import defaultdict, deque
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001"
CHECKS = 0


def check(value: bool, message: str) -> None:
    global CHECKS
    if not value:
        raise AssertionError(message)
    CHECKS += 1


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def sub(x, y):
    return tuple(a - b for a, b in zip(x, y))


ZERO = (0, 0, 0, 0)
PORTS = range(4)
E = tuple(tuple(int(i == a) for i in PORTS) for a in PORTS)
PAIRS = tuple(itertools.combinations(PORTS, 2))


def radius(x):
    check(sum(x) == 0, "radius input outside A3")
    return sum(abs(v) for v in x) // 2


def cell_ball(r):
    out = set()
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            for c in range(-r, r + 1):
                x = (a, b, c, -a - b - c)
                if radius(x) <= r:
                    out.add(x)
    return out


def internal(x, a, b):
    return frozenset(((x, a), (x, b)))


def shared(x, a, b):
    return frozenset(((x, a), (add(x, sub(E[a], E[b])), b)))


def children(site):
    x, a = site
    return add(x, E[a])


def all_edges_anchored(cells):
    edges = set()
    for x in cells:
        for a, b in PAIRS:
            edges.add(internal(x, a, b))
            edges.add(shared(x, a, b))
    return edges


def permute_x(x, sigma):
    y = [0] * 4
    for old, new in enumerate(sigma):
        y[new] = x[old]
    return tuple(y)


def permute_site(site, sigma):
    x, a = site
    return permute_x(x, sigma), sigma[a]


# 1. Canonical edge representatives and exact degree-six incidence.
cells = cell_ball(3)
edges = all_edges_anchored(cells)
for x in cell_ball(2):
    for a in PORTS:
        site = (x, a)
        incident = {edge for edge in edges if site in edge}
        check(len(incident) == 6, "bulk active-link degree is not six")
        same_parent = sum(len({p[0] for p in edge}) == 1 for edge in incident)
        same_child = sum(children(tuple(edge)[0]) == children(tuple(edge)[1]) for edge in incident)
        check((same_parent, same_child) == (3, 3), "degree does not split 3+3")
    touching = {edge for edge in edges if any(site[0] == x for site in edge)}
    check(len(touching) == 18, "cell does not touch 6 internal plus 12 shared terms")

# Each a>b shared partner is the same edge as the re-anchored b<a representative.
for x in cell_ball(2):
    for a in PORTS:
        for b in PORTS:
            if a <= b:
                continue
            y = add(x, sub(E[a], E[b]))
            check(
                frozenset(((x, a), (y, b))) == shared(y, b, a),
                "shared representative re-anchoring failed",
            )

# 2. Translation and S4 covariance, including orientation re-anchoring.
translations = (ZERO, (1, -1, 0, 0), (2, -1, -1, 0))
for x in cell_ball(1):
    for a, b in PAIRS:
        for z in translations:
            for edge in (internal(x, a, b), shared(x, a, b)):
                moved = frozenset((add(site[0], z), site[1]) for site in edge)
                p, q = tuple(moved)
                check(
                    (p[0] == q[0] and p[1] != q[1])
                    or (children(p) == children(q) and p != q),
                    "translation did not preserve either edge relation",
                )
        for sigma in itertools.permutations(PORTS):
            for edge in (internal(x, a, b), shared(x, a, b)):
                moved = frozenset(permute_site(site, sigma) for site in edge)
                p, q = tuple(moved)
                check(
                    (p[0] == q[0] and p[1] != q[1])
                    or (children(p) == children(q) and p != q),
                    "S4 did not preserve either edge relation",
                )

# 3. One common finite FPSS slab embeds a collared finite patch exactly.
for r in range(4):
    patch = cell_ball(r)
    collar = patch | {add(x, sub(E[a], E[b])) for x in patch for a in PORTS for b in PORTS if a != b}
    mins = tuple(min(x[i] for x in collar) for i in PORTS)
    m = tuple(max(1, 1 - mins[i]) for i in PORTS)
    n = sum(m)
    for x in collar:
        mx = add(m, x)
        check(sum(mx) == n and min(mx) >= 1, "collared parent is not strict interior")
    for x in patch:
        for a, b in PAIRS:
            y = add(x, sub(E[a], E[b]))
            check(
                add(add(m, x), E[a]) == add(add(m, y), E[b]),
                "translated shared-child equality failed",
            )

# 4. Exact shell assignment and the coarse 18(2r+1)^3 ceiling.
outer = cell_ball(7)
outer_edges = all_edges_anchored(outer)
assigned = defaultdict(list)
for edge in outer_edges:
    endpoint_cells = sorted({site[0] for site in edge})
    if any(radius(x) > 6 for x in endpoint_cells):
        continue
    owner = min(endpoint_cells, key=lambda x: (radius(x), x))
    assigned[owner].append(edge)
for owner, owned in assigned.items():
    check(len(owned) <= 18, "minimum-radius cell owns more than its 18 touching terms")
for r in range(6):
    nr = sum(len(owned) for owner, owned in assigned.items() if radius(owner) == r)
    shell = len(cell_ball(r) - (cell_ball(r - 1) if r else set()))
    check(nr <= 18 * shell, "interaction shell exceeds 18 times cell shell")
    check(shell <= (2 * r + 1) ** 3, "cell shell exceeds cubic bound")

# 5. Link-distance descent and the dressing-complete matrix filtration.
graph_cells = cell_ball(4)
graph_edges = all_edges_anchored(graph_cells)
adj = defaultdict(set)
for edge in graph_edges:
    p, q = tuple(edge)
    adj[p].add(q)
    adj[q].add(p)
for source_port in PORTS:
    source = (ZERO, source_port)
    distance = {source: 0}
    queue = deque([source])
    while queue:
        p = queue.popleft()
        for q in adj[p]:
            if q not in distance:
                distance[q] = distance[p] + 1
                queue.append(q)
    for target, dl in distance.items():
        if radius(target[0]) <= 3:
            check(dl >= radius(target[0]), "link distance fell below A3 cell distance")

# Powers of the independent influence matrix cannot connect before graph distance.
finite_sites = sorted((x, a) for x in cell_ball(1) for a in PORTS)
index = {p: i for i, p in enumerate(finite_sites)}
nsite = len(finite_sites)
neighbors = [set() for _ in finite_sites]
for edge in all_edges_anchored(cell_ball(2)):
    p, q = tuple(edge)
    if p in index and q in index:
        neighbors[index[p]].add(index[q])
        neighbors[index[q]].add(index[p])
K = [[0 for _ in range(nsite)] for _ in range(nsite)]
for i, ns in enumerate(neighbors):
    K[i][i] = len(ns)
    for j in ns:
        K[i][j] = 1
    check(sum(K[i]) <= 12, "influence row mass exceeds 2 Delta_L")

def multiply(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(nsite)) for j in range(nsite)] for i in range(nsite)]

power = [[int(i == j) for j in range(nsite)] for i in range(nsite)]
for order in range(5):
    for i, source in enumerate(finite_sites):
        # Independent BFS in the restricted graph.
        d = {i: 0}
        queue = deque([i])
        while queue:
            u = queue.popleft()
            for v in neighbors[u]:
                if v not in d:
                    d[v] = d[u] + 1
                    queue.append(v)
        for j, dj in d.items():
            if order < dj:
                check(power[i][j] == 0, "diagonal dressing advanced the endpoint")
    power = multiply(power, K)

# 6. Exact prefactor reduction and integrated factorial-tail identity.
for J in (Fraction(1, 7), Fraction(5, 3), Fraction(11, 2)):
    hbar = Fraction(13, 5)
    lam = 4 * J * 6 / hbar
    check(lam == 24 * J / hbar, "lambda reduction failed")
    check(Fraction(72) * J / hbar / lam == 3, "boundary factor is not three")

def tail(d, x):
    if x == 0:
        return 1.0 if d == 0 else 0.0
    term = math.exp(d * math.log(x) - math.lgamma(d + 1))
    total = term
    for k in range(d + 1, 220):
        term *= x / k
        total += term
        if term <= max(1e-300, total * 1e-17):
            break
    return total

for lam in (0.2, 1.3, 4.7):
    for d in range(8):
        # Simpson quadrature checks integral T_d(lambda u)du=T_(d+1)(lambda t)/lambda.
        steps = 20000
        total = 0.0
        for j in range(steps + 1):
            u = j / steps
            weight = 1 if j in (0, steps) else (4 if j % 2 else 2)
            total += weight * tail(d, lam * u)
        got = total / (3 * steps)
        want = tail(d + 1, lam) / lam
        check(abs(got - want) < 2e-11, "integrated tail identity failed")

for x in (0.5, 2.0, 6.0):
    values = []
    for R in (8, 16, 32, 64):
        values.append(sum((2 * r + 1) ** 3 * tail(r - 2 + 1, x) for r in range(R, 180)))
    check(all(b < a for a, b in zip(values, values[1:])), "boundary tail is not decreasing")
    check(values[-1] < 1e-30, "factorial tail does not beat cubic shell")

# 7. The displayed Z3 cubes are a concrete Folner sequence.
def folner(r):
    return {(a, b, c, -a - b - c) for a in range(-r, r + 1) for b in range(-r, r + 1) for c in range(-r, r + 1)}

for z in ((1, -1, 0, 0), (2, -1, -1, 0), (-3, 1, 1, 1)):
    ratios = []
    for r in (4, 12, 36):
        f = folner(r)
        shifted = {add(x, z) for x in f}
        ratios.append(len(f.symmetric_difference(shifted)) / len(f))
    check(ratios[2] < ratios[1] < ratios[0], "Folner boundary ratio did not shrink")

# 8. Independent exact A1 + E + T2 projectors in pair space.
def mm(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(6)) for j in range(6)] for i in range(6)]

I6 = [[Fraction(int(i == j)) for j in range(6)] for i in range(6)]
A = [[Fraction(int(i != j and bool(set(PAIRS[i]) & set(PAIRS[j])))) for j in range(6)] for i in range(6)]
P1 = [[Fraction(1, 6) for _ in range(6)] for _ in range(6)]
PE = [[v / 12 for v in row] for row in mm(A, [[A[i][j] - 4 * I6[i][j] for j in range(6)] for i in range(6)])]
PT = [[-v / 8 for v in row] for row in mm(
    [[A[i][j] - 4 * I6[i][j] for j in range(6)] for i in range(6)],
    [[A[i][j] + 2 * I6[i][j] for j in range(6)] for i in range(6)],
)]
Z6 = [[Fraction(0) for _ in range(6)] for _ in range(6)]
for P, rank in ((P1, 1), (PE, 2), (PT, 3)):
    check(mm(P, P) == P, "pair-space projector is not idempotent")
    check(sum(P[i][i] for i in range(6)) == rank, "pair-space projector rank failed")
for P, Q in ((P1, PE), (P1, PT), (PE, PT)):
    check(mm(P, Q) == Z6, "pair-space projectors are not orthogonal")
check(
    [[P1[i][j] + PE[i][j] + PT[i][j] for j in range(6)] for i in range(6)] == I6,
    "A1+E+T2 does not resolve pair space",
)

# 9. Reviewed-byte pin and strict claim-ceiling inspection.
theorem = TARGET / "THEOREM.md"
check(
    hashlib.sha256(theorem.read_bytes()).hexdigest()
    == "38cb58ef9fc52e1252e0b0d3415c54488c0471c0ddf25a35b7adf5aba41bccc9",
    "reviewed mutable theorem changed during hostile replay",
)
text = " ".join(theorem.read_text().lower().replace("**", "").split())
for token in (
    "does not turn the infinite net into one infinite record",
    "arbitrarily strong or new boundary laws are outside",
    "not poincaré covariance",
    "not uniqueness, purity, a ground-state or kms property",
    "it is not physical momentum",
    "does not say that any sector is gapless",
    "without importing a graviton, pole, metric, or ricci ansatz",
):
    check(token in text, f"missing strict ceiling: {token}")

print(f"PASS {CHECKS}/{CHECKS}")
