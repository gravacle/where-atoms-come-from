#!/usr/bin/env python3
"""Independent finite-census checks for GL6AI's exact algebra/topology."""

from __future__ import annotations

from collections import deque
from fractions import Fraction
from itertools import combinations, product
from math import exp, factorial


def compositions(n: int):
    return [x for x in product(range(n + 1), repeat=4) if sum(x) == n]


def add_unit(m, a):
    z = list(m)
    z[a] += 1
    return tuple(z)


def sub_unit(c, a):
    z = list(c)
    z[a] -= 1
    return tuple(z)


def make_parent(n: int):
    cells = compositions(n)
    links = [(m, add_unit(m, a)) for m in cells for a in range(4)]
    idx = {e: i for i, e in enumerate(links)}
    adj = [set() for _ in links]
    for i, e in enumerate(links):
        for j, f in enumerate(links):
            if i != j and (e[0] == f[0] or e[1] == f[1]):
                adj[i].add(j)
    cadj = {m: set() for m in cells}
    for m, n2 in combinations(cells, 2):
        if set(add_unit(m, a) for a in range(4)) & set(
            add_unit(n2, b) for b in range(4)
        ):
            cadj[m].add(n2)
            cadj[n2].add(m)
    return cells, links, idx, adj, cadj


def distances(adj, start):
    d = {start: 0}
    q = deque([start])
    while q:
        x = q.popleft()
        for y in adj[x]:
            if y not in d:
                d[y] = d[x] + 1
                q.append(y)
    return d


def degree_energy_direct(bits, links, ud, dstar):
    node_degree = {}
    for bit, (m, c) in zip(bits, links):
        node_degree[m] = node_degree.get(m, 0) + bit
        node_degree[c] = node_degree.get(c, 0) + bit
    nodes = set(x for e in links for x in e)
    return ud * sum((node_degree.get(v, 0) - dstar) ** 2 for v in nodes)


def degree_energy_split(bits, links, adj, ud, dstar):
    # Each physical link has two endpoints.  Each unordered line-graph edge is
    # one same-node pair because the incidence graph is simple.
    constant = ud * len(set(x for e in links for x in e)) * dstar**2
    local = 2 * ud * (1 - 2 * dstar) * sum(bits)
    pair = 0
    for i in range(len(links)):
        for j in adj[i]:
            if i < j:
                pair += 2 * ud * bits[i] * bits[j]
    return constant + local + pair


checks = 0

for n in range(0, 6):
    cells, links, idx, adj, cadj = make_parent(n)
    for i, (_, c) in enumerate(links):
        q = sum(x > 0 for x in c)
        assert len(adj[i]) == q + 2
        assert len(adj[i]) <= 6
        checks += 2

    # Every link path projects to a cell walk, so the exact finite census must
    # satisfy d_link >= d_cell for every pair.  Exhaust all N<=4; N=5 is used
    # for degree saturation only to keep replay fast.
    if n <= 4:
        cd = {m: distances(cadj, m) for m in cells}
        ld = {i: distances({j: adj[j] for j in range(len(links))}, i)
              for i in range(len(links))}
        for i, e in enumerate(links):
            for j, f in enumerate(links):
                assert ld[i][j] >= cd[e[0]][f[0]]
                assert cd[e[0]][f[0]] == sum(
                    abs(a - b) for a, b in zip(e[0], f[0])
                ) // 2
                checks += 2

    # Exhaust the degree-square decomposition for small blocks and a fixed
    # deterministic sample for the larger blocks.
    samples = []
    if len(links) <= 16:
        samples = product((0, 1), repeat=len(links))
    else:
        samples = [
            tuple((17 * i + 13 * k + n) % 2 for i in range(len(links)))
            for k in range(32)
        ]
    for bits in samples:
        bits = tuple(bits)
        assert degree_energy_direct(bits, links, 3, 2) == degree_energy_split(
            bits, links, adj, 3, 2
        )
        checks += 1

# Exact constants: J=2|Ud|, Delta_L=6, conservative influence row has
# diagonal plus off-diagonal mass 2 J Delta_L, and the commutator recursion
# contributes its leading factor 2/hbar.
ud = 7
j = 2 * ud
delta_l = 6
assert 4 * j * delta_l == 48 * ud
checks += 1

# Exact coefficient identity behind AI18:
# integral u^r/r! du = t^(r+1)/(r+1)!.
for r in range(0, 32):
    assert Fraction(1, factorial(r) * (r + 1)) == Fraction(
        1, factorial(r + 1)
    )
    checks += 1

# Marked-tail inequality on a deterministic numerical grid.  This is a
# diagnostic of the displayed elementary step, not floating-point evidence
# for the operator theorem.
for d in range(0, 8):
    for x10 in range(0, 41):
        x = x10 / 10
        tail = sum(x**r / factorial(r) for r in range(d, 80))
        for mu in (0.25, 0.5, 1.0, 2.0):
            assert tail <= exp(x * exp(mu) - mu * d) * (1 + 1e-13)
            checks += 1

print(f"PASS {checks}/{checks}")
