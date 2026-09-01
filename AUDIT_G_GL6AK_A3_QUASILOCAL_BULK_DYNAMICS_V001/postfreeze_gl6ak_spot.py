#!/usr/bin/env python3
"""Distinct post-freeze spot reconstruction for GL6AK.

This does not import either author verifier or the pre-freeze replay.
"""

import cmath
import itertools
import math
from collections import defaultdict
from fractions import Fraction


checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def sub(x, y):
    return tuple(a - b for a, b in zip(x, y))


ports = range(4)
unit = tuple(tuple(int(i == a) for i in ports) for a in ports)
pairs = tuple(itertools.combinations(ports, 2))
zero = (0, 0, 0, 0)


def radius(x):
    check(sum(x) == 0, "A3 point")
    return sum(abs(v) for v in x) // 2


def ball(r):
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
    return frozenset(((x, a), (add(x, sub(unit[a], unit[b])), b)))


def child(site):
    x, a = site
    return add(x, unit[a])


def edges(cells):
    result = set()
    for x in cells:
        for a, b in pairs:
            result.add(internal(x, a, b))
            result.add(shared(x, a, b))
    return result


# Incidence and the 3+3 bulk degree.
edge_set = edges(ball(3))
for x in ball(2):
    touching = {edge for edge in edge_set if any(p[0] == x for p in edge)}
    check(len(touching) == 18, "cell touches 18 terms")
    for a in ports:
        incident = {edge for edge in edge_set if (x, a) in edge}
        check(len(incident) == 6, "site degree six")
        same_parent = sum(len({p[0] for p in edge}) == 1 for edge in incident)
        same_child = sum(child(tuple(edge)[0]) == child(tuple(edge)[1])
                         for edge in incident)
        check((same_parent, same_child) == (3, 3), "degree split")

# A non-symmetric collared finite patch embeds in one strict-interior slab.
patch = {(0, 0, 0, 0), (3, -2, -1, 0), (-2, 1, 2, -1)}
collar = patch | {
    add(x, sub(unit[a], unit[b]))
    for x in patch for a in ports for b in ports if a != b
}
m = tuple(max(1, 1 - min(x[i] for x in collar)) for i in ports)
N = sum(m)
for x in collar:
    image = add(m, x)
    check(min(image) >= 1 and sum(image) == N, "strict-interior embedding")
for x in patch:
    for a, b in pairs:
        y = add(x, sub(unit[a], unit[b]))
        check(add(add(m, x), unit[a]) == add(add(m, y), unit[b]),
              "literal shared child")

# Minimum-radius ownership and boundary constants.
outer_edges = edges(ball(6))
owned = defaultdict(int)
for edge in outer_edges:
    endpoint_cells = {p[0] for p in edge}
    if max(radius(x) for x in endpoint_cells) > 5:
        continue
    owner = min(endpoint_cells, key=lambda x: (radius(x), x))
    owned[owner] += 1
check(max(owned.values()) <= 18, "owner ceiling 18")
for r in range(5):
    shell = ball(r) - (ball(r - 1) if r else set())
    nr = sum(count for x, count in owned.items() if radius(x) == r)
    check(nr <= 18 * len(shell), "shell term ceiling")
    check(len(shell) <= (2 * r + 1) ** 3, "cubic shell ceiling")
for J in (Fraction(1, 9), Fraction(7, 4), Fraction(13, 2)):
    hbar = Fraction(11, 3)
    lam = 24 * J / hbar
    check(72 * J / hbar / lam == 3, "boundary factor three")


def tail(d, x):
    log_term = d * math.log(x) - math.lgamma(d + 1)
    term = 0.0 if log_term < -745 else math.exp(log_term)
    if term == 0.0:
        return 0.0
    total = term
    for k in range(d + 1, 250):
        term *= x / k
        total += term
        if term < max(1e-300, total * 1e-17):
            break
    return total


for x in (0.5, 2.0, 5.0):
    values = [sum((2 * r + 1) ** 3 * tail(r - 1, x)
                  for r in range(R, 180)) for R in (8, 16, 32, 64)]
    check(all(values[i + 1] < values[i] for i in range(3)), "tail decrease")
    check(values[-1] < 1e-30, "tail vanishes")

# Concrete Følner sequence behind translation averaging.
def folner(r):
    return {(a, b, c, -a - b - c)
            for a in range(-r, r + 1)
            for b in range(-r, r + 1)
            for c in range(-r, r + 1)}


for z in ((1, -1, 0, 0), (2, -1, -1, 0), (-3, 1, 1, 1)):
    ratios = []
    for r in (4, 12, 36):
        f = folner(r)
        shifted = {add(x, z) for x in f}
        ratios.append(len(f.symmetric_difference(shifted)) / len(f))
    check(ratios[2] < ratios[1] < ratios[0], "Folner decay")

# Six-pair projector algebra.
def mm(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(6))
             for j in range(6)] for i in range(6)]


I = [[Fraction(int(i == j)) for j in range(6)] for i in range(6)]
A = [[Fraction(int(i != j and bool(set(pairs[i]) & set(pairs[j]))))
      for j in range(6)] for i in range(6)]
P1 = [[Fraction(1, 6) for _ in range(6)] for _ in range(6)]
PE = [[v / 12 for v in row]
      for row in mm(A, [[A[i][j] - 4 * I[i][j] for j in range(6)]
                        for i in range(6)])]
PT = [[-v / 8 for v in row]
      for row in mm([[A[i][j] - 4 * I[i][j] for j in range(6)]
                     for i in range(6)],
                    [[A[i][j] + 2 * I[i][j] for j in range(6)]
                     for i in range(6)])]
Z = [[Fraction(0) for _ in range(6)] for _ in range(6)]
for P, rank in ((P1, 1), (PE, 2), (PT, 3)):
    check(mm(P, P) == P, "projector")
    check(sum(P[i][i] for i in range(6)) == rank, "projector rank")
for P, Q in ((P1, PE), (P1, PT), (PE, PT)):
    check(mm(P, Q) == Z, "orthogonal sectors")
check([[P1[i][j] + PE[i][j] + PT[i][j] for j in range(6)]
       for i in range(6)] == I, "projector resolution")

# Retarded sign for U(t)=exp(itL), using one positive-frequency atom.
nu, t, hbar = 1.7, 0.31, 0.8
weight = 0.43
Fplus = weight * cmath.exp(1j * nu * t)
Fminus = weight * cmath.exp(-1j * nu * t)
chi = -1j / hbar * (Fminus - Fplus)
check(abs(chi.imag) < 1e-14, "retarded reality")
check(abs(chi.real + 2 * weight * math.sin(nu * t) / hbar) < 1e-14,
      "retarded sign")

print(f"PASS__INDEPENDENT_GL6AK_POSTFREEZE_SPOT__{checks}/{checks}")
print("INCIDENCE=DEGREE6_3PLUS3;BOUNDARY=18_72J_FACTOR3_TAIL")
print("STATE=TIME_FOLNER_S4_EXISTENCE_ONLY;SPECTRAL=A1_E_T2_POSITIVE")
print("CEILING=NO_GLOBAL_RECORD_STATE_SELECTION_POLE_MOMENTUM_RICCI_GRAVITY_G")
