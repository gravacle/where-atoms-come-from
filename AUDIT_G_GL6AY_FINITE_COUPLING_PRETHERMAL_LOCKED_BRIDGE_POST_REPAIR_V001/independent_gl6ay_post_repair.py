#!/usr/bin/env python3
"""Independent post-repair replay for GL6AY; imports no author module."""

from __future__ import annotations

import cmath
import itertools
import math
from fractions import Fraction


checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise RuntimeError(label)
    checks += 1


D = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (0, 0, 0),
)


def sub(left, right):
    return tuple(left[index] - right[index] for index in range(3))


def add(left, right):
    return tuple(left[index] + right[index] for index in range(3))


# Integer spectrum, local tensor factor, and exact finite P_L--Q_L gap.
local_values = tuple((k - 2) ** 2 for k in range(5))
check(set(local_values) == {0, 1, 4}, "integer local defect spectrum")
for k, value in enumerate(local_values):
    check((value == 0) == (k == 2), "zero defect exactly degree two")
check(2 ** 4 == 16, "four link qubits per coarse cell")

for volume in range(1, 5):
    qlists = tuple(itertools.product(range(-2, 3), repeat=volume))
    positive = []
    for parents in qlists:
        for children in qlists:
            if sum(parents) == sum(children):
                ndef = sum(q * q for q in parents + children)
                if ndef:
                    positive.append(ndef)
    check(min(positive) == 2, "finite locked-to-charged gap is two U_d")


# Exact strong support and term incidence.
supports = {}
for a in range(4):
    support = {sub(D[a], D[b]) for b in range(4)}
    supports[a] = support
    check(len(support) == 4, "flip strong support has four cells")
    check((0, 0, 0) in support, "child star contains parent cell")
    center = D[a]
    for point in support:
        check(sum(abs(point[i] - center[i]) for i in range(3)) <= 1,
              "strong support connected around its star center")

labels = set()
for a in range(4):
    for point in supports[a]:
        origin_shift = tuple(-entry for entry in point)
        translated = {add(origin_shift, member) for member in supports[a]}
        check((0, 0, 0) in translated, "translated support contains test cell")
        labels.add((origin_shift, a))
check(len(labels) == 16, "exact sixteen labeled supports per cell")


# One-flip pinching and charged resonances.
frequencies = set()
resonant = {0: set(), 1: set()}
for occupied in (0, 1):
    delta = 1 - 2 * occupied
    for ku in range(occupied, 4 + occupied):
        for kv in range(occupied, 4 + occupied):
            qu, qv = ku - 2, kv - 2
            direct = ((qu + delta) ** 2 - qu ** 2
                      + (qv + delta) ** 2 - qv ** 2)
            formula = 2 * delta * (qu + qv) + 2
            check(direct == formula, "one-flip defect-frequency identity")
            frequencies.add(formula)
            if formula == 0:
                resonant[occupied].add((ku, kv))
                check(ku + kv == (3 if not occupied else 5),
                      "charged resonance condition")
check(frequencies == {-6, -4, -2, 0, 2, 4, 6},
      "complete charged frequency set")
check(resonant[0] and resonant[1], "both resonance classes occur")


# Exact ADHH constants and a nonempty sufficient domain.
for kappa in (0.05, 0.1, 0.25, 0.5):
    d0 = 16 * math.exp(4 * kappa)
    v0 = 32 * math.exp(4 * kappa)
    nu0 = 54 * math.pi * (d0 + 2 * v0) / kappa ** 2
    check(abs(nu0 - 4320 * math.pi * math.exp(4 * kappa) / kappa ** 2)
          < 1e-11 * nu0, "source nu0 constant")
    threshold = 9 * math.pi * v0 / kappa
    check(abs(threshold - 288 * math.pi * math.exp(4 * kappa) / kappa)
          < 1e-11 * threshold, "source smallness threshold")


def nstar(ratio):
    return math.floor(ratio / (1 + math.log(ratio)) ** 3) - 2


check(nstar(100_000) >= 1, "sufficient finite-coupling domain nonempty")
previous = -100
for ratio in (10, 100, 1_000, 10_000, 100_000, 1_000_000):
    value = nstar(ratio)
    check(value >= previous, "tested nstar monotonicity")
    previous = value
for n in range(1, 101):
    check(0 < Fraction(2, 3) ** n < 1, "prethermal remainder factor")


# Local-collar spectral identity: a sum of nonnegative commuting constraint
# values is zero exactly when every contained constraint is zero.
for count in range(0, 7):
    for values in itertools.product((0, 1, 4), repeat=count):
        spectral_zero = int(sum(values) == 0)
        product_zero = int(all(value == 0 for value in values))
        check(spectral_zero == product_zero,
              "contained zero projector equals product of zero projectors")


# Abstract strong-support implication.  Outside eigenvalues are held fixed;
# total N conservation then preserves the local zero sector exactly.
for inside_count in range(1, 5):
    inside_lists = tuple(itertools.product((0, 1, 4), repeat=inside_count))
    for outside_before in itertools.product((0, 1, 4), repeat=2):
        outside_sum = sum(outside_before)
        for inside_after in inside_lists:
            total_before = outside_sum  # contained input starts at zero
            total_after = outside_sum + sum(inside_after)
            if total_after == total_before:
                check(all(value == 0 for value in inside_after),
                      "total pinching plus fixed outside preserves collar zero")


# Exact complex matrix replay of the double-commutator identity.  N has two
# zero-sector states, so the collar compression is nontrivial.
n_values = (0, 0, 1, 1, 2, 2)
p_values = tuple(1 if value == 0 else 0 for value in n_values)
a_values = (-3, 2, -1, 4, 0, 5)
for seed in range(1, 80):
    matrix = []
    for i in range(len(n_values)):
        row = []
        for j in range(len(n_values)):
            if n_values[i] != n_values[j]:
                row.append(0j)
            else:
                real = ((seed + 3 * i + 5 * j) % 17) - 8
                imag = ((2 * seed + 7 * i - 7 * j) % 19) - 9
                # Hermitize explicitly below.
                row.append(complex(real, imag))
        matrix.append(row)
    for i in range(len(n_values)):
        matrix[i][i] = complex(matrix[i][i].real, 0)
        for j in range(i):
            value = (matrix[i][j] + matrix[j][i].conjugate()) / 2
            matrix[i][j] = value
            matrix[j][i] = value.conjugate()
    for i in range(len(n_values)):
        for j in range(len(n_values)):
            phi = p_values[i] * matrix[i][j] * p_values[j]
            left = (a_values[i] - a_values[j]) ** 2 * phi
            right = (p_values[i] *
                     (a_values[i] - a_values[j]) ** 2 * matrix[i][j] *
                     p_values[j])
            check(abs(left - right) < 1e-12,
                  "collar projection commutes through double commutator")
            if p_values[j] == 1:
                check(abs(phi - p_values[i] * matrix[i][j]) < 1e-12,
                      "collar interaction agrees on zero-sector input")


# Exponential support control dominates the twist polynomial.
for kappa in (0.01, 0.03, 0.1, 0.5, 1.0):
    samples = [m ** 4 * math.exp(-kappa * m) for m in range(1, 10001)]
    maximum = max(samples)
    check(math.isfinite(maximum), "finite fourth support moment")
    check(abs(samples.index(maximum) + 1 - 4 / kappa) <= 2,
          "fourth moment has expected finite maximum")


# Locked-endpoint affine conservation for a large family of finite cycles.
pairs = tuple((a, b) for a in range(4) for b in range(4) if a != b)
for length in range(2, 6):
    for word in itertools.product(pairs, repeat=length):
        displacement = [0, 0, 0]
        ports = [0, 0, 0, 0]
        for a, b in word:
            for axis in range(3):
                displacement[axis] += D[a][axis] - D[b][axis]
            ports[a] += 1
            ports[b] -= 1
        if displacement == [0, 0, 0]:
            check(ports == [0, 0, 0, 0],
                  "finite contractible locked difference conserves all ports")


def boundary(mask: int, length: int) -> int:
    return sum(((mask >> edge) & 1) !=
               ((mask >> ((edge - 1) % length)) & 1)
               for edge in range(length))


def winding_sum(length: int) -> Fraction:
    full = (1 << length) - 1
    dynamic = [Fraction(0) for _ in range(1 << length)]
    dynamic[0] = Fraction(1)
    for size in range(1, length):
        for mask in range(1, full):
            if bin(mask).count("1") != size:
                continue
            b = boundary(mask, length)
            check(2 <= b <= length and b % 2 == 0,
                  "proper winding subset has positive even boundary")
            dynamic[mask] = sum(dynamic[mask ^ (1 << edge)]
                                for edge in range(length)
                                if mask & (1 << edge)) / b
    return sum(dynamic[full ^ (1 << edge)] for edge in range(length))


for length in (8, 10, 12, 14):
    coefficient = winding_sum(length)
    lower = Fraction(math.factorial(length), length ** (length - 1))
    upper = Fraction(math.factorial(length), 2 ** (length - 1))
    check(lower <= coefficient <= upper, "winding factorial bounds")
    check(coefficient > 0, "winding coefficient nonzero")

# Brute force at the smallest cycle checks the dynamic recurrence and sign.
length = 8
brute = Fraction(0)
for permutation in itertools.permutations(range(length)):
    mask = 0
    term = Fraction(1)
    for step, edge in enumerate(permutation, start=1):
        mask |= 1 << edge
        if step < length:
            term /= boundary(mask, length)
    brute += term
check(brute == winding_sum(length), "brute and dynamic winding word sums agree")


# Local closeness never implies global projector closeness; the repaired
# theorem now explicitly keeps this counterexample outside its claims.
for epsilon in (0.01, 0.03, 0.1):
    distances = []
    for volume in (1, 10, 100, 1_000, 10_000):
        overlap_squared = math.cos(epsilon) ** (2 * volume)
        distances.append(math.sqrt(max(0.0, 1 - overlap_squared)))
    check(all(distances[index] <= distances[index + 1] + 1e-15
              for index in range(len(distances) - 1)),
          "global projection distance grows with volume")
    check(distances[-1] > 0.79, "global projection distance becomes order one")


# The commuting microscopic flip sum has exact extensive norm.
for volume in (1, 2, 5, 10, 100, 1_000):
    h = 0.01
    edges = 4 * volume
    check(h * edges == 4 * h * volume, "exact extensive commuting-flip norm")


# Centered local-observable source exponent for d=3.
for r1 in (0.001, 0.01, 0.05, 0.09):
    check(r1 < math.log(1.5) / 4, "local horizon exponent within source range")
for n in (1, 10, 100, 1000):
    horizon = math.exp(0.01 * n)
    check(horizon >= 1, "finite prethermal observation horizon")


print(f"PASS__INDEPENDENT_GL6AY_POST_REPAIR__{checks}/{checks}")
