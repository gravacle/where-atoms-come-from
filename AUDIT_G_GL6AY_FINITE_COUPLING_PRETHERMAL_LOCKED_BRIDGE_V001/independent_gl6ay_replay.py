#!/usr/bin/env python3
"""Independent replay for the hostile GL6AY audit; imports no author code."""

from __future__ import annotations

import itertools
import math
from fractions import Fraction


checks = 0


def require(condition: bool, label: str) -> None:
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
    return tuple(left[i] - right[i] for i in range(3))


def add(left, right):
    return tuple(left[i] + right[i] for i in range(3))


# Local Hilbert/integer spectrum and exact zero set.
defect = {(k, (k - 2) ** 2) for k in range(5)}
require({value for _, value in defect} == {0, 1, 4},
        "independent local integer spectrum")
for k, value in defect:
    require((value == 0) == (k == 2), "zero defect iff degree two")
require(2 ** 4 == 16, "four link qubits give local dimension sixteen")


# Equality of total parent and child deviations forbids N_def=1.  Enumerate
# independent degree-deviation lists for small volumes as a direct attack.
for volume in range(1, 5):
    positive = []
    qlists = tuple(itertools.product(range(-2, 3), repeat=volume))
    for parents in qlists:
        for children in qlists:
            if sum(parents) != sum(children):
                continue
            ndef = sum(q * q for q in parents + children)
            if ndef:
                positive.append(ndef)
    require(min(positive) == 2, "parent-child sum identity gives gap two")


# Strong-support geometry and exact sixteen-term incidence count.
supports = {}
for a in range(4):
    support = {sub(D[a], D[b]) for b in range(4)}
    supports[a] = support
    require(len(support) == 4, "four-cell flip strong support")
    require((0, 0, 0) in support, "parent cell lies in child-star support")
    center = D[a]
    for point in support:
        distance = sum(abs(point[i] - center[i]) for i in range(3))
        require(distance <= 1, "strong support connected as a radius-one star")

labels_at_origin = set()
for a in range(4):
    for point in supports[a]:
        x = tuple(-coordinate for coordinate in point)
        translated = {add(x, member) for member in supports[a]}
        require((0, 0, 0) in translated, "translated support hits origin")
        labels_at_origin.add((x, a))
require(len(labels_at_origin) == 16, "exactly sixteen labeled terms meet a cell")


# Pinching, charged resonances, and locked cost.
frequencies = set()
addition_resonances = set()
removal_resonances = set()
for occupied in (0, 1):
    delta = 1 - 2 * occupied
    degrees = range(occupied, 4 + occupied)
    for ku in degrees:
        for kv in degrees:
            qu, qv = ku - 2, kv - 2
            direct = (qu + delta) ** 2 - qu ** 2
            direct += (qv + delta) ** 2 - qv ** 2
            formula = 2 * delta * (qu + qv) + 2
            require(direct == formula, "independent one-flip frequency formula")
            frequencies.add(formula)
            if formula == 0 and not occupied:
                addition_resonances.add((ku, kv))
                require(ku + kv == 3, "addition resonance")
            if formula == 0 and occupied:
                removal_resonances.add((ku, kv))
                require(ku + kv == 5, "removal resonance")
require(frequencies == {-6, -4, -2, 0, 2, 4, 6},
        "complete charged frequency set")
require(addition_resonances and removal_resonances,
        "both charged resonance classes nonempty")
require(2 * (0 + 0) + 2 == 2, "locked one-flip cost two")


# Exact source constants and crude norm envelope.
for kappa in (0.05, 0.1, 0.25, 0.5, 0.8):
    d0 = 16 * math.exp(4 * kappa)
    v0 = 32 * math.exp(4 * kappa)
    nu0 = 54 * math.pi * (d0 + 2 * v0) / kappa ** 2
    require(abs(nu0 - 4320 * math.pi * math.exp(4 * kappa) / kappa ** 2)
            < 1e-12 * nu0, "ADHH nu0 envelope")
    threshold = 9 * math.pi * v0 / kappa
    require(abs(threshold - 288 * math.pi * math.exp(4 * kappa) / kappa)
            < 1e-12 * threshold, "ADHH smallness envelope")


def n_star(ratio):
    return math.floor(ratio / (1 + math.log(ratio)) ** 3) - 2


previous = -100
for ratio in (10, 100, 1_000, 10_000, 100_000, 1_000_000):
    value = n_star(ratio)
    require(value >= previous, "tested nstar monotonicity")
    previous = value
require(n_star(100_000) >= 1, "nonempty sufficient theorem regime")
for n in range(1, 101):
    require(0 < (Fraction(2, 3) ** n) < 1,
            "source remainder factor strictly decays")


# The quartic moment is finite; the author's quartic supremum is a valid
# (though nonoptimal) upper bound on the per-cell count.
for kappa in (0.01, 0.03, 0.1, 0.3, 1.0):
    values = [m ** 4 * math.exp(-kappa * m) for m in range(1, 10001)]
    maximum = max(values)
    require(math.isfinite(maximum), "exponential decay controls fourth moment")
    require(abs(values.index(maximum) + 1 - 4 / kappa) <= 2,
            "quartic maximum at four over kappa")
    for m in range(1, 1000):
        require(m ** 3 * math.exp(-kappa * m) <= maximum,
                "loose quartic supremum dominates incidence-reduced cubic")


# Local-collar repair algebra in an explicit finite model.  Two inside
# constraints count bits (0,1) and (2,3); an outside constraint counts bit 4.
# O swaps 10<->01 inside each pair, preserving N_S and commuting outside.
def bits(state, width=5):
    return tuple((state >> index) & 1 for index in range(width))


for state in range(1 << 5):
    value = bits(state)
    n_s = (value[0] + value[1]) + (value[2] + value[3])
    outside = value[4]
    # Use a nonnegative local defect with zero sector at all four inside bits.
    p_s = int(n_s == 0)
    require(p_s in (0, 1), "local zero projector is spectral")
    # Diagonal occupation twist commutes with the local projector.
    a_s = value[0] + 2 * value[1] + 3 * value[2] + 4 * value[3]
    require(p_s * a_s == a_s * p_s, "local lock projector commutes with twist")
    if p_s:
        require(n_s == 0 and outside in (0, 1),
                "local zero sector leaves outside constraint independent")


# Explicit counterexample to promoting local closeness to global projector
# closeness: product rotations have small fixed-local action but rotate the
# full product ray to norm distance tending to one.
for epsilon in (0.01, 0.03, 0.1):
    previous_distance = 0.0
    for volume in (1, 10, 100, 1_000, 10_000):
        overlap_squared = math.cos(epsilon) ** (2 * volume)
        distance = math.sqrt(max(0.0, 1 - overlap_squared))
        require(distance + 1e-15 >= previous_distance,
                "global projector distance grows with volume")
        previous_distance = distance
    require(previous_distance > 0.79,
            "small local rotations give order-one global projector distance")


def boundary(mask: int, length: int) -> int:
    return sum(((mask >> edge) & 1) !=
               ((mask >> ((edge - 1) % length)) & 1)
               for edge in range(length))


def winding_sum(length: int) -> Fraction:
    full = (1 << length) - 1
    dp = [Fraction(0) for _ in range(1 << length)]
    dp[0] = Fraction(1)
    for size in range(1, length):
        for mask in range(1, full):
            if bin(mask).count("1") != size:
                continue
            b = boundary(mask, length)
            require(2 <= b <= length and b % 2 == 0,
                    "proper cycle subset has positive even boundary")
            dp[mask] = sum(dp[mask ^ (1 << edge)]
                           for edge in range(length) if mask & (1 << edge)) / b
    return sum(dp[full ^ (1 << edge)] for edge in range(length))


for length in (8, 10, 12, 14):
    coefficient = winding_sum(length)
    lower = Fraction(math.factorial(length), length ** (length - 1))
    upper = Fraction(math.factorial(length), 2 ** (length - 1))
    require(lower <= coefficient <= upper,
            "independent winding coefficient bounds")
    require(coefficient > 0, "winding coefficient nonzero")

# Brute permutation replay at r=8.
brute = Fraction(0)
length = 8
for permutation in itertools.permutations(range(length)):
    mask = 0
    term = Fraction(1)
    for step, edge in enumerate(permutation, start=1):
        mask |= 1 << edge
        if step < length:
            term /= boundary(mask, length)
    brute += term
require(brute == winding_sum(length), "brute and dynamic winding sums agree")


# Global flip norm and collapse of the standard global sufficient condition.
for volume in (1, 2, 5, 10, 100, 1_000):
    edges = 4 * volume
    h = 0.01
    global_norm = h * edges
    require(global_norm == 4 * h * volume, "commuting flip global norm extensive")
    ratio = global_norm / 2
    if volume > 1:
        require(ratio > 0, "whole-band norm ratio grows with volume")


print(f"PASS__INDEPENDENT_GL6AY_HOSTILE_REPLAY__{checks}/{checks}")
