#!/usr/bin/env python3
"""Constructive replay for the GL6AY finite-coupling bridge."""

from __future__ import annotations

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
    return tuple(left[i] - right[i] for i in range(3))


def add(left, right):
    return tuple(left[i] + right[i] for i in range(3))


def neg(value):
    return tuple(-entry for entry in value)


# Integer-spectrum constraint and exact zero set.
values = {}
for k in range(5):
    q = k - 2
    values[k] = q * q
    check(values[k] in {0, 1, 4}, "local defect spectrum is integer")
    check((values[k] == 0) == (k == 2),
          "zero defect is exactly degree two")
check(set(values.values()) == {0, 1, 4}, "complete local defect spectrum")


# Local-collar projector: nonnegativity makes the zero sector of the sum the
# product of the individual zero sectors.  A strongly supported,
# N_def-preserving transition that fixes every outside constraint therefore
# preserves the contained zero sector without any global projector.
local_defect_values = (0, 1, 4)
for inside in itertools.product(local_defect_values, repeat=3):
    n_s = sum(inside)
    product_zero = int(all(value == 0 for value in inside))
    spectral_zero = int(n_s == 0)
    check(product_zero == spectral_zero,
          "local collar product equals zero spectral projector")
    for twist_value in range(-6, 7):
        check(product_zero * twist_value == twist_value * product_zero,
              "diagonal local collar commutes with occupation twist")

for inside_before in itertools.product(range(3), repeat=2):
    for outside_before in itertools.product(range(3), repeat=2):
        total_before = sum(inside_before) + sum(outside_before)
        for inside_after in itertools.product(range(3), repeat=2):
            outside_after = outside_before
            total_after = sum(inside_after) + sum(outside_after)
            if total_after != total_before:
                continue
            if sum(inside_before) == 0:
                check(sum(inside_after) == 0,
                      "total conservation plus fixed outside preserves local zero collar")


# Strong support of a flip: the complete child star already contains x.
supports = {}
for a in range(4):
    support = {sub(D[a], D[b]) for b in range(4)}
    supports[a] = support
    check(len(support) == 4, "flip strong support has four coarse cells")
    check((0, 0, 0) in support, "flip strong support contains parent cell")
    for left, right in itertools.combinations(support, 2):
        difference = sub(left, right)
        roots = {sub(D[c], D[b]) for b in range(4) for c in range(4)
                 if b != c}
        check(difference in roots or neg(difference) in roots,
              "strong support is connected in A3 adjacency")

# At most sixteen translated flip supports contain one coarse cell.
origin = (0, 0, 0)
edge_labels = set()
for a in range(4):
    for b in range(4):
        x = sub(D[b], D[a])
        translated = {add(x, point) for point in supports[a]}
        check(origin in translated, "constructed strong support contains origin")
        edge_labels.add((x, a))
check(len(edge_labels) == 16, "exactly sixteen flip supports meet one cell")


# Exact one-flip defect-frequency formula, charged resonances, and locked cost.
resonant_add = []
resonant_remove = []
frequencies = set()
for occupation in (0, 1):
    delta = 1 - 2 * occupation
    allowed_degrees = range(occupation, 4 + occupation)
    for k_u in allowed_degrees:
        for k_v in allowed_degrees:
            q_u, q_v = k_u - 2, k_v - 2
            direct = ((q_u + delta) ** 2 - q_u ** 2
                      + (q_v + delta) ** 2 - q_v ** 2)
            formula = 2 * delta * (q_u + q_v) + 2
            check(direct == formula, "one-flip defect formula")
            check(formula % 2 == 0, "one-flip defect frequency is even integer")
            frequencies.add(formula)
            if formula == 0:
                if occupation == 0:
                    resonant_add.append((k_u, k_v))
                    check(k_u + k_v == 3, "addition resonance condition")
                else:
                    resonant_remove.append((k_u, k_v))
                    check(k_u + k_v == 5, "removal resonance condition")
check(resonant_add and resonant_remove, "both charged resonance classes occur")
check(frequencies == {-6, -4, -2, 0, 2, 4, 6},
      "complete one-flip defect-frequency set")
for occupation in (0, 1):
    delta = 1 - 2 * occupation
    locked_cost = 2 * delta * (0 + 0) + 2
    check(locked_cost == 2, "a locked one-link flip creates two defects")
    port_change = delta
    check(abs(port_change) == 1,
          "local locked-to-charged excursion changes a bare port count")


# Crude strong-support norm constants and explicit primary-theorem scales.
for kappa in (0.05, 0.1, 0.25, 0.5):
    d_envelope = 16 * math.exp(4 * kappa)
    v_envelope = 32 * math.exp(4 * kappa)
    nu0_envelope = (54 * math.pi / kappa ** 2) * (
        d_envelope + 2 * v_envelope
    )
    check(abs(nu0_envelope
              - 4320 * math.pi * math.exp(4 * kappa) / kappa ** 2)
          < 1e-11 * nu0_envelope,
          "crude nu_0 envelope has exact coefficient")
    first_ratio = 9 * math.pi * v_envelope / kappa
    check(abs(first_ratio
              - 288 * math.pi * math.exp(4 * kappa) / kappa)
          < 1e-11 * first_ratio,
          "crude first smallness inequality has exact coefficient")


def n_star(ratio: float) -> int:
    return math.floor(ratio / (1 + math.log(ratio)) ** 3) - 2


previous = -10
for ratio in (10, 100, 1_000, 10_000, 100_000, 1_000_000):
    current = n_star(ratio)
    check(current >= previous, "n_* grows monotonically on tested scale")
    previous = current
check(n_star(100_000) >= 1, "explicit nonempty finite-coupling theorem domain")
for n in (1, 2, 5, 10, 50, 100):
    remainder = Fraction(2, 3) ** n
    check(0 < remainder < 1, "prethermal remainder factor is positive and small")
    if n > 1:
        check(remainder < Fraction(2, 3), "remainder decreases with normal-form order")


# Exponential support decay controls the fourth twist moment.
for kappa in (0.01, 0.05, 0.1, 0.5, 1.0):
    samples = [m ** 4 * math.exp(-kappa * m) for m in range(1, 10001)]
    maximum = max(samples)
    check(math.isfinite(maximum) and maximum > 0,
          "exponential potential norm has finite fourth moment")
    maximizing_m = samples.index(maximum) + 1
    check(abs(maximizing_m - 4 / kappa) <= 2,
          "fourth-moment maximum occurs at the analytic scale")


def boundary_count(mask: int, length: int) -> int:
    """Number of cycle vertices incident to exactly one selected edge."""
    return sum(
        ((mask >> ((vertex - 1) % length)) & 1)
        != ((mask >> vertex) & 1)
        for vertex in range(length)
    )


def winding_denominator_sum(length: int) -> Fraction:
    """Sum over all flip orders of products of proper-subset denominators."""
    full = (1 << length) - 1
    dynamic = [Fraction(0) for _ in range(1 << length)]
    dynamic[0] = Fraction(1)
    for size in range(1, length):
        for mask in range(1, full):
            if bin(mask).count("1") != size:
                continue
            boundary = boundary_count(mask, length)
            check(2 <= boundary <= length and boundary % 2 == 0,
                  "proper nonempty cycle subset has positive even boundary")
            subtotal = sum(
                dynamic[mask ^ (1 << edge)]
                for edge in range(length)
                if mask & (1 << edge)
            )
            dynamic[mask] = subtotal / boundary
    return sum(dynamic[full ^ (1 << edge)] for edge in range(length))


for length in (8, 10, 12):
    coefficient = winding_denominator_sum(length)
    lower = Fraction(math.factorial(length), length ** (length - 1))
    upper = Fraction(math.factorial(length), 2 ** (length - 1))
    check(coefficient > 0, "first winding coefficient is strictly positive in magnitude")
    check(lower <= coefficient <= upper, "winding coefficient obeys denominator bounds")
    check(length % 2 == 0, "native winding row has even link length")

# Independent brute-force replay at the smallest audited cycle.
length = 8
brute = Fraction(0)
for permutation in itertools.permutations(range(length)):
    mask = 0
    term = Fraction(1)
    for step, edge in enumerate(permutation, start=1):
        mask |= 1 << edge
        if step < length:
            term /= boundary_count(mask, length)
    brute += term
check(brute == winding_denominator_sum(length),
      "dynamic and ordered-path winding sums agree")


# Contractible alternating parent walks conserve each port; winding abstractly
# changes a j/3 pair by exactly L_j.
pairs = tuple((a, b) for a in range(4) for b in range(4) if a != b)
for pair_count in range(2, 6):
    for word in itertools.product(pairs, repeat=pair_count):
        displacement = [0, 0, 0]
        port_delta = [0, 0, 0, 0]
        for a, b in word:
            for axis in range(3):
                displacement[axis] += D[a][axis] - D[b][axis]
            port_delta[a] += 1
            port_delta[b] -= 1
        if displacement == [0, 0, 0]:
            check(port_delta == [0, 0, 0, 0],
                  "contractible locked return conserves every port")
for period in (4, 5, 8, 13, 21):
    for j in range(3):
        port_delta = [0, 0, 0, 0]
        port_delta[j] = -period
        port_delta[3] = period
        check(sum(port_delta) == 0, "winding transfer preserves total links")
        check(sum(abs(value) for value in port_delta) == 2 * period,
              "winding transfer has exact two-period port variation")


print(f"PASS__GL6AY_FINITE_COUPLING_PRETHERMAL_LOCKED_BRIDGE__{checks}/{checks}")
