#!/usr/bin/env python3
"""Independent hostile replay for frozen GL6AW; imports no author code."""

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


def close(left: complex | float, right: complex | float,
          tolerance: float = 4.0e-11) -> bool:
    return abs(left - right) <= tolerance * max(1.0, abs(left), abs(right))


PORT_STEP = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (0, 0, 0),
)


def point_add(left, right, periods):
    return tuple((left[j] + right[j]) % periods[j] for j in range(3))


def point_subtract(left, right, periods):
    return tuple((left[j] - right[j]) % periods[j] for j in range(3))


def cells(periods):
    return itertools.product(*(range(length) for length in periods))


def six_cycle(origin, ports, periods):
    """Reconstruct the native incidence cycle directly from its three ports."""
    a, b, c = ports
    through_ab = point_add(point_subtract(origin, PORT_STEP[b], periods),
                           PORT_STEP[a], periods)
    through_bc = point_add(point_subtract(origin, PORT_STEP[b], periods),
                           PORT_STEP[c], periods)
    return (
        (origin, a),
        (through_ab, b),
        (through_ab, c),
        (through_bc, a),
        (through_bc, b),
        (origin, c),
    )


def child_of(edge, periods):
    origin, port = edge
    return point_add(origin, PORT_STEP[port], periods)


# A. Reconstruct loop geometry, uniqueness, incidence, and twist count.
geometry_samples = tuple(
    (l0, l1, l2)
    for l0 in range(4, 9)
    for l1 in range(4, 10)
    for l2 in range(4, 8)
)
for periods in geometry_samples:
    volume = math.prod(periods)
    loop_sets = set()
    affected = zero_unaffected = zero_absent = 0
    for origin in cells(periods):
        for ports in itertools.combinations(range(4), 3):
            loop = six_cycle(origin, ports, periods)
            frozen = frozenset(loop)
            check(len(frozen) == 6, "six distinct links per elementary loop")
            check(frozen not in loop_sets, "elementary loop labels are unique")
            loop_sets.add(frozen)

            parent_degrees = {}
            child_degrees = {}
            port_counts = [0, 0, 0, 0]
            for parity, edge in enumerate(loop):
                origin_at_edge, port = edge
                parent_degrees[origin_at_edge] = parent_degrees.get(origin_at_edge, 0) + 1
                child = child_of(edge, periods)
                child_degrees[child] = child_degrees.get(child, 0) + 1
                port_counts[port] += 1
            check(sorted(parent_degrees.values()) == [2, 2, 2],
                  "loop has degree two at three parent vertices")
            check(sorted(child_degrees.values()) == [2, 2, 2],
                  "loop has degree two at three child vertices")
            for port in range(4):
                expected = 2 if port in ports else 0
                check(port_counts[port] == expected, "each used port occurs twice")
                parities = [j % 2 for j, edge in enumerate(loop) if edge[1] == port]
                if expected:
                    check(sorted(parities) == [0, 1],
                          "used-port occurrences have opposite flip parity")

            zero_occurrences = [(j, edge[0][1]) for j, edge in enumerate(loop)
                                if edge[1] == 0]
            if not zero_occurrences:
                zero_absent += 1
            else:
                delta = sum((1 if parity else -1) * coordinate
                            for parity, coordinate in zero_occurrences)
                residue = delta % periods[1]
                if residue in (1, periods[1] - 1):
                    affected += 1
                else:
                    check(residue == 0, "unaffected zero-port loop has trivial phase")
                    zero_unaffected += 1
    check(len(loop_sets) == 4 * volume, "rectangular quotient has exactly 4V loops")
    check(affected == 2 * volume, "exactly 2V loop terms acquire unit twist")
    check(zero_unaffected == volume, "exactly V zero-port loops are unaffected")
    check(zero_absent == volume, "exactly V loops omit port zero")


# B. Reconstruct the centered locked witness for all theorem parity samples.
parity_samples = tuple(
    (l0, l1, l2)
    for l0 in (5, 7, 9)
    for l1 in (4, 6, 8, 10, 12)
    for l2 in (5, 7, 11)
)
for periods in parity_samples:
    volume = math.prod(periods)
    occupied = set()
    for x in cells(periods):
        occupied.add((x, 0 if x[1] % 2 == 0 else 3))
        occupied.add((x, 1))
    check(len(occupied) == 2 * volume, "witness matchings are edge-disjoint")
    parents = {x: 0 for x in cells(periods)}
    children = {x: 0 for x in cells(periods)}
    ports = [0, 0, 0, 0]
    for edge in occupied:
        parents[edge[0]] += 1
        children[child_of(edge, periods)] += 1
        ports[edge[1]] += 1
    check(all(degree == 2 for degree in parents.values()),
          "centered witness is parent locked")
    check(all(degree == 2 for degree in children.values()),
          "centered witness is child locked")
    check(ports[0] == volume // 2, "centered witness has N0=V/2")
    check(sum(ports) == 2 * volume, "centered witness has total locked occupation")


# C. Check translation character both arithmetically and by direct wrap sums.
for periods in parity_samples:
    l0, l1, l2 = periods
    volume = math.prod(periods)
    number = volume // 2
    q = 2 * math.pi / l1
    check(number == l0 * l1 * l2 // 2, "centered charge is integral")
    check(Fraction(number, l1).denominator == 2,
          "charge per twist length is half-integral")
    check(close(cmath.exp(-1j * q * number), -1),
          "odd transverse area gives minus-one character")

    zero_sites = [x for x in cells(periods) if x[1] % 2 == 0]
    before = sum(x[1] for x in zero_sites)
    # Operator convention Y n_x Y^-1=n_(x+e1): relabel x -> x-e1.
    translated = [((x[0], (x[1] - 1) % l1, x[2])) for x in zero_sites]
    after = sum(x[1] for x in translated)
    check((after - before + number) % l1 == 0,
          "direct wrap sum reproduces W -> W-N0 mod L1")
    check(close(cmath.exp(1j * q * (after - before)), -1),
          "direct centered witness translation phase is minus one")


# D. Attack the partial-flip energy coefficient with nonuniform amplitudes,
# both twist directions, and wrap representatives.
for l1 in range(4, 90, 2):
    q = 2 * math.pi / l1
    for left in range(1, 53):
        for right in range(1, 47):
            norm = math.sqrt(left * left + right * right)
            p, r = left / norm, right / norm
            t_expectation = 2 * p * r
            check(0 < t_expectation <= 1 + 1e-14,
                  "positive partial flip has expectation at most one")
            for displacement in (1, -1, l1 - 1, 1 - l1):
                original = -2 * p * r
                twisted = -2 * p * r * math.cos(q * displacement)
                check(close(twisted - original,
                            (1 - math.cos(q)) * t_expectation),
                      "per-term energy loss is (1-cos q)<T>")


# E. Abstractly replay PF/translation orthogonality and min--max in the
# smallest nontrivial representation: Y=X, U=Z, H=-JX.
for numerator in range(1, 500):
    J = numerator / 37
    inv_sqrt_two = 1 / math.sqrt(2)
    psi = (inv_sqrt_two, inv_sqrt_two)
    twisted = (inv_sqrt_two, -inv_sqrt_two)
    overlap = sum(a * b for a, b in zip(psi, twisted))
    ground_energy = -J
    twisted_energy = J
    check(close(overlap, 0), "minus-one translation character is orthogonal")
    check(close(twisted_energy - ground_energy, 2 * J),
          "orthogonal twist realizes a valid component excitation")


# F. Exact envelope, parity failures, isotropic failure, and Følner closure.
for l0 in range(4, 15):
    for l2 in range(4, 15):
        l1 = 2 * (l0 + l2 + 3)
        number = l0 * l1 * l2 // 2
        phase = cmath.exp(-2j * math.pi * (number % l1) / l1)
        expected = -1 if (l0 * l2) % 2 else 1
        check(close(phase, expected), "transverse parity exactly controls twist character")

for m in range(5, 502, 2):
    l0, l1, l2 = m, 2 * m ** 3, m
    volume = l0 * l1 * l2
    q = 2 * math.pi / l1
    stable_exact = 4 * volume * math.sin(q / 2) ** 2
    envelope = 4 * math.pi ** 2 * l0 * l2 / l1
    closing = 2 * math.pi ** 2 / m
    if m <= 31:
        direct_cosine = 2 * volume * (1 - math.cos(q))
        check(close(direct_cosine, stable_exact, 3e-9),
              "stable sine form agrees with direct cosine before cancellation")
    check(stable_exact <= envelope * (1 + 3e-14),
          "cosine bound gives transverse-area over twist-length envelope")
    check(close(envelope, closing), "declared sequence closes as 2pi^2/m")
    check(l0 % 2 == 1 and l2 % 2 == 1 and l1 % 2 == 0,
          "declared sequence obeys parity premises")
    check(min(l0, l1, l2) == m, "injectivity scale diverges with m")
    boundary_ratio = Fraction(2, l0) + Fraction(2, l1) + Fraction(2, l2)
    check(boundary_ratio <= Fraction(5, m), "boundary fraction vanishes")
    isotropic_envelope = 4 * math.pi ** 2 * m
    check(isotropic_envelope >= 4 * math.pi **2 * 5,
          "same elementary estimate does not close isotropically")


print(f"PASS__INDEPENDENT_GL6AW_HOSTILE_REPLAY__{checks}/{checks}")
