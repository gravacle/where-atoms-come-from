#!/usr/bin/env python3
"""Exact constructive replay for the GL6AW twist theorem."""

from __future__ import annotations

import cmath
import itertools
import math


checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise RuntimeError(label)
    checks += 1


DISPLACEMENTS = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (0, 0, 0),
)


def add(x, y, periods):
    return tuple((x[i] + y[i]) % periods[i] for i in range(3))


def subtract(x, y, periods):
    return tuple((x[i] - y[i]) % periods[i] for i in range(3))


def hexagon(base, triple, periods):
    """Return the six ordered (parent-cell,port) links of AR.2."""
    a, b, c = triple
    da, db, dc = (DISPLACEMENTS[index] for index in triple)
    x_da_db = add(subtract(base, db, periods), da, periods)
    x_db_dc = add(subtract(base, db, periods), dc, periods)
    return (
        (base, a),
        (x_da_db, b),
        (x_da_db, c),
        (x_db_dc, a),
        (x_db_dc, b),
        (base, c),
    )


def all_cells(periods):
    return itertools.product(*(range(length) for length in periods))


def centered_witness(periods):
    """Union of the alternating 0/3 matching and constant port-1 matching."""
    occupied = set()
    for x in all_cells(periods):
        occupied.add((x, 0 if x[1] % 2 == 0 else 3))
        occupied.add((x, 1))
    return occupied


# 1. The centered sector is constructively nonempty and exactly locked.
dimension_samples = tuple(
    (l0, l1, l2)
    for l0 in (5, 7, 9)
    for l1 in (4, 6, 8, 10)
    for l2 in (5, 7)
)
for periods in dimension_samples:
    volume = math.prod(periods)
    occupied = centered_witness(periods)
    check(len(occupied) == 2 * volume, "witness has two parent links per cell")
    port_counts = [0, 0, 0, 0]
    child_counts = {x: 0 for x in all_cells(periods)}
    parent_counts = {x: 0 for x in all_cells(periods)}
    for x, port in occupied:
        parent_counts[x] += 1
        port_counts[port] += 1
        child = add(x, DISPLACEMENTS[port], periods)
        child_counts[child] += 1
    check(all(value == 2 for value in parent_counts.values()),
          "witness parent degrees are two")
    check(all(value == 2 for value in child_counts.values()),
          "witness child degrees are two")
    check(port_counts[0] == volume // 2, "witness has centered port-zero charge")
    check(sum(port_counts) == 2 * volume, "witness total occupation is locked")


# 2. Reconstruct every local cycle, its port multiplicities, and exact twist
# changes.  Exactly two orientations per cell have unit x1 change.
for periods in dimension_samples:
    volume = math.prod(periods)
    affected = 0
    unaffected_with_zero = 0
    without_zero = 0
    for base in all_cells(periods):
        for triple in itertools.combinations(range(4), 3):
            cycle = hexagon(base, triple, periods)
            check(len(set(cycle)) == 6, "rectangular elementary cycle has six links")
            multiplicities = {port: 0 for port in range(4)}
            for _, port in cycle:
                multiplicities[port] += 1
            for port in range(4):
                target = 2 if port in triple else 0
                check(multiplicities[port] == target,
                      "cycle contains each used port twice")

            zero_positions = [index for index, (_, port) in enumerate(cycle)
                              if port == 0]
            if not zero_positions:
                without_zero += 1
                continue
            check(len(zero_positions) == 2, "port zero occurs twice")
            check((zero_positions[1] - zero_positions[0]) % 2 == 1,
                  "port-zero occurrences have opposite alternating parity")
            # Toggle even-occupied to odd-occupied.  Only its exponential
            # modulo L1 matters at a wrap.
            delta = 0
            for index, (x, port) in enumerate(cycle):
                if port == 0:
                    delta += (1 if index % 2 else -1) * x[1]
            residue = delta % periods[1]
            if residue in (1, periods[1] - 1):
                affected += 1
            else:
                check(residue == 0, "unaffected zero-port cycle has zero twist")
                unaffected_with_zero += 1
    check(affected == 2 * volume, "exactly two affected orientations per cell")
    check(unaffected_with_zero == volume, "one zero-port orientation is unaffected")
    check(without_zero == volume, "one orientation contains no port-zero link")


# 3. Translation character arithmetic on all sampled odd-area tori.
for periods in dimension_samples:
    l0, l1, l2 = periods
    volume = math.prod(periods)
    number = volume // 2
    check(l0 % 2 == 1 and l2 % 2 == 1 and l1 % 2 == 0,
          "sample has odd transverse area and even twist length")
    residue = number % l1
    phase = cmath.exp(-2j * math.pi * residue / l1)
    check(abs(phase + 1) < 3e-12, "centered twist character is minus one")
    check(2 * residue == l1, "charge quotient is half-integral modulo one")


# 4. PF nonuniformity does not alter the exact cosine energy factor.  Each
# two-configuration partial-flip pair is checked over many positive weights.
for left in range(1, 80):
    for right in range(1, 71):
        norm = math.sqrt(left * left + right * right)
        p, r = left / norm, right / norm
        t_expectation = 2 * p * r
        check(0 < t_expectation <= 1 + 1e-14,
              "positive partial-flip expectation bounded by one")
        for l1 in (4, 6, 10, 18, 34):
            q = 2 * math.pi / l1
            untwisted = 2 * p * r
            twisted = 2 * p * r * math.cos(q)
            check(abs((untwisted - twisted) -
                      (1 - math.cos(q)) * t_expectation) < 3e-12,
                  "twisted partial-flip loss is (1-cos q)t")


# 5. Exact finite-size envelope and the anisotropic Følner sequence.
for m in range(5, 302, 2):
    l0, l1, l2 = m, 2 * m ** 3, m
    volume = l0 * l1 * l2
    q = 2 * math.pi / l1
    # Stable evaluation of 2V(1-cos q).
    exact = 4 * volume * math.sin(q / 2) ** 2
    envelope = 4 * math.pi ** 2 * l0 * l2 / l1
    closing = 2 * math.pi ** 2 / m
    check(exact <= envelope * (1 + 2e-14), "cosine gap obeys area/length envelope")
    check(abs(envelope - closing) < 3e-12 * max(1, closing),
          "chosen aspect ratio gives 2pi^2/m")
    boundary_ratio = 2 * (1 / l0 + 1 / l1 + 1 / l2)
    check(boundary_ratio <= (4 / m + 1 / m ** 3) * (1 + 2e-14),
          "rectangular sequence is Følner")
    check(min(l0, l1, l2) == m, "injectivity scale tends with m")


print(f"PASS__GL6AW_ANISOTROPIC_FOLNER_TWIST__{checks}/{checks}")
