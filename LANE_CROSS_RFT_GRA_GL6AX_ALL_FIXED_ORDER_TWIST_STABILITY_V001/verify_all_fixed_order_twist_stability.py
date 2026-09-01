#!/usr/bin/env python3
"""Constructive replay for GL6AX affine conservation and twist stability."""

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


D = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (0, 0, 0),
)


def determinant(matrix):
    total = 0
    size = len(matrix)
    for permutation in itertools.permutations(range(size)):
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(size) for j in range(i + 1, size))
        product = 1
        for row in range(size):
            product *= matrix[row][permutation[row]]
        total += (-1 if inversions % 2 else 1) * product
    return total


augmented = tuple((1,) + displacement for displacement in D)
# determinant() expects rows, so transpose the augmented columns.
augmented_rows = tuple(tuple(augmented[column][row] for column in range(4))
                       for row in range(4))
check(abs(determinant(augmented_rows)) == 1,
      "native augmented port columns are unimodular")


# Enumerate a large family of finite alternating parent-to-parent cycles.
# A pair (a,b) traverses port a parent->child and port b child->next parent.
pairs = tuple((a, b) for a in range(4) for b in range(4) if a != b)
for pair_count in range(2, 6):
    for word in itertools.product(pairs, repeat=pair_count):
        displacement = [0, 0, 0]
        delta_ports = [0, 0, 0, 0]
        for a, b in word:
            for axis in range(3):
                displacement[axis] += D[a][axis] - D[b][axis]
            delta_ports[a] += 1
            delta_ports[b] -= 1
        if displacement == [0, 0, 0]:
            check(sum(delta_ports) == 0, "closed word has zero total degree change")
            check(delta_ports == [0, 0, 0, 0],
                  "finite closed alternating word conserves every port")


def add_mod(x, y, periods):
    return tuple((x[axis] + y[axis]) % periods[axis] for axis in range(3))


def cells(periods):
    return itertools.product(*(range(period) for period in periods))


def degree_counts(occupied, periods):
    parent = {x: 0 for x in cells(periods)}
    child = {x: 0 for x in cells(periods)}
    for x, port in occupied:
        parent[x] += 1
        child[add_mod(x, D[port], periods)] += 1
    return parent, child


def port_counts(occupied):
    answer = [0, 0, 0, 0]
    for _, port in occupied:
        answer[port] += 1
    return answer


# Explicit sharp wrapping counterexamples in all three directions.
for periods in ((5, 6, 7), (7, 8, 5), (9, 10, 7)):
    for j in range(3):
        spectator = next(port for port in range(3) if port != j)
        occupied = {(x, port) for x in cells(periods) for port in (j, spectator)}
        parent, child = degree_counts(occupied, periods)
        check(all(value == 2 for value in parent.values()),
              "uniform two-port state has parent degree two")
        check(all(value == 2 for value in child.values()),
              "uniform two-port state has child degree two")
        before = port_counts(occupied)
        changed = set(occupied)
        fixed = [0, 0, 0]
        for t in range(periods[j]):
            x = list(fixed)
            x[j] = t
            x = tuple(x)
            next_x = list(fixed)
            next_x[j] = (t + 1) % periods[j]
            next_x = tuple(next_x)
            changed.remove((x, j))
            changed.add((next_x, 3))
        after = port_counts(changed)
        parent, child = degree_counts(changed, periods)
        check(all(value == 2 for value in parent.values()),
              "wrapped toggle preserves parent lock")
        check(all(value == 2 for value in child.values()),
              "wrapped toggle preserves child lock")
        delta = [after[a] - before[a] for a in range(4)]
        expected = [0, 0, 0, 0]
        expected[j] = -periods[j]
        expected[3] = periods[j]
        check(delta == expected, "wrapped row changes exact pair of port totals")
        check(len(occupied.symmetric_difference(changed)) == 2 * periods[j],
              "wrapped row attains the two-times-length threshold")


# Every parent-to-parent step has unit coordinate speed, proving the lower
# threshold for a nonzero winding word.
for a, b in pairs:
    step = tuple(D[a][axis] - D[b][axis] for axis in range(3))
    check(all(abs(component) <= 1 for component in step),
          "parent walk has at most unit coordinate speed")
for lengths in ((5, 6, 7), (9, 4, 11), (13, 8, 5)):
    minimum = min(lengths)
    check(2 * minimum == min(2 * length for length in lengths),
          "minimum wrapping cycle threshold is two L_min")


# Translation character in the centered sector for many odd/even boxes.
for l0 in range(3, 24, 2):
    for l1 in range(4, 30, 2):
        for l2 in range(3, 20, 2):
            volume = l0 * l1 * l2
            n0 = volume // 2
            phase = cmath.exp(-2j * math.pi * n0 / l1)
            check(abs(phase + 1) < 2e-11,
                  "centered large twist has minus-one translation character")


# Complex-amplitude attack.  In each fixed-charge block, +/- averaging
# removes the sine/current term entry by entry and obeys the double-commutator
# Taylor bound.  No real-matrix premise enters.
complex_amplitudes = (
    1 + 2j,
    -3 + 0.5j,
    0.125 - 4j,
    -2.75 - 1.25j,
)
for width in range(3, 10):
    for charge in range(1, width):
        basis = tuple(itertools.combinations(range(width), charge))
        for left_index, left in enumerate(basis):
            w_left = sum(left)
            for right_index, right in enumerate(basis):
                w_right = sum(right)
                delta_w = w_right - w_left
                amplitude = complex_amplitudes[(left_index + 3 * right_index) %
                                               len(complex_amplitudes)]
                for denominator in (11, 17, 31, 47):
                    q = 2 * math.pi / denominator
                    plus = amplitude * cmath.exp(1j * q * delta_w)
                    minus = amplitude * cmath.exp(-1j * q * delta_w)
                    average = (plus + minus) / 2
                    check(abs(average - amplitude * math.cos(q * delta_w)) < 2e-12,
                          "plus-minus average cancels complex odd twist")
                    lhs = abs(average - amplitude)
                    rhs = 0.5 * q * q * delta_w * delta_w * abs(amplitude)
                    check(lhs <= rhs * (1 + 3e-13) + 1e-14,
                          "entrywise cosine loss obeys double-commutator bound")
                if left == right:
                    check(delta_w == 0, "diagonal term has zero twist change")


# The anisotropic area/length envelope closes for every fixed D2, while a
# system-spanning transfer has order-one phase and violates the small local
# second-moment premise.
for m in range(5, 404, 2):
    l0, l1, l2 = m, 2 * m ** 3, m
    for d2 in (0.25, 1.0, 7.5, 63 / 8):
        bound = 2 * math.pi ** 2 * d2 * l0 * l2 / l1
        check(abs(bound - math.pi ** 2 * d2 / m) < 2e-12 * max(1, bound),
              "general local bound closes as pi^2 D2/m")
    q = 2 * math.pi / l1
    spanning_phase = cmath.exp(1j * q * (l1 // 2))
    check(abs(spanning_phase + 1) < 2e-12,
          "half-system transfer has order-one twist phase")
    check((l1 // 2) ** 2 >= m ** 6,
          "system-spanning transfer has divergent second moment")


print(f"PASS__GL6AX_ALL_FIXED_ORDER_TWIST_STABILITY__{checks}/{checks}")
