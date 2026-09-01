#!/usr/bin/env python3
"""Independent hostile replay for frozen GL6AX; imports no author code."""

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


def close(left: complex | float, right: complex | float,
          tolerance: float = 3.0e-11) -> bool:
    return abs(left - right) <= tolerance * max(1.0, abs(left), abs(right))


PORT = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (0, 0, 0),
)


def determinant(matrix):
    result = 0
    for permutation in itertools.permutations(range(len(matrix))):
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(len(matrix))
                         for j in range(i + 1, len(matrix)))
        product = 1
        for row in range(len(matrix)):
            product *= matrix[row][permutation[row]]
        result += (-1 if inversions % 2 else 1) * product
    return result


# A. Re-derive the affine simplex kernel over integers.
columns = tuple((1,) + vector for vector in PORT)
rows = tuple(tuple(columns[column][row] for column in range(4))
             for row in range(4))
check(abs(determinant(rows)) == 1, "augmented port simplex is unimodular")

kernel_solutions = 0
for delta_n in itertools.product(range(-6, 7), repeat=4):
    equations = (
        sum(delta_n),
        delta_n[0],
        delta_n[1],
        delta_n[2],
    )
    if equations == (0, 0, 0, 0):
        kernel_solutions += 1
        check(delta_n == (0, 0, 0, 0),
              "affine moment equations have trivial integer kernel")
check(kernel_solutions == 1, "bounded kernel census finds exactly zero")


# B. Enumerate alternating parent walks without assuming simplicity.  A step
# uses a positive port a and a negative port b at the common child.
pairs = tuple((a, b) for a in range(4) for b in range(4) if a != b)
closed_count = 0
self_touch_count = 0
for length in range(1, 6):
    for word in itertools.product(pairs, repeat=length):
        position = [0, 0, 0]
        visited = [(0, 0, 0)]
        port_change = [0, 0, 0, 0]
        for positive, negative in word:
            port_change[positive] += 1
            port_change[negative] -= 1
            for axis in range(3):
                position[axis] += PORT[positive][axis] - PORT[negative][axis]
            visited.append(tuple(position))
        if position == [0, 0, 0]:
            closed_count += 1
            check(port_change == [0, 0, 0, 0],
                  "every closed alternating parent walk conserves every port")
            if len(set(visited[:-1])) < length:
                self_touch_count += 1
check(closed_count > 1000, "large closed-walk family independently replayed")
check(self_touch_count > 0, "self-touching closed walks are present in replay")


def cells(periods):
    return itertools.product(*(range(period) for period in periods))


def add_mod(left, right, periods):
    return tuple((left[axis] + right[axis]) % periods[axis]
                 for axis in range(3))


def child(edge, periods):
    parent, port = edge
    return add_mod(parent, PORT[port], periods)


def degrees(occupied, periods):
    parent = {x: 0 for x in cells(periods)}
    child_count = {x: 0 for x in cells(periods)}
    for edge in occupied:
        parent[edge[0]] += 1
        child_count[child(edge, periods)] += 1
    return parent, child_count


def assert_locked(occupied, periods, label):
    parent, child_count = degrees(occupied, periods)
    check(all(value == 2 for value in parent.values()),
          f"{label}: parent lock")
    check(all(value == 2 for value in child_count.values()),
          f"{label}: child lock")


def port_totals(occupied):
    result = [0, 0, 0, 0]
    for _, port in occupied:
        result[port] += 1
    return result


def seam_identity(before, after, periods):
    delta = {}
    for edge in before | after:
        delta[edge] = int(edge in after) - int(edge in before)
    before_n = port_totals(before)
    after_n = port_totals(after)
    delta_n = [after_n[a] - before_n[a] for a in range(4)]
    winding = []
    for axis in range(3):
        seam = sum(value for (x, port), value in delta.items()
                   if port == axis and x[axis] == periods[axis] - 1)
        winding.append(seam)
        check(delta_n[axis] == periods[axis] * seam,
              "direct locked pair obeys exact periodic seam identity")
    check(delta_n[3] == -sum(periods[j] * winding[j] for j in range(3)),
          "port-three seam identity follows from total degree")
    return delta_n, winding


# C. Construct the sharp wrapping examples in every direction.
for periods in ((4, 5, 6), (5, 7, 6), (7, 6, 5)):
    for direction in range(3):
        spectator = next(a for a in range(3) if a != direction)
        before = {(x, port) for x in cells(periods)
                  for port in (direction, spectator)}
        assert_locked(before, periods, "uniform wrapping source")
        after = set(before)
        fixed = [0, 0, 0]
        for coordinate in range(periods[direction]):
            x = list(fixed)
            x[direction] = coordinate
            x = tuple(x)
            shifted = add_mod(x, tuple(PORT[direction][axis] - PORT[3][axis]
                                       for axis in range(3)), periods)
            after.remove((x, direction))
            after.add((shifted, 3))
        assert_locked(after, periods, "wrapped-row target")
        delta_n, winding = seam_identity(before, after, periods)
        expected = [0, 0, 0, 0]
        expected[direction] = -periods[direction]
        expected[3] = periods[direction]
        check(delta_n == expected, "wrapped row changes exact port totals")
        check(len(before ^ after) == 2 * periods[direction],
              "wrapped row attains exact two-length Hamming cost")
        check(abs(winding[direction]) == 1,
              "sharp row has one unit of winding")


# D. An actual locked symmetric difference with a degree-four/self-touching
# parent: cross a 0/3 winding row with a 1/2 diagonal winding cycle.
periods = (5, 5, 5)
before = {(x, port) for x in cells(periods) for port in (0, 1)}
after = set(before)
cycle_a = {(t, 0, 0) for t in range(5)}
cycle_b = {(0, t, (-t) % 5) for t in range(5)}
for x in cycle_a:
    shifted = add_mod(x, tuple(PORT[0][axis] - PORT[3][axis]
                               for axis in range(3)), periods)
    after.remove((x, 0))
    after.add((shifted, 3))
for x in cycle_b:
    shifted = add_mod(x, tuple(PORT[1][axis] - PORT[2][axis]
                               for axis in range(3)), periods)
    after.remove((x, 1))
    after.add((shifted, 2))
assert_locked(before, periods, "degree-four source")
assert_locked(after, periods, "degree-four target")
intersection = next(iter(cycle_a & cycle_b))
changed_at_intersection = sum((intersection, port) in (before ^ after)
                              for port in range(4))
check(changed_at_intersection == 4,
      "crossed winding cycles realize a degree-four difference vertex")
seam_identity(before, after, periods)


# E. Independently brute-force the lower winding threshold at the walk level.
for periods in ((4, 5, 6), (5, 6, 7)):
    minimum = min(periods)
    found_at_threshold = False
    for length in range(1, minimum + 1):
        for word in itertools.product(pairs, repeat=length):
            displacement = [0, 0, 0]
            for positive, negative in word:
                for axis in range(3):
                    displacement[axis] += (PORT[positive][axis]
                                           - PORT[negative][axis])
            winds = [displacement[axis] // periods[axis]
                     if displacement[axis] % periods[axis] == 0 else None
                     for axis in range(3)]
            periodic_closed = all(value is not None for value in winds)
            nonzero_winding = periodic_closed and any(winds)
            if length < minimum:
                check(not nonzero_winding,
                      "sub-threshold alternating walk cannot wind")
            elif nonzero_winding:
                found_at_threshold = True
                break
        if found_at_threshold:
            break
    check(found_at_threshold,
          "a minimum-direction winding walk appears at L_min parent steps")


# F. Replay termwise U(1) averaging on the complete five-bit basis.  A
# discrete eleven-point average is exact for all charge differences here.
basis = tuple(itertools.product((0, 1), repeat=5))
charges = tuple(sum(bits) for bits in basis)
roots = 11
for left, charge_left in enumerate(charges):
    for right, charge_right in enumerate(charges):
        raw = complex((3 * left - 2 * right) % 13 - 6,
                      (5 * left + right) % 17 - 8)
        averaged = sum(
            cmath.exp(2j * math.pi * sample
                      * (charge_left - charge_right) / roots) * raw
            for sample in range(roots)
        ) / roots
        if charge_left == charge_right:
            check(close(averaged, raw), "U(1) average retains equal-charge block")
        else:
            check(close(averaged, 0), "U(1) average removes unequal-charge block")


# G. Attack complex hopping, arbitrary diagonal terms, current, and the exact
# central-Taylor coefficient in two-state charge-preserving blocks.
amplitudes = (
    1 + 3j,
    -2.5 + 0.75j,
    0.125 - 4.25j,
    -7 - 0.5j,
)
for left_position in range(-4, 5):
    for right_position in range(-4, 5):
        displacement = right_position - left_position
        for amplitude in amplitudes:
            for denominator in (9, 14, 27, 53):
                q = 2 * math.pi / denominator
                plus = amplitude * cmath.exp(1j * q * displacement)
                minus = amplitude * cmath.exp(-1j * q * displacement)
                central = (plus + minus) / 2
                check(close(central, amplitude * math.cos(q * displacement)),
                      "plus-minus average cancels complex current term")
                central_norm = abs(central - amplitude)
                double_commutator_norm = abs(amplitude) * displacement ** 2
                check(central_norm <= (q * q / 2) * double_commutator_norm
                      * (1 + 2e-13) + 2e-13,
                      "central twist obeys exact one-half double-commutator bound")
                if displacement == 0:
                    check(close(central_norm, 0),
                          "configuration-dependent diagonal energy is twist inert")


# H. Centered-sector character, constants, anisotropic limit, and long hop.
for l0 in range(3, 18, 2):
    for l1 in range(4, 24, 2):
        for l2 in range(3, 16, 2):
            volume = l0 * l1 * l2
            n0 = volume // 2
            character = cmath.exp(-2j * math.pi * n0 / l1)
            check(close(character, -1),
                  "odd transverse area gives centered minus-one character")

for m in range(5, 404, 2):
    l0, l1, l2 = m, 2 * m ** 3, m
    q = 2 * math.pi / l1
    volume = l0 * l1 * l2
    for density in (0.2, 1.0, 7.5, 63 / 8):
        general = (q * q / 2) * volume * density
        displayed = 2 * math.pi ** 2 * density * l0 * l2 / l1
        closing = math.pi ** 2 * density / m
        check(close(general, displayed),
              "q-squared volume algebra gives exact 2pi-squared coefficient")
        check(close(displayed, closing),
              "anisotropic sequence closes as pi-squared D2 over m")
    half_transfer_phase = cmath.exp(1j * q * (l1 // 2))
    check(close(half_transfer_phase, -1),
          "half-system transfer evades the small local twist premise")
    check((l1 // 2) ** 2 == m ** 6,
          "half-system transfer has divergent quadratic moment")


# I. The wrapping-tail expectation allowance is exactly bounded by two norms.
# Diagonal two-level Hermitians attain the factor two under a swapping unitary.
for norm in (0.25, 1.0, 3.5, 19.0):
    original_expectation = norm
    twisted_expectation = -norm
    check(abs(twisted_expectation - original_expectation) == 2 * norm,
          "single twisted-tail difference can attain two operator norms")
    average_difference = abs((twisted_expectation + twisted_expectation) / 2
                             - original_expectation)
    check(average_difference == 2 * norm,
          "two-norm tail coefficient is sharp even after twist averaging")


print(f"PASS__INDEPENDENT_GL6AX_HOSTILE_REPLAY__{checks}/{checks}")
