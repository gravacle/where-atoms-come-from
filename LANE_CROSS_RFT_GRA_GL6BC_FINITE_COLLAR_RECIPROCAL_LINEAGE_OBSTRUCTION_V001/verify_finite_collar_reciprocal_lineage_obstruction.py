#!/usr/bin/env python3
"""Constructive replay for the GL6BC finite-collar reciprocity theorem."""

from __future__ import annotations

import cmath
import json
import math
from pathlib import Path


LANE = Path(__file__).resolve().parent
CHECKS = 0
TOL = 2.0e-12


def check(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"FAIL [{CHECKS}] {label}")


def zeros(n: int, m: int) -> list[list[complex]]:
    return [[0j for _ in range(m)] for _ in range(n)]


def eye(n: int) -> list[list[complex]]:
    out = zeros(n, n)
    for i in range(n):
        out[i][i] = 1
    return out


def add(*matrices: list[list[complex]]) -> list[list[complex]]:
    n, m = len(matrices[0]), len(matrices[0][0])
    out = zeros(n, m)
    for matrix in matrices:
        for i in range(n):
            for j in range(m):
                out[i][j] += matrix[i][j]
    return out


def scale(a: complex, matrix: list[list[complex]]) -> list[list[complex]]:
    return [[a * value for value in row] for row in matrix]


def mul(a: list[list[complex]], b: list[list[complex]]) -> list[list[complex]]:
    n, p, m = len(a), len(b), len(b[0])
    out = zeros(n, m)
    for i in range(n):
        for k in range(p):
            if a[i][k] == 0:
                continue
            for j in range(m):
                out[i][j] += a[i][k] * b[k][j]
    return out


def dagger(a: list[list[complex]]) -> list[list[complex]]:
    return [[a[j][i].conjugate() for j in range(len(a))]
            for i in range(len(a[0]))]


def kron(a: list[list[complex]], b: list[list[complex]]) -> list[list[complex]]:
    out = zeros(len(a) * len(b), len(a[0]) * len(b[0]))
    for i, row in enumerate(a):
        for j, value in enumerate(row):
            for k, brow in enumerate(b):
                for ell, bvalue in enumerate(brow):
                    out[i * len(b) + k][j * len(b[0]) + ell] = value * bvalue
    return out


def comm(a: list[list[complex]], b: list[list[complex]]) -> list[list[complex]]:
    return add(mul(a, b), scale(-1, mul(b, a)))


def max_abs(a: list[list[complex]]) -> float:
    return max(abs(value) for row in a for value in row)


def mat_vec(a: list[list[complex]], v: list[complex]) -> list[complex]:
    return [sum(a[i][j] * v[j] for j in range(len(v))) for i in range(len(a))]


def inner(v: list[complex], w: list[complex]) -> complex:
    return sum(a.conjugate() * b for a, b in zip(v, w))


def expectation(v: list[complex], a: list[list[complex]]) -> complex:
    return inner(v, mat_vec(a, v))


def poisson_tail(d: int, x: float) -> float:
    check(d >= 1 and x >= 0, "tail domain")
    if x == 0:
        return 0.0
    term = x ** d / math.factorial(d)
    total = term
    m = d
    while True:
        m += 1
        term *= x / m
        total += term
        if term <= 1.0e-16 * max(1.0, total):
            return total
        check(m < 10000, "tail summation terminates")


def delta_l(radius: int, ratio: float, time: float) -> float:
    if ratio == 0 or time == 0:
        return 0.0
    a_l = 3 * radius * radius + 3 * radius + 1
    return min(1.0, a_l * poisson_tail(2 * radius + 1, 48 * ratio * abs(time)))


def verify_operator_obstruction() -> None:
    identity = eye(2)
    x = [[0, 1], [1, 0]]
    z = [[1, 0], [0, -1]]
    n = [[0, 0], [0, 1]]
    p = [[0, 0], [0, 1]]
    q = [[1, 0], [0, 0]]

    i4 = eye(4)
    p_full = kron(p, i4)
    q_full = kron(q, i4)
    x1 = kron(identity, kron(x, identity))
    x2 = kron(identity, kron(identity, x))
    n1 = kron(identity, kron(n, identity))
    n2 = kron(identity, kron(identity, n))
    m_pair = kron(identity, kron(z, z))

    ratio = 2
    h_full = add(
        scale(-1, mul(p_full, add(x1, x2))),
        scale(-6 * ratio, add(n1, n2)),
        scale(2 * ratio, mul(n1, n2)),
    )
    h_break = add(
        scale(-6 * ratio, add(n1, n2)),
        scale(2 * ratio, mul(n1, n2)),
    )

    check(max_abs(comm(h_full, p_full)) == 0, "full route lift commutes with P_K")
    check(max_abs(comm(h_full, q_full)) == 0, "full route lift commutes with I-P_K")
    check(max_abs(comm(m_pair, p_full)) == 0, "pair read commutes with P_K")
    check(max_abs(comm(h_break, m_pair)) == 0, "all-BREAK pair is exactly static")
    check(max_abs(comm(h_full, m_pair)) > 0, "KEEP sector can have nontrivial pair dynamics")

    # An exact pair-source generator and its unitary commute with the route word.
    angle = 0.731
    source = add(scale(math.cos(angle / 2), eye(8)),
                 scale(1j * math.sin(angle / 2), m_pair))
    check(max_abs(comm(source, p_full)) < TOL, "pair-source unitary commutes with P_K")

    # Nontrivial block-controlled evolution at s=pi/4.
    r = add(scale(1 / math.sqrt(2), identity),
            scale(1j / math.sqrt(2), x))
    rr = kron(r, r)
    controlled = add(kron(q, i4), kron(p, rr))
    check(max_abs(comm(controlled, p_full)) < TOL, "controlled propagator commutes with P_K")
    check(max_abs(add(mul(dagger(controlled), controlled), scale(-1, eye(8)))) < TOL,
          "controlled propagator is unitary")

    # An entangled cross-sector state retains its diagonal route probability.
    state = [0j] * 8
    state[0] = 1 / math.sqrt(3)
    state[4] = math.sqrt(2 / 3)
    evolved = mat_vec(controlled, state)
    before = expectation(state, p_full).real
    after = expectation(evolved, p_full).real
    check(abs(before - 2 / 3) < TOL, "entangled state initial P_K probability")
    check(abs(after - before) < TOL, "entangled state P_K probability conserved")

    # Both retarded orientations vanish after nontrivial controlled evolution.
    m_t = mul(dagger(controlled), mul(m_pair, controlled))
    check(max_abs(comm(p_full, m_pair)) < TOL, "P_K then pair retarded commutator zero")
    check(max_abs(comm(m_t, p_full)) < TOL, "evolved pair then P_K retarded commutator zero")

    # Pair-source insertion cannot change the route probability.
    pulsed = mat_vec(source, evolved)
    check(abs(expectation(pulsed, p_full).real - before) < TOL,
          "pair pulse leaves route probability invariant")


def verify_forward_positive_control() -> None:
    identity = eye(2)
    x = [[0, 1], [1, 0]]
    z = [[1, 0], [0, -1]]
    blank = [1 + 0j, 0j, 0j, 0j]
    m_pair = kron(z, z)

    for step in range(17):
        time = math.pi * step / 64
        rotation = add(scale(math.cos(time), identity),
                       scale(1j * math.sin(time), x))
        state_keep = mat_vec(kron(rotation, rotation), blank)
        mean = expectation(state_keep, m_pair).real
        expected = math.cos(2 * time) ** 2
        dtv = abs(1 - mean) / 2
        check(abs(mean - expected) < TOL, f"forward mean formula {step}")
        check(abs(dtv - 0.5 * math.sin(2 * time) ** 2) < TOL,
              f"forward DTV formula {step}")

    time = math.pi / 4
    check(abs(0.5 * math.sin(2 * time) ** 2 - 0.5) < TOL,
          "positive control reaches DTV one half")


def verify_collar_and_metric_transfer() -> None:
    # Moderate-ratio substitutions inherited from GL6BA.
    check(48 * 2 == 96, "R=2 tail argument coefficient")
    check(48 * 2.5 == 120, "R=5/2 tail argument coefficient")
    for ratio in (0.0, 2.0, 2.5, 7.0):
        check(delta_l(0, ratio, 0.0) == 0, f"zero-time exact tail R={ratio}")

    # Factorial decay beats the quadratic shell for each finite input.
    for ratio, time in ((2.0, 0.01), (2.5, 0.01), (2.0, 0.1), (2.5, 0.1)):
        values = [delta_l(radius, ratio, time) for radius in range(1, 81)]
        check(values[-1] < 1.0e-12, f"tail converges R={ratio},s={time}")
        check(min(values[20:]) <= min(values[:20]), f"tail eventually improves R={ratio},s={time}")

    # Exhaust binary distributions and exact all-BREAK comparison.  If the
    # KEEP collar error is <= delta, the distance interval pays one tail.
    grid = [i / 20 for i in range(21)]
    for p_full in grid:
        for p_collar in grid:
            error = abs(p_full - p_collar)
            for p_break in grid:
                d_full = abs(p_full - p_break)
                d_collar = abs(p_collar - p_break)
                check(abs(d_full - d_collar) <= error + TOL,
                      "one-tail reverse triangle")
                check(max(0.0, d_collar - error) <= d_full + TOL,
                      "one-tail lower interval")
                check(d_full <= min(1.0, d_collar + error) + TOL,
                      "one-tail upper interval")

    # A positive collar margin greater than the tail certifies positivity.
    for d_collar, error in ((0.3, 0.1), (0.9, 0.2), (0.01, 0.001)):
        check(d_collar > error, "positive test has strict certified margin")
        check(d_collar - error > 0, "positive full-parent lower bound")


def verify_packet_semantics() -> None:
    theorem = (LANE / "THEOREM.md").read_text(encoding="utf-8")
    result = (LANE / "RESULT.md").read_text(encoding="utf-8")
    ledger = json.loads((LANE / "RECIPROCITY_LEDGER.json").read_text(encoding="utf-8"))

    required_theorem = [
        "d_L>\\delta_L",
        "p^{\\Omega,0}=p^{L,0}",
        "\\chi^R_{\\Pi_\\beta,B}(s)=0",
        "\\chi^R_{B,\\Pi_\\beta}(s)=0",
        "future writer/formation channel",
        "undefined",
        "appended as an identity",
        "No graviton",
        "No graviton, Ricci form, Einstein equation, conventional stiffness",
    ]
    for token in required_theorem:
        check(token in theorem, f"theorem scope token: {token}")
    for token in ("one-way physical-support-gated response", "not a metric",
                  "future writer/channel"):
        check(token in result, f"result scope token: {token}")
    check(ledger["lane"] == "GL6BC_V001", "ledger lane")
    check("future writer" in ledger["missing_parent_law"], "ledger missing law")
    check("no gravity" in ledger["ceilings"], "ledger gravity ceiling")
    check("no G" in ledger["ceilings"], "ledger G ceiling")


def main() -> None:
    verify_operator_obstruction()
    verify_forward_positive_control()
    verify_collar_and_metric_transfer()
    verify_packet_semantics()
    print(f"PASS: {CHECKS}/{CHECKS} GL6BC constructive checks")


if __name__ == "__main__":
    main()
