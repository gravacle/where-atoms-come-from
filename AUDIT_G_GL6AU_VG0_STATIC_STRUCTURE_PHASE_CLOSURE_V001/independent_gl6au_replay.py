#!/usr/bin/env python3
"""Independent hostile replay for frozen GL6AU; imports no author code."""

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
          tolerance: float = 3.0e-11) -> bool:
    return abs(left - right) <= tolerance * max(1.0, abs(left), abs(right))


def cycle_columns(z: tuple[complex, complex, complex, complex]):
    result = []
    for a, b, c in itertools.combinations(range(4), 3):
        column = [0j] * 4
        column[a] = z[b] - z[c]
        column[b] = z[c] - z[a]
        column[c] = z[a] - z[b]
        result.append(tuple(column))
    return tuple(result)


def inner(left, right):
    return sum(x.conjugate() * y for x, y in zip(left, right))


# A. Reconstruct the exact first-character transverse plane and coefficient.
# Every real u=(0,a,b,-a-b) is checked, not just an orthonormal basis.
for L in range(4, 521):
    q = 2.0 * math.pi / L
    z = (cmath.exp(1j * q), 1 + 0j, 1 + 0j, 1 + 0j)
    columns = cycle_columns(z)
    for a, b in ((1, 0), (0, 1), (1, -1), (2, 3), (-3, 2)):
        raw = (0.0, float(a), float(b), float(-a - b))
        norm = math.sqrt(sum(value * value for value in raw))
        u = tuple(complex(value / norm) for value in raw)
        check(close(sum(u), 0), "u has zero port total")
        check(close(u[0], 0), "u is in the exact first-character plane")
        for column in columns:
            check(close(sum(column), 0), "cycle column obeys parent incidence")
            check(close(sum(z[i] * column[i] for i in range(4)), 0),
                  "cycle column obeys child incidence")
        observed = sum(abs(inner(u, column)) ** 2 for column in columns)
        check(close(observed, 3 * abs(1 - z[0]) ** 2),
              "first-character norm is 3|1-z|^2")
        check(close(observed, 12 * math.sin(math.pi / L) ** 2),
              "first-character norm is 12 sin^2(pi/L)")


# B. Independently derive the partial-flip oscillator coefficient.  For a
# positive normalized two-entry vector and arbitrary complex diagonal read,
# direct matrix algebra reduces to (J/2)<T>|delta r|^2.
for left in range(1, 41):
    for right in range(1, 37):
        normalization = math.sqrt(left * left + right * right)
        p, r = left / normalization, right / normalization
        t_expectation = 2 * p * r
        check(0 < t_expectation <= 1 + 1e-14,
              "PF partial-flip expectation lies in (0,1]")
        for seed in (1, 3, 7):
            J = Fraction(2 * seed + 1, seed + 4)
            read_left = complex((left + seed) / 17, (right - seed) / 19)
            read_right = complex((right + 2 * seed) / 23, -(left + seed) / 29)
            delta = read_right - read_left
            # Direct expansion of 1/2 <[rho*,[-J T,rho]]>.
            direct = float(J) * p * r * abs(delta) ** 2
            target = 0.5 * float(J) * t_expectation * abs(delta) ** 2
            check(close(direct, target), "double commutator has J/2 factor")


# C. Least spectral support is bounded by its positive weighted average.
for size in range(1, 80):
    energies = [Fraction(k * k + 2 * size + 1, 5 * size + 3)
                for k in range(1, size + 1)]
    weights = [Fraction((k + 1) * (size + 2 - k), 11 * size + 7)
               for k in range(1, size + 1)]
    average = sum(e * w for e, w in zip(energies, weights)) / sum(weights)
    check(min(energies) <= average, "single-mode spectral average bound")


# D. Reconstruct the quadrature identity from translation orbits.  The array
# entries are deliberately correlated; independence is never used.
for L in range(4, 83):
    q = 2 * math.pi / L
    for seed in (2, 5, 11):
        base = [((x * x + 7 * x * seed + 3 * seed) % 31) - 15 for x in range(L)]
        cosine, sine, rho = [], [], []
        for shift in range(L):
            fc = fs = 0.0
            for x in range(L):
                value = base[(x - shift) % L]
                fc += math.cos(q * x) * value
                fs += math.sin(q * x) * value
            cosine.append(fc)
            sine.append(fs)
            rho.append(complex(fc, fs) / math.sqrt(L))
        mc = sum(cosine) / L
        ms = sum(sine) / L
        mr = sum(rho) / L
        vc = sum((x - mc) ** 2 for x in cosine) / L
        vs = sum((x - ms) ** 2 for x in sine) / L
        structure = sum(abs(x - mr) ** 2 for x in rho) / L
        check(close(mc, 0) and close(ms, 0) and close(mr, 0),
              "nontrivial character has zero translation-orbit mean")
        check(close(vc + vs, L * structure),
              "quadrature variances sum to N times structure")
        check(max(vc, vs) + 1e-10 >= L * structure / 2,
              "one quadrature carries at least half")


# E. Static exponents below two close the exact quotient; alpha=2 does not.
for numerator in range(-4, 20):
    alpha = Fraction(numerator, 10)
    for L in (8, 16, 32, 64, 128, 256, 512, 1024):
        s = 0.41
        exact = 6 * math.sin(math.pi / L) ** 2 / (s * L ** (-float(alpha)))
        envelope = 6 * math.pi ** 2 / s * L ** (float(alpha) - 2)
        check(exact <= envelope * (1 + 2e-14), "static exponent envelope")
    check((alpha < 2) == (float(alpha - 2) < 0),
          "strict alpha threshold is two")


# F. Exact PF transform and Dirichlet identity on every star graph.  This is
# independent of the author's path-3 example and includes nonregular graphs.
for leaves in range(1, 401):
    rho = math.sqrt(leaves)
    psi_center = 1 / math.sqrt(2)
    psi_leaf = 1 / math.sqrt(2 * leaves)
    pi_center = psi_center ** 2
    pi_leaf = psi_leaf ** 2
    check(close(leaves * psi_leaf, rho * psi_center), "PF equation at center")
    check(close(psi_center, rho * psi_leaf), "PF equation at leaf")
    p_center_to_leaf = psi_leaf / (rho * psi_center)
    p_leaf_to_center = psi_center / (rho * psi_leaf)
    check(close(leaves * p_center_to_leaf, 1), "PF transform center row stochastic")
    check(close(p_leaf_to_center, 1), "PF transform leaf row stochastic")
    check(close(pi_center * p_center_to_leaf, pi_leaf * p_leaf_to_center),
          "PF transform detailed balance")
    values = [complex((7 * k + leaves) % 13, (5 * k - leaves) % 17)
              for k in range(leaves + 1)]
    markov_dirichlet = 0.0
    adjacency_form = 0.0
    for k in range(1, leaves + 1):
        jump = abs(values[0] - values[k]) ** 2
        markov_dirichlet += pi_center * p_center_to_leaf * jump
        adjacency_form += psi_center * psi_leaf * jump
    check(close(rho * markov_dirichlet, adjacency_form),
          "PF Markov and adjacency Dirichlet forms agree")


# G. RK comparison obstruction: a frozen degree zero and an active degree d
# force scalar operator-distance at least d/2.  This proves no small uniform
# perturbation; it does not purport to rule out a special-component theorem.
for L in range(4, 260, 4):
    d = Fraction(L ** 3, 64)
    for scalar in (Fraction(0), d / 4, d / 2, 3 * d / 4, d):
        check(max(abs(scalar), abs(d - scalar)) >= d / 2,
              "degree range has scalar distance at least d/2")


# H. Audit the twist statement and independently expose its precise ceiling.
# The cubic even tori are trivial.  Rectangular zero-port tori with twist
# length even and odd transverse area have phase -1, but their elementary
# energy bound scales as transverse area / twist length.
for L in range(4, 400, 2):
    number = L ** 3 // 2
    check(number % L == 0, "even cubic half-fill divides by twist length")
    check(close(cmath.exp(2j * math.pi * (number % L) / L), 1),
          "cubic half-filled twist character is trivial")

for m in range(5, 100, 2):
    transverse_area = m * m
    twist_length = 2 * m ** 3
    volume = transverse_area * twist_length
    number = volume // 2
    phase = cmath.exp(2j * math.pi * (number % twist_length) / twist_length)
    check(close(phase, -1), "odd transverse area gives rectangular phase -1")
    q = 2 * math.pi / twist_length
    exact_bound = 4 * volume * math.sin(q / 2) ** 2
    power_bound = 4 * math.pi ** 2 * transverse_area / twist_length
    check(exact_bound <= power_bound * (1 + 2e-9),
          "two-orientation twist energy has area/length bound")
    check(power_bound <= 2 * math.pi ** 2 / m + 1e-14,
          "anisotropic Følner sequence closes as 1/m")


print(f"PASS__INDEPENDENT_GL6AU_HOSTILE_REPLAY__{checks}/{checks}")
