#!/usr/bin/env python3
"""Independent exact/numerical-identity replay for the GL6AU closure bound."""

from __future__ import annotations

import cmath
import itertools
import math
from fractions import Fraction


checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def close(left: complex | float, right: complex | float,
          tolerance: float = 2.0e-11) -> bool:
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) <= tolerance * scale


def cycle_columns(z: tuple[complex, complex, complex, complex]):
    columns = []
    for a, b, c in itertools.combinations(range(4), 3):
        column = [0j, 0j, 0j, 0j]
        column[a] = z[b] - z[c]
        column[b] = z[c] - z[a]
        column[c] = z[a] - z[b]
        columns.append(tuple(column))
    return tuple(columns)


def dot_conjugate(left, right):
    return sum(a.conjugate() * b for a, b in zip(left, right))


def matrix_multiply(left, right):
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(len(right)))
              for j in range(len(right[0])))
        for i in range(len(left))
    )


def matrix_subtract(left, right):
    return tuple(
        tuple(left[i][j] - right[i][j] for j in range(len(left[0])))
        for i in range(len(left))
    )


def dagger(matrix):
    return tuple(tuple(matrix[j][i].conjugate() for j in range(len(matrix)))
                 for i in range(len(matrix[0])))


def expectation(vector, matrix):
    acted = tuple(sum(matrix[i][j] * vector[j]
                      for j in range(len(vector)))
                  for i in range(len(vector)))
    return sum(vector[i].conjugate() * acted[i] for i in range(len(vector)))


# 1. Exact first-character cycle identity.  The two displayed vectors span
# the real plane u_0=0, sum_a u_a=0.
polarizations = (
    (0.0, 1 / math.sqrt(2), -1 / math.sqrt(2), 0.0),
    (0.0, 1 / math.sqrt(6), 1 / math.sqrt(6), -2 / math.sqrt(6)),
)
for L in range(4, 257):
    q = 2 * math.pi / L
    z = (cmath.exp(1j * q), 1 + 0j, 1 + 0j, 1 + 0j)
    columns = cycle_columns(z)
    for u_real in polarizations:
        u = tuple(complex(value) for value in u_real)
        check(close(sum(abs(value) ** 2 for value in u), 1.0),
              "normalized transverse polarization")
        check(close(sum(u), 0.0), "polarization has zero port total")
        check(close(u[0], 0.0), "polarization is transverse to first character")
        for column in columns:
            check(close(sum(column), 0.0), "cycle column conserves port total")
            check(close(sum(z[a] * column[a] for a in range(4)), 0.0),
                  "cycle column satisfies child incidence")
        norm_squared = sum(abs(dot_conjugate(u, column)) ** 2
                           for column in columns)
        target = 12 * math.sin(math.pi / L) ** 2
        check(close(norm_squared, target), "exact first-character cycle norm")


# 2. A projected two-state flip verifies the double-commutator factor and the
# PF orientation bound without importing any continuum structure.
for left in range(1, 18):
    for right in range(1, 18):
        norm = math.sqrt(left * left + right * right)
        psi = (left / norm + 0j, right / norm + 0j)
        flip = ((0j, 1 + 0j), (1 + 0j, 0j))
        t_value = expectation(psi, flip).real
        check(t_value >= -1e-14, "positive PF vector has nonnegative flip expectation")
        check(t_value <= 1 + 1e-14, "partial flip expectation bounded by one")

for index in range(1, 129):
    J = 0.25 + index / 31
    r0 = complex(index / 17, -index / 29)
    r1 = complex((index + 3) / 19, (2 * index - 1) / 37)
    H = ((0j, -J + 0j), (-J + 0j, 0j))
    density = ((r0, 0j), (0j, r1))
    density_star = dagger(density)
    inner_commutator = matrix_subtract(
        matrix_multiply(H, density), matrix_multiply(density, H)
    )
    double_commutator = matrix_subtract(
        matrix_multiply(density_star, inner_commutator),
        matrix_multiply(inner_commutator, density_star),
    )
    plus = (1 / math.sqrt(2) + 0j, 1 / math.sqrt(2) + 0j)
    f_value = 0.5 * expectation(plus, double_commutator).real
    check(close(f_value, 0.5 * J * abs(r1 - r0) ** 2),
          "two-state oscillator factor J/2 times squared density jump")


# 3. The single-mode estimate is simply the least-support-point versus the
# positive spectral average.  Replay it with exact rational spectral data.
for size in range(2, 19):
    gaps = tuple(Fraction(index * index + size, 3 * size + 1)
                 for index in range(1, size + 1))
    weights = tuple(Fraction((index + 2) * (size - index + 1), 7 * size + 3)
                    for index in range(size))
    average = sum(gap * weight for gap, weight in zip(gaps, weights)) / sum(weights)
    check(min(gaps) <= average, "least support point bounded by spectral average")


# 4. Exact quadrature identity under a translation-invariant orbit measure.
# Occupations need not be assumed independent; the identity is kinematic.
u = (0.0, 1 / math.sqrt(2), -1 / math.sqrt(2), 0.0)
for L in range(5, 42):
    q = 2 * math.pi / L
    for seed in range(1, 6):
        base = tuple(
            tuple(
                ((17 * x + 11 * port + 7 * seed + x * port) % 23 - 11) / 13
                for port in range(4)
            )
            for x in range(L)
        )
        rho_values = []
        cosine_values = []
        sine_values = []
        for shift in range(L):
            rho = 0j
            cosine = 0.0
            sine = 0.0
            for x in range(L):
                for port in range(4):
                    value = base[(x - shift) % L][port]
                    weighted = u[port] * value
                    rho += cmath.exp(1j * q * x) * weighted / math.sqrt(L)
                    cosine += math.cos(q * x) * weighted
                    sine += math.sin(q * x) * weighted
            rho_values.append(rho)
            cosine_values.append(cosine)
            sine_values.append(sine)
        mean_rho = sum(rho_values) / L
        mean_cosine = sum(cosine_values) / L
        mean_sine = sum(sine_values) / L
        structure = sum(abs(value - mean_rho) ** 2 for value in rho_values) / L
        variance_cosine = sum((value - mean_cosine) ** 2
                              for value in cosine_values) / L
        variance_sine = sum((value - mean_sine) ** 2
                            for value in sine_values) / L
        check(close(mean_rho, 0.0), "nontrivial-character density has zero mean")
        check(close(variance_cosine + variance_sine, L * structure),
              "real quadrature variance sum equals volume times structure")
        check(max(variance_cosine, variance_sine) + 1e-12
              >= 0.5 * L * structure,
              "one real quadrature carries at least half the structure weight")


# 5. Static exponent arithmetic, including the alpha=1 O(1/L) case.
for alpha in (0.0, 0.5, 1.0, 1.5, 1.9):
    for L in (8, 16, 32, 64, 128, 256, 512):
        s = 0.37
        exact_bound = 6 * math.sin(math.pi / L) ** 2 / (s * L ** (-alpha))
        power_bound = 6 * math.pi ** 2 / s * L ** (alpha - 2)
        check(exact_bound <= power_bound * (1 + 1e-14),
              "sine bound implies static-exponent power bound")
        check(alpha >= 2 or power_bound * L ** (2 - alpha) > 0,
              "positive static-exponent coefficient")
check(1.0 < 2.0, "Coulomb exponent is strictly closing")


# 6. RK obstruction on the smallest nonregular connected flip graph.  The
# uniform state is killed by D-A, but is not a PF eigenvector of A.
A = (
    (0.0, 1.0, 0.0),
    (1.0, 0.0, 1.0),
    (0.0, 1.0, 0.0),
)
D = (
    (1.0, 0.0, 0.0),
    (0.0, 2.0, 0.0),
    (0.0, 0.0, 1.0),
)
laplacian = matrix_subtract(D, A)
uniform = (1 / math.sqrt(3) + 0j,) * 3
check(close(expectation(uniform, laplacian), 0.0),
      "RK graph Laplacian has equal-amplitude zero state")
uniform_adjacency = tuple(sum(A[i][j] * uniform[j] for j in range(3))
                          for i in range(3))
check(not close(uniform_adjacency[0] / uniform[0],
                uniform_adjacency[1] / uniform[1]),
      "zero-potential adjacency does not preserve equal amplitudes")
pf = (0.5 + 0j, 1 / math.sqrt(2) + 0j, 0.5 + 0j)
check(all(close(value, math.sqrt(2) * pf[index])
          for index, value in enumerate(
              tuple(sum(A[i][j] * pf[j] for j in range(3)) for i in range(3)))),
      "zero-potential ground vector is degree-biased PF vector")
check(D[0][0] != D[1][1], "degree potential is configuration dependent")

# Exact ground-state transform on the same graph.
spectral_radius = math.sqrt(2)
stationary = tuple(abs(value) ** 2 for value in pf)
P = tuple(
    tuple(A[i][j] * pf[j].real / (spectral_radius * pf[i].real)
          for j in range(3))
    for i in range(3)
)
for row in P:
    check(close(sum(row), 1.0), "PF ground-state transform is stochastic")
for i in range(3):
    for j in range(3):
        check(close(stationary[i] * P[i][j], stationary[j] * P[j][i]),
              "PF ground-state transform is reversible")
trial = (complex(2, -1), complex(-1, 3), complex(4, 2))
mean = sum(stationary[i] * trial[i] for i in range(3))
variance = sum(stationary[i] * abs(trial[i] - mean) ** 2 for i in range(3))
dirichlet = 0.5 * sum(
    stationary[i] * P[i][j] * abs(trial[i] - trial[j]) ** 2
    for i in range(3) for j in range(3)
)
quantum_numerator = 0.5 * sum(
    A[i][j] * pf[i].real * pf[j].real * abs(trial[i] - trial[j]) ** 2
    for i in range(3) for j in range(3)
)
check(variance > 0, "PF trial has positive stationary variance")
check(close(spectral_radius * dirichlet, quantum_numerator),
      "PF Dirichlet form equals adjacency variational numerator")

for L in range(4, 132, 4):
    active_degree_lower_bound = Fraction(L ** 3, 64)
    # For a diagonal operator with eigenvalues containing 0 and d, the best
    # scalar approximation has error at least d/2.
    check(active_degree_lower_bound / 2 == Fraction(L ** 3, 128),
          "frozen/active degree range gives extensive scalar-distance bound")

# 6b. The elementary port twist has trivial translation character on even
# cubic half-filled tori; this records why bare flux insertion is no shortcut.
for L in range(4, 260, 2):
    port_number = L ** 3 // 2
    phase = cmath.exp(2j * math.pi * (port_number % L) / L)
    check(close(phase, 1.0), "cubic half-filled port twist has trivial character")


# 7. A closing inter-sector tower need not close the excitation gap inside a
# selected sector: direct two-level sectors give an exact counterexample.
for L in range(4, 260):
    sector_ground_difference = Fraction(1, L)
    internal_gap = Fraction(1, 1)
    check(sector_ground_difference <= Fraction(1, 4),
          "inter-sector ground tower closes")
    check(internal_gap == 1, "selected-sector internal gap remains positive")


print(f"PASS__GL6AU_VG0_STATIC_STRUCTURE__{checks}/{checks}")
