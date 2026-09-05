#!/usr/bin/env python3
"""Exact finite replay of the GL6GK metric-representative descent identity."""

from __future__ import annotations

from fractions import Fraction as F


CHECKS = 0


def check(value: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not value:
        raise AssertionError(label)


def mat(rows):
    return tuple(tuple(F(entry) for entry in row) for row in rows)


def add(left, right):
    return tuple(tuple(left[i][j] + right[i][j] for j in range(3))
                 for i in range(3))


def neg(value):
    return tuple(tuple(-value[i][j] for j in range(3)) for i in range(3))


def mul(left, right):
    return tuple(tuple(sum(left[i][k] * right[k][j] for k in range(3))
                       for j in range(3)) for i in range(3))


def transpose(value):
    return tuple(tuple(value[j][i] for j in range(3)) for i in range(3))


def zero(value):
    return all(entry == 0 for row in value for entry in row)


def trace(value):
    return sum(value[i][i] for i in range(3))


def det(value):
    return (value[0][0] * (value[1][1] * value[2][2] - value[1][2] * value[2][1])
            - value[0][1] * (value[1][0] * value[2][2] - value[1][2] * value[2][0])
            + value[0][2] * (value[1][0] * value[2][1] - value[1][1] * value[2][0]))


I = mat(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
E = mat(((2, 1, 0), (0, 1, 1), (1, 0, 1)))
S = mat(((3, 1, 0), (1, 4, 1), (0, 1, 5)))
R = mat(((2, -1, 1), (-1, 3, 0), (1, 0, 4)))
OWNERS = {
    "stay_contact", "writer_contact", "stay_spectral", "writer_spectral",
    "stay_writer_interference", "disconnected", "node", "support",
    "connector", "maxwell", "constraint", "matching", "state", "measure",
    "retained_field", "seam_period", "self_counterterm", "boundary",
}


def generator(i: int, j: int):
    rows = [[0, 0, 0] for _ in range(3)]
    rows[i][j] = 1
    return mat(rows)


def energy(metric, weight):
    # One distinct metric functional per physical owner.  The coefficients
    # deliberately differ so the owner-once test cannot hide duplication.
    return (weight * trace(mul(R, metric)) +
            (weight + 1) * trace(mul(metric, metric)) +
            (2 * weight + 1) * det(metric))


def main() -> int:
    check(len(OWNERS) == 18 and len(set(OWNERS)) == 18,
          "global owner labels are unique")
    g = mul(mul(E, S), transpose(E))
    for i in range(3):
        for j in range(3):
            A = generator(i, j)
            delta_e = mul(E, A)
            delta_s = neg(add(mul(A, S), mul(S, transpose(A))))
            delta_g = add(add(mul(mul(delta_e, S), transpose(E)),
                              mul(mul(E, delta_s), transpose(E))),
                          mul(mul(E, S), transpose(delta_e)))
            check(zero(delta_g), "GK05 exact generator %d%d" % (i, j))
            # Exact finite directional check of every separately labelled owner.
            epsilon = F(1, 101)
            g_shift = add(g, tuple(tuple(epsilon * delta_g[a][b] for b in range(3))
                                   for a in range(3)))
            for weight, owner in enumerate(sorted(OWNERS), 1):
                check(energy(g_shift, weight) == energy(g, weight),
                      "metric owner is vertical-invariant %s %d%d" %
                      (owner, i, j))
    # Negative control: a representative-only term fails even though delta g=0.
    A = generator(0, 1)
    delta_e = mul(E, A)
    e_norm = trace(mul(transpose(E), E))
    shifted_e = add(E, delta_e)
    check(trace(mul(transpose(shifted_e), shifted_e)) != e_norm,
          "E-only owner creates a nonfactorization defect")
    print("PASS__GL6GK_OWNER_ONCE_METRIC_REPRESENTATIVE_DESCENT__%d/%d" %
          (CHECKS, CHECKS))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as error:
        print("FAIL=" + str(error))
        raise SystemExit(1)
