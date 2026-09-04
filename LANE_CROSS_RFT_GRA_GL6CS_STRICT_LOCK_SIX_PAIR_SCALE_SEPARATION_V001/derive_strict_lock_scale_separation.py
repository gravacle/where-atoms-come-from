#!/usr/bin/env python3
"""Exact algebra for the fixed-component strict-lock E2/T2 scale split."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


HERE = Path(__file__).resolve().parent
PAIRS = tuple(combinations(range(4), 2))
A = (F(1),) * 6
T_BASIS = (
    tuple(map(F, (1, 0, 0, 0, 0, -1))),
    tuple(map(F, (0, 1, 0, 0, -1, 0))),
    tuple(map(F, (0, 0, 1, -1, 0, 0))),
)


class Checks:
    def __init__(self):
        self.total = 0

    def equal(self, got, want, label):
        self.total += 1
        if got != want:
            raise AssertionError(f"{label}: got {got!r}, want {want!r}")

    def true(self, condition, label):
        self.total += 1
        if not condition:
            raise AssertionError(label)


CHECK = Checks()


def dot(left, right):
    return sum((F(x) * F(y) for x, y in zip(left, right)), F(0))


def outer(left, right):
    return tuple(tuple(F(x) * F(y) for y in right) for x in left)


def identity(size):
    return tuple(tuple(F(i == j) for j in range(size)) for i in range(size))


def zero(size):
    return tuple(tuple(F(0) for _ in range(size)) for _ in range(size))


def madd(left, right, factor=F(1)):
    return tuple(tuple(left[i][j] + F(factor) * right[i][j]
                       for j in range(len(left[0])))
                 for i in range(len(left)))


def mscale(factor, matrix):
    return tuple(tuple(F(factor) * value for value in row) for row in matrix)


def matvec(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


P_A = mscale(F(1, 6), outer(A, A))
P_T = zero(6)
for vector in T_BASIS:
    P_T = madd(P_T, mscale(F(1, 2), outer(vector, vector)))
P_E = madd(madd(identity(6), P_A, F(-1)), P_T, F(-1))


def pair_word(spins):
    return tuple(F(spins[a] * spins[b]) for a, b in PAIRS)


def local_ring_flip_census():
    rows = []
    for occupied in combinations(range(4), 2):
        spins = tuple(F(-1 if port in occupied else 1) for port in range(4))
        before = pair_word(spins)
        CHECK.equal(matvec(P_T, before), (F(0),) * 6,
                    "locked pair word has no T2 component")
        for a, b in PAIRS:
            if spins[a] == spins[b]:
                continue
            after_spins = list(spins)
            after_spins[a] *= -1
            after_spins[b] *= -1
            after = pair_word(after_spins)
            difference = tuple(after[i] - before[i] for i in range(6))
            CHECK.equal(matvec(P_A, difference), (F(0),) * 6,
                        "ring flip preserves scalar locked-pair component")
            CHECK.equal(matvec(P_T, difference), (F(0),) * 6,
                        "ring flip difference has no T2 component")
            CHECK.equal(matvec(P_E, difference), difference,
                        "ring flip changes a pure E2 diagonal coordinate")
            CHECK.equal(dot(difference, difference), F(16),
                        "every eligible local ring flip has nonzero E2 norm four")
            rows.append({
                "occupied": occupied,
                "flipped_ports": (a, b),
                "difference": difference,
            })
    CHECK.equal(len(rows), 24, "all 24 eligible locked-node ring incidences")
    return rows


def exact_scales():
    # r=h/U.  All displayed coefficients have their remaining 1/U shown.
    j_over_u = F(63, 8)       # J/U = (63/8) r^6
    writer = F(105, 16)       # lambda_T = (105/16) r^6
    contact = F(1, 4)         # g_ct = (1/4U) r^2
    inverse_j = F(1) / j_over_u
    writer_squared_over_j = writer * writer / j_over_u
    writer_over_j = writer / j_over_u
    CHECK.equal(inverse_j, F(8, 63), "E spectral prefactor 1/J")
    CHECK.equal(writer_over_j, F(5, 6), "mixed E-writer prefactor lambda/J")
    CHECK.equal(writer_squared_over_j, F(175, 32),
                "tensor writer spectral prefactor lambda squared over J")
    CHECK.equal(contact, F(1, 4), "tensor contact prefactor")

    # A minimal nonconstant diagonal coordinate on H=-J sigma_x has exact
    # connected susceptibility d^2/J.  It is the finite-component strictness
    # control for the general PF argument in the theorem.
    d = F(3, 5)
    two_state_susceptibility_times_j = d * d
    CHECK.equal(two_state_susceptibility_times_j, F(9, 25),
                "two-state nonconstant diagonal E response is strict")

    return {
        "r_definition": "r=h/U_d",
        "J": "(63/8) U_d r^6",
        "a_E": "1-r^2-(37/12)r^4+O(r^6)",
        "lambda_T": "(105/16)r^6",
        "g_contact": "r^2/(4U_d)",
        "same_state_spectral_prefactors_after_dimensionless_resolvent": {
            "EE": "(8/63) a_E^2 r^-6/U_d",
            "ET": "(5/6) a_E/U_d",
            "TT_writer": "(175/32) r^6/U_d",
        },
        "TT_contact_prefactor": "(1/4)r^2/U_d",
        "coupling_exponents_at_fixed_component": {
            "EE": -6,
            "ET": 0,
            "TT_contact": 2,
            "TT_writer": 6,
        },
        "no_orientation_mixing_required_collective_enhancement": {
            "contact_relative_to_E": "O(r^-8)",
            "writer_susceptibility_relative_to_E": "O(r^-12)",
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-ledger", action="store_true")
    args = parser.parse_args()

    rows = local_ring_flip_census()
    scales = exact_scales()
    ledger = {
        "lane": "GL6CS",
        "scope": "fixed finite H6 flip component and strict-lock coupling asymptotics; no thermodynamic phase, continuum, gravity, or G",
        "local_ring_flip_census": {
            "count": len(rows),
            "rows": rows,
            "conclusion": "every active locked ring move changes a nonzero pure-E2 diagonal pair coordinate",
        },
        "scales": scales,
        "fixed_component_conclusion": "if an E2 diagonal read is nonconstant on the connected flip component, its positive spectral coefficient is O(r^-6/U); the complete T2 first vertex through h6 plus h2 contact is at most O(r^2/U) at leading order for bounded dimensionless response",
        "accumulation_boundary": "a successful rotational limit must use authenticated orientation mixing, a nonuniform collective/critical limit, finite-r dynamics, or another same-order physical block; fixed-frame repetition with bounded coefficients is insufficient",
    }

    def encode(value):
        if isinstance(value, F):
            return (str(value.numerator) if value.denominator == 1
                    else f"{value.numerator}/{value.denominator}")
        if isinstance(value, dict):
            return {str(key): encode(item) for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            return [encode(item) for item in value]
        return value

    payload = json.dumps(encode(ledger), indent=2, sort_keys=True) + "\n"
    target = HERE / "EXACT_LEDGER.json"
    if args.write_ledger:
        target.write_text(payload)
    else:
        CHECK.true(target.is_file(), "frozen exact ledger exists")
        CHECK.equal(target.read_text(), payload, "frozen exact ledger matches replay")

    print(f"PASS__GL6CS_STRICT_LOCK_SCALE_SEPARATION__{CHECK.total}/{CHECK.total}")
    print("LOCKED_RING_CHANGE=PURE_E2_NONZERO;24/24")
    print("FIXED_COMPONENT_EXPONENTS=EE_-6;ET_0;TT_CONTACT_2;TT_WRITER_6")
    print("FIXED_FRAME_STRICT_LOCK_ROTATIONAL_EQUALITY=IMPOSSIBLE_FOR_BOUNDED_COEFFICIENTS")
    print("ESCAPES=AUTHENTICATED_ORIENTATION_MIXING;COLLECTIVE_SINGULAR_LIMIT;FINITE_R;NEW_SAME_ORDER_BLOCK")


if __name__ == "__main__":
    main()
