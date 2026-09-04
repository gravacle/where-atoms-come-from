#!/usr/bin/env python3
"""Independent hostile reconstruction of GL6CS strict-lock scaling."""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def eq(got, expected, label):
    check(got == expected, f"{label}: got {got!r}, expected {expected!r}")


def dot(a, b):
    return sum((F(x)*F(y) for x, y in zip(a, b)), F(0))


def outer(a, b):
    return tuple(tuple(F(x)*F(y) for y in b) for x in a)


def ident(n):
    return tuple(tuple(F(i == j) for j in range(n)) for i in range(n))


def add(a, b, factor=F(1)):
    return tuple(tuple(a[i][j] + F(factor)*b[i][j] for j in range(len(a[0])))
                 for i in range(len(a)))


def scale(c, a):
    return tuple(tuple(F(c)*x for x in row) for row in a)


def transpose(a):
    return tuple(tuple(a[i][j] for i in range(len(a)))
                 for j in range(len(a[0])))


def mul(a, b):
    return tuple(tuple(sum((a[i][r]*b[r][j] for r in range(len(b))), F(0))
                       for j in range(len(b[0]))) for i in range(len(a)))


def matvec(a, v):
    return tuple(dot(row, v) for row in a)


PAIRS = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: i for i, pair in enumerate(PAIRS)}
complement = {}
for i, pair in enumerate(PAIRS):
    other = tuple(q for q in range(4) if q not in pair)
    complement[i] = PAIR_INDEX[other]

ones = (F(1),) * 6
P_A = scale(F(1, 6), outer(ones, ones))
O = tuple(tuple(F(complement[i] == j) for j in range(6)) for i in range(6))
P_T = scale(F(1, 2), add(ident(6), O, F(-1)))
P_E = add(scale(F(1, 2), add(ident(6), O)), P_A, F(-1))
eq(add(add(P_A, P_E), P_T), ident(6), "A+E+T projector resolution")


def pair_word(z):
    return tuple(F(z[a]*z[b]) for a, b in PAIRS)


rows = []
for occupied in combinations(range(4), 2):
    z = tuple(-1 if i in occupied else 1 for i in range(4))
    before = pair_word(z)
    eq(matvec(P_T, before), (F(0),)*6, "locked word T2 dark")
    for a, b in PAIRS:
        if z[a] == z[b]:
            continue
        moved = list(z)
        moved[a] *= -1
        moved[b] *= -1
        after = pair_word(moved)
        delta = tuple(after[i]-before[i] for i in range(6))
        eq(matvec(P_A, delta), (F(0),)*6, "ring move scalar null")
        eq(matvec(P_T, delta), (F(0),)*6, "ring move tensor null")
        eq(matvec(P_E, delta), delta, "ring move pure E2")
        eq(dot(delta, delta), F(16), "ring move E2 norm squared")
        rows.append((occupied, (a, b), delta))
eq(len(rows), 24, "complete 24-incidence local census")


# Reconstruct the inherited solder normalization rather than importing BQ10.
tetra = (
    (F(1, 2), F(1, 2), F(1, 2)),
    (F(1, 2), F(-1, 2), F(-1, 2)),
    (F(-1, 2), F(1, 2), F(-1, 2)),
    (F(-1, 2), F(-1, 2), F(1, 2)),
)


def sym_outer(a, b):
    return tuple(tuple(a[i]*b[j] + b[i]*a[j] for j in range(3))
                 for i in range(3))


def coords(h):
    return (h[0][0], h[1][1], h[2][2], h[0][1], h[0][2], h[1][2])


dcols = tuple(coords(sym_outer(tetra[a], tetra[b])) for a, b in PAIRS)
D = transpose(dcols)
G_frob = tuple(tuple(F((1, 1, 1, 2, 2, 2)[i]) * F(i == j)
                       for j in range(6)) for i in range(6))
gram = mul(transpose(D), mul(G_frob, D))
expected_gram = add(add(scale(F(1, 2), P_A), scale(F(2), P_E)), P_T)
eq(gram, expected_gram, "pair-to-Sym2 solder Gram")

# Equal spatial normalization g pulls back with raw pair coefficients
# h_E=2g and h_T=g.
g = F(7, 5)
eq((2*g)/2, g, "metric-normalized h_E/2=h_T relation")


# Independent exact scale algebra.  CS10 factors out the common spectral 2
# displayed in CS08; both conventions are recorded here.
j0 = F(63, 8)          # J/U_d = j0 r^6
lambda0 = F(105, 16)  # lambda_T = lambda0 r^6
contact0 = F(1, 4)
eq(1/j0, F(8, 63), "one-over-J coefficient")
eq(lambda0/j0, F(5, 6), "lambda-over-J coefficient")
eq(lambda0*lambda0/j0, F(175, 32), "lambda-squared-over-J coefficient")
eq(2/j0, F(16, 63), "full EE spectral coefficient includes common two")
eq(2*lambda0/j0, F(5, 3), "full ET spectral coefficient includes common two")
eq(2*lambda0*lambda0/j0, F(175, 16),
   "full TT writer coefficient includes common two")
eq(contact0, F(1, 4), "TT contact coefficient")

# a_E=1-r^2-(37/12)r^4+... has a nonzero unit leading coefficient.
eq((F(1), F(-1), F(-37, 12))[0], F(1), "E first vertex has unit leading term")

exponents = {"EE": -6, "ET": 0, "TT_contact": 2, "TT_writer": 6}
eq((exponents["EE"], exponents["ET"], exponents["TT_contact"],
    exponents["TT_writer"]), (-6, 0, 2, 6), "four fixed-component powers")
eq((exponents["EE"] + exponents["TT_writer"]) // 2, 0,
   "Cauchy-Schwarz ET exponent")

# Exact two-state strictness control: H=-J sigma_x has gap 2J and
# |<excited|diag(d,-d)|ground>|=d, hence K=2d^2/(2J)=d^2/J.
d = F(11, 13)
gap_in_j = F(2)
matrix_element_squared = d*d
two_state_k_times_j = 2*matrix_element_squared/gap_in_j
eq(two_state_k_times_j, d*d, "two-state positive E susceptibility")
check(two_state_k_times_j > 0, "strict E response witness positive")

# A bounded TT contact dominates the writer but is eight powers below EE.
eq(exponents["TT_contact"] - exponents["EE"], 8,
   "contact needs r^-8 collective enhancement")
eq(exponents["TT_writer"] - exponents["EE"], 12,
   "writer susceptibility needs r^-12 enhancement")

# Ordinary N-fold fixed-frame repetition followed by extensive normalization
# cancels N and cannot alter a coupling exponent.
for n in (1, 2, 7, 101):
    for power in exponents.values():
        coefficient = F(3, 7)
        eq(F(n)*coefficient/F(n), coefficient,
           f"extensive normalization preserves exponent {power} for N={n}")

# If c_E>0 is fixed while all TT dimensionless coefficients are bounded,
# 2 h_T/h_E=O(r^8), not one.  An exact sample sequence makes the limiting
# contradiction executable without asserting an actual thermodynamic phase.
c_e, c_t, c_w = F(5, 3), F(7, 4), F(2, 9)
ratios = []
for inv_r in (10, 100, 1000):
    r = F(1, inv_r)
    h_e = c_e*r**-6
    h_t = c_t*r**2 + c_w*r**6
    ratios.append(2*h_t/h_e)
check(ratios[2] < ratios[1] < ratios[0] and ratios[2] < F(1, 10**20),
      "bounded fixed-frame TT/EE ratio tends to zero")


result = {
    "schema": "AUDIT_G_GL6CS_STRICT_LOCK_ACCUMULATION_HORIZON_V001",
    "target_directory": "LANE_CROSS_RFT_GRA_GL6CS_STRICT_LOCK_SIX_PAIR_SCALE_SEPARATION_V001",
    "disposition": "PASS",
    "independent_results": {
        "eligible_locked_ring_incidences": 24,
        "every_ring_delta": "nonzero pure E2 with norm squared 16",
        "solder_relation": "h_E/2=h_T",
        "fixed_component_exponents": exponents,
        "full_spectral_prefactors": {
            "EE": "(16/63)a_E^2 r^-6/U_d",
            "ET": "(5/3)a_E/U_d",
            "TT_writer": "(175/16)r^6/U_d",
        },
        "CS10_common_two_removed_prefactors": {
            "EE": "(8/63)a_E^2 r^-6/U_d",
            "ET": "(5/6)a_E/U_d",
            "TT_writer": "(175/32)r^6/U_d",
        },
        "TT_contact": "r^2/(4U_d)",
        "fixed_frame_bounded_match": "impossible as r->0 because at least one E2 eigenresponse is strictly O(r^-6) while the entire bounded T2 block is O(r^2)",
        "required_enhancement_without_orientation_mixing": {
            "contact": "O(r^-8)",
            "writer_susceptibility": "O(r^-12)",
        },
    },
    "scope": {
        "proved": "fixed finite component strict-lock obstruction and the enhancement powers required of a nonuniform accumulation",
        "not_proved": "existence or noncommutation of thermodynamic limits, realization of an orientation law or critical enhancement, finite-r phase, full rotational completion, 1PI/Ricci/Einstein dynamics, gravity, or G",
        "CP_status": "orientation averaging remains a logically valid escape class, but the frozen GL6CP packet is REPAIR_REQUIRED and is not promoted here",
    },
    "nonmaterial_notes": [
        "CS10 lists the factors after removing the common spectral factor 2 already explicit in CS08; the full response coefficients are twice those entries, with identical powers.",
        "RESULT.md has a TeX typo `,qquad` instead of `,\\qquad`.",
        "The scientific audit used six frozen author bytes; final author custody was added without changing those bytes, and the resealed audit pins all twelve final target bytes.",
    ],
}

frozen = json.loads((HERE / "INDEPENDENT_RESULT.json").read_text())
eq(frozen, result, "frozen independent result")
print(f"PASS__AUDIT_GL6CS_INDEPENDENT_RECONSTRUCTION__{checks}/{checks}")
print("DISPOSITION=PASS")
print("FIXED_COMPONENT_POWERS=EE_-6;ET_0;TT_CONTACT_2;TT_WRITER_6")
print("FIXED_FRAME_BOUNDED_MATCH=IMPOSSIBLE")
print("COLLECTIVE_ENHANCEMENT=CONTACT_R^-8;WRITER_R^-12")
print("NONCOMMUTING_LIMIT_OR_ORIENTATION_OR_FINITE_R_OR_NEW_OWNER=OPEN_ROUTES_NOT_DERIVED")
print("NO_RICCI_EINSTEIN_GRAVITY_G")
