#!/usr/bin/env python3
"""Exact complete-H6 diagonal source and homogeneous response on FO.

The certificate enumerates every irreducible diagonal closed flip word of
lengths two, four, and six; differentiates the hopping numerators and the
endpoint-referenced Coulomb-DPAR gaps; restores every Brillouin-Wigner /
Feshbach self-consistency fold through order six; and composes the result
with the independently audited FW pair-plus-ring response.

This is a finite, zero-momentum, fixed-through-H6 theorem under FV-PURE.  It
does not infer a convergence radius, a thermodynamic Ward identity, a
massless tensor pole, RGRL-B, gravity, or Newton's constant.
"""

from collections import Counter
from contextlib import redirect_stdout
from fractions import Fraction as F
from hashlib import sha256
from io import StringIO
from itertools import combinations, permutations
from pathlib import Path
import runpy

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
LANE = Path(__file__).resolve().parent
FO_SCRIPT = (ROOT / "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001" /
             "verify_finite_tt_four_point.py")
FW_SCRIPT = (ROOT / "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001" /
             "verify_projected_response.py")

DEPENDENCIES = {
    "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001/THEOREM.md":
        "44fc28edc9820d2b4ea67cef9f83beef60e53bfe6407b582ffbeccfe42f756c5",
    "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001/verify_finite_tt_four_point.py":
        "fb44d45290c0530098c0e8f9593dff1c0f8149d42598f842374b61004a8ff6c2",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/THEOREM.md":
        "6fc221a31151340b91a946d33e442971c1373500e067c354b6c610e3964edb1c",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/RESULT.md":
        "b5d4c3de99aa4e100519c19a9b74de487b47c1a2d3671204e77740bd9094771a",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/verify_projected_source_rank.py":
        "0e93d84f9eb7cf7fdd62b5a14d5c6705c74841899dd1676bdd7e7a41eb971a00",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/MANIFEST.sha256":
        "651a66b9afd7545b04aa80e5f90952fda9327d011ecd19b973aa80ff51a739f3",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/SEAL.sha256":
        "8301da6bbc026d0e14d985592c5dabe3d91072957c4aaa4b1bebf1f45aadd894",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/THEOREM.md":
        "db3e12d50fd1cb41cddc722a0445cdeaef6a52d49704fa6df1028dfd9abcba1b",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/RESULT.md":
        "81d8f732d8395d757c5405c11c093156bf3a1c2dfae4670cda2db061a5c7e262",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/verify_projected_response.py":
        "87152814b07eeef30794626a003adc8b97f16eeb79252ee99a94d307e822aad8",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/MANIFEST.sha256":
        "17499bb43606bd110657218585bb9b8bebe73358cb49da6276ea8ae7d42d0c47",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/SEAL.sha256":
        "1f829571f82aa16eebac53492470b6a0a883b8401a198ee137a34f9a6b052c9b",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/INDEPENDENT_HOSTILE_AUDIT/AUDIT_MANIFEST.sha256":
        "311296bdf3428c77020495e0f92dc29e71ff94427181b00dc5ee08556e37b686",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/INDEPENDENT_HOSTILE_AUDIT/AUDIT_RESULT.json":
        "b33ba43cec5e716dc5d4ef05156f81e1ef43448d1b6ab8ef01a281dfd89eb9a7",
    "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/INDEPENDENT_HOSTILE_AUDIT/AUDIT_SEAL.sha256":
        "fc018a877a9dbcc65233d32d67e6c97e1fb352f4115a5e5f3c1302479ad49018",
}

MANIFEST_FILES = {
    "DEPENDENCIES.sha256", "README.md", "RESULT.md", "SELF_AUDIT.md",
    "THEOREM.md", "VERIFICATION.txt", "derive_complete_h6_source.py",
}

checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


for relative, expected in DEPENDENCIES.items():
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(),
          f"dependency is a regular file: {relative}")
    check(digest(path) == expected, f"dependency hash is pinned: {relative}")


with redirect_stdout(StringIO()):
    fo = runpy.run_path(str(FO_SCRIPT))

states = fo["states"]
state_index = fo["state_index"]
edges = fo["edges"]
edge_labels = fo["edge_labels"]
incidence = fo["incidence"]
orbits = fo["translation_orbits"]
representatives = tuple(orbit[0] for orbit in orbits)
EDGE_COUNT = len(edges)
VERTEX_COUNT = fo["VERTEX_COUNT"]
CELL_COUNT = fo["CELL_COUNT"]

check(len(states) == 180 and len(orbits) == 6 and
      all(len(orbit) == 30 for orbit in orbits),
      "FO interface supplies six free length-thirty translation orbits")

# The cyclic translation shifts only the cell index; it never rotates a link
# label.  This is load bearing because a lab-frame anisotropic tensor source
# would not be orbit-constant under a label-permuting automorphism.
edge_translation = tuple(4*((edge//4 + 1) % CELL_COUNT) + edge % 4
                         for edge in range(EDGE_COUNT))
vertex_translation = tuple(
    (vertex + 1) % CELL_COUNT if vertex < CELL_COUNT else
    CELL_COUNT + ((vertex-CELL_COUNT+1) % CELL_COUNT)
    for vertex in range(VERTEX_COUNT))
check(all(edge_labels[edge_translation[edge]] == edge_labels[edge]
          for edge in range(EDGE_COUNT)),
      "FO cyclic translation preserves every tetrahedral edge label")
check(all(tuple(sorted(edges[edge_translation[edge]])) ==
          tuple(sorted(vertex_translation[v] for v in edges[edge]))
          for edge in range(EDGE_COUNT)),
      "FO cyclic translation is an incidence-preserving graph automorphism")

SIGNS = ((1, 1, 1), (1, -1, -1),
         (-1, 1, -1), (-1, -1, 1))
PAIRS = tuple(combinations(range(4), 2))


def dyad_integer(vector):
    """Dyad as coefficients against (jxx,jyy,jzz,jxy,jxz,jyz).

    The off-diagonal factor two is therefore intentional: ``j:R`` contains
    ``2*jxy*Rxy`` for symmetric tensors.  These are source covector slots,
    not six orthonormalized tensor components.
    """
    x, y, z = vector
    return (x*x, y*y, z*z, 2*x*y, 2*x*z, 2*y*z)


# Eight times dE_C/dj at lambda=-1/2.  Since dE/dj=(1/4)Rhat,
# and Rhat=dyad(sign_b-sign_a)/8, the scaled row is dyad/4.
# Every tetrahedral difference component is 0 or +/-2, so it is integral.
EPRIME8 = {
    (a, b): tuple(value // 4 for value in
                  dyad_integer(tuple(SIGNS[b][axis]-SIGNS[a][axis]
                                     for axis in range(3))))
    for a, b in PAIRS
}

# D_a=n_a n_a^T has coordinate row dyad_integer(SIGNS[a])/3.
EDGE_DYAD_NUMERATOR = tuple(dyad_integer(vector) for vector in SIGNS)


def add_rows(*rows):
    return tuple(sum(values, F(0)) for values in zip(*rows))


def scale_row(scale, row):
    return tuple(F(scale)*value for value in row)


def poly_mul(first, second):
    """Multiply delta polynomials through degree two."""
    return tuple(sum(first[k]*second[n-k] for k in range(n+1))
                 for n in range(3))


def sequential_path_weight(path_data, numerator_prime):
    """Independent dual-polynomial replay of one path through delta^2."""
    base = (F(1), F(0), F(0))
    derivative = [(F(0), F(0), F(0)) for _ in range(6)]
    for gap, gp8 in path_data:
        inverse = F(1, gap)
        resolvent = (-inverse, -inverse*inverse,
                     -inverse*inverse*inverse)
        previous = base
        base = poly_mul(base, resolvent)
        for component in range(6):
            gp = F(gp8[component], 8)
            resolvent_prime = (gp*inverse*inverse,
                               2*gp*inverse**3,
                               3*gp*inverse**4)
            derivative[component] = add_rows(
                poly_mul(derivative[component], resolvent),
                poly_mul(previous, resolvent_prime))
    derivative = tuple(add_rows(row, scale_row(numerator_prime[c], base))
                       for c, row in enumerate(derivative))
    return base, derivative


def vertex_z(state):
    table = []
    for vertex in range(VERTEX_COUNT):
        row = [None]*4
        for _, edge in incidence[vertex]:
            row[edge_labels[edge]] = 1-2*((state >> edge) & 1)
        assert all(value in (-1, 1) for value in row)
        table.append(tuple(row))
    return tuple(table)


# Exact generator covariance supplies the orbit reduction: translation is a
# bijection on every selected edge cluster, preserves its multiplicities and
# labels, and carries every local occupation row to the translated vertex.
translate_state = fo["translate_state"]
check(all(translate_state(state) in state_index for state in states),
      "one-step translation remains inside the FO component")
check(all(
    vertex_z(translate_state(state))[vertex_translation[v]] ==
    vertex_z(state)[v]
    for state in states for v in range(VERTEX_COUNT)),
    "all local labeled occupations are exactly translation covariant")


def local_eprime8(z):
    total = [0]*6
    for a, b in PAIRS:
        coefficient = z[a]*z[b]
        for component, value in enumerate(EPRIME8[(a, b)]):
            total[component] += coefficient*value
    return tuple(total)


def cluster_signature(selected, state, z_initial):
    """All path-relevant gaps and source derivatives for one edge cluster."""
    selected = tuple(selected)
    labels = tuple(edge_labels[edge] for edge in selected)
    subset_data = []
    for mask in range(1, 1 << len(selected)):
        degree_delta = {}
        affected = set()
        toggled = []
        for position, edge in enumerate(selected):
            if not (mask >> position) & 1:
                continue
            toggled.append(edge)
            occupation_delta = -1 if (state >> edge) & 1 else 1
            for vertex in edges[edge]:
                affected.add(vertex)
                degree_delta[vertex] = (degree_delta.get(vertex, 0) +
                                        occupation_delta)
        gap = sum(value*value for value in degree_delta.values())
        gap_prime8 = [0]*6
        for vertex in affected:
            before = z_initial[vertex]
            after = list(before)
            for edge in toggled:
                if vertex in edges[edge]:
                    after[edge_labels[edge]] *= -1
            difference = tuple(new-old for new, old in
                               zip(local_eprime8(tuple(after)),
                                   local_eprime8(before)))
            for component, value in enumerate(difference):
                gap_prime8[component] += value
        assert gap > 0
        subset_data.append((gap, tuple(gap_prime8)))
    return labels, tuple(subset_data)


ORDER_CACHE = {}


def multiset_orders(multiplicities):
    key = tuple(multiplicities)
    if key not in ORDER_CACHE:
        labels = tuple(position for position, multiplicity in
                       enumerate(multiplicities) for _ in range(multiplicity))
        ORDER_CACHE[key] = tuple(sorted(set(permutations(labels))))
    return ORDER_CACHE[key]


def order_is_irreducible(order):
    parity = 0
    for position in order[:-1]:
        parity ^= 1 << position
        if parity == 0:
            return False
    return True


def even_partitions(total, largest=None):
    """Integer partitions of total with positive even parts."""
    if total == 0:
        return ((),)
    if largest is None:
        largest = total
    result = []
    for part in range(min(total, largest), 1, -2):
        if part % 2:
            continue
        for tail in even_partitions(total-part, part):
            result.append((part,) + tail)
    return tuple(result)


check(set(even_partitions(2)) == {(2,)},
      "length-two closed words have the sole even partition (2)")
check(set(even_partitions(4)) == {(4,), (2, 2)},
      "length-four closed words exhaust partitions (4) and (2,2)")
check(set(even_partitions(6)) == {(6,), (4, 2), (2, 2, 2)},
      "length-six closed words exhaust partitions (6), (4,2), and (2,2,2)")
check(len(multiset_orders((2,))) == 1 and
      len(multiset_orders((2, 2))) == 6 and
      len(multiset_orders((4, 2))) == 15 and
      len(multiset_orders((2, 4))) == 15 and
      len(multiset_orders((2, 2, 2))) == 90,
      "unique multiset-order enumerators have exact multinomial sizes")
check(sum(order_is_irreducible(order) for order in multiset_orders((4,))) == 0 and
      sum(order_is_irreducible(order) for order in multiset_orders((6,))) == 0,
      "one-edge (4) and (6) families are wholly reducible at prefix two")


WEIGHT_CACHE = {}
SEQUENTIAL_WITNESSES = set()


def cluster_weight(signature, multiplicities):
    """Return k(delta) and its six source derivatives through delta^2."""
    cache_key = (signature, tuple(multiplicities))
    if cache_key in WEIGHT_CACHE:
        return WEIGHT_CACHE[cache_key]
    labels, subset_data = signature
    base_total = (F(0), F(0), F(0))
    derivative_total = [(F(0), F(0), F(0)) for _ in range(6)]

    # First derivative of the complete numerator product at j=0:
    # -(1/2) sum_occurrence D_label = -sum m_a dyad(sign_a)/6.
    numerator_prime = [F(0)]*6
    for position, multiplicity in enumerate(multiplicities):
        label = labels[position]
        for component, value in enumerate(EDGE_DYAD_NUMERATOR[label]):
            numerator_prime[component] -= F(multiplicity*value, 6)

    for order in multiset_orders(multiplicities):
        parity = 0
        path_data = []
        irreducible = True
        for position in order[:-1]:
            parity ^= 1 << position
            if parity == 0:
                irreducible = False
                break
            path_data.append(subset_data[parity-1])
        if not irreducible:
            continue

        # Closed formulas for the product of resolvents through delta^2.
        # This is algebraically identical to sequential dual-polynomial
        # multiplication but avoids six repeated polynomial convolutions.
        inverses = [F(1, gap) for gap, _ in path_data]
        p0 = F(-1)
        for inverse in inverses:
            p0 *= inverse
        s1 = sum(inverses)
        s2 = sum(inverse*inverse for inverse in inverses)
        base = (p0, p0*s1, p0*(s1*s1+s2)/2)
        derivative = []
        for component in range(6):
            g0 = g1 = g2 = F(0)
            for (gap, gp8), inverse in zip(path_data, inverses):
                gp = F(gp8[component], 8)
                g0 += gp*inverse
                g1 += gp*inverse*inverse
                g2 += gp*inverse*inverse*inverse
            derivative.append((
                numerator_prime[component]*base[0] - base[0]*g0,
                numerator_prime[component]*base[1] -
                    (base[1]*g0 + base[0]*g1),
                numerator_prime[component]*base[2] -
                    (base[2]*g0 + base[1]*g1 + base[0]*g2),
            ))

        # One exact, independent sequential product per multiplicity family
        # guards the closed coefficient formulas without doubling the full
        # enumeration cost.  Equality is Fraction equality, not tolerance.
        family = tuple(multiplicities)
        if family not in SEQUENTIAL_WITNESSES:
            sequential = sequential_path_weight(path_data, numerator_prime)
            check(sequential == (base, tuple(derivative)),
                  f"closed resolvent derivative matches sequential exact replay for {family}")
            SEQUENTIAL_WITNESSES.add(family)
        base_total = add_rows(base_total, base)
        for component in range(6):
            derivative_total[component] = add_rows(
                derivative_total[component], derivative[component])

    result = base_total, tuple(derivative_total)
    WEIGHT_CACHE[cache_key] = result
    return result


def add_weight(accumulator, weight, multiplicity=1):
    base_acc, derivative_acc = accumulator
    base, derivative = weight
    base_acc = add_rows(base_acc, scale_row(multiplicity, base))
    derivative_acc = tuple(add_rows(row, scale_row(multiplicity, new))
                           for row, new in zip(derivative_acc, derivative))
    return base_acc, derivative_acc


ZERO_PACKET = ((F(0), F(0), F(0)),
               tuple((F(0), F(0), F(0)) for _ in range(6)))


def signature_counts(state, size):
    z_initial = vertex_z(state)
    counts = Counter()
    for selected in combinations(range(EDGE_COUNT), size):
        counts[cluster_signature(selected, state, z_initial)] += 1
    return counts


def irreducible_kernels(state):
    z_initial = vertex_z(state)
    k2 = ZERO_PACKET
    for edge in range(EDGE_COUNT):
        signature = cluster_signature((edge,), state, z_initial)
        k2 = add_weight(k2, cluster_weight(signature, (2,)))

    pair_counts = signature_counts(state, 2)
    k4 = ZERO_PACKET
    k6_pair = ZERO_PACKET
    for signature, count in pair_counts.items():
        k4 = add_weight(k4, cluster_weight(signature, (2, 2)), count)
        combined = (
            add_rows(cluster_weight(signature, (4, 2))[0],
                     cluster_weight(signature, (2, 4))[0]),
            tuple(add_rows(first, second) for first, second in zip(
                cluster_weight(signature, (4, 2))[1],
                cluster_weight(signature, (2, 4))[1])))
        k6_pair = add_weight(k6_pair, combined, count)

    triple_counts = signature_counts(state, 3)
    print(f"    raw signature types: pairs={len(pair_counts)} triples={len(triple_counts)}",
          flush=True)
    k6 = k6_pair
    for signature, count in triple_counts.items():
        k6 = add_weight(k6, cluster_weight(signature, (2, 2, 2)), count)
    return k2, k4, k6, len(pair_counts), len(triple_counts)


def folded_coefficients(k2, k4, k6):
    b2, d2 = k2
    b4, d4 = k4
    b6, d6 = k6
    a2 = b2[0]
    da2 = tuple(row[0] for row in d2)
    a4 = b4[0] + a2*b2[1]
    da4 = tuple(d4[c][0] + da2[c]*b2[1] + a2*d2[c][1]
                for c in range(6))
    a6 = (b6[0] + a2*b4[1] + a4*b2[1] +
          a2*a2*b2[2])
    da6 = tuple(
        d6[c][0] + da2[c]*b4[1] + a2*d4[c][1] +
        da4[c]*b2[1] + a4*d2[c][1] +
        2*a2*da2[c]*b2[2] + a2*a2*d2[c][2]
        for c in range(6))
    return (a2, a4, a6), (da2, da4, da6)


def fmt(value):
    return str(value.numerator) if value.denominator == 1 else str(value)


all_results = []
all_type_counts = []
for orbit_index, state in enumerate(representatives):
    print(f"ORBIT {orbit_index}: enumerating diagonal words", flush=True)
    k2, k4, k6, pair_types, triple_types = irreducible_kernels(state)
    coefficients, derivatives = folded_coefficients(k2, k4, k6)
    all_results.append((coefficients, derivatives))
    all_type_counts.append((pair_types, triple_types))
    print(f"ORBIT {orbit_index}: pair_types={pair_types} triple_types={triple_types}")
    print("  scalar", *(fmt(value) for value in coefficients))
    for order, row in zip((2, 4, 6), derivatives):
        print(f"  d{order}", *(fmt(value) for value in row))


reference_scalar = all_results[0][0]
check(all(result[0] == reference_scalar for result in all_results),
      "source-off folded diagonal remains scalar on all six orbits")

EXPECTED_SCALAR = (F(-60), F(-35), F(-893, 9))
EXPECTED_DERIVATIVES = (
    (
        (F(46), F(22), F(22), F(0), F(0), F(0)),
        (F(46), F(22), F(22), F(0), F(0), F(0)),
        (F(46), F(21), F(23), F(0), F(0), F(0)),
        (F(46), F(23), F(21), F(0), F(0), F(0)),
        (F(46), F(22), F(22), F(0), F(0), F(0)),
        (F(46), F(22), F(22), F(0), F(0), F(0)),
    ),
    (
        (F(541, 6), F(97, 6), F(97, 6), F(0), F(0), F(0)),
        (F(541, 6), F(97, 6), F(97, 6), F(0), F(0), F(0)),
        (F(541, 6), F(157, 12), F(77, 4), F(0), F(0), F(0)),
        (F(541, 6), F(77, 4), F(157, 12), F(0), F(0), F(0)),
        (F(541, 6), F(97, 6), F(97, 6), F(0), F(0), F(0)),
        (F(541, 6), F(97, 6), F(97, 6), F(0), F(0), F(0)),
    ),
    (
        (F(635503, 1350), F(50611, 1350), F(50611, 1350), F(0), F(0), F(0)),
        (F(635503, 1350), F(50611, 1350), F(50611, 1350), F(0), F(0), F(0)),
        (F(635503, 1350), F(52481, 2700), F(149963, 2700), F(0), F(0), F(0)),
        (F(635503, 1350), F(149963, 2700), F(52481, 2700), F(0), F(0), F(0)),
        (F(635503, 1350), F(50611, 1350), F(50611, 1350), F(0), F(0), F(0)),
        (F(635503, 1350), F(50611, 1350), F(50611, 1350), F(0), F(0), F(0)),
    ),
)
EXPECTED_TYPE_COUNTS = ((325, 5444), (362, 6021), (275, 4793),
                        (277, 4823), (360, 5968), (316, 5220))

check(reference_scalar == EXPECTED_SCALAR,
      "source-off BW/Feshbach coefficients are exactly (-60,-35,-893/9)")
observed_derivatives = tuple(
    tuple(all_results[orbit][1][order] for orbit in range(6))
    for order in range(3))
check(observed_derivatives == EXPECTED_DERIVATIVES,
      "all eighteen orbit/order derivative rows match the frozen exact table")
check(tuple(all_type_counts) == EXPECTED_TYPE_COUNTS,
      "pair/triple signature multiplicities match exact orbit checksums")
check(len(WEIGHT_CACHE) == 13725 and
      set(ORDER_CACHE) == {(2,), (4,), (6,), (2, 2), (4, 2),
                           (2, 4), (2, 2, 2)} and
      SEQUENTIAL_WITNESSES == {(2,), (2, 2), (4, 2), (2, 4), (2, 2, 2)},
      "all 13725 weights, five kernels, and two reducible families were replayed")

print("SOURCE_OFF_SCALAR", *(fmt(value) for value in reference_scalar))
for order_index, order in enumerate((2, 4, 6)):
    print(f"DERIVATIVE_ORDER_{order}")
    for orbit_index, result in enumerate(all_results):
        print(orbit_index, *(fmt(value) for value in result[1][order_index]))
print(f"CACHE signatures={len(WEIGHT_CACHE)} order_families={len(ORDER_CACHE)}")


# ---------------------------------------------------------------------------
# Exact component reduction and physical scaling.

def direct_pair_row(state):
    """FV-PURE ideal-Coulomb Q_pair/Ud in source covector coordinates."""
    total = [F(0)]*6
    z_table = vertex_z(state)
    for z in z_table:
        for a, b in PAIRS:
            root = tuple(SIGNS[b][axis]-SIGNS[a][axis]
                         for axis in range(3))
            # lambda=-1/2 and Rhat=dyad(root)/8.
            contribution = scale_row(F(-z[a]*z[b], 16),
                                     dyad_integer(root))
            total = list(add_rows(total, contribution))
    return tuple(total)


direct_rows = tuple(direct_pair_row(state) for state in representatives)
EXPECTED_DIRECT_ROWS = (
    (F(52), F(4), F(4), F(0), F(0), F(0)),
    (F(52), F(4), F(4), F(0), F(0), F(0)),
    (F(52), F(2), F(6), F(0), F(0), F(0)),
    (F(52), F(6), F(2), F(0), F(0), F(0)),
    (F(52), F(4), F(4), F(0), F(0), F(0)),
    (F(52), F(4), F(4), F(0), F(0), F(0)),
)
check(direct_rows == EXPECTED_DIRECT_ROWS,
      "exact direct pair source rows reproduce all six FO orbits")

# Q=-2 dH/dj.  The exact derivatives above were computed at Ud=1, so
# Qdiag^(m)=Ud*x^m*qdiag_unit^(m), x=h/Ud, with qdiag_unit=-2*da_m.
qdiag_unit_rows = tuple(
    tuple(scale_row(-2, row) for row in EXPECTED_DERIVATIVES[order])
    for order in range(3))
relative_pair_coefficients = (F(-1), F(-37, 12), F(-16247, 900))
identity_coefficients = (F(-40), F(-20), F(-374, 135))
for order_index, order in enumerate((2, 4, 6)):
    coefficient = relative_pair_coefficients[order_index]
    identity_coefficient = identity_coefficients[order_index]
    expected_identity = (identity_coefficient,)*3 + (F(0),)*3
    check(all(
        add_rows(qdiag_unit_rows[order_index][orbit],
                 scale_row(-coefficient, direct_rows[orbit])) ==
        expected_identity
        for orbit in range(6)),
        f"Qdiag^{order} equals {coefficient} Qpair plus an exact Hilbert identity")


def f_e(x):
    """Through-H6 diagonal renormalization polynomial for active E."""
    x = F(x)
    return (F(1) - x*x - F(37, 12)*x**4 -
            F(16247, 900)*x**6)


FORMAL_F_E = (F(1), F(0), F(-1), F(0), F(-37, 12),
              F(0), F(-16247, 900))
check(FORMAL_F_E[0] == 1,
      "f_E has unit constant term and cannot cancel E as a formal series")

# A finite polynomial evaluation has one positive root.  This is only the
# root of the through-H6 truncation, not a threshold: no convergence radius
# is established and H8+ can move or remove it.
# Put y=x^2.  The exact derivative coefficients of
# p(y)=1-y-(37/12)y^2-(16247/900)y^3 are all strictly negative for y>=0.
# Together with p(0)=1 and the negative leading term, this proves existence
# and uniqueness of the positive root without numerical root finding.
P_Y = (F(1), F(-1), F(-37, 12), F(-16247, 900))
P_Y_DERIVATIVE = (F(-1), F(-37, 6), F(-16247, 300))
check(P_Y[0] == 1 and P_Y[-1] < 0 and
      all(coefficient < 0 for coefficient in P_Y_DERIVATIVE),
      "exact y-polynomial is strictly decreasing from positive to negative")


def p_y(y):
    y = F(y)
    return sum(coefficient*y**power
               for power, coefficient in enumerate(P_Y))


check(p_y(F(1, 4)) == F(15853, 57600) > 0 and
      p_y(F(729, 2500)) == F(-2157513587, 1562500000000) < 0,
      "exact signs bracket the unique root at 1/2 < x < 27/50")
y_roots = np.roots(tuple(float(value) for value in reversed(P_Y)))
positive_y = [root.real for root in y_roots
              if abs(root.imag) < 1e-12 and root.real > 0]
check(len(positive_y) == 1,
      "the finite through-H6 f_E polynomial has one positive-x root")
truncated_root = float(np.sqrt(positive_y[0]))
check(abs(truncated_root - 0.5398271903) < 5e-10,
      "finite through-H6 cancellation root has its numerical checksum")


# ---------------------------------------------------------------------------
# Compose the complete diagonal source with the audited irreducible H6 ring.

with redirect_stdout(StringIO()):
    fw = runpy.run_path(str(FW_SCRIPT))

H = fw["H"]
ground = fw["ground"]
ground_energy = fw["ground_energy"]
eigenvalues = fw["eigenvalues"]
eigenvectors = fw["eigenvectors"]
Q_direct = fw["Q_direct"]
Q_ring = fw["Q_ring"]
BASIS_COORDINATES = fw["BASIS_COORDINATES"]
translation_permutation = fw["translation_permutation"]
U0 = fw["U0"]
dimension = len(states)

check(np.max(np.abs(U0)) == 1/np.sqrt(30) and
      np.linalg.norm(U0.T@U0-np.eye(6)) < 2e-15,
      "equal-size orbit isometry uses the normalized 1/sqrt(30) basis")

orbit_of_state = {}
for orbit_index, orbit in enumerate(orbits):
    for state in orbit:
        orbit_of_state[state] = orbit_index
check(len(orbit_of_state) == dimension,
      "each FO state belongs to exactly one source orbit")

qdiag_unit = []
for order_index in range(3):
    source = np.zeros((6, dimension, dimension), dtype=float)
    for row, state in enumerate(states):
        source[:, row, row] = np.array(
            qdiag_unit_rows[order_index][orbit_of_state[state]], dtype=float)
    qdiag_unit.append(source)
qdiag_unit = tuple(qdiag_unit)

check(all(np.array_equal(source,
                         source[np.ix_(range(6), translation_permutation,
                                       translation_permutation)])
          for source in qdiag_unit),
      "all generated diagonal tensor sources are exactly translation invariant")
check(all(np.linalg.norm(
    np.array([U0.T@component@U0 for component in source]) -
    np.array([np.diag([float(qdiag_unit_rows[order][orbit][component])
                       for orbit in range(6)])
              for component in range(6)])) < 2e-12
          for order, source in enumerate(qdiag_unit)),
      "normalized zero-momentum blocks reproduce the six exact orbit values")
check(np.linalg.norm(np.array([
    Q_direct[:, state_index[state], state_index[state]]
    for state in representatives]) - np.array(direct_rows, dtype=float)) < 2e-12,
      "FW direct source and independently reconstructed exact rows agree")

identity = np.eye(dimension)
isotropic_coordinate_identity = np.zeros((6, dimension, dimension))
for component in range(3):
    isotropic_coordinate_identity[component] = identity
for order, (source, coefficient, scalar) in enumerate(zip(
        qdiag_unit, relative_pair_coefficients, identity_coefficients),
        start=1):
    exact_reduction = (float(coefficient)*Q_direct +
                       float(scalar)*isotropic_coordinate_identity)
    check(np.linalg.norm(source-exact_reduction) < 2e-12,
          f"full 180-state Qdiag order index {order} obeys exact pair-plus-identity reduction")


def matrix_rank_psd(matrix, relative_tolerance=2e-10):
    eigen = np.linalg.eigvalsh((matrix+matrix.T.conj())/2)
    tolerance = relative_tolerance*max(1.0, float(np.max(np.abs(eigen))))
    return int(np.count_nonzero(eigen > tolerance))


def complete_source(x):
    """Complete homogeneous coordinate source in J6 units through H6."""
    x = F(x)
    rho = F(8, 63)/x**6       # Ud/J6, J6=(63/8) Ud x^6
    diagonal = Q_direct.copy()
    for order, source in zip((2, 4, 6), qdiag_unit):
        diagonal += float(x**order)*source
    coordinate_source = float(rho)*diagonal + Q_ring

    rho_e = rho*f_e(x)
    identity_scalar = rho*sum(identity_coefficients[index]*x**order
                              for index, order in enumerate((2, 4, 6)))
    reduced = (float(rho_e)*Q_direct + Q_ring +
               float(identity_scalar)*isotropic_coordinate_identity)
    check(np.linalg.norm(coordinate_source-reduced) < 2e-9,
          f"complete source at x={x} reduces to rho_E Qpair + Qring plus identity")
    return coordinate_source, rho, rho_e, identity_scalar


def response_ranks(coordinate_source, rho_e, label):
    """Independent 180-state commutator and Lehmann replay."""
    Q = np.einsum("ac,cij->aij", BASIS_COORDINATES, coordinate_source)
    centered = np.array([operator-np.trace(operator)/dimension*identity
                         for operator in Q])
    operator_gram = np.einsum("aij,bij->ab", centered, centered)
    commutators = np.array([H@operator-operator@H for operator in Q])
    commutator_gram = np.einsum("aij,bij->ab", commutators, commutators)

    excited_values = eigenvalues[1:]
    excited_vectors = eigenvectors[:, 1:]
    gaps = excited_values-ground_energy
    amplitudes = np.stack([excited_vectors.T@operator@ground
                           for operator in Q], axis=1)
    delta_1 = 2+2*np.sqrt(2.0)
    delta_2 = 4+2*np.sqrt(2.0)
    mask_1 = np.abs(gaps-delta_1) < 2e-10
    mask_2 = np.abs(gaps-delta_2) < 2e-10
    residue_1 = amplitudes[mask_1].T@amplitudes[mask_1]
    residue_2 = amplitudes[mask_2].T@amplitudes[mask_2]
    residue_all = amplitudes.T@amplitudes
    vector_1 = np.array((0, float(rho_e)/np.sqrt(2),
                         -float(rho_e)*np.sqrt(3/2),
                         -3/np.sqrt(2), -3/np.sqrt(2), 0))
    vector_2 = np.array((0, 0, 0, 3/np.sqrt(2), -3/np.sqrt(2), 0))
    exact_1 = np.outer(vector_1, vector_1)
    exact_2 = np.outer(vector_2, vector_2)
    residue_scale = max(1.0, np.linalg.norm(exact_1)+np.linalg.norm(exact_2))
    check(np.linalg.norm(residue_1-exact_1) < 2e-11*residue_scale,
          f"{label}: first exact pole residue uses rho_E")
    check(np.linalg.norm(residue_2-exact_2) < 2e-11*residue_scale,
          f"{label}: second exact pole residue is unchanged")
    check(np.linalg.norm(residue_all-residue_1-residue_2) <
          2e-11*residue_scale,
          f"{label}: complete source has no additional ground-state pole")

    moment_1 = np.empty((6, 6), dtype=float)
    for a in range(6):
        ad = commutators[a]
        for b in range(6):
            moment_1[a, b] = ground@(ad@Q[b]-Q[b]@ad)@ground
    exact_moment_1 = -2*(delta_1*exact_1+delta_2*exact_2)
    check(np.linalg.norm(moment_1-exact_moment_1) < 3e-10*residue_scale,
          f"{label}: M1 is the exact gap-weighted two-pole moment")

    return (matrix_rank_psd(operator_gram),
            matrix_rank_psd(commutator_gram),
            matrix_rank_psd(residue_all),
            matrix_rank_psd(-moment_1))


for sample_x in (F(2, 5), F(1, 2)):
    source, rho, rho_e, identity_scalar = complete_source(sample_x)
    ranks = response_ranks(source, rho_e, f"x={sample_x}")
    check(ranks == (5, 3, 2, 2),
          f"x={sample_x}: complete finite hierarchy remains 5 -> 3 -> 2 -> 2")
    print(f"SCALE x={sample_x} rho={rho} fE={f_e(sample_x)} rhoE={rho_e} "
          f"identity={identity_scalar}")

# At the exact algebraic root of the *truncated polynomial*, all diagonal
# nonidentity response is absent and only the ring source remains.  This is a
# useful algebraic rank stratum, not a physical critical-point prediction.
root_ranks = response_ranks(Q_ring, F(0), "formal finite-H6 root")
check(root_ranks == (4, 2, 2, 2),
      "finite-H6 cancellation stratum has ranks 4 -> 2 -> 2 -> 2")


# Documentary scope becomes an executable gate once the packet exists.
document_paths = tuple(LANE / name for name in
                       ("THEOREM.md", "RESULT.md", "SELF_AUDIT.md"))
if all(path.is_file() for path in document_paths):
    joined = " ".join(" ".join(path.read_text().split())
                      for path in document_paths)
    for phrase in (
        "FV-PURE", "Q_diag^(2,4,6)", "x=h/U_d", "f_E(x)",
        "formal power-series unit", "finite through-H6 cancellation",
        "not a physical threshold", "homogeneous", "Ward", "RGRL-B",
        "thermodynamic", "gravity", "Newton",
    ):
        check(phrase in joined, f"scope text retains: {phrase}")
    for forbidden in (
        "the root is a gravity threshold",
        "the rank-two response proves a graviton",
        "this finite component proves a Ward identity",
        "Newton's constant is calculated",
    ):
        check(forbidden not in joined, f"forbidden promotion absent: {forbidden}")


manifest = LANE / "MANIFEST.sha256"
if manifest.is_file():
    listed = set()
    for line in manifest.read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split("  ", 1)
        listed.add(relative)
        check(digest(LANE / relative) == expected,
              f"manifest custody {relative}")
    check(listed == MANIFEST_FILES,
          "manifest lists exactly the seven frozen builder files")


print(f"SUMMARY {checks}/{checks} complete-H6 source-response checks passed")
print("EXACT source-off a2=-60 a4=-35 a6=-893/9")
print("EXACT fE=1-x^2-(37/12)x^4-(16247/900)x^6")
print("GENERIC_RANKS operator_mod_identity=5 adH=3 retarded=2 M1=2")
print("TRUNCATED_ROOT x=0.5398271903 ranks=4,2,2,2; NOT_A_THRESHOLD")
print("CEILING finite homogeneous fixed-through-H6 FV-PURE; H8+/Ward/RGRLB/gravity/G open")
