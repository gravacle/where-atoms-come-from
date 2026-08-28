#!/usr/bin/env python3
"""Native-support, nonzero-momentum complete-H6 source response on FO.

The default mode runs only low-cost structural and order-two tests.  The
expensive order-four/order-six diagonal-word and ring-history enumeration is
enabled explicitly with ``--full``.

No source-off interaction is added.  The only new object is a ledger which
retains the native location of each already-frozen source insertion before the
Feshbach reduction:

* Coulomb-pair and virtual-gap insertions live on A/B vertices;
* hopping-numerator insertions live on the physical link midpoint;
* differentiated folds inherit the location of the differentiated factor.

The uniform sum of the ledger must reproduce FX exactly before any m=1 result
is scored.  Graph momentum is an exact translation label; it is not a claim of
emergent continuum locality, a Ward identity, or gravity.
"""

from argparse import ArgumentParser
from contextlib import redirect_stdout
from fractions import Fraction as F
from io import StringIO
from itertools import combinations, permutations
import json
from pathlib import Path
import runpy

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FO_SCRIPT = (ROOT / "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001" /
             "verify_finite_tt_four_point.py")
FW_SCRIPT = (ROOT / "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001" /
             "verify_projected_response.py")

parser = ArgumentParser()
parser.add_argument("--full", action="store_true",
                    help="run the expensive H4/H6 and ring enumeration")
args = parser.parse_args()

checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


with redirect_stdout(StringIO()):
    fo = runpy.run_path(str(FO_SCRIPT))

states = fo["states"]
state_index = fo["state_index"]
edges = fo["edges"]
edge_labels = fo["edge_labels"]
incidence = fo["incidence"]
orbits = fo["translation_orbits"]
representatives = tuple(orbit[0] for orbit in orbits)
translate_state = fo["translate_state"]
hexagons = fo["hexagons"]
ring_patterns = fo["ring_patterns"]
CELL_COUNT = fo["CELL_COUNT"]
VERTEX_COUNT = fo["VERTEX_COUNT"]
EDGE_COUNT = len(edges)
dimension = len(states)

check((CELL_COUNT, VERTEX_COUNT, EDGE_COUNT, len(hexagons), dimension) ==
      (30, 60, 120, 120, 180),
      "FY uses exactly the FO Z30 quotient and 180-state component")
check(len(orbits) == 6 and all(len(orbit) == 30 for orbit in orbits),
      "FY uses the six free FO translation orbits")


# ---------------------------------------------------------------------------
# Native support and exact ledger arithmetic.

SUPPORT_A = 0
SUPPORT_B = 1
SUPPORT_EDGE0 = 2
SUPPORT_NAMES = ("A", "B", "e0", "e1", "e2", "e3")
SUPPORT_COUNT = len(SUPPORT_NAMES)
ZERO_ROW = (F(0),) * 6

SIGNS = ((1, 1, 1), (1, -1, -1),
         (-1, 1, -1), (-1, -1, 1))
PAIRS = tuple(combinations(range(4), 2))


def dyad_integer(vector):
    """Covector coefficients against (jxx,jyy,jzz,jxy,jxz,jyz)."""
    x, y, z = vector
    return (x*x, y*y, z*z, 2*x*y, 2*x*z, 2*y*z)


def add_rows(*rows):
    return tuple(sum(values, F(0)) for values in zip(*rows))


def scale_row(scale, row):
    return tuple(F(scale) * value for value in row)


def ledger_add_scaled(target, source, scale=F(1)):
    scale = F(scale)
    if not scale:
        return target
    for key, row in source.items():
        updated = add_rows(target.get(key, ZERO_ROW), scale_row(scale, row))
        if any(updated):
            target[key] = updated
        elif key in target:
            del target[key]
    return target


def ledger_scaled(source, scale):
    result = {}
    ledger_add_scaled(result, source, scale)
    return result


def ledger_linear(*terms):
    """Return an exact linear combination of sparse native ledgers."""
    result = {}
    for scale, source in terms:
        ledger_add_scaled(result, source, scale)
    return result


def ledger_sum_row(source):
    return add_rows(*(source.values())) if source else ZERO_ROW


def freeze_ledger(source):
    return tuple((support, cell, row)
                 for (support, cell), row in sorted(source.items()))


def thaw_ledger(source):
    return {(support, cell): row for support, cell, row in source}


def shift_ledger(source, amount):
    result = {}
    for (support, cell), row in source.items():
        result[(support, (cell + amount) % CELL_COUNT)] = row
    return result


def vertex_support(vertex):
    if vertex < CELL_COUNT:
        return SUPPORT_A, vertex
    return SUPPORT_B, vertex - CELL_COUNT


def edge_support(edge):
    return SUPPORT_EDGE0 + edge_labels[edge], edge // 4


check(all(0 <= support < SUPPORT_COUNT and 0 <= cell < CELL_COUNT
          for edge in range(EDGE_COUNT)
          for support, cell in (edge_support(edge),)),
      "every hopping insertion has one native link species and cell")
check(all(vertex_support(vertex_translation)[1] ==
          (vertex_support(vertex)[1] + 1) % CELL_COUNT and
          vertex_support(vertex_translation)[0] == vertex_support(vertex)[0]
          for vertex in range(VERTEX_COUNT)
          for vertex_translation in
          ((vertex + 1) % CELL_COUNT if vertex < CELL_COUNT else
           CELL_COUNT + (vertex - CELL_COUNT + 1) % CELL_COUNT,)),
      "A/B source supports translate without species rotation")


# FO embeds the cyclic quotient into the diamond primitive coordinates.  The
# phases below are used only after the exact support ledger has been formed.
BASIS_OFFSET = np.array((-0.25, -0.25, -0.25), dtype=float)


def fractional_wavevector(momentum):
    q = momentum * np.array((1, 5, 19), dtype=float) / CELL_COUNT
    return q - np.round(q)


def basis_phase(support, momentum):
    if momentum == 0 or support == SUPPORT_A:
        return 1.0 + 0.0j
    q = fractional_wavevector(momentum)
    if support == SUPPORT_B:
        return np.exp(2j * np.pi * np.dot(q, BASIS_OFFSET))
    label = support - SUPPORT_EDGE0
    displacement = BASIS_OFFSET.copy()
    if label:
        displacement[label - 1] += 1
    # The link source is carried at the same midpoint used by FO's one-link
    # nonzero-momentum observable.
    return np.exp(1j * np.pi * np.dot(q, displacement))


def evaluate_ledger(source, momentum):
    result = np.zeros(6, dtype=complex)
    species_phases = tuple(basis_phase(support, momentum)
                           for support in range(SUPPORT_COUNT))
    for (support, cell), row in source.items():
        phase = (np.exp(2j * np.pi * momentum * cell / CELL_COUNT) *
                 species_phases[support])
        result += phase * np.array(row, dtype=float)
    return result


check(all(abs(basis_phase(support, 0) - 1) < 1e-15
          for support in range(SUPPORT_COUNT)),
      "all native supports coalesce under the uniform source")
check(abs(basis_phase(SUPPORT_B, 1) - np.exp(2j*np.pi/24)) < 2e-15,
      "B-sublattice m=1 phase matches the frozen FO embedding")
SUPPORT_PHASE_EXPONENT_240 = (0, 10, 5, 9, 25, 201)
check(all(abs(basis_phase(support, 1) -
              np.exp(2j*np.pi*exponent/240)) < 3e-15
          for support, exponent in enumerate(SUPPORT_PHASE_EXPONENT_240)),
      "all A/B/link-midpoint m=1 phases match exact zeta_240 exponents")
check(all(abs(basis_phase(support, CELL_COUNT-1) -
              basis_phase(support, 1).conjugate()) < 5e-14
          for support in range(SUPPORT_COUNT)),
      "every native support obeys phase(-m)=conjugate phase(m)")


# All m=1 phases lie in Q(zeta_240).  The following small exact polynomial
# layer is used only after the H6 ledgers exist.  It decides whether a source
# really lifts the homogeneous coefficient; floating residuals do not decide
# that theorem.
def polynomial_trim(polynomial):
    result = list(polynomial)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def polynomial_exact_divide(numerator, denominator):
    numerator = polynomial_trim(tuple(F(value) for value in numerator))
    denominator = polynomial_trim(tuple(F(value) for value in denominator))
    quotient = [F(0)] * max(1, len(numerator)-len(denominator)+1)
    while len(numerator) >= len(denominator) and any(numerator):
        shift = len(numerator)-len(denominator)
        factor = numerator[-1]/denominator[-1]
        quotient[shift] += factor
        for index, value in enumerate(denominator):
            numerator[index+shift] -= factor*value
        numerator = polynomial_trim(numerator)
    if any(numerator):
        raise AssertionError("exact polynomial division left a remainder")
    return tuple(polynomial_trim(quotient))


def cyclotomic_polynomial(n):
    divisors = tuple(value for value in range(1, n+1) if n % value == 0)
    table = {}
    for divisor in divisors:
        polynomial = (F(-1),) + (F(0),)*(divisor-1) + (F(1),)
        for proper in (value for value in divisors
                       if value < divisor and divisor % value == 0):
            polynomial = polynomial_exact_divide(polynomial, table[proper])
        table[divisor] = polynomial
    return table[n]


PHI_240 = cyclotomic_polynomial(240)
check(len(PHI_240)-1 == 64 and PHI_240[-1] == 1,
      "Phi_240 is exact monic of degree phi(240)=64")


def polynomial_reduce_phi240(polynomial):
    remainder = polynomial_trim(polynomial)
    while len(remainder) >= len(PHI_240):
        shift = len(remainder)-len(PHI_240)
        factor = remainder[-1]
        for index, value in enumerate(PHI_240):
            remainder[index+shift] -= factor*value
        remainder = polynomial_trim(remainder)
    return tuple(remainder + [F(0)]*(64-len(remainder)))


def polynomial_multiply(first, second):
    result = [F(0)] * (len(first)+len(second)-1)
    for i, left in enumerate(first):
        if not left:
            continue
        for j, right in enumerate(second):
            if right:
                result[i+j] += left*right
    return result


def polynomial_subtract(first, second):
    size = max(len(first), len(second))
    return [(first[index] if index < len(first) else F(0)) -
            (second[index] if index < len(second) else F(0))
            for index in range(size)]


def ledger_component_polynomial(ledger, component):
    polynomial = [F(0)]*240
    for (support, cell), row in ledger.items():
        exponent = (8*cell + SUPPORT_PHASE_EXPONENT_240[support]) % 240
        polynomial[exponent] += row[component]
    return tuple(polynomial)


def exact_m1_relation(target_ledgers, reference_ledgers, coefficient):
    """Fraction/cyclotomic proof of target(m=1)=coefficient*reference(m=1)."""
    for target, reference in zip(target_ledgers, reference_ledgers):
        residual = ledger_linear((1, target), (-F(coefficient), reference))
        for component in range(6):
            if any(polynomial_reduce_phi240(
                    ledger_component_polynomial(residual, component))):
                return False
    return True


def exact_m1_cross_witness(target_ledgers, reference_ledgers):
    """Return a nonzero exact 2x2 cross-minor, or None."""
    records = []
    for orbit in range(6):
        target_numeric = evaluate_ledger(target_ledgers[orbit], 1)
        reference_numeric = evaluate_ledger(reference_ledgers[orbit], 1)
        for component in range(6):
            records.append((orbit, component, target_numeric[component],
                            reference_numeric[component]))
    candidates = []
    for first_index, first in enumerate(records):
        for second_index in range(first_index+1, len(records)):
            second = records[second_index]
            cross_numeric = (first[2]*second[3] -
                             second[2]*first[3])
            if abs(cross_numeric) > 1e-9:
                candidates.append((first_index, second_index))
    # Numeric ordering only locates likely witnesses.  Exact reduction is the
    # proof; if needed, the exhaustive tail covers an ill-conditioned minor.
    all_pairs = ((first, second) for first in range(len(records))
                 for second in range(first+1, len(records)))
    seen = set(candidates)
    for pair in tuple(candidates) + tuple(pair for pair in all_pairs
                                          if pair not in seen):
        first, second = (records[index] for index in pair)
        target_first = ledger_component_polynomial(
            target_ledgers[first[0]], first[1])
        target_second = ledger_component_polynomial(
            target_ledgers[second[0]], second[1])
        reference_first = ledger_component_polynomial(
            reference_ledgers[first[0]], first[1])
        reference_second = ledger_component_polynomial(
            reference_ledgers[second[0]], second[1])
        cross = polynomial_subtract(
            polynomial_multiply(target_first, reference_second),
            polynomial_multiply(target_second, reference_first))
        remainder = polynomial_reduce_phi240(cross)
        if any(remainder):
            return (first[0], first[1], second[0], second[1],
                    tuple((power, value) for power, value in
                          enumerate(remainder) if value))
    return None


# ---------------------------------------------------------------------------
# Pair energy and gap derivatives, retained at their native vertices.

EPRIME8 = {
    (a, b): tuple(value // 4 for value in
                  dyad_integer(tuple(SIGNS[b][axis] - SIGNS[a][axis]
                                     for axis in range(3))))
    for a, b in PAIRS
}
EDGE_DYAD_NUMERATOR = tuple(dyad_integer(vector) for vector in SIGNS)


def vertex_z(state):
    table = []
    for vertex in range(VERTEX_COUNT):
        row = [None] * 4
        for _, edge in incidence[vertex]:
            row[edge_labels[edge]] = 1 - 2 * ((state >> edge) & 1)
        assert all(value in (-1, 1) for value in row)
        table.append(tuple(row))
    return tuple(table)


def local_eprime8(z):
    total = [0] * 6
    for a, b in PAIRS:
        coefficient = z[a] * z[b]
        for component, value in enumerate(EPRIME8[(a, b)]):
            total[component] += coefficient * value
    return tuple(total)


def direct_pair_ledger(state):
    """FV-PURE Coulomb Q_pair/Ud, before summing native vertices."""
    result = {}
    for vertex, z in enumerate(vertex_z(state)):
        row = [F(0)] * 6
        for a, b in PAIRS:
            root = tuple(SIGNS[b][axis] - SIGNS[a][axis]
                         for axis in range(3))
            row = list(add_rows(row, scale_row(F(-z[a]*z[b], 16),
                                               dyad_integer(root))))
        result[vertex_support(vertex)] = tuple(row)
    return result


EXPECTED_DIRECT_ROWS = (
    (F(52), F(4), F(4), F(0), F(0), F(0)),
    (F(52), F(4), F(4), F(0), F(0), F(0)),
    (F(52), F(2), F(6), F(0), F(0), F(0)),
    (F(52), F(6), F(2), F(0), F(0), F(0)),
    (F(52), F(4), F(4), F(0), F(0), F(0)),
    (F(52), F(4), F(4), F(0), F(0), F(0)),
)

direct_ledgers = tuple(direct_pair_ledger(state) for state in representatives)
check(tuple(ledger_sum_row(ledger) for ledger in direct_ledgers) ==
      EXPECTED_DIRECT_ROWS,
      "native vertex sum exactly recovers all six FX direct-pair rows")
check(direct_pair_ledger(translate_state(representatives[0])) ==
      shift_ledger(direct_ledgers[0], 1),
      "direct-pair native ledger is exactly translation covariant")


def cluster_signature(selected, state, z_initial):
    """Exact gaps plus sparse vertex-resolved gap-source insertions."""
    selected = tuple(selected)
    supports = tuple(edge_support(edge) for edge in selected)
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
        gap = sum(value * value for value in degree_delta.values())
        gap_prime8 = {}
        for vertex in affected:
            before = z_initial[vertex]
            after = list(before)
            for edge in toggled:
                if vertex in edges[edge]:
                    after[edge_labels[edge]] *= -1
            difference = tuple(new - old for new, old in
                               zip(local_eprime8(tuple(after)),
                                   local_eprime8(before)))
            if any(difference):
                gap_prime8[vertex_support(vertex)] = difference
        assert gap > 0
        subset_data.append((gap, freeze_ledger(gap_prime8)))
    return supports, labels, tuple(subset_data)


ORDER_CACHE = {}


def multiset_orders(multiplicities):
    key = tuple(multiplicities)
    if key not in ORDER_CACHE:
        labels = tuple(position for position, multiplicity in
                       enumerate(multiplicities)
                       for _ in range(multiplicity))
        ORDER_CACHE[key] = tuple(sorted(set(permutations(labels))))
    return ORDER_CACHE[key]


COEFFICIENT_CACHE = {}


def path_coefficients(gaps, multiplicities, max_degree=2):
    """Source-independent word coefficients for one exact gap signature.

    The returned gap coefficient multiplies the stored ``8*dDelta/dj`` row.
    Caching this object keeps the expensive multiset-order enumeration separate
    from the native source locations.
    """
    cache_key = (tuple(gaps), tuple(multiplicities), int(max_degree))
    if cache_key in COEFFICIENT_CACHE:
        return COEFFICIENT_CACHE[cache_key]
    base_total = [F(0), F(0), F(0)]
    gap_coefficients = [[F(0) for _ in gaps] for _ in range(3)]

    for order in multiset_orders(multiplicities):
        parity = 0
        path = []
        irreducible = True
        for position in order[:-1]:
            parity ^= 1 << position
            if parity == 0:
                irreducible = False
                break
            path.append((parity - 1, F(1, gaps[parity - 1])))
        if not irreducible:
            continue

        inverses = [inverse for _, inverse in path]
        p0 = F(-1)
        for inverse in inverses:
            p0 *= inverse
        s1 = sum(inverses, F(0))
        s2 = sum((inverse * inverse for inverse in inverses), F(0))
        base = (p0, p0*s1, p0*(s1*s1+s2)/2)
        for degree in range(max_degree + 1):
            base_total[degree] += base[degree]

        for subset_index, inverse in path:
            for degree in range(max_degree + 1):
                # -B(delta) * [dDelta/(Delta-delta)].  The native
                # ledger stores 8*dDelta, hence the factor 1/8.
                coefficient = F(0)
                for gap_degree in range(degree + 1):
                    coefficient -= (base[degree-gap_degree] *
                                    inverse**(gap_degree + 1) / 8)
                gap_coefficients[degree][subset_index] += coefficient

    result = (tuple(base_total),
              tuple(tuple(row) for row in gap_coefficients))
    COEFFICIENT_CACHE[cache_key] = result
    return result


def cluster_weight(signature, multiplicities, max_degree=2):
    """Return k(delta) and tagged d_j k(delta), both through delta^2."""
    supports, labels, subset_data = signature
    gaps = tuple(gap for gap, _ in subset_data)
    base_total, gap_coefficients = path_coefficients(
        gaps, multiplicities, max_degree)
    derivative_total = [{}, {}, {}]

    numerator_prime = {}
    for position, multiplicity in enumerate(multiplicities):
        row = scale_row(F(-multiplicity, 6),
                        EDGE_DYAD_NUMERATOR[labels[position]])
        ledger_add_scaled(numerator_prime, {supports[position]: row})

    for degree in range(max_degree + 1):
        ledger_add_scaled(derivative_total[degree], numerator_prime,
                          base_total[degree])
        for subset_index, (_, frozen_gp8) in enumerate(subset_data):
            coefficient = gap_coefficients[degree][subset_index]
            if coefficient:
                ledger_add_scaled(derivative_total[degree],
                                  thaw_ledger(frozen_gp8), coefficient)

    return base_total, tuple(derivative_total)


def add_weight(accumulator, weight, multiplicity=1):
    base_acc, derivative_acc = accumulator
    base, derivative = weight
    base_acc = tuple(old + F(multiplicity)*new
                     for old, new in zip(base_acc, base))
    for degree in range(3):
        ledger_add_scaled(derivative_acc[degree], derivative[degree],
                          multiplicity)
    return base_acc, derivative_acc


def order_two_kernel(state):
    z_initial = vertex_z(state)
    accumulator = ((F(0), F(0), F(0)), [{}, {}, {}])
    for edge in range(EDGE_COUNT):
        signature = cluster_signature((edge,), state, z_initial)
        accumulator = add_weight(accumulator,
                                 cluster_weight(signature, (2,)))
    return accumulator


EXPECTED_DA2 = (
    (F(46), F(22), F(22), F(0), F(0), F(0)),
    (F(46), F(22), F(22), F(0), F(0), F(0)),
    (F(46), F(21), F(23), F(0), F(0), F(0)),
    (F(46), F(23), F(21), F(0), F(0), F(0)),
    (F(46), F(22), F(22), F(0), F(0), F(0)),
    (F(46), F(22), F(22), F(0), F(0), F(0)),
)

order_two = tuple(order_two_kernel(state) for state in representatives)
check(all(result[0][0] == F(-60) for result in order_two),
      "native order-two source-off sum is the FX scalar a2=-60")
check(tuple(ledger_sum_row(result[1][0]) for result in order_two) ==
      EXPECTED_DA2,
      "native order-two insertion sum exactly recovers every FX da2 row")

translated_order_two = order_two_kernel(translate_state(representatives[0]))
check(translated_order_two[0] == order_two[0][0] and
      translated_order_two[1][0] == shift_ledger(order_two[0][1][0], 1),
      "tagged order-two folded derivative is exactly translation covariant")

qdiag2_ledgers = tuple(ledger_scaled(result[1][0], -2)
                       for result in order_two)
check(tuple(ledger_sum_row(ledger) for ledger in qdiag2_ledgers) ==
      tuple(scale_row(-2, row) for row in EXPECTED_DA2),
      "Q=-2 dH/dj is applied before the order-two Fourier lift")

for orbit, (direct, q2) in enumerate(zip(direct_ledgers, qdiag2_ledgers)):
    direct_m1 = evaluate_ledger(direct, 1) / np.sqrt(VERTEX_COUNT)
    q2_m1 = evaluate_ledger(q2, 1) / np.sqrt(VERTEX_COUNT)
    check(np.all(np.isfinite(direct_m1)) and np.all(np.isfinite(q2_m1)),
          f"orbit {orbit} has finite native-support m=1 pair and Qdiag2 rows")


# ---------------------------------------------------------------------------
# Expensive complete-H6 stages.  These functions are defined in structural
# mode but are called only by ``--full``.  That separation lets the ownership,
# normalization, and H2 recovery gates run before any long enumeration.

EXPECTED_DA4 = (
    (F(541, 6), F(97, 6), F(97, 6), F(0), F(0), F(0)),
    (F(541, 6), F(97, 6), F(97, 6), F(0), F(0), F(0)),
    (F(541, 6), F(157, 12), F(77, 4), F(0), F(0), F(0)),
    (F(541, 6), F(77, 4), F(157, 12), F(0), F(0), F(0)),
    (F(541, 6), F(97, 6), F(97, 6), F(0), F(0), F(0)),
    (F(541, 6), F(97, 6), F(97, 6), F(0), F(0), F(0)),
)
EXPECTED_DA6 = (
    (F(635503, 1350), F(50611, 1350), F(50611, 1350), F(0), F(0), F(0)),
    (F(635503, 1350), F(50611, 1350), F(50611, 1350), F(0), F(0), F(0)),
    (F(635503, 1350), F(52481, 2700), F(149963, 2700), F(0), F(0), F(0)),
    (F(635503, 1350), F(149963, 2700), F(52481, 2700), F(0), F(0), F(0)),
    (F(635503, 1350), F(50611, 1350), F(50611, 1350), F(0), F(0), F(0)),
    (F(635503, 1350), F(50611, 1350), F(50611, 1350), F(0), F(0), F(0)),
)


def zero_packet():
    return (F(0), F(0), F(0)), [{}, {}, {}]


def state_subset_factory(state):
    """Build exact subset data with state-local vertex tables.

    A triple signature contains three singleton and three pair subsets.  Those
    objects are reused across all triples rather than recomputed 118 times.
    The local pair-energy derivative has only sixteen label-toggle masks per
    q4 vertex, so it is also tabulated once.  This changes no summand or
    ownership assignment.
    """
    z_initial = vertex_z(state)
    occupation_delta = tuple(-1 if (state >> edge) & 1 else 1
                             for edge in range(EDGE_COUNT))
    local_delta = []
    for vertex, before in enumerate(z_initial):
        before_prime = local_eprime8(before)
        rows = []
        for label_mask in range(16):
            after = tuple(-value if (label_mask >> label) & 1 else value
                          for label, value in enumerate(before))
            rows.append(tuple(new-old for new, old in
                              zip(local_eprime8(after), before_prime)))
        local_delta.append(tuple(rows))

    def datum(selected):
        degree_delta = {}
        label_masks = {}
        for edge in selected:
            label_bit = 1 << edge_labels[edge]
            for vertex in edges[edge]:
                degree_delta[vertex] = (degree_delta.get(vertex, 0) +
                                        occupation_delta[edge])
                label_masks[vertex] = label_masks.get(vertex, 0) ^ label_bit
        gap = sum(value*value for value in degree_delta.values())
        if gap <= 0:
            raise AssertionError("diagonal proper-prefix gap left Q2")
        gap_prime8 = {}
        for vertex, label_mask in label_masks.items():
            row = local_delta[vertex][label_mask]
            if any(row):
                gap_prime8[vertex_support(vertex)] = row
        return F(gap), freeze_ledger(gap_prime8)

    return datum


def signature_from_subset_data(selected, subset_data):
    return (tuple(edge_support(edge) for edge in selected),
            tuple(edge_labels[edge] for edge in selected),
            tuple(subset_data))


def complete_irreducible_kernels(state, orbit_index):
    """All diagonal irreducible closed words through H6 for one FO orbit."""
    datum = state_subset_factory(state)
    singleton_data = tuple(datum((edge,)) for edge in range(EDGE_COUNT))
    pair_data = {}
    k2 = order_two[orbit_index]
    k4 = zero_packet()
    k6 = zero_packet()

    for selected in combinations(range(EDGE_COUNT), 2):
        pair_data[selected] = datum(selected)
        first, second = selected
        signature = signature_from_subset_data(
            selected, (singleton_data[first], singleton_data[second],
                       pair_data[selected]))
        k4 = add_weight(k4, cluster_weight(signature, (2, 2), 1))
        k6 = add_weight(k6, cluster_weight(signature, (4, 2), 0))
        k6 = add_weight(k6, cluster_weight(signature, (2, 4), 0))

    print(f"ORBIT {orbit_index}: pair clusters complete", flush=True)
    for count, selected in enumerate(combinations(range(EDGE_COUNT), 3), 1):
        first, second, third = selected
        signature = signature_from_subset_data(selected, (
            singleton_data[first], singleton_data[second],
            pair_data[(first, second)], singleton_data[third],
            pair_data[(first, third)], pair_data[(second, third)],
            datum(selected),
        ))
        k6 = add_weight(k6, cluster_weight(signature, (2, 2, 2), 0))
        if count % 50000 == 0:
            print(f"ORBIT {orbit_index}: triple clusters {count}/280840",
                  flush=True)
    return k2, k4, k6


def folded_coefficients_native(k2, k4, k6):
    """BW/Feshbach folds with the differentiated factor retaining support."""
    b2, d2 = k2
    b4, d4 = k4
    b6, d6 = k6
    a2 = b2[0]
    da2 = dict(d2[0])
    a4 = b4[0] + a2*b2[1]
    da4 = ledger_linear((1, d4[0]), (b2[1], da2), (a2, d2[1]))
    a6 = (b6[0] + a2*b4[1] + a4*b2[1] + a2*a2*b2[2])
    da6 = ledger_linear(
        (1, d6[0]),
        (b4[1], da2), (a2, d4[1]),
        (b2[1], da4), (a4, d2[1]),
        (2*a2*b2[2], da2), (a2*a2, d2[2]),
    )
    return (a2, a4, a6), (da2, da4, da6)


def enumerate_diagonal_ledgers():
    results = []
    for orbit_index, state in enumerate(representatives):
        print(f"ORBIT {orbit_index}: native diagonal H6 enumeration",
              flush=True)
        results.append(folded_coefficients_native(
            *complete_irreducible_kernels(state, orbit_index)))

    expected_derivatives = (EXPECTED_DA2, EXPECTED_DA4, EXPECTED_DA6)
    check(all(result[0] == (F(-60), F(-35), F(-893, 9))
              for result in results),
          "all native source-off folds recover FX a2/a4/a6")
    for order_index, order in enumerate((2, 4, 6)):
        check(tuple(ledger_sum_row(results[orbit][1][order_index])
                    for orbit in range(6)) ==
              expected_derivatives[order_index],
              f"native support sum exactly recovers all FX da{order} rows")
    qdiag = tuple(
        tuple(ledger_scaled(results[orbit][1][order_index], -2)
              for orbit in range(6))
        for order_index in range(3)
    )
    check(qdiag[0] == qdiag2_ledgers,
          "full enumeration reproduces the independently gated H2 ledgers")
    return results, qdiag


def ring_cluster_signature(cycle, state):
    """Proper-prefix gaps for one oriented flippable six-link history."""
    cycle = tuple(cycle)
    z_initial = vertex_z(state)
    supports = tuple(edge_support(edge) for edge in cycle)
    labels = tuple(edge_labels[edge] for edge in cycle)
    subset_data = []
    for mask in range(1, 1 << 6):
        # The full mask is the other ice endpoint and never appears in a
        # five-resolvent path.  A positive dummy keeps the parity-index table
        # rectangular; an executable gate below proves it is never sampled.
        if mask == (1 << 6) - 1:
            subset_data.append((F(1), ()))
            continue
        degree_delta = {}
        affected = set()
        toggled = []
        for position, edge in enumerate(cycle):
            if not (mask >> position) & 1:
                continue
            toggled.append(edge)
            occupation_delta = -1 if (state >> edge) & 1 else 1
            for vertex in edges[edge]:
                affected.add(vertex)
                degree_delta[vertex] = (degree_delta.get(vertex, 0) +
                                        occupation_delta)
        gap = sum(value*value for value in degree_delta.values())
        if gap <= 0:
            raise AssertionError("ring proper-prefix gap left Q2")
        gap_prime8 = {}
        for vertex in affected:
            before = z_initial[vertex]
            after = list(before)
            for edge in toggled:
                if vertex in edges[edge]:
                    after[edge_labels[edge]] *= -1
            difference = tuple(new-old for new, old in
                               zip(local_eprime8(tuple(after)),
                                   local_eprime8(before)))
            if any(difference):
                gap_prime8[vertex_support(vertex)] = difference
        subset_data.append((F(gap), freeze_ledger(gap_prime8)))
    return supports, labels, tuple(subset_data)


def oriented_ring_derivative(state, cycle):
    base, derivative = cluster_weight(
        ring_cluster_signature(cycle, state), (1, 1, 1, 1, 1, 1), 0)
    if base[0] != F(-63, 8):
        raise AssertionError("ring source-off coefficient is not -63/8")
    return derivative[0]


def translate_edge(edge, amount=1):
    return 4*((edge//4 + amount) % CELL_COUNT) + edge % 4


def transition_inventory():
    inventory = {}
    for state in states:
        row = state_index[state]
        for (mask, first, second), cycle in zip(ring_patterns, hexagons):
            if (state & mask) not in (first, second):
                continue
            column = state_index[state ^ mask]
            key = tuple(sorted((row, column)))
            previous = inventory.setdefault(key, tuple(cycle))
            if set(previous) != set(cycle):
                raise AssertionError("one FO transition received two hexagons")
    return inventory


def translated_transition(key, amount=1):
    return tuple(sorted(state_index[translate_state(states[index], amount)]
                        for index in key))


def enumerate_ring_ledgers():
    inventory = transition_inventory()
    check(len(inventory) == 420,
          "native ring inventory has the exact 420 undirected FO transitions")
    unassigned = set(inventory)
    transition_orbits = []
    while unassigned:
        representative = min(unassigned)
        orbit = tuple(translated_transition(representative, amount)
                      for amount in range(CELL_COUNT))
        check(len(set(orbit)) == CELL_COUNT,
              "each native ring-entry orbit has free length thirty")
        transition_orbits.append(orbit)
        unassigned.difference_update(orbit)
    check(len(transition_orbits) == 14,
          "the 420 ring entries reduce to fourteen translation orbits")

    ledgers = {}
    for orbit_index, orbit in enumerate(transition_orbits):
        key = orbit[0]
        cycle = inventory[key]
        first_state, second_state = (states[index] for index in key)
        forward = oriented_ring_derivative(first_state, cycle)
        reverse = oriented_ring_derivative(second_state, cycle)
        # q=-2*dH/dj after the standard 1/2 Hermitian endpoint average.
        q_h6 = ledger_linear((-1, forward), (-1, reverse))
        q_j6 = ledger_scaled(q_h6, F(8, 63))
        for amount, translated_key in enumerate(orbit):
            ledgers[translated_key] = shift_ledger(q_j6, amount)
        print(f"RING orbit {orbit_index}: native 720+720 histories complete",
              flush=True)

    check(set(ledgers) == set(inventory),
          "native ring ledgers cover every FO transition exactly once")
    return inventory, transition_orbits, ledgers


def diagonal_source_from_ledgers(orbit_ledgers, momentum):
    source = np.zeros((6, dimension, dimension), dtype=complex)
    normalization = np.sqrt(VERTEX_COUNT)
    for orbit_index, orbit in enumerate(orbits):
        representative_ledger = orbit_ledgers[orbit_index]
        for amount, state in enumerate(orbit):
            shifted = shift_ledger(representative_ledger, amount)
            value = (np.array(ledger_sum_row(shifted), dtype=complex)
                     if momentum == 0 else evaluate_ledger(shifted, momentum))
            source[:, state_index[state], state_index[state]] = (
                value / normalization)
    return source


def ring_source_from_ledgers(transition_ledgers, momentum):
    source = np.zeros((6, dimension, dimension), dtype=complex)
    normalization = np.sqrt(VERTEX_COUNT)
    for (row, column), ledger in transition_ledgers.items():
        value = (np.array(ledger_sum_row(ledger), dtype=complex)
                 if momentum == 0 else evaluate_ledger(ledger, momentum))
        value = value / normalization
        source[:, row, column] = value
        source[:, column, row] = value
    return source


def translation_covariant(source, momentum, tolerance=3e-11):
    permutation = fo["translation_permutation"]
    shifted = source[:, permutation][:, :, permutation]
    phase = np.exp(2j*np.pi*momentum/CELL_COUNT)
    return np.max(np.abs(shifted-phase*source)) < tolerance


def matrix_rank_gram(vectors, relative_tolerance=3e-10):
    """Numerical rank of a positive Gram matrix, returned with eigenvalues."""
    gram = np.array([[np.vdot(first, second) for second in vectors]
                     for first in vectors], dtype=complex)
    gram = (gram + gram.conj().T) / 2
    eigenvalues_gram = np.linalg.eigvalsh(gram)
    tolerance = relative_tolerance * max(1.0, float(np.max(eigenvalues_gram)))
    return int(np.count_nonzero(eigenvalues_gram > tolerance)), gram


def spectral_packet(source, H, ground, ground_energy, eigenvalues,
                    eigenvectors):
    """Basis-invariant positive-frequency residues for six source channels."""
    vectors = tuple(operator @ ground for operator in source)
    amplitudes = np.stack([eigenvectors.conj().T @ vector
                           for vector in vectors], axis=1)
    gaps = eigenvalues-ground_energy
    weight = np.sum(np.abs(amplitudes)**2, axis=1)
    groups = []
    used = np.zeros(len(gaps), dtype=bool)
    for index in np.argsort(gaps):
        if used[index] or weight[index] < 1e-14:
            continue
        mask = np.abs(gaps-gaps[index]) < 2e-10
        used |= mask
        residue = amplitudes[mask].conj().T @ amplitudes[mask]
        groups.append((float(gaps[index]), residue))
    if any(gap <= 2e-10 for gap, _ in groups):
        raise AssertionError("positive-frequency packet contains a zero gap")
    static = sum((residue/gap for gap, residue in groups),
                 np.zeros((6, 6), dtype=complex))
    moment1 = sum((gap*residue for gap, residue in groups),
                  np.zeros((6, 6), dtype=complex))
    return vectors, groups, static, moment1


def source_ranks(source, H, ground, ground_energy, eigenvalues, eigenvectors):
    operator_rank, _ = matrix_rank_gram(tuple(source))
    commutators = tuple(H@operator-operator@H for operator in source)
    commutator_rank, _ = matrix_rank_gram(commutators)
    vectors, groups, static, moment1 = spectral_packet(
        source, H, ground, ground_energy, eigenvalues, eigenvectors)
    spectral_rank, _ = matrix_rank_gram(vectors)
    static_rank = int(np.linalg.matrix_rank(static, tol=3e-10*max(
        1.0, np.linalg.norm(static))))
    moment_rank = int(np.linalg.matrix_rank(moment1, tol=3e-10*max(
        1.0, np.linalg.norm(moment1))))
    return ((operator_rank, commutator_rank, spectral_rank,
             static_rank, moment_rank), groups, vectors, moment1)


def tensor_contract(source, tensor):
    coefficients = np.array((tensor[0, 0], tensor[1, 1], tensor[2, 2],
                             tensor[0, 1], tensor[0, 2], tensor[1, 2]),
                            dtype=complex)
    return np.einsum("a,aij->ij", coefficients, source)


def longitudinal_sources(source, wavevector):
    kx, ky, kz = wavevector
    return np.array((
        kx*source[0] + ky*source[3]/2 + kz*source[4]/2,
        kx*source[3]/2 + ky*source[1] + kz*source[5]/2,
        kx*source[4]/2 + ky*source[5]/2 + kz*source[2],
    ))


# Exact equivalence witness for the optimized full-enumeration kernel.
_factory_witness_state = representatives[0]
_factory_witness_selected = (0, 17, 83)
_factory_witness_datum = state_subset_factory(_factory_witness_state)
_factory_witness_signature = signature_from_subset_data(
    _factory_witness_selected, (
        _factory_witness_datum((_factory_witness_selected[0],)),
        _factory_witness_datum((_factory_witness_selected[1],)),
        _factory_witness_datum(_factory_witness_selected[:2]),
        _factory_witness_datum((_factory_witness_selected[2],)),
        _factory_witness_datum((_factory_witness_selected[0],
                                _factory_witness_selected[2])),
        _factory_witness_datum((_factory_witness_selected[1],
                                _factory_witness_selected[2])),
        _factory_witness_datum(_factory_witness_selected),
    ))
check(_factory_witness_signature == cluster_signature(
          _factory_witness_selected, _factory_witness_state,
          vertex_z(_factory_witness_state)),
      "optimized singleton/pair reuse is exactly the original triple signature")


# One exact, inexpensive ring witness belongs to the pre-enumeration harness.
# It checks that keeping the native locations has not changed the already-
# frozen homogeneous Hermitian answer.
_witness_inventory = transition_inventory()
_witness_key = min(_witness_inventory)
_witness_cycle = _witness_inventory[_witness_key]
_witness_forward = oriented_ring_derivative(
    states[_witness_key[0]], _witness_cycle)
_witness_reverse = oriented_ring_derivative(
    states[_witness_key[1]], _witness_cycle)
_witness_q = ledger_scaled(
    ledger_linear((-1, _witness_forward), (-1, _witness_reverse)), F(8, 63))
_witness_missing = next(iter(set(range(4)) -
                             {edge_labels[edge] for edge in _witness_cycle}))
_witness_expected = add_rows(
    (F(-31, 6), F(-31, 6), F(-31, 6), F(0), F(0), F(0)),
    scale_row(F(3, 2), dyad_integer(SIGNS[_witness_missing])))
check(ledger_sum_row(_witness_q) == _witness_expected,
      "one native 720+720 ring ledger recovers the exact FW Hermitian row")


def run_full():
    """Enumerate the complete fixed-H6 native source and its m=1 response."""
    diagonal_results, qdiag_ledgers = enumerate_diagonal_ledgers()
    _, ring_orbits, ring_ledgers = enumerate_ring_ledgers()

    with redirect_stdout(StringIO()):
        fw = runpy.run_path(str(FW_SCRIPT))
    H = np.array(fo["hamiltonian"], dtype=complex)
    ground = np.array(fo["ground"], dtype=complex)
    ground_energy = float(fo["ground_energy"])
    eigenvalues = np.array(fo["eigenvalues"], dtype=float)
    eigenvectors = np.array(fo["eigenvectors"], dtype=complex)
    Q_direct_fw = np.array(fw["Q_direct"], dtype=complex)
    Q_ring_fw = np.array(fw["Q_ring"], dtype=complex)
    basis_coordinates = np.array(fw["BASIS_COORDINATES"], dtype=float)

    direct_sources = {
        momentum: diagonal_source_from_ledgers(direct_ledgers, momentum)
        for momentum in (0, 1, CELL_COUNT-1)
    }
    qdiag_sources = tuple({
        momentum: diagonal_source_from_ledgers(orbit_ledgers, momentum)
        for momentum in (0, 1, CELL_COUNT-1)
    } for orbit_ledgers in qdiag_ledgers)
    ring_sources = {
        momentum: ring_source_from_ledgers(ring_ledgers, momentum)
        for momentum in (0, 1, CELL_COUNT-1)
    }

    normalization = np.sqrt(VERTEX_COUNT)
    check(np.max(np.abs(normalization*direct_sources[0]-Q_direct_fw)) < 3e-12,
          "summing native supports at m=0 exactly recovers FW Qpair")

    orbit_of_state = {}
    for orbit_index, orbit in enumerate(orbits):
        for state in orbit:
            orbit_of_state[state] = orbit_index
    expected_derivatives = (EXPECTED_DA2, EXPECTED_DA4, EXPECTED_DA6)
    for order_index, order in enumerate((2, 4, 6)):
        expected = np.zeros((6, dimension, dimension), dtype=complex)
        for state in states:
            row = state_index[state]
            expected[:, row, row] = np.array(scale_row(
                -2, expected_derivatives[order_index][orbit_of_state[state]]),
                dtype=float)
        check(np.max(np.abs(normalization*qdiag_sources[order_index][0]-
                            expected)) < 3e-12,
              f"summing native supports at m=0 exactly recovers FX Qdiag^{order}")

    check(np.max(np.abs(normalization*ring_sources[0]-Q_ring_fw)) < 3e-12,
          "summing native supports at m=0 exactly recovers all FW ring entries")
    check(all(translation_covariant(source[1], 1)
              for source in ((direct_sources,) + qdiag_sources +
                             (ring_sources,))),
          "every complete native source term has exact m=1 translation character")
    for label, source in (("pair", direct_sources),
                          ("diag2", qdiag_sources[0]),
                          ("diag4", qdiag_sources[1]),
                          ("diag6", qdiag_sources[2]),
                          ("ring", ring_sources)):
        check(np.max(np.abs(source[CELL_COUNT-1] -
                            source[1].conj().transpose(0, 2, 1))) < 4e-11,
              f"{label} obeys Q(-m)=Q(m)^dagger")

    # Do not assume that the homogeneous f_E reduction survives source
    # resolution.  Decide each coefficient exactly in Q(zeta_240): equality
    # to its homogeneous coefficient is one admissible result; otherwise a
    # nonzero exact cross-minor must prove nonproportionality.
    homogeneous_coefficients = (F(-1), F(-37, 12), F(-16247, 900))
    lift_results = []
    for order, coefficient, ledgers in zip(
            (2, 4, 6), homogeneous_coefficients, qdiag_ledgers):
        lifts = exact_m1_relation(ledgers, direct_ledgers, coefficient)
        witness = None if lifts else exact_m1_cross_witness(
            ledgers, direct_ledgers)
        check(lifts or witness is not None,
              f"m=1 Qdiag^{order} has an exact lift or nonproportional witness")
        if lifts:
            check(np.max(np.abs(qdiag_sources[(order//2)-1][1] -
                                float(coefficient)*direct_sources[1])) < 4e-11,
                  f"m=1 Qdiag^{order} exactly lifts coefficient {coefficient}")
            lift_results.append({"order": order, "outcome": "exact_lift",
                                 "coefficient": str(coefficient)})
        else:
            lift_results.append({
                "order": order, "outcome": "exact_nonproportional",
                "homogeneous_coefficient": str(coefficient),
                "witness": [str(value) for value in witness],
            })

    ring_exact_nonzero = None
    for key, ledger in ring_ledgers.items():
        for component in range(6):
            remainder = polynomial_reduce_phi240(
                ledger_component_polynomial(ledger, component))
            if any(remainder):
                ring_exact_nonzero = (key, component,
                                      tuple((power, value) for power, value in
                                            enumerate(remainder) if value))
                break
        if ring_exact_nonzero is not None:
            break
    check(ring_exact_nonzero is not None,
          "m=1 ring source has an exact nonzero off-diagonal entry")

    # Floating best-fit residuals are retained only as discovery checksums.
    reference = direct_sources[1].reshape(-1)
    reference_norm = np.vdot(reference, reference).real
    residuals = {}
    for order, source in zip((2, 4, 6), qdiag_sources):
        vector = source[1].reshape(-1)
        coefficient = np.vdot(reference, vector)/reference_norm
        residual = np.linalg.norm(vector-coefficient*reference)/np.linalg.norm(vector)
        residuals[f"diag{order}"] = float(residual)
    ring_vector = ring_sources[1].reshape(-1)
    ring_coefficient = np.vdot(reference, ring_vector)/reference_norm
    residuals["ring"] = float(np.linalg.norm(
        ring_vector-ring_coefficient*reference)/np.linalg.norm(ring_vector))
    check(all(np.isfinite(value) for value in residuals.values()),
          "m=1 best-fit residual checksums are finite")

    q, wavevector, _, _, plus, cross = fo["polarization_data"](1)
    check(np.allclose(q, np.array((1, 5, -11))/30),
          "FY uses FO's frozen shortest m=1 reciprocal alias")
    longitudinal_pair = longitudinal_sources(direct_sources[1], wavevector)
    longitudinal_ratio = (np.linalg.norm(longitudinal_pair) /
                          np.linalg.norm(direct_sources[1]))
    check(longitudinal_ratio > 1e-6,
          "native m=1 pair density fails naive spatial transversality")

    sample_results = []
    fixed_gaps = eigenvalues-ground_energy
    for x in (F(2, 5), F(1, 2)):
        rho = F(8, 63)/x**6
        coordinate_source = float(rho) * (
            direct_sources[1] + float(x**2)*qdiag_sources[0][1] +
            float(x**4)*qdiag_sources[1][1] +
            float(x**6)*qdiag_sources[2][1]) + ring_sources[1]
        source = np.einsum("ac,cij->aij", basis_coordinates,
                           coordinate_source)
        ranks, groups, response_vectors, moment1 = source_ranks(
            source, H, ground, ground_energy, eigenvalues, eigenvectors)
        check(groups and all(np.min(np.abs(fixed_gaps-gap)) < 3e-11
                             for gap, _ in groups),
              f"x={x}: every response pole is an unchanged FO Hamiltonian gap")
        direct_moment = np.array([
            np.vdot(first, (H-ground_energy*np.eye(dimension))@second)
            for first in response_vectors for second in response_vectors
        ]).reshape(6, 6)
        check(np.linalg.norm(moment1-direct_moment) <
              4e-9*max(1.0, np.linalg.norm(moment1)),
              f"x={x}: spectral M1 equals the direct gap-weighted moment")

        plus_operator = tensor_contract(coordinate_source, plus)
        cross_operator = tensor_contract(coordinate_source, cross)
        tt_vectors = (plus_operator@ground, cross_operator@ground)
        tt_rank, tt_gram = matrix_rank_gram(tt_vectors)
        tt_weights = []
        for gap, residue in groups:
            coefficients = np.array((
                plus[0, 0], plus[1, 1], plus[2, 2],
                plus[0, 1], plus[0, 2], plus[1, 2],
            ), dtype=complex)
            coefficients_cross = np.array((
                cross[0, 0], cross[1, 1], cross[2, 2],
                cross[0, 1], cross[0, 2], cross[1, 2],
            ), dtype=complex)
            # Convert coordinate covectors to the orthonormal FW channels
            # used by ``residue``.
            channel_plus = np.linalg.solve(basis_coordinates.T, coefficients)
            channel_cross = np.linalg.solve(basis_coordinates.T,
                                             coefficients_cross)
            weight_plus = float(np.real(channel_plus.conj()@residue@channel_plus))
            weight_cross = float(np.real(channel_cross.conj()@residue@channel_cross))
            if weight_plus+weight_cross > 1e-12:
                tt_weights.append((float(gap), weight_plus, weight_cross))

        complete_longitudinal = longitudinal_sources(coordinate_source,
                                                     wavevector)
        longitudinal_complete_ratio = (np.linalg.norm(complete_longitudinal) /
                                       np.linalg.norm(coordinate_source))
        sample = {
            "x": str(x),
            "rho": str(rho),
            "ranks": ranks,
            "pole_gaps": [gap for gap, _ in groups],
            "pole_ranks": [int(np.linalg.matrix_rank(residue,
                                tol=3e-10*max(1.0, np.linalg.norm(residue))))
                           for _, residue in groups],
            "pole_traces": [float(np.trace(residue).real)
                            for _, residue in groups],
            "tt_rank": tt_rank,
            "tt_gram_eigenvalues": [float(value) for value in
                                    np.linalg.eigvalsh(tt_gram)],
            "tt_poles": tt_weights,
            "longitudinal_ratio": float(longitudinal_complete_ratio),
        }
        sample_results.append(sample)
        print("M1_SAMPLE", json.dumps(sample, sort_keys=True), flush=True)

    check(all(sample["pole_gaps"][0] > 1e-8 for sample in sample_results),
          "the two sampled finite m=1 responses have no zero-energy pole")
    check(all(sample["longitudinal_ratio"] > 1e-6
              for sample in sample_results),
          "complete sampled m=1 sources fail naive spatial transversality")

    discovery = {
        "diagonal_source_off": [str(value) for value in diagonal_results[0][0]],
        "ring_translation_orbits": len(ring_orbits),
        "exact_m1_lift_results": lift_results,
        "exact_ring_witness": [str(value) for value in ring_exact_nonzero],
        "m1_best_fit_residuals": residuals,
        "pair_longitudinal_ratio": float(longitudinal_ratio),
        "two_link_threshold": [float(fo["two_link_threshold"][0]),
                               int(fo["two_link_threshold"][1]),
                               int(fo["two_link_threshold"][2])],
        "samples": sample_results,
    }
    print("FY_RESULT_JSON", json.dumps(discovery, sort_keys=True), flush=True)
    print("FY_CEILING finite Z30 graph momentum; no positive Ward/locality/"
          "massless-pole/gravity inference")


if args.full:
    run_full()

print(f"SUMMARY {checks}/{checks} FY native-support structural checks passed")
print("MODE", "FULL" if args.full else "STRUCTURAL_ONLY")
print("CEILING graph-resolved source ledger; no locality/Ward/gravity claim")
