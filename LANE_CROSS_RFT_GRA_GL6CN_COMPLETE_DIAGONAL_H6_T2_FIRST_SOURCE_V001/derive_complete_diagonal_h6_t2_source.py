#!/usr/bin/env python3
"""Exact diagonal order-h6 pure-T2 first-source vertex on the Q4 F3 parent.

The calculation differentiates every Q resolvent in the canonical GL6AO
Kato expression

    K6 = T6 - b X4 + b^2 A3 - d A2,

with the microscopic source inserted before elimination.  Only Python's
standard library and exact fractions are used.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path


HERE = Path(__file__).resolve().parent
PERIOD = 4
STEPS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
PAIR_ORDER = tuple(combinations(range(4), 2))
T_BASIS = (
    (1, 0, 0, 0, 0, -1),
    (0, 1, 0, 0, -1, 0),
    (0, 0, 1, -1, 0, 0),
)
SELECTOR = frozenset(((0, 1), (0, 2), (1, 1),
                      (1, 3), (2, 0), (2, 1)))
CHECKS = []


def check(condition, label):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def ftext(value):
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def canonical_hash(value):
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(raw).hexdigest()


def add(left, right):
    return tuple(left[i] + right[i] for i in range(3))


def sub(left, right):
    return tuple(left[i] - right[i] for i in range(3))


def reduce_cell(cell):
    return tuple(value % PERIOD for value in cell)


def qadd(left, right):
    return reduce_cell(add(left, right))


def qsub(left, right):
    return reduce_cell(sub(left, right))


CELLS = tuple(product(range(PERIOD), repeat=3))
LINKS = tuple((cell, port) for cell in CELLS for port in range(4))
LINK_INDEX = {link: index for index, link in enumerate(LINKS)}
NODES = tuple((kind, cell) for kind in ("C", "P") for cell in CELLS)


def endpoints(link):
    cell, port = link
    return ("P", cell), ("C", qadd(cell, STEPS[port]))


def incident(node):
    kind, cell = node
    if kind == "P":
        return tuple((cell, port) for port in range(4))
    return tuple((qsub(cell, STEPS[port]), port) for port in range(4))


INCIDENT = {node: incident(node) for node in NODES}
ENDS = {link: endpoints(link) for link in LINKS}
ADJACENT_PAIRS = frozenset(
    tuple(sorted((LINK_INDEX[left], LINK_INDEX[right])))
    for node in NODES for left, right in combinations(INCIDENT[node], 2)
)


def selector(a, b):
    return int((a % PERIOD, b % PERIOD) in SELECTOR)


def background_occupation(link):
    cell, port = link
    x, y, z = cell
    a, b = (y - x) % PERIOD, (x + z) % PERIOD
    values = (selector(a - 1, b + 1), selector(a, b),
              1 - selector(a - 1, b + 1), 1 - selector(a, b))
    return values[port]


BASE_OCCUPIED = frozenset(i for i, link in enumerate(LINKS)
                          if background_occupation(link))


def occupation(occupied, edge_index):
    return int(edge_index in occupied)


def is_locked(occupied):
    return all(sum(LINK_INDEX[link] in occupied for link in INCIDENT[node]) == 2
               for node in NODES)


def subset_energy(occupied, edge_indices):
    edge_indices = tuple(edge_indices)
    charges = tuple(1 if edge not in occupied else -1 for edge in edge_indices)
    energy = 2 * len(edge_indices)
    for left, right in combinations(range(len(edge_indices)), 2):
        pair = tuple(sorted((edge_indices[left], edge_indices[right])))
        if pair in ADJACENT_PAIRS:
            energy += 2 * charges[left] * charges[right]
    return F(energy)


def local_pair_word(occupied, node, toggled=frozenset()):
    z = []
    for link in INCIDENT[node]:
        index = LINK_INDEX[link]
        bit = int(index in occupied) ^ int(index in toggled)
        z.append(1 - 2 * bit)
    return tuple(F(z[a] * z[b]) for a, b in PAIR_ORDER)


def local_t_scores(occupied, node, toggled):
    word = local_pair_word(occupied, node, frozenset(toggled))
    return tuple(sum(F(direction[i]) * word[i] for i in range(6))
                 for direction in T_BASIS)


check(len(LINKS) == 256 and len(NODES) == 128,
      "Q4 has 256 links and 128 nodes")
check(all(len(INCIDENT[node]) == 4 and len(set(INCIDENT[node])) == 4
          for node in NODES), "Q4 is simple degree four")
check(len(ADJACENT_PAIRS) == 768,
      "Q4 has 768 owner-once adjacent-link pairs")
check(is_locked(BASE_OCCUPIED) and len(BASE_OCCUPIED) == 128,
      "selected GL6CC background is exactly degree-two locked")
check(all(local_t_scores(BASE_OCCUPIED, node, ()) == (0, 0, 0)
          for node in NODES), "P M_T P vanishes pointwise on the base state")


def unique_permutations(multiset):
    return tuple(sorted(set(permutations(multiset))))


def differentiated_word_kernel(multiset, powers_rows, energy_by_mask):
    """Coefficient of every intermediate source score in a Q-only word sum.

    Each powers row describes the powers of the successive resolvents.  For
    ordinary T6 it is `(1,1,1,1,1)`.  The three X4 terms use `(2,1,1)`,
    `(1,2,1)`, and `(1,1,2)`.  Differentiating
    `prod_r[-1/E(S_r)]^p_r` multiplies it by
    `-sum_r p_r m(S_r)/E(S_r)`.
    """
    answer = defaultdict(F)
    retained = 0
    for word in unique_permutations(multiset):
        masks = []
        parity = 0
        allowed = True
        for label in word[:-1]:
            parity ^= 1 << label
            if not parity or energy_by_mask[parity] <= 0:
                allowed = False
                break
            masks.append(parity)
        if not allowed:
            continue
        retained += 1
        for powers in powers_rows:
            base = F((-1) ** sum(powers))
            for mask, exponent in zip(masks, powers):
                base /= energy_by_mask[mask] ** exponent
            for mask, exponent in zip(masks, powers):
                answer[mask] += base * (-F(exponent) / energy_by_mask[mask])
    return dict(answer), retained


PAIR_T6_KERNEL = {}
PAIR_X4_KERNEL = {}
for pair_energy in (2, 4, 6):
    energy_map = {1: F(2), 2: F(2), 3: F(pair_energy)}
    combined = defaultdict(F)
    retained_rows = []
    for multiset in ((0, 0, 0, 0, 1, 1), (0, 0, 1, 1, 1, 1)):
        row, retained = differentiated_word_kernel(
            multiset, ((1, 1, 1, 1, 1),), energy_map)
        for mask, value in row.items():
            combined[mask] += value
        retained_rows.append(retained)
    PAIR_T6_KERNEL[pair_energy] = dict(combined)
    check(retained_rows == [4, 4], f"repeated-pair T6 retains 4+4 words at p={pair_energy}")
    xrow, xretained = differentiated_word_kernel(
        (0, 0, 1, 1), ((2, 1, 1), (1, 2, 1), (1, 1, 2)), energy_map)
    PAIR_X4_KERNEL[pair_energy] = xrow
    check(xretained == 4, f"X4 retains four words at p={pair_energy}")


TRIPLE_T6_KERNEL = {}
TRIPLE_CLASSES = (
    ((4, 4, 4), 6),
    ((2, 4, 4), 4),
    ((4, 4, 6), 8),
    ((2, 2, 6), 4),
    ((2, 2, 4), 2),
    ((2, 4, 6), 6),
    ((4, 6, 6), 10),
)
def build_triple_kernel(pair_tuple, triple_energy):
    energy_map = {
        1: F(2), 2: F(2), 4: F(2),
        3: F(pair_tuple[0]), 5: F(pair_tuple[1]),
        6: F(pair_tuple[2]), 7: F(triple_energy),
    }
    row, retained = differentiated_word_kernel(
        (0, 0, 1, 1, 2, 2), ((1, 1, 1, 1, 1),), energy_map)
    check(retained == 60, f"triple T6 retains 60 words for {pair_tuple}/{triple_energy}")
    return row


for pair_tuple, triple_energy in TRIPLE_CLASSES:
    TRIPLE_T6_KERNEL[(pair_tuple, triple_energy)] = build_triple_kernel(
        pair_tuple, triple_energy)


def ordered_triple_signature(occupied, triple):
    a, b, c = triple
    p01 = subset_energy(occupied, (a, b))
    p02 = subset_energy(occupied, (a, c))
    p12 = subset_energy(occupied, (b, c))
    total = subset_energy(occupied, triple)
    signature = ((int(p01), int(p02), int(p12)), int(total))
    canonical = (tuple(sorted(signature[0])), signature[1])
    if canonical not in TRIPLE_CLASSES:
        raise AssertionError(f"triple signature outside AO classes: {signature}")
    if signature not in TRIPLE_T6_KERNEL:
        TRIPLE_T6_KERNEL[signature] = build_triple_kernel(*signature)
    return signature


def triple_kernel_for_signature(signature):
    canonical = (tuple(sorted(signature[0])), signature[1])
    check(canonical in TRIPLE_CLASSES,
          f"abstract rooted triple has an AO class: {signature}")
    if signature not in TRIPLE_T6_KERNEL:
        TRIPLE_T6_KERNEL[signature] = build_triple_kernel(*signature)
    return TRIPLE_T6_KERNEL[signature]


def pair_energy_from_bits(left, right):
    return 2 if left != right else 6


def path_energy(first_pair, second_pair):
    equal_joins = int(first_pair == 6) + int(second_pair == 6)
    return 2 + 4 * equal_joins


def root_inclusive_triple_coefficient(signature):
    """Sum coefficients of all parity subsets containing distinguished label 0."""
    kernel = triple_kernel_for_signature(signature)
    return sum((value for mask, value in kernel.items() if mask & 1), F(0))


def rooted_all_triples_coefficient(m_links):
    """Universal coefficient for a source score that tracks one root edge.

    The counts are the edge-rooted refinement of GL6AO's complete matching,
    one-adjacency, star, and path census.  They depend only on degree four,
    girth at least six, and degree-two locking.
    """
    m_links = F(m_links)
    rows = []

    def add_row(label, count, signature):
        rows.append((label, F(count), signature,
                     root_inclusive_triple_coefficient(signature)))

    matching = (m_links*m_links - 21*m_links + 116) / 2
    add_row("matching", matching, ((4, 4, 4), 6))

    add_row("one_adj_opposite_root_adjacent", 4*(m_links-10),
            ((2, 4, 4), 4))
    add_row("one_adj_opposite_root_disjoint", 2*(m_links-10),
            ((4, 4, 2), 4))
    add_row("one_adj_equal_root_adjacent", 2*(m_links-10),
            ((6, 4, 4), 8))
    add_row("one_adj_equal_root_disjoint", m_links-10,
            ((4, 4, 6), 8))

    add_row("star_root_minority", 2, ((2, 2, 6), 4))
    add_row("star_root_majority", 4, ((2, 6, 2), 4))

    add_row("path_root_middle_both_opposite", 4, ((2, 2, 4), 2))
    add_row("path_root_middle_one_equal", 4, ((2, 6, 4), 6))
    add_row("path_root_middle_both_equal", 1, ((6, 6, 4), 10))
    add_row("path_root_end_both_opposite", 8, ((2, 4, 2), 2))
    add_row("path_root_end_first_opposite", 4, ((2, 4, 6), 6))
    add_row("path_root_end_first_equal", 4, ((6, 4, 2), 6))
    add_row("path_root_end_both_equal", 2, ((6, 4, 6), 10))

    check(sum(count for _label, count, _signature, _coefficient in rows) ==
          (m_links-1)*(m_links-2)/2,
          "edge-rooted AO triple census partitions C(M-1,2)")
    total = sum((count * coefficient
                 for _label, count, _signature, coefficient in rows), F(0))
    return total, rows


def universal_singleton_pattern_coefficient(bits, root_port, m_links):
    """Rooted census prediction for the full T6' local singleton pattern."""
    # Repeated-pair histories, rooted at the distinguished local edge.
    repeated = (
        4 * sum(value for mask, value in PAIR_T6_KERNEL[2].items() if mask & 1)
        + 2 * sum(value for mask, value in PAIR_T6_KERNEL[6].items() if mask & 1)
        + (m_links - 7) *
        sum(value for mask, value in PAIR_T6_KERNEL[4].items() if mask & 1)
    )

    triples, _rows = rooted_all_triples_coefficient(m_links)

    # Correct the root-inclusive coefficient whenever another selected edge
    # is also incident at the source node: those parity masks do not have the
    # local singleton pattern {root}.
    correction = F(0)
    other_ports = tuple(port for port in range(4) if port != root_port)
    root_bit = bits[root_port]

    # Pairs containing a second local edge: remove mask {root,other}.
    for other in other_ports:
        p01 = pair_energy_from_bits(root_bit, bits[other])
        correction += PAIR_T6_KERNEL[p01][3]

    # Triples containing exactly one other local edge.  The third edge is
    # either adjacent at the remote end of root/other or disjoint from both.
    # Remove the masks {root,other} and {root,other,third}.
    for other in other_ports:
        p01 = pair_energy_from_bits(root_bit, bits[other])

        def removed(signature):
            kernel = triple_kernel_for_signature(signature)
            return kernel[3] + kernel[7]

        # Third edge adjacent to root at its remote endpoint: two opposite,
        # one equal relative to root.
        for p02, count in ((2, 2), (6, 1)):
            signature = ((p01, p02, 4), path_energy(p01, p02))
            correction += count * removed(signature)
        # Third edge adjacent to the other local edge at its remote endpoint.
        for p12, count in ((2, 2), (6, 1)):
            signature = ((p01, 4, p12), path_energy(p01, p12))
            correction += count * removed(signature)
        # Otherwise it is disjoint from both selected local edges.
        signature = ((p01, 4, 4), 4 if p01 == 2 else 8)
        correction += (m_links - 10) * removed(signature)

    # Triples containing two other local edges: only mask {root} belongs to
    # the singleton pattern, so remove masks 3,5,7 from the rooted total.
    for first, second in combinations(other_ports, 2):
        p01 = pair_energy_from_bits(root_bit, bits[first])
        p02 = pair_energy_from_bits(root_bit, bits[second])
        p12 = pair_energy_from_bits(bits[first], bits[second])
        kernel = triple_kernel_for_signature(((p01, p02, p12), 4))
        correction += kernel[3] + kernel[5] + kernel[7]

    return repeated + triples - correction


def universal_triple_pattern_coefficient(bits, omitted_port):
    ports = tuple(port for port in range(4) if port != omitted_port)
    p01 = pair_energy_from_bits(bits[ports[0]], bits[ports[1]])
    p02 = pair_energy_from_bits(bits[ports[0]], bits[ports[2]])
    p12 = pair_energy_from_bits(bits[ports[1]], bits[ports[2]])
    kernel = triple_kernel_for_signature(((p01, p02, p12), 4))
    return kernel[7]


def section_universal_rooted_census(m_links=256):
    all_rows = []
    for bits in product((0, 1), repeat=4):
        if sum(bits) != 2:
            continue
        singleton = tuple(universal_singleton_pattern_coefficient(
            bits, port, m_links) for port in range(4))
        triple = tuple(universal_triple_pattern_coefficient(bits, port)
                       for port in range(4))
        check(len(set(singleton)) == 1,
              f"all four rooted singleton coefficients agree for locked word {bits}")
        check(len(set(triple)) == 1,
              f"all four rooted triple coefficients agree for locked word {bits}")

        # Complementary local patterns have identical pair words, and the
        # sum over all four one-port defects is T2-dark.
        defect_words = []
        for port in range(4):
            local = list(bits)
            local[port] ^= 1
            z = tuple(1 - 2*bit for bit in local)
            defect_words.append(tuple(F(z[a]*z[b]) for a, b in PAIR_ORDER))
        summed = tuple(sum(row[pair] for row in defect_words)
                       for pair in range(6))
        check(summed == (F(0),) * 6,
              f"four singleton pair words sum to zero for locked word {bits}")
        check(all(sum(F(t[pair]) * summed[pair] for pair in range(6)) == 0
                  for t in T_BASIS),
              f"uniform odd-pattern sum is T2-dark for locked word {bits}")
        all_rows.append({
            "locked_word": list(bits),
            "singleton_coefficients": [ftext(x) for x in singleton],
            "triple_coefficients": [ftext(x) for x in triple],
        })
    check(len(all_rows) == 6, "all six central locked words covered")
    expected_t6_singleton = (F(15, 128)*m_links*m_links
                             + F(3049, 3456)*m_links + F(8653, 4800))
    expected_x4_singleton = -F(5, 16)*m_links - F(487, 432)
    expected_complete_singleton = (F(1, 128)*m_links*m_links
                                   + F(283, 1152)*m_links + F(8653, 4800))
    check(all(F(row["singleton_coefficients"][0]) == expected_t6_singleton
              for row in all_rows), "rooted T6 singleton coefficient polynomial")
    check(all(F(row["triple_coefficients"][0]) == F(49, 576)
              for row in all_rows), "rooted T6 triple coefficient is 49/576")
    complete_from_folds = (expected_t6_singleton
                           + F(m_links, 2)*expected_x4_singleton
                           + F(m_links*m_links, 4)*F(3, 16)
                           + F(7*m_links, 24)*F(-1, 4))
    check(complete_from_folds == expected_complete_singleton,
          "Kato fold signs give complete singleton coefficient polynomial")
    return {
        "locked_word_rows": all_rows,
        "T6_singleton": ("(15/128)M^2+(3049/3456)M+8653/4800"),
        "T6_triple": "49/576",
        "X4_singleton": "-(5/16)M-487/432",
        "A3_singleton": "3/16",
        "A2_singleton": "-1/4",
        "complete_K6_singleton": "(1/128)M^2+(283/1152)M+8653/4800",
        "complete_K6_triple": "49/576",
        "pure_T2_conclusion": "zero pointwise for all six locked words",
    }


@dataclass(frozen=True)
class Dual:
    value: F
    derivative: F = F(0)

    def __add__(self, other):
        other = as_dual(other)
        return Dual(self.value + other.value,
                    self.derivative + other.derivative)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.value, -self.derivative)

    def __sub__(self, other):
        return self + (-as_dual(other))

    def __rsub__(self, other):
        return as_dual(other) - self

    def __mul__(self, other):
        other = as_dual(other)
        return Dual(self.value * other.value,
                    self.derivative * other.value + self.value * other.derivative)

    __rmul__ = __mul__

    def reciprocal(self):
        return Dual(1 / self.value,
                    -self.derivative / (self.value * self.value))

    def __truediv__(self, other):
        return self * as_dual(other).reciprocal()

    def __pow__(self, exponent):
        answer = Dual(F(1))
        for _ in range(exponent):
            answer *= self
        return answer


def as_dual(value):
    return value if isinstance(value, Dual) else Dual(F(value))


def dual_word_sum(multiset, powers_rows, energy_by_mask, score_by_mask):
    total = Dual(F(0))
    retained = 0
    for word in unique_permutations(multiset):
        masks = []
        parity = 0
        allowed = True
        for label in word[:-1]:
            parity ^= 1 << label
            if not parity or energy_by_mask[parity] <= 0:
                allowed = False
                break
            masks.append(parity)
        if not allowed:
            continue
        retained += 1
        for powers in powers_rows:
            value = Dual(F(1))
            for mask, exponent in zip(masks, powers):
                # R(eta)=-1/[E+eta m], differentiated by dual arithmetic.
                value *= (Dual(F(-1)) /
                          Dual(F(energy_by_mask[mask]), F(score_by_mask[mask]))) ** exponent
            total += value
    return total, retained


def section_dual_number_crosscheck():
    rows = []
    for pair_energy in (2, 4, 6):
        energy = {1: F(2), 2: F(2), 3: F(pair_energy)}
        score = {1: F(7, 5), 2: F(-11, 7), 3: F(13, 9)}
        repeated = Dual(F(0))
        retained = 0
        for multiset in ((0, 0, 0, 0, 1, 1), (0, 0, 1, 1, 1, 1)):
            value, count = dual_word_sum(
                multiset, ((1, 1, 1, 1, 1),), energy, score)
            repeated += value
            retained += count
        predicted = sum(PAIR_T6_KERNEL[pair_energy][mask] * score[mask]
                        for mask in score)
        check(retained == 8 and repeated.derivative == predicted,
              f"dual replay matches repeated-pair T6 derivative p={pair_energy}")

        xvalue, xretained = dual_word_sum(
            (0, 0, 1, 1), ((2, 1, 1), (1, 2, 1), (1, 1, 2)),
            energy, score)
        xpredicted = sum(PAIR_X4_KERNEL[pair_energy][mask] * score[mask]
                         for mask in score)
        check(xretained == 4 and xvalue.derivative == xpredicted,
              f"dual replay matches X4 derivative p={pair_energy}")
        rows.append({
            "pair_energy": pair_energy,
            "T6_derivative": ftext(repeated.derivative),
            "X4_derivative": ftext(xvalue.derivative),
        })

    triple_rows = []
    signatures = set()
    for pair_tuple, total in TRIPLE_CLASSES:
        for permuted in set(permutations(pair_tuple)):
            signatures.add((permuted, total))
    for signature in sorted(signatures):
        pair_tuple, total = signature
        energy = {
            1: F(2), 2: F(2), 4: F(2),
            3: F(pair_tuple[0]), 5: F(pair_tuple[1]),
            6: F(pair_tuple[2]), 7: F(total),
        }
        score = {mask: F(3*mask-11, mask+2) for mask in range(1, 8)}
        value, retained = dual_word_sum(
            (0, 0, 1, 1, 2, 2), ((1, 1, 1, 1, 1),), energy, score)
        kernel = triple_kernel_for_signature(signature)
        predicted = sum(kernel[mask] * score[mask] for mask in score)
        check(retained == 60 and value.derivative == predicted,
              f"dual replay matches triple T6 derivative {signature}")
        triple_rows.append({
            "signature": [list(pair_tuple), total],
            "derivative": ftext(value.derivative),
        })

    # The single-intermediate A2/A3 derivatives provide an independent sign
    # anchor for the Kato folds.
    energy_source = Dual(F(2), F(5, 3))
    r = Dual(F(-1)) / energy_source
    check((r**2).derivative == F(-5, 12),
          "dual A2 derivative has negative sign")
    check((r**3).derivative == F(5, 16),
          "dual A3 derivative has positive sign")

    # Synthetic exact assembly checks the literal AO signs
    # T6 - b X4 + b^2 A3 - d A2 with b=-M/2,d=-7M/24.
    m_links = F(52)
    b, d = -m_links/2, -F(7)*m_links/24
    pieces = (Dual(F(-3), F(7, 11)), Dual(F(5), F(-13, 17)),
              r**3, r**2)
    assembled = pieces[0] - b*pieces[1] + b*b*pieces[2] - d*pieces[3]
    predicted = (pieces[0].derivative - b*pieces[1].derivative
                 + b*b*pieces[2].derivative - d*pieces[3].derivative)
    check(assembled.derivative == predicted,
          "dual Kato assembly preserves T6-bX4+b2A3-dA2 signs")
    return {
        "pair_classes": rows,
        "triple_signatures_checked": len(triple_rows),
        "triple_rows_sha256": canonical_hash(triple_rows),
        "resolvent_rule": "d[-1/(E+eta m)]/deta at zero = +m/E^2",
        "fold_formula": "K6'=T6'-b X4'+b^2 A3'-d A2' for pure T2",
        "A2_prime_single_score": "-1/4",
        "A3_prime_single_score": "3/16",
    }


def mask_edges(edges, mask):
    return tuple(edges[label] for label in range(len(edges)) if mask & (1 << label))


def add_scaled(target, vector, scale):
    for i in range(3):
        target[i] += scale * vector[i]


def diagonal_t2_vertex_at_node(occupied, node):
    """Return T6', X4', A3', A2', and the complete K6' T contractions."""
    local = frozenset(LINK_INDEX[link] for link in INCIDENT[node])
    score_cache = {}

    def score(edges):
        key = frozenset(edge for edge in edges if edge in local)
        if key not in score_cache:
            score_cache[key] = local_t_scores(occupied, node, key)
        return score_cache[key]

    pair_sets = 0
    triple_sets = 0
    local_order = tuple(LINK_INDEX[link] for link in INCIDENT[node])
    local_port = {edge: port for port, edge in enumerate(local_order)}
    pattern_coefficients = {
        "T6_prime": defaultdict(F),
        "X4_prime": defaultdict(F),
        "A3_prime": defaultdict(F),
        "A2_prime": defaultdict(F),
    }

    def pattern(edges):
        value = 0
        for edge in edges:
            if edge in local_port:
                value ^= 1 << local_port[edge]
        return value

    def register(label, edges, coefficient):
        local_pattern = pattern(edges)
        # Pure T2 scores vanish for even local parity, so the exact T2
        # calculation need retain only the eight odd patterns.
        if bin(local_pattern).count("1") % 2:
            pattern_coefficients[label][local_pattern] += coefficient

    # Single-edge A_p terms.  Edges outside the source node have zero score.
    for edge in local:
        # R^3 derivative: (-1/E^3)*(-3m/E)=+3m/16 at E=2.
        register("A3_prime", (edge,), F(3, 16))
        # R^2 derivative: (+1/E^2)*(-2m/E)=-m/4 at E=2.
        register("A2_prime", (edge,), F(-1, 4))

    # Every unordered pair is owned once.  Only pairs meeting the source node
    # can contribute to its local source derivative.
    for anchor in sorted(local):
        for other in range(len(LINKS)):
            if other == anchor:
                continue
            left, right = sorted((anchor, other))
            if anchor != min(local.intersection((left, right))):
                continue
            pair_sets += 1
            edges = (left, right)
            p = int(subset_energy(occupied, edges))
            check(p in PAIR_T6_KERNEL, f"two-edge energy {p} is 2,4,or6")
            for mask, coefficient in PAIR_T6_KERNEL[p].items():
                active_edges = mask_edges(edges, mask)
                register("T6_prime", active_edges, coefficient)
            for mask, coefficient in PAIR_X4_KERNEL[p].items():
                active_edges = mask_edges(edges, mask)
                register("X4_prime", active_edges, coefficient)

    # Every unordered triple is owned once.  Again skip triples disjoint from
    # the local source support; all their intermediate scores are zero.
    all_edges = tuple(range(len(LINKS)))
    for anchor in sorted(local):
        others = tuple(edge for edge in all_edges if edge != anchor)
        for first, second in combinations(others, 2):
            triple = tuple(sorted((anchor, first, second)))
            if anchor != min(local.intersection(triple)):
                continue
            triple_sets += 1
            signature = ordered_triple_signature(occupied, triple)
            kernel = TRIPLE_T6_KERNEL[signature]
            for mask, coefficient in kernel.items():
                active_edges = mask_edges(triple, mask)
                register("T6_prime", active_edges, coefficient)

    def contract(label):
        answer = [F(0), F(0), F(0)]
        for local_pattern, coefficient in pattern_coefficients[label].items():
            edges = tuple(local_order[port] for port in range(4)
                          if local_pattern & (1 << port))
            add_scaled(answer, score(edges), coefficient)
        return answer

    direct = contract("T6_prime")
    folded_x = contract("X4_prime")
    a3 = contract("A3_prime")
    a2 = contract("A2_prime")

    m_links = F(len(LINKS))
    b = -m_links / 2
    d = -F(7) * m_links / 24
    complete = [
        direct[i] - b * folded_x[i] + b*b * a3[i] - d * a2[i]
        for i in range(3)
    ]
    bits = tuple(occupation(occupied, edge) for edge in local_order)
    predicted_t6_singleton = universal_singleton_pattern_coefficient(
        bits, 0, len(LINKS))
    predicted_t6_triple = universal_triple_pattern_coefficient(bits, 0)
    predicted_x_singleton = -F(5, 16)*len(LINKS) - F(487, 432)
    for port in range(4):
        check(pattern_coefficients["T6_prime"][1 << port] ==
              predicted_t6_singleton,
              f"literal Q4 T6 singleton pattern port {port} matches rooted census")
        check(pattern_coefficients["T6_prime"][15 ^ (1 << port)] ==
              predicted_t6_triple,
              f"literal Q4 T6 triple pattern port {port} matches rooted census")
        check(pattern_coefficients["X4_prime"][1 << port] ==
              predicted_x_singleton,
              f"literal Q4 X4 singleton pattern port {port} matches rooted census")
    return {
        "node": repr(node),
        "local_locked_word": list(bits),
        "pair_sets_with_local_support": pair_sets,
        "triple_sets_with_local_support": triple_sets,
        "T6_prime": [ftext(x) for x in direct],
        "X4_prime": [ftext(x) for x in folded_x],
        "A3_prime": [ftext(x) for x in a3],
        "A2_prime": [ftext(x) for x in a2],
        "minus_b_X4_prime": [ftext(-b*x) for x in folded_x],
        "b2_A3_prime": [ftext(b*b*x) for x in a3],
        "minus_d_A2_prime": [ftext(-d*x) for x in a2],
        "K6_prime_T_contractions": [ftext(x) for x in complete],
        "odd_local_pattern_coefficients": {
            label: {format(mask, "04b"): ftext(value)
                    for mask, value in sorted(table.items())
                    if bin(mask).count("1") % 2 and value}
            for label, table in pattern_coefficients.items()
        },
        "score_cache_size": len(score_cache),
    }


def canonical_cycle(cell, ports):
    a, b, c = ports
    ab = qadd(qsub(cell, STEPS[b]), STEPS[a])
    cb = qadd(qsub(cell, STEPS[b]), STEPS[c])
    return ((cell, a), (ab, b), (ab, c),
            (cb, a), (cb, b), (cell, c))


def flippable_cycles(occupied):
    rows = []
    seen = set()
    for cell in CELLS:
        for ports in combinations(range(4), 3):
            cycle = canonical_cycle(cell, ports)
            indices = tuple(LINK_INDEX[link] for link in cycle)
            key = frozenset(indices)
            check(key not in seen, "Q4 canonical cycles are unique")
            seen.add(key)
            word = tuple(occupation(occupied, edge) for edge in indices)
            if word in ((1, 0, 1, 0, 1, 0), (0, 1, 0, 1, 0, 1)):
                rows.append(indices)
    check(len(seen) == 256, "Q4 has exactly 256 canonical six-cycles")
    return tuple(rows)


def toggle(occupied, edges):
    answer = set(occupied)
    answer.symmetric_difference_update(edges)
    answer = frozenset(answer)
    check(is_locked(answer), "alternating six-cycle toggle remains locked")
    return answer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", type=int, default=0,
                        help="exploration override: first N Q4 nodes per state")
    parser.add_argument("--states", type=int, default=0,
                        help="exploration override: N deterministic ring-walk states")
    parser.add_argument("--emit", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    universal = section_universal_rooted_census()
    dual = section_dual_number_crosscheck()
    rows = []
    if args.nodes or args.states:
        node_count = args.nodes or 1
        state_count = args.states or 1
        states = []
        occupied = BASE_OCCUPIED
        for step in range(state_count):
            states.append((f"exploration_state_{step}", occupied))
            moves = flippable_cycles(occupied)
            check(bool(moves), f"ring-walk state {step} has a legal move")
            occupied = toggle(occupied, moves[(11 * step + 7) % len(moves)])
        cases = [(label, state, node)
                 for label, state in states for node in NODES[:node_count]]
    else:
        # One literal all-history Q4 representative for each of the six
        # possible central locked words.
        representatives = {}
        for node in NODES:
            bits = tuple(occupation(BASE_OCCUPIED, LINK_INDEX[link])
                         for link in INCIDENT[node])
            representatives.setdefault(bits, node)
        check(len(representatives) == 6,
              "selected Q4 background contains all six central locked words")
        cases = [("base_word_" + "".join(map(str, bits)), BASE_OCCUPIED, node)
                 for bits, node in sorted(representatives.items())]

        # Independently perturb the environment around one fixed source node
        # by legal incident ring toggles and retain three additional local
        # words.  This checks that the result is not an artifact of the base
        # background's periodic surroundings.
        source_node = ("C", (0, 0, 0))
        occupied = BASE_OCCUPIED
        seen_words = {
            tuple(occupation(occupied, LINK_INDEX[link])
                  for link in INCIDENT[source_node])
        }
        source_edges = {LINK_INDEX[link] for link in INCIDENT[source_node]}
        queue = deque(((occupied, 0),))
        visited = {occupied}
        while queue and len(seen_words) < 4:
            current, depth = queue.popleft()
            if depth >= 4:
                continue
            for move in flippable_cycles(current):
                if len(source_edges.intersection(move)) != 2:
                    continue
                trial = toggle(current, move)
                if trial in visited:
                    continue
                visited.add(trial)
                queue.append((trial, depth + 1))
                bits = tuple(occupation(trial, LINK_INDEX[link])
                             for link in INCIDENT[source_node])
                if bits not in seen_words:
                    seen_words.add(bits)
                    cases.append((
                        f"incident_walk_depth_{depth+1}_word_" + "".join(map(str, bits)),
                        trial, source_node))
                    if len(seen_words) >= 4:
                        break
        check(len(seen_words) >= 4,
              "incident ring walk realizes at least four local locked words")

    for case_index, (label, state, node) in enumerate(cases):
        row = diagonal_t2_vertex_at_node(state, node)
        row["case"] = label
        row["case_index"] = case_index
        check(row["T6_prime"] == ["0", "0", "0"],
              f"literal case {label} direct T6 pure-T2 contraction vanishes")
        check(row["X4_prime"] == ["0", "0", "0"],
              f"literal case {label} X4 pure-T2 contraction vanishes")
        check(row["A3_prime"] == ["0", "0", "0"] and
              row["A2_prime"] == ["0", "0", "0"],
              f"literal case {label} A3/A2 pure-T2 contractions vanish")
        check(row["K6_prime_T_contractions"] == ["0", "0", "0"],
              f"literal case {label} complete K6 pure-T2 contraction vanishes")
        rows.append(row)

    pair_kernel_payload = {
        str(energy): {
            "T6": {str(mask): ftext(value) for mask, value in sorted(PAIR_T6_KERNEL[energy].items())},
            "X4": {str(mask): ftext(value) for mask, value in sorted(PAIR_X4_KERNEL[energy].items())},
        }
        for energy in (2, 4, 6)
    }
    triple_kernel_payload = {
        repr((pair_tuple, total)): {
            str(mask): ftext(value) for mask, value in sorted(
                triple_kernel_for_signature((pair_tuple, total)).items())
        }
        for pair_tuple, total in TRIPLE_CLASSES
    }
    result = {
        "schema": "GL6CN_COMPLETE_DIAGONAL_H6_T2_FIRST_SOURCE_V001",
        "domain": "every locked degree-two state of every simple degree-four girth-at-least-six bipartite incidence parent; literal replay on Q4",
        "source_convention": "H(j)=U_d D+hW+sum_v j_v^T M_v; pure T2 has P M_T P=0",
        "kato_formula": "K6=T6-bX4+b^2A3-dA2; b=-M/2; d=-7M/24",
        "differentiated_formula": "K6_T'=T6_T'-bX4_T'+b^2A3_T'-dA2_T'",
        "dual_number_crosscheck": dual,
        "history_derivative_kernels": {
            "pair": pair_kernel_payload,
            "triple_canonical": triple_kernel_payload,
        },
        "universal_rooted_census": universal,
        "literal_Q4_cases": rows,
        "literal_Q4_case_count": len(rows),
        "all_complete_t_contractions_zero": all(
            row["K6_prime_T_contractions"] == ["0", "0", "0"] for row in rows),
        "all_individual_kato_terms_t_contractions_zero": all(
            row[name] == ["0", "0", "0"]
            for row in rows for name in ("T6_prime", "X4_prime", "A3_prime", "A2_prime")),
        "complete_diagonal_h6_T2_vertex": "zero pointwise",
        "integrated_corollary": "through h6, the complete pure-T2 first-source effective operator is the GL6CH off-diagonal six-cycle writer",
        "open": [
            "source-second contacts, including mixed two-site or nonuniform source-source derivatives",
            "higher-order first-source vertices",
            "stationary phase and complete contact-plus-spectral Hessian",
            "record authentication, bulk/refinement, Ricci, gravity, and G",
        ],
        "row_sha256": canonical_hash(rows),
        "checks": len(CHECKS),
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.emit:
        (HERE / "EXACT_LEDGER.json").write_text(rendered)
    elif not (args.nodes or args.states):
        if (HERE / "EXACT_LEDGER.json").read_text() != rendered:
            raise AssertionError("frozen exact ledger does not match replay")
    if args.json or args.nodes or args.states:
        print(rendered, end="")
    else:
        print(f"PASS__GL6CN_COMPLETE_DIAGONAL_H6_T2_FIRST_SOURCE__{len(CHECKS)}/{len(CHECKS)}")
        print("DUAL_RESOLVENT_AND_KATO_SIGNS=PASS;ALL_8_PLUS_60_RETURN_WORDS_DIFFERENTIATED")
        print("ROOTED_CENSUS=ALL_6_LOCKED_WORDS;ODD_PATTERN_COEFFICIENTS_PORT_UNIFORM")
        print(f"Q4_LITERAL_CASES={len(rows)};BASE_ALL6_PLUS_INCIDENT_RING_WALK;ALL_T6_X4_A3_A2_T2_ZERO")
        print("DIAGONAL_H6_T2=ZERO_POINTWISE;COMPLETE_FIRST_T2_THROUGH_H6=GL6CH_WRITER")
        print("CONTACT_HIGHER_ORDER_PHASE_RECORD_BULK_RICCI_GRAVITY_G=OPEN")


if __name__ == "__main__":
    main()
