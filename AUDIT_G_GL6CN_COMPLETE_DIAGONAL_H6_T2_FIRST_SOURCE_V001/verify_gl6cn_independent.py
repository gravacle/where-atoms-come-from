#!/usr/bin/env python3
"""Independent hostile replay for GL6CN.

This implementation does not import or execute the author derivation.  It
rebuilds multiset histories, dual-number signs, rooted local combinatorics,
and four independent degree-four/girth-six circulant parents from definitions.
"""

from __future__ import annotations

import json
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path


HERE = Path(__file__).resolve().parent
PAIR_ORDER = tuple(combinations(range(4), 2))
T_BASIS = (
    (1, 0, 0, 0, 0, -1),
    (0, 1, 0, 0, -1, 0),
    (0, 0, 1, -1, 0, 0),
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


def unique_words(multiset):
    return tuple(sorted(set(permutations(multiset))))


def parity_masks(word):
    mask = frozenset()
    out = []
    for label in word[:-1]:
        mask = mask.symmetric_difference((label,))
        out.append(mask)
    final = mask.symmetric_difference((word[-1],))
    CHECK.equal(final, frozenset(), "history returns diagonally")
    return tuple(out)


def t6_kernel(multiset, energies):
    kernel = {mask: F(0) for mask in energies}
    retained = 0
    for word in unique_words(multiset):
        masks = parity_masks(word)
        if any(not mask for mask in masks):
            continue
        retained += 1
        denominator = F(1)
        for mask in masks:
            denominator /= energies[mask]
        for mask in masks:
            kernel[mask] += denominator / energies[mask]
    return {mask: value for mask, value in kernel.items() if value}, retained


def x4_kernel(pair_energy):
    energies = {
        frozenset((0,)): F(2),
        frozenset((1,)): F(2),
        frozenset((0, 1)): F(pair_energy),
    }
    kernel = {mask: F(0) for mask in energies}
    retained_words = 0
    for word in unique_words((0, 0, 1, 1)):
        masks = parity_masks(word)
        if any(not mask for mask in masks):
            continue
        retained_words += 1
        for powers in ((2, 1, 1), (1, 2, 1), (1, 1, 2)):
            base = F(1)
            for mask, power in zip(masks, powers):
                base *= F((-1) ** power, energies[mask] ** power)
            for mask, power in zip(masks, powers):
                kernel[mask] += -base * power / energies[mask]
    return {mask: value for mask, value in kernel.items() if value}, retained_words


def all_history_kernels():
    expected_pair = {
        2: ((F(3, 16), F(3, 16), F(1, 4)),
            (F(-1, 2), F(-1, 2), F(-1, 2))),
        4: ((F(3, 64), F(3, 64), F(1, 32)),
            (F(-7, 32), F(-7, 32), F(-3, 32))),
        6: ((F(1, 48), F(1, 48), F(1, 108)),
            (F(-5, 36), F(-5, 36), F(-1, 27))),
    }
    pair_kernels = {}
    masks2 = (frozenset((0,)), frozenset((1,)), frozenset((0, 1)))
    for energy, expected in expected_pair.items():
        energies = {masks2[0]: F(2), masks2[1]: F(2), masks2[2]: F(energy)}
        left, retained_a = t6_kernel((0, 0, 0, 0, 1, 1), energies)
        right, retained_b = t6_kernel((0, 0, 1, 1, 1, 1), energies)
        combined = {mask: left.get(mask, 0) + right.get(mask, 0) for mask in masks2}
        x4, retained_x = x4_kernel(energy)
        CHECK.equal((retained_a, retained_b, retained_x), (4, 4, 4),
                    f"pair energy {energy} retained history counts")
        CHECK.equal(tuple(combined[mask] for mask in masks2), expected[0],
                    f"pair energy {energy} T6 derivative kernel")
        CHECK.equal(tuple(x4[mask] for mask in masks2), expected[1],
                    f"pair energy {energy} X4 derivative kernel")
        pair_kernels[energy] = (combined, x4)

    canonical = {
        ((4, 4, 4), 6): (F(1, 8), F(1, 8), F(3, 64), F(1, 8), F(3, 64), F(3, 64), F(1, 64)),
        ((2, 4, 4), 4): (F(17, 64), F(17, 64), F(1, 4), F(3, 16), F(5, 64), F(5, 64), F(1, 16)),
        ((4, 4, 6), 8): (F(5, 48), F(49, 576), F(7, 192), F(49, 576), F(7, 192), F(1, 54), F(1, 144)),
        ((2, 2, 6), 4): (F(7, 16), F(19, 72), F(5, 16), F(19, 72), F(5, 16), F(19, 432), F(49, 576)),
        ((2, 2, 4), 2): (F(5, 8), F(29, 64), F(1, 2), F(29, 64), F(1, 2), F(9, 64), F(25, 64)),
        ((2, 4, 6), 6): (F(41, 192), F(19, 108), F(13, 72), F(217, 1728), F(35, 576), F(19, 648), F(121, 5184)),
        ((4, 6, 6), 10): (F(41, 576), F(41, 576), F(9, 320), F(7, 120), F(2, 135), F(2, 135), F(49, 14400)),
    }
    masks3 = (
        frozenset((0,)), frozenset((1,)), frozenset((0, 1)),
        frozenset((2,)), frozenset((0, 2)), frozenset((1, 2)),
        frozenset((0, 1, 2)),
    )
    triple_kernels = {}
    labelled_signatures = set()
    for (pair_tuple, triple_energy), expected in canonical.items():
        energies = {
            frozenset((0,)): F(2), frozenset((1,)): F(2),
            frozenset((2,)): F(2), frozenset((0, 1)): F(pair_tuple[0]),
            frozenset((0, 2)): F(pair_tuple[1]),
            frozenset((1, 2)): F(pair_tuple[2]),
            frozenset((0, 1, 2)): F(triple_energy),
        }
        kernel, retained = t6_kernel((0, 0, 1, 1, 2, 2), energies)
        CHECK.equal(retained, 60, f"triple {pair_tuple,triple_energy} history count")
        CHECK.equal(tuple(kernel[mask] for mask in masks3), expected,
                    f"triple {pair_tuple,triple_energy} kernel")
        triple_kernels[(pair_tuple, triple_energy)] = kernel
        for relabel in permutations(range(3)):
            pair_values = (
                energies[frozenset((relabel[0], relabel[1]))],
                energies[frozenset((relabel[0], relabel[2]))],
                energies[frozenset((relabel[1], relabel[2]))],
            )
            labelled_signatures.add((pair_values, F(triple_energy)))
    CHECK.equal(len(labelled_signatures), 22, "twenty-two labelled triple signatures")
    return pair_kernels, triple_kernels, labelled_signatures


class Dual:
    def __init__(self, value, derivative=0):
        self.v = F(value)
        self.d = F(derivative)

    def __add__(self, other):
        other = other if isinstance(other, Dual) else Dual(other)
        return Dual(self.v + other.v, self.d + other.d)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v, -self.d)

    def __sub__(self, other):
        return self + (-other if isinstance(other, Dual) else -Dual(other))

    def __mul__(self, other):
        other = other if isinstance(other, Dual) else Dual(other)
        return Dual(self.v * other.v, self.d * other.v + self.v * other.d)

    __rmul__ = __mul__

    def reciprocal(self):
        return Dual(1 / self.v, -self.d / (self.v * self.v))

    def __truediv__(self, other):
        other = other if isinstance(other, Dual) else Dual(other)
        return self * other.reciprocal()

    def __pow__(self, power):
        answer = Dual(1)
        for _ in range(power):
            answer = answer * self
        return answer


def dual_t6(multiset, energies, active_mask):
    answer = Dual(0)
    for word in unique_words(multiset):
        masks = parity_masks(word)
        if any(not mask for mask in masks):
            continue
        term = Dual(1)
        for mask in masks:
            term *= -Dual(1) / Dual(energies[mask], int(mask == active_mask))
        answer += term
    return answer.d


def dual_x4(pair_energy, active_mask):
    energies = {
        frozenset((0,)): F(2), frozenset((1,)): F(2),
        frozenset((0, 1)): F(pair_energy),
    }
    answer = Dual(0)
    for word in unique_words((0, 0, 1, 1)):
        masks = parity_masks(word)
        if any(not mask for mask in masks):
            continue
        for powers in ((2, 1, 1), (1, 2, 1), (1, 1, 2)):
            term = Dual(1)
            for mask, power in zip(masks, powers):
                resolvent = -Dual(1) / Dual(energies[mask], int(mask == active_mask))
                term *= resolvent ** power
            answer += term
    return answer.d


def dual_crosscheck(pair_kernels, labelled_signatures):
    masks2 = (frozenset((0,)), frozenset((1,)), frozenset((0, 1)))
    for energy, (t6, x4) in pair_kernels.items():
        for mask in masks2:
            got = (dual_t6((0, 0, 0, 0, 1, 1),
                           {masks2[0]: F(2), masks2[1]: F(2), masks2[2]: F(energy)}, mask)
                   + dual_t6((0, 0, 1, 1, 1, 1),
                             {masks2[0]: F(2), masks2[1]: F(2), masks2[2]: F(energy)}, mask))
            CHECK.equal(got, t6[mask], f"dual pair T6 {energy} {mask}")
            CHECK.equal(dual_x4(energy, mask), x4[mask],
                        f"dual pair X4 {energy} {mask}")

    masks3 = tuple(frozenset(x) for n in (1, 2, 3)
                   for x in combinations(range(3), n))
    for pair_values, triple_energy in sorted(labelled_signatures):
        energies = {
            frozenset((0,)): F(2), frozenset((1,)): F(2),
            frozenset((2,)): F(2), frozenset((0, 1)): pair_values[0],
            frozenset((0, 2)): pair_values[1],
            frozenset((1, 2)): pair_values[2],
            frozenset((0, 1, 2)): triple_energy,
        }
        analytic, _ = t6_kernel((0, 0, 1, 1, 2, 2), energies)
        for mask in masks3:
            CHECK.equal(dual_t6((0, 0, 1, 1, 2, 2), energies, mask),
                        analytic[mask], f"dual triple {pair_values,triple_energy} {mask}")

    # Independent two-flip fold signs at E=2.
    eta = Dual(2, 1)
    r = -Dual(1) / eta
    CHECK.equal((r ** 2).d, F(-1, 4), "A2 derivative sign")
    CHECK.equal((r ** 3).d, F(3, 16), "A3 derivative sign")


def make_circulant(n):
    shifts = (0, 1, 3, 7)
    edges = tuple((parent, (parent + shifts[port]) % n, port)
                  for parent in range(n) for port in range(4))
    incidence = {("P", x): [] for x in range(n)}
    incidence.update({("C", x): [] for x in range(n)})
    for index, (parent, child, _port) in enumerate(edges):
        incidence[("P", parent)].append(index)
        incidence[("C", child)].append(index)
    CHECK.true(all(len(items) == 4 for items in incidence.values()),
               f"circulant {n} is degree four")
    CHECK.equal(len(set(edges)), 4 * n, f"circulant {n} is simple")
    # No two parents share two children: equivalent to no four-cycle.
    child_sets = [set(edges[e][1] for e in incidence[("P", x)]) for x in range(n)]
    CHECK.true(all(len(child_sets[x] & child_sets[y]) <= 1
                   for x in range(n) for y in range(x + 1, n)),
               f"circulant {n} has girth at least six")
    endpoints = tuple((('P', p), ('C', c)) for p, c, _ in edges)
    return edges, endpoints, incidence


def adjacent(endpoints, left, right):
    return bool(set(endpoints[left]) & set(endpoints[right]))


def shared_vertex(endpoints, left, right):
    shared = set(endpoints[left]) & set(endpoints[right])
    CHECK.equal(len(shared), 1, "adjacent links have one shared vertex")
    return next(iter(shared))


def energy_of_subset(subset, endpoints, signs):
    delta = {}
    for edge in subset:
        for vertex in endpoints[edge]:
            delta[vertex] = delta.get(vertex, 0) + signs[edge]
    return sum(value * value for value in delta.values())


def census_for_root(edges, endpoints, incidence, signs, root):
    others = [edge for edge in range(len(edges)) if edge != root]
    counts = {
        "matching": 0,
        "one_opp_root_adjacent": 0,
        "one_opp_root_disjoint": 0,
        "one_equal_root_adjacent": 0,
        "one_equal_root_disjoint": 0,
        "star_root_minority": 0,
        "star_root_majority": 0,
        "path_middle_both_opp": 0,
        "path_middle_one_equal": 0,
        "path_middle_both_equal": 0,
        "path_end_both_opp": 0,
        "path_end_first_opp": 0,
        "path_end_first_equal": 0,
        "path_end_both_equal": 0,
    }
    for left, right in combinations(others, 2):
        ar_l = adjacent(endpoints, root, left)
        ar_r = adjacent(endpoints, root, right)
        alr = adjacent(endpoints, left, right)
        if not ar_l and not ar_r:
            if not alr:
                counts["matching"] += 1
            else:
                relation = "opp" if signs[left] != signs[right] else "equal"
                counts[f"one_{relation}_root_disjoint"] += 1
        elif ar_l ^ ar_r:
            near, far = (left, right) if ar_l else (right, left)
            first = "opp" if signs[root] != signs[near] else "equal"
            if alr:
                second = "opp" if signs[near] != signs[far] else "equal"
                if first == second == "opp":
                    key = "path_end_both_opp"
                elif first == "opp":
                    key = "path_end_first_opp"
                elif second == "opp":
                    key = "path_end_first_equal"
                else:
                    key = "path_end_both_equal"
                counts[key] += 1
            else:
                counts[f"one_{first}_root_adjacent"] += 1
        else:
            if alr:
                # All three share the same root endpoint in a girth>=6 graph.
                common_a = shared_vertex(endpoints, root, left)
                common_b = shared_vertex(endpoints, root, right)
                CHECK.equal(common_a, common_b, "rooted triangle is a star")
                minority = signs[root] != signs[left] == signs[right]
                counts["star_root_minority" if minority else "star_root_majority"] += 1
            else:
                relations = (signs[root] != signs[left], signs[root] != signs[right])
                counts["path_middle_both_opp" if all(relations) else
                       "path_middle_both_equal" if not any(relations) else
                       "path_middle_one_equal"] += 1
    return counts


def expected_census(m):
    return {
        "matching": (m * m - 21 * m + 116) // 2,
        "one_opp_root_adjacent": 4 * (m - 10),
        "one_opp_root_disjoint": 2 * (m - 10),
        "one_equal_root_adjacent": 2 * (m - 10),
        "one_equal_root_disjoint": m - 10,
        "star_root_minority": 2,
        "star_root_majority": 4,
        "path_middle_both_opp": 4,
        "path_middle_one_equal": 4,
        "path_middle_both_equal": 1,
        "path_end_both_opp": 8,
        "path_end_first_opp": 4,
        "path_end_first_equal": 4,
        "path_end_both_equal": 2,
    }


def word_vector(occupied, flipped_ports):
    z = tuple((-1 if port in occupied else 1) *
              (-1 if port in flipped_ports else 1) for port in range(4))
    return tuple(F(z[a] * z[b]) for a, b in PAIR_ORDER)


def local_mask(subset, support_port):
    return frozenset(support_port[edge] for edge in subset if edge in support_port)


def aggregate_parent(n, occupied, pair_kernels, labelled_signatures):
    edges, endpoints, incidence = make_circulant(n)
    m = len(edges)
    signs = tuple(F(-1 if port in occupied else 1) for _, _, port in edges)
    source = ("P", 0)
    support = tuple(incidence[source])
    support_port = {edge: edges[edge][2] for edge in support}

    # Rooted census independently checked for all four source ports.
    for root in support:
        got = census_for_root(edges, endpoints, incidence, signs, root)
        CHECK.equal(got, expected_census(m),
                    f"rooted census n={n} word={occupied} root={edges[root][2]}")
        CHECK.equal(sum(got.values()), (m - 1) * (m - 2) // 2,
                    "rooted census partitions every pair of other links")

    agg_t6 = {}
    agg_x4 = {}
    all_edges = range(m)
    pair_sets = tuple(combinations(all_edges, 2))
    for left, right in pair_sets:
        if left not in support_port and right not in support_port:
            continue
        p = energy_of_subset((left, right), endpoints, signs)
        CHECK.true(p in pair_kernels, "pair energy is 2,4,6")
        t6, x4 = pair_kernels[p]
        label_to_edge = {0: left, 1: right}
        for selected_mask, coefficient in t6.items():
            mask = local_mask((label_to_edge[x] for x in selected_mask), support_port)
            agg_t6[mask] = agg_t6.get(mask, F(0)) + coefficient
        for selected_mask, coefficient in x4.items():
            mask = local_mask((label_to_edge[x] for x in selected_mask), support_port)
            agg_x4[mask] = agg_x4.get(mask, F(0)) + coefficient

    agg_triple = {}
    triple_count = 0
    triple_cache = {}
    # Every triple meeting support, without materializing all C(M,3).
    for selected in combinations(all_edges, 3):
        if not any(edge in support_port for edge in selected):
            continue
        triple_count += 1
        energies = {}
        for size in (1, 2, 3):
            for labels in combinations(range(3), size):
                subset = tuple(selected[index] for index in labels)
                energies[frozenset(labels)] = F(energy_of_subset(subset, endpoints, signs))
        signature = tuple(energies[mask] for mask in (
            frozenset((0,)), frozenset((1,)), frozenset((0, 1)),
            frozenset((2,)), frozenset((0, 2)), frozenset((1, 2)),
            frozenset((0, 1, 2)),
        ))
        if ((signature[2], signature[4], signature[5]), signature[6]) not in labelled_signatures:
            raise AssertionError(f"unclassified triple energy signature: {signature}")
        if signature not in triple_cache:
            triple_cache[signature] = t6_kernel((0, 0, 1, 1, 2, 2), energies)
        kernel, retained = triple_cache[signature]
        CHECK.equal(retained, 60, "every triple retains sixty Q-only histories")
        for selected_mask, coefficient in kernel.items():
            mask = local_mask((selected[index] for index in selected_mask), support_port)
            agg_triple[mask] = agg_triple.get(mask, F(0)) + coefficient

    expected_pair_count = (m * (m - 1) // 2) - ((m - 4) * (m - 5) // 2)
    expected_triple_count = (m * (m - 1) * (m - 2) // 6) - ((m - 4) * (m - 5) * (m - 6) // 6)
    CHECK.equal(sum(1 for pair in pair_sets if any(x in support_port for x in pair)),
                expected_pair_count, "all source-touching pairs included")
    CHECK.equal(triple_count, expected_triple_count,
                "all source-touching triples included")

    t6_total = {mask: agg_t6.get(mask, 0) + agg_triple.get(mask, 0)
                for mask in set(agg_t6) | set(agg_triple)}
    c_t6_1 = F(15, 128) * m * m + F(3049, 3456) * m + F(8653, 4800)
    c_t6_3 = F(49, 576)
    c_x4_1 = -F(5, 16) * m - F(487, 432)
    singleton_masks = tuple(frozenset((port,)) for port in range(4))
    triple_masks = tuple(frozenset(set(range(4)) - {port}) for port in range(4))
    for mask in singleton_masks:
        CHECK.equal(t6_total[mask], c_t6_1,
                    f"T6 singleton port-uniform n={n} word={occupied} mask={mask}")
        CHECK.equal(agg_x4[mask], c_x4_1,
                    f"X4 singleton port-uniform n={n} word={occupied} mask={mask}")
    for mask in triple_masks:
        CHECK.equal(t6_total[mask], c_t6_3,
                    f"T6 complement-triple uniform n={n} word={occupied} mask={mask}")
        CHECK.equal(agg_x4.get(mask, 0), F(0), "X4 has no triple mask")

    # Even masks are T2 dark.  Odd masks cancel as full vectors.
    words_single = tuple(word_vector(occupied, mask) for mask in singleton_masks)
    words_triple = tuple(word_vector(occupied, mask) for mask in triple_masks)
    CHECK.equal(tuple(sum((word[i] for word in words_single), F(0)) for i in range(6)),
                (F(0),) * 6, "singleton full-vector sum zero")
    CHECK.equal(words_single, words_triple, "complement triples equal singletons")
    for mask in t6_total:
        if len(mask) % 2 == 0:
            vector = word_vector(occupied, mask)
            CHECK.equal(tuple(sum(F(t[i]) * vector[i] for i in range(6))
                              for t in T_BASIS), (F(0),) * 3,
                        "even local mask T2 dark")

    # Fold combination for every port; signs are K6'=T6'-bX4'+b^2A3'-dA2'.
    b = -F(m, 2)
    d = -F(7 * m, 24)
    c_k6_1 = c_t6_1 - b * c_x4_1 + b * b * F(3, 16) - d * F(-1, 4)
    CHECK.equal(c_k6_1,
                F(1, 128) * m * m + F(283, 1152) * m + F(8653, 4800),
                "Kato folded singleton polynomial")
    total_vector = tuple(
        sum((c_k6_1 * words_single[a][i] + c_t6_3 * words_triple[a][i]
             for a in range(4)), F(0)) for i in range(6)
    )
    CHECK.equal(total_vector, (F(0),) * 6,
                "complete diagonal first-source vector cancels pointwise")
    return {
        "n": n, "M": m, "occupied_ports": tuple(sorted(occupied)),
        "pairs_touching_source": expected_pair_count,
        "triples_touching_source": expected_triple_count,
        "T6_singleton": c_t6_1,
        "T6_complement_triple": c_t6_3,
        "X4_singleton": c_x4_1,
        "K6_singleton": c_k6_1,
    }


def qtext(value):
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def encode(value):
    if isinstance(value, F):
        return qtext(value)
    if isinstance(value, dict):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def main():
    pair_kernels, triple_kernels, labelled = all_history_kernels()
    dual_crosscheck(pair_kernels, labelled)

    # Symbolic census identity independent of finite witnesses.
    for m in (40, 68, 76, 92, 116, 256):
        census = expected_census(m)
        CHECK.equal(sum(census.values()), (m - 1) * (m - 2) // 2,
                    f"symbolic rooted census identity M={m}")
        CHECK.equal(census["matching"],
                    (m - 7) * (m - 8) // 2 - (3 * m - 30),
                    f"root-disjoint matching identity M={m}")

    witnesses = []
    for n in (17, 19, 23, 29):
        for occupied in combinations(range(4), 2):
            witnesses.append(aggregate_parent(n, frozenset(occupied), pair_kernels,
                                                labelled))

    # Literal target Q4 support-set counts follow independently from M=256.
    CHECK.equal(256 * 255 // 2 - 252 * 251 // 2, 1014,
                "Q4 pair support count")
    CHECK.equal(256 * 255 * 254 // 6 - 252 * 251 * 250 // 6, 128020,
                "Q4 triple support count")

    result = {
        "verdict": "PASS",
        "independence": "author derivation was neither imported nor executed",
        "pair_kernel_classes": 3,
        "triple_canonical_classes": len(triple_kernels),
        "triple_labelled_signatures": len(labelled),
        "dual_number_rule": "d[-1/(E+eta m)]/deta=+m/E^2",
        "fold_signs": "K6'=T6'-bX4'+b^2A3'-dA2'; A2'=-m/4; A3'=+3m/16",
        "rooted_domain": "simple degree-four bipartite girth-at-least-six; exact incidence classification uses no Q4-specific fact",
        "independent_parent_witnesses": witnesses,
        "pointwise_mechanism": "port-uniform singleton and complementary-triple coefficients multiply identical word families whose four-word full-vector sum is zero",
        "integrated_corollary": "given the separately sealed lower-order and GL6CH offdiagonal completeness claims, the sole pure-T2 first-source operator through h6 is the GL6CH cycle writer",
        "ceiling": "source-second contacts, h8+, stationary phase, record authentication, bulk/refinement, Ricci, gravity, and G remain open",
        "checks": CHECK.total,
    }
    (HERE / "INDEPENDENT_RESULT.json").write_text(
        json.dumps(encode(result), indent=2, sort_keys=True) + "\n")
    print(f"PASS__GL6CN_INDEPENDENT_HOSTILE_SCIENCE__{CHECK.total}/{CHECK.total}")
    print("HISTORY_KERNELS=3_PAIR;7_CANONICAL_TRIPLE;22_LABELLED_TRIPLE")
    print("DUAL_SIGNS=RESOLVENT_PLUS;A2_MINUS_1_4;A3_PLUS_3_16;KATO_SIGNS_PASS")
    print("ROOTED_CENSUS=EXHAUSTIVE_DEGREE4_GIRTH6;24_INDEPENDENT_PARENT_WORD_WITNESSES")
    print("POINTWISE_T2=PORT_UNIFORM_SINGLETON_PLUS_COMPLEMENT_TRIPLE_CANCEL")
    print("INTEGRATED_THROUGH_H6=OFFDIAGONAL_GL6CH_WRITER_ONLY_CONDITIONAL_ON_PINNED_INPUTS")
    print("SOURCE_SECOND_H8_PHASE_RECORD_BULK_RICCI_GRAVITY_G_OPEN")


if __name__ == "__main__":
    main()
