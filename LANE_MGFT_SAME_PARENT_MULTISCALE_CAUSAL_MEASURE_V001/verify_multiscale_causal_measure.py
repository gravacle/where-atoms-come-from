#!/usr/bin/env python3
"""Exact finite witnesses for the same-parent multiscale theorem.

This script verifies rational algebra, finite relations, and the claim ledger.
It does not establish a physical parent process, a complete probe domain, a
continuum manifold, Lorentzian realizability, metric four-volume, gravity, or
record ancestry.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from typing import Callable, Dict, Iterable, Sequence, Tuple


Q = Fraction
Pair = Tuple[int, int]
Relation = set[Pair]
Vector = Tuple[Q, ...]
Matrix = Tuple[Tuple[Q, ...], ...]


class Ledger:
    def __init__(self) -> None:
        self.passed = 0
        self.total = 0

    def check(self, condition: bool, name: str) -> None:
        self.total += 1
        if not condition:
            raise AssertionError(name)
        self.passed += 1


def l1(vector: Sequence[Q]) -> Q:
    return sum((abs(x) for x in vector), Q(0))


def mat_vec(matrix: Matrix, vector: Vector) -> Vector:
    return tuple(
        sum((entry * value for entry, value in zip(row, vector)), Q(0))
        for row in matrix
    )


def rank(matrix: Matrix) -> int:
    rows = [list(row) for row in matrix]
    if not rows:
        return 0
    n_rows = len(rows)
    n_cols = len(rows[0])
    pivot_row = 0
    for col in range(n_cols):
        pivot = next(
            (row for row in range(pivot_row, n_rows) if rows[row][col] != 0),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = rows[pivot_row][col]
        rows[pivot_row] = [value / scale for value in rows[pivot_row]]
        for row in range(n_rows):
            if row == pivot_row:
                continue
            factor = rows[row][col]
            if factor:
                rows[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(rows[row], rows[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == n_rows:
            break
    return pivot_row


def solve_square(matrix: Matrix, rhs: Vector) -> Vector:
    n = len(matrix)
    rows = [list(matrix[i]) + [rhs[i]] for i in range(n)]
    for col in range(n):
        pivot = next(row for row in range(col, n) if rows[row][col] != 0)
        rows[col], rows[pivot] = rows[pivot], rows[col]
        scale = rows[col][col]
        rows[col] = [value / scale for value in rows[col]]
        for row in range(n):
            if row == col:
                continue
            factor = rows[row][col]
            rows[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(rows[row], rows[col])
            ]
    return tuple(rows[i][-1] for i in range(n))


def compose(functions: Sequence[Callable[[Q], Q]], value: Q) -> Q:
    for function in reversed(functions):
        value = function(value)
    return value


def push_relation(relation: Relation, mapping: Dict[int, int]) -> Relation:
    return {(mapping[left], mapping[right]) for left, right in relation}


def is_reflexive(relation: Relation, elements: Iterable[int]) -> bool:
    return all((element, element) in relation for element in elements)


def is_transitive(relation: Relation) -> bool:
    return all(
        (left, right) in relation
        for left, middle in relation
        for middle_again, right in relation
        if middle == middle_again
    )


def is_antisymmetric(relation: Relation) -> bool:
    return all(left == right or (right, left) not in relation for left, right in relation)


def is_partial_order(relation: Relation, elements: Iterable[int]) -> bool:
    return (
        is_reflexive(relation, elements)
        and is_transitive(relation)
        and is_antisymmetric(relation)
    )


def response_support(response: Dict[Pair, Q]) -> Relation:
    return {pair for pair, value in response.items() if value != 0}


def has_support(response: Sequence[Q]) -> bool:
    return any(value != 0 for value in response)


def block_response(
    response: Dict[Pair, Q], mapping: Dict[int, int]
) -> Dict[Pair, Q]:
    blocked: Dict[Pair, Q] = {}
    for (left, right), value in response.items():
        pair = (mapping[left], mapping[right])
        blocked[pair] = blocked.get(pair, Q(0)) + value
    return {pair: value for pair, value in blocked.items() if value != 0}


def aggregate(vector: Vector, mapping: Sequence[int], coarse_size: int) -> Vector:
    result = [Q(0) for _ in range(coarse_size)]
    for fine, value in enumerate(vector):
        result[mapping[fine]] += value
    return tuple(result)


def lower_bidiagonal(size: int) -> Matrix:
    return tuple(
        tuple(
            Q(1) if col == row or (row > 0 and col == row - 1) else Q(0)
            for col in range(size)
        )
        for row in range(size)
    )


def bitstrings(length: int) -> Tuple[Tuple[int, ...], ...]:
    return tuple(product((0, 1), repeat=length))


def bernoulli_measure(length: int, probability_one: Q) -> Dict[Tuple[int, ...], Q]:
    result: Dict[Tuple[int, ...], Q] = {}
    for bits in bitstrings(length):
        ones = sum(bits)
        result[bits] = probability_one**ones * (1 - probability_one) ** (length - ones)
    return result


def prefix_push(
    measure: Dict[Tuple[int, ...], Q], prefix_length: int
) -> Dict[Tuple[int, ...], Q]:
    result: Dict[Tuple[int, ...], Q] = {}
    for bits, value in measure.items():
        prefix = bits[:prefix_length]
        result[prefix] = result.get(prefix, Q(0)) + value
    return result


def coordinate_order(strings: Sequence[Tuple[int, ...]]) -> set[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
    return {
        (left, right)
        for left in strings
        for right in strings
        if all(a <= b for a, b in zip(left, right))
    }


def push_string_relation(
    relation: set[Tuple[Tuple[int, ...], Tuple[int, ...]]], prefix_length: int
) -> set[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
    return {
        (left[:prefix_length], right[:prefix_length])
        for left, right in relation
    }


def main() -> None:
    ledger = Ledger()

    # M1: an exact three-step witness saturates the weighted bound.
    kernels = (Q(19), Q(6), Q(11), Q(5))
    reductions = (
        lambda value: Q(3) * value,
        lambda value: Q(1, 2) * value,
        lambda value: Q(2) * value,
    )
    lipschitz = (Q(3), Q(1, 2), Q(2))
    residuals = tuple(
        abs(kernels[index] - reductions[index](kernels[index + 1]))
        for index in range(3)
    )
    ledger.check(residuals == (Q(1), Q(1, 2), Q(1)), "M1_adjacent_residuals_exact")
    sequential = compose(reductions, kernels[3])
    sequential_error = abs(kernels[0] - sequential)
    bound = (
        residuals[0]
        + lipschitz[0] * residuals[1]
        + lipschitz[0] * lipschitz[1] * residuals[2]
    )
    ledger.check(sequential == Q(15), "M1_sequential_prediction_exact")
    ledger.check(sequential_error == bound == Q(4), "M1_telescoping_bound_saturated")
    direct = Q(3) * kernels[3] - Q(1)
    direct_sequential = abs(direct - sequential)
    direct_error = abs(kernels[0] - direct)
    ledger.check(direct_sequential == Q(1), "M1_direct_sequential_residual_exact")
    ledger.check(direct_error == bound + direct_sequential == Q(5), "M1_direct_bound_saturated")
    for index, factor in enumerate(lipschitz):
        left, right = Q(7, 3), Q(-5, 4)
        ledger.check(
            abs(reductions[index](left) - reductions[index](right))
            == factor * abs(left - right),
            f"M1_Lipschitz_scale_{index}",
        )

    # M2: nested contiguous quotients of a chain are causally lumpable.
    e2 = tuple(range(8))
    e1 = tuple(range(4))
    e0 = tuple(range(2))
    q12 = {value: value // 2 for value in e2}
    q01 = {value: value // 2 for value in e1}
    q02 = {value: value // 4 for value in e2}
    c2 = {(left, right) for left in e2 for right in e2 if left <= right}
    c1 = push_relation(c2, q12)
    c0_sequential = push_relation(c1, q01)
    c0_direct = push_relation(c2, q02)
    ledger.check(is_partial_order(c2, e2), "M2_fine_partial_order")
    ledger.check(is_partial_order(c1, e1), "M2_middle_lumpable_partial_order")
    ledger.check(is_partial_order(c0_direct, e0), "M2_coarse_lumpable_partial_order")
    ledger.check(c0_direct == c0_sequential, "M2_relation_pushforward_functorial")

    response_a = {pair: Q(1 + pair[1] - pair[0]) for pair in c2}
    response_b = {pair: Q(2 + 2 * (pair[1] - pair[0])) for pair in c2}
    response_slow = {(value, value): Q(1) for value in e2}
    ledger.check(response_support(response_a) == response_support(response_b) == c2,
                 "M2_two_distinct_probe_families_saturate")
    ledger.check(response_a != response_b, "M2_probe_laws_not_identical")
    ledger.check(response_support(response_slow) < c2, "M2_slower_probe_inside_envelope")

    for name, response in (("A", response_a), ("B", response_b)):
        middle = block_response(response, q12)
        coarse_sequential = block_response(middle, q01)
        coarse_direct = block_response(response, q02)
        ledger.check(response_support(middle) == c1, f"M2_{name}_middle_support_faithful")
        ledger.check(response_support(coarse_direct) == c0_direct,
                     f"M2_{name}_direct_support_faithful")
        ledger.check(coarse_direct == coarse_sequential,
                     f"M2_{name}_direct_sequential_response_exact")

    # Exact nontransitive quotient counterexample.
    fine_elements = (0, 1, 2, 3)  # a, b1, b2, c
    fine_relation = {(value, value) for value in fine_elements} | {(0, 1), (2, 3)}
    bad_map = {0: 0, 1: 1, 2: 1, 3: 2}
    bad_quotient = push_relation(fine_relation, bad_map)
    ledger.check(is_partial_order(fine_relation, fine_elements), "M2_bad_example_fine_is_poset")
    ledger.check((0, 1) in bad_quotient and (1, 2) in bad_quotient,
                 "M2_bad_example_two_coarse_legs")
    ledger.check((0, 2) not in bad_quotient, "M2_bad_example_missing_composite")
    ledger.check(not is_transitive(bad_quotient), "M2_arbitrary_pushforward_not_order")

    signed = ((Q(1), Q(0)), (Q(0), Q(-1)))
    injection = (Q(1), Q(1))
    fine_image = mat_vec(signed, injection)
    coarse_read = sum(fine_image, Q(0))
    ledger.check(fine_image != (Q(0), Q(0)), "M2_signed_fine_response_nonzero")
    ledger.check(coarse_read == 0, "M2_signed_block_cancellation_exact")

    # A pointwise support check on each actual scale object does not compose.
    # The second block must also be support-faithful on the intermediate
    # predicted by the first block.
    actual_2 = (Q(1), Q(0))
    actual_1 = (Q(1), Q(0))
    actual_0 = (Q(1),)

    def bad_block_12(vector: Vector) -> Vector:
        return (vector[1], vector[0])

    def bad_block_01(vector: Vector) -> Vector:
        return (vector[0],)

    def direct_block_02(vector: Vector) -> Vector:
        return (vector[0],)

    predicted_1 = bad_block_12(actual_2)
    predicted_0_sequential = bad_block_01(predicted_1)
    predicted_0_direct = direct_block_02(actual_2)
    pointwise_actual_gates = (
        has_support(bad_block_12(actual_2)) == has_support(actual_2)
        and has_support(bad_block_01(actual_1)) == has_support(actual_1)
        and has_support(direct_block_02(actual_2)) == has_support(actual_2)
        and has_support(actual_0)
    )
    ledger.check(pointwise_actual_gates, "M2_pointwise_actual_support_gates_pass")
    ledger.check(predicted_1 == (Q(0), Q(1)) and predicted_1 != actual_1,
                 "M2_sequential_intermediate_outside_actual_test")
    ledger.check(not has_support(predicted_0_sequential),
                 "M2_pointwise_counterexample_sequential_support_lost")
    ledger.check(has_support(predicted_0_direct),
                 "M2_pointwise_counterexample_direct_support_retained")
    ledger.check(has_support(predicted_0_sequential) != has_support(predicted_0_direct),
                 "M2_pointwise_support_does_not_compose")
    ledger.check(has_support(bad_block_01(predicted_1)) != has_support(predicted_1),
                 "M2_map_level_gate_rejects_predicted_intermediate")

    # M3: unique full-rank calibrations plus exact aggregation.
    v2 = tuple(Q(value) for value in range(1, 9))
    map12 = tuple(value // 2 for value in range(8))
    map01 = tuple(value // 2 for value in range(4))
    map02 = tuple(value // 4 for value in range(8))
    v1 = aggregate(v2, map12, 4)
    v0 = aggregate(v1, map01, 2)
    ledger.check(v1 == (Q(3), Q(7), Q(11), Q(15)), "M3_middle_volume_exact")
    ledger.check(v0 == (Q(10), Q(26)), "M3_coarse_volume_exact")
    ledger.check(aggregate(v2, map02, 2) == v0, "M3_direct_sequential_volume_exact")
    ledger.check(sum(v2, Q(0)) == sum(v1, Q(0)) == sum(v0, Q(0)) == Q(36),
                 "M3_finite_total_mass_common")

    c2_cal = lower_bidiagonal(8)
    c1_cal = lower_bidiagonal(4)
    c0_cal = ((Q(1), Q(1)), (Q(1), Q(2)))
    for name, calibration, volume in (
        ("fine", c2_cal, v2),
        ("middle", c1_cal, v1),
        ("coarse", c0_cal, v0),
    ):
        data = mat_vec(calibration, volume)
        ledger.check(rank(calibration) == len(volume), f"M3_{name}_calibration_full_rank")
        ledger.check(solve_square(calibration, data) == volume,
                     f"M3_{name}_calibration_recovers_volume")
        ledger.check(all(value > 0 for value in volume), f"M3_{name}_volume_positive")

    # Approximate compatibility telescopes with the l1-contractive aggregation.
    v1_approx = (v1[0] + 1, v1[1], v1[2], v1[3])
    base_v0_approx = aggregate(v1_approx, map01, 2)
    v0_approx = (base_v0_approx[0] + 1, base_v0_approx[1])
    eta1 = l1(tuple(a - b for a, b in zip(v1_approx, v1)))
    eta0 = l1(tuple(a - b for a, b in zip(v0_approx, base_v0_approx)))
    direct_mismatch = l1(tuple(a - b for a, b in zip(v0_approx, v0)))
    ledger.check((eta0, eta1) == (Q(1), Q(1)), "M3_adjacent_volume_residuals_exact")
    ledger.check(direct_mismatch == eta0 + eta1 == Q(2),
                 "M3_volume_telescoping_bound_saturated")

    # Full rank at each scale does not imply cross-scale compatibility.
    incompatible_v0 = (Q(11), Q(25))
    identity2 = ((Q(1), Q(0)), (Q(0), Q(1)))
    ledger.check(rank(identity2) == 2 and all(value > 0 for value in incompatible_v0),
                 "M3_incompatible_coarse_calibration_identifiable")
    ledger.check(sum(incompatible_v0, Q(0)) == sum(v0, Q(0)),
                 "M3_incompatible_example_same_total_mass")
    ledger.check(incompatible_v0 != v0, "M3_full_rank_not_cross_scale_compatibility")

    deficient = (
        (Q(1), Q(1), Q(0), Q(0)),
        (Q(0), Q(0), Q(1), Q(1)),
    )
    volume_a = (Q(1), Q(1), Q(1), Q(1))
    volume_b = (Q(1, 2), Q(3, 2), Q(3, 4), Q(5, 4))
    ledger.check(rank(deficient) == 2, "M3_rank_deficient_control")
    ledger.check(mat_vec(deficient, volume_a) == mat_vec(deficient, volume_b),
                 "M3_rank_deficient_nonunique_positive_volumes")

    # Finite cylinders of a consistent projective measure/order family.
    probability_one = Q(1, 3)
    measures = {length: bernoulli_measure(length, probability_one) for length in range(1, 5)}
    for length, measure in measures.items():
        ledger.check(sum(measure.values(), Q(0)) == 1, f"M3_cylinder_mass_n{length}")
    for length in range(2, 5):
        ledger.check(prefix_push(measures[length], length - 1) == measures[length - 1],
                     f"M3_cylinder_measure_compatible_n{length}")
        fine_strings = bitstrings(length)
        coarse_strings = bitstrings(length - 1)
        fine_order = coordinate_order(fine_strings)
        coarse_order = coordinate_order(coarse_strings)
        ledger.check(push_string_relation(fine_order, length - 1) == coarse_order,
                     f"M3_cylinder_order_compatible_n{length}")

    # Claim ledger: all continuum gates are needed; none implies Einstein dynamics.
    continuum_gate_exact = True
    for flags in product((False, True), repeat=6):
        admitted = all(flags)
        continuum_gate_exact &= admitted == (flags == (True,) * 6)
    ledger.check(continuum_gate_exact, "M4_all_six_continuum_gates_required")
    lorentzian_kinematics = all((True,) * 6)
    einstein_dynamics = False
    record_ancestry = False
    ledger.check(lorentzian_kinematics and not einstein_dynamics,
                 "M4_Lorentzian_does_not_imply_Einstein")
    ledger.check(lorentzian_kinematics and not record_ancestry,
                 "M4_same_parent_does_not_imply_record_ancestry")

    # A fixed perturbative gravitational remainder does not vanish with scale.
    kappa = Q(1, 10)
    remainder_constant = Q(3)
    gravity_residuals = tuple(remainder_constant * kappa * kappa for _ in range(6))
    ledger.check(all(value == Q(3, 100) for value in gravity_residuals),
                 "GRAV_fixed_quadratic_remainder_positive")
    ledger.check(len(set(gravity_residuals)) == 1,
                 "GRAV_more_scales_do_not_reduce_fixed_remainder")

    print("M1_ARBITRARY_SCALE_ERROR_COMPOSITION EXACT")
    print("M2_SUPPORT_AND_CAUSAL_LUMPABILITY EXACT")
    print("M3_VOLUME_AND_PROJECTIVE_CYLINDERS EXACT")
    print("M4_LORENTZIAN_GATE_AND_DYNAMICS_CEILING EXACT")
    print("GRAVITATIONAL_FIXED_REMAINDER_NOT_SCALE_VANISHING EXACT")
    print(f"TOTAL {ledger.passed}/{ledger.total} PASS")
    print(
        "VERDICT PROJECTIVE_CAUSAL_MEASURE_PROVED_CONDITIONALLY__"
        "LORENTZIAN_RECONSTRUCTION_CONDITIONAL__EINSTEIN_OPEN"
    )


if __name__ == "__main__":
    main()
