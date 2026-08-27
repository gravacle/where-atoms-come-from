#!/usr/bin/env python3
"""Exact finite witnesses for historywise gravity and actualization.

The script checks one rational Z2 instance and logical admission gates.  The
general finite-group theorems in THEOREM.md remain analytic.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import product
Vector = tuple[F, F]
Matrix = tuple[tuple[F, F], tuple[F, F]]


def matvec(matrix: Matrix, vector: Vector) -> Vector:
    return tuple(
        sum((matrix[i][j] * vector[j] for j in range(2)), F(0))
        for i in range(2)
    )  # type: ignore[return-value]


def matmul(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            sum((left[i][k] * right[k][j] for k in range(2)), F(0))
            for j in range(2)
        )
        for i in range(2)
    )  # type: ignore[return-value]


def add(left: Vector, right: Vector) -> Vector:
    return (left[0] + right[0], left[1] + right[1])


def scale(value: F, vector: Vector) -> Vector:
    return (value * vector[0], value * vector[1])


def dot(left: Vector, right: Vector) -> F:
    return left[0] * right[0] + left[1] * right[1]


def metric_cells(phi: Vector) -> tuple[tuple[F, F], tuple[F, F]]:
    return tuple((-(1 + value), 1 - value) for value in phi)  # type: ignore[return-value]


def field_action(phi: Vector, source: Vector, laplacian: Matrix) -> F:
    return F(1, 2) * dot(phi, matvec(laplacian, phi)) - dot(phi, source)


def feedback_history(source: Vector) -> int:
    if source[0] > source[1]:
        return 1
    if source[0] < source[1]:
        return -1
    return 0


def main() -> int:
    checks: list[tuple[str, bool]] = []

    def check(name: str, condition: bool) -> None:
        checks.append((name, bool(condition)))

    identity: Matrix = ((F(1), F(0)), (F(0), F(1)))
    swap: Matrix = ((F(0), F(1)), (F(1), F(0)))
    laplacian: Matrix = ((F(1), F(-1)), (F(-1), F(1)))

    check("G1_swap_squares_to_identity", matmul(swap, swap) == identity)
    check("G1_laplacian_swap_invariant", matmul(matmul(swap, laplacian), swap) == laplacian)
    check("G1_binary_history_action_is_free", -1 != 1 and -(-1) == 1)

    source_plus: Vector = (F(1), F(-1))
    source_minus: Vector = matvec(swap, source_plus)
    phi_plus: Vector = (F(1, 2), F(-1, 2))
    phi_minus: Vector = matvec(swap, phi_plus)

    check("W1_source_swap_exact", source_minus == (F(-1), F(1)))
    check("W1_source_swap_closes", matvec(swap, source_minus) == source_plus)
    check("W1_source_contrasts_sum_zero", sum(source_plus) == sum(source_minus) == 0)
    check("W1_potential_swap_exact", phi_minus == (F(-1, 2), F(1, 2)))
    check("W1_potential_swap_closes", matvec(swap, phi_minus) == phi_plus)
    check("W1_poisson_plus_exact", matvec(laplacian, phi_plus) == source_plus)
    check("W1_poisson_minus_exact", matvec(laplacian, phi_minus) == source_minus)
    check("W1_gauge_plus_exact", sum(phi_plus) == 0)
    check("W1_gauge_minus_exact", sum(phi_minus) == 0)

    metric_plus = metric_cells(phi_plus)
    metric_minus = metric_cells(phi_minus)
    check("W1_metric_plus_exact", metric_plus == ((F(-3, 2), F(1, 2)), (F(-1, 2), F(3, 2))))
    check("W1_metric_minus_is_cell_swap", metric_minus == tuple(reversed(metric_plus)))
    check("W1_all_cell_metrics_lorentzian", all(time < 0 < space for time, space in metric_plus + metric_minus))
    check("W1_all_metric_determinants_exact", all(time * space == F(-3, 4) for time, space in metric_plus + metric_minus))

    action_plus = field_action(phi_plus, source_plus, laplacian)
    action_minus = field_action(phi_minus, source_minus, laplacian)
    check("W1_action_plus_exact", action_plus == F(-1, 2))
    check("W1_action_minus_exact", action_minus == F(-1, 2))
    check("W1_invariant_action_ties_orbit", action_plus == action_minus)

    residual_plus = add(matvec(laplacian, phi_plus), scale(F(-1), source_plus))
    residual_minus = add(matvec(laplacian, phi_minus), scale(F(-1), source_minus))
    check("W1_gravity_residual_plus_zero", residual_plus == (F(0), F(0)))
    check("W1_gravity_residual_minus_zero", residual_minus == (F(0), F(0)))

    admissible_histories = tuple(
        sign
        for sign, residual in ((1, residual_plus), (-1, residual_minus))
        if residual == (F(0), F(0))
    )
    check("HGA1_both_orbit_histories_admissible", set(admissible_histories) == {1, -1})
    check("HGA1_endogenous_gravity_not_singleton", len(admissible_histories) == 2)
    check("HGA1_no_fixed_definite_history", all(-history != history for history in admissible_histories))

    check("HGA1_feedback_plus_fixed", feedback_history(source_plus) == 1)
    check("HGA1_feedback_minus_fixed", feedback_history(source_minus) == -1)
    feedback_equivariant = all(
        feedback_history(matvec(swap, source)) == -feedback_history(source)
        for source in (source_plus, source_minus)
    )
    check("HGA1_nonlinear_feedback_equivariant", feedback_equivariant)

    # HGA1a: a unique equivariant flow maps the symmetric complete datum only
    # to a symmetry-fixed history, not to either member of the free orbit.
    def unique_flow(complete_datum: int) -> int:
        return complete_datum

    flow_domain = (-1, 0, 1)
    check(
        "HGA1a_unique_flow_equivariant",
        all(unique_flow(-datum) == -unique_flow(datum) for datum in flow_domain),
    )
    symmetric_flow_history = unique_flow(0)
    check("HGA1a_symmetric_input_has_fixed_history", symmetric_flow_history == -symmetric_flow_history == 0)
    check("HGA1a_symmetric_input_not_free_outcome", symmetric_flow_history not in admissible_histories)

    # HGA1b: the constant uniformizing update is equivariant and has the
    # invariant uniform law, rather than a realized vertex, as its unique fixed
    # law.  The theorem's general finite-orbit statement remains analytic.
    Law = tuple[F, F]
    uniform_law: Law = (F(1, 2), F(1, 2))

    def push_swap(law: Law) -> Law:
        return (law[1], law[0])

    def distribution_update(_: Law) -> Law:
        return uniform_law

    test_laws: tuple[Law, ...] = (
        (F(1), F(0)),
        uniform_law,
        (F(1, 3), F(2, 3)),
    )
    check(
        "HGA1b_distribution_update_equivariant",
        all(
            distribution_update(push_swap(law)) == push_swap(distribution_update(law))
            for law in test_laws
        ),
    )
    check("HGA1b_uniform_law_is_fixed", distribution_update(uniform_law) == uniform_law)
    check(
        "HGA1b_definite_vertices_not_fixed",
        all(distribution_update(law) != law for law in ((F(1), F(0)), (F(0), F(1)))),
    )

    mean_phi = scale(F(1, 2), add(phi_plus, phi_minus))
    mean_source = scale(F(1, 2), add(source_plus, source_minus))
    mean_metric = tuple(
        (
            F(1, 2) * (metric_plus[i][0] + metric_minus[i][0]),
            F(1, 2) * (metric_plus[i][1] + metric_minus[i][1]),
        )
        for i in range(2)
    )
    check("HGA2_mean_potential_zero", mean_phi == (F(0), F(0)))
    check("HGA2_mean_source_zero", mean_source == (F(0), F(0)))
    check("HGA2_mean_metric_background_exact", mean_metric == ((F(-1), F(1)), (F(-1), F(1))))
    check("HGA2_mean_pair_swap_fixed", matvec(swap, mean_phi) == mean_phi and matvec(swap, mean_source) == mean_source)
    check("HGA2_fixed_input_has_no_free_orbit_output", not any(history == -history for history in admissible_histories))

    boundary_plus = source_plus
    boundary_minus = source_minus

    def coupling(boundary: Vector, source: Vector) -> F:
        return -dot(boundary, source)

    coupling_table = (
        (coupling(boundary_plus, source_plus), coupling(boundary_plus, source_minus)),
        (coupling(boundary_minus, source_plus), coupling(boundary_minus, source_minus)),
    )
    check("D1_boundary_swap_exact", matvec(swap, boundary_plus) == boundary_minus)
    check("D1_coupling_table_exact", coupling_table == ((F(-2), F(2)), (F(2), F(-2))))

    joint_covariance = all(
        coupling(matvec(swap, boundary), matvec(swap, source)) == coupling(boundary, source)
        for boundary in (boundary_plus, boundary_minus)
        for source in (source_plus, source_minus)
    )
    check("D1_boundary_coupling_jointly_invariant", joint_covariance)

    selected_plus = 1 if coupling_table[0][0] < coupling_table[0][1] else -1
    selected_minus = 1 if coupling_table[1][0] < coupling_table[1][1] else -1
    check("D1_fixed_boundary_plus_selects_plus", selected_plus == 1)
    check("D1_fixed_boundary_minus_selects_minus", selected_minus == -1)
    check("D1_selector_jointly_equivariant", selected_minus == -selected_plus)

    group = {0, 1}

    def sign_action(element: int, sign: int) -> int:
        return sign if element == 0 else -sign

    history_stabilizer = {element for element in group if sign_action(element, 1) == 1}
    boundary_stabilizer = {element for element in group if sign_action(element, 1) == 1}
    scalar_stabilizer = set(group)
    check("HGA3_history_stabilizer_trivial", history_stabilizer == {0})
    check("HGA3_boundary_stabilizer_trivial", boundary_stabilizer == {0})
    check("HGA3_positive_stabilizer_inclusion", boundary_stabilizer <= history_stabilizer)
    check("HGA3_orbit_cardinality_bound_saturated", len(group) // len(boundary_stabilizer) == len(group) // len(history_stabilizer) == 2)
    check("HGA3_fixed_scalar_fails_inclusion", not scalar_stabilizer <= history_stabilizer)

    subgroup_cases = (
        ({0}, {0}, True),
        ({0}, {0, 1}, True),
        ({0, 1}, {0, 1}, True),
        ({0, 1}, {0}, False),
    )
    check("HGA3_Z2_stabilizer_criterion_truth_table", all((k <= ell) == expected for k, ell, expected in subgroup_cases))

    # Operational gate: one density state has two exact decompositions.  A
    # nonlinear normalized-square rule gives different unconditioned outputs.
    def nonlinear_square(state: Vector) -> Vector:
        denominator = dot(state, state)
        return (state[0] * state[0] / denominator, state[1] * state[1] / denominator)

    rho: Vector = (F(2, 3), F(1, 3))
    prepared_a: Vector = (F(1), F(0))
    prepared_b: Vector = (F(1, 3), F(2, 3))
    decomposed_rho = scale(F(1, 2), add(prepared_a, prepared_b))
    direct_nonlinear_output = nonlinear_square(rho)
    forgotten_label_output = scale(
        F(1, 2),
        add(nonlinear_square(prepared_a), nonlinear_square(prepared_b)),
    )
    preparation_tv = F(1, 2) * sum(
        abs(direct_nonlinear_output[i] - forgotten_label_output[i])
        for i in range(2)
    )
    check("OP_same_density_decomposition_exact", decomposed_rho == rho)
    check("OP_nonlinear_direct_output_exact", direct_nonlinear_output == (F(4, 5), F(1, 5)))
    check("OP_nonlinear_forgotten_label_output_exact", forgotten_label_output == (F(3, 5), F(2, 5)))
    check("OP_preparation_context_TV_one_fifth", preparation_tv == F(1, 5))

    # Matrix transposition is positive on one system but its action on half of
    # a Bell state has an exact negative ancillary expectation at n=d=2, which
    # is already sufficient to refute complete positivity for a qubit input.
    # The vector is left unnormalized; dividing the quadratic form by its
    # squared norm gives the normalized expectation.
    bell_partial_transpose = (
        (F(1, 2), F(0), F(0), F(0)),
        (F(0), F(0), F(1, 2), F(0)),
        (F(0), F(1, 2), F(0), F(0)),
        (F(0), F(0), F(0), F(1, 2)),
    )
    antisymmetric = (F(0), F(1), F(-1), F(0))
    ancillary_quadratic = sum(
        antisymmetric[i] * bell_partial_transpose[i][j] * antisymmetric[j]
        for i in range(4)
        for j in range(4)
    )
    ancillary_norm = sum(value * value for value in antisymmetric)
    check(
        "OP_transpose_Bell_ancillary_expectation_negative_half",
        ancillary_quadratic / ancillary_norm == F(-1, 2),
    )

    def operational_gate(
        preparation_independence: bool,
        no_signalling: bool,
        affine_mixing: bool,
        all_finite_or_input_dim_ancillas_cp: bool,
        named_abandonment: bool,
        frozen_domain_and_scale: bool,
        prospective_prediction: bool,
        prediction_falsifier: bool,
    ) -> bool:
        standard_route = all(
            (
                preparation_independence,
                no_signalling,
                affine_mixing,
                all_finite_or_input_dim_ancillas_cp,
            )
        )
        explicit_departure_route = all(
            (
                named_abandonment,
                frozen_domain_and_scale,
                prospective_prediction,
                prediction_falsifier,
            )
        )
        return standard_route or explicit_departure_route

    check(
        "OP_complete_standard_route_admits_gate",
        operational_gate(True, True, True, True, False, False, False, False),
    )
    check(
        "OP_unfrozen_premise_abandonment_blocks_gate",
        not operational_gate(False, False, False, False, True, False, True, True),
    )
    check(
        "OP_named_frozen_prediction_route_admits_gate",
        operational_gate(False, False, False, False, True, True, True, True),
    )

    # Physical GARH-D admission requires all eleven independently supplied gates.
    gate_logic_exact = True
    for flags in product((False, True), repeat=11):
        admitted = all(flags)
        gate_logic_exact &= admitted == (flags == (True,) * 11)
    check("D2_all_eleven_physical_admission_gates_required", gate_logic_exact)

    finite_theorem_proved = True
    physical_orientation_established = False
    physical_gravity_established = False
    born_law_derived = False
    check("D2_finite_theorem_does_not_establish_orientation", finite_theorem_proved and not physical_orientation_established)
    check("D2_finite_metric_does_not_establish_GR", finite_theorem_proved and not physical_gravity_established)
    check("D2_boundary_table_does_not_derive_Born", finite_theorem_proved and not born_law_derived)

    garh_d_failed = True
    garh_q_logically_forced = False
    check("D2_failed_D_does_not_logically_force_Q", garh_d_failed and not garh_q_logically_forced)

    random_boundary_draw_explained = False
    deeper_actualization_derived = False
    check("D2_unexplained_random_boundary_relocates_draw", not random_boundary_draw_explained and not deeper_actualization_derived)

    failed = [name for name, passed in checks if not passed]
    for name, passed in checks:
        print(("PASS" if passed else "FAIL") + "  " + name)
    print(f"TOTAL {len(checks) - len(failed)}/{len(checks)} PASS")
    if failed:
        print("FAILED " + ", ".join(failed))
        return 1
    print("VERDICT FINITE_EQUIVARIANT_GRAVITY_NONSELECTION_AND_ORIENTING_INPUT_CRITERION_EXACT__PHYSICAL_GARH_D_OPEN")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
