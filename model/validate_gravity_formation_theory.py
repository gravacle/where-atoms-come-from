#!/usr/bin/env python3
"""Bounded custody, zero-input, immutability, and ceiling gate for GFT."""

from __future__ import annotations

import inspect
from pathlib import Path
import shutil
import tempfile

import gravity_formation_theory as gft
from project_model import URM


def _type_error(callable_object) -> bool:
    try:
        callable_object()
    except TypeError:
        return True
    return False


def _refused(callable_object) -> bool:
    try:
        callable_object()
    except gft.GravityFormationTheoryRefusal:
        return True
    return False


def _mutate(mapping, key, value) -> None:
    mapping[key] = value


def main() -> int:
    assert tuple(inspect.signature(gft.gravity_formation_theory).parameters) == ()
    assert tuple(inspect.signature(gft.gravity_formation_theory_certificate).parameters) == ()
    assert _type_error(lambda: gft.gravity_formation_theory({}))
    assert _type_error(lambda: gft.gravity_formation_theory_certificate(data={}))

    theory = gft.gravity_formation_theory()
    certificate = theory.certificate()
    assert theory.closure_sha256 == gft._CLOSURE_PAIR[0][1]
    assert theory.audit_sha256 == gft._CLOSURE_PAIR[1][1]
    assert certificate["schema"] == gft.SCHEMA
    assert certificate["claim_class"] == gft.CLAIM_CLASS
    assert (
        certificate["rgrl"]["status"]
        == "ADOPTED_WORKING_RECORD_GEOMETRY_REALIZATION_POSTULATE"
    )
    assert (
        certificate["scientific_status"]["nature_obeys_RGRL"]
        == "NOT_ESTABLISHED_BY_THIS_CERTIFICATE"
    )
    assert certificate["rgrl"]["RGRL_C_is_onshell_force_law"] is False
    clarification = certificate["adopted_response_clarification"]
    assert clarification["status"] == "ADOPTED_PROSPECTIVE_PROGRAM_POLICY"
    assert clarification["open_type_join"] == "GI21_NOT_SUPPLIED_BY_RGRL_C"
    assert (
        clarification["fully_matched_response"]
        == "ZERO_ON_V002_WELL_POSED_RETARDED_QUOTIENT_WITH_NO_UNRESOLVED_ZERO_MODE_"
        "WHEN_DRESSED_SOURCE_REMAINDER_AND_PHYSICAL_DATA_DERIVATIVES_VANISH"
    )
    assert clarification["fully_matched_response_premises"] == (
        "V002_WELL_POSED_RETARDED_GAUGE_BOUNDARY_QUOTIENT",
        "NO_UNRESOLVED_ZERO_MODE",
    )
    assert (
        clarification["fully_matched_H_R"]
        == "MAY_VANISH_WITHOUT_REDUCING_OFFSHELL_RANK"
    )
    assert (
        certificate["program_authorizations"]["exact_conditional_working_theory_closure"]
        is True
    )
    assert certificate["program_authorizations"]["empirical_RGRL_confirmation"] is False
    assert certificate["program_authorizations"]["gravity_solver"] is False
    assert certificate["program_authorizations"]["caller_data_admitted"] is False
    assert len(certificate["kernel_reduction"]["S4_response_form_factors"]) == 3
    assert len(certificate["kernel_reduction"]["compatible_O3_reduction"]) == 2
    assert certificate["kernel_reduction"]["nonzero_spatial_momentum_classification"] == "OPEN"
    assert certificate["kernel_reduction"]["GI21_compatibility_type_join"] == "OPEN"
    assert (
        certificate["kernel_reduction"]["onshell_H_R_nonzero"]
        == "NOT_IMPLIED_AND_MAY_VANISH_IN_FULLY_MATCHED_LANE"
    )
    assert (
        certificate["exact_results"]["record_front_lorentz_cone_refinement"]
        == "EXACT_CONDITIONAL_ON_AFR_AND_MICROSCOPIC_EQUAL_DEPTH_EQUAL_QJ_NORM_NULL_STEPS"
    )
    assert certificate["exact_results"]["cone_refinement_is_full_RGRL_or_gravity_dynamics"] is False
    assert certificate["exact_results"]["refined_cone_kind"] == "LOCAL_MATHEMATICAL_3PLUS1_LORENTZ_CONE"
    assert certificate["exact_results"]["six_mode_cone_tangent_scope"] == "J0_S4_SYMMETRIC_POINT"
    assert (
        certificate["exact_results"]["q4_common_child_incidence_identity"]
        == "B_DAGGER_B_EQUALS_4I_PLUS_A_EXACT"
    )
    assert (
        certificate["exact_results"]["q4_affine_refinement"]
        == "EXACT_MATHEMATICAL_ATLAS_NOT_PHYSICAL_MANIFOLD"
    )
    assert (
        certificate["exact_results"]["q4_acoustic_cone"]
        == "EXACT_ONLY_FOR_SEPARATELY_SUPPLIED_MASSLESS_ACTION"
    )
    assert (
        certificate["exact_results"]["q4_diamond_u1_inheritance"]
        == "FINITE_LOCAL_EDGE_BINDING_PROGRAMMABLE__GLOBAL_U1_CONDITIONAL_ON_"
        "SUPPLIED_REGULAR_COMPLETION_AND_LEADING_ORDER"
    )
    assert (
        certificate["exact_results"]["q4_programmed_floquet_detuning"]
        == "EXACT_CHILD_ONLY_PHASE_AND_UNIFORM_PARENT_CHILD_QUASIENERGY_SEPARATION_"
        "FOR_SUPPLIED_REPEATABLE_DUAL_FLIP_FREE_SCHEDULE"
    )
    assert certificate["exact_results"]["q4_static_source_off_child_parent_stagger"] == "ABSENT"
    assert certificate["exact_results"]["q4_ice_local_diagonal_module"] == "A1_PLUS_E_PLUS_T2_EXACT"
    assert certificate["exact_results"]["q4_ice_symmetric_fisher_T2_first_derivative"] == "EXACTLY_ZERO_BY_COMPLEMENT_PARITY"
    assert certificate["exact_results"]["q4_ice_fisher_covariance"] == "F_EQUALS_4_DIAG_W_MINUS_M_M_TRANSPOSE_EXACT"
    assert (
        certificate["exact_results"]["q4_ice_odd_T2_fisher_second_jet"]
        == "S4_ISOMORPHISM_SYM2_T2_TO_A1_E_T2__NOT_O3_EQUIVARIANT"
    )
    assert (
        certificate["exact_results"]["q4_ice_complement_preserving_fisher_tangent"]
        == "E_ONLY__A1_PLUS_E_WITH_SEPARATE_SCALAR"
    )
    assert (
        certificate["exact_results"][
            "q4_ice_generic_broken_background_covariance_map_with_separate_scalar_rank"
        ]
        == 6
    )
    assert (
        certificate["exact_results"]["q4_ice_broken_background_T2_origin"]
        == "VECTOR_MEAN_DYAD_NOT_INDEPENDENT_TENSOR"
    )
    assert certificate["exact_results"]["q4_ice_unlabelled_complement_T2"] == "CANCELS_EXACTLY"
    assert (
        certificate["exact_results"]["q4_ice_complement_symmetrized_conditional_T2"]
        == "RETENTION_REQUIRES_CONTROL_SIGN_RECORD"
    )
    assert certificate["exact_results"]["q4_ice_broken_background_preparation_and_stabilization"] == "NOT_DERIVED"
    assert (
        certificate["exact_results"]["q4_ice_order6_interaction"]
        == "COMPACT_HARDCORE_RING_ALREADY_NON_GAUSSIAN"
    )
    assert (
        certificate["exact_results"]["q4_ice_order8_scope"]
        == "SUPPLIED_FINITE_SIMPLE_Z4_BIPARTITE_GIRTH_GE6_SUPPORT__DSTAR2_ER0_"
        "SYMMETRIC_DETUNING__FIXED_FESHBACH_CONVENTION"
    )
    assert (
        certificate["exact_results"]["q4_ice_order8_endpoint_topologies"]
        == "SCALAR_DIAGONAL_PLUS_DRESSED_HEXAGON_PLUS_NEW_ALTERNATING_OCTAGON"
    )
    assert (
        certificate["exact_results"]["q4_ice_order8_octagon_coefficient"]
        == "J8_EQUALS_429_H8_OVER_16_UD7"
    )
    assert certificate["exact_results"]["q4_ice_order8_diagonal_potential"] == "V8_EQUALS_ZERO"
    microscopic = certificate["microscopic_parent_boundary"]
    assert (
        microscopic["scalar_carrier_transfer_on_supplied_saturated_q4_support"]
        == "EXACT_FROM_UNCHANGED_F3_ONE_CARRIER_RESTRICTION"
    )
    assert microscopic["full_hopping_and_d2_ice_same_n_coexistence"] == "EXACTLY_INCOMPATIBLE"
    assert (
        microscopic["finite_programmed_site_edge_solder"]
        == "EXACT_REVERSIBLE_FIXED_ORTHOGONAL_PROGRAM_WITH_SUPPLIED_PHYSICAL_ANTECEDENTS"
    )
    assert microscopic["finite_programmed_solder_is_autonomous_support_selection"] is False
    assert microscopic["interpair_retarded_kernel"] == "EXACTLY_ZERO_UNDER_INHERITED_DYNAMICS"
    assert microscopic["programmed_floquet_detuning_is_autonomous_phase"] is False
    assert (
        microscopic["physical_q4_link_pair_response"]
        == "EXACT_FINITE_LOCAL_AND_SHARED_LINK_RESPONSE_WITH_OPERATOR_SPREADING"
    )
    assert (
        microscopic["physical_q4_link_walsh_pairs"]
        == "EXACT_OPERATOR_REALIZATION_ON_FPMH_QUALIFIED_FINITE_PROGRAMMED_LINK_FACTORS"
    )
    assert (
        microscopic["q4_pair_operator_realization"]
        == "EXACT_FINITE_WALSH_ALGEBRA__FULL_PMMDC_AND_METRIC_SOLDER_OPEN"
    )
    assert microscopic["physical_q4_link_pairs_are_automatically_PMMDC_records"] is False
    assert microscopic["ice_ring_response_pole"] == "FINITE_2J6_NOT_MASSLESS_TENSOR"
    assert microscopic["direct_gaussian_composite_tensor_route"] == "EXACT_CONDITIONAL_NO_GO"
    assert (
        microscopic["ice_T2_fisher_solder_boundary"]
        == "SECOND_JET_OR_VECTOR_BACKGROUND_ONLY__INDEPENDENT_LINEAR_TENSOR_SOLDER_OPEN"
    )
    assert (
        microscopic["inherited_order8_operator_boundary"]
        == "EXACT_J8_AND_V8_ZERO__HEXAGON_DRESSING_COEFFICIENTS_TYPED_NOT_ALL_REDUCED"
    )
    assert microscopic["new_interaction_or_second_field_adopted"] is False
    assert certificate["kernel_reduction"]["gaussian_maxwell_one_link_pole"] == "SPIN_1_PHOTON_CONDITIONAL_ON_MAXWELL_IR"
    assert certificate["kernel_reduction"]["gaussian_maxwell_isolated_helicity2_pole"] == "ABSENT_IN_DIRECT_COMPOSITE_ROUTE"
    assert (
        certificate["kernel_reduction"]["ice_local_fisher_metric_full_linear_tensor_tangent"]
        == "ABSENT_WITH_COMPLEMENT_PRESERVED"
    )
    assert (
        certificate["kernel_reduction"]["ice_broken_background_full_rank_kind"]
        == "BACKGROUND_CONTROL_TO_COVARIANCE_MAP_NOT_SIX_MODE_FISHER_METRIC"
    )
    assert (
        certificate["kernel_reduction"]["inherited_TT_bare_vertex"]
        == "FINITE_RANGE_ANALYTIC__TT_PROJECTION_ALLOWED_NOT_EVALUATED"
    )
    assert certificate["kernel_reduction"]["inherited_TT_strict_single_insertion_new_pole"] == "NOT_ESTABLISHED"
    assert (
        certificate["kernel_reduction"]["inherited_TT_decisive_sequence"]
        == "CONNECTED_FOUR_POINT_THEN_AMPUTATION_THEN_CHANNEL_2PI_THEN_"
        "BETHE_SALPETER_OR_SPECTRAL_THEN_FINITE_VOLUME_HELICITY_WARD_"
        "COMMON_CONE_RESIDUE"
    )
    assert certificate["numerical_G"]["parameter_free_microscopic_record_derivation"] is False
    assert certificate["numerical_G"]["bounded_remainder_output"] == "IDENTIFIED_INTERVAL_OR_SET_NOT_AUTOMATIC_POINT_VALUE"
    assert certificate["numerical_G"]["calibrated_nonsingular_row_identifies"] == "P_EQUALS_G_TIMES_SOURCE_SCALE"
    assert certificate["numerical_G"]["synthetic_validation"] == "PASS_15_OF_15_WITH_ARBITRARY_NONEMPIRICAL_G"
    assert certificate["numerical_G"]["nist_bipm_analysis_observations"] == 8
    assert certificate["numerical_G"]["nist_bipm_G_only_jacobian_rank"] == 1
    assert certificate["numerical_G"]["nist_bipm_full_GC16_fit_ready"] is False
    assert certificate["numerical_G"]["nist_bipm_public_missing_field_count"] == 10
    assert certificate["numerical_G"]["nist_bipm_independent_G_crosscheck_performed"] is False
    assert certificate["SPAG"]["executed_by_this_surface"] is False
    assert certificate["SPAG"]["full_RGRL_confirmation"] is False
    assert certificate["SPAG"]["old_local_RGRL_C_force_column_labels"] == "RETIRED_FOR_FUTURE_CLAIMS"
    assert certificate["SPAG"]["lane_A"] == "COMPLETE_SOURCE_MATCHED_DISCOVERY_WITH_ZERO_PHYSICAL_PREDICTION"
    assert certificate["SPAG"]["lane_B"] == "INDEPENDENT_SOURCE_CALIBRATED_G_CROSS_CHECK"
    assert certificate["SPAG"]["public_data_substitute_executed"] is True
    assert certificate["SPAG"]["public_same_parent_eight_cell_support"] == "ABSENT"
    assert certificate["SPAG"]["public_beta_TM_identifiable"] is False
    assert certificate["SPAG"]["public_second_pass_frozen_query_count"] == 28
    assert certificate["SPAG"]["public_second_pass_new_lineage_roots"] == 0
    assert certificate["SPAG"]["public_second_pass_component_roots"] == 2
    assert certificate["SPAG"]["public_second_pass_result_is_exhaustive_world_search"] is False
    assert certificate["SPAG"]["panda_response_holdout_opened"] is False
    assert certificate["program_authorizations"]["old_SPAG_local_RGRL_C_force_labels_authorized"] is False
    assert (
        certificate["program_authorizations"][
            "autonomous_q4_support_or_full_pair_field_lift_derived"
        ]
        is False
    )
    assert certificate["program_authorizations"]["complement_preserving_local_fisher_T2_solder_derived"] is False
    assert certificate["program_authorizations"]["generic_broken_background_covariance_rank6_with_separate_scalar_derived"] is True
    assert certificate["program_authorizations"]["order8_inherited_loop_operator_boundary_derived"] is True
    assert certificate["program_authorizations"]["inherited_protected_tensor_pole_derived"] is False
    assert (
        certificate["scientific_status"]["inherited_non_gaussian_TT_kernel"]
        == "LEADING_H6_INTERACTION_IDENTIFIED__NORMALIZED_CONNECTED_TT_FOUR_POINT_OPEN"
    )
    assert certificate["scientific_status"]["normalized_connected_TT_four_point"] == "OPEN"
    assert certificate["scientific_status"]["nonperturbative_TT_pole_Ward_common_cone_residue"] == "OPEN"
    assert (
        certificate["scientific_status"]["local_ice_fisher_T2_solder"]
        == "CLOSED_NEGATIVE_FOR_LOCAL_DIAGONAL_COMPLEMENT_PRESERVING_FAMILIES"
    )
    assert certificate["nonpromotion"]["second_fisher_jet_promotes_linear_metric_tangent"] is False
    assert certificate["nonpromotion"]["broken_vector_background_rank_promotes_tensor_mode"] is False
    assert certificate["nonpromotion"]["control_sign_conditioning_promotes_endogenous_metric_solder"] is False
    assert certificate["nonpromotion"]["bare_vertex_analyticity_promotes_dressed_1PI_analyticity"] is False
    assert certificate["nonpromotion"]["strict_single_insertion_boundary_promotes_nonperturbative_no_pole"] is False
    assert certificate["custody"]["artifact_count"] == 50
    assert len(certificate["custody"]["core_source_pairs"]) == 6
    assert (
        certificate["custody"]["core_source_pairs"][4]["label"]
        == "f3_record_front_lorentz_cone_refinement"
    )
    advance_custody = certificate["custody"]["advance_source_pairs"]
    assert len(advance_custody) == 16
    assert tuple(item["label"] for item in advance_custody) == tuple(
        label for label, _, _ in gft._ADVANCE_SOURCE_PAIRS
    )
    assert advance_custody[3]["label"] == "f3_q4_carrier_lift_derivability_boundary"
    assert advance_custody[4]["label"] == "q4_pair_field_lift_derivability_boundary"
    assert advance_custody[5]["label"] == "f3_q4_finite_programmed_support_solder"
    assert advance_custody[6]["label"] == "f3_q4_programmed_floquet_detuning"
    assert advance_custody[7]["label"] == "f3_q4_authenticated_link_pair_response"
    assert advance_custody[8]["label"] == "f3_q4_ice_hybrid_tensor_response"
    assert advance_custody[9]["label"] == "f3_q4_maxwell_composite_pole_screen"
    assert advance_custody[10]["label"] == "f3_q4_inherited_tt_kernel_boundary"
    assert advance_custody[10]["source_sha256"] == (
        "78f0687c9f597c96d235224dad45c204d12d7e6c973c270e3305a901efc75b25"
    )
    assert advance_custody[10]["audit_sha256"] == (
        "53893c7198241f0f8f6aa766f3595fb75b83d208581833c32656b28d7c7f02b9"
    )
    assert advance_custody[11]["label"] == "f3_q4_ice_t2_fisher_solder_boundary"
    assert advance_custody[11]["source_sha256"] == (
        "be69f15d611827db9841bd932042604deb4f82a777ff9da28b80e4493cef7596"
    )
    assert advance_custody[11]["audit_sha256"] == (
        "32297dc0c4b0454c4a4be88d3763eb679b4ca89bb2385010ba8c2b77f2df47d2"
    )
    assert advance_custody[13]["label"] == "spag_public_data_second_pass"
    assert advance_custody[13]["audit_path"] == (
        "LANE_GRA_SPAG_PUBLIC_DATA_SECOND_PASS_V001/POST_INTEGRATION_CUSTODY_REAUDIT.md"
    )
    assert advance_custody[13]["audit_sha256"] == (
        "dcc170c2674fe0020f679e523ab40689be7efd32e05e71fe7600d9bbf4047e97"
    )
    assert advance_custody[15]["label"] == "nist_bipm_g_forward_readiness"
    clarification_custody = certificate["custody"]["adopted_clarification_chain"]
    assert len(clarification_custody) == 4
    assert tuple(item["label"] for item in clarification_custody) == tuple(
        label for label, _ in gft._ADOPTED_CLARIFICATION_CHAIN
    )
    assert clarification_custody[0]["sha256"] == (
        "4959f99898b216edc7da3e212ce2e26422287899fcf8f3b41cd34ef5d8bb3ff8"
    )

    certificate_again = gft.gravity_formation_theory_certificate()
    assert certificate_again == certificate and certificate_again is not certificate
    assert _type_error(lambda: _mutate(certificate, "schema", "tampered"))
    assert _type_error(
        lambda: _mutate(certificate["program_authorizations"], "empirical_RGRL_confirmation", True)
    )

    assert tuple(inspect.signature(URM.gravity_formation_theory).parameters) == ()
    assert tuple(inspect.signature(URM.gravity_formation_theory_certificate).parameters) == ()
    delegated = URM.gravity_formation_theory_certificate()
    assert delegated == certificate and delegated is not certificate
    assert _type_error(lambda: URM.gravity_formation_theory(packet={}))
    gravity_role = URM().roles()["GRAVITY"]
    assert "RGRL is adopted, not empirically confirmed" in gravity_role
    assert "off-shell constitutive ancestry" in gravity_role
    assert "give exactly zero on-shell response" in gravity_role
    assert "SPAG force/common-freefall labels are retired" in gravity_role
    assert "complete local diagonal Fisher boundary" in gravity_role
    assert "J8=429 h^8/(16 U_d^7)" in gravity_role
    assert "normalized connected TT four-point" in gravity_role
    alpha_role = URM().roles()["ALPHA"]
    assert "Alpha is not a standalone theory" in alpha_role
    assert "ACTVIS(r,W_obs) together with SAI1--SAI8 implies" in alpha_role
    assert "cannot choose a private record- or region-level alpha" in alpha_role

    with tempfile.TemporaryDirectory(prefix="wac-gft-tamper-") as temporary:
        copied_root = Path(temporary)
        for relative, _ in gft._EXPECTED_ARTIFACTS:
            copied = copied_root / relative
            copied.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(gft._REPOSITORY_ROOT / relative, copied)
        tamper_targets = (
            gft._CLOSURE_PAIR[0][0],
            gft._ADOPTED_CLARIFICATION_CHAIN[0][1][0],
            gft._ADOPTED_CLARIFICATION_CHAIN[3][1][0],
            gft._CORE_SOURCE_PAIRS[4][1][0],
            gft._ADVANCE_SOURCE_PAIRS[3][1][0],
            gft._ADVANCE_SOURCE_PAIRS[4][2][0],
            gft._ADVANCE_SOURCE_PAIRS[6][1][0],
            gft._ADVANCE_SOURCE_PAIRS[8][2][0],
            gft._ADVANCE_SOURCE_PAIRS[9][1][0],
            gft._ADVANCE_SOURCE_PAIRS[10][1][0],
            gft._ADVANCE_SOURCE_PAIRS[10][2][0],
            gft._ADVANCE_SOURCE_PAIRS[11][1][0],
            gft._ADVANCE_SOURCE_PAIRS[11][2][0],
            gft._ADVANCE_SOURCE_PAIRS[13][1][0],
            gft._ADVANCE_SOURCE_PAIRS[13][2][0],
            gft._ADVANCE_SOURCE_PAIRS[15][1][0],
        )
        for relative in tamper_targets:
            target = copied_root / relative
            target.write_bytes(target.read_bytes() + b"\ntampered\n")
            assert _refused(lambda: gft._verify_custody(copied_root))
            shutil.copy2(gft._REPOSITORY_ROOT / relative, target)

    print("GRAVITY_FORMATION_THEORY_GATE: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
