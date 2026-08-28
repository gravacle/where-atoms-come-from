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
    assert gft.SCHEMA == "WAC_GRAVITY_FORMATION_THEORY_CERTIFICATE_V007"
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
    assert certificate["exact_results"]["q4_finite_tt_periodic_quotient"] == (
        "30_CELLS_60_VERTICES_120_LINKS_120_ELEMENTARY_HEXAGONS"
    )
    assert certificate["exact_results"]["q4_finite_tt_ring_sector"] == (
        "180_STATES_420_TRANSITIONS_TRANSLATION_CLOSED_H6_ONLY"
    )
    assert certificate["exact_results"]["q4_finite_tt_composite_cumulant"] == (
        "W2_1P130847135995723_OVER_J6__W4_MINUS_0P136825085605100_OVER_J6_CUBED__"
        "GAMMA4_COMP_0P083666214307836_J6"
    )
    assert certificate["exact_results"]["q4_finite_tt_lowest_pole_vs_two_link_proxy"] == (
        "3P194109035554332_J6_ABOVE_2P059674505691458_J6"
    )
    assert certificate["exact_results"]["q4_projected_ice_constraint_species"] == (
        "ONE_SCALAR_U1_GAUSS_ONLY"
    )
    assert certificate["exact_results"]["q4_pair_relations_constraint_status"] == (
        "ALGEBRAIC_ZERO_OPERATORS_NOT_NEW_FIRST_CLASS_GENERATORS"
    )
    assert certificate["exact_results"]["q4_A3_static_cometric_tangent_rank"] == 6
    assert certificate["exact_results"]["q4_inherited_even_bulk_Kubo_dynamic_rank"] == 0
    assert certificate["exact_results"]["q4_fixed_parent_collective_metric_origin"] == (
        "NO_PRESENT_OBJECT_JOINTLY_OWNS_SIX_CONFIG_CHANNELS_"
        "INDEPENDENT_CONJUGATES_AND_VECTOR3_PLUS_SCALAR_NULL_PACKET"
    )
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
    assert microscopic["finite_TT_composite_precursor"] == (
        "EXACT_STATIC_TWO_Q_AND_FOUR_Q_CUMULANTS_AND_FINITE_SPECTRUM"
    )
    assert microscopic["connected_four_one_link_channel_2PI"] == "OPEN_NOT_COMPUTED_BY_FO"
    assert microscopic["inherited_projected_constraint_architecture"] == "SCALAR_U1_ONLY"
    assert microscopic["microscopic_RGRLB_from_current_q4_ice"] == "NOT_DERIVED"
    assert microscopic["six_A3_static_deformation_coefficients"] == "EXACT_RANK6"
    assert microscopic["six_A3_dynamic_metric_fields"] == (
        "ABSENT_IN_INHERITED_EVEN_BULK_SCREEN"
    )
    assert microscopic["q4_block_strain_CTP_status"] == (
        "FROZEN_NEXT_CALCULATION_NOT_EXECUTED"
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
    assert certificate["kernel_reduction"]["finite_TT_composite_cumulant_kind"] == (
        "FOUR_Q_EIGHT_ONE_LINK_NOT_CONNECTED_FOUR_ONE_LINK_OR_PHOTON_AMPUTATED"
    )
    assert certificate["kernel_reduction"]["finite_TT_pole_screen"] == (
        "NO_BELOW_PROXY_OR_ENERGY_EXCLUSIVE_CANDIDATE__NOT_A_NO_BOUND_STATE_THEOREM"
    )
    assert certificate["kernel_reduction"]["finite_TT_composite_legendre_quartic"] == (
        "POSITIVE_IN_BOTH_SUSCEPTIBILITY_EIGENCHANNELS"
    )
    assert certificate["kernel_reduction"]["projected_ice_constraint_packet"] == (
        "ONE_SCALAR_U1_NOT_RGRLB"
    )
    assert certificate["kernel_reduction"]["projected_ice_equal_two_polarizations_implies_RGRLB"] is False
    assert certificate["kernel_reduction"]["A3_static_metric_span"] == (
        "RANK6_COEFFICIENT_TANGENT"
    )
    assert certificate["kernel_reduction"]["A3_inherited_even_source_Kubo_rank"] == 0
    assert certificate["kernel_reduction"]["FJ_unprojected_pair_response"] == (
        "EXACT_CONDITIONAL_RANK6_A1_E_T2_RESPONSE_WITH_NEAREST_CELL_SPREADING"
    )
    assert certificate["kernel_reduction"]["current_fixed_parent_metric_packet"] == (
        "NO_JOINT_SIX_CONFIG_PLUS_CONJUGATE_PLUS_VECTOR3_SCALAR_NULL_OBJECT_"
        "IN_CONSTRUCTED_CATALOG__NOT_THERMODYNAMIC_NO_GO"
    )
    assert certificate["kernel_reduction"]["next_collective_metric_calculation"] == (
        "Q4_BLOCK_STRAIN_CTP_SOURCE_BEFORE_FESHBACH_WITH_H_SOURCE_OFF_UNCHANGED"
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
    assert certificate["numerical_G"]["hust_public_forward_status"] == (
        "PROCESSED_DUAL_CHANNEL_FORWARD_CLOSED__FULL_GC16_NOT_READY"
    )
    assert certificate["numerical_G"]["hust_tos_A_B_A_delta_omega2_s_minus2"] == (
        "1.6626945111323172e-6"
    )
    assert certificate["numerical_G"]["hust_tos_quadratic_delta_omega2_s_minus2"] == (
        "1.6626989120180067e-6"
    )
    assert certificate["numerical_G"]["hust_aaf_processed_forward_count"] == 3
    assert certificate["numerical_G"]["hust_aaf_max_rounding_difference_ppm"] == (
        "LESS_THAN_0P2"
    )
    assert certificate["numerical_G"]["hust_cross_method_separation_ppm"] == (
        "44.9483_DESCRIPTIVE_ONLY"
    )
    assert certificate["numerical_G"]["hust_cross_method_z_condition"] == (
        "2P7196_ONLY_IF_CROSS_COVARIANCE_ZERO"
    )
    assert certificate["numerical_G"]["hust_accepted_G_input_used"] is False
    assert certificate["numerical_G"]["hust_full_GC16_fit_ready"] is False
    assert certificate["numerical_G"]["hust_independent_G_crosscheck_performed"] is False
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
    assert certificate["program_authorizations"]["finite_TT_composite_cumulant_and_spectrum_derived"] is True
    assert certificate["program_authorizations"]["connected_four_one_link_channel_2PI_derived"] is False
    assert certificate["program_authorizations"]["q4_scalar_U1_constraint_only_derived"] is True
    assert certificate["program_authorizations"]["microscopic_RGRLB_from_current_q4_ice_derived"] is False
    assert certificate["program_authorizations"]["six_static_cometric_deformations_derived"] is True
    assert certificate["program_authorizations"]["six_dynamic_metric_fields_from_inherited_even_bulk_derived"] is False
    assert certificate["program_authorizations"]["q4_block_strain_CTP_executed"] is False
    assert certificate["program_authorizations"]["hust_processed_dual_channel_forward_executed"] is True
    assert certificate["program_authorizations"]["hust_full_GC16_executed"] is False
    assert certificate["program_authorizations"]["inherited_protected_tensor_pole_derived"] is False
    assert (
        certificate["scientific_status"]["inherited_non_gaussian_TT_kernel"]
        == "FINITE_COMPOSITE_CUMULANT_AND_SPECTRUM_COMPUTED__"
        "CONNECTED_FOUR_ONE_LINK_CHANNEL_2PI_OPEN"
    )
    assert certificate["scientific_status"]["normalized_connected_TT_four_point"] == (
        "OPEN__FO_COMPUTED_FOUR_Q_COMPOSITE_NOT_THIS_OBJECT"
    )
    assert certificate["scientific_status"]["nonperturbative_TT_pole_Ward_common_cone_residue"] == "OPEN"
    assert certificate["scientific_status"]["current_q4_ice_RGRLB_constraint_origin"] == (
        "CLOSED_NEGATIVE_FOR_CURRENT_FINITE_LOCAL_Q4_ICE_BRANCH__"
        "THERMODYNAMIC_COLLECTIVE_ORIGIN_OPEN"
    )
    assert certificate["scientific_status"]["current_fixed_parent_collective_metric_origin"] == (
        "CLOSED_NEGATIVE_FOR_CONSTRUCTED_OBJECT_CATALOG__NOT_THERMODYNAMIC_NO_GO"
    )
    assert certificate["scientific_status"]["next_no_lab_gravity_calculation"] == (
        "Q4_BLOCK_STRAIN_CTP"
    )
    assert certificate["scientific_status"]["hust_public_G_forward"] == (
        "PROCESSED_DUAL_CHANNEL_CLOSED__FULL_GC16_NOT_READY__NO_NEW_G"
    )
    assert (
        certificate["scientific_status"]["local_ice_fisher_T2_solder"]
        == "CLOSED_NEGATIVE_FOR_LOCAL_DIAGONAL_COMPLEMENT_PRESERVING_FAMILIES"
    )
    assert certificate["nonpromotion"]["second_fisher_jet_promotes_linear_metric_tangent"] is False
    assert certificate["nonpromotion"]["broken_vector_background_rank_promotes_tensor_mode"] is False
    assert certificate["nonpromotion"]["control_sign_conditioning_promotes_endogenous_metric_solder"] is False
    assert certificate["nonpromotion"]["bare_vertex_analyticity_promotes_dressed_1PI_analyticity"] is False
    assert certificate["nonpromotion"]["strict_single_insertion_boundary_promotes_nonperturbative_no_pole"] is False
    assert certificate["nonpromotion"]["four_Q_composite_cumulant_promotes_four_link_1PI_or_binding"] is False
    assert certificate["nonpromotion"]["finite_threshold_proxy_promotes_thermodynamic_no_bound_state"] is False
    assert certificate["nonpromotion"]["equal_two_polarizations_promotes_RGRLB"] is False
    assert certificate["nonpromotion"]["six_static_deformation_coefficients_promote_dynamic_metric_fields"] is False
    assert certificate["nonpromotion"]["current_catalog_obstruction_promotes_thermodynamic_no_go"] is False
    assert certificate["nonpromotion"]["hust_processed_forward_promotes_new_G_or_GFT_confirmation"] is False
    assert certificate["custody"]["artifact_count"] == 58
    assert len(certificate["custody"]["core_source_pairs"]) == 6
    assert (
        certificate["custody"]["core_source_pairs"][4]["label"]
        == "f3_record_front_lorentz_cone_refinement"
    )
    advance_custody = certificate["custody"]["advance_source_pairs"]
    assert len(advance_custody) == 20
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
    assert advance_custody[12]["label"] == "f3_q4_finite_tt_composite_cumulant_screen"
    assert advance_custody[12]["source_sha256"] == (
        "44fc28edc9820d2b4ea67cef9f83beef60e53bfe6407b582ffbeccfe42f756c5"
    )
    assert advance_custody[12]["audit_sha256"] == (
        "84d8c02c3e560198f6a9fae04f5ee81bc72c98354e02f1f58d506a8c3171c453"
    )
    assert advance_custody[13]["label"] == "f3_q4_rgrlb_constraint_origin_screen"
    assert advance_custody[13]["audit_sha256"] == (
        "48350331885e06c5b4fbcd4fa21ccaa876b499e6bb0fd27a644a4fef9574c8e0"
    )
    assert advance_custody[14]["label"] == "f3_q4_collective_metric_origin_screen"
    assert advance_custody[14]["audit_sha256"] == (
        "91aa35170432684a47278e46ee2b9d56658a43acc8acbbb84480d047cdbe6dcf"
    )
    assert advance_custody[16]["label"] == "spag_public_data_second_pass"
    assert advance_custody[16]["audit_path"] == (
        "LANE_GRA_SPAG_PUBLIC_DATA_SECOND_PASS_V001/POST_INTEGRATION_CUSTODY_REAUDIT.md"
    )
    assert advance_custody[16]["audit_sha256"] == (
        "dcc170c2674fe0020f679e523ab40689be7efd32e05e71fe7600d9bbf4047e97"
    )
    assert advance_custody[18]["label"] == "nist_bipm_g_forward_readiness"
    assert advance_custody[19]["label"] == "hust_2018_dual_method_g_forward"
    assert advance_custody[19]["source_sha256"] == (
        "44eb8c81a3d84dfa6829bcd6971d0261215877af0529318eab5cedbd3980c340"
    )
    assert advance_custody[19]["audit_sha256"] == (
        "e8625b4cbf67d73927db495e8111e3ffbd4e85f46b80183e55fbf8b4391d0b2e"
    )
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
    assert "local Fisher shortcut fails" in gravity_role
    assert "J8=429 h^8/(16 U_d^7)" in gravity_role
    assert "four bilinear Q insertions are an eight-one-link object" in gravity_role
    assert "one scalar compact-U1 Gauss species" in gravity_role
    assert "retarded root-source rank zero" in gravity_role
    assert "Q4-BLOCK-STRAIN-CTP" in gravity_role
    assert "HUST-2018" in gravity_role and "full GC16 or a new G" in gravity_role
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
            gft._ADVANCE_SOURCE_PAIRS[12][1][0],
            gft._ADVANCE_SOURCE_PAIRS[12][2][0],
            gft._ADVANCE_SOURCE_PAIRS[13][1][0],
            gft._ADVANCE_SOURCE_PAIRS[13][2][0],
            gft._ADVANCE_SOURCE_PAIRS[14][1][0],
            gft._ADVANCE_SOURCE_PAIRS[14][2][0],
            gft._ADVANCE_SOURCE_PAIRS[16][1][0],
            gft._ADVANCE_SOURCE_PAIRS[16][2][0],
            gft._ADVANCE_SOURCE_PAIRS[18][1][0],
            gft._ADVANCE_SOURCE_PAIRS[19][1][0],
            gft._ADVANCE_SOURCE_PAIRS[19][2][0],
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
