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
    assert microscopic["new_interaction_or_second_field_adopted"] is False
    assert certificate["numerical_G"]["parameter_free_microscopic_record_derivation"] is False
    assert certificate["numerical_G"]["bounded_remainder_output"] == "IDENTIFIED_INTERVAL_OR_SET_NOT_AUTOMATIC_POINT_VALUE"
    assert certificate["numerical_G"]["calibrated_nonsingular_row_identifies"] == "P_EQUALS_G_TIMES_SOURCE_SCALE"
    assert certificate["numerical_G"]["synthetic_validation"] == "PASS_15_OF_15_WITH_ARBITRARY_NONEMPIRICAL_G"
    assert certificate["SPAG"]["executed_by_this_surface"] is False
    assert certificate["SPAG"]["full_RGRL_confirmation"] is False
    assert certificate["SPAG"]["old_local_RGRL_C_force_column_labels"] == "RETIRED_FOR_FUTURE_CLAIMS"
    assert certificate["SPAG"]["lane_A"] == "COMPLETE_SOURCE_MATCHED_DISCOVERY_WITH_ZERO_PHYSICAL_PREDICTION"
    assert certificate["SPAG"]["lane_B"] == "INDEPENDENT_SOURCE_CALIBRATED_G_CROSS_CHECK"
    assert certificate["SPAG"]["public_data_substitute_executed"] is True
    assert certificate["SPAG"]["public_same_parent_eight_cell_support"] == "ABSENT"
    assert certificate["SPAG"]["public_beta_TM_identifiable"] is False
    assert certificate["program_authorizations"]["old_SPAG_local_RGRL_C_force_labels_authorized"] is False
    assert certificate["custody"]["artifact_count"] == 34
    assert len(certificate["custody"]["core_source_pairs"]) == 6
    assert (
        certificate["custody"]["core_source_pairs"][4]["label"]
        == "f3_record_front_lorentz_cone_refinement"
    )
    advance_custody = certificate["custody"]["advance_source_pairs"]
    assert len(advance_custody) == 8
    assert tuple(item["label"] for item in advance_custody) == tuple(
        label for label, _, _ in gft._ADVANCE_SOURCE_PAIRS
    )
    assert advance_custody[3]["label"] == "f3_q4_carrier_lift_derivability_boundary"
    assert advance_custody[4]["label"] == "q4_pair_field_lift_derivability_boundary"
    assert advance_custody[5]["label"] == "f3_q4_finite_programmed_support_solder"
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
            gft._ADVANCE_SOURCE_PAIRS[7][1][0],
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
