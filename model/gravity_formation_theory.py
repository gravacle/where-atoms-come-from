"""Zero-input URM certificate for the closed working Gravity Formation Theory.

This surface verifies a narrow set of pinned source, adoption, seal, and audit
bytes and reports their exact claim ceilings, including the adopted
off-shell/on-shell clarification.  It accepts no observations, parameters, or
caller-selected roots; it is not a gravity solver, an RGRL experiment, or a
numerical-G calculation.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping, NoReturn


SCHEMA = "WAC_GRAVITY_FORMATION_THEORY_CERTIFICATE_V002"
CLAIM_CLASS = (
    "ADOPTED_RGRL_OFFSHELL_ANCESTRY_AND_EXACT_CONDITIONAL_WORKING_THEORY_"
    "CLOSURE_WITH_SEPARATE_ONSHELL_RESPONSE_CEILING"
)

_REPOSITORY_ROOT = Path(__file__).resolve().parent.parent

_CLOSURE_PAIR = (
    (
        "GRAVITY_FORMATION_THEORY_CLOSURE_V001.md",
        "63c37c85442fa96739591a1380a41ee29a9f6f66a6ca0afd5b3470d22fdce028",
    ),
    (
        "GRAVITY_FORMATION_THEORY_CLOSURE_V001.AUDIT.md",
        "afb3602c07ca4c15e99d3f1d18432d11629f05184b207d4dc09336ef437657c2",
    ),
)

_ADOPTED_CLARIFICATION_CHAIN = (
    (
        "onshell_offshell_clarification_adoption",
        (
            "GRAVITY_RGRL_ONSHELL_OFFSHELL_CLARIFICATION_ADOPTION_V001.md",
            "4959f99898b216edc7da3e212ce2e26422287899fcf8f3b41cd34ef5d8bb3ff8",
        ),
    ),
    (
        "onshell_offshell_clarification_adoption_seal",
        (
            "GRAVITY_RGRL_ONSHELL_OFFSHELL_CLARIFICATION_ADOPTION_V001.md.sha256",
            "6fda740bf85f79e9015121ec29a4a26d1080bf2fbd77acfafa1d5acdc574cc97",
        ),
    ),
    (
        "onshell_offshell_spag_clarification",
        (
            "GRAVITY_RGRL_ONSHELL_OFFSHELL_AND_SPAG_CLARIFICATION_V001.md",
            "50e2be7f06e79d943eddb6c37c63094dd758bb8b6798987e73b8eacf8ef448df",
        ),
    ),
    (
        "local_response_G_identifiability_ceiling",
        (
            "GRAVITY_RGRL_LOCAL_RESPONSE_G_IDENTIFIABILITY_CEILING_V002.md",
            "140eb379f22c369b0442b5585a4828122ac4f5858b54ae3f2bf6391866ac84a3",
        ),
    ),
)

_CORE_SOURCE_PAIRS = (
    (
        "rgrl_adoption",
        (
            "GRAVITY_RGRL_ADOPTION_V001.md",
            "bca6146dfa2f2a32cea42db43c85c5d5fb1ee7e6114206e321066809e7c0db1f",
        ),
        (
            "GRAVITY_RGRL_ADOPTION_V001.AUDIT.md",
            "4998b2ff2bc65a1d7bc1dc1b7df41848e9f254523d44ab6d3a8ba905b75f8768",
        ),
    ),
    (
        "record_geometry_structural_theorem",
        (
            "GRAVITY_RGRL_POST_ADOPTION_STRUCTURAL_THEOREM_V001.md",
            "733b18ecaa29c7acd755db6947b790a9ae37240a3c74d199752d5e278280783d",
        ),
        (
            "GRAVITY_RGRL_POST_ADOPTION_STRUCTURAL_THEOREM_V001.AUDIT.md",
            "ededc22f61b9e592b6dbb49934216fe82559a5ba534432d5a3d46b5e0ce79b2f",
        ),
    ),
    (
        "rgrl_ir_endpoint_closure",
        (
            "GRAVITY_RGRL_IR_ENDPOINT_CLOSURE_THEOREM_V001.md",
            "c883c4c9f3816e453766846a1691ef27cb50d6ea7e5676bc52ed1928617f82bf",
        ),
        (
            "GRAVITY_RGRL_IR_ENDPOINT_CLOSURE_THEOREM_V001.AUDIT.md",
            "c70fab6a3a6e78a2146b525cd01dad1de5a66f546b78523b69d24343a1585309",
        ),
    ),
    (
        "s4_lineage_metric_kernel",
        (
            "GRAVITY_RGRL_S4_LINEAGE_METRIC_RESPONSE_KERNEL_THEOREM_V001.md",
            "49e97e9cd3c9d8c75c65f3717156071bfcc0d88b3be3118aa442f74fb711f50d",
        ),
        (
            "GRAVITY_RGRL_S4_LINEAGE_METRIC_RESPONSE_KERNEL_THEOREM_V001.AUDIT.md",
            "6ae98e479d172037c35860aa7ed792d14ae1bad6b1beef6c14afb7061441e448",
        ),
    ),
    (
        "f3_record_front_lorentz_cone_refinement",
        (
            "GRAVITY_F3_RECORD_FRONT_LORENTZ_CONE_REFINEMENT_THEOREM_V001.md",
            "2220189ba41fbb137bd0a7be86b86e4c536c2fac7100f420eb45c20229612dfa",
        ),
        (
            "GRAVITY_F3_RECORD_FRONT_LORENTZ_CONE_REFINEMENT_THEOREM_V001.AUDIT.md",
            "33a9409a7aa11734ebe1d021fc447c9c0bf69232642d9536a5f1db1e18bbc70b",
        ),
    ),
    (
        "spag_prospective_protocol",
        (
            "GRAVITY_RGRL_SPAG_PROSPECTIVE_PROTOCOL_V001.md",
            "9495ca2b9edf3ebf1133e077d746e77b78ebda6e0fc061c178c80109506386b9",
        ),
        (
            "GRAVITY_RGRL_SPAG_PROSPECTIVE_PROTOCOL_V001.AUDIT.md",
            "0acb23e0795dd67a5e1d03b1c1771006a867c8b46f46c265c8cd183bede4f1c7",
        ),
    ),
)

_EXPECTED_ARTIFACTS = _CLOSURE_PAIR + tuple(
    artifact
    for _, theorem, audit in _CORE_SOURCE_PAIRS
    for artifact in (theorem, audit)
) + tuple(artifact for _, artifact in _ADOPTED_CLARIFICATION_CHAIN)


class GravityFormationTheoryRefusal(RuntimeError):
    """One or more pinned Gravity Formation Theory custody bytes failed."""


def _refuse(message: str) -> NoReturn:
    raise GravityFormationTheoryRefusal(
        "GRAVITY FORMATION THEORY REFUSES: " + message
    )


def _sha256_file(path: Path) -> str:
    if path.is_symlink() or not path.is_file():
        _refuse(f"custody artifact is absent, non-file, or symlinked: {path.name}")
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as exc:
        _refuse(f"custody artifact is unreadable: {path.name}: {exc}")


@dataclass(frozen=True)
class _Custody:
    artifacts: tuple[tuple[str, str], ...]

    def digest_for(self, relative: str) -> str:
        return dict(self.artifacts)[relative]


def _verify_custody(root: Path = _REPOSITORY_ROOT) -> _Custody:
    """Verify every pinned artifact twice.  The root parameter is private/test-only."""
    root = Path(root)
    observed: list[tuple[str, str]] = []
    for relative, expected in _EXPECTED_ARTIFACTS:
        digest = _sha256_file(root / relative)
        if digest != expected:
            _refuse(f"custody artifact hash mismatch: {relative}")
        observed.append((relative, digest))

    # Close the narrow read-time race without trusting mutable sidecars or a
    # broad repository manifest.
    for relative, expected in _EXPECTED_ARTIFACTS:
        if _sha256_file(root / relative) != expected:
            _refuse(f"custody artifact changed during verification: {relative}")

    return _Custody(tuple(observed))


def _freeze(value: Any) -> Any:
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, tuple):
        return tuple(_freeze(item) for item in value)
    return value


def _pair_certificate(
    custody: _Custody,
    label: str,
    theorem: tuple[str, str],
    audit: tuple[str, str],
) -> dict[str, str]:
    theorem_path, _ = theorem
    audit_path, _ = audit
    return {
        "label": label,
        "source_path": theorem_path,
        "source_sha256": custody.digest_for(theorem_path),
        "audit_path": audit_path,
        "audit_sha256": custody.digest_for(audit_path),
    }


def _artifact_certificate(
    custody: _Custody,
    label: str,
    artifact: tuple[str, str],
) -> dict[str, str]:
    path, _ = artifact
    return {
        "label": label,
        "path": path,
        "sha256": custody.digest_for(path),
    }


def _certificate(custody: _Custody) -> Mapping[str, Any]:
    closure_path, _ = _CLOSURE_PAIR[0]
    closure_audit_path, _ = _CLOSURE_PAIR[1]
    return _freeze(
        {
            "schema": SCHEMA,
            "claim_class": CLAIM_CLASS,
            "rgrl": {
                "adoption_id": "GRAVITY-POSTULATE-RGRL-V001",
                "status": "ADOPTED_WORKING_RECORD_GEOMETRY_REALIZATION_POSTULATE",
                "clauses": ("RGRL-A", "RGRL-B", "RGRL-C"),
                "record_to_mathematical_geometry": "AXIOMATICALLY_AUTHORIZED_IN_THE_WORKING_THEORY",
                "RGRL_C_scope": "OFFSHELL_CONSTITUTIVE_SPATIAL_METRIC_TANGENT_ANCESTRY",
                "RGRL_C_is_onshell_force_law": False,
                "physical_instantiation": "OPEN",
                "empirical_confirmation": "OPEN_AND_FALSIFIABLE",
            },
            "adopted_response_clarification": {
                "adoption_id": (
                    "GRAVITY-RGRL-ONSHELL-OFFSHELL-CLARIFICATION-ADOPTION-V001"
                ),
                "status": "ADOPTED_PROSPECTIVE_PROGRAM_POLICY",
                "offshell_tangent": "D_LCOORD_G_FULL_RANK_OR_DECLARED_DENSE_RANGE",
                "onshell_kernel": "SEPARATELY_ESTABLISHED_OR_MEASURED_H_R",
                "open_type_join": "GI21_NOT_SUPPLIED_BY_RGRL_C",
                "fully_matched_response": (
                    "ZERO_ON_V002_WELL_POSED_RETARDED_QUOTIENT_WITH_NO_UNRESOLVED_ZERO_MODE_"
                    "WHEN_DRESSED_SOURCE_REMAINDER_AND_PHYSICAL_DATA_DERIVATIVES_VANISH"
                ),
                "fully_matched_response_premises": (
                    "V002_WELL_POSED_RETARDED_GAUGE_BOUNDARY_QUOTIENT",
                    "NO_UNRESOLVED_ZERO_MODE",
                ),
                "fully_matched_H_R": "MAY_VANISH_WITHOUT_REDUCING_OFFSHELL_RANK",
                "derivative_ownership": (
                    "IMPLICIT_RESPONSE_IN_DRESSED_OPERATOR_ONCE",
                    "EXPLICIT_SOURCE_OR_REMAINDER_DERIVATIVE_ONCE",
                    "INITIAL_INCOMING_BOUNDARY_OR_BRANCH_DATA_DERIVATIVE_ONCE",
                ),
            },
            "conditional_closure": {
                "closure_id": "GFT-CLOSURE-V001",
                "status": "EXACT_TRANSITIVE_CONDITIONAL_SYNTHESIS",
                "premises": (
                    "AURFT_UDCL",
                    "ADOPTED_RGRL_A_THROUGH_C",
                    "TYPED_HEALTHY_SAME_PARENT_EIR_1_THROUGH_6",
                ),
                "scope": "DECLARED_CONNECTED_INFRARED_DOMAIN",
                "common_physical_metric": "EXACT_CONDITIONAL_ON_EIR_1_NOT_RGRL_ALONE",
                "leading_nonlinear_Einstein_response": "EXACT_CONDITIONAL_ON_EIR_1_THROUGH_6_IN_DECLARED_LOCAL_METRIC_ONLY_CLASS",
            },
            "exact_results": {
                "first_record_gamma_seed": "EXACT_FOR_A_QUALIFIED_POSITIVE_MARGIN_RECORD",
                "causal_volume_metric": "EXACT_INSIDE_ADOPTED_RGRL_A",
                "six_mode_lineage_ancestry": "EXACT_INSIDE_ADOPTED_RGRL_B_C",
                "six_mode_ancestry_kind": "OFFSHELL_FORMAL_SPATIAL_METRIC_TANGENT",
                "common_physical_metric": "EXACT_CONDITIONAL_ON_EIR_1",
                "leading_nonlinear_Einstein_response": "EXACT_CONDITIONAL_ON_FULL_EIR_1_THROUGH_6_IN_DECLARED_IR_CLASS",
                "record_front_lorentz_cone_refinement": (
                    "EXACT_CONDITIONAL_ON_AFR_AND_MICROSCOPIC_EQUAL_DEPTH_EQUAL_QJ_NORM_NULL_STEPS"
                ),
                "refined_cone_kind": "LOCAL_MATHEMATICAL_3PLUS1_LORENTZ_CONE",
                "six_mode_cone_tangent_scope": "J0_S4_SYMMETRIC_POINT",
                "physical_volume_and_common_probe_identification": "OPEN",
                "fixed_finite_additive_direction_set_obstruction": "EXACT",
                "cone_refinement_is_full_RGRL_or_gravity_dynamics": False,
                "scalar_gamma_or_count_determines_kernel_curvature_or_G": False,
            },
            "kernel_reduction": {
                "scope": "Q4_TETRAHEDRAL_S4_FIXED_POINT_LOCAL_OR_K_ZERO",
                "pair_tangent_dimension": 6,
                "S4_sectors": ("A1_TRACE_1", "E2_DIAGONAL_SHEAR_2", "T2_OFF_DIAGONAL_SHEAR_3"),
                "S4_response_form_factors": ("h_A", "h_E", "h_T"),
                "compatible_O3_reduction": ("h_trace", "h_shear"),
                "compatible_source_and_output_O3_required": True,
                "nonzero_spatial_momentum_classification": "OPEN",
                "form_factor_measurement": "OPEN",
                "absolute_lineage_and_length_calibration": "OPEN",
                "offshell_metric_tangent_rank": "FULL_ON_ADMITTED_SIX_DIRECTIONS",
                "onshell_H_R_nonzero": "NOT_IMPLIED_AND_MAY_VANISH_IN_FULLY_MATCHED_LANE",
                "GI21_compatibility_type_join": "OPEN",
            },
            "numerical_G": {
                "positive_total_endpoint_coefficient": "CONDITIONAL_RESULT_IN_THE_TYPED_ENDPOINT",
                "endpoint_matching": "G_EFF_EQUALS_G_END",
                "value_status": "OBSERVED_ENDPOINT_CALIBRATION_REQUIRED",
                "identification_equation": "FULL_IMPLICIT_DRESSED_FINITE_APPARATUS_FORWARD_MODEL",
                "literal_K0_ratio": "NO_UNIQUE_IDENTIFICATION_WITHOUT_WELL_POSED_QUOTIENT",
                "bounded_remainder_output": "IDENTIFIED_INTERVAL_OR_SET_NOT_AUTOMATIC_POINT_VALUE",
                "calibration_performed_by_this_certificate": False,
                "parameter_free_microscopic_record_derivation": False,
                "gamma_or_record_count_calculates_G": False,
            },
            "SPAG": {
                "status": "HISTORICAL_AUDITED_PROTOCOL_SUPERSEDED_FOR_FUTURE_ONSHELL_USE",
                "load_bearing_for_closure": False,
                "historical_protocol_executed": False,
                "old_local_RGRL_C_force_column_labels": "RETIRED_FOR_FUTURE_CLAIMS",
                "offshell_RGRL_C_direct_test": False,
                "lane_A": "COMPLETE_SOURCE_MATCHED_DISCOVERY_WITH_ZERO_PHYSICAL_PREDICTION",
                "lane_B": "INDEPENDENT_SOURCE_CALIBRATED_G_CROSS_CHECK",
                "source_unresolved_positive": "ANCESTRY_CORRELATED_UNCLASSIFIED_ANOMALY",
                "target_local_coefficient": False,
                "six_mode_rank_or_dense_range": False,
                "full_RGRL_confirmation": False,
                "executed_by_this_surface": False,
                "confirmatory_requirement": "NEW_HASHED_ACQUISITION_DISJOINT_RUN_B_AND_ALL_FROZEN_GATES",
            },
            "scientific_status": {
                "nature_obeys_RGRL": "NOT_ESTABLISHED_BY_THIS_CERTIFICATE",
                "EIR_packet": "TYPED_CONDITIONAL_INPUT_NOT_ESTABLISHED_BY_THIS_CERTIFICATE",
                "empirical_RGRL_confirmation": "OPEN_NOT_PERFORMED",
                "empirical_lineage_gravity_confirmation": "OPEN_NOT_PERFORMED",
                "response_form_factor_measurement": "OPEN_NOT_PERFORMED",
                "numerical_G_calibration": "REQUIRES_OBSERVATION_NOT_CALLER_INPUT",
                "deeper_microscopic_RGRL_derivation": "OPEN",
                "AFR_and_null_step_physical_instantiation": "OPEN",
                "GI21_compatibility_type_join": "OPEN",
                "lineage_source_functional": "OPEN",
                "caller_input_scientific_weight": "ZERO",
            },
            "program_authorizations": {
                "working_RGRL_postulate_adopted": True,
                "exact_conditional_working_theory_closure": True,
                "exact_conditional_EIR_composition": True,
                "offshell_onshell_clarification_adopted": True,
                "empirical_RGRL_confirmation": False,
                "empirical_EIR_confirmation": False,
                "SPAG_executed": False,
                "old_SPAG_local_RGRL_C_force_labels_authorized": False,
                "full_RGRL_experimentally_closed": False,
                "numerical_G_derived_from_records": False,
                "gravity_solver": False,
                "caller_data_admitted": False,
            },
            "nonpromotion": {
                "adoption_promotes_empirical_confirmation": False,
                "conditional_EIR_composition_promotes_executed_EIR_packet": False,
                "kernel_reduction_promotes_measured_form_factors": False,
                "positive_G_sign_promotes_numerical_G_derivation": False,
                "SPAG_protocol_promotes_SPAG_observation_or_full_RGRL": False,
                "offshell_rank_promotes_nonzero_onshell_H_R_or_force": False,
                "conditional_cone_refinement_promotes_AFR_or_gravity_dynamics": False,
                "custody_promotes_physical_evidence": False,
                "conditional_closure_promotes_unconditional_actual_world_theorem": False,
                "S4_symmetry_promotes_compatible_O3": False,
            },
            "custody": {
                "verification": "TWO_PASS_EXACT_SHA256_WITH_SYMLINK_REFUSAL",
                "artifact_count": len(_EXPECTED_ARTIFACTS),
                "closure": {
                    "source_path": closure_path,
                    "source_sha256": custody.digest_for(closure_path),
                    "audit_path": closure_audit_path,
                    "audit_sha256": custody.digest_for(closure_audit_path),
                    "audit_verdict": "CLEAN_WITHIN_EXPLICIT_CONDITIONAL_AND_EMPIRICAL_CEILINGS",
                },
                "core_source_pairs": tuple(
                    _pair_certificate(custody, label, theorem, audit)
                    for label, theorem, audit in _CORE_SOURCE_PAIRS
                ),
                "adopted_clarification_chain": tuple(
                    _artifact_certificate(custody, label, artifact)
                    for label, artifact in _ADOPTED_CLARIFICATION_CHAIN
                ),
            },
            "executable_scope": {
                "caller_arguments": 0,
                "solver": False,
                "observation_loader": False,
                "experiment_performed": False,
                "theorem_machine_proved": False,
                "scientific_output": "PINNED_DOCUMENTARY_CUSTODY_AND_STATUS_CERTIFICATE_ONLY",
            },
        }
    )


@dataclass(frozen=True)
class GravityFormationTheory:
    """Immutable handle produced by one fresh verification of all pinned artifacts."""

    _custody: _Custody

    @property
    def claim_class(self) -> str:
        return CLAIM_CLASS

    @property
    def closure_sha256(self) -> str:
        return self._custody.digest_for(_CLOSURE_PAIR[0][0])

    @property
    def audit_sha256(self) -> str:
        return self._custody.digest_for(_CLOSURE_PAIR[1][0])

    def certificate(self) -> Mapping[str, Any]:
        """Return a fresh recursively immutable certificate; run no physics."""
        return _certificate(self._custody)


def gravity_formation_theory() -> GravityFormationTheory:
    """Verify and expose the pinned working theory with zero caller input."""
    return GravityFormationTheory(_verify_custody())


def gravity_formation_theory_certificate() -> Mapping[str, Any]:
    """Return the immutable zero-input status certificate."""
    return gravity_formation_theory().certificate()
