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


SCHEMA = "WAC_GRAVITY_FORMATION_THEORY_CERTIFICATE_V008"
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

_ADVANCE_SOURCE_PAIRS = (
    (
        "f3_q4_clifford_collective_cone",
        (
            "LANE_GRA_FC_F3_Q4_CLIFFORD_COLLECTIVE_CONE_V001/THEOREM.md",
            "28b6319e3187337da8ebef2212b030ff6e5b9f8168d9844ae172d94f3e0641a6",
        ),
        (
            "LANE_GRA_FC_F3_Q4_CLIFFORD_COLLECTIVE_CONE_V001/REAUDIT.md",
            "0dc9acedd88f8d0ff22c94df747f7d0a529fc6920dada73ebb3da3b52d1b3789",
        ),
    ),
    (
        "f3_q4_common_child_acoustic_cone",
        (
            "LANE_GRA_FD_F3_Q4_COMMON_CHILD_ACOUSTIC_CONE_V001/THEOREM.md",
            "60d012766675c12e82dd1731e202a6c0ed48f24e2697f589b63eecc3cb650287",
        ),
        (
            "LANE_GRA_FD_F3_Q4_COMMON_CHILD_ACOUSTIC_CONE_V001/REAUDIT.md",
            "03ba3b2a5d66a40f54b5c2d9e6fe52535a3fad81e196aed2e3bbe7558fe63510",
        ),
    ),
    (
        "f3_q4_diamond_ice_support_join",
        (
            "LANE_GRA_FE_F3_Q4_DIAMOND_ICE_CARRIER_JOIN_V001/THEOREM.md",
            "4cc63e3e5853b4250a2a5b78256d41b83b195cf527901819a62a04ef53f8d932",
        ),
        (
            "LANE_GRA_FE_F3_Q4_DIAMOND_ICE_CARRIER_JOIN_V001/INDEPENDENT_AUDIT.md",
            "9d7ef0419b3022dba0db1add7a46d145ebe4b6ec035f73a9b760e63b978d1b2b",
        ),
    ),
    (
        "f3_q4_carrier_lift_derivability_boundary",
        (
            "LANE_GRA_FF_F3_Q4_CARRIER_LIFT_DERIVABILITY_NO_GO_V001/THEOREM.md",
            "4c5d476e007f36b20f3e34964607c013ab28ae16bf884c063b7f4ac954178e5a",
        ),
        (
            "LANE_GRA_FF_F3_Q4_CARRIER_LIFT_DERIVABILITY_NO_GO_V001/INDEPENDENT_AUDIT.md",
            "655e4a0f90953cc71cd6a12175d9b1d243a6f63d8d1c73140dd6ce3a426a90be",
        ),
    ),
    (
        "q4_pair_field_lift_derivability_boundary",
        (
            "LANE_GRA_FG_Q4_PAIR_FIELD_LIFT_DERIVABILITY_V001/THEOREM.md",
            "fff521ae41e3f8b83a4738ff96a99715e89f90e2d64724786da8a3ed4732e838",
        ),
        (
            "LANE_GRA_FG_Q4_PAIR_FIELD_LIFT_DERIVABILITY_V001/INDEPENDENT_AUDIT.md",
            "ab2d4a0a7e10186973bc63043e8b31a7bc6f351c92cc24543de9be99706f13ef",
        ),
    ),
    (
        "f3_q4_finite_programmed_support_solder",
        (
            "LANE_GRA_FH_F3_Q4_FINITE_PROGRAMMED_SUPPORT_SOLDER_V001/THEOREM.md",
            "2b88febc569efa0de0238e8000d018bf3f798a8ebed2e4ff1327f053d6bd9284",
        ),
        (
            "LANE_GRA_FH_F3_Q4_FINITE_PROGRAMMED_SUPPORT_SOLDER_V001/INDEPENDENT_REAUDIT.md",
            "5c275748d54743ef44098f74c4c5698aead0845d51e6c2dcf32a1bef63f0c7bf",
        ),
    ),
    (
        "f3_q4_programmed_floquet_detuning",
        (
            "LANE_GRA_FI_F3_Q4_PROGRAMMED_FLOQUET_DETUNING_V001/THEOREM.md",
            "09a9e2ee46acf10dbde91e9578576cb537fe5aff4a9dea513d4c1f208e62de4c",
        ),
        (
            "LANE_GRA_FI_F3_Q4_PROGRAMMED_FLOQUET_DETUNING_V001/INDEPENDENT_REAUDIT.md",
            "7ed8f5ddff642fa45c36718be0346ee25280ed20f1e04a2ec2038393654ea244",
        ),
    ),
    (
        "f3_q4_authenticated_link_pair_response",
        (
            "LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/THEOREM.md",
            "05f4a619a6f80aa40c48570ab4035ab874426502a31468a08f435e66610bd769",
        ),
        (
            "LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/INDEPENDENT_AUDIT.md",
            "44690ad431c85af7a4947a431a4c57ad4cd8b19a346e3d26720b341f77256f90",
        ),
    ),
    (
        "f3_q4_ice_hybrid_tensor_response",
        (
            "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/THEOREM.md",
            "cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98",
        ),
        (
            "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/INDEPENDENT_AUDIT.md",
            "c52eab9d701d1c6e82f1d7ec395841f4d2810e96cccbc3e2504760b6742e81e4",
        ),
    ),
    (
        "f3_q4_maxwell_composite_pole_screen",
        (
            "LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/THEOREM.md",
            "98e2b3bc7a1c998d7839dc1a6b435cc1c8ed6d5a622ba45f63571be9ef646452",
        ),
        (
            "LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/INDEPENDENT_AUDIT.md",
            "327bf6a4476c4c6382757dc156a96c6032233d34c25c1f7935e2582acf6c607a",
        ),
    ),
    (
        "f3_q4_inherited_tt_kernel_boundary",
        (
            "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md",
            "78f0687c9f597c96d235224dad45c204d12d7e6c973c270e3305a901efc75b25",
        ),
        (
            "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/INDEPENDENT_AUDIT.md",
            "53893c7198241f0f8f6aa766f3595fb75b83d208581833c32656b28d7c7f02b9",
        ),
    ),
    (
        "f3_q4_ice_t2_fisher_solder_boundary",
        (
            "LANE_GRA_FN_F3_Q4_ICE_T2_FISHER_SOLDER_BOUNDARY_V001/THEOREM.md",
            "be69f15d611827db9841bd932042604deb4f82a777ff9da28b80e4493cef7596",
        ),
        (
            "LANE_GRA_FN_F3_Q4_ICE_T2_FISHER_SOLDER_BOUNDARY_V001/INDEPENDENT_AUDIT.md",
            "32297dc0c4b0454c4a4be88d3763eb679b4ca89bb2385010ba8c2b77f2df47d2",
        ),
    ),
    (
        "f3_q4_finite_tt_composite_cumulant_screen",
        (
            "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001/THEOREM.md",
            "44fc28edc9820d2b4ea67cef9f83beef60e53bfe6407b582ffbeccfe42f756c5",
        ),
        (
            "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001/INDEPENDENT_AUDIT.md",
            "84d8c02c3e560198f6a9fae04f5ee81bc72c98354e02f1f58d506a8c3171c453",
        ),
    ),
    (
        "f3_q4_rgrlb_constraint_origin_screen",
        (
            "LANE_GRA_FP_F3_Q4_RGRLB_CONSTRAINT_ORIGIN_SCREEN_V001/THEOREM.md",
            "8c4cf9d29e48b7116be05d05e6a11513e662443882390aa4f61cd194f5c23dc6",
        ),
        (
            "LANE_GRA_FP_F3_Q4_RGRLB_CONSTRAINT_ORIGIN_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md",
            "48350331885e06c5b4fbcd4fa21ccaa876b499e6bb0fd27a644a4fef9574c8e0",
        ),
    ),
    (
        "f3_q4_collective_metric_origin_screen",
        (
            "LANE_GRA_FQ_F3_Q4_COLLECTIVE_METRIC_ORIGIN_SCREEN_V001/THEOREM.md",
            "07445c035ed4c5167a5a20280c4db69a5101eeb71831cdeb126b29702d04b69d",
        ),
        (
            "LANE_GRA_FQ_F3_Q4_COLLECTIVE_METRIC_ORIGIN_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md",
            "91aa35170432684a47278e46ee2b9d56658a43acc8acbbb84480d047cdbe6dcf",
        ),
    ),
    (
        "spag_public_data_substitute",
        (
            "LANE_GRA_SPAG_PUBLIC_DATA_SUBSTITUTE_V001/PUBLIC_DATA_SUBSTITUTE.md",
            "5b8455c062766acd25fb40f82ade6ded47f4b7c8443bcde178ba71cf2e451c4f",
        ),
        (
            "LANE_GRA_SPAG_PUBLIC_DATA_SUBSTITUTE_V001/AUDIT.md",
            "031af32f79fa93b15e885ba861b29b45500836cdb351b7dfc029070c491233ca",
        ),
    ),
    (
        "spag_public_data_second_pass",
        (
            "LANE_GRA_SPAG_PUBLIC_DATA_SECOND_PASS_V001/PUBLIC_DATA_SECOND_PASS.md",
            "3d4300b9c2998aab4a485771f097f860e570a3931b8a948be9e1b034925931a8",
        ),
        (
            "LANE_GRA_SPAG_PUBLIC_DATA_SECOND_PASS_V001/POST_INTEGRATION_CUSTODY_REAUDIT.md",
            "dcc170c2674fe0020f679e523ab40689be7efd32e05e71fe7600d9bbf4047e97",
        ),
    ),
    (
        "calibrated_finite_apparatus_g_crosscheck",
        (
            "LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/THEOREM.md",
            "cbf0733633ba93756b08dded7486a9be76beb572807693455c761bd36a8f0f5b",
        ),
        (
            "LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/HOSTILE_SELF_AUDIT.md",
            "d00b69cbc9bf2d352a003502976bf6ccb35af0176b71f70b58b0530b174e1315",
        ),
    ),
    (
        "nist_bipm_g_forward_readiness",
        (
            "LANE_GRA_FI_NIST_BIPM_G_FORWARD_READINESS_V001/READINESS.md",
            "824fd6ea9dc62e564f18875f90f460a1358b9d9acb84a99a6a04d984c6a6d0ef",
        ),
        (
            "LANE_GRA_FI_NIST_BIPM_G_FORWARD_READINESS_V001/INDEPENDENT_AUDIT.md",
            "6917298439d708931fa50201459659b10802aa13a36c6f39cf2dfdf80840b90e",
        ),
    ),
    (
        "hust_2018_dual_method_g_forward",
        (
            "LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/THEOREM.md",
            "44eb8c81a3d84dfa6829bcd6971d0261215877af0529318eab5cedbd3980c340",
        ),
        (
            "LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/INDEPENDENT_HOSTILE_AUDIT.md",
            "e8625b4cbf67d73927db495e8111e3ffbd4e85f46b80183e55fbf8b4391d0b2e",
        ),
    ),
    (
        "f3_q4_additive_block_strain_source_rank",
        (
            "LANE_GRA_FR_F3_Q4_BLOCK_STRAIN_SOURCE_RANK_OBSTRUCTION_V001/THEOREM.md",
            "62c7aaee9433a9ffa970ff6e38bac5585200cf40d6fca2cb70477e7e1e7524eb",
        ),
        (
            "LANE_GRA_FR_F3_Q4_BLOCK_STRAIN_SOURCE_RANK_OBSTRUCTION_V001/INDEPENDENT_HOSTILE_AUDIT.md",
            "d2da0796cfec7cff8f1d7da5c9bc449d38acdbae089dd9778fb5f19cb6e42b88",
        ),
    ),
    (
        "f3_q4_complete_reduced_source_rank",
        (
            "LANE_GRA_FS_F3_Q4_COMPLETE_SOURCE_RANK_AUDIT_V001/THEOREM.md",
            "36879f4c18eec83a22bdf9bd161d9d444b72e1dbda1d5eaa0312c6aab3d95724",
        ),
        (
            "LANE_GRA_FS_F3_Q4_COMPLETE_SOURCE_RANK_AUDIT_V001/INDEPENDENT_HOSTILE_AUDIT.md",
            "6e0ecd0febf6364e4122bbf2f65e1feb93c27e960bce30c0097ea0fbe3f58966",
        ),
    ),
    (
        "f3_q4_degree_pair_geometric_strain_boundary",
        (
            "LANE_GRA_FT_F3_Q4_DEGREE_PAIR_GEOMETRIC_STRAIN_BOUNDARY_V001/THEOREM.md",
            "6000f38871a57061b106665a41aca04b5d09f4c8c8f4bdc8132ccd5f3f1fbe39",
        ),
        (
            "LANE_GRA_FT_F3_Q4_DEGREE_PAIR_GEOMETRIC_STRAIN_BOUNDARY_V001/INDEPENDENT_HOSTILE_AUDIT.md",
            "599eec1cde6260be1c9f536274dd8682f77cb45d94e7e3cbc17a28d7552258bd",
        ),
    ),
    (
        "hust_tos_roundtrip_history_residual",
        (
            "LANE_GRA_HUST_TOS_ROUNDTRIP_HISTORY_RESIDUAL_V001/THEOREM.md",
            "8d989a539883274f3af460231bb60bd83355682a1435e57c9bfa0d35cf5ba5d2",
        ),
        (
            "LANE_GRA_HUST_TOS_ROUNDTRIP_HISTORY_RESIDUAL_V001/INDEPENDENT_HOSTILE_AUDIT.md",
            "931431da193ab3f546db25e7c355b2d22d397b8e85cd5a77c92d490545b4ce76",
        ),
    ),
    (
        "hust_nominal_source_kernel_reconstruction",
        (
            "LANE_GRA_HUST_NOMINAL_SOURCE_KERNEL_RECONSTRUCTION_V001/THEOREM.md",
            "cd729640537d66c52da0c9209fb94c1a95ff1a9dc7580ac4b6a9e4a7cea8e67c",
        ),
        (
            "LANE_GRA_HUST_NOMINAL_SOURCE_KERNEL_RECONSTRUCTION_V001/INDEPENDENT_HOSTILE_AUDIT.md",
            "e9669b5e4ef584db75574010dc0480df3afa944701b0bc14cff2ea9884d5fac7",
        ),
    ),
    (
        "clock_k5_common_potential_cycle_closure",
        (
            "LANE_GRA_CLOCK_K5_COMMON_POTENTIAL_CYCLE_CLOSURE_V001/THEOREM.md",
            "9599c3402710e390c03a661a8fbd22860d844f89a79fda017910fec9a021f232",
        ),
        (
            "LANE_GRA_CLOCK_K5_COMMON_POTENTIAL_CYCLE_CLOSURE_V001/INDEPENDENT_HOSTILE_AUDIT.md",
            "187114e102b92e319aa1ddda5571b012165170936bc994f39ecc9adf9cab9da6",
        ),
    ),
    (
        "hust_conditional_homogeneous_g_crosscheck",
        (
            "LANE_GRA_HUST_CONDITIONAL_HOMOGENEOUS_G_CROSSCHECK_V001/THEOREM.md",
            "3ddcdcc8d4ef9f905c9ff3e07e813efc2848317e0f2cde4141798b9143c0e3a8",
        ),
        (
            "LANE_GRA_HUST_CONDITIONAL_HOMOGENEOUS_G_CROSSCHECK_V001/INDEPENDENT_HOSTILE_AUDIT.md",
            "fcd2fa0fcc6dd70cbd1397f3043659627e561606649279deca5d2fb135eb4e9d",
        ),
    ),
)

_EXPECTED_ARTIFACTS = _CLOSURE_PAIR + tuple(
    artifact
    for _, theorem, audit in _CORE_SOURCE_PAIRS
    for artifact in (theorem, audit)
) + tuple(
    artifact
    for _, theorem, audit in _ADVANCE_SOURCE_PAIRS
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
                "q4_common_child_incidence_identity": "B_DAGGER_B_EQUALS_4I_PLUS_A_EXACT",
                "q4_A3_root_second_moment": "SUM_ALPHA_ALPHA_T_EQUALS_16I_OVER_3_EXACT",
                "q4_A3_cell_covolume": "16_A_STAR_CUBED_OVER_3_SQRT_3_EXACT",
                "q4_affine_refinement": "EXACT_MATHEMATICAL_ATLAS_NOT_PHYSICAL_MANIFOLD",
                "q4_acoustic_cone": "EXACT_ONLY_FOR_SEPARATELY_SUPPLIED_MASSLESS_ACTION",
                "q4_clifford_cone": "EXACT_PROSPECTIVE_STENCIL_NOT_CURRENT_F3",
                "q4_diamond_support_completion": "EXACT_TRANSLATION_COMPLETION_AND_DEEP_LOCAL_EXHAUSTION",
                "q4_diamond_u1_inheritance": (
                    "FINITE_LOCAL_EDGE_BINDING_PROGRAMMABLE__GLOBAL_U1_CONDITIONAL_ON_"
                    "SUPPLIED_REGULAR_COMPLETION_AND_LEADING_ORDER"
                ),
                "q4_finite_programmed_support_solder": (
                    "EXACT_FOR_SUPPLIED_FINITE_HARDWARE_ADDRESS_EDGE_LIST_CAP_SCHEDULE_AND_PORTS"
                ),
                "q4_finite_programmed_hold_condition": (
                    "RAW_AND_K_GATED_INCIDENCE_FLIPS_EXACTLY_OFF_OR_CONTINUOUSLY_CANCELLED"
                ),
                "q4_programmed_floquet_detuning": (
                    "EXACT_CHILD_ONLY_PHASE_AND_UNIFORM_PARENT_CHILD_QUASIENERGY_SEPARATION_"
                    "FOR_SUPPLIED_REPEATABLE_DUAL_FLIP_FREE_SCHEDULE"
                ),
                "q4_programmed_floquet_parent_branch": (
                    "EXACT_FUNCTION_OF_B_DAGGER_B_WITH_CONTROLLED_SIBLING_KERNEL"
                ),
                "q4_static_source_off_child_parent_stagger": "ABSENT",
                "q4_raw_finite_slab_global_d2_ice": "EXACTLY_EMPTY",
                "q4_ice_local_diagonal_module": "A1_PLUS_E_PLUS_T2_EXACT",
                "q4_ice_one_link_module": "T2_RANK_3_EXACT",
                "q4_ice_pair_module": "A1_PLUS_E_RANK_3__CENTERED_TANGENT_E_RANK_2",
                "q4_ice_hybrid_tensor_representation": (
                    "INDEPENDENT_A1_PLUS_PAIR_E_PLUS_ONE_LINK_T2_EXACT_ISOMORPHISM_CANDIDATE"
                ),
                "q4_ice_symmetric_fisher_T2_first_derivative": "EXACTLY_ZERO_BY_COMPLEMENT_PARITY",
                "q4_ice_ring_E_T2_dynamics": "EXACT_FINITE_COMMUTATORS_AND_LINKED_RING_RESPONSE",
                "q4_ice_global_ring_domain": "REQUIRES_COMPATIBLE_BOUNDARY_PERIODIC_OR_CONTROLLED_INFINITE_COMPLETION",
                "q4_gaussian_maxwell_pole_screen": (
                    "EXACT_CONDITIONAL_SPIN1_LINK_POLE_PAIR_CONTINUUM_NO_ISOLATED_HELICITY2"
                ),
                "q4_ice_fisher_covariance": "F_EQUALS_4_DIAG_W_MINUS_M_M_TRANSPOSE_EXACT",
                "q4_ice_odd_T2_fisher_second_jet": (
                    "S4_ISOMORPHISM_SYM2_T2_TO_A1_E_T2__NOT_O3_EQUIVARIANT"
                ),
                "q4_ice_complement_preserving_fisher_tangent": "E_ONLY__A1_PLUS_E_WITH_SEPARATE_SCALAR",
                "q4_ice_generic_broken_background_covariance_map_with_separate_scalar_rank": 6,
                "q4_ice_broken_background_T2_origin": "VECTOR_MEAN_DYAD_NOT_INDEPENDENT_TENSOR",
                "q4_ice_unlabelled_complement_T2": "CANCELS_EXACTLY",
                "q4_ice_complement_symmetrized_conditional_T2": (
                    "RETENTION_REQUIRES_CONTROL_SIGN_RECORD"
                ),
                "q4_ice_broken_background_preparation_and_stabilization": "NOT_DERIVED",
                "q4_ice_order6_interaction": "COMPACT_HARDCORE_RING_ALREADY_NON_GAUSSIAN",
                "q4_ice_order8_scope": (
                    "SUPPLIED_FINITE_SIMPLE_Z4_BIPARTITE_GIRTH_GE6_SUPPORT__DSTAR2_ER0_"
                    "SYMMETRIC_DETUNING__FIXED_FESHBACH_CONVENTION"
                ),
                "q4_ice_order8_endpoint_topologies": (
                    "SCALAR_DIAGONAL_PLUS_DRESSED_HEXAGON_PLUS_NEW_ALTERNATING_OCTAGON"
                ),
                "q4_ice_order8_octagon_coefficient": "J8_EQUALS_429_H8_OVER_16_UD7",
                "q4_ice_order8_diagonal_potential": "V8_EQUALS_ZERO",
                "q4_finite_tt_periodic_quotient": (
                    "30_CELLS_60_VERTICES_120_LINKS_120_ELEMENTARY_HEXAGONS"
                ),
                "q4_finite_tt_ring_sector": (
                    "180_STATES_420_TRANSITIONS_TRANSLATION_CLOSED_H6_ONLY"
                ),
                "q4_finite_tt_composite_cumulant": (
                    "W2_1P130847135995723_OVER_J6__W4_MINUS_0P136825085605100_OVER_J6_CUBED__"
                    "GAMMA4_COMP_0P083666214307836_J6"
                ),
                "q4_finite_tt_lowest_pole_vs_two_link_proxy": (
                    "3P194109035554332_J6_ABOVE_2P059674505691458_J6"
                ),
                "q4_projected_ice_constraint_species": "ONE_SCALAR_U1_GAUSS_ONLY",
                "q4_pair_relations_constraint_status": (
                    "ALGEBRAIC_ZERO_OPERATORS_NOT_NEW_FIRST_CLASS_GENERATORS"
                ),
                "q4_A3_static_cometric_tangent_rank": 6,
                "q4_inherited_even_bulk_Kubo_dynamic_rank": 0,
                "q4_fixed_parent_collective_metric_origin": (
                    "NO_PRESENT_OBJECT_JOINTLY_OWNS_SIX_CONFIG_CHANNELS_"
                    "INDEPENDENT_CONJUGATES_AND_VECTOR3_PLUS_SCALAR_NULL_PACKET"
                ),
                "q4_additive_block_strain_source_rank": 4,
                "q4_additive_block_strain_source_sectors": "A1_PLUS_T2",
                "q4_additive_block_strain_E_null_dimension": 2,
                "q4_complete_reduced_microscopic_source_rank": 4,
                "q4_complete_reduced_effective_rank_through_order8": "AT_MOST_4",
                "q4_same_H0_determines_source_derivative": False,
                "q4_DPAR_conditional_microscopic_source_rank": 6,
                "q4_DPAR_direct_ice_projected_image": "RANK3_A1_PLUS_E",
                "q4_DPAR_local_H6_E_response": "EXACT_NONZERO_ON_EXPLICIT_PERIODIC_WITNESSES",
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
                "finite_physical_link_pair_response": (
                    "EXACT_A1_E_T2_RESPONSE_ON_UNCONSTRAINED_FOUR_LINK_BLOCK"
                ),
                "finite_adjacent_cell_pair_kernel": "EXACT_NONZERO_SHARED_LINK_RESPONSE",
                "ice_projected_pair_tangent": "E_ONLY_RANK_2",
                "ice_hybrid_symmetric_tensor_representation": "A1_PLUS_E_PLUS_T2_CANDIDATE_EXACT",
                "ice_physical_metric_T2_linear_response": "OPEN__DIRECT_SYMMETRIC_QUERY_ZERO",
                "gaussian_maxwell_one_link_pole": "SPIN_1_PHOTON_CONDITIONAL_ON_MAXWELL_IR",
                "gaussian_maxwell_centered_pair_spectrum": (
                    "TWO_PHOTON_CONTINUUM_PLUS_CONTACTS_CONDITIONAL_ON_MAXWELL_IR"
                ),
                "gaussian_maxwell_isolated_helicity2_pole": "ABSENT_IN_DIRECT_COMPOSITE_ROUTE",
                "finite_group_S4_label_determines_continuum_spin": False,
                "ice_local_fisher_metric_full_linear_tensor_tangent": "ABSENT_WITH_COMPLEMENT_PRESERVED",
                "ice_broken_background_full_rank_kind": "BACKGROUND_CONTROL_TO_COVARIANCE_MAP_NOT_SIX_MODE_FISHER_METRIC",
                "inherited_TT_bare_vertex": (
                    "FINITE_RANGE_ANALYTIC__TT_PROJECTION_ALLOWED_NOT_EVALUATED"
                ),
                "inherited_TT_strict_single_insertion_new_pole": "NOT_ESTABLISHED",
                "inherited_TT_decisive_sequence": (
                    "CONNECTED_FOUR_POINT_THEN_AMPUTATION_THEN_CHANNEL_2PI_THEN_"
                    "BETHE_SALPETER_OR_SPECTRAL_THEN_FINITE_VOLUME_HELICITY_WARD_"
                    "COMMON_CONE_RESIDUE"
                ),
                "finite_TT_composite_cumulant_kind": (
                    "FOUR_Q_EIGHT_ONE_LINK_NOT_CONNECTED_FOUR_ONE_LINK_OR_PHOTON_AMPUTATED"
                ),
                "finite_TT_composite_legendre_quartic": (
                    "POSITIVE_IN_BOTH_SUSCEPTIBILITY_EIGENCHANNELS"
                ),
                "finite_TT_pole_screen": (
                    "NO_BELOW_PROXY_OR_ENERGY_EXCLUSIVE_CANDIDATE__NOT_A_NO_BOUND_STATE_THEOREM"
                ),
                "projected_ice_constraint_packet": "ONE_SCALAR_U1_NOT_RGRLB",
                "projected_ice_equal_two_polarizations_implies_RGRLB": False,
                "A3_static_metric_span": "RANK6_COEFFICIENT_TANGENT",
                "A3_inherited_even_source_Kubo_rank": 0,
                "FJ_unprojected_pair_response": (
                    "EXACT_CONDITIONAL_RANK6_A1_E_T2_RESPONSE_WITH_NEAREST_CELL_SPREADING"
                ),
                "current_fixed_parent_metric_packet": (
                    "NO_JOINT_SIX_CONFIG_PLUS_CONJUGATE_PLUS_VECTOR3_SCALAR_NULL_OBJECT_"
                    "IN_CONSTRUCTED_CATALOG__NOT_THERMODYNAMIC_NO_GO"
                ),
                "additive_block_strain_source_rank": 4,
                "additive_block_strain_source_null": "E_RANK2",
                "complete_reduced_CWFM_source_rank": 4,
                "complete_reduced_effective_rank_through_order8": "AT_MOST_4",
                "same_H0_source_derivative_uniqueness": "EXACTLY_FALSE",
                "DPAR_conditional_microscopic_source_rank": 6,
                "DPAR_direct_ice_projected_image": "RANK3_A1_PLUS_E",
                "DPAR_full_CTP_rank": "OPEN_NOT_COMPUTED",
                "next_collective_metric_calculation": (
                    "DERIVE_DPAR_PHYSICAL_REALIZATION_THEN_FULL_CTP_WARD_RESPONSE"
                ),
            },
            "microscopic_parent_boundary": {
                "scalar_carrier_transfer_on_supplied_saturated_q4_support": (
                    "EXACT_FROM_UNCHANGED_F3_ONE_CARRIER_RESTRICTION"
                ),
                "q4_label_to_coexisting_F3_site_edge_solder": "NOT_DERIVED_BY_CURRENT_PARENT",
                "finite_programmed_site_edge_solder": (
                    "EXACT_REVERSIBLE_FIXED_ORTHOGONAL_PROGRAM_WITH_SUPPLIED_PHYSICAL_ANTECEDENTS"
                ),
                "finite_programmed_solder_is_autonomous_support_selection": False,
                "finite_programmed_solder_is_scalable_thermodynamic_instantiation": False,
                "positive_child_parent_carrier_detuning": "NOT_OWNED_BY_CURRENT_STATIC_SOURCE_OFF_PARENT",
                "positive_programmed_floquet_carrier_detuning": (
                    "EXACT_ON_SUPPLIED_REPEATABLE_DUAL_FLIP_FREE_CONTROLLER_SCHEDULE"
                ),
                "programmed_floquet_detuning_is_autonomous_phase": False,
                "full_hopping_and_d2_ice_same_n_coexistence": "EXACTLY_INCOMPATIBLE",
                "six_static_pair_registers": "EXACT_FINITE_S4_REPRESENTATION",
                "inherited_pair_projector_dynamics": "CONSERVE_EVERY_PAIR_PROJECTOR",
                "interpair_retarded_kernel": "EXACTLY_ZERO_UNDER_INHERITED_DYNAMICS",
                "physical_q4_link_walsh_pairs": (
                    "EXACT_OPERATOR_REALIZATION_ON_FPMH_QUALIFIED_FINITE_PROGRAMMED_LINK_FACTORS"
                ),
                "physical_q4_link_pair_response": (
                    "EXACT_FINITE_LOCAL_AND_SHARED_LINK_RESPONSE_WITH_OPERATOR_SPREADING"
                ),
                "physical_q4_link_pairs_are_automatically_PMMDC_records": False,
                "q4_pair_operator_realization": (
                    "EXACT_FINITE_WALSH_ALGEBRA__FULL_PMMDC_AND_METRIC_SOLDER_OPEN"
                ),
                "noncommuting_pair_field_dynamics": (
                    "EXACT_FINITE_LINK_AND_ICE_RING_DYNAMICS__THERMODYNAMIC_TENSOR_FIELD_OPEN"
                ),
                "ice_hybrid_tensor_solder": (
                    "EXACT_REPRESENTATION_CANDIDATE__PHYSICAL_METRIC_T2_RESPONSE_AND_CALIBRATION_OPEN"
                ),
                "ice_ring_response_pole": "FINITE_2J6_NOT_MASSLESS_TENSOR",
                "direct_gaussian_composite_tensor_route": "EXACT_CONDITIONAL_NO_GO",
                "next_tensor_route": "INHERITED_NON_GAUSSIAN_TT_KERNEL_OR_DISTINCT_RANK2_CONSTRAINED_PHASE",
                "ice_T2_fisher_solder_boundary": (
                    "SECOND_JET_OR_VECTOR_BACKGROUND_ONLY__INDEPENDENT_LINEAR_TENSOR_SOLDER_OPEN"
                ),
                "inherited_order8_operator_boundary": (
                    "EXACT_J8_AND_V8_ZERO__HEXAGON_DRESSING_COEFFICIENTS_TYPED_NOT_ALL_REDUCED"
                ),
                "finite_TT_composite_precursor": (
                    "EXACT_STATIC_TWO_Q_AND_FOUR_Q_CUMULANTS_AND_FINITE_SPECTRUM"
                ),
                "connected_four_one_link_channel_2PI": "OPEN_NOT_COMPUTED_BY_FO",
                "inherited_projected_constraint_architecture": "SCALAR_U1_ONLY",
                "microscopic_RGRLB_from_current_q4_ice": "NOT_DERIVED",
                "six_A3_static_deformation_coefficients": "EXACT_RANK6",
                "six_A3_dynamic_metric_fields": "ABSENT_IN_INHERITED_EVEN_BULK_SCREEN",
                "q4_additive_block_strain_source_status": (
                    "EXECUTED_EXACT_RANK4_A1_T2_WITH_E_NULL2"
                ),
                "q4_complete_reduced_source_status": (
                    "EXECUTED_EXACT_MICRO_RANK4_EFFECTIVE_RANK_AT_MOST4_THROUGH_ORDER8"
                ),
                "unreduced_physical_BS_source": "UNDERDETERMINED_BY_CURRENT_PARENT",
                "degree_pair_E_query": "AVAILABLE_AND_LOCALLY_H6_DYNAMICAL",
                "DPAR_status": "SUFFICIENT_CONDITIONAL_LAW__NEITHER_INHERITED_NOR_ADOPTED",
                "DPAR_conditional_microscopic_rank": 6,
                "DPAR_direct_ice_projected_rank": 3,
                "full_state_dependent_CTP_Ward_response": "OPEN_NOT_COMPUTED",
                "new_interaction_or_second_field_adopted": False,
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
                "finite_apparatus_forward_model": "EXACT_DECLARED_EXTENDED_SOURCE_AND_TWO_MODE_SCHUR_RESPONSE",
                "same_source_remainder_data_response": "EXACTLY_ZERO",
                "calibrated_nonsingular_row_identifies": "P_EQUALS_G_TIMES_SOURCE_SCALE",
                "free_global_source_scale_degeneracy": "F_G_S_EQUALS_F_GQ_S_OVER_Q",
                "independent_source_interval_maps_to_G_interval": "G_IN_P_MINUS_OVER_S_PLUS_TO_P_PLUS_OVER_S_MINUS",
                "synthetic_validation": "PASS_15_OF_15_WITH_ARBITRARY_NONEMPIRICAL_G",
                "nist_bipm_public_reduced_forward": "DELTA_N_J_EQUALS_G_A_J_PLUS_R_J_EXACT",
                "nist_bipm_analysis_observations": 8,
                "nist_bipm_nominal_source_stiffness": "ZERO_AT_TORQUE_EXTREMA",
                "nist_bipm_G_only_jacobian_rank": 1,
                "nist_bipm_full_GC16_fit_ready": False,
                "nist_bipm_public_missing_field_count": 10,
                "nist_bipm_independent_G_crosscheck_performed": False,
                "hust_public_forward_status": (
                    "PROCESSED_DUAL_CHANNEL_FORWARD_CLOSED__FULL_GC16_NOT_READY"
                ),
                "hust_tos_A_B_A_delta_omega2_s_minus2": "1.6626945111323172e-6",
                "hust_tos_quadratic_delta_omega2_s_minus2": "1.6626989120180067e-6",
                "hust_aaf_processed_forward_count": 3,
                "hust_aaf_max_rounding_difference_ppm": "LESS_THAN_0P2",
                "hust_cross_method_separation_ppm": "44.9483_DESCRIPTIVE_ONLY",
                "hust_cross_method_z_condition": "2P7196_ONLY_IF_CROSS_COVARIANCE_ZERO",
                "hust_accepted_G_input_used": False,
                "hust_full_GC16_fit_ready": False,
                "hust_independent_G_crosscheck_performed": False,
                "hust_roundtrip_history_return_rank": 36,
                "hust_roundtrip_history_differential_rank": 18,
                "hust_roundtrip_history_status": (
                    "DESCRIPTIVE_HISTORY_CONFOUND__NO_MATCHED_NO_EXCURSION_ARM"
                ),
                "hust_nominal_source_kernel_status": (
                    "CONDITIONAL_HOMOGENEOUS_FUNCTIONALS_RECONSTRUCTED__FULL_MASS_STRESS_REMAINDER_OPEN"
                ),
                "hust_nominal_AAF_kernel_kg_m_minus3": "6926P660438859097_TO_6926P700007763433",
                "hust_nominal_ToS_kernel_kg_m_minus3": "24914_TO_25005_APPROXIMATE",
                "clock_K5_cut_rank": 4,
                "clock_K5_cycle_dimension": 6,
                "clock_K5_marginal_box_rho_star": "27_OVER_82",
                "clock_K5_status": "COMMON_NODE_SCALAR_COMPATIBLE__NOT_INDEPENDENT_GRAVITY_TEST",
                "hust_conditional_AAF_G_SI": (
                    "6P674235591785795E_MINUS11",
                    "6P674022699178638E_MINUS11",
                    "6P674260454934994E_MINUS11",
                ),
                "hust_conditional_ToS_zero_correction_anchor_range_SI": (
                    "6P673451644467455E_MINUS11_TO_6P673682353003419E_MINUS11"
                ),
                "hust_conditional_G_source_premise": "R_NORM_EQUALS_ZERO_HOMOGENEOUS_APPARATUS",
                "hust_conditional_G_accepted_or_CODATA_input_used": False,
                "hust_public_physical_G_point_or_compact_interval_identified": False,
                "hust_AAF_normalization_collision_ppm": "1631P537953",
                "hust_ToS_normalization_collision_ppm": "152P258414_TO_153P469498",
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
                "public_data_substitute_executed": True,
                "public_packets_admitted": ("PAGE_GEILKER", "FUCHS", "NIST_BIPM_2026", "PANDA_HOLDOUT_PRESERVED"),
                "public_same_parent_eight_cell_support": "ABSENT",
                "public_beta_TM_identifiable": False,
                "page_geilker_generous_proxy_rank": "2_OF_8",
                "nist_bipm_planning_envelope_nN_m": "0.00019286904345467502_TO_0.000578607130364025",
                "public_result_is_lineage_null": False,
                "public_second_pass_executed": True,
                "public_second_pass_frozen_query_count": 28,
                "public_second_pass_new_lineage_roots": 0,
                "public_second_pass_component_roots": 2,
                "public_second_pass_result_is_exhaustive_world_search": False,
                "panda_response_holdout_opened": False,
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
                "public_same_parent_lineage_estimand": "NOT_IDENTIFIED_IN_ADMITTED_DATA",
                "finite_apparatus_G_measurement": "NOT_PERFORMED",
                "autonomous_q4_support_selection": "OPEN",
                "autonomous_q4_detuning_phase": "OPEN",
                "q4_pair_field_lift": "PARTIAL_FINITE_OPERATOR_AND_RESPONSE_REALIZATION",
                "physical_metric_pair_solder": "OPEN",
                "thermodynamic_massless_tensor_phase": "OPEN",
                "compatible_global_ice_completion": "SUPPLIED_NOT_DERIVED",
                "direct_gaussian_composite_helicity2_route": "CLOSED_NEGATIVE",
                "inherited_non_gaussian_TT_kernel": (
                    "FINITE_COMPOSITE_CUMULANT_AND_SPECTRUM_COMPUTED__"
                    "CONNECTED_FOUR_ONE_LINK_CHANNEL_2PI_OPEN"
                ),
                "normalized_connected_TT_four_point": (
                    "OPEN__FO_COMPUTED_FOUR_Q_COMPOSITE_NOT_THIS_OBJECT"
                ),
                "nonperturbative_TT_pole_Ward_common_cone_residue": "OPEN",
                "local_ice_fisher_T2_solder": (
                    "CLOSED_NEGATIVE_FOR_LOCAL_DIAGONAL_COMPLEMENT_PRESERVING_FAMILIES"
                ),
                "independent_even_rank2_collective_variable": "OPEN",
                "current_q4_ice_RGRLB_constraint_origin": (
                    "CLOSED_NEGATIVE_FOR_CURRENT_FINITE_LOCAL_Q4_ICE_BRANCH__"
                    "THERMODYNAMIC_COLLECTIVE_ORIGIN_OPEN"
                ),
                "current_fixed_parent_collective_metric_origin": (
                    "CLOSED_NEGATIVE_FOR_CONSTRUCTED_OBJECT_CATALOG__NOT_THERMODYNAMIC_NO_GO"
                ),
                "q4_source_rank_boundary": (
                    "ADDITIVE_AND_COMPLETE_REDUCED_SOURCE_RANK4_WITH_E_NULL2__"
                    "DPAR_CONDITIONALLY_REPAIRS_MICROSCOPIC_RANK_TO6"
                ),
                "next_no_lab_gravity_calculation": (
                    "DERIVE_DPAR_PHYSICAL_REALIZATION_OR_RETAIN_EXPLICIT_CONDITIONAL_PREMISE_"
                    "THEN_TEST_FULL_CTP_WARD_RESPONSE"
                ),
                "hust_public_G_forward": (
                    "PROCESSED_DUAL_CHANNEL_CLOSED__FULL_GC16_NOT_READY__NO_NEW_G"
                ),
                "hust_public_source_reconstruction": (
                    "CONDITIONAL_HOMOGENEOUS_KERNELS_COMPUTED__MATCHED_FULL_NUMERATOR_OPEN"
                ),
                "hust_public_conditional_G_quotient": (
                    "THREE_AAF_QUOTIENTS_AND_SEVEN_TOS_AFFINE_FAMILIES_COMPUTED__"
                    "NO_PUBLIC_PHYSICAL_G_POINT_OR_COMPACT_INTERVAL__NO_NEW_G"
                ),
                "public_history_and_clock_diagnostics": (
                    "HISTORY_CONFOUND_MAPPED_AND_K5_NODE_SCALAR_COMPATIBILITY_PASSED__"
                    "NO_LINEAGE_OR_INDEPENDENT_GRAVITY_RESULT"
                ),
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
                "public_data_lineage_null_authorized": False,
                "finite_apparatus_G_measurement_authorized": False,
                "autonomous_q4_support_or_full_pair_field_lift_derived": False,
                "programmed_floquet_detuning_derived": True,
                "autonomous_q4_detuning_phase_derived": False,
                "finite_physical_link_pair_response_derived": True,
                "physical_metric_pair_solder_derived": False,
                "ice_hybrid_tensor_representation_derived": True,
                "massless_tensor_from_ice_ring_response_derived": False,
                "gaussian_composite_helicity2_pole_derived": False,
                "complement_preserving_local_fisher_T2_solder_derived": False,
                "generic_broken_background_covariance_rank6_with_separate_scalar_derived": True,
                "order8_inherited_loop_operator_boundary_derived": True,
                "finite_TT_composite_cumulant_and_spectrum_derived": True,
                "connected_four_one_link_channel_2PI_derived": False,
                "q4_scalar_U1_constraint_only_derived": True,
                "microscopic_RGRLB_from_current_q4_ice_derived": False,
                "six_static_cometric_deformations_derived": True,
                "six_dynamic_metric_fields_from_inherited_even_bulk_derived": False,
                "q4_block_strain_source_rank_screen_executed": True,
                "q4_complete_reduced_source_rank_screen_executed": True,
                "q4_degree_pair_DPAR_boundary_executed": True,
                "q4_DPAR_adopted": False,
                "q4_DPAR_inherited_from_current_parent": False,
                "q4_DPAR_conditional_microscopic_rank6_derived": True,
                "q4_physical_metric_source_derived": False,
                "q4_full_CTP_Ward_response_executed": False,
                "hust_processed_dual_channel_forward_executed": True,
                "hust_roundtrip_history_diagnostic_executed": True,
                "hust_nominal_source_kernel_reconstruction_executed": True,
                "hust_conditional_homogeneous_G_crosscheck_executed": True,
                "clock_K5_common_node_compatibility_executed": True,
                "hust_full_GC16_executed": False,
                "inherited_protected_tensor_pole_derived": False,
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
                "mathematical_q4_refinement_promotes_physical_spacetime": False,
                "supplied_massless_action_promotes_current_F3_massless_phase": False,
                "q4_diamond_shape_promotes_autonomous_support_or_gravity": False,
                "static_pair_register_promotes_propagating_pair_field": False,
                "public_data_design_ceiling_promotes_empirical_null": False,
                "synthetic_finite_apparatus_G_fit_promotes_measurement": False,
                "programmed_floquet_gap_promotes_static_or_autonomous_detuning": False,
                "finite_link_pair_response_promotes_PMMDC_solder_or_metric": False,
                "ice_representation_isomorphism_promotes_physical_metric": False,
                "finite_ring_response_promotes_massless_tensor": False,
                "spin1_photon_or_two_photon_continuum_promotes_graviton": False,
                "second_fisher_jet_promotes_linear_metric_tangent": False,
                "broken_vector_background_rank_promotes_tensor_mode": False,
                "control_sign_conditioning_promotes_endogenous_metric_solder": False,
                "bare_vertex_analyticity_promotes_dressed_1PI_analyticity": False,
                "strict_single_insertion_boundary_promotes_nonperturbative_no_pole": False,
                "four_Q_composite_cumulant_promotes_four_link_1PI_or_binding": False,
                "finite_threshold_proxy_promotes_thermodynamic_no_bound_state": False,
                "equal_two_polarizations_promotes_RGRLB": False,
                "six_static_deformation_coefficients_promote_dynamic_metric_fields": False,
                "current_catalog_obstruction_promotes_thermodynamic_no_go": False,
                "public_component_dataset_promotes_lineage_estimand": False,
                "nist_summary_reduction_promotes_independent_G_measurement": False,
                "hust_processed_forward_promotes_new_G_or_GFT_confirmation": False,
                "rank4_additive_source_promotes_all_physical_source_no_go": False,
                "DPAR_conditional_rank6_promotes_inherited_physical_metric_source": False,
                "HUST_history_residual_promotes_lineage_gravity": False,
                "HUST_nominal_kernel_promotes_full_apparatus_G": False,
                "HUST_conditional_homogeneous_quotient_promotes_new_or_independent_G": False,
                "clock_K5_compatibility_promotes_independent_GR_or_common_metric": False,
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
                "advance_source_pairs": tuple(
                    _pair_certificate(custody, label, theorem, audit)
                    for label, theorem, audit in _ADVANCE_SOURCE_PAIRS
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
