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


SCHEMA = "WAC_GRAVITY_FORMATION_THEORY_CERTIFICATE_V013"
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
    (
        "f3_q4_pair_resolved_coulomb_dpar_realization",
        (
            "LANE_GRA_FU_F3_Q4_PAIR_RESOLVED_MAXWELL_DPAR_DERIVATION_V001/THEOREM.md",
            "f088346f72861b3b11ae737fe6b882d43da9e747fc1d1d1f6bd446a7fd2b6272",
        ),
        (
            "LANE_GRA_FU_F3_Q4_PAIR_RESOLVED_MAXWELL_DPAR_DERIVATION_V001/INDEPENDENT_HOSTILE_AUDIT.md",
            "6f859779566177b5999cfe02c01cd569c5bd7b0b4ec2b21b0b3e79ebf26f9277",
        ),
    ),
    (
        "f3_q4_coulomb_projected_source_rank",
        (
            "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/THEOREM.md",
            "6fc221a31151340b91a946d33e442971c1373500e067c354b6c610e3964edb1c",
        ),
        (
            "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md",
            "3801fd9ba6ba3c0fe80c9f4792abfdeb6dd7c37c7145663be05b4d56f8160723",
        ),
    ),
    (
        "f3_q4_coulomb_fv_witness_finite_response",
        (
            "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/THEOREM.md",
            "db3e12d50fd1cb41cddc722a0445cdeaef6a52d49704fa6df1028dfd9abcba1b",
        ),
        (
            "LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md",
            "77ecf5f906c09fc171c20e8b1c431ac76a8db7d19369a74506d49877428ba463",
        ),
    ),
    (
        "bense_public_memory_history_frontier",
        (
            "LANE_GRA_BENSE_PUBLIC_MEMORY_HISTORY_FRONTIER_V001/THEOREM.md",
            "9eef510038cf541f70b22082613231024e42d81e17818fb0cf8b866df540f75c",
        ),
        (
            "LANE_GRA_BENSE_PUBLIC_MEMORY_HISTORY_FRONTIER_V001/INDEPENDENT_HOSTILE_AUDIT.md",
            "d4592d19ca0219a601b45074d30fcae6c9d9319a01881d27e98d523b617e3689",
        ),
    ),
    (
        "hust_public_calibrated_source_identifiability",
        (
            "LANE_GRA_HUST_PUBLIC_CALIBRATED_SOURCE_IDENTIFIABILITY_V001/THEOREM.md",
            "ce2177498cdae8a4fe611deb8a84526554177bd33b9fa59eff48f7b5f7391aa2",
        ),
        (
            "LANE_GRA_HUST_PUBLIC_CALIBRATED_SOURCE_IDENTIFIABILITY_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md",
            "95693c55b483243034358d29dfe6dbe0b4ac4a4e96a906a0c0e72c4b1db68232",
        ),
    ),
    (
        "hust_public_bounded_completion_search",
        (
            "LANE_GRA_HUST_PUBLIC_COMPLETION_SEARCH_V001/THEOREM.md",
            "b763e8d313aa034368f297ffd49bbad016209f4ad70f84a33082c3b049ec733a",
        ),
        (
            "LANE_GRA_HUST_PUBLIC_COMPLETION_SEARCH_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md",
            "faf1b012878b3b4c7daa335a9a9520e23eaa7cc3688f53d7e2a950a7c971c02c",
        ),
    ),
    (
        "f3_q4_coulomb_complete_h6_homogeneous_response",
        (
            "LANE_GRA_FX_F3_Q4_COULOMB_COMPLETE_H6_SOURCE_RESPONSE_V001/THEOREM.md",
            "2bf65e602dfbb5cf8cad7b69d5f22aa8ae01904924e320006322e251fc9ca5a4",
        ),
        (
            "LANE_GRA_FX_F3_Q4_COULOMB_COMPLETE_H6_SOURCE_RESPONSE_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md",
            "57ccae45274bdb6f118e8bf98c4be1b2cea279658f96ad004380b19383a4dac8",
        ),
    ),
    (
        "f3_q4_native_support_m1_complete_h6_response",
        (
            "LANE_GRA_FY_F3_Q4_NATIVE_SUPPORT_M1_COMPLETE_H6_RESPONSE_V001/THEOREM.md",
            "8db3dd16c36e0205b5c98fc3154e8a2f1876d243c3c1d2068424c1276ee68f28",
        ),
        (
            "LANE_GRA_FY_F3_Q4_NATIVE_SUPPORT_M1_COMPLETE_H6_RESPONSE_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md",
            "f9baf2b21ead24d947866192eb7f0cb6d4e353ffb2b1107569cc442564804f21",
        ),
    ),
    (
        "f3_q4_m1_continuity_contact_ward_boundary",
        (
            "LANE_GRA_FZ_F3_Q4_M1_CONTINUITY_CONTACT_WARD_BOUNDARY_V001/THEOREM.md",
            "c0750b7d8a6a7f1b12d3ef76e8d5a6a3754a86e714f75b7efc203a56c7cfeaf9",
        ),
        (
            "LANE_GRA_FZ_F3_Q4_M1_CONTINUITY_CONTACT_WARD_BOUNDARY_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md",
            "a6ad33451d21f94a530a8b06a89b12c4e1a085de284e353fdaa9893f5748b63e",
        ),
    ),
    (
        "f3_q4_fu09b_encoded_charge_current_lift",
        (
            "LANE_GRA_GA_F3_Q4_FU09B_ENCODED_CHARGE_CURRENT_LIFT_V001/THEOREM.md",
            "374168c75f928cdbce55ca790157ba00df2d3a7f2045611d805d1b1b3c090336",
        ),
        (
            "LANE_GRA_GA_F3_Q4_FU09B_ENCODED_CHARGE_CURRENT_LIFT_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md",
            "651ba2ef0cfba5d87faaa370030df3efca8a08413931236f1142b72289eaa83b",
        ),
    ),
    (
        "f3_q4_fixed_support_energy_momentum_ward_boundary",
        (
            "LANE_GRA_GB_F3_Q4_FIXED_SUPPORT_ENERGY_MOMENTUM_WARD_BOUNDARY_V001/THEOREM.md",
            "60555987258ac723ac86a1f879dbea3222fc60496350390d880de5d346dff1fa",
        ),
        (
            "LANE_GRA_GB_F3_Q4_FIXED_SUPPORT_ENERGY_MOMENTUM_WARD_BOUNDARY_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md",
            "ccedbf8067d15cb716a6871a1f5e07c3098fecad9cd84e9e0b6c568cdbee388d",
        ),
    ),
    (
        "f3_q4_finite_family_common_cone_boundary",
        (
            "LANE_GRA_GC_F3_Q4_FINITE_FAMILY_COMMON_CONE_BOUNDARY_V001/THEOREM.md",
            "f74b8e5e8f28643a207e0ad26b378c1989d1fbe8007cca2044661951395dc98f",
        ),
        (
            "LANE_GRA_GC_F3_Q4_FINITE_FAMILY_COMMON_CONE_BOUNDARY_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md",
            "2588d4b11bf83b6db78e0ee86a1d9edb6767063ae2735f4b2a3b2a3a6b771585",
        ),
    ),
    (
        "f3_q4_translation_owning_recoil_parent",
        (
            "LANE_GRA_GD_F3_Q4_TRANSLATION_OWNING_RECOIL_PARENT_V001/THEOREM.md",
            "c1fde7eca80d8d555e42556b3573c0444e6ac9ce5001a26b15c0eb04e125a4a6",
        ),
        (
            "LANE_GRA_GD_F3_Q4_TRANSLATION_OWNING_RECOIL_PARENT_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md",
            "d1d74fe66e5f2d888a18348e2e347a8abe32260ba50f94562422a8a688166341",
        ),
    ),
    (
        "f3_q4_public_spectral_payload_search",
        (
            "LANE_GRA_GD_F3_Q4_PUBLIC_SPECTRAL_PAYLOAD_SEARCH_V001/RESULT.md",
            "a750daa417a1ffc70aed572fb3e4127617904dcad2c682df8cf1ea5cc65c07d4",
        ),
        (
            "LANE_GRA_GD_F3_Q4_PUBLIC_SPECTRAL_PAYLOAD_SEARCH_V001/VERIFICATION.txt",
            "c7d4af7e5454ba17ba2430fa59f018fc4402196699181599700d9983b0c5dfa3",
        ),
    ),
    (
        "f3_q4_canonical_ir_ancestry_observable_contract",
        (
            "LANE_GRA_GF_F3_Q4_CANONICAL_IR_ANCESTRY_OBSERVABLE_V001/THEOREM.md",
            "928541b88f2c306905d715ea1e7b7ea0ab0f59d405134fe768f0d3f126c62e91",
        ),
        (
            "LANE_GRA_GF_F3_Q4_CANONICAL_IR_ANCESTRY_OBSERVABLE_V001/INDEPENDENT_HOSTILE_REAUDIT_V005/INDEPENDENT_HOSTILE_REAUDIT_V005.md",
            "63fdc94ba5e726f6b2bccccd2c85d9481ad93e66b99b1e8e0f8c68ced3f06025",
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
                "q4_lumped_tetrahedral_capacitance_first_variation": (
                    "A1_ONLY__E_NULLITY_2_EXACT"
                ),
                "q4_pair_resolved_grounded_elastance_match": (
                    "EXACT_UD_D_MINUS_2_SQUARED_ON_ALL_16_Q4_STATES"
                ),
                "q4_central_kernel_DPAR_slope": "LAMBDA_EQUALS_R0_VPRIME_OVER_2V",
                "q4_ideal_fixed_coupling_coulomb_lambda": "MINUS_ONE_HALF",
                "q4_ideal_coulomb_Ud_relation": (
                    "2_KAPPA_SQUARED_ALPHA_MU0_HBAR_C_OVER_EPSILON_R_R0"
                ),
                "q4_pair_resolved_DPAR_realization_scope": (
                    "EXACT_CONDITIONAL_LOCAL_PHYSICAL_COMPLETION__NOT_CURRENT_F3_"
                    "GLOBAL_SHARED_LINK_VISIBLE_U1_OR_TRANSVERSE_MAXWELL_STIFFNESS"
                ),
                "q4_FV_PURE_status": (
                    "EXPLICIT_S10_COMPLETE_NONIDENTITY_SOURCE_PREMISE__STRONGER_THAN_FU_S1_TO_S9"
                ),
                "q4_FV_direct_projected_nonidentity_E_rank": 2,
                "q4_FV_H6_Hermitian_ring_nonidentity_rank": 4,
                "q4_FV_H6_Hermitian_ring_sectors": "A1_PLUS_T2",
                "q4_FV_projected_nonidentity_source_rank": 6,
                "q4_FV_operator_witness_determinant": "MINUS_4678629417_OVER_256",
                "q4_FV_formal_through_H8_rank": 6,
                "q4_FV_complete_CTP_Ward_tensor_response": "OPEN_NOT_COMPUTED",
                "q4_FW_source_scope": "FV_WITNESS_PAIR0_PLUS_RING6_IRREDUCIBLE_ONLY",
                "q4_FW_FV_family_offshell_rank": 6,
                "q4_FW_component_operator_rank_mod_identity": 5,
                "q4_FW_commutator_rank": 3,
                "q4_FW_ground_retarded_rank": 2,
                "q4_FW_first_nonzero_moment_rank": 2,
                "q4_FW_exact_gap_1_J6": "2_PLUS_2_SQRT2",
                "q4_FW_exact_gap_2_J6": "4_PLUS_2_SQRT2",
                "q4_FW_residue_ranks": (1, 1),
                "q4_FW_generated_diagonal_and_fold_scope": "OMITTED_NOT_UPPER_BOUND",
                "q4_FX_scope": (
                    "COMPLETE_HOMOGENEOUS_K0_SELECTED_FO180_FV_PURE_THROUGH_H6_MOD_IDENTITY"
                ),
                "q4_FX_source_off_coefficients": ("MINUS60", "MINUS35", "MINUS893_OVER9"),
                "q4_FX_fE": "1_MINUS_X2_MINUS_37_OVER12_X4_MINUS_16247_OVER900_X6",
                "q4_FX_component_operator_rank_mod_identity": 5,
                "q4_FX_commutator_rank": 3,
                "q4_FX_ground_retarded_rank": 2,
                "q4_FX_first_nonzero_moment_rank": 2,
                "q4_FX_finite_root_is_physical_threshold": False,
                "q4_FX_nonzero_momentum_source": "SUPERSEDED_BY_EXACT_FY_M1_REFINEMENT",
                "q4_FY_scope": (
                    "NATIVE_SUPPORT_M1_SELECTED_FO180_FV_PURE_THROUGH_H6__"
                    "SAMPLED_FINITE_RESPONSE"
                ),
                "q4_FY_support_species": ("A", "B", "E0", "E1", "E2", "E3"),
                "q4_FY_exact_m0_recovery": True,
                "q4_FY_exact_field": "Q_ZETA240_MOD_PHI240",
                "q4_FY_conjugate_mode": "Q29_EQUALS_Q1_DAGGER",
                "q4_FY_diagonal_lift_coefficients": (
                    "MINUS1",
                    "MINUS37_OVER12",
                    "MINUS16247_OVER900",
                ),
                "q4_FY_same_fE_as_homogeneous": True,
                "q4_FY_ring_offdiagonal_independent": True,
                "q4_FY_response_samples_x": ("2_OVER5", "1_OVER2"),
                "q4_FY_sampled_rank_hierarchy": "6_TO_6_TO_6_TO_6_TO_6",
                "q4_FY_TT_ground_image_rank": 2,
                "q4_FY_response_gap_count": 4,
                "q4_FY_residue_ranks": (1, 3, 1, 1),
                "q4_FY_naive_spatial_transversality": False,
                "q4_FY_temporal_density_current_contact": "OPEN_NOT_CONSTRUCTED",
                "q4_FZ_projected_incidence_charge": (
                    "ZERO_OPERATOR_ON_SELECTED_ICE_COMPONENT"
                ),
                "q4_FZ_H6_incidence_continuity": (
                    "EXACT_TRIVIAL_ZERO_CHARGE_CONTINUITY"
                ),
                "q4_FZ_supplied_embedding_source_contraction": "NONZERO_EXACT",
                "q4_FZ_TT_projector_rank": 2,
                "q4_FZ_full_spacetime_Ward": (
                    "OPEN_UNDECIDED_WITHOUT_PHYSICAL_DIVERGENCE_TEMPORAL_CURRENT_CONTACT"
                ),
                "q4_GA_fixed_total_charge_code_lift": "EXACT_UNITARY_EQUIVALENCE",
                "q4_GA_internal_charge_continuity": (
                    "QDOT_L_PLUS_I_LR_ZERO__QDOT_R_MINUS_I_LR_ZERO"
                ),
                "q4_GA_outer_port_continuity": (
                    "QDOT_R_MINUS_I_LR_PLUS_I_RBOUNDARY_ZERO"
                ),
                "q4_GA_FY_nonidentity_preservation": (
                    "EXACT_ONLY_ON_FULL_CODE_SOURCE_INDEPENDENT_SCALAR_HOLD"
                ),
                "q4_GA_spatial_bond_current_or_T0J": "OPEN_NOT_DERIVED",
                "q4_GB_local_ring_energy_continuity": "EXACT",
                "q4_GB_translation_kind": "GLOBAL_Z30_REPRESENTATION_LABEL_ONLY",
                "q4_GB_local_momentum_density_and_Ward": "OPEN_NOT_DERIVED",
                "q4_GC_family": "G_L_WITH_L_EQUALS_5_TIMES_2_TO_R",
                "q4_GC_vertices_links_native_supports_H6": (
                    "2L3__4L3__6L3__4L3"
                ),
                "q4_GC_kmin": "3PI_OVER_2_L_ASTAR",
                "q4_GC_exact_TT_projector_rank": 2,
                "q4_GC_count_normalization_is_physical_canonical": False,
                "q4_GC_massless_common_cone_pole": "OPEN_NOT_ESTABLISHED",
                "q4_GD_recoil_gate": (
                    "B1_ALGEBRAIC_EXISTENCE_ON_COMMON_AUXILIARY_MECHANICAL_TORUS"
                ),
                "q4_GD_total_momentum_ownership": (
                    "EXACT_EQUAL_AND_OPPOSITE_FACTOR_EDGE_IMPULSE_WITH_ADMITTED_BOUNDARY"
                ),
                "q4_GD_FY_source_preservation": (
                    "EXACT_ON_SOURCE_INDEPENDENT_FULL_CODE_SCALAR_HOLD"
                ),
                "q4_GD_physical_diamond_placement": "OPEN_NOT_DERIVED",
                "q4_public_GC_FY_spectral_payload_found": False,
                "q4_public_spectral_enabling_root": (
                    "ZHOU_2026_QFI_TABLES_AND_PYROCHLORE_QMC_ED_CODE_ONLY"
                ),
                "q4_GF_G2_observable_contract": (
                    "V005_NARROW_DESIGN_CONTRACT_INDEPENDENTLY_SEALED"
                ),
                "q4_GF_amplitude_classifier": (
                    "DISJOINT_TOTAL_LOWER_AND_UPPER_BOUNDED__ANCESTRY_SEPARATE"
                ),
                "q4_GF_helicity_gate": (
                    "POINCARE_LITTLE_GROUP_PLUS_MINUS2_REQUIRED__SCALAR_DOUBLET_REJECTED"
                ),
                "q4_GF_matched_spectral_physics": "OPEN_UNEXECUTED",
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
                "pair_resolved_central_kernel_DPAR": (
                    "EXACT_CONDITIONAL_LOCAL_PHYSICAL_REALIZATION"
                ),
                "ideal_coulomb_DPAR_lambda": "MINUS_ONE_HALF",
                "FV_PURE_complete_source_premise": "EXPLICIT_NOT_DERIVED_FROM_FU",
                "FV_projected_offshell_operator_rank": 6,
                "FV_projected_rank_sectors": "A1_1_PLUS_E_2_PLUS_T2_3",
                "FW_finite_homogeneous_response_hierarchy": "6_TO_5_TO_3_TO_2_TO_2",
                "FW_scope": "FV_WITNESS_SUBOPERATOR_NOT_COMPLETE_H6_SOURCE_OR_CTP",
                "FW_component_conservation_is_Ward_identity": False,
                "FX_complete_homogeneous_H6_response_hierarchy": "5_TO_3_TO_2_TO_2",
                "FX_generated_diagonal_source_direction": "PAIR_SOURCE_PLUS_HILBERT_IDENTITY_ONLY",
                "FX_finite_polynomial_zero_is_threshold": False,
                "next_collective_metric_calculation": (
                    "SAME_ANCESTRY_CANONICALLY_NORMALIZED_MASSLESS_HELICITY2_POLE_"
                    "WITH_COMMON_CONE_AND_FACTORIZATION"
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
                "DPAR_pair_resolved_coulomb_realization": (
                    "EXACT_CONDITIONAL_LOCAL_COMPLETION__GLOBAL_F3_GAUGE_SOLDER_OPEN"
                ),
                "FV_PURE_projected_offshell_rank6": "EXACT_CONDITIONAL_ON_EXPLICIT_S10",
                "FV_WITNESS_finite_homogeneous_response": (
                    "EXACT_6_TO_5_TO_3_TO_2_TO_2__NOT_COMPLETE_H6_SOURCE"
                ),
                "complete_generated_H6_diagonal_and_fold_response": (
                    "EXACT_PAIR_SOURCE_RENORMALIZATION_PLUS_IDENTITIES_ON_SELECTED_K0_COMPONENT"
                ),
                "native_support_nonzero_momentum_H6_response": (
                    "EXACT_M1_DIAGONAL_LIFT_PLUS_INDEPENDENT_RING__SAMPLED_FINITE_RESPONSE"
                ),
                "typed_m1_continuity_contact_Ward_boundary": (
                    "EXACT_ZERO_PROJECTED_CHARGE_AND_NONZERO_SUPPLIED_EMBEDDING_"
                    "CONTRACTION__FULL_WARD_UNDECIDED"
                ),
                "encoded_charge_current_lift": (
                    "EXACT_FIXED_CHARGE_CODE_AND_INTERNAL_OUTER_CONTINUITY__"
                    "SPATIAL_CURRENT_OPEN"
                ),
                "fixed_support_energy_momentum_boundary": (
                    "EXACT_RING_ENERGY_CONTINUITY__LOCAL_MOMENTUM_OPEN"
                ),
                "finite_family_common_cone_boundary": (
                    "EXACT_AFFINE_SOURCE_MOMENTUM_TT_KINEMATICS__"
                    "POLE_AND_CANONICAL_RESIDUE_OPEN"
                ),
                "translation_owning_recoil_parent": (
                    "EXACT_B1_AUXILIARY_TORUS_TOTAL_MOMENTUM_OWNERSHIP__"
                    "PHYSICAL_PLACEMENT_OPEN"
                ),
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
                "hust_calibrated_partial_AAF_kernel_kg_m_minus3": (
                    "6926P5053_TO_6926P5456_APPROXIMATE"
                ),
                "hust_calibrated_partial_ToS_kernel_kg_m_minus3": (
                    "24911P97_TO_25003P12_APPROXIMATE"
                ),
                "hust_public_signed_ToS_anelastic_corrections_ppm": (
                    "MINUS6P01",
                    "MINUS8P38",
                    "MINUS5P68",
                    "MINUS6P92",
                ),
                "hust_calibrated_partial_AAF_zero_remainder_G_SI": (
                    "6P6741755E_MINUS11_TO_6P6744092E_MINUS11"
                ),
                "hust_calibrated_partial_ToS_zero_remainder_G_SI": (
                    "6P6739974E_MINUS11_TO_6P6742407E_MINUS11"
                ),
                "hust_processed_kernel_comparator_gap_AAF_ppm": "18P86_TO_29P95_APPROXIMATE",
                "hust_processed_kernel_comparator_gap_ToS_ppm": "0P57_TO_12P11_APPROXIMATE",
                "hust_author_rule_covariance_ppm": ("11P616", "11P637"),
                "hust_independently_owned_physical_harmonic_remainder": "NOT_PUBLICLY_OWNED",
                "hust_raw_event_design_covariance_refit_ready": False,
                "hust_public_completion_search_status": (
                    "BOUNDED_NO_QUALIFYING_ROOT_ON_DECLARED_SURFACES"
                ),
                "hust_public_completion_search_world_exhaustive": False,
                "hust_public_completion_search_acquisition_leads": 2,
                "hust_public_completion_search_confirmed_dissertation_leads": 1,
                "hust_public_completion_search_numerical_G_advanced": False,
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
                "bense_public_history_samples": 101628,
                "bense_public_history_branches": 38,
                "bense_public_history_within_file_comparisons": 97,
                "bense_public_history_force_event_proxies": 34,
                "bense_matched_endpoint_lineage_estimand_identified": False,
                "bense_gravity_observable_present": False,
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
                    "PAIR_RESOLVED_COULOMB_DPAR_CONDITIONALLY_REALIZED__"
                    "FV_PURE_PROJECTED_OFFSHELL_RANK6__"
                    "COMPLETE_HOMOGENEOUS_H6_RESPONSE_5_TO_3_TO_2_TO_2__"
                    "NATIVE_M1_H6_SAMPLED_RESPONSE_6_TO_6_TO_6_TO_6_TO_6"
                ),
                "phase_A_candidate_status": (
                    "FULLY_SEALED_FY_NATIVE_RESPONSE_AND_GC_FINITE_FAMILY_KINEMATICS"
                ),
                "direct_shortcut_G1_status": (
                    "SEALED_B1_TOTAL_MOMENTUM_OWNERSHIP_ON_AUXILIARY_TORUS"
                ),
                "direct_shortcut_G2_status": (
                    "V005_OBSERVABLE_CONTRACT_SEALED__MATCHED_MASSLESS_HELICITY2_PHYSICS_OPEN"
                ),
                "direct_shortcut_G3_status": (
                    "SOFT_FACTOR_DEPENDENCY_THEOREM_IN_PROGRESS__"
                    "POSITIVE_CLOSURE_REQUIRES_G2"
                ),
                "next_no_lab_gravity_calculation": (
                    "ALL_GL_NATIVE_H6_LEDGER__THEN_MATCHED_MASSLESS_HELICITY2_POLE__"
                    "THEN_SOFT_UNIVERSALITY__THEN_EINSTEIN_SELF_COUPLING"
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
                "hust_public_calibrated_source": (
                    "CENTRAL_CORRECTIONS_AND_SIGNED_TOS_ANAELASTICITY_RECOVERED__"
                    "ONE_INDEPENDENT_ROW_HARMONIC_REMAINDER_STILL_REQUIRED"
                ),
                "hust_public_completion_search": (
                    "NO_QUALIFYING_ROOT_ON_DECLARED_BOUNDED_SURFACE__"
                    "TARGETED_ACQUISITION_OR_AUTHOR_REQUEST_NEXT"
                ),
                "public_history_and_clock_diagnostics": (
                    "HUST_AND_BENSE_HISTORY_CONFOUNDS_MAPPED_AND_K5_NODE_SCALAR_COMPATIBILITY_PASSED__"
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
                "q4_pair_resolved_coulomb_DPAR_realization_executed": True,
                "q4_FV_projected_offshell_rank6_executed": True,
                "q4_FW_witness_finite_response_executed": True,
                "q4_complete_generated_H6_source_response_executed": True,
                "q4_native_support_nonzero_momentum_response_executed": True,
                "q4_m1_continuity_contact_Ward_boundary_executed": True,
                "q4_encoded_charge_current_lift_executed": True,
                "q4_fixed_support_energy_continuity_executed": True,
                "q4_local_momentum_Ward_derived": False,
                "q4_finite_family_common_cone_boundary_executed": True,
                "q4_same_ancestry_massless_helicity2_pole_derived": False,
                "q4_translation_owning_recoil_B1_executed": True,
                "q4_physical_diamond_recoil_placement_derived": False,
                "q4_qualifying_public_spectral_payload_found": False,
                "q4_physical_metric_source_derived": False,
                "q4_full_CTP_Ward_response_executed": False,
                "hust_processed_dual_channel_forward_executed": True,
                "hust_roundtrip_history_diagnostic_executed": True,
                "hust_nominal_source_kernel_reconstruction_executed": True,
                "hust_conditional_homogeneous_G_crosscheck_executed": True,
                "hust_calibrated_partial_source_reconstruction_executed": True,
                "hust_bounded_public_completion_search_executed": True,
                "bense_public_history_development_analysis_executed": True,
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
                "pair_resolved_coulomb_DPAR_promotes_visible_EM_or_global_F3_solder": False,
                "FV_projected_rank6_promotes_complete_CTP_or_tensor_response": False,
                "FW_finite_witness_rank2_promotes_complete_response_or_Ward_tensor_pole": False,
                "FX_complete_homogeneous_rank2_promotes_local_Ward_tensor_pole_or_gravity": False,
                "FX_finite_polynomial_zero_promotes_physical_threshold": False,
                "FY_native_m1_response_promotes_continuum_locality_Ward_or_gravity": False,
                "FY_nonzero_spatial_contraction_promotes_complete_Ward_failure": False,
                "FY_sampled_rank6_promotes_generic_in_x_rank6": False,
                "FZ_zero_projected_charge_promotes_complete_Ward_packet": False,
                "GA_internal_charge_current_promotes_spatial_momentum_current": False,
                "GB_ring_energy_continuity_promotes_local_stress_Ward": False,
                "GC_affine_TT_kinematics_promotes_massless_common_cone_pole": False,
                "GD_auxiliary_recoil_promotes_physical_spacetime_placement": False,
                "public_QFI_QMC_root_promotes_GC_FY_spectral_evidence": False,
                "HUST_history_residual_promotes_lineage_gravity": False,
                "Bense_path_dependence_promotes_causal_lineage_or_gravity": False,
                "HUST_nominal_kernel_promotes_full_apparatus_G": False,
                "HUST_conditional_homogeneous_quotient_promotes_new_or_independent_G": False,
                "HUST_processed_comparator_remainder_promotes_independent_source_ownership": False,
                "HUST_bounded_no_root_search_promotes_global_data_nonexistence": False,
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
