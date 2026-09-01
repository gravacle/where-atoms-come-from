"""Zero-input URM checkpoint for sealed microscopic-gravity progress through GL6AW.

This additive surface verifies the frozen theorem/audit packets that establish
record-gated finite response, interaction-owned nonfactorization, an
authenticated finite relational atlas, quasi-local bulk dynamics, native
degree locking, the complete order-six locked Hamiltonian, thermodynamic
locked sectors, collective-response constraints, and the exact quantum-ice
comparison point; the first-character and anisotropic finite-size closure
theorems; and the exact record-conditioned pure-loop clock/typed-atlas bridge.
It deliberately does not alter the V014 working-theory certificate and does
not promote these results to an all-orders phase, selected GNS mode, physical
metric, gravity, Einstein dynamics, or a numerical value of G.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path, PurePosixPath
import posixpath
import re
from types import MappingProxyType
from typing import Any, Mapping, NoReturn


SCHEMA = "WAC_GRAVITY_MICROSCOPIC_PROGRESS_CERTIFICATE_V003"
CLAIM_CLASS = "SEALED_MICROSCOPIC_F3_PROGRESS_THROUGH_GL6AW_WITH_STRICT_GNS_METRIC_AND_GRAVITY_CEILINGS"

_REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
_HASH_ROW = re.compile(r"^([0-9a-f]{64})  ([^\r\n]+)$")


@dataclass(frozen=True)
class _PacketPin:
    label: str
    author_dir: str
    author_claim_sha256: str
    author_manifest_sha256: str
    author_seal_sha256: str | None
    audit_dir: str
    audit_sha256: str
    audit_manifest_sha256: str
    audit_seal_sha256: str
    disposition: str
    author_claim_file: str = "THEOREM.md"
    audit_claim_file: str = "AUDIT.md"


_PACKETS = (
    _PacketPin(
        "GL6T",
        "LANE_CROSS_RFT_GRA_GL6T_F3_LINEAGE_GATED_Q4_E2_RESPONSE_V001",
        "c610f0dc4092e8ca738c1b1a32a5e034f1ccbaafb032d3c1b5afdf3efaca539b",
        "c64f60789f51247a70cace0ee3626bfb7d3c5eaad36c318f7d674b5cc4a797d9",
        None,
        "AUDIT_G_GL6T_F3_LINEAGE_GATED_Q4_E2_RESPONSE_V001",
        "cefa812da0cc0b75fd689097326f74f30d679571ff5001d0fa00d09632f944b1",
        "eef159a264065610ad7ded02f1eb828e0334c9e135b13e66b59609466be72718",
        "bd14425c1016e98fe0c96427b97a3898537465e490c695fb999bf5b4f3f03a6e",
        "PASS_AT_EXACT_FINITE_BRANCHWISE_PAIR_RESPONSE_CAPACITY_SCOPE__NO_COLLECTIVE_STIFFNESS_RICCI_GRAVITY_OR_G_PROMOTION",
    ),
    _PacketPin(
        "GL6U",
        "LANE_CROSS_RFT_GRA_GL6U_F3_DEGREE_INTERACTION_Q4_RESPONSE_V001",
        "3f23084cd7e5fa500d331c6e9739b6f0b3e7875bdce3b33f5369eafd54bc7965",
        "7297d5b6456c106a7267fcd8fa4d2aa09d72a19f23605e0becde1ad14b2dba04",
        None,
        "AUDIT_G_GL6U_F3_DEGREE_INTERACTION_Q4_RESPONSE_V001",
        "3d0ccff68227c0bbc5db4ae32b00b714e8b2cd7892a5f4dcd5284671ad8cd052",
        "6b74caff9f35958c27204e78ab464a38fb8483963a26818891047d19a6d1ebba",
        "b5d47d99414d851f61963254b3146289f64c7c6bd3fad416662eba63cf36767b",
        "PASS_AT_EXACT_FINITE_INTERACTION_OWNED_FACTORIZATION_DEFECT_AND_FULL_PAIR_RESPONSE_SCOPE__NO_COLLECTIVE_STIFFNESS_RICCI_GRAVITY_OR_G_PROMOTION",
    ),
    _PacketPin(
        "GL6AA",
        "LANE_CROSS_RFT_GRA_GL6AA_RECORD_AUTHENTICATED_SHARED_CHILD_ATLAS_V001",
        "faea49e3dcd5f2b4d5b3bab9026432d741192339dd789a939ef2318236848c0e",
        "9bebed96f2b864738639571f870ae7a34dc440dea5ab4418bc4e7cc9c2eb2a63",
        "e7d6b64a713fa37170fa70ecd7978a93028fc1594b4e2bd6e8befffad61d3004",
        "AUDIT_G_GL6AA_RECORD_AUTHENTICATED_SHARED_CHILD_ATLAS_V001",
        "cc34f4f8ea3824f6788209db4e3b9a1e03f034d38048b1098b221f5d741b0e0a",
        "1287cfd1f3755b3995bc1aa023b093f9532b15dfeb727e55737718ed852160d2",
        "4204c93bc373c4551660eb391c4a139979afe90f20f261c6abf3f234e85f4499",
        "PASS__OLD_PAIRWISE_QUERY_NONIDENTIFYING__SELECTED_INDEPENDENT_NODE_ID_AND_PORT_FRAME_COMPLETION_EXACT__LITERAL_SHARED_CHILD_AND_COCYCLES_AUTHENTICATED__AUTONOMOUS_SUPPORT_LENGTH_CONE_RICCI_GRAVITY_AND_G_OPEN",
    ),
    _PacketPin(
        "GL6AF",
        "LANE_CROSS_RFT_GRA_GL6AF_FORMATION_PATTERN_E2_SOURCE_V001",
        "6e20f188aa84a369d90eb67d741eb0c8aa0eb89013598b57021fa7464249fc71",
        "f2f4b416db1caa3a935e641cf145a488728b7281325a0e2a4861d6181e48677f",
        "7d440fb501063acd960658bb5318f46483568599e5a272669e48b1992c8a9f86",
        "AUDIT_G_GL6AF_FORMATION_PATTERN_E2_SOURCE_V001",
        "5566872b05d8ee04c67901f2f92199d129618bfd884df57f3156bb8355da6d81",
        "3a0dae720356080d151328259ffabc170414e30b8e8b46770feff76d5a1d2308",
        "a3c0e0803651bff50a764aa291b8be24034baa0448f36cdb51d6fb9cf7a8afd8",
        "PASS__ALL_16_FORMATION_PATTERN_RESPONSE_EXACT__FIXED_E_RESTRICTION_RANK_THRESHOLD_EXACT__BROKEN_S4_MIXING_CEILING_EXACT__K_BRANCH_PROJECTOR_AND_CTP_NORMALIZATION_LAWFUL__PHYSICAL_K_NOT_SEMANTIC_REC__NO_COLLECTIVE_STIFFNESS_OR_GRAVITY_PROMOTION",
    ),
    _PacketPin(
        "GL6AG",
        "LANE_CROSS_RFT_GRA_GL6AG_N1_FORMATION_NEIGHBOR_PROPAGATION_V001",
        "8551a4dc37183b8ab83ac48a0774dc0b82bd280faa5f9e469272f8c7ef898dea",
        "f69fd10402b1ec428b95f826ba19ec5fb9635a2079b05c0b49514b3084979b6e",
        "12b3e091624a8be4233b342a5303ec09439e140003ea5b770a7c7323e91d55c7",
        "AUDIT_G_GL6AG_N1_FORMATION_NEIGHBOR_PROPAGATION_V001",
        "9b313ccc6b15161e088531d51585479831ad4d0995e537f6353a04581d9a8112",
        "b2e1e59425fb050cfb3fb21ff0a4d25057009fded74181d81b080b01aeac19d1",
        "210b84d267511698df94c2b7e3980d7057eb356dbbff252dd3322924b3c91d63",
        "PASS__FULL16_COMMON_PARENT_AND_K_BRANCH_RESTRICTION_LAWFUL__MATCHED_NOT_ABSOLUTE_RECEIVER_CONTRAST__FIXED_E_COLUMNS_TYPED_AND_BROKEN_S4_CEILING_KEPT__ALL16_Q12_AND_PAIR_MOBIUS_Q16_EXACT__BRIDGE_ABLATION_FACTORIZATION_DIAGNOSTIC_ONLY__K_NOT_REC__NO_BULK_RICCI_GRAVITY_OR_G_PROMOTION",
    ),
    _PacketPin(
        "GL6AH",
        "LANE_CROSS_RFT_GRA_GL6AH_BULK_BOUNDARY_IRREP_SEPARATION_V001",
        "79b04596c6df950a86bfe25fb02f8cf7822d65f6eb6101615deba7d7f88b58eb",
        "f83a7524592be535a8e21e8fc94fdbcc0ba0f6c6d874779bf904769f3c115367",
        "5a94269aca3092c4be83c96738cf977e37369150ea7097a8aed6f000a0c61ef8",
        "AUDIT_G_GL6AH_BULK_BOUNDARY_IRREP_SEPARATION_V001",
        "c393b2c98304d3b83dd4b02d5de0cff8f4879001644590365bad44c17a3122f0",
        "a4a1c5aa81adce48f775ba192ea6207f498c7d6e1bd4356c7e3e9bebe57a0433",
        "0915014ca7cf8468225c42d65e467ec3ad8eb5950e229e6c2dcf95dc9c897a26",
        "PASS__DIRECT_Q6_FULL_SIX_A1_PLUS_T2_AND_E_NULL_EXACT__GENERAL_EDGE_WORD_BOUND_LAWFUL__Q12_Q16_RECEIVER_HELPER_SIGNS_AND_SUPPORT_OWNERSHIP_EXACT__ETA_HOMOGENEOUS_CANCELLATION_EXACT__DELTA0_IDENTITIES_AND_CONTINUITY_SCOPE_SEPARATED__NO_BULK_SHEAR_CONE_RICCI_GRAVITY_OR_G_PROMOTION",
    ),
    _PacketPin(
        "GL6AI",
        "LANE_CROSS_RFT_GRA_GL6AI_F3_RELATIONAL_INFLUENCE_ENVELOPE_V001",
        "a51e802f6ba148e5f9848e95f41a80073795b24b7eaf87e36c0766b0856aa494",
        "fc50cad54dca00aab1c30d7c12ef07147df1242f94483f63955185695073f706",
        "12daa03d45cd653db24622ae8b3d8166015291534b3295e3c426eb37180fc918",
        "AUDIT_G_GL6AI_F3_RELATIONAL_INFLUENCE_ENVELOPE_V001",
        "5734ad57122c64e3174aa7706b0e7aa86102b3a18a3b868aca20af0997ab462a",
        "b6f2b4c98f7747f074f1793efc3ccaee15593667d500d892fa384ec5678e8705",
        "b1f80a18353a647c31fd5d58bbbacb72c1977ed2c3d9f43f81b33020dbe28c5e",
        "PASS__EXACT_DEGREE_SQUARE_SPLIT__LINK_DEGREE_Q_PLUS2_LE6__DRESSING_COMPLETE_LAMBDA48_EXACT__RETAINED_SOURCE_K_DUHAMEL_TYPED__LINK_TO_CELL_DISTANCE_DESCENT__UNIFORM_ANALYTIC_EXPONENTIAL_TAIL__NO_EXACT_SPEED_LORENTZ_RICCI_GRAVITY_OR_G_PROMOTION",
    ),
    _PacketPin(
        "GL6AK",
        "LANE_CROSS_RFT_GRA_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001",
        "083d5fbb8a48e27e365167075da132ffa23e395587a4c0e40cc572d8b761ad30",
        "d38f89c618ea6f77c7b399b005ad0f0abe04d3865e06921f8c765feb44f40620",
        "322bf51a00f8fea3f36a09656dda4ebf89ba56b9a88d60b50e9cc7ab33223987",
        "AUDIT_G_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001",
        "ec1e452e51fd381ad23bc24f49d4613132bcc570b778697e72ac63fd006cd22e",
        "01d9dd3e2a7cc3247ea11719b5d69c1e0c82f3370770fef1b91fe98452d12839",
        "3137e843b9abc98ecdececeb67a204c56993b673e79562dde87f9c2588f1fe7f",
        "PASS__FINITE_AUTHENTICATED_ANCESTRY_AND_A3_DEGREE6_EXACT__BOUNDARY_FACTOR3_AND_CAUCHY_TAIL_EXACT__JOINT_INVARIANT_STATE_EXISTENCE_SOUND__POSITIVE_A1_E_T2_LIOUVILLIAN_MEASURES_SOUND__NO_GLOBAL_RECORD_STATE_SELECTION_POLE_PHYSICAL_MOMENTUM_RICCI_GRAVITY_OR_G",
        audit_claim_file="POSTFREEZE_AUDIT.md",
    ),
    _PacketPin(
        "GL6AM",
        "LANE_CROSS_RFT_GRA_GL6AM_AUTHENTICATED_BULK_RESPONSE_FUNCTIONAL_V001",
        "8407cee5196bfa4240f02159a5f59f941903dcf7a10e2baa18cf52a01ac8f743",
        "8b4ac6f6ceda2acb117480201ee96ce22be97fe0a99d8c097d8267100efa8c44",
        None,
        "AUDIT_G_GL6AM_AUTHENTICATED_BULK_RESPONSE_FUNCTIONAL_V001",
        "2659392e6f2fe3c0062068426faf5d516cceaf3a106017742b8ea88c21517b00",
        "d35a2a9b2e581db963ca0513a26ffbcbbbbae28efd7c3a019dfa4a9f0db50301",
        "d1d70e294e4be73f8efdc26301a37bf11808f9f248b433b3ec5aaf7f87e324c7",
        "PASS__FINITE_AUTHENTICATED_PULSE_AND_DEFECT_WINDOWS_HAVE_BOUNDARY_INDEPENDENT_BULK_FUNCTIONALS__RETARDED_THETA_AND_FACTORIAL_TAIL_EXACT_NOT_STRICT_CONE__FINITE_WINDOW_CORRELATION_MEASURE_POSITIVE__A1_E_T2_ONLY_S4_CLOSED__K_DEFECT_COCYCLE_AND_STATE_INDEPENDENT_TAIL_SOUND__NO_STATE_SELECTION_BULK_COEFFICIENT_MOMENTUM_CONE_GRAVITY_OR_G",
    ),
    _PacketPin(
        "GL6AN",
        "LANE_CROSS_RFT_GRA_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001",
        "32f597edc51a609a37b86144487cd7db3bd2f14a65adb754a893d47ef6807e81",
        "24a71c01ed1b7a92830e92ec7682882c892667289e2794dafb4af5905ad71b2e",
        "a946902f027c555f91cd1f2e9ce93e3182f8edeca319955cd691a0bc929fba51",
        "AUDIT_G_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001",
        "8754210b2ff0077e8cec4c5ce0f771ba40cfee2b8957d37f5be64a50ba49d0b4",
        "73b618e88b96d40ca40e32bd33ab578df98ed20af9fc359d90bfe4ae75f5c91b",
        "d0a50cda599842d8854db0bc2ab9e665f823e2ab7e2048a0eafa672ae3ad7b7e",
        "PASS__NATIVE_DEGREE_LOCK_AND_LINEAR_WARD_NO_GO_EXACT__LOCAL_LOCKED_PAIR_VARIATIONS_EXACTLY_E__GRAM_EIGENVALUE_QUADRATIC_SINGULAR_VALUE_LINEAR__DECLARED_Q4_GIRTH6__CANONICAL_H2_H4_SCALAR__HEXAGON_MINUS63_OVER8__NO_RECORD_POLE_CONE_GRAVITY_OR_G_PROMOTION",
    ),
    _PacketPin(
        "GL6AO",
        "LANE_CROSS_RFT_GRA_GL6AO_COMPLETE_SIXTH_ORDER_LOCKED_HAMILTONIAN_V001",
        "f75edcb115c3f7c86c6598f4597366b36e363df2d03ad919cc607b57dfb6b20c",
        "c690665043fbbb277aae307a4308e8d30a41f0fbf87be8b1501d0ba86874a494",
        "9df82d1cdc53822bb88b1d419f67db367caf872e2168a9bbeb69d1a6acc9f0ae",
        "AUDIT_G_GL6AO_COMPLETE_SIXTH_ORDER_LOCKED_HAMILTONIAN_V001",
        "84dc32f946e6f60fc4822f1b5cf39b0f1a73c857911ab87301670f8adf4d174c",
        "8a3f438f8f894e6996cf54b2dcd994838bfb325e99f09fac99da0136af45ade5",
        "1870aefe176aeacaa127d4afd10458c4565caa25abc49c91096da5129f01040b",
        "PASS__CANONICAL_KATO_ORDER6_FORMULA_EXACT__DIRECT_AND_FOLDED_WORD_CENSUS_COMPLETE__M3_M2_CANCEL__COMMON_DIAGONAL_MINUS893_OVER1080_M__ONLY_ALTERNATING_HEXAGON_OFFDIAGONAL_MINUS63_OVER8__Q4_AND_FORMAL_LINKED_SCOPE_SOUND__NO_PHASE_POLE_MOMENTUM_CONE_GRAVITY_OR_G",
    ),
    _PacketPin(
        "GL6AP",
        "LANE_CROSS_RFT_GRA_GL6AP_LOCKED_IR_CONDITIONAL_RESPONSE_V001",
        "aab21a99ecd0c6696084dd0a22ed2489533588f2171bf0c319ffc00a4922033e",
        "1c6a3db34bda121f7b7fdc64a85aaed0d146dcb12d8f960449d10d873ff94e1f",
        "f56c24594e041485edc6d14cff91a3317bc400fc92ece595ef590291331a8288",
        "AUDIT_G_GL6AP_LOCKED_IR_CONDITIONAL_RESPONSE_V001",
        "d89f7a8ebfbf561254ccb924134328f3ff7f536f4f4bb4256a9721c2812fc3e3",
        "45e5a88bd01fe5d03ea9197f02c35437776f6684682d91d1e9711d6f28e1e814",
        "372223ff43c4e4a610a66a5d93ca56955cb18678b7602b172462e6b2026a602a",
        "PASS__LINK_KERNEL_T2_NOT_PAIR_E__HOM_S4_ZERO__ALL24_Q4_AUTOMORPHISMS__NATIVE_LOOP_BREAKS_E_COUNT__MASS_AND_CUBIC_ALLOWED__RECIPROCAL_QUADRATIC_CLASSIFICATION_EXACT__SPECTRAL_AND_PHYSICAL_MOMENTUM_CEILINGS_SOUND",
    ),
    _PacketPin(
        "GL6AQ",
        "LANE_CROSS_RFT_GRA_GL6AQ_AUTHENTICATED_E_LOOP_SELECTION_OBSTRUCTION_V001",
        "1d1b01380ec8fd7ce83c69d45b68d9bde36bbe1dacdd32e3a5909ee6723a5ace",
        "1adb601b5c0a957e9182ecfb89fd24fc2e63a318a7006759826911592d979afd",
        None,
        "AUDIT_G_GL6AQ_AUTHENTICATED_E_LOOP_SELECTION_OBSTRUCTION_V001",
        "8d0a451c1400244337f2c2fa08395fe3879f2e811a09249d8d1a194c31b5ea94",
        "f740eb62935dfd7b981b0252a4802d9ae2e5e5122c11563ebfa89dcba8bb4af9",
        "e795534605156e510daaee1cfa502d298186d2e57d09efcd60531e99167bb91f",
        "PASS__AUTHENTICATED_PAIR_QUERY_HAS_EXACT_NONZERO_LOCAL_LOCKED_E_OVERLAP__TRANSVERSE_K_HAS_ZERO_DIRECT_LOCKED_AND_ONE_CELL_S4_LINEAR_E_PROJECTION__SIX_RETAINED_SUPPORTS_GATE_MINUS63_OVER8_E_CHANGING_LOOP__TRACE_REFUTES_ONLY_UNIVERSAL_NONZERO_STATIONARY_CONTRAST__EXISTENTIAL_SELECTED_LOCKED_BULK_WITNESS_OPEN__NO_STATE_POLE_CONE_GRAVITY_OR_G",
    ),
    _PacketPin(
        "GL6AR",
        "LANE_CROSS_RFT_GRA_GL6AR_LOCKED_HEXAGON_THERMODYNAMIC_SECTOR_V001",
        "960af10c683496ac921c3371a8182db7f018cf788f38dc486ecb0af95c089555",
        "cdb1acc430b9c928807d8daabbf85848c40aeb1469da9a8fa4ec052ce34042cc",
        "42860871df4d86af4ce554b35460f928aade1f470de27fd617ebe9fc3dd15466",
        "AUDIT_G_GL6AR_LOCKED_HEXAGON_THERMODYNAMIC_SECTOR_V001",
        "f0fd37305d9278b29814816d338f574c9a51c2001011715e92ea085253cd9a04",
        "67cbd61e89ea1ec10d8cf3c9be8c87f06204b00f61e0e9517e6d2b4ff0cdfbf3",
        "9a84a348f4034f69762d61590c2630272c1787aff9647b541240b6326b77c935",
        "PASS__FINITE_AND_INFINITE_LOCKED_DYNAMICS_WELL_DEFINED__FLUX_AND_COMPONENT_SCOPE_EXACT__PF_AND_DIRICHLET_IDENTITIES_EXACT__FINITE_PERIODIC_ACTIVE_ENERGY_BOUND_ONLY__DELTA_LE_18_T_L_OVER_VARIANCE_EXACT_FOR_POSITIVE_VARIANCE__FINITE_SIZE_CLOSURE_CONDITIONAL__FULL_GNS_BRIDGE_STILL_MISSING__NO_MOMENTUM_PHOTON_CONE_GRAVITY_OR_G",
    ),
    _PacketPin(
        "GL6AS",
        "LANE_CROSS_RFT_GRA_GL6AS_NATIVE_HEXAGON_COLLECTIVE_RESPONSE_V001",
        "bfe36071a24ccc7d6d7a16afeeea1b5554a95562ae91ac59c709db478000db9f",
        "0b7c12e51ff2892cc44e0e4e39d68b6939c5fdb47fce8c018db1ad09191e0f0f",
        "9e7f058f68e989a00e16a17611246392fb4cdc9957a24a0f770ac530d029a957",
        "AUDIT_G_GL6AS_NATIVE_HEXAGON_COLLECTIVE_RESPONSE_V001",
        "eb4b17b3bf0429f41f5d09e5af9d395cb427223298be541b60551815b8c00fe1",
        "15c0498389fbd0424ead7ee17d0cfe4b95f119bdc17d36816c811121c844084d",
        "fda5c90d62a3b6bec0fcd75914f0cc5d9c7d066e970da601d26dac3cbf622c91",
        "PASS__PURE_HEXAGON_PORT_T2_CONSERVATION_AND_CYCLE_CONTINUITY_EXACT__SMA_REQUIRES_STRUCTURE_FACTOR__HARMONIC_TWO_CHARACTER_MODES_CONDITIONAL__PAIR_E_IS_COMPLEMENT_EVEN_QUADRATIC_WITH_CONDITIONAL_TWO_T2_CHANNEL__RETAINED_SOURCE_A1_PLUS_T2_ZERO_CHARACTER_E_CROSS__SYM2_TRACELESS_T2_EQUALS_E_PLUS_T2_ALGEBRA_ONLY__NO_POLE_PHYSICAL_CONE_STRESS_GRAVITY_OR_G",
    ),
    _PacketPin(
        "GL6AT",
        "LANE_CROSS_RFT_GRA_GL6AT_PRIMARY_QUANTUM_ICE_CROSSWALK_V001",
        "1c2b48c40a88000a88f8446fa2aea5116acb0896d94848274c87a22003bbcad9",
        "9acb465079e8dbf5fa15ffc6b08d8c4d34c65ee086f1f822e8276156f484ce61",
        "2a668784df7a4f3702e3e27c076ca8c9b2d15394e3f0e46135d6e440649ee575",
        "AUDIT_G_GL6AT_PRIMARY_QUANTUM_ICE_CROSSWALK_V001",
        "ebb24afb71cff956a70baa5ad787f983ce2a863a4e5742890a3255aa3cc648fe",
        "fa474b179408ec7852d613aea07d9c44e38133b79035d428166ef0f19134d1fa",
        "2ac91c1b48541e6f525fab905853d4e919ff293a00967341987d1819d64c0dd6",
        "PASS__EXACT_ORDER6_FPL_POINT_V_OVER_G_ZERO__PRIMARY_SOURCE_LABELS_SOUND__NO_SPECTRAL_OR_GRAVITY_PROMOTION",
        author_claim_file="RESULT.md",
    ),
    _PacketPin(
        "GL6AU",
        "LANE_CROSS_RFT_GRA_GL6AU_VG0_STATIC_STRUCTURE_PHASE_CLOSURE_V001",
        "ce2cdf1394053778300c71f3b2d25f79914efe82390b5b66f84fde807e0f2612",
        "347bfade725e1cbefbe6d205f596ea059e550d98f08aa452591bb47a45ae9705",
        "3a55ca000cc044305d5daae932a0c935cc73131f4c3f55db716d82753d1f7a4f",
        "AUDIT_G_GL6AU_VG0_STATIC_STRUCTURE_PHASE_CLOSURE_V001",
        "40fc6dabad36cbba2dbf5d1eb6a09b06d1d04577edf3e355e64177f094804a13",
        "fe726ece60a1dc386607221736f9d2092d681cb0673396e5eaeb305325ad8ef3",
        "c0ff5575e5927c7e74496b1515c487b7c5b88348bd782413928db29cb8f3200a",
        "PASS__FIRST_CHARACTER_STATIC_EXPONENT_CLOSURE_EXACT__ICE0_STATIC_AND_GNS_BRIDGE_EXPLICITLY_UNPROVED__NO_PHYSICAL_CONE_STRESS_GRAVITY_OR_G",
    ),
    _PacketPin(
        "GL6AV",
        "LANE_CROSS_RFT_GRA_GL6AV_RECORD_CONDITIONED_COLLECTIVE_METRIC_BRIDGE_V001",
        "d448cd4667f39405b0173185646e3db1f65435fd6eff60d2f35485e4d7958102",
        "283e1010be399b5e31cf93da34dee4c075dc29f9a98d9b4804ba6ea2e411073b",
        "c2ea7f69d86542583f58dffa5b8f6207c1b28b77c941ee016e1ac5f69c76b191",
        "AUDIT_G_GL6AV_RECORD_CONDITIONED_COLLECTIVE_METRIC_BRIDGE_V001",
        "e88fd0eab85b8f1a1a59d40d2f8f61b34bbcb03e087f57126c7f8921f12ba86a",
        "a3af036aa3de2af9fa9f3c0a3f0042e509978d18c3a026126c0cfcb438acf49f",
        "2a8276b6b2cb1f450072ed9a6d8d70c58bbe75b12e79418e5fbdc80fb6fcb6cc",
        "PASS__EXACT_PURE_LOOP_CLOCK_AND_TYPED_ATLAS__OPERATIONAL_METRIC_STRICTLY_CONDITIONAL__NO_GRAVITY_PROMOTION",
    ),
    _PacketPin(
        "GL6AW",
        "LANE_CROSS_RFT_GRA_GL6AW_ANISOTROPIC_FOLNER_TWIST_CLOSURE_V001",
        "ad315f2434b1042f183eeeb0244f9d323214717bebf01643c46eccd110373f9d",
        "73b0055d4fca77c044282fe80c9de78e3976c9bc77b37c31037d75f8e1ecb234",
        "382b10d28301c536f007fe79241167a1a830b68ce47135dbd5dfca37850b93e0",
        "AUDIT_G_GL6AW_ANISOTROPIC_FOLNER_TWIST_CLOSURE_V001",
        "e3a43679e686ddde45b40bdd737a610b2fad608a0c76495c7016df11e6224c89",
        "d26e96d417795664c9bbea492e216d9c3e4ebbeecaaad5a0478a280f6965ac6d",
        "42eab5756a691b89e5d46cb88d4669f95d0b35f2c21dd59982c804506593684a",
        "PASS__EXACT_FINITE_CENTERED_SECTOR_DICHOTOMY__ANISOTROPIC_FOLNER_CLOSURE__GNS_AND_PHYSICAL_PROMOTION_EXPLICITLY_OPEN",
    ),
)


class GravityMicroscopicProgressRefusal(RuntimeError):
    """A pinned microscopic-progress artifact failed closed custody checks."""


def _refuse(message: str) -> NoReturn:
    raise GravityMicroscopicProgressRefusal(
        "GRAVITY MICROSCOPIC PROGRESS REFUSES: " + message
    )


def _sha256(path: Path) -> str:
    try:
        if not path.is_file() or path.is_symlink():
            _refuse(f"custody object is absent, non-file, or symlinked: {path}")
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as exc:
        _refuse(f"custody object is unreadable: {path}: {exc}")


def _strict_text(path: Path) -> str:
    try:
        text = path.read_bytes().decode("utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        _refuse(f"text custody object is unreadable or not UTF-8: {path}: {exc}")
    if "\r" in text or "\x00" in text or not text.endswith("\n"):
        _refuse(f"text custody object has forbidden bytes or lacks final newline: {path}")
    return text


def _root_path(relative: str) -> Path:
    posix = PurePosixPath(relative)
    if posix.is_absolute() or ".." in posix.parts or not posix.parts:
        _refuse(f"unsafe custody path: {relative}")
    path = _REPOSITORY_ROOT.joinpath(*posix.parts)
    try:
        path.relative_to(_REPOSITORY_ROOT)
    except ValueError:
        _refuse(f"custody path escapes repository: {relative}")
    return path


def _declared_target(declaring_relative: str, target_relative: str) -> tuple[str, Path]:
    """Resolve both repository-rooted and packet-local historical hash rows."""
    target_posix = PurePosixPath(target_relative)
    if target_posix.is_absolute():
        _refuse(f"absolute hash target in {declaring_relative}: {target_relative}")
    root_candidate = None
    if ".." not in target_posix.parts:
        root_candidate = _root_path(target_relative)
    parent = PurePosixPath(declaring_relative).parent
    local_relative = posixpath.normpath(str(parent / target_posix))
    if local_relative == "." or local_relative.startswith("../"):
        _refuse(f"hash target escapes repository in {declaring_relative}: {target_relative}")
    local_candidate = _root_path(local_relative)
    root_exists = bool(
        root_candidate is not None
        and (root_candidate.exists() or root_candidate.is_symlink())
    )
    local_exists = local_candidate.exists() or local_candidate.is_symlink()
    if root_candidate is not None and root_candidate != local_candidate and root_exists and local_exists:
        _refuse(
            f"ambiguous root-versus-packet-local hash target in {declaring_relative}: "
            f"{target_relative}"
        )
    if local_exists and not root_exists:
        return local_relative, local_candidate
    if root_candidate is None:
        return local_relative, local_candidate
    return target_relative, root_candidate


def _verify_hash_list(relative: str, visited: set[str]) -> int:
    """Verify one declared hash list and recursively follow declared hash lists."""
    if relative in visited:
        return 0
    visited.add(relative)
    path = _root_path(relative)
    rows: set[str] = set()
    count = 0
    for line in _strict_text(path).splitlines():
        match = _HASH_ROW.fullmatch(line)
        if match is None:
            _refuse(f"malformed hash row in {relative}")
        expected, target_relative = match.groups()
        if target_relative in rows:
            _refuse(f"duplicate hash target in {relative}: {target_relative}")
        rows.add(target_relative)
        canonical_relative, target = _declared_target(relative, target_relative)
        if _sha256(target) != expected:
            _refuse(f"hash mismatch for {target_relative} declared by {relative}")
        count += 1
        same_packet = (
            PurePosixPath(canonical_relative).parent
            == PurePosixPath(relative).parent
        )
        if same_packet and target.name in {
            "MANIFEST.sha256",
            "DEPENDENCIES.sha256",
            "AUDITED_TARGETS.sha256",
        }:
            count += _verify_hash_list(canonical_relative, visited)
    if not rows:
        _refuse(f"empty hash list: {relative}")
    return count


def _verify_seal(directory: str, expected_seal: str, expected_manifest: str) -> None:
    relative = f"{directory}/SEAL.sha256"
    path = _root_path(relative)
    if _sha256(path) != expected_seal:
        _refuse(f"seal-file hash mismatch: {relative}")
    expected_text = f"{expected_manifest}  {directory}/MANIFEST.sha256\n"
    if _strict_text(path) != expected_text:
        _refuse(f"seal does not exactly bind the pinned manifest: {relative}")


@dataclass(frozen=True)
class _Custody:
    packet_rows: tuple[Mapping[str, Any], ...]
    declared_hash_rows_checked: int


def _verify_custody() -> _Custody:
    visited: set[str] = set()
    rows = []
    checked = 0
    for pin in _PACKETS:
        author_claim_relative = f"{pin.author_dir}/{pin.author_claim_file}"
        author_manifest_relative = f"{pin.author_dir}/MANIFEST.sha256"
        audit_claim_relative = f"{pin.audit_dir}/{pin.audit_claim_file}"
        audit_manifest_relative = f"{pin.audit_dir}/MANIFEST.sha256"

        if _sha256(_root_path(author_claim_relative)) != pin.author_claim_sha256:
            _refuse(f"{pin.label} author-claim hash mismatch")
        if _sha256(_root_path(author_manifest_relative)) != pin.author_manifest_sha256:
            _refuse(f"{pin.label} author-manifest hash mismatch")
        if _sha256(_root_path(audit_claim_relative)) != pin.audit_sha256:
            _refuse(f"{pin.label} audit-claim hash mismatch")
        if _sha256(_root_path(audit_manifest_relative)) != pin.audit_manifest_sha256:
            _refuse(f"{pin.label} audit-manifest hash mismatch")

        checked += _verify_hash_list(author_manifest_relative, visited)
        checked += _verify_hash_list(audit_manifest_relative, visited)
        if pin.author_seal_sha256 is not None:
            _verify_seal(pin.author_dir, pin.author_seal_sha256, pin.author_manifest_sha256)
        _verify_seal(pin.audit_dir, pin.audit_seal_sha256, pin.audit_manifest_sha256)

        audit_text = _strict_text(_root_path(audit_claim_relative))
        marker = f"**Disposition:** `{pin.disposition}`"
        if marker not in audit_text:
            _refuse(f"{pin.label} accepted audit disposition is absent")
        rows.append(
            MappingProxyType(
                {
                    "gate": pin.label,
                    "author_directory": pin.author_dir,
                    "author_claim_file": pin.author_claim_file,
                    "author_claim_sha256": pin.author_claim_sha256,
                    "theorem_sha256": (
                        pin.author_claim_sha256
                        if pin.author_claim_file == "THEOREM.md"
                        else None
                    ),
                    "author_manifest_sha256": pin.author_manifest_sha256,
                    "author_seal_sha256": pin.author_seal_sha256,
                    "audit_directory": pin.audit_dir,
                    "audit_claim_file": pin.audit_claim_file,
                    "audit_sha256": pin.audit_sha256,
                    "audit_manifest_sha256": pin.audit_manifest_sha256,
                    "audit_seal_sha256": pin.audit_seal_sha256,
                    "audit_disposition": pin.disposition,
                }
            )
        )

    # Close the most exposed read-time race over the primary claim bytes.
    for pin in _PACKETS:
        if _sha256(
            _root_path(f"{pin.author_dir}/{pin.author_claim_file}")
        ) != pin.author_claim_sha256:
            _refuse(f"{pin.label} author claim changed during verification")
        if _sha256(
            _root_path(f"{pin.audit_dir}/{pin.audit_claim_file}")
        ) != pin.audit_sha256:
            _refuse(f"{pin.label} audit claim changed during verification")
    return _Custody(tuple(rows), checked)


def _freeze(value: Any) -> Any:
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    return value


def _certificate(custody: _Custody) -> Mapping[str, Any]:
    return _freeze(
        {
            "schema": SCHEMA,
            "claim_class": CLAIM_CLASS,
            "relationship_to_V014": "ADDITIVE_PROGRESS_SURFACE__V014_MEANING_UNCHANGED",
            "exact_results": {
                "record_gated_local_pair_response": {
                    "premise": "N=0;U_d=0;h>0;Delta>0;0<vartheta<2*pi",
                    "KEEP": "D=-8*h*x*I6-4*h*x*z^2*A_L",
                    "sectors": (
                        "D_A1=-8*h*x*(1+2*z^2)",
                        "D_E2=-8*h*x*(1-z^2)",
                        "D_T2=-8*h*x",
                    ),
                    "open_domain_rank": 6,
                    "matched_BREAK": "D_BREAK=0",
                },
                "interaction_owned_nonfactorization": (
                    "y-x*z^2=-(16/3)*h^3*U_d*(tau/hbar)^4+O((tau/hbar)^6)"
                ),
                "authenticated_finite_atlas": {
                    "shared_child": "J_e^1=J_f^1 iff c=c_prime",
                    "edge_distance": "d_G(m,n)=one_half*norm_1(n-m)=1",
                    "translation_cocycle": "sum_closed_path(delta)=0",
                    "frame_cocycle": "T_ln*T_nm=T_lm",
                    "premise": "INDEPENDENTLY_FORMED_NODE_ID_AND_PORT_LABEL_RECORDS",
                },
                "formation_pattern_source_threshold": {
                    "premise": "N=0;U_d=0;h>0;Delta>0;0<vartheta<2*pi",
                    "fixed_E_restriction": (
                        "ZERO_OR_ONE_FORMED_LINK_GIVES_ZERO__"
                        "EXACTLY_TWO_GIVES_RANK_ONE"
                    ),
                    "two_formed": "E^T*D^(ab)*E=-16*h*x*(1-z)*w_ab^T*w_ab",
                    "scope": "FIXED_E_X_E_RESTRICTION_NOT_BROKEN_S4_INVARIANT_BLOCK",
                },
                "matched_neighbor_response": {
                    "premise": (
                        "N=1;h=U_d=E0>0;Delta=6*E0;d_star=2;delta=0;"
                        "receiver_c_in_{1,2,3}"
                    ),
                    "reference": "MATCHED_TO_SOURCE_PATTERN_0000_WITH_RECEIVERS_FIXED_FORMED_KEEP",
                    "expansion_variable": "u=E0*tau/hbar",
                    "expansion": (
                        "Delta_mu_c^E=-(5626/42525)*u^12*kappa_c*w_0c+O(u^14)"
                    ),
                },
                "homogeneous_direct_edge_propagation": {
                    "premise": "GENERAL_N_DIRECT_EDGE;delta=0",
                    "raw": "Delta_Q_n6=-128*h^4*U_d^2*u_b_in",
                    "normalized_witness": "Delta_mu_n=(8/45)*u^6*u_b_in+O(u^8)",
                    "irrep": "u_b_in=one_half*one+one_half*s_b in A1_plus_T2",
                    "direct_E": "E_TRANSPOSE*x_q=0_AT_EVERY_ORDER",
                    "displayed_helper_E": (
                        "NONZERO_ON_INCOMPLETE_ENDPOINT_INCIDENCE__CANCELS_ON_COMPLETE_"
                        "HOMOGENEOUS_INTERIOR_ENDPOINT_INCIDENCE"
                    ),
                },
                "uniform_quasilocal_influence_envelope": {
                    "active_link_degree": "degree_L(e)=q(child)+2<=6",
                    "lambda_F3": "48*abs(U_d)/hbar",
                    "tail": (
                        "norm_commutator<=2*norm(A)*norm(B)*T_distance(lambda_F3*abs(t))"
                    ),
                    "exponential": (
                        "norm_commutator<=2*norm(A)*norm(B)*exp(-mu*(d_cell-v_mu*abs(t)))"
                    ),
                    "v_mu": "lambda_F3*exp(mu)/mu",
                    "v_1": "48*e*abs(U_d)/hbar_CELL_STEPS_PER_PARENT_TIME_UPPER_ENVELOPE_ONLY",
                },
                "authenticated_a3_bulk_dynamics": {
                    "premise": "SELECTED_HOMOGENEOUS_ALL_FORMED_F3_MEMBER_ON_A3_X_FOUR_PORTS",
                    "site_set": "L=A3_x_{1,2,3,4}",
                    "pair_interaction_degree": 6,
                    "lambda_F3": "48*abs(U_d)/hbar",
                    "boundary_comparison": (
                        "norm(tau_t^S(A)-tau_t^R(A))<=3*norm(A)*abs(X)*"
                        "sum_{r=R}^infinity(2r+1)^3*T_{r-r_X+1}(lambda_F3*abs(t))"
                    ),
                    "bulk_limit": "BOUNDARY_INDEPENDENT_STRONGLY_CONTINUOUS_QUASILOCAL_DYNAMICS",
                    "stationary_state": "AT_LEAST_ONE_JOINT_TIME_TRANSLATION_S4_INVARIANT_STATE_EXISTS_NOT_SELECTED",
                    "spectral_measure": (
                        "mu_AB(B)=mu_A1(B)*P_A1+mu_E(B)*P_E+mu_T2(B)*P_T2"
                    ),
                },
                "authenticated_finite_window_bulk_response": {
                    "premise": "FIXED_FINITE_AUTHENTICATED_CELLS_PORTS_PULSES_AND_READS_IN_A_CHOSEN_GL6AK_STATE",
                    "pulse_boundary": (
                        "norm(O_S(B)-O_R(B))<=E_R(B,t)+2*norm(B)*sum_k E_R(V_k,s_k)"
                    ),
                    "retarded_kernel": (
                        "G^R_{beta,alpha}(t)=i*E_star^2/(2*hbar)*Theta(t)*"
                        "omega([tau_t(M_beta),M_alpha])"
                    ),
                    "correlation_measure": "FINITE_WINDOW_LIOUVILLIAN_CORRELATION_MEASURE_POSITIVE",
                    "commutator_measure": "NOT_POSITIVE_WITHOUT_PASSIVITY_OR_KMS",
                    "formation_defect": "V_kappa=h*sum_{p_in_D}(1-kappa_p)*X_p",
                    "defect_scope": "BOUNDARY_INDEPENDENT_NONEQUILIBRIUM_CONTRAST_NOT_STATIONARY_SECTOR_COEFFICIENT",
                },
                "native_degree_lock_sector": {
                    "premise": (
                        "r_epsilon=-6*r_U;U_d>0;Delta=4*U_d*(d_star-2)>0;d_star>2"
                    ),
                    "strict_witness": "d_star=3;Delta=4*U_d",
                    "hamiltonian": "H=-h*sum_e X_e+U_d*sum_v(k_v-2)^2+C",
                    "linear_Ward_no_go": "[H,k_v-2]=-i*h*sum_{e_incident_v}Y_e",
                    "locked_pair_sector": "A1_FIXED__PAIR_VARIATIONS_EXACTLY_E__T2_ABSENT",
                    "uniform_locked_covariance": "connected_covariance_(A1,E,T2)=(0,8/3,0)",
                    "constraint_Gram": (
                        "TWO_GENERIC_FLAT_DIRECTIONS_PLUS_ONE_EIGENVALUE_QUADRATIC_"
                        "AT_TRIVIAL_CHARACTER__SINGULAR_VALUE_LINEAR"
                    ),
                    "declared_Q4_lower_orders": (
                        "H2=-(M/2)*(h^2/U_d)*P;H4=-(7M/24)*(h^4/U_d^3)*P;M=256"
                    ),
                    "first_off_diagonal": "-(63/8)*(h^6/U_d^5)_ALTERNATING_HEXAGON",
                },
                "complete_order6_locked_hamiltonian": {
                    "premise": "DECLARED_PERIOD4_GIRTH_AT_LEAST6_Q4;M=256;CONTROLLED_AS_h_over_U_d_TO_ZERO",
                    "formula": (
                        "H_eff=C*P-(M/2)*(h^2/U_d)*P-(7M/24)*(h^4/U_d^3)*P-"
                        "(h^6/U_d^5)*[(893M/1080)*P+(63/8)*sum_{c_in_Hex(Q4)}T_c]"
                        "+O(h^8/U_d^7)"
                    ),
                    "configuration_changing_term": "ONLY_ONE_ALTERNATING_HEXAGON_AT_ORDER6",
                    "hexagon_amplitude": "-63/8",
                    "order6_diagonal": "COMMON_SCALAR_-(893/1080)*M__NO_FLIPPABILITY_POTENTIAL",
                    "infinite_scope": "FORMAL_UNIFORMLY_FINITE_RANGE_LINKED_INTERACTION",
                },
                "locked_ir_representation_and_response_boundary": {
                    "trivial_character_link_kernel": "ker(B(1))=T2;dim=3",
                    "locked_pair_fluctuations": "E;dim=2",
                    "mismatch": "Hom_S4(T2,E)=0",
                    "uniform_pair_E_conservation": "BROKEN_BY_NATIVE_ALTERNATING_HEXAGON_AT_ORDER6",
                    "reciprocal_quadratic_inverse": (
                        "Gamma_E^R=[a0^R(omega)+c0^R(omega)*I2(theta)]*I2+"
                        "c2^R(omega)*T(Q_E(theta))+O(abs(theta)^4)"
                    ),
                    "phase_scope": "E_MASS_AND_CUBIC_ALLOWED__POLE_GAP_AND_DISPERSION_STATE_CONDITIONAL",
                },
                "authenticated_E_loop_selection_boundary": {
                    "premise": "c_in_E=ker(R);c_nonzero;LOCAL_k_x=2",
                    "pair_read": "O_x(c)=sum_{a<b}c_ab*Z_(x,a)*Z_(x,b)",
                    "uniform_locked_variance": "(8/3)*norm(c)^2",
                    "linear_K_obstruction": "P_lock*X_e*P_lock=0__FOUR_ONE_CELL_PORT_SOURCES_HAVE_NO_E",
                    "hexagon_E_change": "norm(delta_M)^2=16_AT_EACH_OF_THREE_PARENT_NODES",
                    "retained_support_gate": (
                        "matrix_element=-(63/8)*(h^6/U_d^5)*product_{e_in_C}kappa_e"
                    ),
                    "stationary_scope": "PRODUCT_TRACE_REFUTES_UNIVERSAL_NONZERO_DEFECT_RESPONSE",
                },
                "locked_hexagon_thermodynamic_sector": {
                    "premise": "PURE_ORDER6_DEGREE2_HEXAGON_INTERACTION;t=(63/8)*h^6/U_d^5",
                    "finite_component": "H_C=-t*A_C;Delta_C=t*(rho_C-lambda_2(C))",
                    "native_flux": "THREE_INDEPENDENT_COORDINATE_CUT_FLUXES_PLUS_DEPENDENT_FOURTH_PORT_COUNT",
                    "active_periodic_bound": "abs(E0(Q_L))/L^3>=t/64_FOR_4_DIVIDES_L__NO_LIMIT_CLAIM",
                    "quasilocal_interaction": "SUPPORT18;sup_e sum_{c:e_in_supp(tau_c)}norm(Phi(c))<=18*t",
                    "variance_bound": "Delta_L<=18*t*norm(w)_infinity^2*L/Var(F_L)",
                    "closure_scope": "FINITE_COMPONENT_GAP_CLOSURE_CONDITIONAL__SELECTED_GNS_BRIDGE_OPEN",
                },
                "native_hexagon_collective_response": {
                    "conserved_density": "CENTERED_PORT_T2",
                    "cycle_symbol": (
                        "C1(theta)*C1(theta)^T=4*(abs(theta)^2*P_T2-theta*theta^T)"
                    ),
                    "oscillator_strength": (
                        "f_u(chi)=(J*t_hex/2)*u_dagger*C(chi)*C(chi)^dagger*u"
                    ),
                    "single_mode_bound": "Delta_T2(chi;u)<=f_u(chi)/S_u_plus(chi)",
                    "harmonic_premise": "POSITIVE_COHERENT_K0_G0_HESSIANS_AND_SELECTED_SYMMETRIC_PHASE",
                    "conditional_isotropic_character_modes": (
                        "omega_1^2=omega_2^2=4*g*kappa*I2(theta)+O(abs(theta)^3)"
                    ),
                    "authenticated_pair_read": "O_x(c)=4*sum_{a<b}c_ab*e_a(x)*e_b(x)",
                    "pair_channel": "COMPLEMENT_EVEN__ONE_T2_ODD_OVERLAP_ZERO_IF_COMPLEMENT_SYMMETRIC",
                    "retained_source": "A1_PLUS_T2_AT_ZERO_CHARACTER__E_CROSS_PROJECTION_ZERO",
                    "algebraic_composite": "Sym2_0(T2)=E_PLUS_T2__NO_PHYSICAL_ROTATION_OR_STRESS_PROMOTION",
                },
                "order6_quantum_ice_crosswalk": {
                    "premise": "DISPLAYED_ORDER6_LOCKED_INTERACTION_AFTER_COMMON_SCALAR_REMOVAL",
                    "coupling": "g=(63/8)*h^6/U_d^5>0",
                    "diamond_embedding": "t_a=e_a-one/4;dot(t_a,t_b)=delta_ab-one/4",
                    "Hilbert_map": "DEGREE2_EQUALS_TWO_DIMERS_PER_DIAMOND_SITE_EQUALS_FULLY_PACKED_LOOPS",
                    "ring_operator": "T_c=P*product_{e_in_c}X_e*P",
                    "comparison_family": "H(v)=-g*sum_c T_c+v*sum_c F_c",
                    "exact_parameter": "v/g=0",
                    "distinct_RK_point": "v/g=1",
                    "order_scope": "NOT_AN_ALL_ORDERS_FINITE_h_over_U_d_IDENTIFICATION",
                },
                "vg0_first_character_static_closure": {
                    "premise": (
                        "PURE_ORDER6_v_over_g_ZERO;PERIODIC_Q_L;L>=4;"
                        "TRANSLATION_STABLE_CONNECTED_COMPONENT"
                    ),
                    "coupling": "J=(63/8)*h^6/U_d^5>0",
                    "first_character": "q_L=2*pi/L;z=(exp(i*q_L),1,1,1)",
                    "transverse_plane": "u=(0,u1,u2,u3);u1+u2+u3=0;norm(u)=1",
                    "cycle_norm": "norm(C(z)^dagger*u)^2=12*sin(pi/L)^2",
                    "oscillator_bound": "f_u(q_L)<=6*J*sin(pi/L)^2",
                    "component_gap_bound": (
                        "Delta_C(L)<=6*J*sin(pi/L)^2/S_u,L(q_L)"
                    ),
                    "static_exponent_premise": "S_u,L(q_L)>=s*L^(-alpha);s>0;alpha<2",
                    "static_exponent_closure": (
                        "Delta_C(L)<=(6*pi^2*J/s)*L^(alpha-2)->0"
                    ),
                    "quadrature_identity": "Var(F_c)+Var(F_s)=L^3*S_u,L(q_L)",
                    "alpha_one_consequence": (
                        "S>=s/L_IMPLIES_ONE_REAL_QUADRATURE_VARIANCE>=s*L^2/2_"
                        "AND_Delta_C=O(J/L)"
                    ),
                    "PF_transform": (
                        "P_C(n,m)=A_C(n,m)*psi(m)/(rho_C*psi(n));"
                        "pi_C=psi^2;Delta_C=J*rho_C*gap(P_C)"
                    ),
                    "closure_scope": "STATIC_LOWER_BOUND_AND_SELECTED_GNS_BRIDGE_UNPROVED",
                },
                "record_conditioned_collective_clock_and_typed_atlas": {
                    "premise": "ISOLATED_ORDER6_OFF_DIAGONAL_PURE_LOOP_COMPONENT",
                    "record_conditioned_loop": (
                        "H_hex(kappa)=-sum_c J_c(kappa)*tau_c;"
                        "J_c=J*product_{e_in_c}kappa_e;J=(63/8)*h^6/U_d^5"
                    ),
                    "nonuniform_ceiling": (
                        "LOWER_ORDER_DIAGONAL_SHIFTS_MAY_BE_CONFIGURATION_DEPENDENT__"
                        "NOT_COMPLETE_CONDITIONED_HAMILTONIAN"
                    ),
                    "homogeneous_formal_family": (
                        "H_hex(q)=q^6*H_hex(1);alpha_t^q=alpha_{q^6*t}^1"
                    ),
                    "fixed_state_spectrum": (
                        "L_q=q^6*L_1;P_q(B)=P_1(q^-6*B);"
                        "chi_q^R(omega)=q^-6*chi_1^R(omega/q^6);q>0"
                    ),
                    "time_dependent_formal_family": (
                        "sigma(t2,t1)=integral_{t1}^{t2}q(s)^6*ds"
                    ),
                    "orientation_coupling": (
                        "J_d=J*product_{a_not_equal_d}kappa_a^2;"
                        "delta_log_J_d=6*s-2*(r_c)_d"
                    ),
                    "tetrahedral_evaluation": (
                        "rank(E)=4;im(E)=A1_plus_T2;ker(E)=E"
                    ),
                    "typed_atlas": (
                        "(A1_plus_T2)_FORMAL_GENERATOR_COEFFICIENTS_PLUS_"
                        "E_AUTHENTICATED_PAIR_QUERY_ISOMORPHIC_TO_Sym2_R3_"
                        "BUT_NOT_ONE_PHYSICAL_METRIC_TANGENT"
                    ),
                    "conditional_bridge": (
                        "AV_PHASE_PLUS_AV_CONTINUUM_PLUS_AV_CLOCK_PLUS_"
                        "AV_CONSTITUTIVE_PLUS_AV_UPDATE"
                    ),
                    "exact_disposition": (
                        "FINITE_RETAINED_PATTERN_CONDITIONS_FUTURE_LEADING_LOOP_"
                        "GENERATOR__TWO_WAY_BACK_REACTION_OPEN"
                    ),
                },
                "anisotropic_folner_twist_closure": {
                    "premise": (
                        "PURE_ORDER6_RECTANGULAR_TORUS;L0,L2_ODD;L1_EVEN;"
                        "ALL_PERIODS>=4;CENTERED_PORT0_N0=V/2"
                    ),
                    "nonempty_sector": (
                        "ALTERNATING_PORT0_PORT3_MATCHING_UNION_CONSTANT_PORT1_MATCHING"
                    ),
                    "twist": "U0=exp(i*(2*pi/L1)*sum_x x1*n_(x,0))",
                    "translation_character": "Y*U0*Y^-1=-U0",
                    "finite_sector_dichotomy": (
                        "TRANSLATION_RELATED_EXACT_GROUND_COMPONENT_DEGENERACY_OR_"
                        "SAME_COMPONENT_ORTHOGONAL_EXCITATION"
                    ),
                    "gap_bound": (
                        "Delta_C<=2*J*V*(1-cos(2*pi/L1))"
                        "<=4*pi^2*J*L0*L2/L1"
                    ),
                    "folner_sequence": "(L0,L1,L2)=(m,2*m^3,m);m>=5_ODD",
                    "folner_closure": "Delta_C(m)<=2*pi^2*J/m->0_IN_STABLE_COMPONENT_BRANCH",
                    "closure_scope": (
                        "UNCONDITIONAL_FINITE_CENTERED_SECTOR_DICHOTOMY__"
                        "ISOTROPIC_SCALING_SELECTED_GNS_AND_HIGHER_ORDER_STABILITY_OPEN"
                    ),
                },
            },
            "controlled_evidence": {
                "quantum_ice_v_over_g_zero": {
                    "Shannon": "ZERO_TEMPERATURE_GFMC_PLUS_FINITE_ED__NUMERICAL_U1_LIQUID_PHASE_EVIDENCE",
                    "Benton": (
                        "QMC_CALIBRATED_GAUSSIAN_THEORY__omega(k)_approximately_c_abs(k)__"
                        "c=(0.6_plus_or_minus_0.1)*g*a0/hbar"
                    ),
                    "visibility": "EFFECTIVE_SINGLE_LINK_POLE_WEIGHT_PROPORTIONAL_TO_omega(k)",
                    "claim_class": "NUMERICAL_AND_EFFECTIVE_EVIDENCE_NOT_A_PHASE_GAP_OR_POLE_THEOREM",
                },
            },
            "open_gates": (
                "ISOTROPIC_V_OVER_G_ZERO_PHASE_CONTROL_AND_SELECTED_FINITE_TO_GNS_BRIDGE",
                "CONTROLLED_ALL_ORDERS_FINITE_H_OVER_U_D_SURVIVAL_OF_THE_ORDER6_LOCKED_PHASE",
                "PHYSICAL_PREPARATION_HOMOGENIZATION_AND_CAUSAL_UPDATE_OF_RETAINED_COEFFICIENT_FIELD",
                "NONZERO_SELECTED_PROPAGATION_AND_RETAINED_LINEAGE_VISIBLE_STATIONARY_RESPONSE",
                "CALIBRATED_PHYSICAL_MOMENTUM_LENGTH_TIME_SINGLE_CONE_AND_RANK6_CONSTITUTIVE_JOIN",
                "COMMON_COUPLING_COMPLETE_STRESS_WARD_CONTACT_OWNERSHIP_AND_RGRL_B",
                "NATIVE_CONTINUUM_OPERATOR_IDENTIFICATION_WITHOUT_RICCI_ANSATZ",
                "HELD_OUT_RICCI_EINSTEIN_COMPARISON_AFTER_NATIVE_OPERATOR_IS_EARNED",
                "MICROSCOPIC_RICCI_COEFFICIENT_AND_G_MODEL",
            ),
            "ceilings": {
                "physical_K_is_semantic_REC": False,
                "finite_atlas_is_autonomous_physical_space": False,
                "analytic_tail_is_exact_finite_speed": False,
                "quasilocal_envelope_is_Lorentz_or_common_physical_cone": False,
                "direct_A1_T2_signal_is_bulk_shear": False,
                "invariant_bulk_state_is_selected_physical_state": False,
                "positive_correlation_measure_implies_nonzero_retarded_response": False,
                "degree_lock_supplies_exact_Ward_identity": False,
                "order6_linked_hamiltonian_is_all_orders_parent": False,
                "finite_component_gap_closure_is_GNS_gaplessness": False,
                "conditional_T2_character_mode_is_physical_photon_or_cone": False,
                "pair_E_or_Sym2_T2_is_complete_stress_or_spin2": False,
                "quantum_ice_numerics_are_mathematical_phase_proof": False,
                "diamond_embedding_supplies_physical_scale": False,
                "static_exponent_premise_is_proved_from_first_principles": False,
                "anisotropic_finite_twist_closure_is_isotropic_or_GNS_gaplessness": False,
                "formal_global_or_fractional_q_is_authenticated_record_field": False,
                "homogeneous_q6_clock_rescaling_supplies_physical_cone_or_curvature": False,
                "typed_Sym2_source_read_atlas_is_rank6_metric_constitutive_law": False,
                "record_conditioned_future_generator_is_two_way_back_reaction": False,
                "microscopic_J_is_Newton_G": False,
                "native_operator_is_Ricci": False,
                "Einstein_equations_derived_here": False,
                "gravity_derived_here": False,
                "G_calculated_here": False,
                "empirical_confirmation": False,
            },
            "custody": {
                "packet_count": len(custody.packet_rows),
                "packets": custody.packet_rows,
                "declared_hash_rows_checked": custody.declared_hash_rows_checked,
                "verification": (
                    "PINNED_PRIMARY_HASHES__SEALED_MANIFESTS__INDEPENDENT_AUDIT_"
                    "DISPOSITIONS__TOP_PACKET_DECLARED_DEPENDENCY_ROWS"
                ),
            },
            "executable_scope": {
                "caller_arguments": 0,
                "physics_recalculated": False,
                "observation_loaded": False,
                "gravity_solver": False,
                "scientific_output": "PINNED_EXACT_PROGRESS_AND_OPEN_GATE_CERTIFICATE",
            },
        }
    )


@dataclass(frozen=True)
class GravityMicroscopicProgress:
    """Immutable handle produced by fresh sealed-packet custody verification."""

    _custody: _Custody

    @property
    def claim_class(self) -> str:
        return CLAIM_CLASS

    def certificate(self) -> Mapping[str, Any]:
        return _certificate(self._custody)


def gravity_microscopic_progress() -> GravityMicroscopicProgress:
    """Verify and expose sealed microscopic progress with no caller input."""
    return GravityMicroscopicProgress(_verify_custody())


def gravity_microscopic_progress_certificate() -> Mapping[str, Any]:
    """Return the immutable zero-input progress certificate."""
    return gravity_microscopic_progress().certificate()
