"""Zero-input URM checkpoint for sealed microscopic-gravity progress through GL6AI.

This additive surface verifies the frozen theorem/audit packets that establish
record-gated finite response, interaction-owned nonfactorization, an
authenticated finite relational atlas, homogeneous direct-edge propagation,
and a uniform quasi-local influence envelope.  It deliberately does not alter
the V014 working-theory certificate and does not promote these results to a
Lorentz cone, gravity, Einstein dynamics, or a numerical value of G.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path, PurePosixPath
import posixpath
import re
from types import MappingProxyType
from typing import Any, Mapping, NoReturn


SCHEMA = "WAC_GRAVITY_MICROSCOPIC_PROGRESS_CERTIFICATE_V001"
CLAIM_CLASS = "SEALED_MICROSCOPIC_F3_PROGRESS_THROUGH_GL6AI_WITH_STRICT_IR_CEILING"

_REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
_HASH_ROW = re.compile(r"^([0-9a-f]{64})  ([^\r\n]+)$")


@dataclass(frozen=True)
class _PacketPin:
    label: str
    author_dir: str
    theorem_sha256: str
    author_manifest_sha256: str
    author_seal_sha256: str | None
    audit_dir: str
    audit_sha256: str
    audit_manifest_sha256: str
    audit_seal_sha256: str
    disposition: str


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
        theorem_relative = f"{pin.author_dir}/THEOREM.md"
        author_manifest_relative = f"{pin.author_dir}/MANIFEST.sha256"
        audit_relative = f"{pin.audit_dir}/AUDIT.md"
        audit_manifest_relative = f"{pin.audit_dir}/MANIFEST.sha256"

        if _sha256(_root_path(theorem_relative)) != pin.theorem_sha256:
            _refuse(f"{pin.label} theorem hash mismatch")
        if _sha256(_root_path(author_manifest_relative)) != pin.author_manifest_sha256:
            _refuse(f"{pin.label} author-manifest hash mismatch")
        if _sha256(_root_path(audit_relative)) != pin.audit_sha256:
            _refuse(f"{pin.label} audit hash mismatch")
        if _sha256(_root_path(audit_manifest_relative)) != pin.audit_manifest_sha256:
            _refuse(f"{pin.label} audit-manifest hash mismatch")

        checked += _verify_hash_list(author_manifest_relative, visited)
        checked += _verify_hash_list(audit_manifest_relative, visited)
        if pin.author_seal_sha256 is not None:
            _verify_seal(pin.author_dir, pin.author_seal_sha256, pin.author_manifest_sha256)
        _verify_seal(pin.audit_dir, pin.audit_seal_sha256, pin.audit_manifest_sha256)

        audit_text = _strict_text(_root_path(audit_relative))
        marker = f"**Disposition:** `{pin.disposition}`"
        if marker not in audit_text:
            _refuse(f"{pin.label} accepted audit disposition is absent")
        rows.append(
            MappingProxyType(
                {
                    "gate": pin.label,
                    "author_directory": pin.author_dir,
                    "theorem_sha256": pin.theorem_sha256,
                    "author_manifest_sha256": pin.author_manifest_sha256,
                    "author_seal_sha256": pin.author_seal_sha256,
                    "audit_directory": pin.audit_dir,
                    "audit_sha256": pin.audit_sha256,
                    "audit_manifest_sha256": pin.audit_manifest_sha256,
                    "audit_seal_sha256": pin.audit_seal_sha256,
                    "audit_disposition": pin.disposition,
                }
            )
        )

    # Close the most exposed read-time race over the fourteen primary bytes.
    for pin in _PACKETS:
        if _sha256(_root_path(f"{pin.author_dir}/THEOREM.md")) != pin.theorem_sha256:
            _refuse(f"{pin.label} theorem changed during verification")
        if _sha256(_root_path(f"{pin.audit_dir}/AUDIT.md")) != pin.audit_sha256:
            _refuse(f"{pin.label} audit changed during verification")
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
            },
            "open_gates": (
                "STATIONARY_BULK_IR_LAW_AND_NONZERO_SOURCE_VISIBLE_LOW_FREQUENCY_WEIGHT",
                "COMMON_PHYSICAL_CONE_SHARED_BY_RECORD_MATTER_EM_AND_CLOCK_PROBES",
                "CONSERVATION_AND_COMPLETE_STRESS_WARD_CONTACT_OWNERSHIP",
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
