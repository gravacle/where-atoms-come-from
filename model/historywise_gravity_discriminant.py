"""Pinned URM exposure of the historywise-gravity formal discriminant.

This module exposes one sealed finite-group theorem lane.  It is deliberately
not a gravity solver, an outcome selector, a physical GARH-D admission engine,
or empirical evidence.  No caller input participates in custody or in the
scientific result.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
from types import MappingProxyType
from typing import Any, Mapping, NoReturn


SCHEMA = "WAC_HISTORYWISE_GRAVITY_DISCRIMINANT_CERTIFICATE_V001"
CLAIM_CLASS = "FORMAL_FINITE_GROUP_DISCRIMINANT_ONLY"
LANE_ID = "LANE_CROSS_RFT_MGFT_HISTORYWISE_GRAVITY_ACTUALIZATION_V001"
MANIFEST_SHA256 = "de2d0483821aa3983a401c8617c9684b115c6517c0d7493ced59acdee9ee53b8"
AUDIT_DISPOSITION = (
    "POST_REPAIR_INDEPENDENT_MATHEMATICAL_AND_PHYSICAL_SCOPE_AUDITS_ACCEPT__"
    "EXACT_FINITE_GROUP_DISCRIMINANT_ONLY__PHYSICAL_GARH_D_Q_AND_GR_OPEN"
)
INDEPENDENT_AUDIT_VERDICT = (
    "ACCEPT_EXACT_FINITE_GROUP_NONSELECTION_AND_STABILIZER_CRITERION__"
    "FINITE_WITNESSES_EXACT__PHYSICAL_GARH_D_GARH_Q_ACTUALIZATION_BORN_GR_"
    "AND_RECORD_CAUSED_GRAVITY_OPEN"
)
VERIFIER_TOTAL = "TOTAL 64/64 PASS"
VERIFIER_VERDICT = (
    "VERDICT FINITE_EQUIVARIANT_GRAVITY_NONSELECTION_AND_ORIENTING_INPUT_"
    "CRITERION_EXACT__PHYSICAL_GARH_D_OPEN"
)

_ARTIFACT_SHA256 = {
    "README.md": "6d701f71e334bd41eb35238b645dfd56c5a274b34ee80c899ca393fb5e36f3a0",
    "THEOREM.md": "94d14119aa1c5e01f7009801f869a38ddaf536b861c4df1ca9acf0c30a4340e3",
    "COUNTEREXAMPLES.md": "e6b1010d4c103c47f05c118a9fa02619b0e87bd2885cb7082de744a16df77b2d",
    "RESULT.md": "4877810bbaba54be25186fa71c57a7e9a3f1f87af8b366d096f8f533c9860d66",
    "AUDIT.md": "3e7870e18580705971b00538fb08ab8eef8033c100926cfd3441ce1a48e51b12",
    "PRIMARY_SOURCES.md": "a1c02c46115097f8717560b65b045e99c6dfcf34c5ec8935ac7f0ade6d869437",
    "URM_INTEGRATION.md": "fa875d7ed0f2271a2ece3a1e922a11f8681f0345673954faab373759982b6334",
    "verify_historywise_gravity_actualization.py": (
        "5175652d4d426753dd0ba6f9ddca8a9c692c90cdac5b387539b907833667cfbf"
    ),
    "VERIFICATION.txt": "06f52ccacbdf656b2a69c3b3adb2aea1da65dba58a24e0469d8466e252e9cd5f",
}
_MANIFEST_LINE = re.compile(r"^([0-9a-f]{64})  ([^\r\n]+)$")
_REPOSITORY_ROOT = Path(__file__).absolute().parent.parent
_LANE_ROOT = _REPOSITORY_ROOT / LANE_ID


class HistorywiseGravityDiscriminantRefusal(RuntimeError):
    """The pinned formal lane failed custody or exact reproduction."""


def _refuse(message: str) -> NoReturn:
    raise HistorywiseGravityDiscriminantRefusal(
        "HISTORYWISE GRAVITY DISCRIMINANT REFUSES: " + message
    )


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    try:
        return _sha256_bytes(path.read_bytes())
    except OSError as exc:
        _refuse(f"custody object is unreadable: {path.name}: {exc}")


def _strict_text(path: Path, label: str) -> str:
    try:
        payload = path.read_bytes()
        text = payload.decode("utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        _refuse(f"{label} is unreadable or not UTF-8: {exc}")
    if "\r" in text or "\x00" in text:
        _refuse(f"{label} contains forbidden carriage-return or NUL bytes")
    return text


def _parse_manifest(text: str) -> dict[str, str]:
    rows: dict[str, str] = {}
    if not text.endswith("\n"):
        _refuse("manifest is not newline terminated")
    for line in text.splitlines():
        match = _MANIFEST_LINE.fullmatch(line)
        if match is None:
            _refuse("manifest has a malformed row")
        digest, relative = match.groups()
        posix = PurePosixPath(relative)
        if posix.is_absolute() or ".." in posix.parts:
            _refuse("manifest contains an absolute or escaping path")
        expected_prefix = LANE_ID + "/"
        if not relative.startswith(expected_prefix):
            _refuse("manifest row is outside the pinned lane")
        name = relative[len(expected_prefix):]
        if not name or "/" in name or name in rows:
            _refuse("manifest contains a nested or duplicate artifact")
        rows[name] = digest
    if rows != _ARTIFACT_SHA256:
        _refuse("manifest is not the exact closed nine-artifact inventory")
    return rows


def _validate_audit_text(text: str) -> None:
    if AUDIT_DISPOSITION not in text:
        _refuse("accepted audit disposition is absent")
    final_marker = "## Final audit verdict\n\n" + INDEPENDENT_AUDIT_VERDICT
    if final_marker not in text:
        _refuse("independent final audit ACCEPT verdict is absent")
    if "Three independent post-repair reviews" not in text:
        _refuse("independent post-repair audit statement is absent")


def _sealed_fresh_stdout(verification_text: str) -> str:
    start = "EXECUTED CHECKS\n\n"
    end = "\n\nREPRODUCTION AND SCOPE CHECKS"
    if verification_text.count(start) != 1 or verification_text.count(end) != 1:
        _refuse("verification transcript sections are malformed")
    body = verification_text.split(start, 1)[1].split(end, 1)[0]
    if body.count(VERIFIER_TOTAL) != 1 or body.count(VERIFIER_VERDICT) != 1:
        _refuse("verification transcript lacks the exact 64/64 total or verdict")
    if INDEPENDENT_AUDIT_VERDICT not in verification_text:
        _refuse("verification transcript lacks the accepted independent verdict")
    ceiling = (
        "They do not machine-prove the general finite-group theorem, establish a physical "
        "orienting input, choose GARH-D or GARH-Q, derive a Born law, reconstruct a physical "
        "metric, derive general relativity, or prove that records cause gravity."
    )
    if ceiling not in verification_text:
        _refuse("verification transcript lacks the executable claim ceiling")
    return body + "\n"


def _run_fresh_verifier(lane_root: Path, expected_stdout: str) -> str:
    verifier = lane_root / "verify_historywise_gravity_actualization.py"
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment["PYTHONPYCACHEPREFIX"] = "/private/tmp/wac-historywise-gravity-pycache"
    try:
        result = subprocess.run(
            [sys.executable, "-I", "-B", str(verifier)],
            cwd=lane_root,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=30,
            check=False,
        )
    except (OSError, UnicodeError, subprocess.SubprocessError) as exc:
        _refuse(f"fresh sealed-verifier execution failed: {exc}")
    if result.returncode != 0:
        _refuse(f"fresh sealed verifier returned {result.returncode}")
    if result.stderr:
        _refuse("fresh sealed verifier emitted stderr")
    if result.stdout != expected_stdout:
        _refuse("fresh verifier output differs from the sealed 64/64 transcript and verdict")
    return _sha256_bytes(result.stdout.encode("utf-8"))


@dataclass(frozen=True)
class _Custody:
    manifest_sha256: str
    audit_sha256: str
    verification_sha256: str
    verifier_sha256: str
    fresh_verifier_stdout_sha256: str


def _verify_custody(lane_root: Path = _LANE_ROOT) -> _Custody:
    """Verify the pinned closed lane.  The parameter is private and test-only."""
    if lane_root.name != LANE_ID:
        _refuse("lane identity differs from the pinned lane")
    if not lane_root.is_dir() or lane_root.is_symlink():
        _refuse("pinned lane is absent, not a directory, or symlinked")

    expected_names = set(_ARTIFACT_SHA256) | {"MANIFEST.sha256"}
    try:
        children = tuple(lane_root.iterdir())
    except OSError as exc:
        _refuse(f"pinned lane cannot be enumerated: {exc}")
    if {child.name for child in children} != expected_names:
        _refuse("lane is not the exact closed manifest-plus-artifacts inventory")
    if any(child.is_symlink() for child in children):
        _refuse("lane contains a symlinked custody object")
    if any(not child.is_file() for child in children):
        _refuse("lane contains a non-file custody object")

    manifest_path = lane_root / "MANIFEST.sha256"
    if _sha256_file(manifest_path) != MANIFEST_SHA256:
        _refuse("closed manifest hash mismatches the pinned seal")
    manifest_text = _strict_text(manifest_path, "manifest")
    rows = _parse_manifest(manifest_text)
    for name, expected_digest in rows.items():
        artifact = lane_root / name
        if not artifact.is_file() or artifact.is_symlink():
            _refuse(f"manifest artifact is absent, non-file, or symlinked: {name}")
        if _sha256_file(artifact) != expected_digest:
            _refuse(f"manifest artifact hash mismatch: {name}")

    audit_text = _strict_text(lane_root / "AUDIT.md", "accepted audit")
    verification_text = _strict_text(
        lane_root / "VERIFICATION.txt", "verification transcript"
    )
    _validate_audit_text(audit_text)
    expected_stdout = _sealed_fresh_stdout(verification_text)
    fresh_digest = _run_fresh_verifier(lane_root, expected_stdout)

    # Close the narrow execution-time race over the three executable custody
    # objects without rewriting or trusting any broad model manifest.
    if _sha256_file(manifest_path) != MANIFEST_SHA256:
        _refuse("closed manifest changed during verification")
    for name in ("AUDIT.md", "VERIFICATION.txt", "verify_historywise_gravity_actualization.py"):
        if _sha256_file(lane_root / name) != _ARTIFACT_SHA256[name]:
            _refuse(f"custody object changed during verification: {name}")

    return _Custody(
        manifest_sha256=MANIFEST_SHA256,
        audit_sha256=_ARTIFACT_SHA256["AUDIT.md"],
        verification_sha256=_ARTIFACT_SHA256["VERIFICATION.txt"],
        verifier_sha256=_ARTIFACT_SHA256["verify_historywise_gravity_actualization.py"],
        fresh_verifier_stdout_sha256=fresh_digest,
    )


def _freeze(value: Any) -> Any:
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, tuple):
        return tuple(_freeze(item) for item in value)
    return value


def _certificate(custody: _Custody) -> Mapping[str, Any]:
    return _freeze(
        {
            "schema": SCHEMA,
            "claim_class": CLAIM_CLASS,
            "formal_results": {
                "HGA1": "EXACT_FINITE_GROUP_ORBIT_NONSELECTION",
                "HGA1_feedback": "EXACT_EQUIVARIANT_FIXED_POINT_ORBIT_CLOSURE",
                "HGA1a": "UNIQUE_EQUIVARIANT_FLOW_FIXED_INPUT_IMPLIES_FIXED_HISTORY",
                "HGA1b": "UNIQUE_FIXED_LAW_ON_FINITE_TRANSITIVE_ORBIT_IS_UNIFORM_NOT_A_SAMPLE",
                "HGA2": "CONDITIONAL_FINITE_AFFINE_MEAN_FIELD_NONSELECTION",
                "HGA3": "EXACT_TRANSITIVE_FINITE_G_SET_STABILIZER_CRITERION_UP_TO_CONJUGACY",
            },
            "discriminant": (
                "ENDOGENOUS_EQUIVARIANT_SINGLETON_EXCLUDED_ON_FIXED_POINT_FREE_ORBIT"
            ),
            "orienting_input": "K_SUBSET_L_MATHEMATICAL_CAPACITY_ONLY",
            "witnesses": {
                "negative": "EXACT_FINITE_TWO_CELL_ALGEBRAIC_PROXY_NOT_PHYSICAL_GR",
                "positive": "NOMINATED_ORIENTING_INPUT_CAPACITY_ONLY_NOT_A_PHYSICAL_FIELD",
            },
            "statuses": {
                "physical_GARH_D": "NOT_ADMITTED_BY_THIS_DISCRIMINANT",
                "GARH_Q": "NOT_DERIVED_NOT_FORCED_BY_THIS_DISCRIMINANT",
                "GARH_D_Q_DECISION": "NOT_MADE",
                "objective_actualization": "OPEN_IN_THIS_LANE",
                "physical_gravity": "NO_PROOF_OUTPUT",
                "record_causes_gravity": "NO_PROOF_OUTPUT",
                "Born_law": "NO_PROOF_OUTPUT",
                "general_relativity": "NO_PROOF_OUTPUT",
                "empirical_validation": "NONE",
                "caller_input_scientific_weight": "ZERO",
            },
            "nonpromotion": {
                "failed_GARH_D_promotes_GARH_Q": False,
                "positive_boundary_promotes_physical_orientation": False,
                "finite_metric_proxy_promotes_physical_metric_or_GR": False,
                "custody_or_execution_promotes_empirical_validation": False,
            },
            "authorizations": {
                "physical_GARH_D": False,
                "GARH_Q": False,
                "GARH_D_Q_decision": False,
                "objective_actualization": False,
                "physical_orienting_input": False,
                "physical_gravity": False,
                "record_causes_gravity": False,
                "Born_law": False,
                "general_relativity": False,
                "empirical_validation": False,
                "scientific_readiness": False,
                "gravity_solver": False,
                "outcome_selector": False,
            },
            "custody": {
                "lane_id": LANE_ID,
                "manifest_sha256": custody.manifest_sha256,
                "manifest_artifact_count": len(_ARTIFACT_SHA256),
                "audit_sha256": custody.audit_sha256,
                "audit_disposition": AUDIT_DISPOSITION,
                "independent_audit_verdict": INDEPENDENT_AUDIT_VERDICT,
                "verification_sha256": custody.verification_sha256,
                "verifier_sha256": custody.verifier_sha256,
                "fresh_verifier_stdout_sha256": custody.fresh_verifier_stdout_sha256,
                "sealed_verifier_total": VERIFIER_TOTAL,
                "sealed_verifier_verdict": VERIFIER_VERDICT,
            },
            "executable_scope": {
                "finite_witness_checks": "64/64_REPRODUCED",
                "general_finite_group_theorem_machine_proved": False,
                "statement": "EXECUTABLE_DOES_NOT_PROVE_GENERAL_FINITE_GROUP_THEOREM",
                "physical_or_empirical_proof_weight": "ZERO",
            },
        }
    )


@dataclass(frozen=True)
class HistorywiseGravityDiscriminant:
    """Immutable result of one fresh, pinned custody and reproduction check."""

    _custody: _Custody

    @property
    def claim_class(self) -> str:
        return CLAIM_CLASS

    @property
    def manifest_sha256(self) -> str:
        return self._custody.manifest_sha256

    def certificate(self) -> Mapping[str, Any]:
        """Return a fresh, recursively immutable formal-only certificate."""
        return _certificate(self._custody)


def historywise_gravity_discriminant() -> HistorywiseGravityDiscriminant:
    """Reproduce and expose the pinned formal discriminant with no caller input."""
    return HistorywiseGravityDiscriminant(_verify_custody())


def historywise_gravity_discriminant_certificate() -> Mapping[str, Any]:
    """Return the fresh immutable certificate; authorize no physical claim."""
    return historywise_gravity_discriminant().certificate()
