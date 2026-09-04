#!/usr/bin/env python3
"""Fail-closed verifier for the independent GL6CS hostile audit."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6CS_STRICT_LOCK_SIX_PAIR_SCALE_SEPARATION_V001"
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def unique_object(pairs):
    out = {}
    for key, value in pairs:
        if key in out:
            raise ValueError(f"duplicate JSON key: {key}")
        out[key] = value
    return out


required = {
    "README.md", "AUDIT_REPORT.md", "INDEPENDENT_RESULT.json",
    "TARGET.sha256", "CONTEXT.sha256", "VERIFICATION.txt",
    "verify_gl6cs_independent.py", "verify_packet.py",
    "MANIFEST.sha256", "SEAL.sha256",
}
for name in sorted(required):
    check((HERE / name).is_file(), f"required audit file {name}")

# Freeze and replay the final twelve author bytes.  Six science bytes were
# audited first; the author then added six custody bytes without changing the
# science bytes.
target_lines = [line for line in (HERE / "TARGET.sha256").read_text().splitlines()
                if line.strip()]
check(len(target_lines) == 12, "twelve final target bytes pinned")
target_names = set()
for line in target_lines:
    expected, relative = line.split("  ", 1)
    check(relative not in target_names, f"unique target path {relative}")
    target_names.add(relative)
    path = ROOT / relative
    check(path.is_file(), f"target exists {relative}")
    check(digest(path) == expected, f"target hash {relative}")

target_required = {
    "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md",
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "EXACT_LEDGER.json",
    "VERIFICATION.txt", "derive_strict_lock_scale_separation.py",
    "verify_packet.py", "MANIFEST.sha256", "SEAL.sha256",
}
check(target_names == {f"{TARGET.name}/{name}" for name in target_required},
      "target pins every final packet byte exactly once")
target_manifest_lines = [
    line for line in (TARGET / "MANIFEST.sha256").read_text().splitlines()
    if line.strip()
]
target_manifest_names = set()
for line in target_manifest_lines:
    expected, name = line.split("  ", 1)
    check(Path(name).parent == Path("."), f"target manifest local path {name}")
    check(name not in target_manifest_names, f"target manifest unique path {name}")
    target_manifest_names.add(name)
    check(digest(TARGET / name) == expected, f"target manifest hash {name}")
check(target_manifest_names == target_required - {"MANIFEST.sha256", "SEAL.sha256"},
      "target manifest covers all non-custody bytes")
target_seal = (TARGET / "SEAL.sha256").read_text().strip().split("  ", 1)
check(target_seal == [digest(TARGET / "MANIFEST.sha256"), "MANIFEST.sha256"],
      "target seal names and hashes target manifest")

# Pin the exact theorem/custody/audit context used in the scope assessment.
context_lines = [line for line in (HERE / "CONTEXT.sha256").read_text().splitlines()
                 if line.strip()]
check(len(context_lines) == 36, "thirty-six context bytes pinned")
context = {}
for line in context_lines:
    expected, relative = line.split("  ", 1)
    check(relative not in context, f"unique context path {relative}")
    context[relative] = expected
    path = ROOT / relative
    check(path.is_file(), f"context exists {relative}")
    check(digest(path) == expected, f"context hash {relative}")

# Every context seal must name and hash the manifest pinned beside it.
seal_paths = sorted(relative for relative in context if relative.endswith("/SEAL.sha256"))
check(len(seal_paths) == 12, "twelve upstream custody seals pinned")
for relative in seal_paths:
    directory = str(Path(relative).parent)
    manifest_relative = f"{directory}/MANIFEST.sha256"
    check(manifest_relative in context, f"manifest paired with {relative}")
    fields = (ROOT / relative).read_text().strip().split("  ", 1)
    check(fields[0] == context[manifest_relative] and
          fields[1] in {"MANIFEST.sha256", manifest_relative},
          f"valid pinned custody seal {relative}")

result = json.loads((HERE / "INDEPENDENT_RESULT.json").read_text(),
                    object_pairs_hook=unique_object)
check(result["schema"] == "AUDIT_G_GL6CS_STRICT_LOCK_ACCUMULATION_HORIZON_V001",
      "audit result schema")
check(result["target_directory"] == TARGET.name, "exact target directory")
check(result["disposition"] == "PASS", "pass disposition")
core = result["independent_results"]
check(core["eligible_locked_ring_incidences"] == 24, "24 local incidences")
check(core["every_ring_delta"] == "nonzero pure E2 with norm squared 16",
      "pure E2 local change")
check(core["solder_relation"] == "h_E/2=h_T", "solder normalization")
check(core["fixed_component_exponents"] == {
    "EE": -6, "ET": 0, "TT_contact": 2, "TT_writer": 6,
}, "four fixed-component exponents")
check(core["CS10_common_two_removed_prefactors"] == {
    "EE": "(8/63)a_E^2 r^-6/U_d",
    "ET": "(5/6)a_E/U_d",
    "TT_writer": "(175/32)r^6/U_d",
}, "CS10 factored coefficients")
check(core["full_spectral_prefactors"] == {
    "EE": "(16/63)a_E^2 r^-6/U_d",
    "ET": "(5/3)a_E/U_d",
    "TT_writer": "(175/16)r^6/U_d",
}, "full spectral coefficients")
check(core["TT_contact"] == "r^2/(4U_d)", "contact coefficient")
check(core["required_enhancement_without_orientation_mixing"] == {
    "contact": "O(r^-8)", "writer_susceptibility": "O(r^-12)",
}, "collective enhancement powers")
check("impossible as r->0" in core["fixed_frame_bounded_match"],
      "fixed-frame bounded obstruction")
scope = result["scope"]
check("fixed finite component" in scope["proved"], "proved scope")
for phrase in ("noncommutation", "Ricci/Einstein", "gravity", "G"):
    check(phrase in scope["not_proved"], f"open scope {phrase}")
check("REPAIR_REQUIRED" in scope["CP_status"], "CP not promoted")

independent = subprocess.run(
    [sys.executable, "-B", str(HERE / "verify_gl6cs_independent.py")],
    cwd=ROOT, capture_output=True, text=True, check=False,
)
check(independent.returncode == 0, "independent reconstruction exits zero")
for token in (
    "PASS__AUDIT_GL6CS_INDEPENDENT_RECONSTRUCTION__138/138",
    "DISPOSITION=PASS",
    "FIXED_COMPONENT_POWERS=EE_-6;ET_0;TT_CONTACT_2;TT_WRITER_6",
    "FIXED_FRAME_BOUNDED_MATCH=IMPOSSIBLE",
    "COLLECTIVE_ENHANCEMENT=CONTACT_R^-8;WRITER_R^-12",
    "NONCOMMUTING_LIMIT_OR_ORIENTATION_OR_FINITE_R_OR_NEW_OWNER=OPEN_ROUTES_NOT_DERIVED",
    "NO_RICCI_EINSTEIN_GRAVITY_G",
):
    check(token in independent.stdout, f"independent replay token {token}")

target_replay = subprocess.run(
    [sys.executable, "-B", str(TARGET / "derive_strict_lock_scale_separation.py")],
    cwd=ROOT, capture_output=True, text=True, check=False,
)
check(target_replay.returncode == 0, "frozen target replay exits zero")
for token in (
    "PASS__GL6CS_STRICT_LOCK_SCALE_SEPARATION__110/110",
    "LOCKED_RING_CHANGE=PURE_E2_NONZERO;24/24",
    "FIXED_COMPONENT_EXPONENTS=EE_-6;ET_0;TT_CONTACT_2;TT_WRITER_6",
    "FIXED_FRAME_STRICT_LOCK_ROTATIONAL_EQUALITY=IMPOSSIBLE_FOR_BOUNDED_COEFFICIENTS",
    "ESCAPES=AUTHENTICATED_ORIENTATION_MIXING;COLLECTIVE_SINGULAR_LIMIT;FINITE_R;NEW_SAME_ORDER_BLOCK",
):
    check(token in target_replay.stdout, f"frozen target replay token {token}")

target_packet = subprocess.run(
    [sys.executable, "-B", str(TARGET / "verify_packet.py")],
    cwd=ROOT, capture_output=True, text=True, check=False,
)
check(target_packet.returncode == 0, "final target packet exits zero")
check("PASS__GL6CS_PACKET__174/174" in target_packet.stdout,
      "final target packet replay token")

report = " ".join((HERE / "AUDIT_REPORT.md").read_text().split())
for token in (
    "Disposition: PASS",
    "24",
    "h_E/2=h_T",
    "fixed finite component",
    "r^{-8}",
    "r^{-12}",
    "does not prove that either iterated limit exists",
    "No continuum",
    "pins all twelve final target bytes",
    "174/174",
):
    check(token in report, f"report token {token}")

manifest_lines = [line for line in (HERE / "MANIFEST.sha256").read_text().splitlines()
                  if line.strip()]
manifest_names = set()
for line in manifest_lines:
    expected, name = line.split("  ", 1)
    check(Path(name).parent == Path("."), f"manifest local path {name}")
    check(name not in manifest_names, f"manifest unique path {name}")
    manifest_names.add(name)
    check(digest(HERE / name) == expected, f"manifest hash {name}")
check(manifest_names == required - {"MANIFEST.sha256", "SEAL.sha256"},
      "manifest covers all non-custody audit bytes")
seal = (HERE / "SEAL.sha256").read_text().strip().split("  ", 1)
check(seal == [digest(HERE / "MANIFEST.sha256"), "MANIFEST.sha256"],
      "seal hashes and names audit manifest")

print(f"PASS__AUDIT_GL6CS_PACKET__{checks}/{checks}")
