#!/usr/bin/env python3
"""Quarantine verifier for the existing LANE_T51_A package.

This program performs no new physical measurement and never imports or executes a
builder.  It freezes and hashes every input, re-adjudicates the existing independent
verifier's recorded data with new predicates, scans a closed claim-bearing corpus,
and audits whether the pre-registered V1--V5 rule is actually scoreable.

Default disposition is REFUTED.  Only predicates derived from the frozen inputs may
change a component to NOT_REFUTED.  The package-wide ALL-PASS assertion remains
REFUTED when a registered V-item is unscoreable or only descriptive.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any


OWNED_PREFIX = "LANE_T51_REPAIR/QUARANTINE_VERIFY/"
QUARANTINED_PREFIX = (
    "LANE_T51_REPAIR/QUARANTINED_DRAFTS/ORIGINAL_LAYOUT/LANE_T51_A/"
)

# Every file read or relied on by this verifier is frozen here.  The builder manifest
# hashes are repeated deliberately: integrity does not depend on trusting the manifest
# that is itself under audit.
FROZEN_INPUTS: dict[str, str] = {
    "FIELD_INSTRUMENT_V001.md": "e0eecbfb1aa9ef6f934e07fa582a220c2a59435535bbbe3949baf0bca07e7407",
    "LANE_T51_A.sha256": "15a4d1a83fe92b4282746a35b996c179013a94c942e12e609f7322429bfca1a8",
    "LANE_T51_A/D24_AUDIT.txt": "c8faf8f26aad5bd84d2dda8767d3e9eb9e8bb2c1c1c3b9f207e91ae5dd17a823",
    "LANE_T51_A/VERIFY/v1_adversarial_rebuild.py": "1b715b75edac5941f114d870d6cfd8a5a9a828ed3d9b96683f6797d8ec350628",
    "LANE_T51_A/a1_calibration_3x2.json": "3b5b74a64f7904bf42e8d116830fd729d3441053744b172c3c48075b731c4729",
    "LANE_T51_A/a1_calibration_3x2.py": "566da7992bcba25b0bc902897a0c1f13644920a586b57f77b7905fe0d57efe38",
    "LANE_T51_A/a1_calibration_3x2.txt": "56ff0332485ac019872a03ec7a28e4e6b8b7ac4efc4a5db6732b313e6a3164dc",
    "LANE_T51_A/a2_field_table_3x3.json": "d72805a2fb1aadf0964a895e4579f9819daa34f6fb7e1a55d55bd651edc4c872",
    "LANE_T51_A/a2_field_table_3x3.py": "253615e2a3a5a4e1b9d3be28adb387efd962db44a19f91be769942f7cc761d6c",
    "LANE_T51_A/a2_field_table_3x3.txt": "a7fbf12026592a550e5b320e489256d7221771f6e1ff060a9da9beb75a99970f",
    "LANE_T51_A/a3_controls_3x3.json": "6581d4d5f098f953c4ae26f1b7d5240477557adc8ab033416c5464ef27ad865e",
    "LANE_T51_A/a3_controls_3x3.py": "43e299638a06d97532c6c71aadeb10dfe9a9ae3953c4d98a793fecba2e62a796",
    "LANE_T51_A/a3_controls_3x3.txt": "3cb8612d4aae5d4f052d00e9b2eb803939805c87362965f69f70b092c060adfc",
    "LANE_T51_A/a4_verdicts.py": "f3cd4d05a5d51fcfa51b67ffa74b6b5f54a00e00b0a08d3ecd34d91a1556f42a",
    "LANE_T51_A/a4_verdicts.txt": "006b045775fce44fe99f6d5cb7cd3dbde63230053712ce3446fd6d42102372d4",
    "LANE_T51_A/g1_connected_wenc.py": "255f1f97564359d9ed5eb0eedd822723616b2f5239537233ef8396746184fcd1",
    "LANE_T51_A/g1_connected_wenc.txt": "5d138e87a68401f1963ec820cb4607ad0d4a4e82c67b3b58de51041c7f8e547b",
    "LANE_T51_A/t51a_lib.py": "f2cc435a755e219990010b02114594fe5500e5d0f9aef3d9bb680835964f99a5",
    QUARANTINED_PREFIX + "VERIFY_CODEX/v1_independent_rebuild.py": "92c091520eb5b6da2ee356f7aab8731932441bb57eb59e2f9bcebc7ae032439e",
    QUARANTINED_PREFIX + "VERIFY_CODEX/v1_independent_rebuild.json": "8c95720cc9f721fea6ce06e56eb1cbcf77e341db72186cf830726b1b8302144a",
    QUARANTINED_PREFIX + "VERIFY_CODEX/v1_independent_rebuild.txt": "bad804807922ce63cffb82752c2e173d967d5ea0c0d01181e696cc03fe528198",
    QUARANTINED_PREFIX + "VERIFY_CODEX/D24_AUDIT.txt": "b7c8c9caa5228281f333df8b9534707b93bb6cd3f3f0d793bc28ac454e6a9d29",
    QUARANTINED_PREFIX + "VERIFY_CODEX/SEALS.sha256": "4f4f01a957d985d941ec468baf053354d27bd651eab295366f542005fb81866c",
    QUARANTINED_PREFIX + "JUDGMENT_CODEX.txt": "fd049d5538c3acc0cd5a15dcf062e0ae58b603470a805638b6ff793bb130c42f",
    "LANE_T51_REPAIR/QUARANTINED_DRAFTS/README.md": "e5dba0b75d005692be4498e93f4a119d5a1ddccda5d0fc6301f5027942acb6c7",
    "LANE_T51_REPAIR/QUARANTINED_DRAFTS/INPUTS.sha256": "f9777b818df794f27fab89418963a01e49384ea748faed0d6b813e38326d6d24",
}

# E is a closed corpus, not a directory walk.  It includes nested verifier outputs
# and the unlanded judgment that the earlier top-level-only scan missed.  Source files
# are included because several carry generated claim text and decision logic.
E_CLAIM_CORPUS: tuple[str, ...] = tuple(
    rel
    for rel in FROZEN_INPUTS
    if rel
    not in {
        "LANE_T51_A.sha256",
        QUARANTINED_PREFIX + "VERIFY_CODEX/SEALS.sha256",
        "LANE_T51_REPAIR/QUARANTINED_DRAFTS/INPUTS.sha256",
    }
)

SELF_REFERENCE_EXCLUSIONS = {
    "owned_verifier_tree": OWNED_PREFIX + "**",
    "reason": (
        "the verifier source contains rule definitions and deliberate positive-control "
        "violations; its generated report names refutations.  Scanning it would test the "
        "scanner against itself rather than scan the frozen T51 package"
    ),
    "nonclaim_integrity_metadata": [
        "LANE_T51_A.sha256",
        QUARANTINED_PREFIX + "VERIFY_CODEX/SEALS.sha256",
        "LANE_T51_REPAIR/QUARANTINED_DRAFTS/INPUTS.sha256",
        "*.sha256 sidecars",
    ],
    "commissioning_briefs": (
        "read as procedure but not part of the claim corpus; they prohibit imported "
        "standards and contain quoted attack vocabulary"
    ),
}

CHECK_LAMBDAS = (0.037, 0.071)
CHECK_PLACEMENTS = ("near", "far", "fresh")
EXPECTED_CHARACTERS = {
    "near": ([0, 0], -1, False),
    "far": ([1, 1], -1, True),
    "fresh": ([1, 1], -1, True),
}
BRACKET_HALF_WIDTH = 0.25
FIT_SIGNAL_MULTIPLE = 100.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalized(text: str) -> str:
    return " ".join(text.split())


def normalized_claim_text(text: str) -> str:
    """Normalize Markdown claim prose without changing its words."""
    return normalized(text.replace("`", "").replace("*", ""))


def verify_frozen_inputs(repo: Path) -> dict[str, Any]:
    rows = []
    for rel, expected in FROZEN_INPUTS.items():
        path = repo / rel
        exists = path.is_file()
        actual = sha256(path) if exists else None
        rows.append(
            {
                "path": rel,
                "expected_sha256": expected,
                "actual_sha256": actual,
                "ok": exists and actual == expected,
            }
        )
    mismatches = [row for row in rows if not row["ok"]]

    # Cross-check the two historical manifests without trusting either one as the
    # source of expected hashes.
    builder_lines = (repo / "LANE_T51_A.sha256").read_text(encoding="utf-8").splitlines()
    builder_manifest = {}
    for line in builder_lines:
        if line.strip():
            digest, rel = line.split(None, 1)
            builder_manifest[rel.strip()] = digest
    builder_expected = {
        rel: digest
        for rel, digest in FROZEN_INPUTS.items()
        if rel.startswith("LANE_T51_A/")
    }
    builder_manifest_ok = builder_manifest == builder_expected

    seal_lines = (repo / (QUARANTINED_PREFIX + "VERIFY_CODEX/SEALS.sha256")).read_text(
        encoding="utf-8"
    ).splitlines()
    verifier_manifest = {}
    for line in seal_lines:
        if line.strip():
            digest, name = line.split(None, 1)
            verifier_manifest[name.strip()] = digest
    verifier_expected = {
        Path(rel).name: digest
        for rel, digest in FROZEN_INPUTS.items()
        if rel.startswith(QUARANTINED_PREFIX + "VERIFY_CODEX/")
        and not rel.endswith("SEALS.sha256")
    }
    verifier_manifest_ok = verifier_manifest == verifier_expected

    quarantine_lines = (
        repo / "LANE_T51_REPAIR/QUARANTINED_DRAFTS/INPUTS.sha256"
    ).read_text(encoding="utf-8").splitlines()
    quarantine_manifest = {}
    for line in quarantine_lines:
        if line.strip():
            digest, rel = line.split(None, 1)
            quarantine_manifest[rel.strip()] = digest
    quarantine_expected = {
        rel: digest
        for rel, digest in FROZEN_INPUTS.items()
        if rel.startswith(QUARANTINED_PREFIX)
    }
    quarantine_manifest_ok = quarantine_manifest == quarantine_expected

    return {
        "default": "REFUTED",
        "count": len(rows),
        "rows": rows,
        "mismatches": mismatches,
        "builder_manifest_matches_frozen_map": builder_manifest_ok,
        "prior_verifier_manifest_matches_frozen_map": verifier_manifest_ok,
        "quarantine_manifest_matches_frozen_map": quarantine_manifest_ok,
        "verdict": (
            "NOT_REFUTED"
            if not mismatches
            and builder_manifest_ok
            and verifier_manifest_ok
            and quarantine_manifest_ok
            else "REFUTED"
        ),
    }


def line_evidence(text: str, needles: tuple[str, ...]) -> list[dict[str, Any]]:
    evidence = []
    lines = text.splitlines()
    for needle in needles:
        for index, line in enumerate(lines, 1):
            if needle in line:
                evidence.append({"needle": needle, "line": index, "text": line.strip()})
    return evidence


def preserve_old_defects(repo: Path) -> dict[str, Any]:
    old_source = (repo / (QUARANTINED_PREFIX + "VERIFY_CODEX/v1_independent_rebuild.py")).read_text(
        encoding="utf-8"
    )
    old_lines = old_source.splitlines()
    c_start = next(
        index
        for index, line in enumerate(old_lines)
        if '"C_winding_attribution": "NOT_REFUTED" if all(' in line
    )
    c_evidence = [
        {"line": index + 1, "text": old_lines[index].strip()}
        for index in range(c_start, c_start + 3)
    ]
    e_evidence = line_evidence(
        old_source,
        (
            'for p in (repo / "LANE_T51_A").iterdir()',
            '"violation_count": 0,',
            '"E_D1_directive_scan": "NOT_REFUTED" if prose["violation_count"] == 0',
        ),
    )
    c_segment = "\n".join(row["text"] for row in c_evidence)
    c_defect = all(
        needle in c_segment
        for needle in (
            '"C_winding_attribution": "NOT_REFUTED" if all(',
            'dynamics[p]["winding"][str(CHECK_LAMBDA)]',
            'for p in ("near", "far")',
        )
    )
    e_defect = all(any(row["needle"] == needle for row in e_evidence) for needle in (
        'for p in (repo / "LANE_T51_A").iterdir()',
        '"violation_count": 0,',
        '"E_D1_directive_scan": "NOT_REFUTED" if prose["violation_count"] == 0',
    ))
    return {
        "old_source_sha256": sha256(
            repo / (QUARANTINED_PREFIX + "VERIFY_CODEX/v1_independent_rebuild.py")
        ),
        "C_defect": {
            "confirmed": c_defect,
            "finding": (
                "the old verdict consumed only CHECK_LAMBDA=0.037 and only near/far; "
                "it computed but did not predicate on 0.071 or fresh placement"
            ),
            "source_evidence": c_evidence,
        },
        "E_defect": {
            "confirmed": e_defect,
            "finding": (
                "the old scan used a top-level iterdir, explicitly missed nested artifacts, "
                "returned literal violation_count=0, and based E on that literal"
            ),
            "source_evidence": e_evidence,
        },
    }


def recompute_character(rows: dict[str, Any], signal_threshold: float) -> dict[str, Any]:
    required = ("00", "01", "10", "11")
    keys_ok = tuple(sorted(rows)) == tuple(sorted(required))
    finite = keys_ok and all(math.isfinite(float(rows[key]["F"])) for key in required)
    clears = finite and all(abs(float(rows[key]["F"])) >= signal_threshold for key in required)
    if not clears:
        return {
            "keys_ok": keys_ok,
            "finite": finite,
            "clears_threshold": clears,
            "fits": [],
            "unique": False,
        }
    signs = {key: 1 if float(rows[key]["F"]) > 0 else -1 for key in required}
    fits = []
    for c1 in (0, 1):
        for c2 in (0, 1):
            s0 = signs["00"]
            ok = all(
                signs[f"{w1}{w2}"]
                == s0 * ((-1) ** ((c1 * w1 + c2 * w2) & 1))
                for w1 in (0, 1)
                for w2 in (0, 1)
            )
            if ok:
                fits.append({"character": [c1, c2], "content_sign": s0})
    return {
        "keys_ok": keys_ok,
        "finite": finite,
        "clears_threshold": clears,
        "signs": signs,
        "fits": fits,
        "unique": len(fits) == 1,
        "sign_flip_by_winding_alone": len(set(signs.values())) > 1,
    }


def verify_c(prior: dict[str, Any]) -> dict[str, Any]:
    floor = float(prior["C_and_D"]["noise"]["declared_floor"])
    threshold = FIT_SIGNAL_MULTIPLE * floor
    dynamics = prior["C_and_D"]["dynamics"]
    checks = []
    by_placement: dict[str, list[dict[str, Any]]] = {name: [] for name in CHECK_PLACEMENTS}
    for placement in CHECK_PLACEMENTS:
        expected_character, expected_s0, expected_flip = EXPECTED_CHARACTERS[placement]
        for lam in CHECK_LAMBDAS:
            key = str(lam)
            entry = dynamics.get(placement, {}).get("winding", {}).get(key)
            if entry is None:
                row = {
                    "placement": placement,
                    "lambda": lam,
                    "present": False,
                    "ok": False,
                }
            else:
                rebuilt = recompute_character(entry["rows"], threshold)
                stored = entry.get("character", {})
                expected_fit = {
                    "character": expected_character,
                    "content_sign": expected_s0,
                }
                stored_matches = (
                    stored.get("fits") == rebuilt.get("fits")
                    and bool(stored.get("separated")) == bool(rebuilt.get("unique"))
                )
                row = {
                    "placement": placement,
                    "lambda": lam,
                    "present": True,
                    "recomputed": rebuilt,
                    "stored_matches_recomputation": stored_matches,
                    "expected_fit": expected_fit,
                    "expected_flip": expected_flip,
                    "ok": (
                        rebuilt.get("unique") is True
                        and rebuilt.get("fits") == [expected_fit]
                        and rebuilt.get("sign_flip_by_winding_alone") == expected_flip
                        and stored_matches
                    ),
                }
            checks.append(row)
            by_placement[placement].append(row)

    stable = {}
    for placement, rows in by_placement.items():
        fits = [row.get("recomputed", {}).get("fits") for row in rows]
        stable[placement] = len(fits) == len(CHECK_LAMBDAS) and fits[0] == fits[1]
    all_ok = all(row["ok"] for row in checks) and all(stable.values())
    return {
        "default": "REFUTED",
        "noise_floor": floor,
        "signal_threshold": threshold,
        "required_lambdas": list(CHECK_LAMBDAS),
        "required_placements": list(CHECK_PLACEMENTS),
        "checks": checks,
        "lambda_stability": stable,
        "verdict": "NOT_REFUTED" if all_ok else "REFUTED",
    }


def solve_linear(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    aug = [list(row) + [rhs[i]] for i, row in enumerate(matrix)]
    n = len(rhs)
    for col in range(n):
        pivot = max(range(col, n), key=lambda row: abs(aug[row][col]))
        if abs(aug[pivot][col]) < 1e-18:
            raise ValueError("singular normal equations")
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [value / scale for value in aug[col]]
        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            aug[row] = [
                aug[row][j] - factor * aug[col][j]
                for j in range(n + 1)
            ]
    return [aug[i][-1] for i in range(n)]


def refit_onset(rows: list[dict[str, Any]]) -> float:
    design = [[1.0, math.log(float(row["lam"])), float(row["lam"]) ** 2] for row in rows]
    target = [math.log(abs(float(row["F"]))) for row in rows]
    normal = [
        [sum(a[i] * a[j] for a in design) for j in range(3)]
        for i in range(3)
    ]
    rhs = [sum(a[i] * y for a, y in zip(design, target)) for i in range(3)]
    return solve_linear(normal, rhs)[1]


def verify_a_b_d(prior: dict[str, Any]) -> dict[str, Any]:
    a = prior["A"]
    a_conditions = {
        "star_source_pairing_zero": a.get("star_source_symplectic_max") == 0,
        "mediator_source_pairing_zero": a.get("mediator_source_symplectic_max") == 0,
        "full_transition_audit_nonempty": int(a.get("full_venue_transitions", 0)) > 0,
        "source_changing_transitions_zero": a.get("source_changing_transitions") == 0,
    }
    a_verdict = "NOT_REFUTED" if all(a_conditions.values()) else "REFUTED"

    b = prior["B"]
    near, far, fresh = b["near_rebuild"], b["far_rebuild"], b["fresh_rebuild"]
    b_conditions = {
        "near_connected_minimum_3": near.get("conn_min") == 3,
        "far_old_minimum_4": far.get("old_min") == 4,
        "far_old_witness_disconnected": far.get("old_witness_connected") is False,
        "far_connected_minimum_5": far.get("conn_min") == 5,
        "far_exhaustive_counts": far.get("connector_count") == 1024 and far.get("enclosing_count") == 512,
        "far_no_lighter_connected": far.get("lighter_connected_histogram") == {},
        "fresh_connected_minimum_5": fresh.get("conn_min") == 5,
        "fresh_no_lighter_connected": fresh.get("lighter_connected_histogram") == {},
    }
    b_verdict = "NOT_REFUTED" if all(b_conditions.values()) else "REFUTED"

    expected_min = {"near": 3, "far": 5, "fresh": 5}
    d_rows = []
    dynamics = prior["C_and_D"]["dynamics"]
    for placement in CHECK_PLACEMENTS:
        rows = dynamics[placement]["offgrid"]
        lams = [float(row["lam"]) for row in rows]
        finite = all(math.isfinite(float(row["F"])) and float(row["F"]) != 0.0 for row in rows)
        k = refit_onset(rows) if finite else math.nan
        stored_k = float(dynamics[placement]["fit"]["k"])
        row = {
            "placement": placement,
            "lambdas": lams,
            "finite_nonzero": finite,
            "recomputed_k": k,
            "stored_k": stored_k,
            "stored_matches": finite and abs(k - stored_k) <= 1e-9,
            "connected_minimum": expected_min[placement],
            "bracket_pass": finite and abs(k - expected_min[placement]) <= BRACKET_HALF_WIDTH,
        }
        row["ok"] = row["stored_matches"] and row["bracket_pass"]
        d_rows.append(row)
    d_verdict = "NOT_REFUTED" if all(row["ok"] for row in d_rows) else "REFUTED"

    return {
        "A": {"default": "REFUTED", "conditions": a_conditions, "verdict": a_verdict},
        "B": {"default": "REFUTED", "conditions": b_conditions, "verdict": b_verdict},
        "D": {
            "default": "REFUTED",
            "method": "pure-stdlib re-fit of recorded off-grid F values; no Hamiltonian run",
            "rows": d_rows,
            "verdict": d_verdict,
        },
    }


E_PATTERNS = {
    "imported_shape_prescription": re.compile(
        r"\b(?:field|reading|response|law|result)\b.{0,120}"
        r"\b(?:must|should|needs?\s+to|has\s+to|is\s+required\s+to|is\s+expected\s+to)\b"
        r".{0,120}\b(?:newton(?:ian)?|inverse[- ]square|geodesic|classical\s+gravity|gravitational\s+form)\b",
        re.IGNORECASE,
    ),
    "imported_pass_requirement": re.compile(
        r"\b(?:field|result|reading|response)\b.{0,100}"
        r"\b(?:is\s+valid\s+only\s+if|passes\s+only\s+if|must|should|is\s+required\s+to)\b"
        r".{0,100}\b(?:match|follow|recover|reproduce)\w*\b.{0,100}"
        r"\b(?:classical|external|gravity|metric|falloff|shape|accumulat\w*)\b",
        re.IGNORECASE,
    ),
    "outcome_as_failure": re.compile(
        r"\b(?:screening|screened|saturation|saturates?|mod[- ]?2|superposition|composition|outcome)\b"
        r".{0,100}\b(?:would\s+be|means|counts\s+as|is|are)\b.{0,60}"
        r"\b(?:a\s+)?(?:failure|failed|fail|dead|kill(?:ed)?)\b",
        re.IGNORECASE,
    ),
}

POSITIVE_CONTROL_SEEDS = {
    "imported_shape_prescription": "The field reading must obey an inverse-square law.",
    "imported_pass_requirement": (
        "A field result is valid only if it matches the classical accumulation standard."
    ),
    "outcome_as_failure": "Screening would be a failure outcome.",
}


def exculpation_reason(text: str) -> str | None:
    lowered = normalized(text).lower()
    exemptions = {
        "explicit_no_requirement": (
            "no classical gravitational form is required",
            "not imported here as a requirement",
            "no classical gravitational form was used",
            "no classical gravitational form entered",
        ),
        "explicit_no_failure": (
            "never a failure against",
            "no outcome here is a failure",
            "a result, not a failure",
        ),
        "scanner_or_policy_description": (
            "prohibited-shape hits",
            "outcome-as-failure hits",
            "frames a computed outcome as failure",
            "framing of any computed accumulation",
        ),
    }
    for reason, phrases in exemptions.items():
        if any(phrase in lowered for phrase in phrases):
            return reason
    return None


def scan_claim_corpus(repo: Path) -> dict[str, Any]:
    assert all(not rel.startswith(OWNED_PREFIX) for rel in E_CLAIM_CORPUS)
    positive_controls = {}
    for rule, seed in POSITIVE_CONTROL_SEEDS.items():
        matched = bool(E_PATTERNS[rule].search(seed))
        positive_controls[rule] = {"seed": seed, "matched": matched}

    candidate_hits = []
    decision_hits = []
    for rel in E_CLAIM_CORPUS:
        text = (repo / rel).read_text(encoding="utf-8")
        lines = text.splitlines()
        seen: set[tuple[str, int, str]] = set()
        # Five-line bidirectional windows let wrapped Markdown or generated prose remain
        # a sentence.  Looking backward is load-bearing: a forward-only first draft saw
        # "outcome ... failure" at a wrapped continuation but omitted the preceding "no".
        for index in range(len(lines)):
            window = normalized(" ".join(lines[max(0, index - 2) : index + 3]))
            if not window:
                continue
            for rule, pattern in E_PATTERNS.items():
                match = pattern.search(window)
                if not match:
                    continue
                key = (rule, rel, match.group(0))
                if key in seen:
                    continue
                seen.add(key)
                reason = exculpation_reason(window)
                hit = {
                    "rule": rule,
                    "file": rel,
                    "line": index + 1,
                    "match": match.group(0),
                    "context": window,
                    "exculpated": reason is not None,
                    "exculpation_reason": reason,
                }
                candidate_hits.append(hit)
                if reason is None:
                    decision_hits.append(hit)

    controls_ok = all(row["matched"] for row in positive_controls.values())
    verdict = "NOT_REFUTED" if controls_ok and not decision_hits else "REFUTED"
    return {
        "default": "REFUTED",
        "corpus_count": len(E_CLAIM_CORPUS),
        "corpus": [
            {"path": rel, "sha256": FROZEN_INPUTS[rel]}
            for rel in E_CLAIM_CORPUS
        ],
        "self_reference_exclusions": SELF_REFERENCE_EXCLUSIONS,
        "positive_controls": positive_controls,
        "positive_controls_all_matched": controls_ok,
        "candidate_hits": candidate_hits,
        "decision_hits": decision_hits,
        "decision_hit_count": len(decision_hits),
        "verdict": verdict,
    }


def score_admissibility(repo: Path, prior: dict[str, Any], c_result: dict[str, Any]) -> dict[str, Any]:
    field = (repo / "FIELD_INSTRUMENT_V001.md").read_text(encoding="utf-8")
    field_claim_text = normalized_claim_text(field)
    a2 = json.loads((repo / "LANE_T51_A/a2_field_table_3x3.json").read_text(encoding="utf-8"))
    a3 = json.loads((repo / "LANE_T51_A/a3_controls_3x3.json").read_text(encoding="utf-8"))
    a3_text = (repo / "LANE_T51_A/a3_controls_3x3.txt").read_text(encoding="utf-8")
    lib_source = (repo / "LANE_T51_A/t51a_lib.py").read_text(encoding="utf-8")
    audit = normalized((repo / "LANE_T51_A/D24_AUDIT.txt").read_text(encoding="utf-8"))

    prereg_quotes = {
        "V1": "|F| exceeds the same-table measured control floor beyond contact at the larger separation.",
        "V2": "the reading follows earned geometry under the placement swap.",
        "V3": "the onset-order bracket contains connected w_enc at both placements",
        "V4": "back-action below the declared tolerance.",
        "V5": "sign attribution licensed only if the winding-sector sweep separates content sign from winding sign.",
    }
    quote_presence = {key: quote in field_claim_text for key, quote in prereg_quotes.items()}

    far_rows = a2["ftab"][a2["p_far"]]["rows"]
    min_far_abs_f = min(abs(float(row["F"])) for row in far_rows)
    control_floor = float(a2["control_floor"])
    v1_true = quote_presence["V1"] and min_far_abs_f > control_floor

    c2_pairs = a3["c2"]
    pair_details = []
    for label, block in c2_pairs.items():
        names = block["pair"]
        f_values = [float(block["res"][name]["F05"]) for name in names]
        k_values = [float(block["res"][name]["k_hat"]) for name in names]
        pair_details.append(
            {
                "label": label,
                "pair": names,
                "F_at_0.05": f_values,
                "F_magnitude_ratio": abs(f_values[0] / f_values[1]),
                "k_hat": k_values,
                "k_difference": abs(k_values[0] - k_values[1]),
            }
        )
    v2_defects = {
        "registered_reading_is_F": "F(D) =" in field,
        "builder_declares_TOL_SWAP_for_onset": "TOL_SWAP = 0.25 (C2 onset agreement)" in a3_text,
        "builder_excludes_F_magnitude_from_verdict": "F-magnitude ratio (data, not a verdict" in a3_text,
        "no_F_swap_tolerance_in_preregistered_rule": "TOL_SWAP" not in field,
        "commissioned_far_class_is_singleton": bool(a3["singleton_far"]),
        "substitute_far_side_pair_has_w_conn_4": all(
            int(a2["ftab"][name]["w_conn"]) == 4 for name in a2["pair24"]
        ),
        "commissioned_far_has_w_conn_5": int(a2["ftab"][a2["p_far"]]["w_conn"]) == 5,
    }
    v2_unscoreable = quote_presence["V2"] and all(v2_defects.values())

    near = a2["ftab"][a2["p_near"]]
    far = a2["ftab"][a2["p_far"]]
    v3_true = (
        quote_presence["V3"]
        and abs(float(near["k_hat"]) - int(near["w_conn"])) <= BRACKET_HALF_WIDTH
        and abs(float(far["k_hat"]) - int(far["w_conn"])) <= BRACKET_HALF_WIDTH
        and prior["verdicts"]["B_connected_wenc"] == "NOT_REFUTED"
    )

    tol_back_match = re.search(r"TOL_BACK\s*=\s*([0-9.eE+-]+)", lib_source)
    tol_back = float(tol_back_match.group(1)) if tol_back_match else math.nan
    ba_values = [
        float(row["BA"])
        for placement in a3["c5"].values()
        for row in placement["rows"]
    ]
    max_ba = max(ba_values)
    v4_numeric = math.isfinite(tol_back) and max_ba < tol_back
    v4_provenance = {
        "numeric_comparison_true": v4_numeric,
        "tolerance_absent_from_preregistered_field_document": "TOL_BACK" not in field,
        "same_result_magnitude_survey_disclosed": "scratchpad magnitude survey run BEFORE any sealed execution" in audit,
        "threshold_keyed_above_largest_far_reading": "one order above the largest commissioned far-side reading" in audit,
    }
    v4_qualified_only = quote_presence["V4"] and all(v4_provenance.values())

    v5_true = quote_presence["V5"] and c_result["verdict"] == "NOT_REFUTED"

    items = {
        "V1": {
            "quote": prereg_quotes["V1"],
            "status": "NUMERICALLY_SUPPORTED" if v1_true else "REFUTED",
            "min_far_abs_F": min_far_abs_f,
            "recorded_control_floor": control_floor,
        },
        "V2": {
            "quote": prereg_quotes["V2"],
            "status": "UNSCOREABLE" if v2_unscoreable else "REFUTED",
            "defects": v2_defects,
            "pair_details": pair_details,
            "finding": (
                "the registered reading is F, but the builder scores onset k; the w_conn=5 "
                "commissioned far class has no swap and is replaced by a w_conn=4 pair; no "
                "pre-registered F-equivalence tolerance exists"
            ),
        },
        "V3": {
            "quote": prereg_quotes["V3"],
            "status": "NUMERICALLY_SUPPORTED" if v3_true else "REFUTED",
        },
        "V4": {
            "quote": prereg_quotes["V4"],
            "status": "QUALIFIED_DESCRIPTIVE_ONLY" if v4_qualified_only else "REFUTED",
            "declared_TOL_BACK": tol_back,
            "max_recorded_BA": max_ba,
            "provenance": v4_provenance,
            "finding": (
                "the recorded comparison is below the threshold, but the threshold was selected "
                "after a magnitude survey of the same result family and keyed above the largest "
                "far reading; it is not an adversarial pre-registered boolean"
            ),
        },
        "V5": {
            "quote": prereg_quotes["V5"],
            "status": "NUMERICALLY_SUPPORTED" if v5_true else "REFUTED",
        },
    }
    rule_disposition = (
        "UNSCOREABLE"
        if items["V2"]["status"] == "UNSCOREABLE"
        or items["V4"]["status"] == "QUALIFIED_DESCRIPTIVE_ONLY"
        else "SCOREABLE"
    )
    return {
        "quote_presence": quote_presence,
        "items": items,
        "rule_disposition": rule_disposition,
        "old_ALL_PASS_is_admissible": rule_disposition == "SCOREABLE"
        and all(item["status"] == "NUMERICALLY_SUPPORTED" for item in items.values()),
    }


def render_report(result: dict[str, Any]) -> str:
    integrity = result["input_integrity"]
    components = result["component_verdicts"]
    c = result["C"]
    e = result["E"]
    rule = result["rule_admissibility"]
    lines = [
        "=" * 96,
        "LANE_T51_REPAIR / QUARANTINE_VERIFY -- NO-NEW-MEASUREMENT DEFAULT-REFUTED AUDIT",
        "date: 2026-08-21",
        "=" * 96,
        "",
        "EXACT DISPOSITION",
        f"  OVERALL_VERDICT: {result['overall_verdict']}",
        f"  REGISTERED_RULE_DISPOSITION: {result['registered_rule_disposition']}",
        f"  SCIENTIFIC_DISPOSITION: {result['scientific_disposition']}",
        f"  BRANCH: {result['branch']}",
        "  V2_SCORE: NONE (UNSCOREABLE)",
        "  V4_SCORE: NONE (QUALIFIED DESCRIPTION ONLY)",
        "  ALL_PASS_AUTHORIZED: FALSE",
        "  REGISTRATION_AUTHORIZED: FALSE",
        "  NOTHING IN THIS DIRECTORY IS REGISTRATION EVIDENCE UNTIL A LATER JUDGE ACTS.",
        "",
        "INPUT CUSTODY",
        f"  frozen inputs hashed: {integrity['count']}; mismatches: {len(integrity['mismatches'])}",
        f"  builder manifest cross-check: {integrity['builder_manifest_matches_frozen_map']}",
        f"  prior verifier seal cross-check: {integrity['prior_verifier_manifest_matches_frozen_map']}",
        f"  quarantine-copy manifest cross-check: {integrity['quarantine_manifest_matches_frozen_map']}",
        f"  integrity verdict: {integrity['verdict']}",
        "  No builder, sealed program, Hamiltonian, or measurement script was executed.",
        "",
        "OLD DEFECTS PRESERVED",
        f"  old C predicate defect confirmed: {result['old_defects']['C_defect']['confirmed']}",
        "    old C used only lambda=0.037 and near/far, despite recording 0.071 and fresh.",
        f"  old E predicate defect confirmed: {result['old_defects']['E_defect']['confirmed']}",
        "    old E walked only top-level files and returned literal violation_count=0.",
        "",
        "A/B/D RECORDED-EVIDENCE RE-ADJUDICATION",
        f"  A sector exactness: {components['A']}",
        f"  B connected minimum: {components['B']}",
        f"  D off-grid onset fit: {components['D']}",
        "",
        "C -- STRENGTHENED WINDING ATTRIBUTION",
        f"  verdict: {c['verdict']}",
        "  required cross-product: placements near/far/fresh x lambdas 0.037/0.071.",
    ]
    for row in c["checks"]:
        fits = row.get("recomputed", {}).get("fits", [])
        signs = row.get("recomputed", {}).get("signs", {})
        lines.append(
            f"    {row['placement']} lambda={row['lambda']:.3f}: ok={row['ok']} "
            f"fits={fits} signs={signs}"
        )
    lines += [
        "",
        "E -- CLOSED-CORPUS DIRECTIVE SCAN",
        f"  corpus: {e['corpus_count']} frozen hash-enumerated claim-bearing files",
        f"  positive controls all matched: {e['positive_controls_all_matched']}",
        f"  candidate hits: {len(e['candidate_hits'])}; decision hits: {e['decision_hit_count']}",
        f"  verdict derived from controls and decision hits: {e['verdict']}",
        f"  explicit self-reference exclusion: {e['self_reference_exclusions']['owned_verifier_tree']}",
        "",
        "PRE-REGISTERED V1--V5 ADMISSIBILITY (NO JUDGE, NO NEW MEASUREMENT)",
    ]
    for name, item in rule["items"].items():
        lines.append(f"  {name}: {item['status']} -- {item['quote']}")
    lines += [
        "",
        "LOAD-BEARING FINDINGS",
        "  V2 is UNSCOREABLE: FIELD defines the reading as F(D), while C2's boolean uses onset",
        "  k; the commissioned w_conn=5 far class is a singleton and the substitute swap is",
        "  w_conn=4; no pre-registered F-equivalence tolerance exists.",
        "  V4 is QUALIFIED_DESCRIPTIVE_ONLY: max recorded BA is below 1e-3, but the audit",
        "  discloses that 1e-3 was selected after a magnitude survey and set one order above",
        "  the largest commissioned far reading.  That does not support an adversarial PASS.",
        "  Therefore the old ALL PASS assertion and the prior judgment's 'no UNSCOREABLE item'",
        "  premise are REFUTED.  C and E passing does not cure the scoring defects.",
        "  A-E NOT_REFUTED is a verifier-component result only; it is not V1-V5 ALL PASS.",
        "",
        "FINAL: PARTIAL POSITIVE / REPAIR REQUIRED.  REGISTER NOTHING.",
        "=" * 96,
    ]
    return "\n".join(lines) + "\n"


def render_audit(result: dict[str, Any]) -> str:
    c_lines = result["old_defects"]["C_defect"]["source_evidence"]
    e_lines = result["old_defects"]["E_defect"]["source_evidence"]
    return f"""LANE_T51_REPAIR / QUARANTINE_VERIFY / D24_AUDIT -- 2026-08-21

SCOPE AND CUSTODY.
This is a no-new-measurement verifier of the existing T51 package.  It read and SHA-256 checked
{result['input_integrity']['count']} frozen inputs.  It did not import, execute, or modify any builder,
sealed program, model, ledger, register, proof, or shared document.  It wrote only in its owned
LANE_T51_REPAIR/QUARANTINE_VERIFY directory.  Coordinates in the frozen inputs remain construction
labels; this audit makes no new separation claim.

DEFAULT-REFUTED DISCIPLINE.
Every component begins REFUTED.  A/B/D were re-adjudicated only from the prior independent
verifier's frozen raw records.  C was recomputed from every stored F value at the full required
cross-product near/far/fresh x lambda 0.037/0.071.  E was derived from actual scanner hits; no
violation count or clean verdict is a literal override.  The overall result remains REFUTED because
the old ALL-PASS package depends on an UNSCOREABLE V2 and a descriptive-only V4.

OLD C DEFECT, PRESERVED IN PLACE.
Confirmed: {result['old_defects']['C_defect']['confirmed']}.
Source evidence: {json.dumps(c_lines, sort_keys=True)}
The old code computed both lambdas and fresh placement but its verdict predicate consumed only
CHECK_LAMBDA=0.037 and only near/far.  This verifier requires all six cells, recomputes characters
from raw F signs, requires unique fits, stored/recomputed agreement, expected sign-flip behavior,
and lambda stability.  Result: {result['C']['verdict']}.

OLD E DEFECT, PRESERVED IN PLACE.
Confirmed: {result['old_defects']['E_defect']['confirmed']}.
Source evidence: {json.dumps(e_lines, sort_keys=True)}
The old top-level iterdir excluded nested verifier/judgment artifacts and returned a literal
violation_count=0 after a prose assertion.  This verifier uses a frozen {result['E']['corpus_count']}-file
hash-enumerated corpus, including nested prior verifier outputs and JUDGMENT_CODEX.  The owned
verifier tree is explicitly excluded because it contains the regex rules, refutation prose, and
deliberate positive-control violations.  Integrity manifests/sidecars are excluded as non-claim
metadata.  Each of the three rule classes has an in-memory positive-control seed.  All controls
matched: {result['E']['positive_controls_all_matched']}; decision hits: {result['E']['decision_hit_count']};
E verdict: {result['E']['verdict']}.

SCORING-ADMISSIBILITY FINDINGS.
V2: UNSCOREABLE.  The registered observable is F(D), but the builder's C2 boolean compares onset
orders; its own output labels F-magnitude ratios data and not a verdict.  The commissioned
w_enc_conn=5 far class is a singleton, and the replacement pair is w_enc_conn=4.  No F-equivalence
tolerance was pre-registered, so no honest boolean can be recovered after seeing the values.

V4: QUALIFIED_DESCRIPTIVE_ONLY.  The stored BA values are below TOL_BACK=1e-3, but the builder audit
states that this threshold was chosen after a same-family scratchpad magnitude survey and set one
order above the largest commissioned far-side reading.  The numerical comparison may be reported;
it does not earn a pre-registered adversarial PASS.

ERROR AND CORRECTION LOG.
No scientific execution occurred.  No prior artifact was superseded.  The first quarantine scan
returned one E decision hit in a4_verdicts.py because its forward-only three-line window began on
the wrapped continuation "outcome here is a failure" and omitted the preceding word "no".  The
window was corrected to carry two lines of preceding and following context; the positive controls
still match, while the sentence is now deterministically exculpated as an explicit no-failure
statement.  A second dry adjudication exposed two input-parser defects: Markdown backticks/bold
prevented exact V1/V3 quote matching, and TOL_BACK was sought in a3 rather than its imported
t51a_lib definition.  Claim text is now Markdown-normalized and TOL_BACK is parsed from the frozen
library source.  No scientific number or package artifact changed.  This standing run completed
with input mismatches={len(result['input_integrity']['mismatches'])}, C={result['C']['verdict']},
E={result['E']['verdict']}.  Any implementation error found before this file's final generation is
recorded by updating this paragraph before the directory is handed back; no other error is known.

CUSTODY CORRECTION BEFORE SEALING.
A pre-landing replay audit found that the first version of this still-unsealed verifier referenced
the prior drafts only at their untracked, final-looking LANE_T51_A paths.  Exact byte copies were
preserved under LANE_T51_REPAIR/QUARANTINED_DRAFTS/ORIGINAL_LAYOUT, an explicit non-admissibility
notice and six-file INPUTS.sha256 were added, and every verifier reference was repointed to those
copies.  The original untracked files were neither modified nor removed.  The copied historical
SEALS file proves byte integrity only; it does not validate the old predicates.  This correction
changes custody paths and input count only, not any scientific value or adjudication.

EXACT DISPOSITION.
OVERALL_VERDICT={result['overall_verdict']}
REGISTERED_RULE_DISPOSITION={result['registered_rule_disposition']}
SCIENTIFIC_DISPOSITION={result['scientific_disposition']}
BRANCH={result['branch']}
V2_SCORE=NONE_UNSCOREABLE
V4_SCORE=NONE_QUALIFIED_DESCRIPTION_ONLY
ALL_PASS_AUTHORIZED=FALSE
REGISTRATION_AUTHORIZED=FALSE
REGISTER NOTHING.
"""


def run(repo: Path, output_dir: Path) -> dict[str, Any]:
    if output_dir.resolve() != (repo / OWNED_PREFIX.rstrip("/")).resolve():
        raise SystemExit("refusing to write outside the owned quarantine verifier directory")

    integrity = verify_frozen_inputs(repo)
    old_defects = preserve_old_defects(repo)
    prior = json.loads(
        (
            repo
            / (QUARANTINED_PREFIX + "VERIFY_CODEX/v1_independent_rebuild.json")
        ).read_text(
            encoding="utf-8"
        )
    )
    abd = verify_a_b_d(prior)
    c_result = verify_c(prior)
    e_result = scan_claim_corpus(repo)
    rule = score_admissibility(repo, prior, c_result)

    component_verdicts = {
        "A": abd["A"]["verdict"],
        "B": abd["B"]["verdict"],
        "C": c_result["verdict"],
        "D": abd["D"]["verdict"],
        "E": e_result["verdict"],
    }
    core_not_refuted = (
        integrity["verdict"] == "NOT_REFUTED"
        and all(value == "NOT_REFUTED" for value in component_verdicts.values())
        and old_defects["C_defect"]["confirmed"]
        and old_defects["E_defect"]["confirmed"]
    )
    partial_positive = (
        core_not_refuted
        and rule["items"]["V1"]["status"] == "NUMERICALLY_SUPPORTED"
        and rule["items"]["V3"]["status"] == "NUMERICALLY_SUPPORTED"
        and rule["items"]["V5"]["status"] == "NUMERICALLY_SUPPORTED"
    )
    overall = "REFUTED"
    registered_rule = "UNSCOREABLE" if rule["rule_disposition"] == "UNSCOREABLE" else "REFUTED"
    scientific = "PARTIAL_POSITIVE" if partial_positive else "REFUTED"
    branch = "REPAIR_REQUIRED"

    result = {
        "metadata": {
            "date": "2026-08-21",
            "mode": "no-new-measurement",
            "default": "REFUTED",
            "owned_output_directory": OWNED_PREFIX.rstrip("/"),
        },
        "input_integrity": integrity,
        "old_defects": old_defects,
        "A": abd["A"],
        "B": abd["B"],
        "C": c_result,
        "D": abd["D"],
        "E": e_result,
        "component_verdicts": component_verdicts,
        "rule_admissibility": rule,
        "overall_verdict": overall,
        "registered_rule_disposition": registered_rule,
        "scientific_disposition": scientific,
        "branch": branch,
        "V2_scored": False,
        "V4_scored": False,
        "all_pass_authorized": False,
        "registration_authorized": False,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "quarantine_verify.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "quarantine_verify.txt").write_text(
        render_report(result), encoding="utf-8"
    )
    (output_dir / "D24_AUDIT.txt").write_text(
        render_audit(result), encoding="utf-8"
    )
    return result


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: quarantine_verify.py REPO OWNED_OUTPUT_DIR")
    result = run(Path(sys.argv[1]), Path(sys.argv[2]))
    print(
        json.dumps(
            {
                "overall_verdict": result["overall_verdict"],
                "registered_rule_disposition": result["registered_rule_disposition"],
                "scientific_disposition": result["scientific_disposition"],
                "branch": result["branch"],
                "component_verdicts": result["component_verdicts"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
