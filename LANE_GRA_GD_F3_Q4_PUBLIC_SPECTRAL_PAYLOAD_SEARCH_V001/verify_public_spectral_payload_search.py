#!/usr/bin/env python3
"""Static replay for the bounded GD public spectral-payload search."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

LANE = Path(__file__).resolve().parent
ROOT = LANE.parent
checks: list[str] = []


def check(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)
    checks.append(label)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


required = {
    "README.md",
    "RESULT.md",
    "RESULT.json",
    "PAYLOAD_BINDING_MATRIX.md",
    "SEARCH_LEDGER.md",
    "SOURCE_CUSTODY.json",
    "verify_public_spectral_payload_search.py",
}
for name in sorted(required):
    check((LANE / name).is_file(), f"required file: {name}")

custody = json.loads((LANE / "SOURCE_CUSTODY.json").read_text())
result = json.loads((LANE / "RESULT.json").read_text())

check(custody["lane"] == result["lane"], "lane identifiers agree")
check(custody["retrieved_utc_date"] == "2026-08-28", "retrieval date pinned")

for dep in custody["local_dependencies"]:
    path = (LANE / dep["path"]).resolve()
    check(path.is_file(), f"dependency exists: {path.name}")
    check(sha256(path) == dep["sha256"], f"dependency hash: {path.name}")

roots = {item["id"]: item for item in custody["primary_roots"]}
check(len(roots) == 5, "five primary roots inventoried")
check(roots["SHANNON-2012-ARXIV-SOURCE"]["raw_numeric_arrays"] is False,
      "Shannon source has no raw arrays")
check(roots["SIKORA-2011-ARXIV-SOURCE"]["program_source"] is False,
      "Sikora source has no program")
check(roots["HUANG-2018-ARXIV-SOURCE"]["raw_numeric_arrays"] is False,
      "Huang source has no raw arrays")
check(roots["HKU-QFI-2026-DATA-METADATA"]["doi"] ==
      "10.25442/hku.32404548.v1", "HKU data DOI pinned")
check(roots["HKU-QFI-2026-CODE-METADATA"]["doi"] ==
      "10.25442/hku.32412273.v1", "HKU code DOI pinned")

check(len(custody["hku_data_files"]) == 10, "ten HKU data files hashed")
check(len(custody["hku_code_files"]) == 4, "four HKU code files hashed")
check(len({name for name, _ in custody["hku_data_files"]}) == 10,
      "HKU data filenames unique")
check(len({digest for _, digest in custody["hku_code_files"]}) == 4,
      "HKU code hashes unique")
check(all(len(digest) == 64 for _, digest in custody["hku_data_files"]),
      "HKU data SHA-256 shape")
check(all(len(digest) == 64 for _, digest in custody["hku_code_files"]),
      "HKU code SHA-256 shape")

notes = custody["inspection_notes"]
check(notes["qed_native_max_body"] == 3, "QED body-rank ceiling pinned")
check("runtime rejection" in notes["qed_ring_exchange"],
      "QED ring rejection pinned")
check("full XXZ" in notes["qmc_released_operator"],
      "QMC Hamiltonian/operator mismatch pinned")
check("not released" in notes["qmc_covariance"],
      "covariance gap pinned")

check(len(custody["repository_searches"]) == 3,
      "three repository-native searches pinned")
check({entry["provider"] for entry in custody["repository_searches"]} ==
      {"Figshare", "Zenodo"}, "repository providers pinned")

check(result["qualifying_public_packet_found"] is False,
      "no qualifying packet result")
check(result["new_executable_root_found"] is True,
      "new executable root result")
check(result["new_root"]["direct_gc_fy_admissibility"] is False,
      "new root not overpromoted")
check(len(result["missing_minimal_payload"]) == 9,
      "nine exact payload gaps")
check("adapt" in result["strongest_next_step"],
      "next step is source adaptation")
check("non-exhaustive" in result["ceiling"],
      "non-exhaustive ceiling explicit")

docs = "\n".join((LANE / name).read_text() for name in
                 ("README.md", "RESULT.md", "PAYLOAD_BINDING_MATRIX.md",
                  "SEARCH_LEDGER.md"))
for phrase in ("no such data exist", "proves absence", "gravity is proved",
               "common cone is proved"):
    check(phrase not in docs.lower(), f"forbidden overclaim absent: {phrase}")
check("GC16--GC19" in docs, "exact target named")
check("implementation substrate" in docs, "substrate/evidence distinction named")

print(f"PASS {len(checks)}/{len(checks)}")

