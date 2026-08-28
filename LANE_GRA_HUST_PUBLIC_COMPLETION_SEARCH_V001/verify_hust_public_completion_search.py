#!/usr/bin/env python3
"""Verify the static custody and claim discipline of the HUST completion search."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name: str):
    return json.loads((HERE / name).read_text(encoding="utf-8"))


checks: list[tuple[str, bool]] = []


def check(label: str, condition: bool) -> None:
    checks.append((label, bool(condition)))


ledger = load("SEARCH_LEDGER.json")
result = load("RESULT.json")
custody = load("SOURCE_CUSTODY.json")
theorem = (HERE / "THEOREM.md").read_text(encoding="utf-8")

check("ledger schema", ledger["schema"] == "WAC_HUST_PUBLIC_COMPLETION_SEARCH_LEDGER_V001")
check("result schema", result["schema"] == "WAC_HUST_PUBLIC_COMPLETION_SEARCH_RESULT_V001")
check("custody schema", custody["schema"] == "WAC_HUST_PUBLIC_COMPLETION_SEARCH_CUSTODY_V001")
check("single query date", {ledger["date"], result["date"], custody["date"]} == {"2026-08-27"})
check("target DOI", ledger["target_doi"] == "10.1038/s41586-018-0431-5")

rule = ledger["qualification_rule"]
check("four-part admission rule", len(rule) == 4 and all(rule.values()))
surfaces = ledger["surfaces"]
check("declared surface count", len(surfaces) == ledger["counts"]["declared_surface_entries"] == 11)
check("unique surface IDs", len({s["id"] for s in surfaces}) == len(surfaces))
scope = ledger["search_evidence_scope"]
reproducible_ids = set(scope["reproducible_literal_endpoint_entries"])
curator_ids = set(scope["curator_recorded_nonexecutably_frozen_entries"])
check("search evidence classes partition all entries",
      not (reproducible_ids & curator_ids)
      and reproducible_ids | curator_ids == {s["id"] for s in surfaces})
check("aggregate entries are curator-recorded not executable",
      all(s.get("reproducibility") == "CURATOR_RECORDED_NONEXECUTABLE"
          for s in surfaces if s["id"] in curator_ids))
check("zero qualifying roots in every surface", all(s["qualifying_completion_roots"] == 0 for s in surfaces))
check("zero qualifying roots total", ledger["counts"]["qualifying_completion_roots"] == 0)
check("two acquisition leads", ledger["counts"]["new_targeted_acquisition_leads"] == 2)
check("publisher inventory searched", any(s["id"] == "nature_official_release_inventory" for s in surfaces))
check("DOI registries searched", {"crossref_doi_relation", "datacite_exact_related_identifier"}.issubset({s["id"] for s in surfaces}))
check("data repositories searched", {"figshare_exact_resource_doi", "zenodo_exact_doi", "dryad_osf_targeted"}.issubset({s["id"] for s in surfaces}))
check("institutional surfaces searched", {"hust_institutional_pages", "campaign_dissertation_catalogues"}.issubset({s["id"] for s in surfaces}))

r = result["result"]
check("no root claimed", r["qualifying_public_completion_root_located"] is False)
check("no numerical G advance claimed", r["numerical_G_crosscheck_advanced"] is False)
check("parent family unchanged", r["prior_calibrated_forward_family_changed"] is False)
check("world exhaustiveness denied", r["search_surface_exhaustive_of_world"] is False)
check("seven exact missing-object groups", len(result["missing_objects"]) == 7)
check("all missing statuses bounded", all(x["status"] == "NOT_LOCATED_IN_BOUNDED_PUBLIC_SEARCH" for x in result["missing_objects"]))
required_ids = {
    "row_harmonic_remainders",
    "measured_planar_density_map",
    "source_sphere_density_orientation_maps",
    "cmm_residual_coordinates",
    "attachment_coating_maps",
    "aaf_shelf_deformation_compensation_maps",
    "campaign_raw_correction_covariance_packet",
}
check("missing-object IDs complete", {x["id"] for x in result["missing_objects"]} == required_ids)
check("two leads not treated as payload", len(result["targeted_acquisition_leads"]) == 2 and all(not x["payload_inspected"] for x in result["targeted_acquisition_leads"]))
check("all leads acquisition-only", all(x["claim"] == "ACQUISITION_TARGET_ONLY" for x in result["targeted_acquisition_leads"]))
leads = result["targeted_acquisition_leads"]
check("one confirmed dissertation lead",
      leads[0]["record_type"] == "CONFIRMED_HUST_2021_DOCTORAL_DISSERTATION_LEAD")
check("second lead is unverified title only",
      leads[1]["record_type"] == "UNVERIFIED_TITLE_ONLY_BIBLIOGRAPHIC_LEAD"
      and leads[1]["author"] is None and leads[1]["institution"] is None
      and leads[1]["year"] is None)

check("no new binary admitted", custody["new_binary_sources_acquired"] == 0)
check("no completion root admitted", custody["new_qualifying_completion_roots"] == 0)
dependency_paths = {dep["path"] for dep in custody["hash_pinned_dependencies"]}
check("seven-object inventory dependency is explicitly pinned",
      "../LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/SOURCE_CUSTODY.json"
      in dependency_paths)
for dep in custody["hash_pinned_dependencies"]:
    path = (HERE / dep["path"]).resolve()
    check(f"dependency exists: {dep['path']}", path.is_file())
    check(f"dependency hash: {dep['path']}", path.is_file() and digest(path) == dep["sha256"])

lower = theorem.lower()
check("theorem says bounded", "bounded" in lower)
check("theorem disclaims global absence", "makes no claim that the requested data do not exist in the world" in lower)
check("theorem preserves G ceiling", "adds no numerical value of" in lower)
check("theorem keeps author request path", "corresponding authors" in lower)
check("theorem distinguishes curator search from reproducible endpoints",
      "curator-recorded dated searches" in lower
      and "not executable completeness certificates" in lower)
check("theorem distinguishes prior inventory from local admission",
      "does not mean every remote object was locally downloaded and admitted" in lower)
check("ledger has non-exhaustive ceiling", "does not establish" in ledger["non_exhaustive_ceiling"].lower())

failed = [label for label, ok in checks if not ok]
for label, ok in checks:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
print(f"SUMMARY {len(checks) - len(failed)}/{len(checks)} checks passed")
if failed:
    raise SystemExit(1)
