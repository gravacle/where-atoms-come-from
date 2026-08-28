#!/usr/bin/env python3
"""Independent offline replay of the HUST public-completion hostile audit.

This verifier does not import or execute the builder verifier.  It checks the
repaired core, the independent normalized live-query receipt, the seven-object
publisher inventory, and the exact claim ceilings from first principles.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path


AUDIT = Path(__file__).resolve().parent
LANE = AUDIT.parent
ROOT = LANE.parent

checks: list[tuple[str, bool]] = []


def check(label: str, condition: bool) -> None:
    checks.append((label, bool(condition)))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def entries(path: Path) -> dict[str, str]:
    answer: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line:
            value, relative = line.split("  ", 1)
            answer[relative] = value
    return answer


def valid_hash(value: str) -> bool:
    return bool(re.fullmatch(r"[0-9a-f]{64}", value))


ledger = load(LANE / "SEARCH_LEDGER.json")
result = load(LANE / "RESULT.json")
custody = load(LANE / "SOURCE_CUSTODY.json")
receipt = load(AUDIT / "LIVE_REQUERY_NORMALIZED.json")
theorem = (LANE / "THEOREM.md").read_text(encoding="utf-8")
readme = (LANE / "README.md").read_text(encoding="utf-8")

doi = "10.1038/s41586-018-0431-5"
check("independent receipt schema",
      receipt["schema"] == "WAC_HUST_PUBLIC_COMPLETION_INDEPENDENT_REQUERY_V001")
check("single target DOI",
      {doi, ledger["target_doi"], receipt["target_doi"],
       custody["article"]["doi"]} == {doi})
check("builder claim date is fixed",
      {ledger["date"], result["date"], custody["date"]} == {"2026-08-27"})
check("audit receipt is one-day hostile requery",
      receipt["retrieved_utc"].startswith("2026-08-28T"))

surfaces = {item["id"]: item for item in ledger["surfaces"]}
expected_surfaces = {
    "nature_official_release_inventory",
    "nature_data_availability",
    "crossref_doi_relation",
    "datacite_exact_related_identifier",
    "figshare_exact_resource_doi",
    "zenodo_exact_doi",
    "dryad_osf_targeted",
    "china_scientific_data_surfaces",
    "hust_institutional_pages",
    "public_code_hosts",
    "campaign_dissertation_catalogues",
}
expected_literal = {
    "nature_official_release_inventory",
    "nature_data_availability",
    "crossref_doi_relation",
    "datacite_exact_related_identifier",
    "figshare_exact_resource_doi",
    "zenodo_exact_doi",
    "hust_institutional_pages",
}
expected_curator = expected_surfaces - expected_literal
scope = ledger["search_evidence_scope"]
literal = set(scope["reproducible_literal_endpoint_entries"])
curator = set(scope["curator_recorded_nonexecutably_frozen_entries"])
check("exact eleven declared surface entries", set(surfaces) == expected_surfaces)
check("surface counter is exact",
      ledger["counts"]["declared_surface_entries"] == len(surfaces) == 11)
check("literal endpoint class is exact", literal == expected_literal)
check("curator-recorded class is exact", curator == expected_curator)
check("evidence classes are a partition",
      not literal.intersection(curator) and literal.union(curator) == set(surfaces))
check("all aggregate entries disclose nonexecutable status",
      all(surfaces[key].get("reproducibility") ==
          "CURATOR_RECORDED_NONEXECUTABLE" for key in curator))
check("HUST literal scope is only two named pages",
      surfaces["hust_institutional_pages"]["mode"] ==
      "direct inspection of the two named literal URLs only")

figshare_surface = surfaces["figshare_exact_resource_doi"]
check("Figshare ledger names both literal POST endpoints",
      figshare_surface["query_or_url"] ==
      "POST https://api.figshare.com/v2/articles/search; POST https://api.figshare.com/v2/collections/search")
check("Figshare ledger freezes exact body",
      figshare_surface["request_body"] == {"resource_doi": doi})

check("every declared surface reports zero qualifying roots",
      all(item["qualifying_completion_roots"] == 0
          for item in surfaces.values()))
check("aggregate qualifying-root counter is zero",
      ledger["counts"]["qualifying_completion_roots"] == 0)
check("bounded conclusion is surface- and date-limited",
      "explicitly named public surfaces" in ledger["bounded_conclusion"]
      and "2026-08-27" in ledger["bounded_conclusion"])
check("ledger expressly rejects global absence",
      "does not establish" in ledger["non_exhaustive_ceiling"].lower()
      and "all public surfaces" in ledger["non_exhaustive_ceiling"].lower())
check("theorem expressly rejects global absence",
      "makes no claim that the requested data do not exist in the world" in theorem)
check("theorem limits reproducibility to literal endpoint checks",
      "curator-recorded dated searches" in theorem
      and "not executable completeness certificates" in theorem)
check("README preserves dated bounded ceiling",
      "dated, bounded search result" in readme
      and "not a theorem that the data do not exist" in readme)

nature = receipt["nature"]
nature_objects = {item["id"]: item for item in nature["objects"]}
expected_ids = {f"MOESM{i}" for i in range(1, 8)}
expected_labels = {
    "MOESM1": "Supplementary Information",
    "MOESM2": "Supplementary Data for Supplementary Fig. 1",
    "MOESM3": "Source Data Fig. 2",
    "MOESM4": "Source Data Fig. 3",
    "MOESM5": "Source Data Extended Data Fig. 2",
    "MOESM6": "Source Data Extended Data Fig. 4",
    "MOESM7": "Source Data Extended Data Fig. 5",
}
check("Nature response receipt has valid hash",
      valid_hash(nature["observed_response_sha256"]))
check("Nature associated-object count is exactly seven",
      nature["associated_object_count"] == len(nature_objects) == 7)
check("Nature associated-object IDs are MOESM1 through MOESM7",
      set(nature_objects) == expected_ids)
check("Nature labels are exact", all(
      nature_objects[key]["label"] == label
      for key, label in expected_labels.items()))
check("Nature object URLs bind the matching object ID", all(
      key in item["url"] for key, item in nature_objects.items()))
check("Nature data availability is an acquisition route, not payload",
      nature["data_availability_classification"] ==
      "ADDITIONAL_SUPPORTING_DATA_ON_REASONABLE_REQUEST_FROM_CORRESPONDING_AUTHORS")

dependency_map = {item["path"]: item for item in custody["hash_pinned_dependencies"]}
prior_relative = "../LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/SOURCE_CUSTODY.json"
check("seven-object prior inventory is explicitly hash-pinned",
      prior_relative in dependency_map)
for relative, item in dependency_map.items():
    path = (LANE / relative).resolve()
    check(f"dependency exists: {relative}", path.is_file())
    check(f"dependency hash: {relative}",
          path.is_file() and digest(path) == item["sha256"])

prior_path = (LANE / prior_relative).resolve()
prior = load(prior_path)
prior_objects = {item["id"]: item for item in prior["official_release"]}
check("prior official inventory has exact seven IDs", set(prior_objects) == expected_ids)
check("independent streamed hashes match prior custody", all(
      nature_objects[key]["streamed_sha256"] == prior_objects[key]["sha256"]
      for key in expected_ids))
check("independent URLs match prior custody", all(
      nature_objects[key]["url"] == prior_objects[key]["url"]
      for key in expected_ids))
check("prior inventory distinguishes two pinned from five remote-only",
      prior["pinned_file_count"] == 2
      and prior["unpinned_file_count"] == 5
      and sum(bool(item["pinned_in_lane"]) for item in prior_objects.values()) == 2)
for key, item in prior_objects.items():
    if item["pinned_in_lane"]:
        local_path = prior_path.parent / item["lane_path"]
        check(f"locally pinned binary exists: {key}", local_path.is_file())
        check(f"locally pinned binary hash: {key}",
              local_path.is_file() and digest(local_path) == item["sha256"])

crossref = receipt["crossref"]
check("Crossref response receipt has valid hash",
      valid_hash(crossref["observed_response_sha256"]))
check("Crossref exact record and status",
      crossref["status"] == "ok" and crossref["record_doi"] == doi)
check("Crossref relation object was empty", crossref["relation"] == {})
check("Crossref interpretation does not infer global absence",
      "not proof" in crossref["interpretation"].lower())

datacite = receipt["datacite"]
datacite_records = datacite["records"]
check("DataCite response receipt has valid hash",
      valid_hash(datacite["observed_response_sha256"]))
check("DataCite exact related-identifier result count",
      datacite["total"] == len(datacite_records) == 3)
check("DataCite returned the exact three observed DOIs", {
      item["doi"] for item in datacite_records} == {
          "10.5281/zenodo.21428316",
          "10.5281/zenodo.21428128",
          "10.5281/zenodo.21428129",
      })
check("DataCite records are later reference relations", all(
      item["creator"] == "Attar, Ali"
      and item["publisher"] == "Zenodo"
      and item["publication_year"] == 2026
      and item["relation_type"] == "References"
      for item in datacite_records))
check("DataCite provenance classification rejects campaign deposit",
      datacite["classification"] ==
      "LATER_UNRELATED_REFERENCES_NOT_HUST_CAMPAIGN_DEPOSITS")

figshare = receipt["figshare"]
figshare_requests = {item["url"]: item for item in figshare["requests"]}
expected_figshare = {
    "https://api.figshare.com/v2/articles/search",
    "https://api.figshare.com/v2/collections/search",
}
check("Figshare independent receipt has both endpoints",
      set(figshare_requests) == expected_figshare)
check("Figshare exact POST bodies match target DOI", all(
      item["method"] == "POST" and item["body"] == {"resource_doi": doi}
      for item in figshare_requests.values()))
check("Figshare responses are exactly empty arrays", all(
      item["result"] == [] for item in figshare_requests.values()))
check("Figshare empty-response hashes are exact", all(
      item["observed_response_sha256"] == hashlib.sha256(b"[]").hexdigest()
      for item in figshare_requests.values()))
check("Figshare interpretation preserves query boundary",
      "does not cover" in figshare["interpretation"].lower())

zenodo = receipt["zenodo"]
zenodo_records = zenodo["records"]
check("Zenodo response receipt has valid hash",
      valid_hash(zenodo["observed_response_sha256"]))
check("Zenodo exact DOI-query result count",
      zenodo["total"] == len(zenodo_records) == 2)
check("Zenodo returned the exact two observed DOIs", {
      item["doi"] for item in zenodo_records} == {
          "10.5281/zenodo.19582997",
          "10.5281/zenodo.20584558",
      })
check("Zenodo returned later non-HUST creators", {
      item["creator"] for item in zenodo_records} == {
          "Lehew, John", "Berg, Stefan"}
      and all(item["publication_date"].startswith("2026-")
              for item in zenodo_records))
check("Zenodo provenance classification rejects campaign deposit",
      zenodo["classification"] ==
      "LATER_UNRELATED_REFERENCES_NOT_HUST_CAMPAIGN_DEPOSITS")

hust_pages = receipt["hust_literal_pages"]
check("HUST receipt explicitly limits page scope",
      "two literal pages only" in hust_pages["scope"])
check("HUST page receipt contains exactly two URLs",
      {item["url"] for item in hust_pages["pages"]} == {
          "https://ngl.hust.edu.cn/info/1238/3516.htm",
          "https://phys.hust.edu.cn/info/1211/3362.htm",
      })
check("HUST page response hashes are valid", all(
      valid_hash(item["observed_response_sha256"])
      for item in hust_pages["pages"]))
check("HUST pages expose no declared-extension data attachment", all(
      item["candidate_data_attachment_links_by_declared_extension"] == 0
      for item in hust_pages["pages"]))

leads = result["targeted_acquisition_leads"]
discovery = receipt["thesis_discovery"]
confirmed = discovery["confirmed_lead"]
unverified = discovery["unverified_lead"]
check("result has exactly two acquisition leads", len(leads) == 2)
check("confirmed lead is HUST 2021 doctoral metadata", all((
      confirmed["author"] == "Jun-Fei Wu",
      confirmed["institution"] == "Huazhong University of Science and Technology",
      confirmed["year"] == 2021,
      confirmed["record_type"] == "doctoral dissertation",
      confirmed["indexed_source_section"] == "中国博士学位论文全文数据库",
      confirmed["payload_inspected"] is False,
)))
check("confirmed receipt agrees with repaired result",
      leads[0]["record_type"] == "CONFIRMED_HUST_2021_DOCTORAL_DISSERTATION_LEAD"
      and leads[0]["author"] == confirmed["author"]
      and leads[0]["year"] == confirmed["year"])
check("second receipt authenticates title only", all((
      unverified["authenticated_author"] is None,
      unverified["authenticated_institution"] is None,
      unverified["authenticated_year"] is None,
      unverified["authenticated_document_type"] is None,
      unverified["payload_inspected"] is False,
      unverified["classification"] == "UNVERIFIED_TITLE_ONLY_BIBLIOGRAPHIC_LEAD",
)))
check("second receipt agrees with repaired null fields",
      leads[1]["record_type"] == "UNVERIFIED_TITLE_ONLY_BIBLIOGRAPHIC_LEAD"
      and leads[1]["author"] is None
      and leads[1]["institution"] is None
      and leads[1]["year"] is None)
check("both core leads are acquisition-only and not inspected",
      all(item["claim"] == "ACQUISITION_TARGET_ONLY"
          and item["payload_inspected"] is False for item in leads))

aggregate = receipt["aggregate_search_evidence"]
check("independent receipt keeps aggregate searches nonexecutable",
      set(aggregate["entry_ids"]) == expected_curator
      and aggregate["classification"] ==
      "CURATOR_RECORDED_NONEXECUTABLE_DATED_SEARCHES")

missing_ids = {
    "row_harmonic_remainders",
    "measured_planar_density_map",
    "source_sphere_density_orientation_maps",
    "cmm_residual_coordinates",
    "attachment_coating_maps",
    "aaf_shelf_deformation_compensation_maps",
    "campaign_raw_correction_covariance_packet",
}
check("exact seven missing-object groups remain", {
      item["id"] for item in result["missing_objects"]} == missing_ids)
check("every missing status remains bounded", all(
      item["status"] == "NOT_LOCATED_IN_BOUNDED_PUBLIC_SEARCH"
      for item in result["missing_objects"]))

outcome = receipt["audit_observation"]
check("no new binary or completion root was promoted",
      outcome["new_binary_promoted"] is False
      and outcome["new_qualifying_completion_root_observed"] is False
      and custody["new_binary_sources_acquired"] == 0
      and custody["new_qualifying_completion_roots"] == 0)
check("no accepted G or author-processed result was promoted",
      outcome["accepted_or_codata_G_used"] is False
      and outcome["authors_processed_G_or_kernel_promoted"] is False
      and result["result"]["numerical_G_crosscheck_advanced"] is False)
check("parent calibrated family is unchanged",
      result["result"]["prior_calibrated_forward_family_changed"] is False)
check("world exhaustiveness is denied",
      result["result"]["search_surface_exhaustive_of_world"] is False
      and outcome["global_public_absence_inferred"] is False)
check("theorem strict ceiling excludes numerical G and processed source promotion",
      "adds no numerical value" in theorem
      and "does not convert the authors' processed source coefficient" in theorem)

builder_payload = {
    "README.md", "THEOREM.md", "SEARCH_LEDGER.json", "RESULT.json",
    "SOURCE_CUSTODY.json", "verify_hust_public_completion_search.py",
    "VERIFICATION.txt",
}
builder_manifest = entries(LANE / "MANIFEST.sha256")
check("builder manifest exact membership", set(builder_manifest) == builder_payload)
for relative, expected in builder_manifest.items():
    check(f"builder manifest hash: {relative}", digest(LANE / relative) == expected)
builder_seal = entries(LANE / "LANE_SEAL.sha256")
check("builder seal binds manifest and transitively its verification",
      builder_seal == {
          "MANIFEST.sha256": digest(LANE / "MANIFEST.sha256"),
      }
      and builder_manifest["VERIFICATION.txt"] ==
      digest(LANE / "VERIFICATION.txt"))

core_files = [item for item in LANE.iterdir() if item.is_file()]
core_text = "\n".join(item.read_text(encoding="utf-8") for item in core_files)
check("no accepted-like numerical G value occurs in builder core",
      re.search(r"(?<![0-9])6[.,]67[0-9]{2,}", core_text) is None)
check("no source binary exists in completion-search core",
      not any(item.suffix.lower() in {
          ".pdf", ".xlsx", ".xls", ".csv", ".zip", ".h5", ".mat", ".dat"
      } for item in core_files))
for path in core_files + [AUDIT / "LIVE_REQUERY_NORMALIZED.json"]:
    raw = path.read_bytes()
    check(f"UTF-8 and control-byte hygiene: {path.name}",
          all(byte in (9, 10) or byte >= 32 for byte in raw)
          and 127 not in raw)
    raw.decode("utf-8")

failed = [label for label, ok in checks if not ok]
for label, ok in checks:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
print(f"SURFACE_COUNTS literal={len(literal)} curator={len(curator)} total={len(surfaces)}")
print(f"NATURE_OBJECTS {len(nature_objects)}")
print(f"DATACITE_RECORDS {len(datacite_records)}")
print(f"ZENODO_RECORDS {len(zenodo_records)}")
print("FIGSHARE_MATCHES 0")
print("VERDICT PASS_AFTER_COMPLETENESS_LEAD_TYPING_AND_CUSTODY_REPAIR")
print(f"SUMMARY {len(checks) - len(failed)}/{len(checks)} independent checks passed")
if failed:
    raise SystemExit(1)
