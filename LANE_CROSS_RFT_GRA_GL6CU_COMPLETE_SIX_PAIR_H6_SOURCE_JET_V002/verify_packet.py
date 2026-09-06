#!/usr/bin/env python3
"""Fail-closed custody, exact replay, and scope verifier for GL6CU V002."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def unique_object(pairs):
    answer = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


required = {
    "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md",
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "EXACT_LEDGER.json",
    "VERIFICATION.txt", "derive_complete_six_pair_h6_source_jet.py",
    "verify_packet.py", "MANIFEST.sha256", "SEAL.sha256",
}
for name in sorted(required):
    check((HERE / name).is_file(), f"required packet file {name}")


# Exact frozen dependencies and repair provenance.
dependency_lines = [
    line for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines()
    if line.strip()
]
check(len(dependency_lines) == 16, "sixteen exact dependency bytes pinned")
dependency_names = set()
for line in dependency_lines:
    expected, relative = line.split("  ", 1)
    check(relative not in dependency_names, f"unique dependency {relative}")
    dependency_names.add(relative)
    candidate = ROOT / relative
    check(candidate.is_file(), f"dependency exists {relative}")
    check(digest(candidate) == expected, f"dependency hash {relative}")
for token in ("GL6CF", "GL6CG", "GL6CH", "GL6CN", "GL6BV", "GL6AO",
              "GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V001",
              "AUDIT_G_GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V001"):
    check(any(token in name for name in dependency_names),
          f"dependency family pinned {token}")
check(any(name.endswith(
      "GL6BV_STRICT_LOCK_T2_SOURCE_FESHBACH_CONTACT_V001/THEOREM.md")
      for name in dependency_names), "defining BV CTP theorem pinned")
check(any(name.endswith(
      "AUDIT_G_GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V001/AUDIT_REPORT.md")
      for name in dependency_names), "sealed V001 audit report pinned")


ledger = json.loads((HERE / "EXACT_LEDGER.json").read_text(),
                    object_pairs_hook=unique_object)
check(ledger["schema"] == "GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V002" and
      ledger["status"] == "PASS", "ledger V002 schema and status")
check(ledger["checks"] == 21849, "exact internal check census")

source = ledger["source_convention"]
check(source["hamiltonian"].startswith("H=U_d D+hW+sum") and
      source["fixed_objects"].startswith("P, z, M"),
      "plus-source convention and fixed objects")
check("selected Q4 locked branch" in source["retained_space"] and
      "not the full global locked-manifold matrix" in source["retained_space"],
      "retained Hilbert/state scope")
check(source["downstream_writer_sign"].startswith("for H=H0-JY, Y_A=-V_A"),
      "downstream writer sign bridge")
check("D2F[G_A,G_B]+DF[G_AB]" in source["canonical_chain_rule"] and
      "bare G_AB=0" in source["canonical_chain_rule"],
      "canonical second-chain-rule scope")
check(source["ctp_definition"].startswith(
      "Z[j+,j-]=Tr(T exp[-i integral H(j+)dt]") and
      "W=-i log Z" in source["ctp_definition"] and
      "sigma=+1 forward and -1 backward" in source["ctp_definition"],
      "defining CTP functional and branches")
check(source["ctp_branch_contact"].startswith(
      "direct same-branch connected contact is -sigma K''_AB delta(t-s)") and
      "forward ++ is -K''" in source["ctp_branch_contact"] and
      "backward -- is +K''" in source["ctp_branch_contact"] and
      "mixed direct contacts vanish" in source["ctp_branch_contact"],
      "branchwise CTP contact signs")

census = ledger["owner_census"]
check(census == {
    "active_alternating_cycle_transitions_in_selected_Q4_state": 64,
    "bare_node_sources": 128,
    "diagonal_nonbare_operator_jets_computed": 3840,
    "elementary_six_cycle": 256,
    "graph_wide_geometric_supports_enumerated": 4096,
    "inactive_six_cycle_supports_enumerated_not_differentiated": 192,
    "one_edge": 256,
    "ownership": "Möbius inversion by distinct flipped-edge support; one undirected cycle owns both Hermitian directions",
    "selected_row_nonbare_operator_jets_computed": 3904,
    "three_edge_path": 2304,
    "three_edge_star": 512,
    "total_nonbare_geometric_supports": 4096,
    "two_edge_wedge": 768,
}, "4096 enumerated supports versus 3904 computed selected-row Jets")

bare = ledger["bare_endpoint_source_selected_Q4"]
check(bare["value"] == "0" and bare["gradient_entries"] == 768 and
      bare["hessian_entries"] == 0 and
      sum(row["count"] for row in bare["local_pair_word_census"]) == 128,
      "complete source-independent-P endpoint source")
diagonal = ledger["diagonal_selected_Q4"]
check([diagonal[str(order)]["value"] for order in (2, 4, 6)] ==
      ["-128", "-224/3", "-28576/135"],
      "AO source-off diagonal values")

stencils = ledger["exact_local_operator_stencils"]
expected_owner_counts = {
    "one_edge_h2": 256, "one_edge_h4": 256, "one_edge_h6": 256,
    "wedge_h4": 768, "wedge_h6": 768,
    "star_h6": 512, "path_h6": 2304,
}
check(set(stencils) == set(expected_owner_counts),
      "every diagonal family/order record present")
for name, owner_count in expected_owner_counts.items():
    row = stencils[name]
    check(row["owner_count"] == owner_count and
          sum(int(multiplicity) * class_count
              for multiplicity, class_count in
              row["class_multiplicity_histogram"].items()) == owner_count and
          len(row["all_exact_stencils_sha256"]) == 64,
          f"complete class coverage {name}")
    check(sum(row["class_multiplicity_histogram"].values()) ==
          row["class_count"] and len(row["all_class_records_sha256"]) == 64 and
          all(len(entry["stencil_sha256"]) == 64 and entry["count"] > 0
              for entry in row["representative_class_records"]),
          f"class census and representative hashes {name}")
    check("first_exact_stencil" in row and
          len(row["representative_class_records"]) <= 4,
          f"compact literal/representative storage {name}")

cycle = ledger["cycle_writer_selected_Q4"]
check(cycle["active_count"] == 64 and
      cycle["source_free_coefficient"] == "-63/8" and
      "105/8" in cycle["raw_first_coefficient_per_cycle_node"],
      "active CH writer regressions")
cycle_stencil = cycle["exact_local_operator_stencils"][
    "alternating_cycle_offdiag_h6"]
check(cycle_stencil["owner_count"] == 64 and
      sum(cycle_stencil["class_multiplicity_histogram"].values()) ==
      cycle_stencil["class_count"], "active-cycle source-second classes")
check(min(cycle_stencil["hessian_entry_count_classes"]) > 0 and
      all(entry["hessian_entries"] > 0
          for entry in cycle_stencil["representative_class_records"]),
      "every cycle class has nonzero mixed second source")
check("first_exact_stencil" in cycle_stencil and
      len(cycle_stencil["representative_class_records"]) <= 4,
      "compact cycle literal/representative storage")

storage = ledger["stencil_storage_contract"]
check(storage == {
    "all_reconstructed_local_stencil_classes_committed": 1654,
    "custody": "compact ledger stores counts aggregate hashes representative metadata and one literal first exact stencil per family/order; sealed executable reconstructs all 1654 classes and verifies aggregate hashes",
    "family_order_records": 8,
    "literal_first_exact_stencils_stored": 8,
    "representative_records_per_family_order_maximum": 4,
}, "compact ledger storage contract")
check(sum(row["class_count"] for row in stencils.values()) +
      cycle_stencil["class_count"] == 1654,
      "all 1654 class commitments counted")

check(ledger["disconnected_owner_tests"] == {
    "three_edge_matching": "zero complete Jet",
    "two_edge_matching": "zero complete Jet",
    "wedge_plus_isolated": "zero complete Jet",
}, "disconnected linked-owner regressions")
check(any("not a stationary-state expectation" in row
          for row in ledger["ceilings"]) and
      any("no connected-to-1PI" in row for row in ledger["ceilings"]),
      "operator versus connected/1PI ceiling")


# Exact replay of the repaired target engine.
replay = subprocess.run(
    [sys.executable, "-B", str(HERE / "derive_complete_six_pair_h6_source_jet.py")],
    cwd=ROOT, capture_output=True, text=True, check=False,
)
check(replay.returncode == 0, "exact replay exits zero")
for token in (
    "PASS__GL6CU_V002_COMPLETE_SIX_PAIR_H6_SOURCE_JET__21849/21849",
    "CENSUS=4096_GEOMETRIC_SUPPORTS_ENUMERATED;3904_SELECTED_ROW_NONBARE_JETS_COMPUTED;192_INACTIVE_CYCLES",
    "OWNERS=128_BARE+256_LINK+768_WEDGE+512_STAR+2304_PATH+64_ACTIVE_CYCLE_JETS",
    "STENCILS=1654_CLASSES_HASH_COMMITTED_AND_EXECUTABLY_RECONSTRUCTED;8_LITERAL_FIRST_STENCILS_STORED",
    "DIAGONAL_VALUES=h2:-128,h4:-224/3,h6:-28576/135",
    "ACTIVE_CYCLE_WRITERS=64;H6=-63/8;FIRST_RAW=105/8;SECOND=NONZERO",
    "CTP=DEFINED_Z_PLUS_MINUS_AND_W_MINUS_I_LOG_Z;DIRECT_CONTACT=-SIGMA_ENERGY_HESSIAN_DELTA",
):
    check(token in replay.stdout, f"exact replay token {token}")


# Human-readable repair and ceiling custody.
theorem_raw = (HERE / "THEOREM.md").read_text()
result_raw = (HERE / "RESULT.md").read_text()
readme_raw = (HERE / "README.md").read_text()
self_raw = (HERE / "SELF_AUDIT.md").read_text()
theorem = " ".join(theorem_raw.split())
result = " ".join(result_raw.split())
readme = " ".join(readme_raw.split())
self_audit = " ".join(self_raw.split())
for token in (
    "effective-Hamiltonian operator theorem", "source-before-projection",
    "D_C^{-1}(k\\odot\\xi)", "des-Cloizeaux", "energy Hessian",
    "4096 **enumerated supports**", "=3904", "192 graph-wide cycle supports",
    "commits 1654 reconstructed local stencil classes",
    "same-branch coincident connected contact operator",
    "Z[j^+,j^-]", "W[j^+,j^-]=-i\\log Z[j^+,j^-]",
    "-\\sigma K''_{AB}\\,\\delta(t-s)", "gravity, and `G`",
):
    check(token in theorem, f"theorem repair/scope token {token}")
check("4096 graph-wide nonbare geometric supports" in result and
      "3904 selected-row nonbare operator Jets" in result and
      "192 six-cycle supports" in result,
      "result distinguishes supports and computed Jets")
check("commits 1654 exact local stencil classes" in result and
      "does not literally store every stencil coefficient" in result and
      "reconstructs all 1654 classes" in result,
      "result states compact ledger truthfully")
check("not a stationary connected response or a `1PI` kernel" in result and
      "`delta j` is not a Ward-field" in result,
      "result response and source/field ceilings")
check("4096 graph-wide geometric supports" in readme and
      "computes 3904 selected-row nonbare source Jets" in readme and
      "does not literally store every stencil coefficient" in readme,
      "readme census and compact-ledger repair")
check("4096 geometric supports" in self_audit and
      "3904 nonbare operator Jets" in self_audit and
      "does not literally store all 1654 coefficient tables" in self_audit,
      "self-audit repair boundaries")
check("does not replace the fresh independent hostile audit" in self_audit,
      "self-audit independence ceiling")
for forbidden in (
    "evaluates all 4096 nonbare geometric owners",
    "and every exact local nonuniform source stencil class",
    "connected local CTP contact is minus the displayed energy Hessian",
):
    check(forbidden not in result and forbidden not in readme and
          forbidden not in json.dumps(ledger),
          f"V001 overclaim absent {forbidden}")

verification = (HERE / "VERIFICATION.txt").read_text()
check("PASS__GL6CU_V002_COMPLETE_SIX_PAIR_H6_SOURCE_JET__21849/21849" in
      verification, "verification exact target count")
check("CENSUS=4096_GEOMETRIC_SUPPORTS_ENUMERATED;3904_SELECTED_ROW_NONBARE_JETS_COMPUTED;192_INACTIVE_CYCLES" in
      verification, "verification census repair")
check("HOSTILE_AUDIT=NOT_INCLUDED" in verification,
      "fresh hostile audit still required")


# Local immutable packet custody.
manifest_lines = [
    line for line in (HERE / "MANIFEST.sha256").read_text().splitlines()
    if line.strip()
]
manifest_names = set()
for line in manifest_lines:
    expected, name = line.split("  ", 1)
    check(Path(name).parent == Path("."), f"manifest local path {name}")
    check(name not in manifest_names, f"manifest unique path {name}")
    manifest_names.add(name)
    check((HERE / name).is_file(), f"manifest target exists {name}")
    check(digest(HERE / name) == expected, f"manifest hash {name}")
check(manifest_names == required - {"MANIFEST.sha256", "SEAL.sha256"},
      "manifest covers every non-custody packet byte")
seal = (HERE / "SEAL.sha256").read_text().strip().split("  ", 1)
check(seal == [digest(HERE / "MANIFEST.sha256"), "MANIFEST.sha256"],
      "seal names and hashes packet manifest")

print(f"PASS__GL6CU_V002_PACKET__{checks}/{checks}")
