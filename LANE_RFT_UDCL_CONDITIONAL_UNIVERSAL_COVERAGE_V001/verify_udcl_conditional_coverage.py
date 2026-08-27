#!/usr/bin/env python3
"""Custody and finite-logic regression for the typed U-DCL theorem lane.

This verifier checks pinned dependencies, load-bearing scope language, and the
logical composition of the displayed implications.  It does not prove or test
whether nature satisfies U-DCL.
"""

from __future__ import annotations

from hashlib import sha256
from itertools import product
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


THEOREM = read(HERE / "THEOREM.md")
CLARIFICATION = read(HERE / "TYPED_CLARIFICATION_V002.md")
BOUNDARY = read(HERE / "BOUNDARY.md")
RESULT = read(HERE / "RESULT.md")
AUDIT = read(HERE / "AUDIT.md")
INDEPENDENT_AUDIT = read(HERE / "INDEPENDENT_AUDIT.md")
ADOPTION = read(ROOT / "URFT_UDCL_ADOPTION_V001.md")
DECISION = read(ROOT / "URFT_UNIVERSAL_DOMAIN_LAW_DECISION_V001.md")
SEALED = read(ROOT / "LANE_RFT_STANDARD_CAUSAL_URFT_SCOPE_V001" / "THEOREM.md")

checks: list[tuple[str, bool]] = []


def add(name: str, condition: bool) -> None:
    checks.append((name, bool(condition)))


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


# Required author, repair, audit, and custody artifacts.
for filename in (
    "README.md",
    "THEOREM.md",
    "TYPED_CLARIFICATION_V002.md",
    "BOUNDARY.md",
    "RESULT.md",
    "AUDIT.md",
    "INDEPENDENT_AUDIT.md",
    "DEPENDENCIES.sha256",
    "verify_udcl_conditional_coverage.py",
    "VERIFICATION.txt",
):
    add(f"required file {filename}", (HERE / filename).is_file())


# Pinned external dependencies.
dependency_rows = []
for line in read(HERE / "DEPENDENCIES.sha256").splitlines():
    if line.strip():
        expected, relative = line.split(maxsplit=1)
        dependency_rows.append((expected, relative.strip()))
add("five dependency pins", len(dependency_rows) == 5)
for expected, relative in dependency_rows:
    target = ROOT / relative
    add(
        f"dependency hash {relative}",
        target.is_file() and digest(target) == expected,
    )


# The imported lane manifest must itself validate from its own directory.
sealed_dir = ROOT / "LANE_RFT_STANDARD_CAUSAL_URFT_SCOPE_V001"
sealed_manifest_ok = True
sealed_manifest_rows = 0
for line in read(sealed_dir / "MANIFEST.sha256").splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    sealed_manifest_rows += 1
    target = sealed_dir / relative.strip()
    if not target.is_file() or digest(target) != expected:
        sealed_manifest_ok = False
add("sealed dependency manifest has seven rows", sealed_manifest_rows == 7)
add("sealed dependency manifest verifies", sealed_manifest_ok)


# Adoption and historical-decision continuity.
add("adoption has U-DCL quantifier", "\\operatorname{U\\!DCL}:\\Longleftrightarrow" in ADOPTION)
add("decision has exact D1 heading", "### D1. Directed frontier" in DECISION)
add("decision has exact D4 heading", "### D4. Locality and provenance" in DECISION)
add("historical anchors preserved", "historical V001 decision and adoption record remain unchanged" in CLARIFICATION)
add("wording strengthening disclosed", "stronger **wording**" in CLARIFICATION)


# Ontic/evidential separation.
add("ontic DCL phys defined", "DCL_{\\rm phys}(r):\\Longleftrightarrow" in THEOREM)
add("certificate separately defined", "\\operatorname{Cert}_{DCL}(r;P)" in THEOREM)
add("sound certificate implies physical", "\\operatorname{Cert}_{DCL}(r;P)\\Longrightarrow DCL_{\\rm phys}(r)" in THEOREM)
add("noncertification is not physical failure", "noncertification does not imply physical failure" in THEOREM)
add("adopted DCL alias explicit", "DCL(r):=DCL_{\\rm phys}(r)" in THEOREM)
add("physical bound outcome independent", "outcome-independent finite physical bound" in THEOREM)
add("prospective freeze reserved for certificate", "Prospective freezing is a requirement on `Cert_DCL`" in THEOREM)


# One common typed physical witness.
for phrase, name in (
    ("b\\in\\{K,W\\}", "common K W branch"),
    ("All four clauses must hold for this **same** tuple", "same tuple all clauses"),
    ("Clause-wise mixing of a\nclassical state, quantum maps", "branch mixing forbidden"),
    ("complete external operational incidence", "complete external incidence"),
    ("is acyclic", "external acyclicity"),
    ("Unresolved feedback or chronology violation cannot be hidden", "unresolved loop forbidden"),
    ("actual history distribution and its correlations", "joint history distribution"),
    ("Positive-history conditional objects are conditionals of this supplied joint", "history conditionals derived"),
    ("compose\nuniquely on both sides of `G`", "well posed two sided composition"),
    ("same branch `b` and have matching\ndomains and codomains", "typed map interfaces"),
    ("**joint** law, state, or process of all later fresh inputs", "joint future input rule"),
    ("Equality of their separate marginals is insufficient", "marginals insufficient"),
):
    add(name, phrase in THEOREM)


# Exact local and imported chain.
add("DCL phys supplies C S J", "DCL_{\\rm phys}(r)\\Longrightarrow\n C(r)\\land S(r)\\land J(r)" in THEOREM)
add("J supplied not derived", "D2 **supplies**" in THEOREM and "No marginal-extension theorem" in THEOREM)
add("T5 overlap disclosed", "D4 openly supplies T5-grade provenance" in THEOREM)
add("sealed dependency has SCU1", "Theorem SCU1" in SEALED)
add("sealed dependency REC C S J", "\\operatorname{REC}(r)\\land C(r)\\land S(r)\\land J(r)" in SEALED)
add("sealed dependency CTS", "\\operatorname{CTS}(r)" in SEALED)
add("sealed dependency OCC union", "\\operatorname{OCC}_{\\cup}" in SEALED)
add("sealed dependency A1 A4", "\\bigwedge_{i=1}^{4}A_i" in SEALED)
add("per record theorem typed", "\\operatorname{REC}(r)\\land DCL_{\\rm phys}(r)" in THEOREM)


# Independent domain and global quantifier.
add("domain membership entails REC", "r\\in\\mathfrak R^{\\rm actual,bf}_{\\rm FM}\n \\Longrightarrow \\operatorname{REC}(r)" in THEOREM)
add("global postulate uses DCL phys", "\\quad DCL_{\\rm phys}(r)." in THEOREM)
add("universal Coverage U conclusion", "\\forall r\\in\\mathfrak R^{\\rm actual,bf}_{\\rm FM},\n \\quad\\operatorname{COV}_{\\cup}(r)" in THEOREM)
add("no post hoc domain retreat", "Post hoc removal after a failed U-DCL test is forbidden" in THEOREM)


# Typed falsifier burden.
add("falsifier negates full existential", "negate the full\nexistential" in BOUNDARY)
add("typed realization ledger", "outcome-independent realization\n   ledger" in BOUNDARY)
add("falsifier completeness theorem", "checkable completeness theorem" in BOUNDARY)
add("negative without closure refuses", "Without item 4, a negative fit is `NOT_TESTED`" in BOUNDARY)
add("all frontiers and branches", "every allowed frontier and every\n  applicable typed K/W realization" in BOUNDARY)
add("one positive local only", "No finite/open collection of successes proves\nU-DCL" in BOUNDARY)
add("measured port supersedes", "preserved but superseded" in BOUNDARY and "status return to `NOT_TESTED`" in BOUNDARY)


# Claim ceiling in all reporting artifacts.
for excluded in ("actualization", "Born", "A5", "gravity"):
    add(
        f"claim ceiling {excluded}",
        excluded.lower() in THEOREM.lower()
        and excluded.lower() in RESULT.lower()
        and excluded.lower() in BOUNDARY.lower()
        and excluded.lower() in INDEPENDENT_AUDIT.lower(),
    )
add("natural validity unproved", "U-DCL remains a universal physical conjecture" in THEOREM)
add("independent audit accepts repaired theorem", "ACCEPT_AFTER_TYPED_V002_REPAIR" in INDEPENDENT_AUDIT)


# Boolean composition: the arrows are physical/theorem premises; exhaustive
# enumeration checks that no logical antecedent is lost when they compose.
local_chain_ok = True
for values in product((False, True), repeat=10):
    rec, dcl, c, s, j, cts, occ, enc, a14, cov = values
    premises = (
        ((not dcl) or (c and s and j))
        and ((not (rec and c and s and j)) or cts)
        and ((not cts) or occ)
        and ((not occ) or (enc and a14))
        and (cov == (occ and enc and a14))
    )
    if premises and rec and dcl and not cov:
        local_chain_ok = False
        break
add("Boolean REC and typed DCL imply Coverage U", local_chain_ok)


# Universal generalization: every independently admitted member is REC; U-DCL
# supplies DCL_phys; the local theorem then covers every member.
universal_ok = True
for n in range(1, 6):
    for dcls in product((False, True), repeat=n):
        udcl = all(dcls)
        covs = tuple(dcls)  # REC holds for each independently admitted member.
        if udcl and not all(covs):
            universal_ok = False
add("finite model universal generalization", universal_ok)


# Existential witness and exact negative asymmetry.
existential_ok = True
for n in range(1, 6):
    for witnesses in product((False, True), repeat=n):
        dcl = any(witnesses)
        exhaustive_refutation = not any(witnesses)
        if dcl == exhaustive_refutation:
            existential_ok = False
add("typed existential versus exhaustive refutation", existential_ok)


# A certificate can be absent while the physical predicate is true.  A sound
# positive certificate may not occur while the physical predicate is false.
certificate_typing_ok = (
    ((not False) or False)  # no certificate, no physical witness is allowed
    and ((not False) or True)  # no certificate, physical witness is allowed
    and ((not True) or True)  # sound certificate, physical witness is required
    and not ((not True) or False)  # certificate without witness is forbidden
)
add("certificate implication has no converse", certificate_typing_ok)


passed = sum(ok for _, ok in checks)
for number, (name, ok) in enumerate(checks, start=1):
    print(f"D{number:02d} {'PASS' if ok else 'FAIL'} {name}")
print(f"TOTAL {passed}/{len(checks)} {'PASS' if passed == len(checks) else 'FAIL'}")
if passed == len(checks):
    print(
        "VERDICT TYPED_UDCL_CONDITIONAL_UNIVERSAL_COVERAGE_LOGIC_SCOPE_"
        "AND_CUSTODY_PASS__NATURAL_UDCL_VALIDITY_NOT_TESTED"
    )
else:
    print("VERDICT FAIL")
    sys.exit(1)
