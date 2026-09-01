#!/usr/bin/env python3
"""Fail-closed verifier for the independent repaired-byte GL6AY audit."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6AY_FINITE_COUPLING_PRETHERMAL_LOCKED_BRIDGE_V001"
PRIOR = ROOT / "AUDIT_G_GL6AY_FINITE_COUPLING_PRETHERMAL_LOCKED_BRIDGE_V001"
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise RuntimeError(label)
    checks += 1


required = {
    "README.md",
    "AUDIT.md",
    "SOURCE_REPLAY.md",
    "AUDITED_TARGETS.sha256",
    "PRIOR_AUDIT_CUSTODY.sha256",
    "VERIFICATION.txt",
    "independent_gl6ay_post_repair.py",
    "verify_audit_packet.py",
    "MANIFEST.sha256",
    "SEAL.sha256",
}
check({path.name for path in HERE.iterdir() if path.is_file()} == required,
      "post-repair audit inventory exact")


def verify_rows(path: Path, allowed_parent: Path | None = None):
    answer = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split(maxsplit=1)
        check(len(fields) == 2, f"two hash-row fields: {line}")
        expected, relative = fields
        check(len(expected) == 64 and
              all(char in "0123456789abcdef" for char in expected),
              f"lowercase sha256 syntax: {relative}")
        check(relative not in answer, f"unique hash row: {relative}")
        target = ROOT / relative
        check(target.is_file() and not target.is_symlink(),
              f"regular nonsymlink target: {relative}")
        if allowed_parent is not None:
            check(target.parent == allowed_parent,
                  f"hash target confined: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"hash resolves: {relative}")
        answer.append(relative)
    return answer


# Repaired author custody is exact and exhaustive.
targets = verify_rows(HERE / "AUDITED_TARGETS.sha256", AUTHOR)
check(len(targets) == 12, "all twelve repaired author files pinned")
check({Path(item).name for item in targets} ==
      {path.name for path in AUTHOR.iterdir() if path.is_file()},
      "pinned author set equals frozen author inventory")
check(hashlib.sha256((AUTHOR / "THEOREM.md").read_bytes()).hexdigest() ==
      "5b86ab5eb2998eb719dffd09e05add131863fd2a3290d87fb749dc8aebc1891c",
      "final repaired theorem hash")
author_manifest_hash = hashlib.sha256(
    (AUTHOR / "MANIFEST.sha256").read_bytes()).hexdigest()
check(author_manifest_hash ==
      "e81ec1cfd4bdcdc43b4709b8f90f9eceac3dfba82be80701dd4a2a7e08de089b",
      "final repaired author manifest hash")
check(hashlib.sha256((AUTHOR / "SEAL.sha256").read_bytes()).hexdigest() ==
      "740f051b3347d7387e481a9991f536bd61a7e47ad51d80b680475dff394e5cbb",
      "final repaired author seal-file hash")
author_seal = [line for line in (AUTHOR / "SEAL.sha256").read_text().splitlines()
               if line.strip()]
check(len(author_seal) == 1, "one repaired author seal row")
author_seal_fields = author_seal[0].split(maxsplit=1)
check(author_seal_fields == [author_manifest_hash,
      f"{AUTHOR.name}/MANIFEST.sha256"], "repaired author seal resolves")

# The original FAIL audit is preserved by exact additive custody.
prior = verify_rows(HERE / "PRIOR_AUDIT_CUSTODY.sha256", PRIOR)
check(len(prior) == 3, "three original FAIL-audit custody objects")
check({Path(item).name for item in prior} ==
      {"AUDIT.md", "MANIFEST.sha256", "SEAL.sha256"},
      "original FAIL-audit custody set exact")
check(hashlib.sha256((PRIOR / "MANIFEST.sha256").read_bytes()).hexdigest() ==
      "687755718fa49553f33a738f7313d6ad89ae62b026be761f47860317ce9497e8",
      "original FAIL-audit manifest frozen")
prior_seal = [line for line in (PRIOR / "SEAL.sha256").read_text().splitlines()
              if line.strip()]
check(len(prior_seal) == 1, "one original FAIL-audit seal row")
prior_fields = prior_seal[0].split(maxsplit=1)
check(prior_fields == [
    "687755718fa49553f33a738f7313d6ad89ae62b026be761f47860317ce9497e8",
    f"{PRIOR.name}/MANIFEST.sha256",
], "original FAIL-audit seal resolves")

# This additive audit packet is itself exhaustive and sealed.
manifest = verify_rows(HERE / "MANIFEST.sha256", HERE)
check(len(manifest) == 8, "eight post-repair audit content rows")
check({Path(item).name for item in manifest} == required - {
    "MANIFEST.sha256", "SEAL.sha256"
}, "post-repair audit manifest inventory exact")
seal_rows = [line for line in (HERE / "SEAL.sha256").read_text().splitlines()
             if line.strip()]
check(len(seal_rows) == 1, "one post-repair audit seal row")
seal_fields = seal_rows[0].split(maxsplit=1)
check(len(seal_fields) == 2, "two post-repair seal fields")
check(seal_fields[1] == f"{HERE.name}/MANIFEST.sha256",
      "post-repair seal targets own manifest")
check(seal_fields[0] ==
      hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest(),
      "post-repair audit seal resolves")

audit = " ".join((HERE / "AUDIT.md").read_text().split())
readme = " ".join((HERE / "README.md").read_text().split())
source = " ".join((HERE / "SOURCE_REPLAY.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())
theorem = " ".join((AUTHOR / "THEOREM.md").read_text().split())
result = " ".join((AUTHOR / "RESULT.md").read_text().split())
self_audit = " ".join((AUTHOR / "SELF_AUDIT.md").read_text().split())
author_docs = " ".join((theorem, result, self_audit,
                        " ".join((AUTHOR / "README.md").read_text().split())))

disposition = (
    "PASS__FINITE_COUPLING_PRETHERMAL_BRIDGE__LOCAL_COLLAR_TYPED__"
    "GLOBAL_PROJECTOR_CLAIMS_REMOVED"
)
check(disposition in audit and disposition in readme,
      "post-repair PASS disposition explicit")
for token in (
    "F3 strong-support mapping — pass",
    "Finite local collar — pass",
    "[D_hat(S),N_S]=0",
    "without any global infinite-volume projector",
    "Locked-endpoint port conservation — pass",
    "Second twist moment — pass",
    "Finite-order contact and finite-horizon dynamics — pass",
    "Leakage, winding, and whole-band boundary — pass",
    "Hostile verdict: PASS",
):
    check(token in audit, f"hostile-audit finding present: {token}")

for token in (
    "1509.05386v3", "Theorem 3.1", "Theorem 3.3",
    "1704.08703v2", "Appendix A",
    "1105.0675v1", "Section 4",
):
    check(token in source, f"exact source replay token: {token}")

for token in (
    "There is no corresponding global spectral projection",
    "N_S=sum_(v: supp(q_v^2) subset S)q_v^2",
    "P_S^0=chi(N_S=0)",
    "Phi(S)=P_S^0 D_hat(S)P_S^0",
    "[D_hat(S),N_S]=0",
    "n,n' globally locked",
    "P_L D_hat_L P_L",
    "exact termwise port U(1)^4 there",
    "P_S^0[A_S,[A_S,D_hat_L(S)]]P_S^0",
    "`||Y_L^*P_LY_L-P_L||` uniformly in volume",
    "no global dressed spectral subspace exists",
    "intermediate return to `P_L`",
    "gap from `P_L` to `Q_L`",
    "Q_L(E-Q_L H_L Q_L)^(-1)Q_L",
    "descendant of `P_L`",
):
    check(token in theorem, f"repaired theorem scope token: {token}")

lower_author = author_docs.lower()
for forbidden in (
    "phi_l(s)=p d_hat_l(s)p",
    "boxed: p d_hat p",
    "physical dressed space `y^*p`",
    "intermediate return to `p`,",
    "gap from `p` to `q`",
    "q(e-qh_lq)^(-1)q",
    "descendant of `p`,",
    "proves exact finite-coupling locked phase",
    "proves selected-gns gaplessness",
    "proves an isotropic cone",
    "is a graviton",
    "is gravity",
    "derives the einstein equation",
    "derives newton's constant",
):
    check(forbidden not in lower_author,
          f"forbidden repaired-author regression absent: {forbidden}")

for token in (
    "PASS__GL6AY_FINITE_COUPLING_PRETHERMAL_LOCKED_BRIDGE__11145/11145",
    "PASS__GL6AY_PACKET__219/219",
    "PASS__INDEPENDENT_GL6AY_POST_REPAIR__32095/32095",
    "PASS__GL6AY_POST_REPAIR_HOSTILE_AUDIT_PACKET__",
):
    check(token in verification, f"verification replay recorded: {token}")

print(f"PASS__GL6AY_POST_REPAIR_HOSTILE_AUDIT_PACKET__{checks}/{checks}")
