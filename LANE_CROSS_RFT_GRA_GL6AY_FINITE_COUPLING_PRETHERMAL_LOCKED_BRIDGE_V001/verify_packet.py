#!/usr/bin/env python3
"""Fail-closed custody, inventory, source, and scope verifier for GL6AY."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise RuntimeError(label)
    checks += 1


required = {
    "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md",
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "PRIMARY_SOURCES.md",
    "VERIFICATION.txt", "verify_finite_coupling_prethermal_bridge.py",
    "verify_packet.py", "MANIFEST.sha256", "SEAL.sha256",
}
check({path.name for path in HERE.iterdir() if path.is_file()} == required,
      "author file set exact")


def verify_rows(path: Path, allowed_parent: Path | None = None):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        check(len(expected) == 64 and
              all(char in "0123456789abcdef" for char in expected),
              f"hash syntax: {relative}")
        check(relative not in rows, f"unique row: {relative}")
        target = ROOT / relative
        check(target.is_file() and not target.is_symlink(),
              f"regular target: {relative}")
        if allowed_parent is not None:
            check(target.parent == allowed_parent, f"target confined: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"hash matches: {relative}")
        rows.append(relative)
    return rows


dependencies = verify_rows(HERE / "DEPENDENCIES.sha256")
check(len(dependencies) == 18, "eighteen exact local dependencies")
for marker in ("GL6AN", "GL6AO", "GL6AX"):
    check(sum(marker in row for row in dependencies) == 6,
          f"six {marker} author/audit objects")
    check(sum(f"AUDIT_G_{marker}" in row for row in dependencies) == 3,
          f"three distinct {marker} audit objects")
for seal_relative in (row for row in dependencies if row.endswith("/SEAL.sha256")):
    seal_path = ROOT / seal_relative
    rows = [line for line in seal_path.read_text().splitlines() if line.strip()]
    check(len(rows) == 1, f"one dependency seal row: {seal_relative}")
    expected, manifest_relative = rows[0].split(maxsplit=1)
    manifest_path = ROOT / manifest_relative
    check(manifest_path.is_file() and not manifest_path.is_symlink(),
          f"dependency manifest regular: {manifest_relative}")
    check(hashlib.sha256(manifest_path.read_bytes()).hexdigest() == expected,
          f"dependency seal resolves: {seal_relative}")

manifest = verify_rows(HERE / "MANIFEST.sha256", HERE)
check(len(manifest) == 10, "author manifest has ten rows")
check({Path(row).name for row in manifest} == required - {
    "MANIFEST.sha256", "SEAL.sha256"
}, "author manifest file set exact")
seal = [line for line in (HERE / "SEAL.sha256").read_text().splitlines()
        if line.strip()]
check(len(seal) == 1, "one author seal row")
expected, relative = seal[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "seal targets manifest")
check(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest()
      == expected, "seal resolves")

theorem = " ".join((HERE / "THEOREM.md").read_text().split())
readme = " ".join((HERE / "README.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
sources = " ".join((HERE / "PRIMARY_SOURCES.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())

for document in (theorem, readme):
    check("author frozen and sealed" in document, "frozen author status")
check("distinct independent hostile audit" in theorem.lower() or
      "distinct independent hostile audit" in readme.lower(),
      "distinct hostile audit required")

for token in (
    "spec(q_v^2)={0,1,4} subset Z",
    "A potential term `Z_S` is **strongly supported**",
    "||Z||^str_kappa",
    "U_d >=9pi v_0/kappa_0",
    "n_* =floor{(U_d/nu_0)/[1+ln(U_d/nu_0)]^3}-2",
    "Y_L H_L Y_L^*",
    "||V_hat_L||^str_(kappa_*)",
    "An absent-edge addition is resonant",
    "There is no corresponding global spectral projection",
    "N_S=sum_(v: supp(q_v^2) subset S)q_v^2",
    "P_S^0=chi(N_S=0)",
    "Phi(S)=P_S^0 D_hat(S)P_S^0",
    "[D_hat(S),N_S]=0",
    "n,n' globally locked",
    "P_L D_hat_L P_L",
    "exact termwise port U(1)^4 there",
    "P_S^0[A_S,[A_S,D_hat_L(S)]]P_S^0",
    "sup_(m>=1){m^4 exp(-kappa_*m)}",
    "-(63/8)(h^6/U_d^5)sum_c T_c",
    "0<r_1<ln(3/2)/4",
    "t_obs<=t_*",
    "local P_L -> Q_L leakage",
    "port-changing P_L -> P_L return",
    "Neither statement bounds",
    "no global dressed spectral subspace exists",
    "A_C^(r)!=0",
    "intermediate return to `P_L`",
    "gap from `P_L` to `Q_L`",
    "Q_L(E-Q_L H_L Q_L)^(-1)Q_L",
    "descendant of `P_L`",
    "||W_L||=h|E_L|",
    "failure of the standard whole-band proof route",
    "record-authenticated finite-horizon application theorem",
):
    check(token in theorem, f"theorem scope token: {token}")

for token in (
    "arXiv:`1509.05386v3`",
    "https://arxiv.org/html/1509.05386v3",
    "arXiv:`1704.08703v2`",
    "https://arxiv.org/html/1704.08703v2",
    "arXiv:`1105.0675v1`",
    "https://arxiv.org/html/1105.0675v1",
    "strong support",
    "fixed-order linked truncation",
):
    check(token in sources, f"primary source pin: {token}")

check("excursion `P_L -> Q_L`" in result and
      "returning `P_L -> P_L` process" in result,
      "result keeps leakage/winding distinction")
check("Phi_S=P_S^0 D_hat(S)P_S^0" in result,
      "result uses local collar interaction")
check("not the norm distance of a global" in result,
      "result rejects global dressed-subspace inference")
check("physical clock calibration open" in self_audit,
      "self-audit keeps clock calibration open")
check("PASS__GL6AY_FINITE_COUPLING_PRETHERMAL_LOCKED_BRIDGE__" in verification,
      "constructive replay recorded")

aggregate = " ".join((theorem, readme, result, self_audit)).lower()
for forbidden in (
    "intermediate return to `p`,",
    "gap from `p` to `q`",
    "q(e-qh_lq)^(-1)q",
    "descendant of `p`",
    "phi_l(s)=p d_hat_l(s)p",
    "boxed: p d_hat p",
    "physical dressed space `y^*p`",
    "its deviation from the bare description is controlled",
    "proves exact finite-coupling locked phase",
    "proves the remainder vanishes",
    "proves selected-gns gaplessness",
    "proves an isotropic cone",
    "is a physical photon",
    "is a graviton",
    "is gravity",
    "derives the einstein equation",
    "derives newton's constant",
):
    check(forbidden not in aggregate, f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6AY_PACKET__{checks}/{checks}")
