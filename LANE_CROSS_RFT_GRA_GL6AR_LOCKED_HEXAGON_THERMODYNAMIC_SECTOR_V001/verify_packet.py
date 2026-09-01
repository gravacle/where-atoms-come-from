#!/usr/bin/env python3
"""Fail-closed custody and scope verifier for the GL6AR author packet."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


required = {
    "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md",
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "VERIFICATION.txt",
    "verify_locked_hexagon_thermodynamics.py", "verify_packet.py",
    "MANIFEST.sha256", "SEAL.sha256",
}
check({path.name for path in HERE.iterdir() if path.is_file()} == required,
      "packet directory file set exact")


def verify_rows(path, allowed_parent=None):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        check(len(expected) == 64 and all(char in "0123456789abcdef" for char in expected),
              f"hash syntax: {relative}")
        check(relative not in rows, f"unique row: {relative}")
        target = ROOT / relative
        check(target.is_file(), f"target exists: {relative}")
        if allowed_parent is not None:
            check(target.parent == allowed_parent, f"target confined: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"hash matches: {relative}")
        rows.append(relative)
    return rows


dependencies = verify_rows(HERE / "DEPENDENCIES.sha256")
check(len(dependencies) == 12, "exact dependency row count")
check(all("GL6AO" in row or "GL6AN" in row for row in dependencies),
      "dependencies confined to sealed AO/AN")
for lane in ("GL6AO", "GL6AN"):
    check(sum(lane in row for row in dependencies) == 6,
          f"six author/audit custody rows for {lane}")
check(any(row.endswith("AUDIT_G_GL6AO_COMPLETE_SIXTH_ORDER_LOCKED_HAMILTONIAN_V001/SEAL.sha256")
          for row in dependencies), "AO audit seal pinned")
check(any(row.endswith("AUDIT_G_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/SEAL.sha256")
          for row in dependencies), "AN audit seal pinned")

manifest = verify_rows(HERE / "MANIFEST.sha256", HERE)
check(len(manifest) == 9, "exact author manifest count")
check({Path(row).name for row in manifest} == required - {"MANIFEST.sha256", "SEAL.sha256"},
      "manifest file set exact")

seal_rows = [line for line in (HERE / "SEAL.sha256").read_text().splitlines()
             if line.strip()]
check(len(seal_rows) == 1, "one seal row")
expected, relative = seal_rows[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "seal targets manifest")
check(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest() == expected,
      "seal hash")

theorem = " ".join((HERE / "THEOREM.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())
for token in (
    "Omega(D,eta)={n locked",
    "sup_e sum_{c:e in supp(tau_c)} ||Phi(c)|| <= 18t",
    "sum_{i=0}^2 [s_(x,i)-s_(x-e_i,i)] = 0",
    "Every connected component of the finite flip graph",
    "Delta_C=t[rho_C-lambda_2(C)]>0",
    "Fixed-boundary exhaustion theorem",
    "does not prove retention of a nonlocal asymptotic sector",
    "-4tL^3 <= E_0(Q_L) <= -(t/64)L^3",
    "convergence of the energy density or existence of a thermodynamic energy-density limit is not proved",
    "Delta_L <= 18t ||w||_infinity^2 L / V_L",
    "If `V_L>0`",
    "unconditional selected-component finite-size gap closure is **not yet proved or refuted**",
    "no global infinite locked projector is inserted",
    "No completeness claim is made",
    "No fourth geometric cut has been derived",
    "projection onto the full zero-energy subspace",
    "mathcal E_omega(A_R)=<xi_R,H_omega xi_R> -> 0",
    "Merely requiring `omega(A_R)=0` would be insufficient",
):
    check(token in theorem, f"theorem scope token: {token}")
for token in (
    "conditional on an extensive variance lower bound",
    "not physical position or momentum",
    "not called a photon or graviton",
    "Finite-size component-gap closure is not promoted to thermodynamic/GNS spectral gaplessness",
):
    check(token in theorem + " " + self_audit, f"promotion ceiling: {token}")
check("finite-size component-gap closure is neither proved nor refuted" in result,
      "result preserves finite-size closure obstruction")
check("would still not by itself prove that a selected infinite-volume GNS ground state is spectrally gapless" in result,
      "result preserves infinite-volume spectral bridge")
for token in (
    "PASS__GL6AR_LOCKED_HEXAGON_THERMODYNAMICS__72161/72161",
    "SOFT=DELTA_LE_18_T_L_OVER_VARIANCE;EXTENSIVE_VARIANCE_IMPLIES_SELECTED_COMPONENT_GAP_L_MINUS2",
    "OBSTRUCTION=NO_DERIVED_VARIANCE_OR_GROUND_SECTOR;NO_GNS_SPECTRAL_BRIDGE",
    "SECTORS=THREE_COORDINATE_CUT_FLUXES_PLUS_DEPENDENT_FOURTH_PORT_COUNT;COMPONENTS_REFINE_INVARIANTS",
    "GROUND=FIXED_BOUNDARY_EXHAUSTION_WEAKSTAR_LIMIT;FINITE_PERIODIC_ACTIVE_ENERGY_BOUND_4_DIVIDES_L",
    "PASS__GL6AR_PACKET__",
):
    check(token in verification, f"verification token: {token}")

combined = (theorem + " " + result + " " + self_audit).lower()
for forbidden in (
    "the flux is electromagnetism",
    "the mode is a photon",
    "the mode is a graviton",
    "proves unconditional gaplessness",
    "proves gns gaplessness",
    "proves thermodynamic gaplessness",
    "four independent cut fluxes",
    "derives gravity",
    "derives gravitational coupling",
    "derives newton's constant",
):
    check(forbidden not in combined, f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6AR_PACKET__{checks}/{checks}")
