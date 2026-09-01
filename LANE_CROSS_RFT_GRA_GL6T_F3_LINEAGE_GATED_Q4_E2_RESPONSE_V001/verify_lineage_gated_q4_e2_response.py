#!/usr/bin/env python3
"""Exact algebra replay for GL6T's lineage-gated q4 response."""

from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parent
checks: list[str] = []


def check(label: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(label)
    checks.append(label)


edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
line = [[F(0) for _ in edges] for _ in edges]
for i, e in enumerate(edges):
    for j, f in enumerate(edges):
        if i != j and len(set(e) & set(f)) == 1:
            line[i][j] = F(1)

check("line degree four", all(sum(row) == 4 for row in line))

# Exact rational witness h=2, Delta=3, theta=pi/3.
h = F(2)
x = F(6, 25)
z = F(17, 25)
diag = -8 * h * x
adj = -4 * h * x * z * z
D = [[diag if i == j else adj * line[i][j] for j in range(6)] for i in range(6)]

check("diag witness", diag == F(-96, 25))
check("adj witness", adj == F(-13872, 15625))

# S4 edge-representation probes.
v_A = [F(1)] * 6
v_E1 = [F(1), F(-1), F(0), F(0), F(-1), F(1)]
v_E2 = [F(1), F(1), F(-2), F(-2), F(1), F(1)]
v_T1 = [F(1), F(0), F(0), F(0), F(0), F(-1)]
v_T2 = [F(0), F(1), F(0), F(0), F(-1), F(0)]
v_T3 = [F(0), F(0), F(1), F(-1), F(0), F(0)]


def matvec(M, v):
    return [sum(M[i][j] * v[j] for j in range(len(v))) for i in range(len(M))]


def eigen_check(label, v, lam):
    check(label, matvec(D, v) == [lam * y for y in v])


lam_A = F(-115488, 15625)
lam_E = F(-32256, 15625)
lam_T = F(-96, 25)
eigen_check("A1 eigen", v_A, lam_A)
eigen_check("E1 eigen", v_E1, lam_E)
eigen_check("E2 eigen", v_E2, lam_E)
eigen_check("T1 eigen", v_T1, lam_T)
eigen_check("T2 eigen", v_T2, lam_T)
eigen_check("T3 eigen", v_T3, lam_T)
check("all negative", lam_A < 0 and lam_E < 0 and lam_T < 0)

det_sector = lam_A * lam_E**2 * lam_T**3
det_formula = (8 * h * x) ** 6 * (1 + 2 * z * z) * (1 - z * z) ** 2
check("det sector formula", det_sector == det_formula and det_formula > 0)

# Direct exact determinant by fraction-preserving elimination.
def determinant(M):
    A = [row[:] for row in M]
    out = F(1)
    for i in range(len(A)):
        pivot = next((r for r in range(i, len(A)) if A[r][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            A[i], A[pivot] = A[pivot], A[i]
            out = -out
        q = A[i][i]
        out *= q
        for r in range(i + 1, len(A)):
            f = A[r][i] / q
            for c in range(i + 1, len(A)):
                A[r][c] -= f * A[i][c]
    return out


check("direct determinant", determinant(D) == det_formula)
check("break zero", determinant([[F(0)] * 6 for _ in range(6)]) == 0)

ledger = json.loads((HERE / "LINEAGE_Q4_RESPONSE_LEDGER.json").read_text())
check("four keep records", ledger["matched_routes"]["KEEP"]["K"] == 4)
check("four break quarantine", ledger["matched_routes"]["BREAK"]["G"] == 4)
check("N0 parent", "FPSS_N0" in ledger["selected_parent"])
check("complete instrument", "no_success_filter" in ledger["instrument_scope"])
check("blank chronology", "prewait_is_first" in ledger["prewait_chronology"])
check("K inclusive query", "K_and_n" in ledger["record_query"])
check("factorized dynamics ceiling", "no_interlink_interaction" in ledger["dynamical_scope"])
check("E2 ledger", ledger["s4_eigenvalues"]["E2"] == "-8*h*x*(1-z^2)")
check("no G", "no_G" in ledger["ceilings"])

theorem = (HERE / "THEOREM.md").read_text()
for phrase in (
    "D^{\\rm KEEP}(\\tau)",
    "-8hx\\,I_6-4hxz^2A_L",
    "D_{E_2}&=-8hx(1-z^2)",
    "\\det D^{\\rm KEEP}",
    "D^{\\rm BREAK}(\\tau)=0",
    "D_{E_2}&=-{32h^4\\Delta\\tau^4",
    "Specialize FPSS to `N=0`",
    "twelve explicit",
    "not a success-filtered unconditional experiment",
    "omit both the optional FPSS saturation pulse `U_KX`",
    "H_{\\rm full}=H_\\star",
    "complete common projective query",
    "Q^{K_a}_{F_a,u_{-a},{\\rm KEEP}}",
    "No such product is asserted after mixing",
    "does not prove an inter-link interaction",
    "It does not say that the active pair operators are records",
    "Nor is this a provenance-only force",
    "not called positive Ricci stiffness",
):
    check(f"theorem phrase {phrase}", phrase in theorem)

for line_hash in (HERE / "DEPENDENCIES.sha256").read_text().splitlines():
    expected, rel = line_hash.split("  ", 1)
    target = REPO / rel
    check(f"dependency exists {rel}", target.is_file())
    actual = hashlib.sha256(target.read_bytes()).hexdigest()
    check(f"dependency hash {rel}", actual == expected)

print(f"PASS__GL6T_LINEAGE_GATED_Q4_E2_RESPONSE__{len(checks)}/{len(checks)}")
