#!/usr/bin/env python3
"""Exact replay for the bounded RMLB implication/countermodel packet."""

from fractions import Fraction as Q
from pathlib import Path


ROOT = Path(__file__).resolve().parent
passed = 0
total = 0
failures: list[str] = []


def check(condition: bool, name: str) -> None:
    global passed, total
    total += 1
    if condition:
        passed += 1
    else:
        failures.append(name)


def mm(a, b):
    return [
        [sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2)]
        for i in range(2)
    ]


def madd(a, b):
    return [[a[i][j] + b[i][j] for j in range(2)] for i in range(2)]


def msub(a, b):
    return [[a[i][j] - b[i][j] for j in range(2)] for i in range(2)]


def mscale(c, a):
    return [[c * a[i][j] for j in range(2)] for i in range(2)]


I = [[1, 0], [0, 1]]
sx = [[0, 1], [1, 0]]
sy = [[0, -1j], [1j, 0]]
sz = [[1, 0], [0, -1]]
Pxp = mscale(Q(1, 2), madd(I, sx))
Pzp = mscale(Q(1, 2), madd(I, sz))

check(mm(Pxp, Pxp) == Pxp, "X effect is a projection")
check(mm(Pzp, Pzp) == Pzp, "Z effect is a projection")
comm = msub(mm(Pzp, Pxp), mm(Pxp, Pzp))
expected_comm = mscale(0.5j, sy)
check(comm == expected_comm, "sharp X/Z commutator exact")
check(comm != [[0, 0], [0, 0]], "sharp X/Z incompatible")


# Exact noisy-record Fisher/QFI separation.
eta = Q(3, 5)
for t in [Q(0), Q(1, 3), Q(-1, 2), Q(4, 5)]:
    f_query = eta * eta / (1 - eta * eta * t * t)
    f_sld = 1 / (1 - t * t)
    gap = (1 - eta * eta) / ((1 - t * t) * (1 - eta * eta * t * t))
    check(f_sld - f_query == gap, f"Fisher gap identity t={t}")
    check(gap > 0, f"Fisher gap positive t={t}")
    check(f_query > 0, f"noisy record remains informative t={t}")

# Exact squared-fidelity witness: u=3/5 has state square-root term 4/5;
# choosing eta=25/39 makes eta*u=5/13 with query square-root term 12/13.
fid_u = Q(3, 5)
fid_eta = Q(25, 39)
gamma_state = (1 + Q(4, 5)) / 2
gamma_query = (1 + Q(12, 13)) / 2
check(fid_eta * fid_u == Q(5, 13), "fidelity witness eta*u exact")
check(gamma_state == Q(9, 10), "state squared fidelity exact")
check(gamma_query == Q(25, 26), "query squared fidelity exact")
check(gamma_query > gamma_state, "noisy query fidelity strictly exceeds state fidelity")


# A concrete Markov read kernel and data-processing check.
# Raw q in +/- has p(+|x)=(1+x)/2.  A binary symmetric read has correlation eta.
for x in [Q(-2, 3), Q(0), Q(1, 4), Q(3, 5)]:
    p_plus = (1 + x) / 2
    p_minus = (1 - x) / 2
    keep = (1 + eta) / 2
    flip = (1 - eta) / 2
    y_plus = keep * p_plus + flip * p_minus
    check(y_plus == (1 + eta * x) / 2, f"Markov factorization x={x}")
    f_raw = 1 / (1 - x * x)
    f_score = eta * eta / (1 - eta * eta * x * x)
    check(f_score <= f_raw, f"Fisher data processing x={x}")


# No one constant scale solders the binary Fisher metric to a constant g_xx.
f0 = Q(1)
fhalf = Q(4, 3)
check(f0 != fhalf, "binary Fisher metric varies with x")
ell2_from_zero = Q(1) / f0
ell2_from_half = Q(1) / fhalf
check(ell2_from_zero != ell2_from_half, "no common soldering scale")


theorem = (ROOT / "THEOREM.md").read_text()
normalized_theorem = " ".join(theorem.split())
required = [
    "DCL **completeness**",
    "Statistical **sufficiency**",
    "common **separator experiment**",
    "tomographically complete",
    "standard-Borel query POVM",
    "parameter-independent POVM",
    "F^{Y_m}\\preceq F^{Q_m}",
    "F^{Q_m}\\preceq F^{\\rm SLD}_{\\rho}",
    "[P^Z_+,P^X_+]",
    "F^{\\rm SLD}(t)-F^Q(t)",
    "\\gamma_Q(0,u)",
    "FERS-P3/P4",
    "no one constant information-to-length scale",
    "RLS -- record-localization soldering law",
    "logically independent",
    "not a definition of spacetime by information",
]
for snippet in required:
    check(snippet in normalized_theorem, f"scope snippet {snippet[:28]}")

check(theorem.count("\\[") == theorem.count("\\]"), "display math balanced")
check("DOES_NOT_PROVE_ONE_PHYSICALLY_READABLE_BLACKWELL_COMPLETE_CLASSICAL_QUERY" in theorem,
      "disposition retains query ceiling")
check("DOES_NOT_PROVE_INFORMATION_METRIC_EQUALS_CAUSAL_PROPER_METRIC" in theorem,
      "disposition retains soldering ceiling")

if failures:
    print(f"FAIL {passed}/{total}")
    for failure in failures:
        print(f" - {failure}")
    raise SystemExit(1)

print(f"PASS {passed}/{total}")
print("PER_INSTRUMENT_RECORD_MEDIATION_EXACT")
print("INCOMPATIBLE_RECORDED_QUERIES_REFUTE_UNIVERSAL_CLASSICAL_SUFFICIENCY")
print("NOISY_BONA_FIDE_RECORD_STRICTLY_FAILS_QFI_SATURATION")
print("METRIC_SPECTATOR_REFUTES_PROPER_GEOMETRY_DERIVATION")
print("RLS_COARSE_INTERFACE_TWO_CLAUSE_RESIDUAL_ISOLATED")
