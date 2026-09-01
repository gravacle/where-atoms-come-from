#!/usr/bin/env python3
"""Exact finite-series replay for GL6U's inherited degree interaction."""

from fractions import Fraction as F
from math import factorial
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


def mv(A, v):
    return [sum(A[i][j] * v[j] for j in range(len(v))) for i in range(len(A))]


def mm(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]


def gconj(z):
    return (z[0], -z[1])


def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1],
            a[0] * b[1] + a[1] * b[0])


def gadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


# Exact rational witness U21 on the complete 16-state active-star space.
h, Delta, U, dstar = F(1), F(13), F(1), F(2)
dimension = 16
H = [[F(0) for _ in range(dimension)] for _ in range(dimension)]
for state in range(dimension):
    bits = [(state >> a) & 1 for a in range(4)]
    degree = sum(bits)
    H[state][state] = (
        Delta * degree
        + U * ((degree - dstar) ** 2
               + sum((bit - dstar) ** 2 for bit in bits))
    )
    for a in range(4):
        H[state][state ^ (1 << a)] -= h

v0 = [F(0)] * dimension
v0[0] = F(1)
powers = [v0]
for _ in range(6):
    powers.append(mv(H, powers[-1]))

# c_r=(-i)^r H^r/r!, represented as Gaussian rationals (real,imag).
minus_i_power = [(F(1), F(0)), (F(0), F(-1)),
                 (F(-1), F(0)), (F(0), F(1))]
coefficients = []
for order in range(7):
    re, im = minus_i_power[order % 4]
    coefficients.append([
        (re * value / factorial(order), im * value / factorial(order))
        for value in powers[order]
    ])


def op_x(a):
    out = [[F(0) for _ in range(dimension)] for _ in range(dimension)]
    for state in range(dimension):
        out[state ^ (1 << a)][state] = F(1)
    return out


def op_z(a):
    out = [[F(0) for _ in range(dimension)] for _ in range(dimension)]
    for state in range(dimension):
        out[state][state] = F(1) if ((state >> a) & 1) == 0 else F(-1)
    return out


def expectation_series(operator, max_order=6):
    out = []
    for order in range(max_order + 1):
        total = (F(0), F(0))
        for left_order in range(order + 1):
            right_order = order - left_order
            for i in range(dimension):
                left = gconj(coefficients[left_order][i])
                if left == (F(0), F(0)):
                    continue
                for j in range(dimension):
                    weight = operator[i][j]
                    if weight == 0:
                        continue
                    right = coefficients[right_order][j]
                    total = gadd(total, gmul(left, (weight * right[0],
                                                     weight * right[1])))
        check(f"real expectation order {order}", total[1] == 0)
        out.append(total[0])
    return out


def convolve(a, b, max_order=6):
    return [sum(a[p] * b[order - p] for p in range(order + 1))
            for order in range(max_order + 1)]


X0 = op_x(0)
Z0, Z1, Z2 = op_z(0), op_z(1), op_z(2)
XZZ = mm(mm(X0, Z1), Z2)
x_series = expectation_series(X0)
z_series = expectation_series(Z0)
y_series = expectation_series(XZZ)
z2_series = convolve(z_series, z_series)
xz2_series = convolve(x_series, z2_series)
x_minus_y = [x_series[r] - y_series[r] for r in range(7)]
factor_defect = [y_series[r] - xz2_series[r] for r in range(7)]

delta1 = Delta + 2 * U * (1 - 2 * dstar)
delta2 = 2 * Delta + 2 * U * (3 - 4 * dstar)
check("delta1 witness", delta1 == 7)
check("delta2 witness", delta2 == 16)
check("gap identity", 2 * delta1 - delta2 == -2 * U)
check("x order two", x_series[2] == h * delta1)
check("z order two", z_series[2] == -2 * h * h)
check("x-y lower orders zero", x_minus_y[:4] == [F(0)] * 4)
check("x-y order four",
      x_minus_y[4] == F(4, 3) * h**3 * (2 * delta2 - delta1))
check("factor defect lower orders zero", factor_defect[:4] == [F(0)] * 4)
check("factor defect order four", factor_defect[4] == -F(16, 3) * h**3 * U)
check("A1 leading", -24 * h * h * delta1 == -168)
check("T2 leading", -8 * h * h * delta1 == -56)
check("E2 leading", -F(32, 3) * h**4 * (2 * delta2 - delta1) == -F(800, 3))

# Exact L(K4) sector census for U12-U13.
edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
line = [[F(int(i != j and len(set(e) & set(f)) == 1))
         for j, f in enumerate(edges)] for i, e in enumerate(edges)]
check("line degree", all(sum(row) == 4 for row in line))
vA = [F(1)] * 6
vE = [F(1), F(-1), F(0), F(0), F(-1), F(1)]
vT = [F(1), F(0), F(0), F(0), F(0), F(-1)]
check("A eigen four", mv(line, vA) == [4 * q for q in vA])
check("E eigen minus two", mv(line, vE) == [-2 * q for q in vE])
check("T eigen zero", mv(line, vT) == [F(0)] * 6)

ledger = json.loads((HERE / "INTERACTING_RESPONSE_LEDGER.json").read_text())
check("native degree term", "degree" in ledger["restored_native_term"])
check("E2 ledger", ledger["s4_eigenvalues"]["E2"] == "-8*h*(x-y)")
check("no G", "no_G" in ledger["ceilings"])

theorem = (HERE / "THEOREM.md").read_text()
for phrase in (
    "delta_r=r\\Delta+U_dr(r+1-4d_\\star)",
    "D^{\\rm KEEP}=-8hxI_6-4hyA_L",
    "D_{E_2}=-8h(x-y)",
    "\\mathfrak C_{XZZ}:=y-xz^2",
    "-{16\\over3}h^3U_ds^4",
    "2\\delta_1-\\delta_2=-2U_d",
    "not called a full connected cumulant",
    "not yet a collective thermodynamic phase",
):
    check(f"theorem phrase {phrase}", phrase in theorem)

for line_hash in (HERE / "DEPENDENCIES.sha256").read_text().splitlines():
    expected, rel = line_hash.split("  ", 1)
    target = REPO / rel
    check(f"dependency exists {rel}", target.is_file())
    actual = hashlib.sha256(target.read_bytes()).hexdigest()
    check(f"dependency hash {rel}", actual == expected)

print(f"PASS__GL6U_DEGREE_INTERACTION_Q4_RESPONSE__{len(checks)}/{len(checks)}")
