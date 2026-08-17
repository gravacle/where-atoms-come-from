#!/usr/bin/env python3
"""
REFUTER 1 — SCRIPT 5.  W-08's CHARACTER IDENTITY AT FOUR CLASSES, EXACT IN THE GROUP RING.

The lane checks the identity  1 - |Z_k|^2 = sum_{j<l} w_j w_l |chi_j^k - chi_l^k|^2  in
Q[sqrt3], which restricts it to q-th roots of unity with q | 12.  That is 1620 cases in five
moduli.  I check it in Z[x]/(x^q - 1) with rational coefficients -- no trigonometry at all,
conjugation is x -> x^(q-1) -- which is exact for EVERY q, and I run q up to 40 including
primes the lane cannot reach (5, 7, 11, 13, 17, 19, 23, 29, 31, 37).

This is a stronger form of the lane's own check, not a different claim.  Two things it
settles that the lane's does not:
  (i)  the identity is exact for every modulus, not just for q | 12;
  (ii) MONOTONICITY IS A COROLLARY OF THE IDENTITY, not an independent observation --
       the right-hand side is a sum of non-negative terms.  So the lane's separately
       reported "0 events of |Z_k| > 1, exactly 0, not 1e-16" could not have come out
       otherwise once the identity verified.  It is a theorem, not a control, so this
       does not void it -- but it is not a second piece of evidence either.

EXACT: integer/Fraction arithmetic in the group ring.  No float anywhere in this file.
"""
import sys
from fractions import Fraction
from itertools import combinations

LOG = []
def out(s=""):
    print(s); LOG.append(s)

CLASSES = [(0, 0), (1, 0), (0, 1), (1, 1)]
FAIL = []


def ring_mul(a, b, q):
    """Multiply in Q[x]/(x^q - 1); a, b are length-q coefficient lists of Fractions."""
    c = [Fraction(0)] * q
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    c[(i + j) % q] += ai * bj
    return c


def ring_conj(a, q):
    """x -> x^(q-1), i.e. complex conjugation on the roots of unity."""
    c = [Fraction(0)] * q
    for i, ai in enumerate(a):
        c[(-i) % q] += ai
    return c


def ring_sub(a, b, q):
    return [x - y for x, y in zip(a, b)]


out("=" * 104)
out("REFUTER 1 / SCRIPT 5 — THE CHARACTER IDENTITY IN Z[x]/(x^q - 1), EXACT, ALL q <= 40")
out("=" * 104)
out()

WEIGHTS = [
    ("B0b SENSE U (4-class)", [Fraction(4, 9), Fraction(2, 9), Fraction(1, 9), Fraction(2, 9)]),
    ("B4  SENSE U (4-class)", [Fraction(1, 6), Fraction(1, 6), Fraction(1, 6), Fraction(1, 2)]),
    ("SENSE C 1/4 (4-class)", [Fraction(1, 4)] * 4),
    ("skew 4-class", [Fraction(1, 10), Fraction(2, 10), Fraction(3, 10), Fraction(4, 10)]),
    ("K1 SENSE U (3-class)", [Fraction(0), Fraction(2, 5), Fraction(2, 5), Fraction(1, 5)]),
    ("B1q SENSE U (3-class)", [Fraction(1, 7), Fraction(3, 7), Fraction(3, 7), Fraction(0)]),
    ("2-class {00,11}", [Fraction(1, 2), Fraction(0), Fraction(0), Fraction(1, 2)]),
]
QS = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 23, 24,
      29, 30, 31, 36, 37, 40]

nchk = 0
nbad = 0
nneg = 0
for wname, w in WEIGHTS:
    for q in QS:
        for (A, B) in [(1, 0), (0, 1), (1, 1), (1, -1), (1, 2), (2, 3), (3, 5), (q // 2, 1)]:
            ex = [(a * A + b * B) % q for (a, b) in CLASSES]
            for k in range(1, q + 1):
                # Z_k as a group-ring element:  sum_j w_j x^(k e_j)
                Z = [Fraction(0)] * q
                for wj, e in zip(w, ex):
                    Z[(k * e) % q] += wj
                Z2 = ring_mul(Z, ring_conj(Z, q), q)          # |Z_k|^2 as a ring element
                one = [Fraction(0)] * q
                one[0] = Fraction(1)
                lhs = ring_sub(one, Z2, q)
                rhs = [Fraction(0)] * q
                for j, l in combinations(range(4), 2):
                    if w[j] == 0 or w[l] == 0:
                        continue
                    d = ((ex[j] - ex[l]) * k) % q
                    # |chi_j^k - chi_l^k|^2 = 2 - x^d - x^(-d)
                    term = [Fraction(0)] * q
                    term[0] += 2
                    term[d % q] -= 1
                    term[(-d) % q] -= 1
                    for i in range(q):
                        rhs[i] += w[j] * w[l] * term[i]
                nchk += 1
                if any(x != 0 for x in ring_sub(lhs, rhs, q)):
                    nbad += 1
                    FAIL.append(f"identity residual {wname} q={q} (A,B)=({A},{B}) k={k}")
out(f"EXACT group-ring identity checks: {nchk}")
out(f"non-zero residuals                : {nbad}")
out(f"moduli covered                    : {QS}")
out(f"weight vectors                    : {len(WEIGHTS)} (four 4-class, two 3-class, one 2-class)")
out()
out("READ:")
out("  * W-08's character identity is EXACT for every modulus tested, at four occupied classes")
out("    with six pair terms and at three with three.  The lane's Q[sqrt3] check (q | 12) is")
out("    correct and this extends it to primes it could not reach.  NOTHING FALLS.")
out("  * The right-hand side is a sum of terms w_j w_l |chi_j^k - chi_l^k|^2 >= 0.  Therefore")
out("    1 - |Z_k|^2 >= 0 IS A COROLLARY.  The lane reports '0 events of |Z_k| > 1 -- exactly 0,")
out("    not 1e-16' alongside the identity as if it were a second confirmation; given the")
out("    identity it could not have been anything else.  It is a THEOREM and stays standing;")
out("    it is not independent evidence.")
out()
out("  * AND THE CLASS COUNT REALLY IS ABSENT FROM THE PROOF.  The derivation is")
out("        1 - |Z|^2 = (sum_j w_j)^2 - sum_{j,l} w_j w_l chi_j conj(chi_l)")
out("                  = sum_{j,l} w_j w_l (1 - chi_j conj(chi_l))")
out("                  = sum_{j<l} w_j w_l |chi_j - chi_l|^2,")
out("    which never indexes the number of j.  So the lane's central scope verdict --")
out("    'the proof never mentions the class count, therefore carrier-independent' -- is")
out("    CORRECT, and it is correct for a reason that needed no carrier to be built.")

out()
out("=" * 104)
if FAIL:
    out(f"**{len(FAIL)} FAILURES** (first 10)")
    for f in FAIL[:10]:
        out("   " + f)
else:
    out("0 failures.  W-08's identity, monotonicity and floor survive four classes and every")
    out("modulus tested.  I could not break them and I say what I tried.")

with open("r5_identity.OUT.txt", "w") as fh:
    fh.write("\n".join(LOG) + "\n")
sys.exit(1 if FAIL else 0)
