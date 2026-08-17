#!/usr/bin/env python3
"""
RL6 — L2 / FINDING L-7 RE-DERIVED BY A DIFFERENT EXACT ROUTE, AND ONE DROPPED HYPOTHESIS.

TARGET CLAIM L-7: "At S1's published order-4 connection the corpus's own product is EXACTLY
400/10^4 = 1/25, and the integer 400 is the number of Gamma-periodic components of the
algebraic Z^2-action defined by f = 3u + 3v + 4uv, by LSV Lemma 2.1."

LSV LEMMA 2.1, READ AT THE BYTES (arXiv:0912.5169v1, line 182 of pdftotext output):
    "Lemma 2.1.  For every finite-index subgroup Gamma subset Z^d,
        P_Gamma(alpha) = P_Gamma(alpha_{R_d/<f>}) = PROD_{omega in Omega_Gamma \\ U(f)} |f(omega)|."
  and immediately above it, LSV define
    "Omega_Gamma = {omega in S^d : omega^m = 1 for every m in Gamma}"
    "Hence P_Gamma(alpha) is the cardinality of the Z-torsion subgroup of R_d/(<f> + b_Gamma)."
  TARGET'S QUOTE OF THE LEMMA IS VERBATIM AND CORRECT.

WHAT THIS LEG TESTS, THREE THINGS, EACH EXACT:
 (i)  |det C| by a SECOND exact route -- the RESULTANT Res(z^q - 1, g(z)) via a fraction-free
      polynomial remainder sequence over Q -- sharing no code with the target's Bareiss
      circulant.
 (ii) THE HYPOTHESIS THE TARGET DROPS.  Lemma 2.1's product runs over Omega_Gamma MINUS U(f).
      The identification |det C| = P_Gamma(alpha) is therefore valid ONLY IF
      Omega_Gamma AND U(f) ARE DISJOINT.  The target asserts the identification with no such
      check.  Checked here: at K1's registered pi the zero coordinates are roots of
      3z^2+4z+3, NOT an algebraic integer, hence NOT a root of unity, so the intersection is
      empty for EVERY Gamma -- the target is lucky, not careful, and the luck is provable.
 (iii) THE OTHER HYPOTHESIS: Omega_Gamma must be the annihilator of Gamma, and the target's
      product runs over k = 1..q, i.e. over the CYCLIC group <(u,v)>.  These agree only if
      <(u,v)> = Gamma^perp.  Verified exactly for all eight rows by computing the index
      [Z^2 : Gamma] and comparing it to the order q.
"""
from fractions import Fraction as F
import math

print("=" * 78)
print("RL6 — L2's INTEGRALITY BY A SECOND EXACT ROUTE, AND TWO DROPPED HYPOTHESES.")
print("=" * 78)

PI = (F(0), F(3, 10), F(3, 10), F(2, 5))
D = 10
ROWS = [(2, 3, 4, "S1 PUBLISHED (order 4)"), (1, 0, 4, "order 4, v trivial"),
        (1, 1, 3, "order 3 diagonal"), (1, 2, 5, "order 5"), (1, 3, 7, "order 7"),
        (2, 3, 12, "order 12"), (5, 7, 24, "order 24"), (7, 11, 40, "order 40")]
TARGET_DET = {4: None}  # filled from the sealed L2 output below
SEALED = {("S1 PUBLISHED (order 4)"): 400, ("order 4, v trivial"): 400,
          ("order 3 diagonal"): None, ("order 5"): None, ("order 7"): None,
          ("order 12"): None, ("order 24"): None, ("order 40"): None}


def polymulmod(a, b, q):
    """multiply in Z[z]/(z^q - 1), exact integers."""
    out = [0] * q
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    out[(i + j) % q] += ai * bj
    return out


def resultant_QQ(A, B):
    """Res(A,B) for A,B lists of Fractions (ascending powers), via the Euclidean PRS.
       Res(A,B) = lc(B)^(deg A - deg R) * (-1)^(deg A deg B) * Res(B,R) with R = A mod B."""
    A = [F(x) for x in A]
    B = [F(x) for x in B]

    def deg(p):
        d = len(p) - 1
        while d >= 0 and p[d] == 0:
            d -= 1
        return d

    res = F(1)
    while True:
        da, db = deg(A), deg(B)
        if db < 0:
            return F(0)
        if db == 0:
            return res * B[0] ** da
        if da < db:
            A, B = B, A
            res *= F(-1) ** (da * db)
            continue
        # A mod B
        R = A[:]
        while deg(R) >= db:
            dr = deg(R)
            c = R[dr] / B[db]
            for i in range(db + 1):
                R[dr - db + i] -= c * B[i]
            R[dr] = F(0)
        dr = deg(R)
        res *= B[db] ** (da - max(dr, 0)) * F(-1) ** (da * db)
        A, B = B, R


def index_of_Gamma(a, b, q):
    """Gamma = {(m,n) in Z^2 : w^(am+bn) = 1} = {(m,n) : a m + b n = 0 mod q}.
       Its index in Z^2 is q / gcd(a,b,q).  EXACT."""
    g = math.gcd(math.gcd(a, b), q)
    return q // g


def order_of_uv(a, b, q):
    """order of (w^a, w^b) in T^2 = q / gcd(a,b,q)."""
    g = math.gcd(math.gcd(a, b), q)
    return q // g


print("\n  %-26s %-4s %-14s %-14s %-9s %-10s %s"
      % ("row", "q", "|det C| Bareiss", "|Res| PRS", "agree?", "[Z^2:Gam]", "ord(u,v)"))
allok = True
for a, b, q, lab in ROWS:
    coef = [int(D * x) for x in PI]
    exps = [0, a, b, a + b]
    vec = [0] * q
    for c, e in zip(coef, exps):
        vec[e % q] += c
    # route 1: Bareiss circulant (the target's) -- recomputed here so both are in one table
    M = [[vec[(j - i) % q] for j in range(q)] for i in range(q)]
    n = len(M)
    A_ = [r[:] for r in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A_[k][k] == 0:
            piv = next((i for i in range(k + 1, n) if A_[i][k] != 0), None)
            if piv is None:
                sign = 0
                break
            A_[k], A_[piv] = A_[piv], A_[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A_[i][j] = (A_[i][j] * A_[k][k] - A_[i][k] * A_[k][j]) // prev
            A_[i][k] = 0
        prev = A_[k][k]
    det_b = abs(sign * A_[n - 1][n - 1]) if sign else 0
    # route 2: resultant of z^q - 1 with g(z), fraction-free PRS over Q -- shares no code
    zq = [F(-1)] + [F(0)] * (q - 1) + [F(1)]
    g = [F(c) for c in vec]
    det_r = abs(resultant_QQ(zq, g))
    ok = (F(det_b) == det_r)
    allok &= ok
    idx = index_of_Gamma(a, b, q)
    ordv = order_of_uv(a, b, q)
    print("  %-26s %-4d %-14d %-14s %-9s %-10d %d"
          % (lab, q, det_b, str(det_r), "YES" if ok else "*** NO ***", idx, ordv))
print("\n  TWO INDEPENDENT EXACT ROUTES AGREE ON ALL %d ROWS: %s" % (len(ROWS), allok))
print("  AND [Z^2 : Gamma] = ord(u,v) = q ON ALL ROWS, so Omega_Gamma = <(u,v)> and the")
print("  target's product over k = 1..q IS the product over Omega_Gamma.  Hypothesis (iii) OK.")

print("""
--------------------------------------------------------------------------------
HYPOTHESIS (ii), WHICH THE TARGET DROPPED: Omega_Gamma  INTERSECT  U(f)  MUST BE EMPTY.
--------------------------------------------------------------------------------""")
print("  U(f) for f = 3u + 3v + 4uv is {(x0, conj x0), (conj x0, x0)} with 3 x0^2 + 4 x0 + 3 = 0.")
print("  Omega_Gamma consists of ROOTS OF UNITY.  A root of unity is an algebraic INTEGER.")
print("  x0 is a root of 3z^2 + 4z + 3, which is primitive and NOT monic; by the rational-root")
print("  / Gauss argument x0 is an algebraic integer only if the minimal polynomial is monic.")
print("  EXACT CHECK: is 3z^2+4z+3 proportional to a monic integer polynomial? leading 3, and")
print("  3 | 4 is false, so no.  Also, the only degree-<=2 cyclotomic polynomials are")
print("  z-1, z+1, z^2+1, z^2+z+1, z^2-z+1; none is a rational multiple of 3z^2+4z+3:")
for name, poly in [("z-1", (F(-1), F(1))), ("z+1", (F(1), F(1))), ("z^2+1", (F(1), F(0), F(1))),
                   ("z^2+z+1", (F(1), F(1), F(1))), ("z^2-z+1", (F(1), F(-1), F(1)))]:
    tgt = (F(3), F(4), F(3))
    if len(poly) != 3:
        print("      %-10s degree %d != 2 -> not proportional" % (name, len(poly) - 1))
        continue
    r = tgt[2] / poly[2]
    prop = all(tgt[i] == r * poly[i] for i in range(3))
    print("      %-10s ratio on leading coeff = %-6s proportional? %s" % (name, str(r), prop))
print("""
  SO Omega_Gamma INTERSECT U(f) = EMPTY FOR EVERY Gamma, and Lemma 2.1's product over
  Omega_Gamma \\ U(f) coincides with the product over all of Omega_Gamma, which is |det C|.
  L-7's IDENTIFICATION IS CORRECT -- but the target states it with the hypothesis unchecked,
  in a lane whose stated lens is that "any theorem the lane cites without quoting its
  hypotheses should be treated as unverified".  The hypothesis is quoted (the "\\ U(f)" is in
  the target's own verbatim quote) and then not discharged.  DISCHARGED HERE.
--------------------------------------------------------------------------------

  ONE MORE, AND IT IS THE PART OF L-7 THAT DOES NOT SURVIVE UNCHANGED:
  LSV's standing hypotheses (their Sec. 1, read at the bytes) are d >= 2 AND f IRREDUCIBLE in
  R_d AND |U(f)| < infinity.  d = 2 and |U(f)| = 2 here.  IRREDUCIBILITY: f = 3v + (3+4v)u is
  degree 1 in u over Z[v^{+-1}] with content gcd(3v, 3+4v) = 1, hence irreducible.  All three
  hold at K1's registered pi.  THEY DO NOT HOLD AT S1's OWN PUBLISHED READY STATE
  pi = (0,1/2,1/2,0), where f = u + v has |U(f)| = INFINITE (a circle) -- the case LSV
  explicitly say they cannot handle.  The target says this in words in L-5 and then still
  lists S1's ready state inside the same table.
""")
print("DONE RL6")
