#!/usr/bin/env python3
"""
M1_04 — THE DURABILITY CRITERION, STATED EXACTLY, AND ITS HYPOTHESES.

PART 1.  THE PRODUCT LEMMA, AND THE HYPOTHESIS THE BRIEF'S FORM IS MISSING.
  Claimed:  |Omega_N| = prod_{k<=N} |Z_k| -> 0  <=>  sum_k (1-|Z_k|) = infinity,
            "hypotheses: 0 <= |Z_k| <= 1".
  0 <= a_k <= 1 IS NOT SUFFICIENT.  Take a_1 = 0, a_k = 1 for k >= 2: the product is 0 for
  every N >= 1 while sum (1-a_k) = 1 < infinity.  The correct lemma is

      LEMMA.  Let a_k in [0,1].  Then prod_{k<=N} a_k -> 0  iff  ( some a_k = 0 )  or
              ( sum_k (1-a_k) = infinity ).  Equivalently, ASSUMING a_k > 0 for all k,
              prod a_k -> 0 iff sum (1-a_k) = infinity.
      PROOF.  If some a_k = 0 the product is eventually 0.  Otherwise all a_k in (0,1] and
      log prod = sum log a_k.  For a in (0,1]:  -log a >= 1-a, so sum(1-a_k) = infinity
      forces sum(-log a_k) = infinity and the product -> 0.  Conversely if sum(1-a_k) <
      infinity then a_k -> 1, so a_k >= 1/2 eventually, and on [1/2,1] one has
      -log a <= (1-a)/a <= 2(1-a); hence sum(-log a_k) < infinity and the product tends to a
      strictly positive limit.  []
  THE MISSING HYPOTHESIS IS NOT ACADEMIC ON THIS CARRIER: S1's published connection with
  S1's published ready state has Z_k = 0 for every odd k (M1_03 exhibit D).

PART 2.  ON K1 THE EXCLUDED BRANCH CANNOT DECIDE ANYTHING, AND THE CRITERION IS EXACT.
  Because Z_k = P(u^k, v^k) with P a 3-term character sum with NON-NEGATIVE coefficients
  summing to 1, the sequence 1-|Z_k| is the orbit sampling of a CONTINUOUS function on the
  closed subgroup H = closure{ (u^k,v^k) : k in Z } <= T^2.  Weyl on the monothetic group H
  gives  (1/N) sum_{k<=N} (1-|Z_k|)  ->  int_H (1 - |P|) dHaar_H  =: c(H,pi) >= 0,
  WITH NO DIOPHANTINE HYPOTHESIS (1-|P| is continuous; no singularity to dodge).  And
      c = 0  <=>  |P| = 1 identically on H
            <=>  all characters in the support of pi agree on H  (strict triangle inequality)
            <=>  G := < chi_a/chi_b : a,b in supp(pi) > = {1}.
  Hence on K1:

      |Omega_N| -> 0   <=>   G != {1},
      and when it holds the divergence of sum (1-|Z_k|) is LINEAR with density c > 0,
      and |Omega_N| -> 0 EXPONENTIALLY with rate lambda_H = int_H log|P| < 0.
      When G = {1}, |Z_k| = 1 for every k and |Omega_N| = 1 for every N.
  There is no intermediate case: the sum cannot converge with infinitely many terms non-zero.
  G is exactly W-01/W-02's corrected criterion; here it drops out as a statement about the
  RATE, not about a binary firing.

PART 3.  lambda_H VERSUS m(P).  lambda_H = m(P) iff H = T^2 iff u,v have no multiplicative
  relation.  m(P) < 0 does NOT imply lambda < 0: exhibited below on a NON-TRIVIAL connection.

Precision: float64.  The four-class characters are exact where they are roots of unity.
"""
import numpy as np
from M1_02_mahler_machinery import m_R1

print("=" * 78)
print("M1_04 — THE DURABILITY CRITERION ON K1")
print("=" * 78)

# ---------------------------------------------------------------- Part 1 checks
a = np.linspace(1e-12, 1.0, 2000001)
print("\nPART 1.  The two inequalities the lemma runs on (checked on a 2e6 grid, float64):")
print("   min over (0,1] of ( -log a - (1-a) )        = %.3e   (must be >= 0)"
      % float(np.min(-np.log(a) - (1 - a))))
b = np.linspace(0.5, 1.0, 1000001)
print("   min over [1/2,1] of ( 2(1-b) + log b )      = %.3e   (must be >= 0)"
      % float(np.min(2 * (1 - b) + np.log(b))))
print("   COUNTEREXAMPLE to the lemma as stated with only 0<=a_k<=1:")
print("     a_1 = 0, a_k = 1 (k>=2):  prod = 0 for all N;  sum (1-a_k) = 1 < infinity.")

# ---------------------------------------------------------------- Part 2 machinery
def Zseq(u, v, p10, p01, p11, K):
    k = np.arange(1, K + 1)
    return np.abs(p10 * u ** k + p01 * v ** k + p11 * (u * v) ** k)

def G_trivial(u, v, p10, p01, p11, tol=1e-12):
    """G = {1} iff all characters present in the support agree at (u,v)."""
    ch = []
    if p10 > 0: ch.append(u)
    if p01 > 0: ch.append(v)
    if p11 > 0: ch.append(u * v)
    return all(abs(ch[0] - c) < tol for c in ch)

CASES = [
    ("full support, generic connection",      np.exp(-2j), np.exp(1.1j), 0.3, 0.3, 0.4),
    ("full support, S1 published (ord 4)",    -1 + 0j,     -1j,          0.3, 0.3, 0.4),
    ("full support, TRIVIAL connection",      1 + 0j,      1 + 0j,       0.3, 0.3, 0.4),
    ("S={10,01} and u = v  (W_F W_C = 1)",    np.exp(0.7j), np.exp(0.7j), 0.5, 0.5, 0.0),
    ("S={10,01}, u != v",                     np.exp(0.7j), np.exp(1.3j), 0.5, 0.5, 0.0),
    ("S={10,11} and W_C = 1",                 np.exp(0.7j), 1 + 0j,      0.5, 0.0, 0.5),
    ("S={01,11} and W_F = 1",                 1 + 0j,      np.exp(0.7j), 0.0, 0.5, 0.5),
    ("S={11} only  (the ROOT alone)",         np.exp(-2j), np.exp(1.1j), 0.0, 0.0, 1.0),
    ("S1 published conn + S1 published state", -1 + 0j,    -1j,          0.0, 0.5, 0.5),
]
K = 200000
print("\nPART 2.  For each case: G trivial?  Cesaro density of (1-|Z_k|) at N=2e5;")
print("         sum_{k<=N}(1-|Z_k|); min|Z_k|; log|Omega_N|/N; and m(P) for comparison.")
print("  %-42s %-6s %11s %11s %11s %11s %11s" %
      ("case", "G={1}", "density", "min|Z|", "logOmega/N", "m(P)", "#Z_k=0"))
for (name, u, v, p10, p01, p11) in CASES:
    az = Zseq(u, v, p10, p01, p11, K)
    dens = float(np.mean(1 - az))
    nz = int(np.sum(az == 0.0))
    with np.errstate(divide='ignore'):
        lo = float(np.mean(np.log(az)))
    print("  %-42s %-6s %11.3e %11.3e %11.4f %11.6f %11d" %
          (name, str(G_trivial(u, v, p10, p01, p11)), dens, float(az.min()),
           lo, m_R1(0.0, p10, p01, p11), nz))

print("\n  READ-OFF:")
print("   * Every G={1} row has density EXACTLY 0 and log|Omega_N|/N = 0: the record is written")
print("     and un-written; |Omega_N| = 1 for all N.  Three of those rows are NON-TRIVIAL")
print("     connections (W_F != 1 or W_C != 1) — S3-audit COR-B's 'four families that never")
print("     form', re-derived here as the G={1} locus of the four possible supports.")
print("   * Every G!={1} row has density bounded away from 0 -> sum (1-|Z_k|) diverges LINEARLY.")
print("   * ROW 4 IS THE HYPOTHESIS SEPARATION: u = v, non-trivial connection, m(P) = -0.693147")
print("     yet lambda = 0 and NOTHING DECAYS.  m(P) < 0 does not imply durability.")
print("   * The last row has Z_k = 0 on odd k: durable via the branch the brief's lemma excludes.")

# ---------------------------------------------------------------- density is bounded below
print("\nPART 2b.  Is the Cesaro density bounded away from 0 off the G={1} locus?  1000 random")
print("          (connection, ready state) draws, N = 20000, float64, seed 20260816+40.")
rng = np.random.default_rng(20260816 + 40)
dens_min, worst = np.inf, None
for _ in range(1000):
    f, c = rng.uniform(0, 2 * np.pi, 2)
    u, v = np.exp(-1j * f), np.exp(1j * c)
    w = rng.dirichlet([1, 1, 1])
    az = Zseq(u, v, w[0], w[1], w[2], 20000)
    d = float(np.mean(1 - az))
    if d < dens_min:
        dens_min, worst = d, (f, c, w)
print("   min density over 1000 draws = %.6e   at f=%.4f c=%.4f pi=(%.3f,%.3f,%.3f)"
      % (dens_min, worst[0], worst[1], *worst[2]))
print("   (it is small only when the draw is NEAR the G={1} locus; it is 0 only ON it.)")

# ---------------------------------------------------------------- Part 3: rate vs m(P)
print("\nPART 3.  lambda_H vs m(P): the three regimes, one line each (from M1_03).")
print("   rank L = 0 (Diophantine)   lambda = m(P)          = -0.767507880   [= the N1 claim]")
print("   rank L = 1 (11,20)         lambda = subtorus mean = -0.767014993   != m(P)")
print("   rank L = 2 (order 4)       lambda = 4-point mean  = -0.804718956   != m(P)")
print("   and on the G={1} locus     lambda = 0                              != m(P) < 0")
print("\n   So N1's identification lambda = m(P) is TRUE EXACTLY ON THE COMPLEMENT of a dense,")
print("   Haar-null set of connections, and NOT on either of the two connections this corpus")
print("   has actually published.")

# ---------------------------------------------------------------- exponential decay check
print("\nPART 4.  |Omega_N| itself, on the Diophantine connection (float64 underflows; we track")
print("         log|Omega_N| = sum log|Z_k|, which does not).")
alpha = -(2.0 ** (1.0 / 3.0)) % 1.0
beta = (4.0 ** (1.0 / 3.0)) % 1.0
D = 2 ** 39
A = int(np.floor(alpha * D)); B = int(np.floor(beta * D))
dA = alpha - A / D; dB = beta - B / D
k = np.arange(1, 1000001, dtype=np.int64)
fa = np.mod(((k * A) % D) / D + k * dA, 1.0)
fb = np.mod(((k * B) % D) / D + k * dB, 1.0)
x = np.exp(2j * np.pi * fa); y = np.exp(2j * np.pi * fb)
az = np.abs(0.3 * x + 0.3 * y + 0.4 * x * y)
cl = np.cumsum(np.log(az))
cs = np.cumsum(1 - az)
print("      N        log|Omega_N|      log|Omega_N|/N     sum(1-|Z_k|)   first N with |Omega|<1e-300")
first = int(np.argmax(cl < np.log(1e-300))) + 1
for N in (10, 100, 1000, 10000, 100000, 1000000):
    print("   %8d   %15.4f   %15.9f   %14.2f" % (N, cl[N - 1], cl[N - 1] / N, cs[N - 1]))
print("   float64 underflow of |Omega_N| (|Omega_N| < 1e-300) first at N = %d" % first)
print("\nDONE M1_04")
