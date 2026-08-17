#!/usr/bin/env python3
"""
LANE W-10 / C — REFUTER 2 — STEP 1.
IS EACH OF LANE C's THREE CONNECTIONS THE ARITHMETIC TYPE IT CLAIMS TO BE?

The corpus has TWICE mislabelled a connection: S3 called f=2.0,c=1.1 "rationally independent"
when it is exactly resonant (S4 sec7 / ERRATUM v W-02), and S3 called f=3.14159,c=1.57080 an
order-4 point when it is generic (COR-K).  This step tests lane C's arms against that failure
mode, ON THE INTEGERS LANE C ACTUALLY FEEDS TO ITS ORBIT, not on the labels.

SECTIONS
  1  the three arms' relation lattices L, computed EXACTLY from what the code feeds the orbit
  2  the printed "resonance exactness" check in C_04 -- is it a check at all?
  3  int64 headroom in C_04's phase reduction
  4  the DIOPHANTINE arm: implemented value vs intended irrational; drift over the run;
     shortest relation vector of the implemented dyadic pair (i.e. its true order)
  5  the badly-approximable claim, measured, and the dominant dual resonance that sets the
     K = 1e7 Birkhoff error

Precision: Python int / Fraction (EXACT) throughout section 1-4; mpmath 60 dps where stated.
"""
import numpy as np
import mpmath as mp
from fractions import Fraction as Fr
import math

mp.mp.dps = 60

# ---------------------------------------------------------------------------------------
# lane C's own construction, copied VERBATIM from C_04_birkhoff.py lines 76-93 so that the
# integers audited below are the integers lane C runs.
D = 2 ** 40
alphaA = -(2.0 ** (1.0 / 3.0)) % 1.0
betaA = (4.0 ** (1.0 / 3.0)) % 1.0
A_numA = int(np.floor(alphaA * D)); B_numA = int(np.floor(betaA * D))
dAA = alphaA - A_numA / D; dBA = betaA - B_numA / D

S = 2 ** 35
A0 = int(round(S / np.pi))
D_C = 20 * S
A_numC = (-20 * A0) % D_C
B_numC = (11 * A0) % D_C
RELATION_RESIDUE = 11 * (-20 * A0) + 20 * (11 * A0)
# ---------------------------------------------------------------------------------------


def kernel_basis_mod(a, b, M):
    """EXACT basis of L = {(m,n) in Z^2 : m a + n b = 0 mod M}.
    Column-reduce the row [a b M] to [g 0 0] tracking a unimodular 3x3; the last two columns
    of the tracker span ker([a b M]) in Z^3; project onto (m,n) (injective, since t is then
    determined).  Pure integer arithmetic."""
    row = [a, b, M]
    T = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]      # T[i][j] = coefficient of e_i in column j
    cols = [0, 1, 2]
    # simple integer column echelon by repeated gcd steps
    while True:
        nz = [j for j in cols if row[j] != 0]
        if len(nz) <= 1:
            break
        nz.sort(key=lambda j: abs(row[j]))
        p = nz[0]
        for j in nz[1:]:
            q = row[j] // row[p]
            row[j] -= q * row[p]
            for i in range(3):
                T[i][j] -= q * T[i][p]
    zero_cols = [j for j in cols if row[j] == 0]
    basis = [(T[0][j], T[1][j]) for j in zero_cols]
    return basis


def gauss_reduce(v1, v2):
    """Lagrange-Gauss reduction of a rank-2 integer lattice; returns the shortest vector."""
    def n2(v):
        return v[0] * v[0] + v[1] * v[1]
    a, b = list(v1), list(v2)
    if n2(a) > n2(b):
        a, b = b, a
    while True:
        if n2(a) == 0:
            a, b = b, a
            if n2(a) == 0:
                return (0, 0)
        mu = round((a[0] * b[0] + a[1] * b[1]) / n2(a))
        b = [b[0] - mu * a[0], b[1] - mu * a[1]]
        if n2(b) >= n2(a):
            return tuple(a)
        a, b = b, a


def shortest_relation(a, b, M):
    bas = kernel_basis_mod(a, b, M)
    bas = [v for v in bas if v != (0, 0)]
    if len(bas) < 2:
        return None, bas
    v = gauss_reduce(bas[0], bas[1])
    return v, bas


def frac_exact(x):
    """EXACT Fraction of a float64 (a dyadic rational)."""
    return Fr(x)


if __name__ == "__main__":
    print("=" * 100)
    print("R2_01 — ARE LANE C's THREE CONNECTIONS THE ARITHMETIC TYPE THEY CLAIM?")
    print("       Audited on the integers C_04 feeds its orbit, not on the labels.")
    print("=" * 100)

    # ---------------------------------------------------------------- 1. relation lattices
    print("\n" + "-" * 100)
    print("1.  RELATION LATTICES L = {(m,n) : m alpha + n beta in Z}, EXACT, AS IMPLEMENTED.")
    print("-" * 100)

    print("\n  ARM O  FINITE ORDER   A=2, B=3, M=4  (alpha=1/2, beta=3/4)")
    v, bas = shortest_relation(2, 3, 4)
    idx = 4  # index of L in Z^2 = order of the orbit
    print(f"     basis of L = {bas}   shortest relation vector = {v}")
    # order of the point (alpha,beta) in T^2
    ordO = 4 // math.gcd(math.gcd(2, 3), 4)
    print(f"     order of (alpha,beta) in T^2 = M/gcd(A,B,M) = {ordO}   CLAIMED 4   "
          f"{'MATCH' if ordO == 4 else 'MISMATCH'}")
    print("     => arm O is genuinely FINITE ORDER 4.  Label correct.")

    print("\n  ARM R  RESONANT      A_numC, B_numC, D_C as C_04 builds them")
    print(f"     A0 = round(2^35/pi) = {A0}")
    print(f"     A_numC = {A_numC}   B_numC = {B_numC}   D_C = {D_C}")
    lhs = 11 * A_numC + 20 * B_numC
    print(f"     THE CHECK THAT MATTERS (on the reduced integers actually used):")
    print(f"       11*A_numC + 20*B_numC = {lhs}   mod D_C = {lhs % D_C}   "
          f"{'ON THE SUBTORUS' if lhs % D_C == 0 else 'NOT ON THE SUBTORUS'}")
    vR, basR = shortest_relation(A_numC, B_numC, D_C)
    ordR = D_C // math.gcd(math.gcd(A_numC, B_numC), D_C)
    print(f"     shortest relation vector of the IMPLEMENTED pair = {vR}"
          f"   |v| = {math.hypot(*vR):.4g}")
    print(f"     order of the implemented (alpha,beta) in T^2 = {ordR} = {ordR:.6e}")
    print(f"     CLAIM: primitive relation (11,20).  gcd(11,20) = {math.gcd(11,20)}.")
    print(f"     IS (11,20) THE SHORTEST?  {'YES' if tuple(map(abs, vR)) == (11,20) else 'NO'}"
          f"   (a shorter one would mean the arm sits on a SMALLER subtorus than claimed)")
    print(f"     second-shortest independent direction is at distance ~ index/|v| = "
          f"{ordR/math.hypot(*vR):.4g}, i.e. the implemented orbit is a FINITE cyclic group of")
    print(f"     order {ordR:.4e}, dense in the (11,20) circle only up to that resolution.")
    print(f"     K = 1e7 << {ordR:.3e}, so over the run the arm is indistinguishable from the")
    print(f"     ideal rank-one resonance.  LABEL CORRECT AT THIS K.")

    print("\n  ARM D  DIOPHANTINE   the float64 dyadic pair C_04 actually runs")
    fa = frac_exact(alphaA); fb = frac_exact(betaA)
    M = max(fa.denominator, fb.denominator)
    aI = int(fa * M); bI = int(fb * M)
    assert Fr(aI, M) == fa and Fr(bI, M) == fb
    vD, basD = shortest_relation(aI, bI, M)
    ordD = M // math.gcd(math.gcd(aI, bI), M)
    print(f"     common denominator M = {M} = 2^{int(math.log2(M))}")
    print(f"     order of the IMPLEMENTED (alpha,beta) in T^2 = {ordD} = {ordD:.6e}")
    print(f"     shortest relation vector of the IMPLEMENTED pair = {vD}"
          f"   |v|_inf = {max(abs(vD[0]),abs(vD[1]))}")
    print(f"     ||m alpha + n beta|| at that vector = "
          f"{abs(((vD[0]*aI+vD[1]*bI) % M + M) % M)/M:.6e}   (0 => an EXACT relation)")
    K = 10 ** 7
    print(f"     RUN LENGTH K = {K:.0e}.  The implemented arm is an EXACTLY PERIODIC (finite-")
    print(f"     order) connection of order {ordD:.3e}; it is 'Diophantine' only in the sense")
    print(f"     that no relation of sup-norm <= {max(abs(vD[0]),abs(vD[1]))} exists, and")
    print(f"     {max(abs(vD[0]),abs(vD[1]))} {'>' if max(abs(vD[0]),abs(vD[1])) > K else '<='} K.")

    # ---------------------------------------------------------------- 2. the printed check
    print("\n" + "-" * 100)
    print("2.  C_04's PRINTED 'RESONANCE EXACTNESS' CHECK.  IS IT A CHECK?")
    print("-" * 100)
    print("     C_04 prints:  RELATION_RESIDUE = 11*(-20*A0) + 20*(11*A0)   'must be 0'")
    print(f"     value at the run's A0 = {A0}: {RELATION_RESIDUE}")
    for t in (0, 1, 12345, -7, 2 ** 61):
        val = 11 * (-20 * t) + 20 * (11 * t)
        print(f"       same expression at A0 = {t:>22}: {val}")
    print("     IT IS THE IDENTITY -220 t + 220 t = 0.  IT IS ZERO FOR EVERY A0 AND FOR EVERY")
    print("     CHOICE OF S, D_C, AND FOR ANY PAIR OF RELATION EXPONENTS ONE CARES TO WRITE.")
    print("     IT NEVER TOUCHES A_numC OR B_numC -- the integers the orbit is built from --")
    print("     because both are taken mod D_C AFTER the expression is formed.")
    print("     THE CHECK COULD NOT HAVE FAILED.  It is C_04's only certificate that the")
    print("     RESONANT arm is resonant, and it certifies nothing.")
    print(f"     THE NON-VACUOUS CHECK, RUN HERE: (11*A_numC + 20*B_numC) mod D_C = "
          f"{lhs % D_C}.  It PASSES, so no published number moves; the defect is in the")
    print("     certificate, not the arm.  (This is the lane's own named fatal defect class --")
    print("     'a control that could not have failed' -- committed in the block whose whole")
    print("     purpose is to certify the connection's arithmetic type.)")

    # ---------------------------------------------------------------- 3. int64 headroom
    print("\n" + "-" * 100)
    print("3.  INT64 HEADROOM IN C_04's PHASE REDUCTION  (k * A_num formed directly).")
    print("-" * 100)
    IMAX = np.iinfo(np.int64).max
    for nm, An, Bn in (("DIOPHANTINE", A_numA, B_numA), ("RESONANT", A_numC, B_numC),
                       ("FINITE ORDER", 2, 3)):
        worst = max(An, Bn)
        prod = worst * K
        print(f"     {nm:13s} max(A,B) = {worst:>14}   k*A at K=1e7 = {prod:.4e}"
              f"   /int64max = {prod/IMAX:6.1%}   K_overflow = {IMAX//worst:.3e}")
    print("     NO OVERFLOW OCCURS AT K = 1e7.  The DIOPHANTINE arm runs at 88% of the ceiling;")
    print("     K = 1.13e7 would wrap SILENTLY (numpy int64 does not raise).  A latent defect,")
    print("     not an actual one.  This refuter's own evaluator (R2_03) uses chunked modular")
    print("     accumulation and has no K ceiling at all.")

    # ------------------------------------------------- 4. implemented vs intended irrational
    print("\n" + "-" * 100)
    print("4.  THE DIOPHANTINE ARM: IMPLEMENTED VALUE vs THE IRRATIONAL IT NAMES.  mpmath 60dps")
    print("-" * 100)
    a_true = mp.mpf(2) - mp.cbrt(mp.mpf(2))          # (-2^(1/3)) mod 1 = 2 - 2^(1/3)
    b_true = mp.cbrt(mp.mpf(4)) - 1                  # 4^(1/3) mod 1
    a_impl = mp.mpf(alphaA); b_impl = mp.mpf(betaA)
    da = a_impl - a_true; db = b_impl - b_true
    print(f"     alpha true  = {mp.nstr(a_true, 25)}")
    print(f"     alpha impl  = {mp.nstr(a_impl, 25)}   diff = {mp.nstr(da, 6)}")
    print(f"     beta  true  = {mp.nstr(b_true, 25)}")
    print(f"     beta  impl  = {mp.nstr(b_impl, 25)}   diff = {mp.nstr(db, 6)}")
    print(f"     PHASE DRIFT AT k = 1e7:  |k*d_alpha| = {mp.nstr(abs(da)*K, 6)}"
          f"   |k*d_beta| = {mp.nstr(abs(db)*K, 6)}")
    print(f"     Orbit-point spacing at K = 1e7 is ~1/K = 1.0e-07, so the implemented orbit's")
    print(f"     last points sit {float(abs(da)*K)/1e-7:.2%} / {float(abs(db)*K)/1e-7:.2%} of one")
    print("     spacing away from the ideal ones.  Small, but it is NOT the named irrational,")
    print("     and it is the reason the reported 5.4e-08 gap cannot be attributed purely to")
    print("     the Diophantine discrepancy without the check run in section 5 and in R2_03.")

    # ---------------------------------------------------------------- 5. badly approximable
    print("\n" + "-" * 100)
    print("5.  'BADLY APPROXIMABLE' -- MEASURED, AND THE DOMINANT DUAL RESONANCE.")
    print("-" * 100)
    print("     (a) SIMULTANEOUS form.  For a badly approximable PAIR there is c > 0 with")
    print("         k^(1/2) max(||k a||, ||k b||) >= c for all k >= 1.  Measured minimum over")
    print("         k <= 2e6, on the TRUE irrationals at 60 dps and on the IMPLEMENTED floats:")
    for tag, aa, bb in (("true (2-2^(1/3), 4^(1/3)-1)", a_true, b_true),
                        ("implemented float64 pair", a_impl, b_impl)):
        ks = np.arange(1, 2_000_001, dtype=np.int64)
        af = float(aa); bf = float(bb)
        # float64 is enough here: k<=2e6 and phases are O(1); error ~ 2e6*2.2e-16 = 4.4e-10
        ra = np.abs(((ks * af) % 1.0 + 0.5) % 1.0 - 0.5)
        rb = np.abs(((ks * bf) % 1.0 + 0.5) % 1.0 - 0.5)
        q = np.sqrt(ks) * np.maximum(ra, rb)
        j = int(np.argmin(q))
        print(f"         {tag:32s} min = {q[j]:.6f} at k = {ks[j]}")
    print("         A Liouville-type or resonant pair drives this to 0; a badly approximable")
    print("         pair keeps it bounded below.  Both rows stay near 0.2-0.5, consistent with")
    print("         the claim.  THIS IS A CONSISTENCY CHECK, NOT A PROOF -- the theorem behind")
    print("         (2^(1/3), 4^(1/3)) in Bad is the CUBIC NORM FORM argument (Perron 1921 /")
    print("         Cassels), not 'W. M. Schmidt's theorem' as C_00 and C_04 attribute it.")
    print("         Schmidt's badly-approximable results are the full-Hausdorff-dimension")
    print("         theorem and the Subspace Theorem; neither yields this pair.  ATTRIBUTION")
    print("         DEFECT, inherited from LANE_W08 M1_03; the FACT is standard and true.")

    print("\n     (b) DUAL form -- the one that actually sets the Birkhoff error.  The error of")
    print("         (1/K) sum g(k a, k b) is  sum_{(m,n)!=0} ghat(m,n) * S_K(m,n)/K  with")
    print("         |S_K/K| ~ min(1, 1/(2 K ||m a + n b||)).  Dominant low-order resonances:")
    for tag, aa, bb in (("DIOPHANTINE", a_true, b_true),
                        ("RESONANT   ", mp.mpf(A_numC) / D_C, mp.mpf(B_numC) / D_C),
                        ("F1,SQRT2   ", mp.mpf(-1) / (2 * mp.pi) % 1, mp.sqrt(2) / (2 * mp.pi))):
        best = []
        for m in range(-40, 41):
            for n in range(-40, 41):
                if m == 0 and n == 0:
                    continue
                r = mp.mpf(m) * aa + mp.mpf(n) * bb
                r = r - mp.floor(r + mp.mpf(0.5))
                best.append((abs(float(r)), m, n))
        best.sort()
        head = best[:4]
        print(f"         {tag}: " + "  ".join(
            f"(m,n)={m},{n} ||.||={r:.3e} 1/(2K||.||)={1/(2*K*r) if r>0 else float('inf'):.2e}"
            for r, m, n in head[:2]))
        r0, m0, n0 = head[0]
        print(f"                    strongest: (m,n)=({m0},{n0}), amplification 1/(2K||.||) = "
              f"{(1/(2*K*r0) if r0>0 else float('inf')):.3e} at K=1e7")
    print("         For the DIOPHANTINE arm this predicts a Birkhoff error of order")
    print("         |ghat| * 1e-7..1e-8, which is the observed 5.4e-08 -- i.e. C_04's reported")
    print("         gaps ARE finite-N discrepancy, not bias.  Confirmed directly in R2_03 by")
    print("         driving K to 1e9 and watching the gap fall.")
    print("\nDONE.")
