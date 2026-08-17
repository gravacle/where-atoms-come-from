#!/usr/bin/env python3
"""
LANE W-10 / C — REFUTER 2 — STEP 5.  THE ISOLATION VERDICT, ASSEMBLED.

1  DIFF THE THREE ARMS OF C_04 EXPLICITLY (the ledger says the arms were diffed; here is the
   diff, on every field of the call, not just on (alpha,beta)).
2  THE RESOLUTION TABLE: effect size (R2_02, exact) against ruler (R2_03, measured) at each K.
3  A SPOT-CHECK OF C-9's q_1 BOUND, the one rigorous limit-failure in the lane.
4  THIS REFUTER'S OWN DEFECTS, recorded rather than patched.
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 60

# --- C_04's three arms, reconstructed exactly as C_04 builds them -------------------------
D = 2 ** 40
alphaA = -(2.0 ** (1.0 / 3.0)) % 1.0
betaA = (4.0 ** (1.0 / 3.0)) % 1.0
A_numA = int(np.floor(alphaA * D)); B_numA = int(np.floor(betaA * D))
dAA = alphaA - A_numA / D; dBA = betaA - B_numA / D
S = 2 ** 35; A0 = int(round(S / np.pi)); D_C = 20 * S
ARMS = [
    ("DIOPHANTINE", A_numA, B_numA, D, dAA, dBA, "m(P)", "Jensen trapezoid n = 2^20"),
    ("RESONANT", (-20 * A0) % D_C, (11 * A0) % D_C, D_C, 0.0, 0.0, "subtorus",
     "trapezoid on log|Q| n = 2^22"),
    ("FINITE ORDER", 2, 3, 4, 0.0, 0.0, "finite", "mpmath exact 4-point average"),
]

# --- effect sizes (R2_02, roots/Jensen, 50 dps) and rulers (R2_03, this refuter's runs) ---
EFFECT = {                       # |m(Q_(11,20)) - m(P)|, EXACT
    "B0b": mp.mpf('1.16812798517e-8'),
    "B4": mp.mpf('4.84421064413e-10'),
    "K1": mp.mpf('2.36576616316e-4'),
    "SENSE C": mp.mpf(0),
}
RULER = {                        # |avg - the limit it converges to|, at each K (R2_03 sec 1)
    "B0b": {10 ** 7: 2.068e-08, 10 ** 8: 3.526e-09, 3 * 10 ** 8: 1.393e-09, 10 ** 9: 2.502e-10},
    "B4": {10 ** 7: 7.267e-09, 10 ** 8: 5.270e-09, 3 * 10 ** 8: 1.836e-10, 10 ** 9: 8.108e-11},
    "K1": {10 ** 7: 7.904e-08, 10 ** 8: 2.775e-08, 3 * 10 ** 8: 1.148e-08, 10 ** 9: 1.063e-10},
    "SENSE C": {10 ** 7: 8.917e-07, 10 ** 8: 1.260e-07, 3 * 10 ** 8: 1.314e-07, 10 ** 9: 4.469e-08},
}

if __name__ == "__main__":
    print("=" * 110)
    print("R2_05 (1) — THE THREE ARMS OF C_04, DIFFED ON EVERY FIELD OF THE CALL")
    print("=" * 110)
    hdr = f"  {'field':16s} " + "".join(f"{a[0]:>26s}" for a in ARMS)
    print(hdr)
    fields = [("A_num", 1), ("B_num", 2), ("Dn", 3), ("dA", 4), ("dB", 5)]
    for nm, i in fields:
        print(f"  {nm:16s} " + "".join(f"{str(a[i]):>26s}" for a in ARMS))
    print(f"  {'alpha':16s} " + "".join(f"{a[1]/a[3]+a[4]:26.18f}" for a in ARMS))
    print(f"  {'beta':16s} " + "".join(f"{a[2]/a[3]+a[5]:26.18f}" for a in ARMS))
    print(f"  {'TARGET evaluator':16s} " + "".join(f"{a[7]:>30s}" for a in ARMS))
    print("\n  THE ORBIT SIDE IS ISOLATED: same orbit_average(), same K = 1e7, same CHUNK = 1e6,")
    print("  same checkpoints, same ready state, and the three (alpha,beta) are distinct.")
    print("  TWO THINGS MOVE WITH THE CONNECTION AND ARE NOT NAMED IN THE LEDGER:")
    print("    (i)  the phase REPRESENTATION -- arm 1 carries a nonzero float64 correction")
    print("         (dA,dB); arms 2 and 3 are exact rationals with dA = dB = 0.  FORCED by the")
    print("         arithmetic type, so not a defect, but it is a second moving coordinate.")
    print("         MEASURED IMMATERIAL: R2_03 re-ran arm 1 with exact 80-bit phases on the TRUE")
    print("         irrational and got -0.8109302706 at K = 1e7 against C_04's -0.810930271.")
    print("         THIS ATTACK FAILED.  Named because it was tried, not because it landed.")
    print("    (ii) the TARGET evaluator -- three different quadratures of three different")
    print("         accuracies (row above).  This one is NOT immaterial: see (2).")

    print("\n" + "=" * 110)
    print("R2_05 (2) — RESOLUTION.  CAN THE ARM SEE THE EFFECT IT IS CLAIMED TO SEE?")
    print("            effect = |m(Q_(11,20)) - m(P)|, EXACT (R2_02, roots + trapezoid agree")
    print("            to 5e-50).  ruler = |avg - limit| at that K (R2_03, exact 80-bit phases).")
    print("=" * 110)
    print(f"  {'case':10s} {'effect (exact)':>18s} " +
          "".join(f"{'K=%.0e' % k:>13s}" for k in [10**7, 10**8, 3*10**8, 10**9]))
    for c in ("B0b", "B4", "K1", "SENSE C"):
        line = f"  {c:10s} {mp.nstr(EFFECT[c], 6):>18s} "
        for k in (10 ** 7, 10 ** 8, 3 * 10 ** 8, 10 ** 9):
            r = RULER[c][k]
            snr = float(EFFECT[c]) / r if r > 0 else float('inf')
            line += f"{snr:13.2f}"
        print(line)
    print("  (entries are effect/ruler.  < 3 = the arm cannot tell the two candidate limits")
    print("   apart at that K.  0.00 for SENSE C means the two limits ARE THE SAME NUMBER.)")
    print("\n  AT LANE C's K = 1e7 THE RESONANT ARM RESOLVES ON EXACTLY ONE OF THE FOUR ROWS,")
    print("  AND IT IS K1 -- THE THREE-CLASS CONTROL THIS ROUND EXISTS TO MOVE OFF.")
    print("  On B0b it reads 0.56 and the measured value lands NEARER m(P) than nearer the")
    print("  subtorus (8.78e-09 vs 2.05e-08).  On B4 it reads 0.07.  On SENSE C it is 0.00")
    print("  BY EXACT IDENTITY: Q = (1+z^11)(1+z^20)/4, so m(Q) = -log 4 = m(P).")
    print("  THE LIMIT CLAIM IS NEVERTHELESS TRUE ON B0b AND B4 -- at K = 1e8 .. 1e9 the")
    print("  average leaves m(P) and converges on m(Q) (R2_03 sec 1).  The THEOREM survives;")
    print("  the MEASUREMENT that was offered as its evidence does not reach it.")

    print("\n" + "=" * 110)
    print("R2_05 (3) — SPOT-CHECK OF C-9's q_1 BOUND (the lane's one rigorous limit failure).")
    print("=" * 110)
    # alpha = gamma/2, gamma = 2^-3 + 2^-256 + ... ;  q = 8 -> 8 alpha = 1/2 + 2^-254 + tiny
    mp.mp.dps = 300          # eps = 2^-254 ~ 3.4e-77 is BELOW 60-dps resolution of 1/2 + eps
    eps = mp.mpf(2) ** -254
    u8 = mp.expjpi(2 * (mp.mpf('0.5') + eps))
    v = mp.expjpi(2 * (8 * (mp.cbrt(4) - 1)))
    Z8 = (1 + u8) * (1 + v) / 4
    bound = mp.log(abs(Z8)) / 8
    print(f"     |1 + u^8| = 2 sin(pi eps) = {mp.nstr(2*mp.sin(mp.pi*eps), 6)}"
          f"   (C-9 states 2.17e-76)")
    print(f"     |1 + v^8|                 = {mp.nstr(abs(1+v), 8)}")
    print(f"     log|Z_8| / 8              = {mp.nstr(bound, 8)}   (C-9 states <= -21.864)")
    print(f"     m(P) for SENSE C          = {mp.nstr(-mp.log(4), 8)}"
          f"   ratio = {mp.nstr(bound/(-mp.log(4)), 6)}   (C-9 states factor 15.8)")
    print("     CONFIRMED.  Note this bound uses beta = 4^(1/3), the SAME beta as the")
    print("     DIOPHANTINE arm; only alpha changes.  That is a clean one-variable move and")
    print("     it is the strongest thing in the lane.")

    print("\n" + "=" * 110)
    print("R2_05 (4) — THIS REFUTER'S OWN DEFECTS, RECORDED RATHER THAN PATCHED")
    print("=" * 110)
    print("""  R-1  R2_03 section 4 was written and RUN with the header 'f = 1.0, c = 3.0 -- a second
       published connection [S4:330].  Also generic.'  IT IS NOT GENERIC.  Both f and c are
       rational, so -3f + c = 0 and the primitive relation is (3,1): the connection is
       exactly resonant.  The run's flat errors (7.9e-02, 2.4e-02, 1.5e-02, not falling with
       K) are the correct behaviour of a subtorus limit, not a failure of convergence.  I
       caught this by RUNNING it, not by reading it.  The header is left in R2_03's output
       with this correction beside it rather than edited away, because the defect is the
       exact one this lens was commissioned to hunt -- and I committed it in the act of
       hunting it.  R2_04 then types every published connection by computation instead.
  R-2  R2_02 computed m(P) with an UNSPLIT 2^16 Jensen trapezoid.  That is O(n^-2) at a
       branch crossing and O(n^-1) at a log singularity, so R2_02's m(P) is wrong by 8.4e-11
       (K1), 1.7e-10 (B0b*) and 1.06e-05 (SENSE C).  Consequences: R2_02's 'effect size' row
       for SENSE C reads 1.06e-05 when the true effect is EXACTLY 0, and its K1 row is off in
       the eleventh place.  R2_04(B) recomputes both by kink-split adaptive quadrature and
       CONFIRMS lane C's values to every printed place; R2_03 and R2_05 use the exact values.
       R2_02's SENSE C effect-size line is SUPERSEDED, not deleted.
  R-3  The resolution table above uses, as its 'ruler', the distance from the average to the
       limit it actually converged to.  That is the right ruler for a power statement but it
       is measured from THIS refuter's runs, not from lane C's; lane C's own K = 1e7 numbers
       are reproduced to 4e-10 by an independent evaluator (R2_03 sec 1 vs C_04.out), so the
       substitution is sound, but it IS a substitution.
  R-4  I did not re-derive C-2's reversal locus, C-3's zero-set criterion beyond a hand check,
       C-8's Boyd-Lawton ladder, or C-1's edge-transport arm.  Anything I say about those is
       a spot-check, not an audit, and is graded as such.
  R-5  NOT LINEAGE-INDEPENDENT.  Opus 5 refuting Opus 5, sharing the phase-reduction device
       and the Mahler-measure machinery.  Discount as one block with lane C and with
       W-07/W-08/W-09.""")
    print("\nDONE.")
