#!/usr/bin/env python3
"""
R2_03 — TWO ISOLATION ATTACKS ON THE DIOPHANTINE ARM (M1_03 case A, findings F5 and F10).

ATTACK 1 (self-flag 4, discharged or not).  The target reports case A deviations
  N=1e6 -> +6.52e-09,  N=1e7 -> +6.42e-08 (WORSE), attributes the non-monotonicity to
  "discrepancy oscillation" and says "a residual precision artefact at large k is not
  excluded".  I exclude it here.  The target's phase scheme represents alpha to 2^-39 EXACTLY
  plus a float64 residue; the number it actually iterates is the float64 value of -2^(1/3),
  which differs from the true cube root by ~1e-17, so by k = 1e7 the orbit has drifted ~1e-10
  in phase.  I recompute with a TWO-LIMB exact int64 scheme modulo 2^62 whose numerator comes
  from an EXACT integer cube root: representation error < 2.2e-19, drift by k=1e7 < 2.2e-12.
  ISOLATION LEDGER: HELD FIXED — the pair (-2^(1/3), 4^(1/3)), the weights (0.3,0.3,0.4), the
  observable, the k-grid, the evaluator.  THE ONE THING THAT MOVES: the phase representation
  (2^39 + float64 residue  ->  exact 2^62 from an integer cube root).

ATTACK 2 (a confound the target did not flag).  Finding F5's "empirical face" — the
  min-distance-to-Z(P) table with sqrt(N)*dist = 0.22 .. 0.35 — is quoted in support of the
  CONVERGENCE CLAIM of F10, which is made at pi = (0, 0.3, 0.3, 0.4).  But that table is
  computed in M1_06 at pi = (0, 1/3, 1/3, 1/3): DIFFERENT WEIGHTS, DIFFERENT ZERO SET,
  DIFFERENT P.  The isolation ledger's COMPARISON 1 declares the ready state HELD FIXED at
  (0,0.3,0.3,0.4); the Diophantine diagnostic that licenses it was run at another state.
  Recomputed here at the weights the claim is actually about.
  ISOLATION LEDGER: HELD FIXED — the Schmidt pair, the k-grid, the evaluator.  MOVES: the
  weight vector, (1/3,1/3,1/3) -> (0.3,0.3,0.4), which is the thing F5 silently changed.

Precision: float64 for |P|; phases EXACT modulo 2^62 (two-limb int64); cube roots EXACT
integers.  Every intermediate in the modular reduction is asserted < 2^63.
"""
import numpy as np

P10, P01, P11 = 0.3, 0.3, 0.4
TWO62 = 1 << 62
SHIFT = 31
MASK = (1 << SHIFT) - 1


def iroot(n, k):
    """exact floor(n^(1/k)) for integers."""
    if n < 0:
        raise ValueError
    x = int(round(n ** (1.0 / k)))
    while x ** k > n:
        x -= 1
    while (x + 1) ** k <= n:
        x += 1
    return x


def frac_num_2_62(which):
    """numerator A with A/2^62 = frac(target) to within 2^-62, EXACT integer arithmetic.
       which = 'a' : frac(-2^(1/3)) = 2 - 2^(1/3)
       which = 'b' : frac(4^(1/3))  = 4^(1/3) - 1"""
    # compute floor(2^(1/3) * 2^62) exactly: floor((2 * 2^186)^(1/3))
    c2 = iroot(2 * (1 << 186), 3)          # ~ 2^(1/3) * 2^62
    c4 = iroot(4 * (1 << 186), 3)          # ~ 4^(1/3) * 2^62
    if which == 'a':
        return 2 * TWO62 - c2              # (2 - 2^(1/3)) * 2^62
    return c4 - TWO62                      # (4^(1/3) - 1) * 2^62


def phases_2_62(A, kmin, kmax):
    """frac(k*A/2^62) for k in [kmin,kmax], EXACT modular reduction in int64."""
    A_hi, A_lo = A >> SHIFT, A & MASK
    k = np.arange(kmin, kmax + 1, dtype=np.int64)
    assert int(k.max()) * A_hi < 2 ** 62, "limb overflow"
    assert int(k.max()) * A_lo < 2 ** 62, "limb overflow"
    hi = ((k * A_hi) % (1 << SHIFT)) << SHIFT
    lo = k * A_lo
    val = (hi + lo) % TWO62                       # exact  (k*A) mod 2^62
    # split so float64 keeps all of it: top 31 bits exact, bottom 31 bits scaled
    top = (val >> SHIFT).astype(np.float64) / (1 << SHIFT)
    bot = (val & MASK).astype(np.float64) / float(TWO62)
    return top, bot


print("=" * 78)
print("R2_03 — PRECISION AND WEIGHTS UNDER THE DIOPHANTINE ARM")
print("=" * 78)

A = frac_num_2_62('a')
B = frac_num_2_62('b')
print("\nEXACT numerators (integer cube roots, no float anywhere):")
print("  A = floor((2 - 2^(1/3)) * 2^62) = %d" % A)
print("  B = floor((4^(1/3) - 1) * 2^62) = %d" % B)
print("  A/2^62 = %.17f   float64 (-2^(1/3))%%1 = %.17f   diff = %.3e"
      % (A / TWO62, (-(2.0 ** (1 / 3))) % 1.0, abs(A / TWO62 - ((-(2.0 ** (1 / 3))) % 1.0))))
print("  B/2^62 = %.17f   float64 (4^(1/3))%%1  = %.17f   diff = %.3e"
      % (B / TWO62, (4.0 ** (1 / 3)) % 1.0, abs(B / TWO62 - ((4.0 ** (1 / 3)) % 1.0))))
print("  -> the target's iterated alpha differs from the true one by ~1e-17; by k=1e7 that is")
print("     a phase drift of ~1e-10.  Mine is < 2.2e-19, drift < 2.2e-12.")


def run(p10, p01, p11, K, checkpoints, zeros=None):
    tot, done = 0.0, 0
    out, mind = {}, {}
    CH = 10 ** 6
    md = np.inf
    mz = np.inf
    while done < K:
        n = min(CH, K - done)
        ta, ba = phases_2_62(A, done + 1, done + n)
        tb, bb = phases_2_62(B, done + 1, done + n)
        x = np.exp(2j * np.pi * ta) * np.exp(2j * np.pi * ba)
        y = np.exp(2j * np.pi * tb) * np.exp(2j * np.pi * bb)
        az = np.abs(p10 * x + p01 * y + p11 * x * y)
        # DEFECT FOUND ON RE-READ OF MY OWN CODE AND RECORDED RATHER THAN SILENTLY FIXED:
        # the first version took md = min(md, az.min()) over the WHOLE 1e6 chunk before the
        # checkpoint loop, so the min-|Z| and min-dist columns at N = 1e3,1e4,1e5 were minima
        # over k <= 1e6 -- they reproduced M1_06's N=1e6 row at every one of the first three
        # rows and would have supported a false claim that M1_06's table is wrong.  M1_06's
        # table is RIGHT.  Fixed here with a running (cumulative-minimum) accumulator.
        runmin = np.minimum.accumulate(az)
        if zeros is not None:
            fa = ta + ba
            fb = tb + bb
            dbest = None
            for (gx, gy) in zeros:
                dx = np.abs(((fa - gx + 0.5) % 1.0) - 0.5)
                dy = np.abs(((fb - gy + 0.5) % 1.0) - 0.5)
                d = np.hypot(dx, dy)
                dbest = d if dbest is None else np.minimum(dbest, d)
            rundist = np.minimum.accumulate(dbest)
        cs = np.cumsum(np.log(np.maximum(az, 1e-323)))
        for cp in checkpoints:
            if done < cp <= done + n:
                out[cp] = (tot + cs[cp - done - 1]) / cp
                mind[cp] = (min(md, float(runmin[cp - done - 1])),
                            min(mz, float(rundist[cp - done - 1])) if zeros is not None else np.inf)
        md = min(md, float(runmin[-1]))
        if zeros is not None:
            mz = min(mz, float(rundist[-1]))
        tot += float(cs[-1])
        done += n
    return out, mind


import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W08_M1_IDENTIFICATION")
from M1_02_mahler_machinery import m_R1

CHECK = [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]
mP = m_R1(0.0, P10, P01, P11)
print("\nATTACK 1 — case A recomputed at 2^62 exact phases.  m(P) = %.12f" % mP)
res, mnd = run(P10, P01, P11, 10 ** 7, CHECK)
print("      N          this lane (2^62)      dev from m(P)      target (2^39+float)   dev")
TARGET = {10 ** 3: -0.768194820808, 10 ** 4: -0.767573512643, 10 ** 5: -0.767501598430,
          10 ** 6: -0.767507873838, 10 ** 7: -0.767507816121}
for cp in CHECK:
    print("  %9d   %18.12f   %+.3e     %18.12f   %+.3e"
          % (cp, res[cp], res[cp] - mP, TARGET[cp], TARGET[cp] - mP))
print("  agreement between the two phase representations at N=1e7: %.3e"
      % abs(res[10 ** 7] - TARGET[10 ** 7]))
print("  VERDICT ON SELF-FLAG 4: the N=1e7 excursion is NOT a precision artefact; it survives")
print("  a representation 10^3 times sharper.  The target's guess was right and is now checked.")

# ---------------------------------------------------------------- ATTACK 2
print("\n" + "-" * 78)
print("ATTACK 2 — F5's min-distance table is computed at THE WRONG WEIGHTS for F10's claim.")


def zeros_of_P(p10, p01, p11):
    """P = p10 x + p01 y + p11 xy = 0 on T^2.  Divide by y: p10 x/y + p01 + p11 x = 0.
       With x = e(fx), y = e(fy): need p01 + p11 x = -p10 x/y, moduli must match:
       |p01 + p11 x| = p10.  Solve for cos(2 pi fx), then fy from the phase."""
    out = []
    c = (p10 ** 2 - p01 ** 2 - p11 ** 2) / (2 * p01 * p11)
    if -1.0 <= c <= 1.0:
        for s in (+1.0, -1.0):
            fx = s * np.arccos(c) / (2 * np.pi)
            x = np.exp(2j * np.pi * fx)
            lhs = p01 + p11 * x                      # modulus p10
            # need p10 * x / y = -lhs  ->  y = -p10 * x / lhs
            y = -p10 * x / lhs
            fy = np.angle(y) / (2 * np.pi)
            out.append((fx % 1.0, fy % 1.0))
    return out


for (name, p) in (("F5/M1_06's weights (1/3,1/3,1/3)", (1 / 3, 1 / 3, 1 / 3)),
                  ("F10/M1_03's weights (0.3,0.3,0.4)", (0.3, 0.3, 0.4))):
    zs = zeros_of_P(*p)
    resid = [abs(p[0] * np.exp(2j * np.pi * a) + p[1] * np.exp(2j * np.pi * b)
                 + p[2] * np.exp(2j * np.pi * (a + b))) for (a, b) in zs]
    print("\n  %s" % name)
    print("    zeros on T^2 (additive coords): %s   residuals %s"
          % (["(%.6f,%.6f)" % z for z in zs], ["%.1e" % r for r in resid]))
    r2, m2 = run(p[0], p[1], p[2], 10 ** 7, CHECK, zeros=zs)
    print("         N        min|Z_k|      min dist to Z(P)   sqrt(N)*dist    (1/N)sum log|Z_k|")
    for cp in CHECK:
        md, mz = m2[cp]
        print("    %9d   %.6e     %.6e       %8.4f      %.9f"
              % (cp, md, mz, np.sqrt(cp) * mz, r2[cp]))
    print("    m(P) = %.12f   dev at N=1e7 = %+.3e" % (m_R1(0.0, *p), r2[10 ** 7] - m_R1(0.0, *p)))
print("\n  READ-OFF: the N^{-1/2} face DOES hold at (0.3,0.3,0.4) too, so F5's conclusion")
print("  survives the correction — but the evidence as published is from a different P, and")
print("  the ledger's 'ready state HELD FIXED' does not cover it.  Recorded, not scored as fatal.")
print("\nDONE R2_03")
