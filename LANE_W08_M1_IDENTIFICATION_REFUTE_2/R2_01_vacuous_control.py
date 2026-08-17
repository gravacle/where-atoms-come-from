#!/usr/bin/env python3
"""
R2_01 — ATTACK ON M1_03's "DIFFERENT POINT ON THE RESONANCE LINE" CONTROL.

THE TARGET'S CLAIM (finding F10): the (11,20) subtorus limit is "reproduced from a DIFFERENT
point on the same resonance line (-0.767016018 at N = 1e6), confirming that the limit depends
only on the primitive relation vector."

THE CODE (M1_03, lines near the end):
    A0 = int(round(INVPI * S))                 # INVPI = 1.0/np.pi ,  S = 2**35
    ...
    A2 = int(round(0.31830988618379 * S))      # "a different irrational point on the same line"

ISOLATION LEDGER FOR THIS BLOCK
  HELD FIXED : relation vector (11,20); D_C = 20*2^35; weights (0.3,0.3,0.4); K = 1e6; code path.
  THE ONE THING THAT IS SUPPOSED TO MOVE : the point on the resonance circle, i.e. A0 -> A2.
  WHAT ACTUALLY MOVES : nothing.  0.31830988618379 is 1/pi truncated at 14 digits; the
  truncation error 6.7e-16 times S = 3.44e10 is 2.3e-5, far below the 0.5 needed to change the
  rounded integer.  A2 == A0 EXACTLY, so the "control" re-runs the identical orbit.
Part 2 supplies the control that was intended: genuinely different points on the same circle.
"""
import numpy as np

S = 2 ** 35
D_C = 20 * S
P10, P01, P11 = 0.3, 0.3, 0.4
CHUNK = 10 ** 6


def orbit_average(A_num, B_num, D, K, checkpoints):
    tot, done, out = 0.0, 0, {}
    minabs = np.inf
    while done < K:
        n = min(CHUNK, K - done)
        k = np.arange(done + 1, done + n + 1, dtype=np.int64)
        assert k.max() * max(A_num, B_num) < 2 ** 62, "int64 overflow guard"
        fa = ((k * A_num) % D).astype(np.float64) / D
        fb = ((k * B_num) % D).astype(np.float64) / D
        x = np.exp(2j * np.pi * fa)
        y = np.exp(2j * np.pi * fb)
        a = np.abs(P10 * x + P01 * y + P11 * x * y)
        minabs = min(minabs, float(a.min()))
        cs = np.cumsum(np.log(np.maximum(a, 1e-323)))
        for cp in checkpoints:
            if done < cp <= done + n:
                out[cp] = (tot + cs[cp - done - 1]) / cp
        tot += float(cs[-1])
        done += n
    return out, minabs


print("=" * 78)
print("R2_01 — THE (11,20) 'DIFFERENT POINT' CONTROL IS THE SAME POINT")
print("=" * 78)

A0 = int(round((1.0 / np.pi) * S))
A2 = int(round(0.31830988618379 * S))
print("\n  1/pi as float64                    = %.17f" % (1.0 / np.pi))
print("  the literal used for the 'other' pt = %.17f" % 0.31830988618379)
print("  difference                          = %.3e" % abs(1.0 / np.pi - 0.31830988618379))
print("  S = 2^35 = %d ;  difference * S = %.4e   (needs > 0.5 to change round())"
      % (S, abs(1.0 / np.pi - 0.31830988618379) * S))
print("\n  A0 = round(1/pi * S)                = %d" % A0)
print("  A2 = round(0.31830988618379 * S)    = %d" % A2)
print("  A2 - A0                             = %d" % (A2 - A0))
print("  ==> IDENTICAL.  THE CONTROL IS THE ORIGINAL ORBIT RE-RUN.  IT COULD NOT HAVE FAILED.")

# reproduce both, to show the printed numbers coincide to the last digit
for tag, A in (("A0 (case C)", A0), ("A2 ('different point')", A2)):
    An = (-20 * A) % D_C
    Bn = (11 * A) % D_C
    res, mn = orbit_average(An, Bn, D_C, 10 ** 6, [10 ** 6])
    print("  %-24s  (1/N)sum log|Z_k| at N=1e6 = %.12f" % (tag, res[10 ** 6]))

# ---------------------------------------------------------------- the control that was intended
print("\n" + "-" * 78)
print("PART 2 — THE CONTROL THAT WAS INTENDED: GENUINELY DIFFERENT POINTS ON THE SAME CIRCLE.")
print("  H = {(z^20, z^-11)}.  Parametrise by alpha = -t, beta = 11 t / 20 for t in (0,1),")
print("  i.e. A_num = (-20*T) % D, B_num = (11*T) % D with T = round(t*S) and t VARIED.")
print("  HELD FIXED: relation (11,20), weights, D, K, code path.  MOVES: t only.")
print("  Exact subtorus value (Jensen, computed in part 3) is the prediction.")
ts = [1.0 / np.pi, 0.123456789012345, 0.5555555555555555 * 0.7, np.sqrt(2) - 1.0,
      np.exp(1.0) - 2.0, 0.8090169943749474]
print("\n     t                     N=1e5           N=1e6           min|Z_k|")
for t in ts:
    T = int(round(t * S))
    An = (-20 * T) % D_C
    Bn = (11 * T) % D_C
    res, mn = orbit_average(An, Bn, D_C, 10 ** 6, [10 ** 5, 10 ** 6])
    print("   %.15f   %14.9f  %14.9f   %.3e" % (t, res[10 ** 5], res[10 ** 6], mn))

# ---------------------------------------------------------------- part 3: the exact prediction
coef = np.zeros(32)
coef[31] = P10
coef[20] = P11
coef[0] = P01
rts = np.roots(coef[::-1])
mH = np.log(abs(P10)) + float(np.sum(np.log(np.maximum(np.abs(rts), 1.0))))
print("\n  EXACT subtorus value m(p10 z^31 + p11 z^20 + p01) by Jensen = %.12f" % mH)
print("  (this is the prediction every row of Part 2 must approach)")
print("\nDONE R2_01")
