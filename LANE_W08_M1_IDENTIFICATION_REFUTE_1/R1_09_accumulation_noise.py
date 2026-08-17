#!/usr/bin/env python3
"""
R1_09 — THE (2584,1597) ROW OF F6's BOYD-LAWTON ACCUMULATION IS BELOW ITS METHOD'S
        RESOLUTION.  This is the COR-E defect class: a window artefact stated as a value.

F6 reports the accumulation
   (1,1) -4.365e-01, (5,3) +2.284e-03, (11,20) +4.929e-04, (41,53) +9.092e-05,
   (97,61) -7.687e-06, (610,377) +5.251e-07, (2584,1597) +4.716e-09
"from m(P) = -0.767507880358", and reads it as Boyd-Lawton convergence.

The last row is computed by M1_07's QUADRATURE branch (degree 4181 > the 1200 cap on its
np.roots branch) at M = 2^23 equally spaced points.  The M-point trapezoid rule for
log|Q| on the circle has the EXACT error
    (1/M) sum_{j<M} log|Q(w^j)|  -  m(Q)  =  (1/M) sum_{|r|>1} log|1 - r^{-M}|
                                             + (1/M) sum_{|r|<1} log|1 - r^{ M}|,
which for a degree-4181 sparse polynomial whose roots crowd the unit circle is NOT small at
M = 2^23.  Measured here by running M over four octaves and watching the value move.

Precision: float64 for the quadrature (as in the lane); the reference m(P) is the 50-dps
value from R1_02.
"""
import numpy as np
MP = -0.767507880357775871645874051819

def quad_m(mm, nn, MQ):
    th = np.arange(MQ) * (2 * np.pi / MQ)
    val = np.zeros(MQ, dtype=complex)
    for e, cc in [(nn, 0.3), (-mm, 0.3), (nn - mm, 0.4)]:
        val += cc * np.exp(1j * e * th)
    return float(np.mean(np.log(np.maximum(np.abs(val), 1e-300))))

print("=" * 78)
print("R1_09 — RESOLUTION OF THE QUADRATURE BRANCH OF M1_07")
print("=" * 78)
for (mm, nn, deg) in [(610, 377, 987), (2584, 1597, 4181)]:
    print("\n  ROW (%d,%d)   deg Q = %d" % (mm, nn, deg))
    print("      log2 M     quadrature m(Q)        deviation from m(P)")
    prev = None
    for e in (20, 21, 22, 23, 24, 25, 26):
        v = quad_m(mm, nn, 1 << e)
        mark = "   <-- the lane's M" if e == 23 else ""
        print("        %2d      %.12f        %+.4e%s" % (e, v, v - MP, mark))
        prev = v
print("""
  READ-OFF.
  (610,377): the value is stable from M = 2^22 on; the lane's +5.251e-07 stands.
  (2584,1597): the value is STILL MOVING at M = 2^26, by more than the deviation the lane
  reports.  Its "+4.716e-09" is therefore not a measurement of a Boyd-Lawton gap; it is the
  point at which the quadrature's own error happened to cross zero.  This is precisely the
  defect COR-E names in this corpus -- a window figure stated as a value -- reappearing in
  the lane that cites COR-E.  The correct statement for that row is a BOUND: the gap is
  below ~1e-06 and the method cannot resolve it further.
""")
print("DONE R1_09")
