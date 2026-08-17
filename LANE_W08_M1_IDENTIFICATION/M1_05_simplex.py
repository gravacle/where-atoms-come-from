#!/usr/bin/env python3
"""
M1_05 — m(P) OVER K1'S OWN SIMPLEX  { p00 = 0, p10 + p01 + p11 = 1, p >= 0 }.
WHERE DOES IT VANISH, AND WHERE IS IT BOUNDED AWAY FROM 0?

p00 = 0 is forced by K1's incidence, so the parameter space is the 2-SIMPLEX in
(p10, p01, p11), not the 3-simplex.  Its three vertices have carrier meanings:
      p11 = 1   the whole ready state sits on v0, THE ROOT   (in both loops)
      p10 = 1   the whole ready state sits on {v1,v2}        (face triangle only)
      p01 = 1   the whole ready state sits on {v3,v4}        (unfilled triangle only)

Established here (grid + closed form, float64; the closed forms are exact):
  S1  m(P) = log( max(p10,p01,p11) )  EXACTLY, whenever max > 1/2 (the non-triangle region).
  S2  m(P) <= -log 2  on the whole triangle region max <= 1/2, with equality only on its
      boundary; the minimum over the simplex is at the centroid,
      m = log(1/3) + m(1+x+y) = -1.0986122886681098 + 0.3230659472... = -0.7755463414...
  S3  m(P) = 0  <=>  max = 1  <=>  ONE CLASS CARRIES THE WHOLE READY STATE.  Three points.
  S4  UNIFORM BOUND:  m(P) <= -(1 - max(p10,p01,p11)).  So lambda is at least the weight
      that is NOT in the heaviest class.  There is NO uniform gap: m -> 0 linearly as the
      state concentrates on one class.
  S5  Kronecker / Lehmer, placed correctly: they concern INTEGER coefficients.  Here the
      coefficients are probabilities and the vanishing question is settled elementarily by
      Jensen (S3).  Lehmer's problem is NOT the relevant background for "is lambda bounded
      away from 0" over the simplex — the answer there is NO, and no gap theorem can help.
      It IS relevant in one specialisation: for a ready state UNIFORM on n cells, n*P has
      non-negative integer coefficients, m(nP) >= 0 (Kronecker), and max <= (n-1)/n unless
      the state is confined to one class, so lambda <= -1/n.  On K1 with n = 5: lambda <= -1/5.
"""
import numpy as np
from M1_02_mahler_machinery import m_R1, m_CM

print("=" * 78)
print("M1_05 — m(P) OVER THE SIMPLEX p00 = 0")
print("=" * 78)

NG = 400
rows = []
best = (1e9, None)
worst_S1 = 0.0
worst_S4 = -1e9
n_nontri = 0
zeros = []
for i in range(NG + 1):
    for j in range(NG + 1 - i):
        p10 = i / NG
        p01 = j / NG
        p11 = 1.0 - p10 - p01
        if p11 < -1e-15:
            continue
        p11 = max(p11, 0.0)
        m = m_R1(0.0, p10, p01, p11)
        mx = max(p10, p01, p11)
        if m < best[0]:
            best = (m, (p10, p01, p11))
        if mx > 0.5:
            n_nontri += 1
            worst_S1 = max(worst_S1, abs(m - np.log(mx)))
        worst_S4 = max(worst_S4, m + (1.0 - mx))          # must stay <= 0
        if m > -1e-9:
            zeros.append((p10, p01, p11, m))

print("\nS1  max > 1/2  =>  m(P) = log(max) exactly.")
print("    worst |m(P) - log(max)| over the %d grid points with max > 1/2 : %.3e"
      % (n_nontri, worst_S1))
print("S4  worst value of  m(P) + (1 - max)   over the whole grid : %.3e   (must be <= 0)"
      % worst_S4)
print("S3  grid points with m(P) > -1e-9 :")
for z in zeros:
    print("      (p10,p01,p11) = (%.3f, %.3f, %.3f)   m = %+.3e" % z)
print("    -> exactly the three vertices.  m(P) = 0 iff one class carries all the weight.")
print("S2  minimum on the grid: m = %.12f at (p10,p01,p11) = (%.4f, %.4f, %.4f)"
      % (best[0], *best[1]))
cent = m_R1(0.0, 1 / 3, 1 / 3, 1 / 3)
print("    centroid closed form  log(1/3) + m(1+x+y) = %.12f ; computed %.12f ; diff %.2e"
      % (np.log(1 / 3) + m_CM(1, 1, 1), cent, abs(np.log(1 / 3) + m_CM(1, 1, 1) - cent)))
print("    boundary of the triangle region (max = 1/2): m = log(1/2) = %.12f" % np.log(0.5))
print("    so on the whole triangle region  m(P) in [%.9f, %.9f]" % (best[0], np.log(0.5)))

print("\nTHE THREE VERTICES, IN CARRIER TERMS")
print("    p11 = 1  : ready state on v0 alone, THE ROOT.  m = 0.  This is W-01's registered")
print("               'the root can never fire', recovered as a Mahler-measure statement.")
print("    p10 = 1  : ready state on {v1,v2} (filled triangle, not the cycle).  m = 0.")
print("    p01 = 1  : ready state on {v3,v4} (unfilled cycle, not the face).  m = 0.")
print("    Any other ready state whatsoever has m(P) < 0 STRICTLY.")

print("\nCARRIER-SPECIFIC VALUES")
for name, pv in [("S1 published p=(1/2,0,0,1/4,1/4) -> (0, 1/2, 1/2)", (0.0, 0.5, 0.5)),
                 ("uniform on all 5 vertices        -> (2/5, 2/5, 1/5)", (0.4, 0.4, 0.2)),
                 ("S3/S4 generic                    -> (3/10, 3/10, 2/5)", (0.3, 0.3, 0.4)),
                 ("uniform on the 4 non-root vtcs   -> (1/2, 1/2, 0)", (0.5, 0.5, 0.0)),
                 ("root only                        -> (0, 0, 1)", (0.0, 0.0, 1.0))]:
    m = m_R1(0.0, *pv)
    print("   %-52s m(P) = %+.12f   bound -(1-max) = %+.6f"
          % (name, m, -(1 - max(pv))))

print("\nNO UNIFORM GAP — m(P) near a vertex (approach along p11 -> 1):")
for e in (1e-1, 1e-2, 1e-3, 1e-4, 1e-6):
    m = m_R1(0.0, e / 2, e / 2, 1 - e)
    print("    1 - p11 = %8.1e   m(P) = %+.10e   ratio m/(1-p11) = %.6f" % (e, m, m / e))
print("    m(P) ~ -(1-p11) as the state concentrates: the rate degrades CONTINUOUSLY to 0.")
print("    A Lehmer-type gap cannot exist over a simplex of real weights.")

print("\nDONE M1_05")
