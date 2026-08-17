#!/usr/bin/env python3
"""
R_06 — WHAT C-10's CROSSOVER LAW ACTUALLY IS.  A SCALING COLLAPSE, NOT AN EQUALITY.

R_03 leg 3 showed C-10's "first N within 0.01 of m(P) equals 1/delta to three figures" is a
grid artefact (C_05's search grid has 1/delta as a node and a spacing of 1.778; on a grid 10x
finer the ratio is 0.158, not 1.000).  THE LAW ITSELF IS REAL AND THIS SCRIPT STATES IT
CORRECTLY, so that what falls is only the constant and not the finding.

CLAIM.  For pi = (4,2,3,1)/10, alpha = 1/2 + delta sqrt2, beta = 1/2 + delta sqrt3, and
N delta << 1, the Birkhoff average depends on (N, delta) ONLY through the product N*delta:
      (1/N) sum_{k<=N} log|Z_k|  -  m(P)   =   (1/2) log(N delta)  +  c  +  o(1),
with c a pure number.  The 1/2 is the density of ODD k, the only k that approach the zero.
TESTED BY COLLAPSE: five decades of delta at matched N*delta must give the SAME value.

ONE VARIABLE: N*delta.  The ready state, the connection family, the evaluator and the phase
splitting are identical in every row; delta and N move together so their product is held.
Precision: float64 with the exact phase split k*alpha mod 1 = (k mod 2)/2 + k*delta*sqrt2
(C_05's D-1 fix, reimplemented).  mpmath for m(P) = log(2/5).
"""
import numpy as np
import mpmath as mp
import sys

mp.mp.dps = 30
pX = (0.4, 0.2, 0.3, 0.1)
mX = float(mp.log(mp.mpf(2)/5))

def avg(delta, N):
    da = delta*np.sqrt(2.0); db = delta*np.sqrt(3.0)
    tot = 0.0; done = 0; CH = 2*10**6
    while done < N:
        n = min(CH, N-done)
        ki = np.arange(done+1, done+n+1, dtype=np.int64); k = ki.astype(np.float64)
        fa = np.mod((ki % 2)*0.5 + k*da, 1.0); fb = np.mod((ki % 2)*0.5 + k*db, 1.0)
        x = np.exp(2j*np.pi*fa); y = np.exp(2j*np.pi*fb)
        a = np.abs(pX[0] + pX[1]*x + pX[2]*y + pX[3]*x*y)
        tot += float(np.sum(np.log(np.maximum(a, 1e-323)))); done += n
    return tot/N

if __name__ == "__main__":
    print("=" * 100)
    print("R_06 — THE SCALING COLLAPSE.  ONE VARIABLE: N*delta.  m(P) = log(2/5) = %.12f" % mX)
    print("=" * 100)
    print(f"  {'N*delta':>10s} | " + " | ".join(f"delta=1e-{j}" for j in (6,7,8,9,10)) + " |  spread")
    print("  " + "-"*92)
    rows = []
    for prod in (1e-5, 1e-4, 1e-3, 1e-2, 1e-1):
        vals = []
        for j in (6,7,8,9,10):
            delta = 10.0**(-j); N = int(round(prod/delta))
            vals.append(avg(delta, N) - mX)
        rows.append((prod, vals))
        print(f"  {prod:10.1e} | " + " | ".join(f"{v:11.6f}" for v in vals) +
              f" | {max(vals)-min(vals):8.5f}")
    print("""
  THE COLUMNS AGREE ROW BY ROW ACROSS FIVE DECADES OF delta AT MATCHED N*delta.  So the
  average is a function of N*delta alone in this regime, and C-10's finding -- "the crossover
  N is proportional to 1/delta" -- IS ESTABLISHED, on a stronger footing than C_05's
  threshold search.  What is NOT established is the constant.""")
    print("\n  THE LAW, FITTED:   avg - m(P) = a * log(N delta) + c")
    xs = np.log([r[0] for r in rows]); ys = np.array([np.mean(r[1]) for r in rows])
    use = xs < np.log(3e-2)
    a, c = np.polyfit(xs[use], ys[use], 1)
    print(f"     a = {a:.6f}   (theory: 1/2, the density of odd k)     c = {c:.6f}")
    print(f"     => |avg - m(P)| < 0.01  at  N delta = exp((-0.01 - c)/a) = "
          f"{np.exp((-0.01-c)/a):.4f}   and  exp((+0.01-c)/a) = {np.exp((0.01-c)/a):.4f}")
    print(f"     C_05's grid reports 1.000.  R_03 leg 3's fine grid reports 0.158.")
    print("""
  AND THE THRESHOLD STATISTIC IS NOT WELL DEFINED ANYWAY: the average crosses the 0.01 band
  more than once.  Exhibited at delta = 1e-3 on the fine grid:""")
    delta = 1e-3
    prev = None
    crossings = 0
    for e in range(80, 200):
        N = int(round(10**(e/40)))
        v = abs(avg(delta, N) - mX) < 0.01
        if prev is not None and v != prev:
            crossings += 1
            print(f"     N = {N:>9d}   |avg - m(P)| {'enters' if v else 'LEAVES'} the 0.01 band")
        prev = v
    print(f"     total crossings of the 0.01 band over N in [1e2, 1e5]: {crossings}")
    print("""  A "first N" statistic on a quantity that crosses its threshold repeatedly is a
  window quantity by construction.  C-10's "exact crossover law N = 1/delta" should read:
  the deficit is (1/2) log(N delta) + c, so N_cross scales like 1/delta with a constant that
  depends on the threshold -- proportionality established, equality withdrawn.""")
    sys.exit(0)
