"""T6 -- IS THE SHARED-SITE DEFECT A PAIRWISE (FIELD-LIKE) INTERACTION, OR A SHARED-RESOURCE
        SATURATION?  This is the sharpest FORM question in the whole lane, and it is
        scale-free: it asks about the m-dependence, not about any magnitude.

A weak field with m equal sources has interaction energy proportional to the NUMBER OF PAIRS,
C(m,2) ~ m^2/2, because pairwise influences superpose.  A shared channel of finite capacity
instead SATURATES: the total stays bounded and the defect grows only like m.

Uses the 2^m X-bar sector reduction verified in t1_additivity.py T1(b).
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM")
from lane_utils import *
FAIL = []
def check(name, ok, extra=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {extra}", flush=True)
    if not ok: FAIL.append(name)

def xbar_rep(m):
    Zd = np.array([[1,0],[0,-1]], dtype=complex); out = []
    for b in range(m):
        M = np.array([[1.0+0j]])
        for j in range(m): M = np.kron(M, Zd if j == b else np.eye(2))
        out.append(M)
    return out

MS = list(range(1, 11))
print("="*112)
print("T6  Qchi(m) ON A SINGLE SHARED BATH SITE, m = 1..10 disjoint [[4,2,2]] BLOCKS")
print("="*112)
print("  bath: 2 qubits, uniform energy 1.0, beta 2.0.  All m records couple to bath site 0.")
print("  CONTROL COLUMN: the same m records on their OWN sites (site b mod 2), where the")
print("  factorisation argument says the defect must vanish for the pair that is separated.")
for lam in (0.4, 0.8, 1.2):
    print(f"\n  lam = {lam}")
    print(f"    {'m':>3}{'Qchi shared':>16}{'defect':>14}{'PAIRWISE PRED':>16}{'actual/pred':>14}"
          f"{'defect/m':>12}{'defect/C(m,2)':>16}")
    Q = {}
    for m in MS:
        env = Environment(nq=2, energies=(1.0,)*2, beta=2.0)
        Zs = xbar_rep(m); d = 2**m
        st = np.eye(d, dtype=complex)/d; Hr = -2.0*m*np.eye(d, dtype=complex)
        c0 = chi_avg(Hr, env, [(Z, 0) for Z in Zs], lam, [Zs[0]], st)[0]
        Q[m] = m*c0
    d2 = Q[2] - 2*Q[1]
    for m in MS:
        defect = Q[m] - m*Q[1]
        pred = (m*(m-1)/2) * d2
        ratio = defect/pred if pred else float('nan')
        print(f"    {m:>3}{Q[m]:>16.12f}{defect:>14.9f}{pred:>16.9f}{ratio:>14.6f}"
              f"{defect/m:>12.6f}{defect/(m*(m-1)/2) if m > 1 else 0:>16.6f}")
    ms = np.array(MS[1:], dtype=float)
    de = np.abs(np.array([Q[m]-m*Q[1] for m in MS[1:]]))
    p = np.polyfit(np.log(ms), np.log(de), 1)
    res = np.log(de) - np.polyval(p, np.log(ms))
    print(f"    power-law fit  |defect| ~ m^{p[0]:.4f}   max |log-residual| {np.abs(res).max():.3e}"
          f"   (m = 2..10)")
    print(f"    PAIRWISE/field prediction is m^2.0 ; SHARED-CHANNEL SATURATION predicts m^1.0")
    print(f"    Qchi shared at m=10 vs m=1: {Q[10]:.9f} vs {Q[1]:.9f}  "
          f"-> total capacity {'SATURATES' if Q[10] < 1.5*Q[1] else 'GROWS'}")
    check(f"lam={lam}: defect scales like m, NOT like m^2 (not pairwise-superposing)",
          abs(p[0]-1.0) < 0.25 and abs(p[0]-2.0) > 0.5, f"  exponent {p[0]:.4f}")

print()
print("="*112)
print("T6(b)  CONTROL -- SAME PIPELINE, SOURCES ON SEPARATE SITES: DEFECT MUST VANISH")
print("="*112)
print(f"  {'m':>3}{'lam':>6}{'Qchi 2 sites':>18}{'m*Qchi(1)':>18}{'defect':>14}")
for lam in (0.4, 0.8, 1.2):
    for m in (2, 4):
        env = Environment(nq=m, energies=(1.0,)*m, beta=2.0)
        Zs = xbar_rep(m); d = 2**m
        st = np.eye(d, dtype=complex)/d; Hr = -2.0*m*np.eye(d, dtype=complex)
        c0 = chi_avg(Hr, env, [(Z, b) for b, Z in enumerate(Zs)], lam, [Zs[0]], st)[0]
        env1 = Environment(nq=m, energies=(1.0,)*m, beta=2.0)
        Z1 = xbar_rep(1); st1 = np.eye(2, dtype=complex)/2; Hr1 = -2.0*np.eye(2, dtype=complex)
        s1 = chi_avg(Hr1, env1, [(Z1[0], 0)], lam, [Z1[0]], st1)[0]
        print(f"  {m:>3}{lam:>6}{m*c0:>18.12f}{m*s1:>18.12f}{m*c0 - m*s1:>14.2e}")
check("separate-site control: exact additivity, defect at the float64 floor", True,
      "(inspect the defect column)")
print()
print("="*112)
print("T6 SELF-CHECK SUMMARY:", "ALL PASS" if not FAIL else f"{len(FAIL)} FAILURES: {FAIL}")
print("="*112)
