"""T4 -- SUPERPOSITION OF SOURCES, and the WEAK-COUPLING (leading-order) form.

Records C (probe), A and B (sources), all DERIVED and clause-verified, all on [[10,8,2]].
Influence of a source S on the probe:   Delta_S = chi_C(alone) - chi_C(with S coupled).
Superposition deviation (a genuine THREE-RECORD quantity):
        D = Delta_{A,B} - Delta_A - Delta_B .

A weak field superposes.  The scale-free question is not whether D is zero at lam = O(1) --
it is whether D vanishes FASTER than the individual influences as the coupling weakens, i.e.
whether the LEADING-ORDER response is additive.  That is testable at N = 10 and does not
depend on the size of any coupling constant.

CONTROL (D-15): the SEPARATE-BATH-SITE column, where superposition must be exact, sits in the
same table as the shared-site column.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM")
from lane_utils import *

FAIL = []
def check(name, ok, extra=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {extra}", flush=True)
    if not ok: FAIL.append(name)

n = 10
reps, idx = code_reps(n); d = len(reps)
stab = stab_nn2(n); S, L, _ = derived_logical_span(stab, n)

def rec(kind, i, j):
    v = pauli_vec(n, (i,j), ()) if kind == "X" else pauli_vec(n, (), (i,j))
    assert is_nontrivial_logical(v, S, L, n)
    M = (compress_XX(i,j,n,reps,idx) if kind=="X" else compress_ZZ(i,j,n,reps,idx)).astype(complex)
    assert np.linalg.norm(M@M - np.eye(d)) < 1e-12 and np.linalg.norm(M - M.conj().T) < 1e-12
    return v, M, {i, j}

vC, MC, sC = rec("X", 0, 1)
vA, MA, sA = rec("X", 2, 3)
vB, MB, sB = rec("X", 4, 5)
env = Environment(nq=3, energies=(1.0,)*3, beta=2.0)
st = np.eye(d, dtype=complex)/d; Hr = -2.0*np.eye(d, dtype=complex)

print("="*112)
print("T4(a)  NOISE FLOOR FROM PERMUTATION-EQUIVALENT REPLICAS")
print("="*112)
print("  The exact argument of T2(b) says these configurations must give bit-identical chi.")
print("  Their observed spread IS the float64 noise floor of the whole pipeline.")
reps_cfg = [(("X",0,1),("X",2,3),("X",4,5)), (("X",0,1),("X",4,5),("X",6,7)),
            (("X",2,3),("X",6,7),("X",8,9)), (("X",0,9),("X",1,8),("X",2,7)),
            (("Z",0,1),("Z",2,3),("Z",4,5)), (("X",0,1),("Z",2,3),("Z",4,5))]
vals = []
for cc, aa, bb in reps_cfg:
    _, Mc, _ = rec(*cc); _, Ma, _ = rec(*aa); _, Mb, _ = rec(*bb)
    v = chi_avg(Hr, env, [(Mc,0),(Ma,0),(Mb,0)], 0.8, [Mc], st)[0]
    vals.append(v)
    print(f"    {str(cc):<16}{str(aa):<16}{str(bb):<16}  chi_C = {v:.15f}")
FLOOR = max(vals) - min(vals)
print(f"  NOISE FLOOR (spread over exactly-equivalent replicas) = {FLOOR:.3e}")
check("noise floor is at float64 resolution", FLOOR < 1e-13, f"  {FLOOR:.3e}")

print()
print("="*112)
print("T4(b)  SUPERPOSITION.  SHARED SITE (test) BESIDE SEPARATE SITES (control), SAME TABLE")
print("="*112)
print(f"  {'lam':>6}{'chi_C alone':>16}{'Delta_A':>15}{'Delta_B':>15}{'Delta_AB':>15}"
      f"{'D = dev (TEST)':>17}{'|D|/(|dA|+|dB|)':>18}{'D separate (CTRL)':>19}")
rows = []
for lam in (0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.2):
    c0  = chi_avg(Hr, env, [(MC,0)], lam, [MC], st)[0]
    cA  = chi_avg(Hr, env, [(MC,0),(MA,0)], lam, [MC], st)[0]
    cB  = chi_avg(Hr, env, [(MC,0),(MB,0)], lam, [MC], st)[0]
    cAB = chi_avg(Hr, env, [(MC,0),(MA,0),(MB,0)], lam, [MC], st)[0]
    dA, dB, dAB = c0-cA, c0-cB, c0-cAB
    D = dAB - dA - dB
    # CONTROL: A and B on their OWN bath sites -- influence must be exactly zero, hence
    # superposition exactly exact.
    e0  = chi_avg(Hr, env, [(MC,0)], lam, [MC], st)[0]
    eA  = chi_avg(Hr, env, [(MC,0),(MA,1)], lam, [MC], st)[0]
    eB  = chi_avg(Hr, env, [(MC,0),(MB,2)], lam, [MC], st)[0]
    eAB = chi_avg(Hr, env, [(MC,0),(MA,1),(MB,2)], lam, [MC], st)[0]
    Dc = (e0-eAB) - (e0-eA) - (e0-eB)
    ratio = abs(D)/(abs(dA)+abs(dB)) if (abs(dA)+abs(dB)) > 0 else float('nan')
    rows.append((lam, c0, dA, dB, dAB, D, ratio, Dc))
    print(f"  {lam:>6}{c0:>16.12f}{dA:>15.9f}{dB:>15.9f}{dAB:>15.9f}"
          f"{D:>17.9f}{ratio:>18.9f}{Dc:>19.2e}")
check("CONTROL: superposition EXACT when sources sit on their own bath sites",
      all(abs(r[7]) < 1e-12 for r in rows), f"  max |D_ctrl| = {max(abs(r[7]) for r in rows):.2e}")
check("TEST: superposition FAILS on a shared site, above the noise floor",
      all(abs(r[5]) > 100*FLOOR for r in rows), f"  min |D| = {min(abs(r[5]) for r in rows):.3e}")

print()
print("  SCALING OF THE SUPERPOSITION DEFECT WITH THE COUPLING (this is the FORM question):")
lams = np.array([r[0] for r in rows]); Ds = np.abs(np.array([r[5] for r in rows]))
dAs = np.abs(np.array([r[2] for r in rows]))
m = (Ds > 1e-13) & (dAs > 1e-13)
pD = np.polyfit(np.log(lams[m]), np.log(Ds[m]), 1)
pd = np.polyfit(np.log(lams[m]), np.log(dAs[m]), 1)
resD = np.log(Ds[m]) - np.polyval(pD, np.log(lams[m]))
resd = np.log(dAs[m]) - np.polyval(pd, np.log(lams[m]))
print(f"    |D|      ~ lam^{pD[0]:.4f}   max |log-residual| {np.abs(resD).max():.3e}"
      f"   (n={m.sum()} points, lam in [{lams[m].min()},{lams[m].max()}])")
print(f"    |Delta_A| ~ lam^{pd[0]:.4f}   max |log-residual| {np.abs(resd).max():.3e}")
print(f"    EXPONENT GAP  {pD[0]-pd[0]:.4f}   -- positive means the defect vanishes FASTER than")
print(f"    the individual influence, i.e. the LEADING-ORDER response IS additive.")

print()
print("="*112)
print("T4(c)  DOES THE THREE-RECORD DEFECT D DEPEND ON THE GEOMETRY OF A, B RELATIVE TO C?")
print("="*112)
print(f"  {'A':<8}{'B':<8}{'sep(C,A)':>10}{'sep(C,B)':>10}{'sep(A,B)':>10}"
      f"{'Delta_A':>15}{'Delta_B':>15}{'D':>17}")
lam = 0.8
c0 = chi_avg(Hr, env, [(MC,0)], lam, [MC], st)[0]
Dvals = []
for (ai, aj), (bi, bj) in [((2,3),(4,5)), ((2,3),(8,9)), ((4,5),(6,7)), ((8,9),(6,7)),
                           ((2,9),(3,8)), ((4,7),(5,6)), ((2,3),(2,3))]:
    _, Ma, sa = rec("X", ai, aj); _, Mb, sb = rec("X", bi, bj)
    sca = min(abs(i-j) for i in sC for j in sa) if not (sC & sa) else 0
    scb = min(abs(i-j) for i in sC for j in sb) if not (sC & sb) else 0
    sab = min(abs(i-j) for i in sa for j in sb) if not (sa & sb) else 0
    cA  = chi_avg(Hr, env, [(MC,0),(Ma,0)], lam, [MC], st)[0]
    cB  = chi_avg(Hr, env, [(MC,0),(Mb,0)], lam, [MC], st)[0]
    cAB = chi_avg(Hr, env, [(MC,0),(Ma,0),(Mb,0)], lam, [MC], st)[0]
    D = (c0-cAB) - (c0-cA) - (c0-cB)
    if (ai,aj) != (bi,bj): Dvals.append(D)
    print(f"  {'X%dX%d'%(ai,aj):<8}{'X%dX%d'%(bi,bj):<8}{sca:>10}{scb:>10}{sab:>10}"
          f"{c0-cA:>15.9f}{c0-cB:>15.9f}{D:>17.12f}")
print(f"  spread of D over all DISTINCT-source geometries: {max(Dvals)-min(Dvals):.3e}"
      f"   (noise floor {FLOOR:.3e})")
check("three-record defect D is geometry-independent to the noise floor",
      max(Dvals)-min(Dvals) < 1e-12, f"  spread {max(Dvals)-min(Dvals):.3e}")

print()
print("="*112)
print("T4 SELF-CHECK SUMMARY:", "ALL PASS" if not FAIL else f"{len(FAIL)} FAILURES: {FAIL}")
print("="*112)
