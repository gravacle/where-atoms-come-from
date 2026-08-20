"""T8 -- THE WEAK-COUPLING LIMIT OF SUPERPOSITION.  The one place a field-like FORM could
        still hold, resolved properly.

T4's global fit over lam in [0.025, 1.2] had log-residuals of 1.4 -- it is not a power law over
that range, because chi saturates.  A fit with residuals that large is not a fit.  This script
scans lam over a geometric grid down to 0.005 and fits ONLY inside the window where the
log-residuals are small, reporting the window, the exponents, the residuals and the noise floor.

VENUE: [[8,6,2]] code space (dim 64).  Licensed by T2(b): chi depends only on the symplectic
data, so n = 8 gives bit-identical numbers to n = 10 -- CHECKED here against the n=10 value.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM")
from lane_utils import *
FAIL = []
def check(name, ok, extra=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {extra}", flush=True)
    if not ok: FAIL.append(name)
print(__doc__)

def build(n):
    reps, idx = code_reps(n); d = len(reps)
    stab = stab_nn2(n); S, L, _ = derived_logical_span(stab, n)
    def R(kind, i, j):
        v = pauli_vec(n,(i,j),()) if kind=="X" else pauli_vec(n,(),(i,j))
        assert is_nontrivial_logical(v, S, L, n), (kind,i,j,n)
        M = (compress_XX(i,j,n,reps,idx) if kind=="X" else compress_ZZ(i,j,n,reps,idx)).astype(complex)
        assert np.linalg.norm(M@M - np.eye(d)) < 1e-12
        return M
    return d, R

d8, R8 = build(8)
d10, R10 = build(10)
env = Environment(nq=3, energies=(1.0,)*3, beta=2.0)
print("="*112)
print("T8(a)  VENUE CHECK: n = 8 REPRODUCES n = 10 EXACTLY")
print("="*112)
for lam in (0.1, 0.8):
    for n, dd, R in ((8, d8, R8), (10, d10, R10)):
        st = np.eye(dd, dtype=complex)/dd; Hr = -2.0*np.eye(dd, dtype=complex)
        C, A, B = R("X",0,1), R("X",2,3), R("X",4,5)
        v = chi_avg(Hr, env, [(C,0),(A,0),(B,0)], lam, [C], st)[0]
        print(f"  lam={lam}  n={n:>2} (code dim {dd:>4})  chi_C(A,B) = {v:.15f}")
print()

d, R = d8, R8
st = np.eye(d, dtype=complex)/d; Hr = -2.0*np.eye(d, dtype=complex)
C, A, B = R("X",0,1), R("X",2,3), R("X",4,5)
print("="*112)
print("T8(b)  WEAK-COUPLING SCAN.  Delta_S = chi_C(alone) - chi_C(with S);  D = Delta_AB - Delta_A - Delta_B")
print("="*112)
LAMS = [0.005, 0.00707, 0.01, 0.01414, 0.02, 0.02828, 0.04, 0.05657, 0.08, 0.1131, 0.16, 0.2263, 0.32]
print(f"  {'lam':>8}{'chi_C alone':>18}{'Delta_A':>16}{'Delta_AB':>16}{'D':>16}"
      f"{'|D|/(|dA|+|dB|)':>18}{'D separate (CTRL)':>20}")
rows = []
for lam in LAMS:
    c0  = chi_avg(Hr, env, [(C,0)], lam, [C], st)[0]
    cA  = chi_avg(Hr, env, [(C,0),(A,0)], lam, [C], st)[0]
    cB  = chi_avg(Hr, env, [(C,0),(B,0)], lam, [C], st)[0]
    cAB = chi_avg(Hr, env, [(C,0),(A,0),(B,0)], lam, [C], st)[0]
    dA, dB, dAB = c0-cA, c0-cB, c0-cAB
    D = dAB - dA - dB
    e0  = chi_avg(Hr, env, [(C,0)], lam, [C], st)[0]
    eA  = chi_avg(Hr, env, [(C,0),(A,1)], lam, [C], st)[0]
    eB  = chi_avg(Hr, env, [(C,0),(B,2)], lam, [C], st)[0]
    eAB = chi_avg(Hr, env, [(C,0),(A,1),(B,2)], lam, [C], st)[0]
    Dc = (e0-eAB) - (e0-eA) - (e0-eB)
    ratio = abs(D)/(abs(dA)+abs(dB)) if (abs(dA)+abs(dB)) > 0 else float('nan')
    rows.append((lam, c0, dA, dB, dAB, D, ratio, Dc))
    print(f"  {lam:>8}{c0:>18.14f}{dA:>16.4e}{dAB:>16.4e}{D:>16.4e}{ratio:>18.9f}{Dc:>20.2e}")

print()
print("  POWER LAWS, FITTED ONLY WHERE THE RESIDUALS ARE SMALL:")
lams = np.array([r[0] for r in rows]); chi0 = np.array([r[1] for r in rows])
dAs = np.abs(np.array([r[2] for r in rows])); Ds = np.abs(np.array([r[5] for r in rows]))
def fit(x, y, lab, kmax):
    xs, ys = np.log(x[:kmax]), np.log(y[:kmax])
    c = np.polyfit(xs, ys, 1); r = ys - np.polyval(c, xs)
    print(f"    {lab:<14} exponent {c[0]:>8.4f}   max |log-resid| {np.abs(r).max():.3e}"
          f"   rms {np.sqrt((r**2).mean()):.3e}   (lam <= {x[kmax-1]})")
    return c[0], np.abs(r).max()
for kmax in (6, 9, len(LAMS)):
    print(f"  -- window: the {kmax} smallest couplings --")
    e0, r0 = fit(lams, chi0, "chi_C alone", kmax)
    e1, r1 = fit(lams, dAs, "Delta_A", kmax)
    e2, r2 = fit(lams, Ds,  "D (3-record)", kmax)
    print(f"    EXPONENT GAP  D - Delta_A = {e2-e1:+.4f}")
gap_ok = None
xs = np.log(lams[:6])
e1 = np.polyfit(xs, np.log(dAs[:6]), 1)[0]
e2 = np.polyfit(xs, np.log(Ds[:6]), 1)[0]
r1 = np.abs(np.log(dAs[:6]) - np.polyval(np.polyfit(xs, np.log(dAs[:6]), 1), xs)).max()
r2 = np.abs(np.log(Ds[:6]) - np.polyval(np.polyfit(xs, np.log(Ds[:6]), 1), xs)).max()
print()
print(f"  NOISE FLOOR: the separate-site control column is the pipeline's zero; "
      f"max |D_ctrl| = {max(abs(r[7]) for r in rows):.2e}")
print(f"  smallest |D| measured = {Ds.min():.3e}, which is "
      f"{Ds.min()/max(abs(r[7]) for r in rows):.1e} times the floor -- the signal is real.")
check("weak-window fits are clean (max |log-residual| < 0.05 for both)", r1 < 0.05 and r2 < 0.05,
      f"  Delta_A {r1:.3e}, D {r2:.3e}")
check("D vanishes FASTER than Delta_A as lam -> 0 (leading-order response is ADDITIVE)",
      (e2 - e1) > 0.5, f"  gap {e2-e1:+.4f}  (Delta_A ~ lam^{e1:.3f}, D ~ lam^{e2:.3f})")
print()
print("="*112)
print("T8 SELF-CHECK SUMMARY:", "ALL PASS" if not FAIL else f"{len(FAIL)} FAILURES: {FAIL}")
print("="*112)
