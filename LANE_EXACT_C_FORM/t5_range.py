"""T5 -- WHERE COULD A SEPARATION DEPENDENCE POSSIBLY COME FROM?

Two exact structural facts, each verified numerically beside a positive control that WOULD
have registered the opposite.

  FACT A (from clause (ii) alone).  Every record commutes with H.  If every coupling operator
     is also a record, then the total Hamiltonian
           H_tot = H (x) I + I (x) H_B + lam * sum_a R_a (x) B_{j(a)}
     commutes with every R_a.  Each record's VALUE is an exact constant of motion.  There is
     therefore NO dynamical channel in the carrier through which one record can influence
     another: any record-record interaction must be mediated ENTIRELY by a shared environment
     degree of freedom.  Durability forbids carrier-mediated interaction.  Exact, all N.

  FACT B (from the environment's structure).  The model's Environment has H_B = sum_j e_j Z_j
     -- a sum of SINGLE-qubit terms with NO bath-bath coupling.  Hence the joint evolution
     factorises over bath qubits within each joint record sector, and a record coupled to bath
     qubit j is invisible to a record coupled to bath qubit j' != j.  The record-record
     interaction is therefore a STEP FUNCTION of the bath-site assignment: exactly its alone
     value off-site, a fixed O(1) reduction on-site.  Range zero.  Exact.

  CONSEQUENCE.  Gravity's form requires an influence that FALLS OFF over a distance.  Here
  there is no distance for it to fall off over: the carrier channel is forbidden by clause (ii)
  and the environment channel is a zero-range contact term.  Any separation dependence must be
  IMPORTED as carrier/environment data, exactly as clause (v) must be.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM")
from lane_utils import *

FAIL = []
def check(name, ok, extra=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {extra}", flush=True)
    if not ok: FAIL.append(name)

print(__doc__)

# ------------------------------------------------------------------ FACT A, on the FULL carrier
print("="*112)
print("T5(a)  RECORD VALUES ARE EXACT CONSTANTS OF MOTION -- with a POSITIVE CONTROL that moves them")
print("="*112)
n = 8
Xall = xz_to_matrix([1]*n+[0]*n, n); Zall = xz_to_matrix([0]*n+[1]*n, n)
H = -(Xall + Zall)
M = RecordModel(H)
stab = stab_nn2(n); S, L, _ = derived_logical_span(stab, n)
env = Environment(nq=3, energies=(1.0, 1.4, 0.7), beta=2.0)
Pg, kdim = M.ground_space()
RA = xz_to_matrix(pauli_vec(n,(0,1),()), n)
RB = xz_to_matrix(pauli_vec(n,(),(4,5)), n)
ERR = xz_to_matrix(pauli_vec(n,(),(0,)), n)          # single-qubit Z_0 -- NOT a record (clause ii fails)
print(f"  carrier [[{n},{n-2},2]] FULL space dim {2**n}, bath 3 qubits, ground space dim {kdim}")
print(f"  Z_0 is a record: {is_nontrivial_logical(pauli_vec(n,(),(0,)), S, L, n)}   "
      f"(it must be False -- it is the POSITIVE CONTROL)")
# THE STATE MUST BE POLARISED IN THE RECORDS OR THE TEST IS VACUOUS.  A first version used
# the maximally mixed code state, where <R> = 0 identically and NOTHING can move; the positive
# control failed and caught it.  Here rho0 is the maximally mixed state of the joint
# R_A = +1, R_B = +1 sector, so <R_A> = <R_B> = 1 at t = 0 and motion is detectable.
Q = Pg @ (np.eye(2**n) + RA) @ (np.eye(2**n) + RB) / 4
rho0 = Q / np.real(np.trace(Q))
assert abs(np.real(np.trace(rho0@RA)) - 1) < 1e-9 and abs(np.real(np.trace(rho0@RB)) - 1) < 1e-9
print(f"  initial state: maximally mixed on the code sector R_A=+1, R_B=+1; "
      f"<R_A>(0) = {np.real(np.trace(rho0@RA)):.12f}, <R_B>(0) = {np.real(np.trace(rho0@RB)):.12f}")
print(f"  {'coupling':<28}{'coupling is a record':>22}{'lam':>6}{'t':>5}"
      f"{'|<R_A>(t) - <R_A>(0)|':>24}{'|<R_B>(t)-<R_B>(0)|':>22}")
rows = []
for cname, Cop, isrec in (("R_A = X0X1", RA, True), ("R_B = Z4Z5", RB, True),
                          ("R_A + R_B", RA+RB, True), ("Z_0  (NOT a record)", ERR, False)):
    for lam in (0.4, 1.2):
        for t in (3.0, 9.0):
            r = M.evolve([(Cop, 0)], env, lam=lam, t=t, state0=rho0)
            rS = r.reshape(2**n, env.dim, 2**n, env.dim).trace(axis1=1, axis2=3)
            mA = abs(float(np.real(np.trace(rS@RA))) - float(np.real(np.trace(rho0@RA))))
            mB = abs(float(np.real(np.trace(rS@RB))) - float(np.real(np.trace(rho0@RB))))
            rows.append((cname, isrec, lam, t, mA, mB))
            print(f"  {cname:<28}{str(isrec):>22}{lam:>6}{t:>5}{mA:>24.3e}{mB:>22.3e}")
check("record values EXACTLY conserved under every record coupling",
      all(r[4] < 1e-12 and r[5] < 1e-12 for r in rows if r[1]),
      f"  max motion {max(max(r[4],r[5]) for r in rows if r[1]):.2e}")
check("POSITIVE CONTROL: a non-record coupling DOES move a record value",
      any(max(r[4], r[5]) > 1e-6 for r in rows if not r[1]),
      f"  max motion {max(max(r[4],r[5]) for r in rows if not r[1]):.2e}")

# ------------------------------------------------------------------ FACT B, range of the interaction
print()
print("="*112)
print("T5(b)  RANGE OF THE RECORD-RECORD INTERACTION IN THE ENVIRONMENT")
print("="*112)
nn = 10
reps, idx = code_reps(nn); d = len(reps)
stab = stab_nn2(nn); S, L, _ = derived_logical_span(stab, nn)
def R2(kind, i, j):
    v = pauli_vec(nn,(i,j),()) if kind=="X" else pauli_vec(nn,(),(i,j))
    assert is_nontrivial_logical(v, S, L, nn)
    return (compress_XX(i,j,nn,reps,idx) if kind=="X" else compress_ZZ(i,j,nn,reps,idx)).astype(complex)
MA = R2("X",0,1)
st = np.eye(d, dtype=complex)/d; Hr = -2.0*np.eye(d, dtype=complex)
print(f"  carrier [[{nn},{nn-2},2]] code space dim {d}; bath 2 qubits, uniform energies")
print(f"  A = X0X1 always at bath site 0.  Partner B scanned over CARRIER separation AND")
print(f"  over BATH SITE.  ratio = chi_A(with B) / chi_A(alone).")
env4 = Environment(nq=2, energies=(1.0,)*2, beta=2.0)
for lam in (0.4, 0.8, 1.2):
    alone = chi_avg(Hr, env4, [(MA,0)], lam, [MA], st)[0]
    print(f"\n  lam={lam}   chi_A alone = {alone:.15f}")
    print(f"    {'B':<8}{'carrier sep':>12}{'B bath site':>13}{'chi_A':>20}"
          f"{'ratio':>20}{'|ratio-1|':>14}")
    off, on = [], []
    for (p, q) in ((2,3), (4,5), (6,7), (8,9)):
        MB = R2("X", p, q)
        sepv = p - 1
        for site in (0, 1):
            c = chi_avg(Hr, env4, [(MA,0),(MB,site)], lam, [MA], st)[0]
            ratio = c/alone
            (on if site == 0 else off).append(ratio)
            print(f"    {'X%dX%d'%(p,q):<8}{sepv:>12}{site:>13}{c:>20.15f}"
                  f"{ratio:>20.15f}{abs(ratio-1):>14.2e}")
    print(f"    OFF-SITE  ratios: max |ratio-1| = {max(abs(r-1) for r in off):.3e}"
          f"   (n={len(off)}, every carrier separation)")
    print(f"    ON-SITE   ratios: min {min(on):.12f}  max {max(on):.12f}  "
          f"spread over carrier separation {max(on)-min(on):.3e}")
    check(f"lam={lam}: OFF-SITE interaction is EXACTLY zero (ratio 1) at every separation",
          max(abs(r-1) for r in off) < 1e-13, f"  max |ratio-1| {max(abs(r-1) for r in off):.2e}")
    check(f"lam={lam}: ON-SITE interaction is O(1) and separation-independent",
          (1-max(on)) > 1e-3 and (max(on)-min(on)) < 1e-13,
          f"  1-ratio = {1-np.mean(on):.6f}, spread {max(on)-min(on):.2e}")

print()
print("="*112)
print("T5(c)  UPPER BOUND ON ANY HIDDEN SEPARATION-DEPENDENT TERM")
print("="*112)
print("""  The OFF-SITE column above is exactly 1 for the exact reason given in FACT B, and the
  ON-SITE column is exactly constant for the exact reason given in T2(b).  Both are EXACT
  ARGUMENTS, so the bound on a hidden separation-dependent term in these quantities is
  IDENTICALLY ZERO at every N -- not a numerical bound, and therefore not defeated by the
  weakness objection.

  Where only numbers are available (the m-dependence of the shared-bath defect in T1(c)), the
  bound is the measured noise floor: chi values are O(0.1) and equivalent replicas agree to
  ~1e-15 absolute, i.e. ~1e-14 relative.  That is FOURTEEN orders of magnitude, not thirty-six:
  a gravity-strength residual would be invisible there and the numerical columns CANNOT exclude
  one.  Only the exact columns can, and where they apply they exclude it completely.""")

print()
print("="*112)
print("T5 SELF-CHECK SUMMARY:", "ALL PASS" if not FAIL else f"{len(FAIL)} FAILURES: {FAIL}")
print("="*112)
