"""O-5 A.  THE CARRIER, AND WHAT EXACT CLAUSE (ii) ACTUALLY COVERS.

Before asking what the RIGHT approximate clause (ii) is, establish two things with numbers:

  1. the carrier is what we say it is  (dim 256, ground degeneracy 4, distance 2, logicals verified);
  2. what exact clause (ii) already excludes.

On (2) there is a specific thing to check that the program has not checked.  P1's verification of the
five clauses used jump operators L_k = Z_l -- Z-type local noise.  The record R is a Z-type logical,
so [Z_l, R] = 0 identically and clause (ii) passes for free.  What happens with X-type local noise,
which is the OTHER half of a generic local bath?

SELF-CHECKS, each with an independently known answer:
  S1  the 8 stabilizers commute pairwise                                  -> 0
  S2  spectrum of H is contained in {-8,-6,...,8}                         -> known
  S3  ground energy -8 with degeneracy 4 = 2^{2g}, g=1                    -> known
  S4  logicals commute with H; Z_A/X_A and Z_B/X_B anticommute; cross pairs commute
  S5  code distance: no weight-1 Pauli acts non-trivially inside the code space, some weight-2 does
POSITIVE CONTROL for every reported zero: an operator that WOULD register if a nonzero existed.
"""
import numpy as np
from o5_common import (DIM, NQ, STARS, PLAQS, Zop, Xop, Yop, toric_H, sym_H,
                       Z_A_SUP, X_A_SUP, Z_B_SUP, X_B_SUP)

np.set_printoptions(precision=4, suppress=True)
print("=" * 100)
print("O-5 A.  CARRIER + WHAT EXACT CLAUSE (ii) COVERS")
print("=" * 100)

Av = [Xop(s) for s in STARS]
Bp = [Zop(p) for p in PLAQS]
H = toric_H()
Z_A, X_A, Z_B, X_B = Zop(Z_A_SUP), Xop(X_A_SUP), Zop(Z_B_SUP), Xop(X_B_SUP)

print(f"\n  lattice 2x2 torus, {NQ} edges, dim = {DIM}")
print(f"  stars      {STARS}")
print(f"  plaquettes {PLAQS}")
print(f"  RECORD  R = Z on {Z_A_SUP}        WRITER W = X on {X_A_SUP}")

# ---- S1 -----------------------------------------------------------------------------------------
S = Av + Bp
c = max(np.linalg.norm(a @ b - b @ a, 2) for a in S for b in S)
print(f"\n  S1  max ||[S_i,S_j]|| over the 8 stabilizers = {c:.2e}   {'PASS' if c < 1e-12 else 'FAIL'}")
# positive control: two stabilizers replaced by an anticommuting pair must register
ctrl = np.linalg.norm(Xop([0]) @ Zop([0]) - Zop([0]) @ Xop([0]), 2)
print(f"      POSITIVE CONTROL ||[X_0,Z_0]||               = {ctrl:.2e}   "
      f"{'PASS (the measure registers when there IS non-commutation)' if ctrl > 1 else 'FAIL'}")

# ---- S2, S3 -------------------------------------------------------------------------------------
ev = np.linalg.eigvalsh(H)
uniq = np.unique(np.round(ev, 8))
ok2 = set(uniq.tolist()) <= set(float(x) for x in range(-8, 9, 2))
g0 = int(np.sum(np.abs(ev - ev[0]) < 1e-9))
print(f"\n  S2  distinct eigenvalues of H = {uniq}")
print(f"      contained in {{-8,-6,...,8}}                   {'PASS' if ok2 else 'FAIL'}")
print(f"  S3  ground energy {ev[0]:+.6f}, degeneracy {g0}   (expected -8, 4 = 2^(2g), g=1)   "
      f"{'PASS' if abs(ev[0] + 8) < 1e-9 and g0 == 4 else 'FAIL'}")
print(f"      gap to first excited level = {ev[g0] - ev[0]:.6f}   (expected 4)")

# ---- S4 -----------------------------------------------------------------------------------------
def com(A, B): return np.linalg.norm(A @ B - B @ A, 2)
def acom(A, B): return np.linalg.norm(A @ B + B @ A, 2)
print("\n  S4  logical algebra")
for nm, O in [("Z_A", Z_A), ("X_A", X_A), ("Z_B", Z_B), ("X_B", X_B)]:
    print(f"      ||[H,{nm}]|| = {com(H, O):.2e}   {'PASS' if com(H, O) < 1e-12 else 'FAIL'}")
pairs = [("Z_A", "X_A", "anti"), ("Z_B", "X_B", "anti"),
         ("Z_A", "X_B", "comm"), ("Z_B", "X_A", "comm"), ("Z_A", "Z_B", "comm")]
d = dict(Z_A=Z_A, X_A=X_A, Z_B=Z_B, X_B=X_B)
for a, b, kind in pairs:
    x = acom(d[a], d[b]) if kind == "anti" else com(d[a], d[b])
    print(f"      {a},{b} {kind:4s}: ||{'{,}' if kind == 'anti' else '[,]'}|| = {x:.2e}   "
          f"{'PASS' if x < 1e-12 else 'FAIL'}")

# ---- S5  distance -------------------------------------------------------------------------------
e, U = np.linalg.eigh(H)
P = U[:, :4] @ U[:, :4].conj().T          # code-space projector
def code_action(O):
    """how far P O P is from a multiple of the identity ON the code space"""
    M = U[:, :4].conj().T @ O @ U[:, :4]
    return float(np.linalg.norm(M - np.trace(M) / 4 * np.eye(4), 2))
w1 = []
for l in range(NQ):
    for nm, O in [("X", Xop([l])), ("Y", Yop(l)), ("Z", Zop([l]))]:
        w1.append((code_action(O), f"{nm}_{l}"))
w1.sort(reverse=True)
best1 = w1[0]
print(f"\n  S5  DISTANCE.  max over ALL 24 weight-1 Paulis of ||P O P - scalar|| = {best1[0]:.2e} ({best1[1]})")
print(f"      {'PASS -- no weight-1 operator acts non-trivially in the code space' if best1[0] < 1e-12 else 'FAIL'}")
print(f"      POSITIVE CONTROL, weight-2 logical Z_A: ||P Z_A P - scalar|| = {code_action(Z_A):.4f}   "
      f"{'PASS (a weight-2 operator DOES)' if code_action(Z_A) > 0.5 else 'FAIL'}")
print("      => code distance d = 2.")

# =================================================================================================
print("\n" + "=" * 100)
print("  WHAT EXACT CLAUSE (ii) COVERS:  [L_k, R] = 0 DEPENDS ENTIRELY ON THE NOISE MODEL")
print("=" * 100)
print("""
  P1 verified the five clauses with L_k = Z_l.  R is a Z-type logical, so those commutators vanish
  identically.  Here are BOTH halves of a generic local bath, on the same record.
""")
print(f"    {'jump operator family':38s} {'max_k ||[L_k,R]||':>20s}   clause (ii)")
for nm, fam in [("L_k = Z_l   (dephasing, P1's choice)", [Zop([l]) for l in range(NQ)]),
                ("L_k = X_l   (bit flip)", [Xop([l]) for l in range(NQ)]),
                ("L_k = Y_l", [Yop(l) for l in range(NQ)]),
                ("L_k = {X_l, Y_l, Z_l}  (generic local)",
                 [f([l]) if f is not Yop else Yop(l) for l in range(NQ) for f in (Xop, Yop, Zop)])]:
    m = max(np.linalg.norm(L @ Z_A - Z_A @ L, 2) for L in fam)
    print(f"    {nm:38s} {m:20.4f}   {'HOLDS EXACTLY' if m < 1e-12 else 'FAILS -- and by O(1), not by epsilon'}")

print("""
  READ THIS PLAINLY.  Under a generic local bath the exact clause (ii) does not fail by a small
  amount for the program's own record.  It fails MAXIMALLY: ||[X_l,R]|| = 2 = 2||R||, the largest
  value a commutator of two involutions can take.  Exact clause (ii) as stated is satisfied by the
  toric-code record only for a bath chosen to commute with it.

  THIS IS THE REAL CONTENT OF O-5, AND IT IS NOT AN EPSILON PROBLEM.  Any relaxation of the form
  ||[L_k,R]|| <= epsilon must take epsilon >= 2 to admit the carrier's own record -- and epsilon = 2
  is no constraint at all, since ||[L,R]|| <= 2 always for involutions.  Quantified below.
""")
mx = 0.0
rng = np.random.default_rng(11)
for _ in range(200):
    A = rng.normal(size=(DIM, DIM)) + 1j * rng.normal(size=(DIM, DIM))
    A = (A + A.conj().T) / 2
    A = A / np.linalg.norm(A, 2)
    mx = max(mx, np.linalg.norm(A @ Z_A - Z_A @ A, 2))
print(f"    sup over 200 random unit-norm Hermitian L of ||[L,R]||   = {mx:.4f}   (bound is 2)")
print(f"    => the tolerance epsilon = 2 admits EVERY operator.  The relaxed L-clause is vacuous at")
print(f"       the tolerance needed to admit the program's own record.")

# --- the same statement for the symmetry carrier, for later contrast ------------------------------
Hs = sym_H()
Rs = Zop([0])
es = np.linalg.eigvalsh(Hs)
gs = int(np.sum(np.abs(es - es[0]) < 1e-9))
print("\n" + "-" * 100)
print("  SYMMETRY CARRIER (same 256-dim space, H_sym = -sum_{l>=1} X_l)")
print(f"    ground energy {es[0]:+.4f}, degeneracy {gs}  (expected -7, 2)   "
      f"{'PASS' if abs(es[0] + 7) < 1e-9 and gs == 2 else 'FAIL'}")
print(f"    gap = {es[gs] - es[0]:.4f}  (expected 2)")
print(f"    record R_sym = Z_0 :  ||[H_sym,R_sym]|| = {com(Hs, Rs):.2e}   "
      f"{'PASS (exact record of H_sym)' if com(Hs, Rs) < 1e-12 else 'FAIL'}")
Usx = np.linalg.eigh(Hs)[1]
Ms = Usx[:, :gs].conj().T @ Rs @ Usx[:, :gs]
nsc = np.linalg.norm(Ms - np.trace(Ms) / gs * np.eye(gs), 2)
print(f"    non-trivial on the ground space: ||PR_symP - scalar|| = {nsc:.4f}   "
      f"{'PASS (clause iii holds)' if nsc > 0.5 else 'FAIL'}")
print(f"    ||[X_0, R_sym]|| = {np.linalg.norm(Xop([0]) @ Rs - Rs @ Xop([0]), 2):.4f}  -- IDENTICAL to")
print(f"       the topological carrier's ||[X_l,R]|| = 2.  Remember this number.  A local X flips")
print(f"       R_sym at zero energy cost here, and at cost 4 on the toric code, and the commutator")
print(f"       norm CANNOT TELL THE TWO APART.  Measured in o5_d_theorems.py.")
print("\n  A DONE")
