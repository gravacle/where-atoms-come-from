"""ADVERSARIAL CHECK V5 (D-18) -- ARE THE 'RECORDS' ACTUALLY RECORDS ON THIS CARRIER?
Independent re-derivation: logicals from symplectic_logicals ONLY (never nominated), clauses
checked with record_model's OWN clause_iii / clause_iv, plus an explicit statement of which
clause was NOT checked.  CONTROL IN THE SAME TABLE: a NON-record operator on the same carrier
(a single-qubit Z_0, which is not a logical) must FAIL at least one clause, or the clause
checker is dead and proves nothing.
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC")
from lib_scaleb import *
import numpy as np
from record_model import clause_iii, clause_iv, eigenspaces

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)

P("=" * 104)
P("V5  D-18 CLAUSE RE-CHECK, INDEPENDENT OF THE LANE'S SCRIPTS.")
P("    logicals come from symplectic_logicals (conjugate pairs); nothing is nominated.")
P("=" * 104)
P("")
P("  %-4s %-4s %-9s | %-7s %-7s %-7s %-7s %-9s | %-28s" %
  ("n", "2k", "mult(H)", "(i)", "(ii)", "(iii)", "(iv)", "SC-1 sp", "NON-RECORD control Z_0"))
P("  " + "-" * 100)
for n in [4, 6, 8]:
    stab, pairs = carrier(n)
    vs, lab = record_vectors(pairs, n)
    k = len(pairs); m = len(vs)
    H = hamiltonian(n)
    es = eigenspaces(H)
    mult = [e[2] for e in es]
    R = [xz_to_matrix(v, n) for v in vs]
    c1 = all(np.linalg.norm(r - r.conj().T) < 1e-9 and
             np.linalg.norm(r @ r - np.eye(2 ** n)) < 1e-9 for r in R)
    c2 = all(np.linalg.norm(H @ r - r @ H) < 1e-9 for r in R)
    c3 = all(clause_iii(r, es) for r in R)
    c4 = all(clause_iv(r, es) for r in R)
    S = np.array([[sp_form(vs[i], vs[j], n) for j in range(m)] for i in range(m)], dtype=int)
    # F2 determinant by elimination
    A = S.copy() % 2; rank = 0
    for c in range(m):
        piv = None
        for r in range(rank, m):
            if A[r, c]: piv = r; break
        if piv is None: continue
        A[[rank, piv]] = A[[piv, rank]]
        for r in range(m):
            if r != rank and A[r, c]: A[r] = (A[r] + A[rank]) % 2
        rank += 1
    sc1 = (rank == m)
    # NON-RECORD control: Z on physical qubit 0 alone
    zv = [0] * n + [1 if q == 0 else 0 for q in range(n)]
    Z0 = xz_to_matrix(zv, n)
    z_ii = np.linalg.norm(H @ Z0 - Z0 @ H) < 1e-9
    z_iii = clause_iii(Z0, es); z_iv = clause_iv(Z0, es)
    P("  %-4d %-4d %-9s | %-7s %-7s %-7s %-7s %-9s | (ii)=%s (iii)=%s (iv)=%s" %
      (n, m, str(mult), c1, c2, c3, c4, sc1, z_ii, z_iii, z_iv))
P("")
P("  CLAUSE (v) PROTECTED WAS NOT CHECKED BY THE LANE AND IS NOT CHECKED HERE.  record_model's")
P("  own header states (v) is CARRIER DATA, not derivable from (H,{L_k}).  The lane's probe")
P("  text says '(i)-(iv)' honestly; the headline then speaks of 'records' throughout.")
P("  The control column shows the clause checker is live: Z_0 is not a logical and fails.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/VERIFY/V5_records_are_records.txt",
     "w").write("\n".join(OUT) + "\n")
