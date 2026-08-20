"""V3 -- (a) IS THE SEPARATION NULL A STATEMENT ABOUT RECORDS, OR JUST THE CARRIER'S SYMMETRY?
        (b) IS 'chi DEPENDS ON THE GRAM MATRIX ALONE' TRUE OFF THE CODE SPACE?

(a)  t7_assumptions.py states the exact fact itself: H = -(X^n + Z^n) is invariant under EVERY
     qubit permutation, so any state that is a function of H is too.  A permutation that fixes
     supp(A) and carries supp(B) to supp(B') then conjugates the whole joint evolution and
     leaves chi EXACTLY invariant.  That argument never mentions records.  So the null must be
     equally true for operators that are NOT records.  Tested here.

(b)  The lane's load-bearing theorem is proved on the CODE SPACE, where Hs = -2*I and rho = I/d
     and two stabiliser-equivalent representatives compress to the SAME matrix.  Off the code
     space they do not.  X0X1 and X2X3X4X5X6X7 = (X^8)*(X0X1) are the SAME logical class, same
     symplectic Gram matrix with every partner, same bath site -- but different weight and
     different support size.  If chi differs between them on the full carrier, then "chi is a
     function of the Gram matrix and the bath-site assignment ALONE" is FALSE outside the code
     space, and the lane's no-metric-form conclusion is scoped to a venue in which geometry was
     projected out before the first float was computed.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM")
from lane_utils import *

n = 8; N = 2**n
Xall = xz_to_matrix([1]*n+[0]*n, n); Zall = xz_to_matrix([0]*n+[1]*n, n)
H = -(Xall + Zall)
stab = stab_nn2(n); S, L, _ = derived_logical_span(stab, n)
w, V = np.linalg.eigh(H)
def thermal(beta):
    p = np.exp(-beta*(w - w.min())); p /= p.sum()
    return (V*p) @ V.conj().T
env = Environment(nq=2, energies=(1.0,)*2, beta=2.0)

def perm_matrix(p):
    P = np.zeros((N, N))
    for v in range(N):
        u = 0
        for i in range(n):
            if (v >> (n-1-i)) & 1: u |= 1 << (n-1-p[i])
        P[u, v] = 1
    return P

print("="*110)
print("V3(a)  THE NULL IS A SYMMETRY IDENTITY, AND IT DOES NOT MENTION RECORDS")
print("="*110)
pi = list(range(n)); pi[2], pi[6] = 6, 2; pi[3], pi[7] = 7, 3
P = perm_matrix(pi)
st2 = thermal(2.0)
A = xz_to_matrix(pauli_vec(n,(0,1),()), n)
B1 = xz_to_matrix(pauli_vec(n,(2,3),()), n)     # sep 1
B2 = xz_to_matrix(pauli_vec(n,(6,7),()), n)     # sep 5
print(f"  permutation pi = (2 6)(3 7)")
print(f"    || P H P^T - H ||           = {np.linalg.norm(P@H@P.T - H):.3e}")
print(f"    || P rho P^T - rho ||       = {np.linalg.norm(P@st2@P.T - st2):.3e}   (beta = 2 thermal)")
print(f"    || P A P^T - A ||           = {np.linalg.norm(P@A@P.T - A):.3e}   (A = X0X1 is FIXED)")
print(f"    || P B(sep1) P^T - B(sep5)||= {np.linalg.norm(P@B1@P.T - B2):.3e}")
print("  => chi(A with B at sep 1) and chi(A with B at sep 5) are the SAME NUMBER by conjugation.")
print("     No clause of the record definition is used anywhere in that sentence.")
print()
print("  CONSEQUENCE, TESTED: the same flatness for a pair that are NOT records at all.")
print("  X_p is X-type of ODD weight, so it ANTICOMMUTES with Z^(x)8: it fails clause (ii).")
Anr = xz_to_matrix(pauli_vec(n,(0,),()), n)
print(f"    X_0 is a derived record? {is_nontrivial_logical(pauli_vec(n,(0,),()), S, L, n)}")
print(f"  {'partner':<8}{'sep':>5}{'chi alone':>18}{'chi crowded':>18}{'interaction':>18}")
Inr = []
a0 = chi_avg(H, env, [(Anr,0)], 0.8, [Anr], st2)[0]
for p in range(2, n):
    Bnr = xz_to_matrix(pauli_vec(n,(p,),()), n)
    c = chi_avg(H, env, [(Anr,0),(Bnr,0)], 0.8, [Anr], st2)[0]
    Inr.append(a0-c)
    print(f"  {'X%d'%p:<8}{p:>5}{a0:>18.12f}{c:>18.12f}{a0-c:>18.12f}")
print(f"  SPREAD over separation for the NON-RECORD pair: {max(Inr)-min(Inr):.3e}")
print("  If this spread is at the floor, 'record-record interactions are separation-blind' is")
print("  not a fact about records: EVERYTHING in this carrier is separation-blind, by symmetry.")

print()
print("="*110)
print("V3(b)  SAME LOGICAL CLASS, SAME GRAM MATRIX, DIFFERENT WEIGHT -- ON THE FULL CARRIER")
print("="*110)
vA  = pauli_vec(n,(0,1),())
vAp = [(x+y)%2 for x, y in zip(vA, [1]*n+[0]*n)]          # X^8 * X0X1 = X2X3X4X5X6X7
print(f"  A  = X0X1                supp size 2  derived record? {is_nontrivial_logical(vA, S, L, n)}")
print(f"  A' = X2X3X4X5X6X7        supp size 6  derived record? {is_nontrivial_logical(vAp, S, L, n)}")
print(f"  same logical class (A+A' in stabiliser span)? "
      f"{in_span([(x+y)%2 for x,y in zip(vA,vAp)], S, 2*n)}")
Ap = xz_to_matrix(vAp, n)
print(f"  {'beta':>6}{'chi(A) alone':>22}{'chi(A_prime) alone':>22}{'difference':>18}")
for beta in (0.0, 0.5, 2.0):
    st = thermal(beta)
    ca  = chi_avg(H, env, [(A,0)],  0.8, [A],  st)[0]
    cap = chi_avg(H, env, [(Ap,0)], 0.8, [Ap], st)[0]
    print(f"  {beta:>6}{ca:>22.12f}{cap:>22.12f}{ca-cap:>18.3e}")
print()
print("  and the two-record interaction with the partner's representative changed:")
print(f"  {'beta':>6}{'I with B=X2X3 (wt2)':>24}{'I with B=X0X1X4..X7 (wt6)':>28}{'difference':>16}")
vB  = pauli_vec(n,(2,3),())
vBp = [(x+y)%2 for x, y in zip(vB, [1]*n+[0]*n)]          # X^8 * X2X3 = X0X1X4X5X6X7
Bw2 = xz_to_matrix(vB, n); Bw6 = xz_to_matrix(vBp, n)
for beta in (0.0, 0.5, 2.0):
    st = thermal(beta)
    al = chi_avg(H, env, [(A,0)], 0.8, [A], st)[0]
    c2 = chi_avg(H, env, [(A,0),(Bw2,0)], 0.8, [A], st)[0]
    c6 = chi_avg(H, env, [(A,0),(Bw6,0)], 0.8, [A], st)[0]
    print(f"  {beta:>6}{al-c2:>24.12f}{al-c6:>28.12f}{(al-c2)-(al-c6):>16.3e}")
print()
print("  A NON-ZERO difference here refutes 'chi is a function of the Gram matrix ALONE' as a")
print("  statement about the carrier; it would hold only after the code-space reduction.")
