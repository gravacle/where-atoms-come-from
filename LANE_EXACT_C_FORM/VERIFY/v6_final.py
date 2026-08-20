"""V6 -- corrected exact F_2 recount, plus the anticommuting-partner number from a 1-qubit model.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM")
from lane_utils import *

print("="*100)
print("V6(a)  EXACT F_2 RECOUNT OVER m DISJOINT [[4,2,2]] BLOCKS (corrected labelling)")
print("="*100)
print(f"  {'m':>4}{'n=4m':>6}{'Qk = #logical qubits':>22}{'F_2 rank of full Gram':>24}{'defect vs m*Qk(1)':>20}")
ok = True
for m in (1,2,3,4,6,8,10,12):
    n = 4*m
    S, L, pairs = derived_logical_span(stab_blocks(m,4), n)
    k = len(pairs)
    gens = [p for pr in pairs for p in pr]
    r = f2_rank([[sp(x,y,n) for y in gens] for x in gens], 2*k)
    print(f"  {m:>4}{n:>6}{k:>22}{r:>24}{k - m*2:>20}")
    ok &= (k == 2*m) and (r == 2*k)
print(f"  Qk = 2m and the Gram matrix is FULL RANK 2k at every m: {ok}   -> lane CONFIRMED exactly")
print("  But: the stabiliser group of a direct sum of blocks IS the direct sum of the groups, so")
print("  k, rank and log2(protected dim) are additive BY CONSTRUCTION.  Confirming them is a")
print("  consistency check on the code, not evidence that a physical source is additive.")

print()
print("="*100)
print("V6(b)  THE 'ANISOTROPY' NUMBER 0.254089251974 FROM A ONE-QUBIT SYSTEM, NO CODE, NO RECORDS")
print("="*100)
Xp = np.array([[0,1],[1,0]], complex); Zp = np.array([[1,0],[0,-1]], complex)
env = Environment(nq=2, energies=(1.0,)*2, beta=2.0)
st = np.eye(2, dtype=complex)/2; Hr = -2.0*np.eye(2, dtype=complex)
for lam, ref in ((0.4, 0.196063), (0.8, 0.254089251974), (1.2, None)):
    v = chi_avg(Hr, env, [(Zp,0),(Xp,0)], lam, [Zp], st)[0]
    print(f"  lam={lam}: readout Z, partner X (they ANTICOMMUTE, sp=1), same bath site: "
          f"chi = {v:.12f}   lane sp=1 value: {ref}")
print("  and the commuting case, for the same one-qubit reduction:")
Z1 = np.kron(Zp, np.eye(2)); Z2 = np.kron(np.eye(2), Zp)
st4 = np.eye(4, dtype=complex)/4; Hr4 = -2.0*np.eye(4, dtype=complex)
print(f"  lam=0.8: two COMMUTING +-1 operators on one bath site: "
      f"chi = {chi_avg(Hr4, env, [(Z1,0),(Z2,0)], 0.8, [Z1], st4)[0]:.12f}   lane sp=0 value: 0.136408688972")
print("  => table G's 'ANISOTROPIC IN THE F_2 PAIRING, ISOTROPIC IN GEOMETRY' is the statement")
print("     that two commuting +-1 observables share a channel differently from two")
print("     anticommuting ones.  One qubit reproduces it.  No carrier, no code, no geometry.")
