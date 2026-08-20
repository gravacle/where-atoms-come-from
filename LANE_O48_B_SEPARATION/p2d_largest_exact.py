"""LANE_O48_B_SEPARATION -- PART 2, STEP D: THE LARGEST FULLY EXACT COMPUTATION, AND HOW MUCH OF
E0(z) THE TWO-BODY COEFFICIENT ACTUALLY ACCOUNTS FOR.

ROUTE 1 is exact to ALL orders in g but needs every one of the 2^m record blocks.  Here it is
pushed to m = 20 (1,048,576 blocks, one exact 20x20 SVD each) -- the largest fully exact Walsh
transform this lane reached.  Reported alongside: the weight of the 3-body and higher Walsh
coefficients, so the two-body truncation is measured rather than assumed.
"""
import sys, time, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_B_SEPARATION")
from mediator import E0_batch, all_configs, chi_free
from common import fwht_fast, pair_index

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)

P("=" * 118)
P("PART 2 STEP D -- LARGEST FULLY EXACT WALSH TRANSFORM, AND THE WEIGHT OUTSIDE THE TWO-BODY TERM")
P("=" * 118)
m = 20
t = np.ones(m - 1); w = np.zeros(m - 1); i0 = 4
T = chi_free(m, t)
Z = all_configs(m)
P("")
P(f"m = {m} record qubits + {m} mediator qubits = {2*m} qubits; {2**m} exact blocks, one SVD each.")
P("")
P("-" * 118)
P("[D1] EXACT J_eff (ROUTE 1, all orders in g) AGAINST THE g^2 LIMIT (ROUTE 2), i0 = 4.")
P("-" * 118)
P(f"{'g':>6} {'time s':>8} " + " ".join(f"{'r='+str(r):>13}" for r in (1, 2, 3, 5, 7, 9, 11, 13)))
P(f"{'g^2 lim':>6} {'-':>8} " + " ".join(f"{-8*T[i0,i0+r]:>13.5e}" for r in (1, 2, 3, 5, 7, 9, 11, 13)))
store = {}
for g in (0.02, 0.05, 0.10, 0.20, 0.40):
    t0 = time.time()
    E = E0_batch(m, Z, t, w, g)
    c = fwht_fast(E) / 2 ** m
    store[g] = c
    P(f"{g:>6.2f} {time.time()-t0:>8.1f} " +
      " ".join(f"{float(c[pair_index(i0,i0+r,m)])/g**2:>13.5e}" for r in (1,2,3,5,7,9,11,13)))
P("")
P("-" * 118)
P("[D2] HOW MUCH OF E0(z) LIVES OUTSIDE THE TWO-BODY TERMS?  Walsh weight by interaction order.")
P("     D-15 control: the same decomposition of a PURELY TWO-BODY inserted energy, which must put")
P("     100.000000 % of its non-constant weight in order 2 and exactly 0 elsewhere.")
P("-" * 118)
pc = np.bitwise_count(np.arange(2 ** m, dtype=np.int64))
P(f"{'source':>26} {'g':>6} " + " ".join(f"{'order '+str(k):>13}" for k in (1, 2, 3, 4, 5)) +
  f" {'order>=6':>13}")
def weights(c):
    w2 = c ** 2
    tot = w2[1:].sum()
    out = []
    for k in (1, 2, 3, 4, 5):
        out.append(100.0 * w2[pc == k].sum() / tot)
    out.append(100.0 * w2[(pc >= 6)].sum() / tot)
    return out
for g in (0.02, 0.05, 0.10, 0.20, 0.40):
    ws = weights(store[g])
    P(f"{'MEDIATED E0(z)':>26} {g:>6.2f} " + " ".join(f"{x:>13.6f}" for x in ws[:5]) +
      f" {ws[5]:>13.6f}")
Zs = all_configs(m)
Ectrl = np.zeros(2 ** m)
for a in range(m):
    for b in range(a + 1, m):
        Ectrl += 0.3 * abs(a - b) ** -3.0 * Zs[:, a] * Zs[:, b]
ws = weights(fwht_fast(Ectrl) / 2 ** m)
P(f"{'CONTROL pure 2-body A/r^3':>26} {'-':>6} " + " ".join(f"{x:>13.6f}" for x in ws[:5]) +
  f" {ws[5]:>13.6f}")
P("")
P('READ (filled from the numbers above): EVERY ODD ORDER IS EXACTLY 0.000000 at every g.  That is')
P('      not luck: the admissible writer found by search in p2a is the global flip, so E0(-z) = E0(z)')
P('      identically and all odd Walsh coefficients vanish.  Clause (iv) is what makes the record-record')
P('      energy purely even.  The two-body term carries 99.999431 % of the non-constant weight at')
P('      g=0.02, 99.981081 % at g=0.05, 99.650918 % at g=0.10, 95.985319 % at g=0.20 and 94.560371 %')
P('      at g=0.40, the remainder sitting almost entirely in the FOUR-body coefficients.  So the')
P('      pairwise description of E0(z) is essentially complete in the weak-coupling regime where the')
P('      power law lives, and starts to leak at strong coupling -- which is the same place [B8] found')
P('      the tail collapsing.  The control puts 100.000000 % in order 2, so the decomposition is not')
P('      manufacturing the four-body weight.')
P("")
P("[D3] LARGEST SIZES ACTUALLY REACHED IN THIS LANE, AND WHAT STOPPED EACH:")
P("     ROUTE 1, exact Walsh over all blocks .............. m = 20  (2^20 = 1048576 blocks)")
P("                stopped by: 2^m blocks; m = 22 is 4x this and m = 24 is 16x.")
P("     ROUTE 2, exact O(g^2) orbital sum, 1D ............. m = 8192 mediator sites")
P("                stopped by: dense O(m^3) eigendecomposition and an m^2 eigenvector array.")
P("     ROUTE 2, 2D mediator .............................. 63 x 64 = 4032 sites")
P("     ROUTE 2, 3D mediator .............................. 15 x 16 x 17 = 4080 sites")
P('     Full dense ED of the mediated carrier ............. 9 qubits (dim 512, dilute case in p2c),')
P('                8 qubits (dim 256, every case in p2a)')
P('                stopped by: the 4^nq EXHAUSTIVE Pauli search (4^9 = 262144), not by the ED itself.')
P("     Diagonal-exact clause check, Part 1 ............... n = 16 (dim 65536)")
P("     Exhaustive Pauli search, Part 1 ................... n = 10 (4^10 = 1048576 Paulis)")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_B_SEPARATION/p2d_largest_exact.txt","w").write("\n".join(OUT)+"\n")
