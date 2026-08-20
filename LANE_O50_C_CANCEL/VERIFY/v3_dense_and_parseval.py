"""
VERIFY 3.
 (A) Reproduce the lane's DENSE L=2 clause numbers and its step-6 zeros independently.
 (B) Test the lane's escape #2 ("non-linear functionals CLOSED EXACTLY by Parseval").
     Parseval bounds nothing about the SIZE of Var relative to m.  If a functional exists whose
     RESPONSIVE part is EXTENSIVE at a typical configuration, the lane's clause [C5](c)
     ("the responsive part is O(sqrt(m)) ... SUB-EXTENSIVE") is false outside linear functionals.
"""
import sys, math
from fractions import Fraction
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
import numpy as np
from record_model import (RecordModel, symplectic_logicals, xz_to_matrix, eigenspaces,
                          clause_iii, clause_iv)

def line(s=""): print(s, flush=True)

L = 2; n = 2 * L * L
def h(i, j): return (i % L) * L + (j % L)
def v(i, j): return L * L + (i % L) * L + (j % L)
S = []
for i in range(L):
    for j in range(L):
        r = [0] * (2 * n)
        for e in (h(i, j), h(i, j - 1), v(i, j), v(i - 1, j)): r[e] ^= 1
        S.append(r)
        r = [0] * (2 * n)
        for e in (h(i, j), h(i + 1, j), v(i, j), v(i, j + 1)): r[n + e] ^= 1
        S.append(r)
H = np.zeros((256, 256), dtype=complex)
for row in S: H -= xz_to_matrix(row, n)

line("=" * 100)
line("V3-A.  DENSE L=2 TORIC CARRIER, REBUILT.  THE LANE'S CLAUSE TABLE, INDEPENDENTLY.")
line("=" * 100)
es = eigenspaces(H)
w = np.linalg.eigvalsh(H)
vals, counts = np.unique(np.round(w, 6), return_counts=True)
line(f"  dim 256, eigenvalue multiplicities = {list(counts)}   (lane: [4, 48, 152, 48, 4])")
line(f"  ground energy = {vals[0]:.6f}  ground-space dim = {counts[0]}   (lane: -8.0, 4)")
pairs = symplectic_logicals(S, n)
line(f"  conjugate pairs returned = {len(pairs)}   (lane: 2)")
line()
line(f"  {'operator':<26} {'||[R,H]||':>12} {'||R^2-I||':>12} {'(iii)':>7} {'(iv)':>7} {'max|Tr(P_E R)|':>16}")
ops = []
for idx, pr in enumerate(pairs):
    for half, name in ((0, 'Zbar'), (1, 'Xbar')):
        M = xz_to_matrix(pr[half], n)
        ops.append((f"{name}_{idx} (computed)", M))
zc = [0] * (2 * n); zc[0] = 0; zc[n + 0] = 1
ops.append(("CONTROL single-qubit Z_0", xz_to_matrix(zc, n)))
ops.append(("CONTROL identity I", np.eye(256, dtype=complex)))
bp = [0] * (2 * n)
for e in (h(0, 0), h(1, 0), v(0, 0), v(0, 1)): bp[n + e] ^= 1
ops.append(("CONTROL stabiliser B_p", xz_to_matrix(bp, n)))
for name, M in ops:
    cH = np.linalg.norm(M @ H - H @ M)
    sq = np.linalg.norm(M @ M - np.eye(256))
    mx = max(abs(np.trace(P @ M)) for _, P, _ in es)
    line(f"  {name:<26} {cH:>12.2e} {sq:>12.2e} {str(clause_iii(M, es)):>7} "
         f"{str(clause_iv(M, es)):>7} {mx:>16.2e}")
line("  (lane: computed logicals 0.00e+00 / 0.00e+00 / True / True / 6.74e-16..8.92e-15;")
line("   Z_0 fails (ii) at 4.53e+01; identity 1.52e+02; B_p 2.40e+01)")
line()
R1 = xz_to_matrix(pairs[0][0], n); R2 = xz_to_matrix(pairs[1][0], n)
rm = RecordModel(H)
fam = rm.commuting_family([R1, R2])
line(f"  commuting_family -> {len(fam)} members ; independently_writable -> "
     f"{rm.independently_writable(fam)}   (lane: 2 members, [0,1])")
line()
line("  STEP-6 T-1 REPRODUCED: <R_1> in Gibbs states (must be 0 by clause (iv)).")
for beta in (0.0, 0.25, 1.0, 4.0, 40.0):
    ev, U = np.linalg.eigh(H)
    p = np.exp(-beta * (ev - ev.min())); p /= p.sum()
    rho = (U * p) @ U.conj().T
    line(f"    beta = {beta:>6}:  <R_1> = {np.real(np.trace(rho @ R1)):>10.2e}   "
         f"<R_2> = {np.real(np.trace(rho @ R2)):>10.2e}")
line("  ADVERSARIAL STATE the lane did NOT test: rho ~ e^{-bH}(1+R_1)/2 -- STATIONARY ([rho,H]=0),")
line("  perfectly admissible as a state, and NOT of the form f(H):")
ev, U = np.linalg.eigh(H)
p = np.exp(-4.0 * (ev - ev.min())); p /= p.sum()
rho0 = (U * p) @ U.conj().T
proj = (np.eye(256) + R1) / 2
rho1 = proj @ rho0 @ proj
rho1 = rho1 / np.trace(rho1)
line(f"    <R_1> = {np.real(np.trace(rho1 @ R1)):.6f} ,  ||[rho,H]|| = {np.linalg.norm(rho1@H-H@rho1):.2e}")
line("    -> a stationary state with record magnetisation 1.  T-1's zero is a zero of f(H) states")
line("       ONLY; it is clause (iv) restated, not an exclusion of ordered stationary states.")

# ---------------------------------------------------------------- B. Parseval does not close it
line()
line("=" * 100)
line("V3-B.  IS ESCAPE #2 ('non-linear functionals, CLOSED EXACTLY by Parseval') ACTUALLY CLOSED?")
line("=" * 100)
line("  Lane's [C5](c): 'the responsive part is O(sqrt(m)) at a typical configuration -- SUB-EXTENSIVE'.")
line("  Counter-functional, entirely inside the record configuration space:  f(s) = (sum_i s_i)^2 .")
line(f"  {'m':>7} {'MEAN=fhat(0)':>13} {'sqrt(VAR) exact':>18} {'sqrt(VAR)/m':>12} {'linear sqrt(VAR)/m':>19}")
for m in (16, 64, 256, 1024, 4096):
    mean = m                                     # E[(sum s)^2] = m
    var = 2 * m * (m - 1)                        # Var[(sum s)^2] = 2m(m-1)
    line(f"  {m:>7} {mean:>13} {math.sqrt(var):>18.4f} {math.sqrt(var)/m:>12.6f} "
         f"{math.sqrt(m)/m:>19.6f}")
line("  BRUTE-FORCE CHECK of the two moments by full enumeration of all 2^m configurations:")
for m in (8, 12, 14):
    tot = 0; tot2 = 0
    for x in range(1 << m):
        s = sum(1 - 2 * ((x >> i) & 1) for i in range(m))
        f = s * s
        tot += f; tot2 += f * f
    mean = Fraction(tot, 1 << m); var = Fraction(tot2, 1 << m) - mean * mean
    line(f"    m={m:>3}  mean = {mean}  (predicted {m})   Var = {var}  (predicted {2*m*(m-1)})  "
         f"match={mean == m and var == 2*m*(m-1)}")
line()
line("  READ: for f = (sum s)^2 the RECORD-DEPENDENT part has sqrt(VAR) = sqrt(2m(m-1)) ~ 1.41*m --")
line("  EXTENSIVE, not O(sqrt(m)).  Parseval gives MEAN = fhat(empty) and VAR = sum_{S!=0} fhat(S)^2")
line("  but places NO bound on VAR relative to m.  So the lane's Cauchy-Schwarz sub-extensivity")
line("  argument is a statement about LINEAR functionals with equal weights ONLY, and escape #2 is")
line("  NOT closed by Parseval.  What still kills f is sign-definiteness of the RESPONSE -- and that")
line("  is the bijection identity, which holds for any permutation of any finite set (see v1-F).")
line("  CONTROL: the response of (sum s)^2 to a write of record 1, summed over all configurations:")
for m in (8, 12):
    tot = 0
    for x in range(1 << m):
        s = sum(1 - 2 * ((x >> i) & 1) for i in range(m))
        y = x ^ 1
        s2 = sum(1 - 2 * ((y >> i) & 1) for i in range(m))
        tot += s2 * s2 - s * s
    line(f"    m={m}: sum over configurations of delta f = {tot}  (0 => response takes both signs)")
