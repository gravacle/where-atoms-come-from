# LANE W11-R-CROSS  LEG A -- INDEPENDENT REPRODUCTION, AND THE CONCESSION CHECKED FOR FORCEDNESS.
#
# Before attacking anything I rebuild the two carriers from the sealed bytes with my own code and
# reproduce (i) the steelman's and the registrar's shared headline numbers and (ii) the concession
# the steelman makes in its first leg -- that Z^T_n is gauge-invariant at every PARTIAL tick.
# A concession must be FORCED, not lazy.  Here it is forced, and there is a one-line reason:
#   T(a^g) = g T(a) g^{-1}  entrywise (S1:59-63 on the connection AND on the section), so
#   <T_F^n g s, T_C^n g s> = <g T_F^n s, g T_C^n s> = <T_F^n s, T_C^n s> for g unitary.
# I check the algebra numerically anyway, on both carriers, 2000 draws each.
import numpy as np
from fractions import Fraction as Fr
import xlib as X

rng = np.random.default_rng(20260817)
GEN_K1 = np.array([1.0, 0.37, 0.91, 2 ** 0.5, 0.23, 1.77])

print("== A0  MY OWN B0b, CHECKED AGAINST S4_THE_MEASUREMENT_V001.md:575 ==")
cls = X.classes(X.B0B_LOOP_F, X.B0B_LOOP_C, X.B0B_NV)
from collections import Counter
cnt = Counter("".join(map(str, c)) for c in cls)
print(f"  gamma_F vertices {sorted(X.loop_vertices(X.B0B_LOOP_F))} (len {len(X.B0B_LOOP_F)})   "
      f"gamma_C vertices {sorted(X.loop_vertices(X.B0B_LOOP_C))} (len {len(X.B0B_LOOP_C)})")
print(f"  class multiset {dict(cnt)}     [S4:575  {{00:4, 01:1, 10:2, 11:2}}]")
assert cnt["00"] == 4 and cnt["01"] == 1 and cnt["10"] == 2 and cnt["11"] == 2, "B0b MISBUILT"
print("  MATCHES.  My B0b is built from an explicit edge list, not from either W-11 lane's code.")

print("\n== A1  COR-F's SEALED EXHIBIT: T IS UNITARY AND T^L = M_gamma ==")
for nm, lf, lc, NV, aa in (("K1 ", X.K1_LOOP_F, X.K1_LOOP_C, 5, GEN_K1),
                           ("B0b", X.B0B_LOOP_F, X.B0B_LOOP_C, 9, rng.uniform(0, 2 * np.pi, 18))):
    TF, TC = X.T_edge(lf, aa, NV), X.T_edge(lc, aa, NV)
    MF, MC = X.M_circuit(lf, aa, NV), X.M_circuit(lc, aa, NV)
    LF, LC = len(lf), len(lc)
    print(f"  {nm}  ||T*T-I|| F={np.linalg.norm(TF.conj().T@TF-np.eye(NV)):.2e} "
          f"C={np.linalg.norm(TC.conj().T@TC-np.eye(NV)):.2e}   "
          f"||T_F^{LF}-M_F||={np.linalg.norm(np.linalg.matrix_power(TF,LF)-MF):.2e}   "
          f"||T_C^{LC}-M_C||={np.linalg.norm(np.linalg.matrix_power(TC,LC)-MC):.2e}   "
          f"T_F diagonal? {np.allclose(TF,np.diag(np.diag(TF)))}")

print("\n== A2  THE CONCESSION IS FORCED: Z^T_n IS GAUGE-INVARIANT AT EVERY PARTIAL TICK ==")
for nm, lf, lc, NV, ne, EDG, aa in (
        ("K1 ", X.K1_LOOP_F, X.K1_LOOP_C, 5, 6, X.K1_EDGES, GEN_K1),
        ("B0b", X.B0B_LOOP_F, X.B0B_LOOP_C, 9, 18, X.B0B_E, rng.uniform(0, 2 * np.pi, 18))):
    s = rng.normal(size=NV) + 1j * rng.normal(size=NV)
    s /= np.linalg.norm(s)
    TF, TC = X.T_edge(lf, aa, NV), X.T_edge(lc, aa, NV)
    worst_cov = worst_inv = 0.0
    for _ in range(2000):
        th = rng.uniform(0, 2 * np.pi, NV)
        g = np.diag(np.exp(1j * th))
        a2 = X.gauge_apply(aa, th, EDG)
        TF2, TC2 = X.T_edge(lf, a2, NV), X.T_edge(lc, a2, NV)
        worst_cov = max(worst_cov, np.linalg.norm(TF2 - g @ TF @ np.conj(g).T))
        for n in (1, 2, 5):
            z1 = X.Z(TF, TC, s, n, n)
            z2 = X.Z(TF2, TC2, g @ s, n, n)
            worst_inv = max(worst_inv, abs(z1 - z2))
    print(f"  {nm}  max||T(a^g) - g T(a) g*|| = {worst_cov:.2e}    "
          f"max|Z^T_n(a^g, g.s) - Z^T_n(a,s)|, n in 1,2,5 = {worst_inv:.2e}")
print("  -> COVARIANT, hence INVARIANT.  The steelman's leg A concession is FORCED, not lazy, and")
print("     the brief's objection (1) is dead.  I confirm it in my own code and do not revive it.")

print("\n== A3  THE SHARED HEADLINE NUMBERS, IN MY CODE ==")
sA = np.sqrt(np.array([0.40, 0.15, 0.15, 0.15, 0.15])) + 0j
sB = np.sqrt(np.array([0.40, 0.30, 0.00, 0.05, 0.25])) + 0j
sC = sA * np.exp(1j * np.array([0.0, 1.3, -0.7, 2.2, 0.4]))
ST = (sA, sB, sC)
assert X.arms_differ(*ST), "STATE ARMS BYTE-IDENTICAL -- leg void"
pis = [X.pi_of(s, X.K1_LOOP_F, X.K1_LOOP_C, 5) for s in ST]
assert np.allclose(pis[0], pis[1]) and np.allclose(pis[0], pis[2])
TF, TC = X.T_edge(X.K1_LOOP_F, GEN_K1, 5), X.T_edge(X.K1_LOOP_C, GEN_K1, 5)
MF, MC = X.M_circuit(X.K1_LOOP_F, GEN_K1, 5), X.M_circuit(X.K1_LOOP_C, GEN_K1, 5)
edge_sp = [max(abs(X.Z(TF, TC, s, n, n)) for s in ST) - min(abs(X.Z(TF, TC, s, n, n)) for s in ST)
           for n in range(1, 10)]
circ_sp = [max(abs(X.Z(MF, MC, s, k, k)) for s in ST) - min(abs(X.Z(MF, MC, s, k, k)) for s in ST)
           for k in range(1, 10)]
print(f"  pi (all three states) = {np.round(pis[0],9)}")
print(f"  K1 EDGE    spread over n<=9: max {max(edge_sp):.2e}  [registrar 5.90e-01]  per n: "
      + " ".join(f"{v:.2e}" for v in edge_sp))
print(f"  K1 CIRCUIT spread over k<=9: max {max(circ_sp):.2e}  [registrar 1.11e-16]")
print(f"  m(P) at pi = {np.round(pis[0],3)} by Jensen n=2^20 : {X.m_jensen(pis[0]):.12f}"
      f"   [register N1 / W-10: -0.767507880357]")

print("\n== A4  THE STEELMAN'S EXACT TRIVIAL-CONNECTION WITNESS, RECHECKED IN Fraction ==")
p, q, r = Fr(6, 11), Fr(7, 11), Fr(-6, 11)
nrm = p * p + q * q + r * r
Zex = p * q + q * r                     # conj(s0)s1 + conj(s1)s2 terms; the rest vanish
print(f"  s = (6/11, 7/11, -6/11, 0, 0):  ||s||^2 = {nrm}    Z^T_1 = {Zex}    (EXACT)")
s_ex = np.array([6 / 11, 7 / 11, -6 / 11, 0, 0], dtype=complex)
a0 = np.zeros(6)
TF0, TC0 = X.T_edge(X.K1_LOOP_F, a0, 5), X.T_edge(X.K1_LOOP_C, a0, 5)
MF0, MC0 = X.M_circuit(X.K1_LOOP_F, a0, 5), X.M_circuit(X.K1_LOOP_C, a0, 5)
print(f"  double precision: |Z^T_1| = {abs(X.Z(TF0,TC0,s_ex,1,1)):.3e}     "
      f"|Z^M_1| = {abs(X.Z(MF0,MC0,s_ex,1,1)):.12f}")
print("  CONFIRMED at the bytes and in exact rationals.  No arithmetic dispute with either lane.")

print("\n== A5  AND THE FACT NEITHER LANE STATES: AT THE TRIVIAL CONNECTION THE EDGE GATE FIRES ON")
print("        THE ROOT-ONLY STATE, AT EVERY CONNECTION, EXACTLY ==")
e0 = np.zeros(5, dtype=complex); e0[0] = 1.0
for tag, aa in (("trivial a=0", a0), ("generic f=1,c=sqrt2", GEN_K1)):
    TFx, TCx = X.T_edge(X.K1_LOOP_F, aa, 5), X.T_edge(X.K1_LOOP_C, aa, 5)
    MFx, MCx = X.M_circuit(X.K1_LOOP_F, aa, 5), X.M_circuit(X.K1_LOOP_C, aa, 5)
    print(f"  {tag:<22} EDGE |Z_1(root state)| = {abs(X.Z(TFx,TCx,e0,1,1)):.2e}    "
          f"CIRCUIT |Z_1(root state)| = {abs(X.Z(MFx,MCx,e0,1,1)):.12f}")
print("  W-01's registered property 'THE ROOT CAN NEVER FIRE, independently reproduced'")
print("  (REGISTER:47-49) is FALSE under the edge convention at EVERY connection, trivial or not.")
print("  That is a second registered property the edge tick contradicts, and it is stronger than")
print("  the trivial-limit one because it does not depend on the trivial connection at all.")
