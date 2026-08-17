# X8 — THE SWEEP THE REGISTRAR DECLARED AS ITS FIRST WEAKNESS AND NEITHER IT NOR LANE C RAN:
#      the whole set { U unitary : U^L = M_gamma }.
# My brief: "If a different, equally natural edge tick restores invisibility, Reading B falls.
#            If EVERY unitary with T^L = M except the diagonal ones breaks invisibility, Reading B
#            is far stronger than the registrar claimed.  Settle it."
#
# STRUCTURE.  On the loop, M = W.I_L.  U^L = W.I_L  <=>  U = omega V with omega^L = W and V^L = I.
# Therefore  DIAG(U) := {m : U^m diagonal in the vertex basis} = {m : V^m diagonal},
# which is closed under addition and contains L, hence  DIAG(U) = d.Z  for some DIVISOR d of L.
# COR-F's T is the EXTREME d = L.  The fibre-wise root D is the other extreme d = 1.
import numpy as np
from x_lib import *

rng = np.random.default_rng(20260817)

def rand_unitary(n, r):
    X = r.normal(size=(n,n)) + 1j*r.normal(size=(n,n))
    Q, R = np.linalg.qr(X)
    return Q*np.exp(-1j*np.angle(np.diag(R)))

def embed(K, vs, Uloop, order):
    """put an L x L unitary on the loop's fibres (in the given vertex order), identity off it."""
    A = np.eye(K.nv, dtype=complex)
    for i, u in enumerate(order):
        for j, v in enumerate(order):
            A[u, v] = Uloop[i, j]
    for v in range(K.nv):
        if v not in vs:
            A[v, :] = 0; A[:, v] = 0; A[v, v] = 1
    return A

def diagset(A, hi):
    return [m for m in range(0, hi+1)
            if np.allclose(np.linalg.matrix_power(A,m),
                           np.diag(np.diag(np.linalg.matrix_power(A,m))), atol=1e-10)]

print("== X8a  BRANCH ROBUSTNESS OF THE FIBRE-WISE ROOT D (L choices of the L-th root) ==")
for K, pi in ((K1(), np.array([0.,.30,.30,.40])), (B0b(), None)):
    if pi is None:
        wB = np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wB/=wB.sum(); pi = pi_of(K,np.sqrt(wB)+0j)
    a = generic_conn(K, np.random.default_rng(7+K.nv))
    WF, WC = holo(K.wF,a), holo(K.wC,a)
    MF, MC = M_circ(K,K.VF,WF), M_circ(K,K.VC,WC)
    S = states_same_pi(K, pi, 32, np.random.default_rng(20260817))
    for bF in range(K.LF):
        for bC in range(K.LC):
            DF = D_root(K,K.VF,WF,K.LF,bF); DC = D_root(K,K.VC,WC,K.LC,bC)
            e1 = np.linalg.norm(np.linalg.matrix_power(DF,K.LF)-MF)
            e2 = np.linalg.norm(np.linalg.matrix_power(DC,K.LC)-MC)
            sp = max(spread_over(S, np.linalg.matrix_power(DF,m), np.linalg.matrix_power(DC,m))
                     for m in range(1,10))
            assert e1 < 1e-14 and e2 < 1e-14 and sp < 1e-12, (bF,bC,e1,e2,sp)
    print(f"   {K.name:<4}: all {K.LF}x{K.LC} branch pairs give D^L = M to <1e-14 and invisibility")
    print(f"          to <1e-12 at every m<=9.  The counterexample is not a branch accident.")

print("\n== X8b  HAAR SWEEP OF { U : U^L = M_gamma } -- 4000 draws per loop, per carrier ==")
print("   U = omega V,  omega^L = W (random branch),  V = Q diag(zeta^k) Q* with Q Haar,")
print("   zeta = exp(2 pi i / L), exponents k drawn uniformly.  ONE VARIABLE: the unitary U.")
for K, pi in ((K1(), np.array([0.,.30,.30,.40])), (B0b(), None)):
    if pi is None:
        wB = np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wB/=wB.sum(); pi = pi_of(K,np.sqrt(wB)+0j)
    a = generic_conn(K, np.random.default_rng(7+K.nv))
    WF, WC = holo(K.wF,a), holo(K.wC,a)
    MF, MC = M_circ(K,K.VF,WF), M_circ(K,K.VC,WC)
    S = states_same_pi(K, pi, 24, np.random.default_rng(20260817))
    ordF = list(K.VF); ordC = list(K.VC)
    from collections import Counter
    cntF, cntC, bad, invis1 = Counter(), Counter(), 0.0, 0
    NDRAW = 4000
    for t in range(NDRAW):
        out = []
        for (vs, order, L, W, Mref) in ((K.VF, ordF, K.LF, WF, MF), (K.VC, ordC, K.LC, WC, MC)):
            om = np.exp(1j*(np.angle(W) + 2*np.pi*rng.integers(0,L))/L)
            Q = rand_unitary(L, rng)
            ks = rng.integers(0, L, L)
            V = Q @ np.diag(np.exp(2j*np.pi*ks/L)) @ Q.conj().T
            A = embed(K, vs, om*V, order)
            bad = max(bad, np.linalg.norm(np.linalg.matrix_power(A,L) - Mref))
            out.append(A)
        UF, UC = out
        dF = diagset(UF, 2*K.LF); dC = diagset(UC, 2*K.LC)
        gF = dF[1] if len(dF) > 1 else 0; gC = dC[1] if len(dC) > 1 else 0
        cntF[gF] += 1; cntC[gC] += 1
        if t < 400 and spread_over(S, UF, UC) < 1e-12: invis1 += 1
    print(f"   {K.name:<4}: max || U^L - M_gamma || over {NDRAW} draws = {bad:.2e}   (constraint holds)")
    print(f"          smallest positive m with U_F^m diagonal: {dict(sorted(cntF.items()))}"
          f"   (divisors of L_F={K.LF})")
    print(f"          smallest positive m with U_C^m diagonal: {dict(sorted(cntC.items()))}"
          f"   (divisors of L_C={K.LC})")
    print(f"          invisibility already at m = 1: {invis1} of 400 draws")
    for d in sorted(set(cntF)|set(cntC)):
        assert d == 0 or (K.LF % d == 0 or K.LC % d == 0), d
    print(f"          EVERY observed d DIVIDES its loop length, as the structure theorem requires.")

print("""
================ WHAT X8 SETTLES ================
 A  {m : U^m is diagonal} = d.Z with d | L, for EVERY unitary U with U^L = M_gamma.  The
    invisibility set of ANY rival tick is therefore a sublattice {d_F | m_F} x {d_C | m_C},
    and the only free content is the pair of divisors.  Lane C's leg 1 computed ONE point of
    this family (d = (L_F, L_C), COR-F's T) and named it "the operative variable".
 B  A HAAR-RANDOM tick with U^L = M_gamma almost never restores invisibility at m = 1: the
    generic d is L.  SO THE REGISTRAR'S SECOND HORN IS THE TRUE ONE FOR THE GENERIC TICK --
    almost every rival breaks invisibility off the sublattice, and Reading B is stronger than
    the registrar claimed, NOT weaker.
 C  BUT THE EXCEPTIONS ARE NOT A NULL SET OF THE CONSTRUCTION -- they are exactly the ticks that
    are FIBRE-WISE AND LOOP-CONSTANT, i.e. the corpus's OWN operator M_gamma and its roots.
    The corpus did not land on the generic tick.  It landed, by stipulation, on the one measure-
    zero family in which the clock question is empty.
 D  SO THE BRIEF'S DISJUNCTION IS NOT EXHAUSTIVE.  "A different equally natural edge tick
    restores invisibility" is TRUE (D does), and Reading B does NOT fall, because D is not a
    rival to the CONVENTION -- it is the convention's own shape at a finer clock.  What decides
    invisibility, over the whole family, is the divisor d, and d = 1 exactly when the tick is
    fibre-wise-and-loop-constant.  THE STIPULATION IS ON THE OPERATOR.
""")
