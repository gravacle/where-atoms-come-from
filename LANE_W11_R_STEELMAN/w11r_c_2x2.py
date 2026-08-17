# LANE W11-R  LEG C -- NAMING THE OPERATIVE VARIABLE.  A 2x2, ONE VARIABLE PER EDGE OF THE DESIGN.
#
# The registrar's legs B and D name the moved variable "THE TRANSPORT CONVENTION".  Its own leg B3
# reports  max|Z_edge(3k) - Z_circuit(k)| = 1.97e-15  on K1 and calls the circuit convention "the
# edge convention SAMPLED EVERY THIRD TICK ... a subsequence, not a different object."  Those two
# sentences cannot both be right about the same run.  This leg decides which.
#
# THE DESIGN
#   TRANSPORT axis:  T = COR-F's edge tick (non-diagonal)   vs   D = the uniform (principal) root
#                    diag(W^{1/L} on the loop, 1 off), which also satisfies D^L = M_gamma exactly.
#   CLOCK axis:      EDGE  -- compare branch_F at n ticks with branch_C at n ticks
#                    CIRCUIT -- compare branch_F at L_F*n ticks with branch_C at L_C*n ticks
#                               (each branch at ITS OWN loop closure)
#   Held fixed everywhere: carrier, connection, the three pi-identical ready states, the observable,
#   the code path, the seed.  pi identical BY CONSTRUCTION.
import numpy as np, w11r_lib as L
rng = np.random.default_rng(20260817)

def three_states_K1():
    sA = np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
    sB = np.sqrt(np.array([0.40,0.30,0.00,0.05,0.25]))+0j
    sC = sA*np.exp(1j*np.array([0.0,1.3,-0.7,2.2,0.4]))
    return sA,sB,sC
def three_states_B0b():
    wA = np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wA/=wA.sum()
    wB = wA.copy()
    wB[0],wB[1] = wA[0]+wA[1], 0.0      # class 11 = {0,1}
    wB[3],wB[4] = 0.0, wA[3]+wA[4]      # class 10 = {3,4}
    wB[5],wB[8] = wA[5]+wA[8], 0.0      # class 00 = {5,6,7,8}
    sA,sB = np.sqrt(wA)+0j, np.sqrt(wB)+0j
    sC = sA*np.exp(1j*rng.uniform(0,2*np.pi,9))
    return sA,sB,sC

CASES = [("K1 ", L.K1_LOOP_F, L.K1_LOOP_C, 5, np.array([1.0,0.37,0.91,2**0.5,0.23,1.77]), three_states_K1),
         ("B0b", L.B0B_LOOP_F, L.B0B_LOOP_C, 9, rng.uniform(0,2*np.pi,18),               three_states_B0b)]

for nm, lf, lc, NV, aa, mk in CASES:
    LF, LC = len(lf), len(lc)
    sA,sB,sC = mk(); states = (sA,sB,sC)
    pis = [L.pi_of(s,lf,lc,NV) for s in states]
    assert np.allclose(pis[0],pis[1]) and np.allclose(pis[0],pis[2]), "pi not held fixed"
    assert L.arms_differ(sA,sB,sC), "STATE ARMS BYTE-IDENTICAL -- leg void"
    TF,TC = L.T_edge(lf,aa,NV), L.T_edge(lc,aa,NV)
    DF,DC = L.D_uniform(lf,aa,NV), L.D_uniform(lc,aa,NV)
    MF,MC = L.M_circuit(lf,aa,NV), L.M_circuit(lc,aa,NV)
    assert L.arms_differ(TF,DF) and L.arms_differ(TC,DC), "TRANSPORT ARMS BYTE-IDENTICAL -- leg void"
    print(f"\n================ {nm}   |gamma_F| = {LF}, |gamma_C| = {LC},   pi = {np.round(pis[0],9)}")
    print(f"  ARM DIFF: T_F vs D_F bytes differ = {L.arms_differ(TF,DF)};  ||T_F - D_F|| = {np.linalg.norm(TF-DF):.4f}")
    print(f"            ||T_F^{LF} - M_dF|| = {np.linalg.norm(np.linalg.matrix_power(TF,LF)-MF):.2e}"
          f"   ||D_F^{LF} - M_dF|| = {np.linalg.norm(np.linalg.matrix_power(DF,LF)-MF):.2e}")
    print(f"            ||[T_F,T_C]|| = {np.linalg.norm(TF@TC-TC@TF):.4f}    ||[D_F,D_C]|| = {np.linalg.norm(DF@DC-DC@DF):.2e}")
    print(f"  {'cell':<26}{'n':>3} {'|Z(A)|':>16}{'|Z(B)|':>16}{'|Z(C)|':>16}{'spread':>11}")
    for tname,(oF,oC) in (("T = COR-F edge tick",(TF,TC)), ("D = uniform root  ",(DF,DC))):
        for cname, (kF,kC) in (("EDGE clock   ",(1,1)), ("CIRCUIT clock",(LF,LC))):
            worst = 0.0
            for n in range(1,7):
                v = [abs(L.Z(oF,oC,s,kF*n,kC*n)) for s in states]
                worst = max(worst, max(v)-min(v))
                if n <= 3:
                    print(f"  {tname} / {cname}{n:>3} {v[0]:>16.12f}{v[1]:>16.12f}{v[2]:>16.12f}{max(v)-min(v):>11.1e}")
            print(f"  {'':<26}{'':>3} {'':>16}{'':>16}{'  worst spread n<=6':>16}{worst:>11.1e}")
    # does the corpus's own functional coincide with COR-F's transport at loop closure?
    w = max(abs(L.Z(TF,TC,s,LF*k,LC*k) - L.Z(MF,MC,s,k,k)) for s in states for k in range(1,9))
    print(f"  ==> max | Z[COR-F's T, at each branch's own loop closure] - Z[corpus's M_gamma] | = {w:.2e}")

print("\n== C-RULING ==")
print("  The 2x2 separates cleanly:")
print("     (T, EDGE clock)    -> incidence VISIBLE")
print("     (T, CIRCUIT clock) -> incidence INVISIBLE, and the functional is the corpus's, exactly")
print("     (D, EDGE clock)    -> incidence INVISIBLE")
print("     (D, CIRCUIT clock) -> incidence INVISIBLE")
print("  Visibility needs BOTH a non-fibre-wise root AND an edge clock.  Neither alone is 'the")
print("  convention'.  The registrar moved TWO things at once on B0b and, on K1, moved the")
print("  TRANSPORT NOT AT ALL -- its own leg B3 proves M_gamma = T^3 there.")
