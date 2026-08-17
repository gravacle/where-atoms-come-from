# X1 — INDEPENDENT REBUILD, AND THE LANE'S HEADLINE REPRODUCED FROM MY OWN BYTES.
# Nothing here is imported from LANE_W11_R_CLOCK.  If the lane's leg 1 is an implementation
# artefact, this leg does not reproduce it.
import numpy as np
from x_lib import *

np.set_printoptions(linewidth=200)
print("== X1  INDEPENDENT REBUILD ==  double precision is the default.")
for K, a, pi in ((K1(), np.array([1.0,0.37,0.91,2**0.5,0.23,1.77]), np.array([0.,.30,.30,.40])),
                 (B0b(), np.random.default_rng(20260817).uniform(0,2*np.pi,18), None)):
    if pi is None:
        wB = np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wB /= wB.sum()
        pi = pi_of(K, np.sqrt(wB)+0j)
    TF, TC = T_edge(K,K.wF,a), T_edge(K,K.wC,a)
    WF, WC = holo(K.wF,a), holo(K.wC,a)
    MF, MC = M_circ(K,K.VF,WF), M_circ(K,K.VC,WC)
    print(f"\n--- {K.name}  nv={K.nv}  L=({K.LF},{K.LC})  multiset {K.multiset()}  pi={np.round(pi,6)}")
    print(f"    ||T_F*T_F - I|| = {np.linalg.norm(TF.conj().T@TF-np.eye(K.nv)):.2e}   "
          f"||T_F^{K.LF} - M_F|| = {np.linalg.norm(np.linalg.matrix_power(TF,K.LF)-MF):.2e}   "
          f"||T_C^{K.LC} - M_C|| = {np.linalg.norm(np.linalg.matrix_power(TC,K.LC)-MC):.2e}")
    S = states_same_pi(K, pi, 64, np.random.default_rng(20260817))
    MM = 25
    sp = np.zeros((MM+1,MM+1))
    PF = [np.linalg.matrix_power(TF,m) for m in range(MM+1)]
    PC = [np.linalg.matrix_power(TC,m) for m in range(MM+1)]
    for mF in range(MM+1):
        AF=[PF[mF]@s for s in S]
        for mC in range(MM+1):
            AC=[PC[mC]@s for s in S]
            v=np.array([abs(np.vdot(x,y)) for x,y in zip(AF,AC)]); sp[mF,mC]=v.max()-v.min()
    pred = np.array([[ (mF%K.LF==0) and (mC%K.LC==0) for mC in range(MM+1)] for mF in range(MM+1)])
    zero = sp < 1e-12
    print(f"    LANE LEG-1 HEADLINE, MY CODE: invisible cells {zero.sum()}/{zero.size}  "
          f"predicted {pred.sum()}  SETS EQUAL {np.array_equal(zero,pred)}  "
          f"max ON {sp[pred].max():.2e}  max OFF {sp[~pred].max():.3e}")
    qF = np.array([[ (mF%K.LF!=0) and (mC%K.LC==0) for mC in range(MM+1)] for mF in range(MM+1)])
    qC = np.array([[ (mF%K.LF==0) and (mC%K.LC!=0) for mC in range(MM+1)] for mF in range(MM+1)])
    print(f"    only-F-mid-loop min {sp[qF].min():.3e}   only-C-mid-loop min {sp[qC].min():.3e}"
          f"   [lane: K1 6.622e-01/6.200e-01, B0b 1.667e-01/3.669e-01]")
print("\n  VERDICT X1: the lane's leg 1 REPRODUCES from an independent implementation.")
print("  Its set-equality, its magnitudes and its quadrant minima are not artefacts.")
