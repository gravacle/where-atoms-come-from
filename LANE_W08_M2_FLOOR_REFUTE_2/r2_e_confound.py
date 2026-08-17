# LANE W-08 / M2 REFUTER 2 — leg E: THE LANE'S DECLARED CONFOUND IS MIS-STATED ON THE CORPUS'S
# OWN READY STATE, AND THE ISOLATION IT DECLARED IMPOSSIBLE IS AVAILABLE THERE.
# Lane, leg D header and ISOLATION LEDGER ROW 5, stated with no ready-state restriction:
#   "|Z_k| = 1 for some k >= 1  <=>  alpha and beta are both rational  <=>  H = closure{(x^k,y^k)}
#    is FINITE.  So attainment cannot be moved without moving H."
# That biconditional holds on RS-G (all three weights positive).  On RS-P (K1's PUBLISHED ready
# state, p10 = 0) it is FALSE in the second half: attainment depends on alpha ALONE.
# ISOLATION: carrier, ready state RS-P, observable, k-grid, code path held.  MOVES: beta only.
import numpy as np, mpmath as mp
mp.mp.dps=40
FACE_V=[0,1,2]; CYC_V=[0,3,4]
def Z_ops(f,c,p,K):
    WF=np.exp(1j*f); WC=np.exp(1j*c); p=np.asarray(p,float)
    inF=np.array([1.0 if v in FACE_V else 0.0 for v in range(5)])
    inC=np.array([1.0 if v in CYC_V else 0.0 for v in range(5)])
    return np.array([np.sum(p*np.conj(WF**(k*inF))*(WC**(k*inC))) for k in range(1,K+1)])
RSP=[0.5,0,0,0.25,0.25]
astar=float(np.mod(2*np.cos(2*np.pi/7),1))
print("== E1  A CONNECTION THAT IS ATTAINED AND HAS INFINITE H, ON K1'S PUBLISHED READY STATE ==")
print("   alpha = 1/2 (f = pi, S1's own curvature), beta = the cubic irrational (c = 2pi beta*).")
print("   x = conj(W_F) = -1 has order 2; y = W_C has infinite order; so H = Z/2 x T is INFINITE,")
print("   while |Z_k| = 1 EXACTLY on every even k.  Attained AND H infinite, simultaneously.")
f=np.pi; c=2*np.pi*astar
Z=Z_ops(f,c,RSP,20)
print("   k      : " + " ".join(f"{k:>7d}" for k in range(1,11)))
print("   |Z_k|  : " + " ".join(f"{abs(z):>7.4f}" for z in Z[:10]))
print(f"   cells with |Z_k| = 1 to 1e-15, k <= 20: {int(np.sum(np.abs(np.abs(Z)-1)<1e-15))}"
      f"   cells with Z_k = 0 exactly: {int(np.sum(np.abs(Z)<1e-15))}")
print("   CONSEQUENCE FOR THE LANE'S LEDGER ROW 5: on RS-P attainment CAN be moved while H stays")
print("   infinite (move alpha 1/2 -> irrational with beta held irrational), so the confound the")
print("   lane declared unavoidable is an artefact of RS-G.  On RS-P the operative variable is")
print("   the closure of {x^k} ALONE -- W_C is invisible, which is the lane's own M2-9.\n")

print("== E2  AND THE SAME ROW IS THE CORPUS'S PUBLISHED POINT, WHERE THE 'COSMETIC' CLAIM DIES ==")
print("   S1's published connection AND S1's published ready state (W-01's registered firing):")
Z0=Z_ops(np.pi,3*np.pi/2,RSP,8)
print("   |Z_k| for k=1..8 : " + " ".join(f"{abs(z):.3f}" for z in Z0))
print("   Z_1 = 0 EXACTLY -> Omega_N = 0 for every N >= 1: the record is COMPLETE AFTER ONE")
print("   CIRCUIT.  An APPROACHED connection on the same ready state gives |Omega_N| = 2^{-N}(1+o(1)).")
print("   The honest-schedule distinction at the corpus's own point is therefore: complete at")
print("   N = 1 versus N ~ 46 circuits to reach the same 1e-14.  Not 4.8%, and not on the")
print("   schedule axis -- on the honest one.")
N=np.arange(1,60)
print(f"   approached: N needed for |Omega_N| < 1e-14 at rate -log2 : "
      f"{int(np.ceil(14*np.log(10)/np.log(2)))} circuits;  attained: 1 circuit.")
