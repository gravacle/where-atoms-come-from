# W-11 leg B — THE DECISIVE TEST.
# The corpus's whole carrier-independent layer rests on: Z_k depends on the ready state ONLY
# through pi = (p00,p10,p01,p11). That is what "the incidence labels are invisible" (W-03/N2)
# and "lambda = m(P)" (N1) MEAN. Under M_gamma it is true because both branch operators are
# DIAGONAL, so only |s_v|^2 can enter, and |s_v|^2 enters only through its class sum.
#
# ISOLATION LEDGER
#   Held fixed in every comparison: carrier K1, connection, the observable <branch_F, branch_C>,
#   the number of CIRCUITS, the code path, the seed. And pi is held fixed BY CONSTRUCTION.
#   Moved, one at a time:
#     (i)  the WITHIN-CLASS distribution of |s_v|^2, with every class sum identical
#     (ii) the PHASES of s, with every |s_v| identical
#     (iii) the TRANSPORT CONVENTION, M_gamma versus COR-F's T
#   pi cannot see (i) or (ii). The question is whether the FUNCTIONAL can.
import numpy as np
rng=np.random.default_rng(20260817)
LOOP_F=[(0,1,0),(1,2,1),(2,0,2)]; LOOP_C=[(0,3,3),(3,4,4),(4,0,5)]
FACE_V={0,1,2}; CYC_V={0,3,4}
def Top(loop,a):
    U=np.exp(1j*np.asarray(a)); T=np.zeros((5,5),dtype=complex); on={v for v,_,_ in loop}
    for v in range(5):
        if v not in on: T[v,v]=1.0
    for (s_,d_,e) in loop: T[d_,s_]=U[e]
    return T
def Mop(vs,W):
    M=np.eye(5,dtype=complex)
    for v in vs: M[v,v]=W
    return M
def hol(a): return np.exp(1j*(a[0]+a[1]+a[2])), np.exp(1j*(a[3]+a[4]+a[5]))
def pi_of(s):
    w=np.abs(s)**2
    return np.array([0.0, w[1]+w[2], w[3]+w[4], w[0]])      # p00,p10,p01,p11 on K1

a=np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])                # generic: c involves sqrt(2)
WF,WC=hol(a); TF,TC=Top(LOOP_F,a),Top(LOOP_C,a); MF,MC=Mop(FACE_V,WF),Mop(CYC_V,WC)

# three ready states: SAME pi, differing only within classes / in phase
sA=np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
sB=np.sqrt(np.array([0.40,0.30,0.00,0.05,0.25]))+0j          # same class sums, different split
sC=sA*np.exp(1j*np.array([0.0,1.3,-0.7,2.2,0.4]))            # same moduli, different phases
for nm,s in (("A",sA),("B",sB),("C",sC)):
    print(f"  pi({nm}) = {np.round(pi_of(s),12)}   |s|^2 = {np.round(np.abs(s)**2,4)}")
assert np.allclose(pi_of(sA),pi_of(sB)) and np.allclose(pi_of(sA),pi_of(sC))
print("  -> all three have IDENTICAL pi. Under the corpus's own theorems they are indistinguishable.\n")

def Zc(s,k,op_F,op_C):                                       # k CIRCUITS
    return np.vdot(np.linalg.matrix_power(op_F,k)@s, np.linalg.matrix_power(op_C,k)@s)
def Ze(s,n,):                                                # n EDGE TICKS
    return np.vdot(np.linalg.matrix_power(TF,n)@s, np.linalg.matrix_power(TC,n)@s)

print("== B1  UNDER THE CORPUS'S CONVENTION (M_gamma), CIRCUIT k ==")
print(f"  {'k':>3} {'|Z(A)|':>16} {'|Z(B)|':>16} {'|Z(C)|':>16} {'max spread':>13}")
for k in range(1,7):
    v=[abs(Zc(s,k,MF,MC)) for s in (sA,sB,sC)]
    print(f"  {k:>3} {v[0]:>16.12f} {v[1]:>16.12f} {v[2]:>16.12f} {max(v)-min(v):>13.2e}")
print("  -> identical. The invisibility theorem, working exactly as registered.\n")

print("== B2  UNDER COR-F's CONVENTION (T), EDGE TICK n.  ONE VARIABLE MOVED: THE CONVENTION ==")
print(f"  {'n':>3} {'|Z(A)|':>16} {'|Z(B)|':>16} {'|Z(C)|':>16} {'max spread':>13}  circuit?")
for n in range(1,10):
    v=[abs(Ze(s,n)) for s in (sA,sB,sC)]
    tag = f"= circuit {n//3}" if n%3==0 else ""
    print(f"  {n:>3} {v[0]:>16.12f} {v[1]:>16.12f} {v[2]:>16.12f} {max(v)-min(v):>13.2e}  {tag}")
print()
print("== B3  DO THE TWO CONVENTIONS AGREE WHERE THEY OVERLAP? ==")
w=max(abs(Ze(s,3*k)-Zc(s,k,MF,MC)) for s in (sA,sB,sC) for k in range(1,8))
print(f"  max | Z_edge(3k) - Z_circuit(k) | over 3 states, k<=7 = {w:.2e}")
print("  -> on K1 both loops have length 3, so the circuit convention is EXACTLY the edge")
print("     convention SAMPLED EVERY THIRD TICK. It is a subsequence, not a different object.")
