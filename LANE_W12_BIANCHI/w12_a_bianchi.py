# W-12 leg A — does the Bianchi constraint bite on a carrier with b2 >= 1?
# The proposal (registrar, from the DNA linking-number analogy): N3 kills every ABSOLUTELY
# CONTINUOUS connection measure, so no Wilson action moves the rate. A topological constraint is
# SINGULAR, so it escapes N3's wording. K1 has b2 = 0 and S1:102-104 says Bianchi is vacuous there.
# S4 owns five carriers with b2 >= 1 and never imposed it.
#
# ISOLATION: same connection sampler, same rate estimator, same designation. The ONE thing that
# moves is whether the carrier has closed 2-cycles.
import numpy as np, itertools
from math import gcd
np.set_printoptions(precision=6, suppress=True)

def smith2(M):
    """invariant factors of an integer matrix of rank <= 2"""
    ent = [abs(int(x)) for x in M.ravel() if x]
    s1 = 0
    for x in ent: s1 = gcd(s1, x)
    minors = []
    for i, j in itertools.combinations(range(M.shape[0]), 2):
        for k, l in itertools.combinations(range(M.shape[1]), 2):
            minors.append(abs(int(round(M[i,k]*M[j,l] - M[i,l]*M[j,k]))))
    g2 = 0
    for x in minors: g2 = gcd(g2, x)
    return s1, (g2 // s1 if s1 else 0)

# ---- B0b, incidence verbatim from LANE_W10_A_CARRIERS/w10a_lib.py:182-213
def V(i,j): return 3*(j%3)+(i%3)
def H(i,j): return 3*(j%3)+(i%3)
def Wg(i,j): return 9+3*(j%3)+(i%3)
E_b0b=[None]*18
for j in range(3):
    for i in range(3):
        E_b0b[H(i,j)]=(V(i,j),V(i+1,j)); E_b0b[Wg(i,j)]=(V(i,j),V(i,j+1))
FACES_b0b=[]
for j in range(3):
    for i in range(3):
        FACES_b0b.append([(H(i,j),1),(Wg(i+1,j),1),(H(i,j+1),-1),(Wg(i,j),-1)])
gF_b0b=[(H(0,0),1),(Wg(1,0),1),(H(0,1),-1),(Wg(0,0),-1)]
gC_b0b=[(H(0,0),1),(H(1,0),1),(H(2,0),1)]

def d2_matrix(E,FACES):
    M=np.zeros((len(E),len(FACES)))
    for fi,f in enumerate(FACES):
        for (e,s) in f: M[e,fi]+=s
    return M

print("== A1  B0b: b2, the 2-cycles, and what Bianchi actually says ==")
d2=d2_matrix(E_b0b,FACES_b0b)
# 2-cycles = kernel of d2 (a 2-chain with zero boundary)
u_,s_,vt_=np.linalg.svd(d2)
ker=vt_[np.sum(s_>1e-9):]
print(f"  E={len(E_b0b)} F={len(FACES_b0b)}  rank(d2)={int(np.sum(s_>1e-9))}  dim ker(d2) = b2 = {ker.shape[0]}")
z=ker[0]/np.abs(ker[0][np.argmax(np.abs(ker[0]))])
print(f"  the 2-cycle (up to scale): {np.round(z,6)}   -> it is the sum of ALL NINE faces")
print()
print("  BIANCHI for a 2-cycle z:  sum_F z_F * f(F) = a(boundary of z) = a(0) = 0.")
print("  With the connection given by EDGE PHASES a_e -- which is how S1 sec3 defines it --")
print("  this holds IDENTICALLY for every a. Verified over 2000 random connections:")
worst=0.0
rng=np.random.default_rng(20260820)
for _ in range(2000):
    a=rng.uniform(0,2*np.pi,18)
    f=np.array([sum(s*a[e] for (e,s) in F) for F in FACES_b0b])
    worst=max(worst,abs(float(z@f)))
print(f"    max | sum_F z_F f(F) |  = {worst:.2e}    (exactly 0 up to float)")
print()
print("  ==> BIANCHI IS NOT A CONSTRAINT ON THE CONNECTION. It is an IDENTITY.")
print("      It restricts which CURVATURE ASSIGNMENTS are realizable, not which connections exist.")
print("      MY PROPOSAL WAS WRONG AS STATED, and this is the reason.\n")

print("== A2  AND IT NEVER PINS A SINGLE DESIGNATED HOLONOMY ==")
print("  W_F is ONE face curvature among nine; one linear relation among nine quantities")
print("  leaves any single one of them free. Verified: sweep the connection and see whether")
print("  (W_F, W_C) covers the torus.")
def hol(g,a):
    return float(sum(s*a[e] for (e,s) in g))
pts=[]
for _ in range(20000):
    a=rng.uniform(0,2*np.pi,18)
    pts.append((hol(gF_b0b,a)%(2*np.pi), hol(gC_b0b,a)%(2*np.pi)))
pts=np.array(pts)
Hgrid,_,_=np.histogram2d(pts[:,0],pts[:,1],bins=20,range=[[0,2*np.pi]]*2)
print(f"  20000 draws into a 20x20 torus grid: cells occupied {int((Hgrid>0).sum())} of 400,"
      f"  min {int(Hgrid.min())}, max {int(Hgrid.max())}")
print("  -> (W_F, W_C) covers T^2. Bianchi constrains it nowhere.")
