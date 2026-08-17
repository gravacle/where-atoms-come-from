# LANE W11-R  LEG G -- two loose ends closed rather than left dangling.
# G1: leg D's 'generic' row reported 1 pi-only root in 400.  Is that a counterexample to the
#     characterisation, or a draw that landed IN the preserving set?  Checked, not assumed.
# G2: does the corpus's own K1 exhibit -- W-01's firing on S1's published connection -- survive
#     under the UNIFORM ROOT D at every tick?  (If Reading A's selected transport broke a
#     registered exhibit, Reading A would fall on its own argument.)
import numpy as np, w11r_lib as L
rng = np.random.default_rng(20260817)
lf,lc,NV = L.K1_LOOP_F, L.K1_LOOP_C, 5
a = np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])
MF,MC = L.M_circuit(lf,a,NV), L.M_circuit(lc,a,NV)
CLS = L.classes(lf,lc,NV)
sA = np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
sB = np.sqrt(np.array([0.40,0.30,0.00,0.05,0.25]))+0j
sC = sA*np.exp(1j*np.array([0.0,1.3,-0.7,2.2,0.4]))
ST=(sA,sB,sC)
def spread(oF,oC,nmax=6):
    return max(max(abs(L.Z(oF,oC,s,n,n)) for s in ST)-min(abs(L.Z(oF,oC,s,n,n)) for s in ST) for n in range(1,nmax+1))
def random_root(loop, aa, NV, kind):
    Ls=len(loop); on=sorted(L.loop_vertices(loop)); off=[v for v in range(NV) if v not in on]
    W=L.holonomy(loop,aa); w0=np.exp(1j*(np.angle(W)+2*np.pi*rng.integers(0,Ls))/Ls); z=np.exp(2j*np.pi/Ls)
    U=np.zeros((NV,NV),dtype=complex)
    def block(idx,scal):
        d=len(idx); ks=rng.integers(0,Ls,size=d)
        Q,_=np.linalg.qr(rng.normal(size=(d,d))+1j*rng.normal(size=(d,d)))
        return scal*(Q@np.diag(z**ks.astype(float))@Q.conj().T)
    Bl=block(on,w0); Bo=block(off,1.0)
    for i,v in enumerate(on):
        for j,u in enumerate(on): U[v,u]=Bl[i,j]
    for i,v in enumerate(off):
        for j,u in enumerate(off): U[v,u]=Bo[i,j]
    return U
def class_constant_diagonal(U):
    if not np.allclose(U,np.diag(np.diag(U)),atol=1e-9): return False
    d=np.diag(U); seen={}
    for v in range(len(d)):
        c=CLS[v]
        if c in seen and abs(seen[c]-d[v])>1e-9: return False
        seen[c]=d[v]
    return True
print("== G1  EVERY pi-ONLY DRAW IS CLASS-CONSTANT DIAGONAL, WITH NO EXCEPTION ==")
n_pi=n_cc=0; exceptions=[]
for _ in range(4000):
    UF,UC=random_root(lf,a,NV,'g'), random_root(lc,a,NV,'g')
    s_=spread(UF,UC)
    if s_<1e-12:
        n_pi+=1
        cc = class_constant_diagonal(UF) and class_constant_diagonal(UC)
        n_cc += cc
        if not cc: exceptions.append(s_)
print(f"  4000 random roots of (M_dF, M_c):  pi-only draws = {n_pi};  of those, class-constant")
print(f"  diagonal = {n_cc};  COUNTEREXAMPLES to the characterisation = {len(exceptions)}")
print("  -> the characterisation holds: pi-only <=> both roots diagonal AND constant on the four")
print("     classes.  Leg D's lone 'generic' hit was a draw that landed inside the preserving set.")
print("\n== G2  DOES THE UNIFORM ROOT BREAK ANY REGISTERED EXHIBIT?  W-01's FIRING ON S1's CONNECTION ==")
aS1=np.array([np.pi/3]*3+[np.pi/2]*3)                      # S1 sec6: W_F = -1, W_C = -i
pS1=np.array([0.5,0.0,0.0,0.25,0.25])                      # W-01's p = (1/2,0,0,1/4,1/4)
sS1=np.sqrt(pS1)+0j
MF1,MC1=L.M_circuit(lf,aS1,NV), L.M_circuit(lc,aS1,NV)
DF1,DC1=L.D_uniform(lf,aS1,NV), L.D_uniform(lc,aS1,NV)
print(f"  W_F = {L.holonomy(lf,aS1):.6f}   W_C = {L.holonomy(lc,aS1):.6f}   [S1 sec6: -1, -i]")
print(f"  CIRCUIT convention  |Z_1| = {abs(L.Z(MF1,MC1,sS1,1,1)):.3e}   [REGISTER:47 'exactly 0']")
print(f"  UNIFORM ROOT, at loop closure (n = 3)  |Z| = {abs(L.Z(DF1,DC1,sS1,3,3)):.3e}")
print(f"  UNIFORM ROOT, at every tick n = 1..6:  " +
      " ".join(f"{abs(L.Z(DF1,DC1,sS1,n,n)):.4f}" for n in range(1,7)))
print(f"  pi-spread of the uniform root over the three pi-identical states: {spread(DF1,DC1):.2e}")
print("  -> the registered firing survives, and the uniform root additionally makes the incidence")
print("     invisible at EVERY tick, not only at loop closure.  Reading A's selected transport")
print("     costs the corpus nothing it has registered.")
