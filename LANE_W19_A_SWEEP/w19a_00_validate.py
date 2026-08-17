# W19-A step 0. VALIDATE THE INSTRUMENT BEFORE USING IT.
# (1) reproduce the sealed T1 numbers exactly, with a DIFFERENT method
#     (T1: 8x8 full space + 1e6 penalty projection.  here: exact 4-dim physical sector, gauge-fixed)
# (2) check the Gauss law holds on the reconstructed full-space vector
# (3) check the claimed exact structure: rho_link is X-diagonal, and X_l = X(C_i)^{-1} on physical states
import numpy as np, sys
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_A_SWEEP")
from zn_gauge import *
from carriers import *

np.set_printoptions(precision=6, suppress=True)
print("="*100)
print("W19-A / 00 VALIDATE.  Instrument = exact physical-sector Z_N gauge theory (gauge-fixed orbit basis).")
print("="*100)

# ---------------------------------------------------------------- (1) T1 reproduction
print("\n[1] T1 REPRODUCTION.  Sealed T1 (LANE_T1_NEW_PROGRAM/PUBLISHED_CONVENTIONS.txt) reported")
print("    I(1:2)=0.690763  I(1:3)=0.384496  I(1:{2,3})=1.075259  for")
print("    H = -0.7*(X1+X2+X3) - (Z1Z2 + Z2Z3) on the Z_2 theta graph, physical sector.")
V,E = theta()
g = ZNGauge("theta", V, E, 2)
print(f"    theta: V={g.V} L={g.L} C={g.C} dim_phys={g.dimP}   (T1 said physical dim 4)")
# T1's plaquettes were specifically W12 = Z1Z2 and W23 = Z2Z3.  Override the MCB choice to match.
g.plaq  = [np.array([1,-1,0]), np.array([0,1,-1])]
g.avec  = [tuple(int(p[g.chords[c]])%2 for c in range(g.C)) for p in g.plaq]
psi,E0,gap = g.ground(mag=1.0, elec=0.7)          # H = -1*sum_p W_p - 0.7*sum_e X_e  for Z_2
Psi = g.full_vector(psi)
for F,lab in [([1],"I(1:2)"),([2],"I(1:3)"),([1,2],"I(1:{2,3})")]:
    print(f"    {lab:<12} = {mutual_information(Psi,g.L,2,[0],F):.6f}")
print(f"    S(link 1) = {S_of(Psi,g.L,2,[0]):.6f} bits;  E0 = {E0:.6f};  gap = {gap:.6f}")

# ---------------------------------------------------------------- (2) Gauss law on the full vector
print("\n[2] GAUSS LAW CHECK on the reconstructed full-space vector (no penalty used anywhere).")
def gauss_residual(gg, Psi):
    N,L,V_ = gg.N, gg.L, gg.V
    T = Psi.reshape([N]*L); worst = 0.0
    for v in range(V_):
        sh = [ (1 if b==v else 0) - (1 if a==v else 0) for (a,b) in gg.edges ]
        Tv = T
        for e,s in enumerate(sh):
            if s % N: Tv = np.roll(Tv, s % N, axis=e)     # X_e^{s}: |j> -> |j+s>
        worst = max(worst, float(np.abs(Tv-T).max()))
    return worst
for nm,(V_,E_),N in [("theta",theta(),2),("ladder_2sq",ladder(2),2),("cube_Q3",cube(),2),
                     ("ladder_2sq",ladder(2),3)]:
    gg = ZNGauge(nm,V_,E_,N); p,_,_ = gg.ground(2.0,2.0); P = gg.full_vector(p)
    print(f"    {nm:<12} Z_{N}: L={gg.L} C={gg.C} dimP={gg.dimP:<6} ||G_v Psi - Psi||_inf = {gauss_residual(gg,P):.3e}"
          f"   norm={np.linalg.norm(P):.12f}")

# ---------------------------------------------------------------- (3) exact structure claims
print("\n[3] EXACT STRUCTURE CLAIMS (these are theorems; the numbers only confirm the code).")
print("    (3a) rho_link is diagonal in the X (electric) basis, so H(S) = H of the electric-flux")
print("         distribution.  Test: <Z_l> and all X-off-diagonal weight.")
gg = ZNGauge("cube_Q3",*cube(),) if False else ZNGauge("cube_Q3",cube()[0],cube()[1],2)
p,_,_ = gg.ground(2.0,2.0); P = gg.full_vector(p)
r = rdm(P,gg.L,2,[0])
# in the Z basis rho = [[a,b],[b,c]]; X-diagonal means a==c  (i.e. <Z>=0)
print(f"         rho_link(Z basis) = {r.ravel()}   <Z_l> = {r[0,0]-r[1,1]:+.3e}  (0 => X-diagonal)")
Xb = np.array([[1,1],[1,-1]])/np.sqrt(2)
print(f"         rho_link(X basis) off-diagonal = {abs((Xb@r@Xb)[0,1]):.3e}")
print("    (3b) X_l = X(C_i)^dag EXACTLY on the physical sector for every cut C_i through l, at every")
print("         coupling.  Test: <X_l X(C_i)> = 1 to machine precision.")
def expect_X_product(gg,Psi,links):
    N,L = gg.N, gg.L; T = Psi.reshape([N]*L); Tv = T
    for e in links: Tv = np.roll(Tv, 1, axis=e)
    return float((T.ravel()*Tv.ravel()).sum())
cuts,d = level_cuts(gg.V,gg.edges,0)
print(f"         cube_Q3, link 0: d = dist_(G-l)(u,v) = {d}; cuts = {cuts}")
for i,C in enumerate(cuts):
    print(f"         <X_0 * X(C_{i+1})> = {expect_X_product(gg,P,[0]+C):+.12f}")
print("    (3c) girth-through-link and d for every carrier in the ladder:")
for nm,(V_,E_) in LADDER:
    gs=[]; ds=[]
    for l in range(len(E_)):
        gi,di = girth_through(V_,E_,l)
        if gi is not None: gs.append(gi); ds.append(di)
    print(f"         {nm:<20} V={V_:<3} L={len(E_):<3} C={len(E_)-V_+1:<3} "
          f"girth_through(min,max)=({min(gs)},{max(gs)})  d(min,max)=({min(ds)},{max(ds)})")
print("\nDONE 00.")
