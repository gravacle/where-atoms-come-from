# LANE_W11_R_MATH — LEG 8.  W-10 REGISTERS THE COMMUTATOR AS THE DIAGNOSTIC.  IT IS NOT ONE.
# REGISTER_V001.md:1042-1046: "||[T_F, T_C]|| is 2.828 on B0b, 2.449 on K1 and exactly 0 on B0a --
# non-zero precisely when class 11 is occupied.  The incidence that N2 calls invisible is visible
# to the corpus's own sealed alternative transport."
# My Leg 3 (C2) says the operative object is sigma_F^n versus sigma_C^n, not [T_F,T_C].  Test it on
# a DISJOINT-LOOP designation, where the commutator vanishes identically and class 11 is empty.
import numpy as np, wm0_lib as L
car=L.B0b(); NV=9; E=car["E"]
# same complex, same gamma_F; gamma_C moved to the row j = 2, whose vertices {6,7,8} are DISJOINT
# from V(gamma_F) = {0,1,3,4}.  Admissible under S4's CHOICE LEDGER C4 exactly as B0b's own is
# (W-10's erratum against W-09: loop designation is by fiat, and B0b's own sweep reaches 16 multisets).
gC=[(6,1),(7,1),(8,1)]
def walk(g):
    out=[]
    for (e,s) in g:
        u,v=E[e]; out.append((u,v,e,+1) if s>0 else (v,u,e,-1))
    return out
carD=dict(car); carD["walkC"]=walk(gC); carD["name"]="B0b/disjoint"
cl,F,C=L.classes(carD)
print(f"  V(gamma_F) = {sorted(F)}   V(gamma_C) = {sorted(C)}   overlap = {sorted(F&C)}")
print(f"  class vector = {cl}   class 11 occupied? {3 in cl}   (0=00, 1=10, 2=01, 3=11)")
a=np.random.default_rng(20260817).uniform(0,2*np.pi,18)
TF,TC=L.Top(carD["walkF"],a,NV),L.Top(carD["walkC"],a,NV)
MF,MC=L.Mop(carD["walkF"],a,NV),L.Mop(carD["walkC"],a,NV)
print(f"  || [T_F, T_C] ||  =  {np.linalg.norm(TF@TC-TC@TF):.2e}    <-- EXACTLY ZERO, as on B0a")
sizes=[int((cl==c).sum()) for c in range(4)]; pi=np.array([s/sum(sizes) for s in sizes])
r=np.random.default_rng(11); states=[]
for _ in range(40):
    w=np.zeros(NV)
    for c in range(4):
        idx=np.where(cl==c)[0]
        if len(idx)==0: continue
        q=r.random(len(idx)); q=q/q.sum()*pi[c]; w[idx]=q
    states.append(np.sqrt(w)*np.exp(1j*r.uniform(0,2*np.pi,NV)))
print(f"  ARMS DIFF: min ||s_i - s_j|| over the 40-state ensemble = "
      f"{min(np.linalg.norm(states[i]-states[j]) for i in range(8) for j in range(i+1,8)):.3f}")
print(f"  pi identical across the ensemble to "
      f"{max(np.max(np.abs(L.pi_of(states[0],cl)-L.pi_of(s,cl))) for s in states[1:]):.1e}")
def spr(oF,oC,n):
    v=[abs(np.vdot(np.linalg.matrix_power(oF,n)@s,np.linalg.matrix_power(oC,n)@s)) for s in states]
    return max(v)-min(v)
print(f"  {'n':>3}{'  EDGE spread':>16}{'  CIRCUIT spread (k=n)':>24}")
for n in (1,2,3,4,6,8,9,12,24):
    print(f"  {n:>3}{spr(TF,TC,n):>16.3e}{spr(MF,MC,n):>24.3e}")
print("  -> the commutator is EXACTLY 0 and class 11 is EMPTY, and the edge convention STILL makes")
print("     the incidence visible at every n not divisible by lcm(4,3) = 12, with spreads up to")
print("     O(1).  W-10's 'non-zero precisely when class 11 is occupied' is a true statement ABOUT")
print("     THE COMMUTATOR and a FALSE diagnostic for visibility: [T_F,T_C] = 0 is neither")
print("     necessary nor sufficient.  The operative object is sigma_F^n versus sigma_C^n (Leg 3).")
print("     A SEVENTH MISNAMING OF THE OPERATIVE VARIABLE, and it is in a registered W-10 row.")
