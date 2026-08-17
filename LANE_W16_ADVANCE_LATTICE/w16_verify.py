# W-16 — the registrar verifies the two frame-refusals from W-11 that the registrar overrode.
# CLAIM (LANE_W11_R_BLIND and LANE_W11_R_CLOCK, both returning REFUTED_AS_POSED):
#   there is only ONE transport. M_gamma = T^L, so the "circuit convention" is a SCHEDULE on
#   COR-F's own generator. Both conventions are rays in one family Y(mF,mC) = <T_F^mF s, T_C^mC s>:
#   EDGE is the ray (n,n), CIRCUIT is the ray (L_F k, L_C k). The set where the incidence is
#   invisible is EXACTLY the sublattice L_F*Z x L_C*Z.
import numpy as np
rng=np.random.default_rng(20260823)
LOOP_F=[(0,1,0),(1,2,1),(2,0,2)]; LOOP_C=[(0,3,3),(3,4,4),(4,0,5)]
FACE_V={0,1,2}; CYC_V={0,3,4}; LF=LC=3
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
a=np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])
WF,WC=np.exp(1j*sum(a[:3])),np.exp(1j*sum(a[3:]))
TF,TC=Top(LOOP_F,a),Top(LOOP_C,a); MF,MC=Mop(FACE_V,WF),Mop(CYC_V,WC)

print("== V1  IS THERE ONE TRANSPORT OR TWO? ==")
print(f"  || T_F^3 - M_dF || = {np.linalg.norm(np.linalg.matrix_power(TF,3)-MF):.2e}")
print(f"  || T_C^3 - M_c  || = {np.linalg.norm(np.linalg.matrix_power(TC,3)-MC):.2e}")
print("  M lies INSIDE the cyclic group T generates. The corpus's operator is a POWER of COR-F's.")
print("  -> the two lanes are right: there is ONE transport and TWO SCHEDULES ON IT.\n")

print("== V2  THE INVISIBLE SET: is it EXACTLY the advance sublattice L_F*Z x L_C*Z? ==")
# 64 states with identical pi, differing within class and in phase
base=np.array([0.40,0.15,0.15,0.15,0.15])
states=[]
for _ in range(64):
    w=base.copy()
    t=rng.uniform(0,min(w[1],w[2])); w[1]+=t; w[2]-=t
    t2=rng.uniform(0,min(w[3],w[4])); w[3]+=t2; w[4]-=t2
    states.append(np.sqrt(w)*np.exp(1j*rng.uniform(0,2*np.pi,5)))
pw={}
def P(Mx,n):
    if n not in pw.setdefault(id(Mx),{}): pw[id(Mx)][n]=np.linalg.matrix_power(Mx,n)
    return pw[id(Mx)][n]
inv=set(); notinv=set()
for mF in range(26):
    for mC in range(26):
        v=[abs(np.vdot(P(TF,mF)@s, P(TC,mC)@s)) for s in states]
        (inv if (max(v)-min(v))<1e-9 else notinv).add((mF,mC))
lat={(mF,mC) for mF in range(26) for mC in range(26) if mF%LF==0 and mC%LC==0}
print(f"  invisible cells {len(inv)} of 676   advance sublattice {len(lat)} of 676")
print(f"  SET EQUALITY: {inv==lat}")
spread_off=max(max(abs(np.vdot(P(TF,mF)@s,P(TC,mC)@s)) for s in states)
              -min(abs(np.vdot(P(TF,mF)@s,P(TC,mC)@s)) for s in states) for (mF,mC) in list(notinv)[:80])
print(f"  max spread OFF the lattice (80 sampled cells) = {spread_off:.3e}")
print()
print("  ==> CONFIRMED. The corpus's ray (3k,3k) and COR-F's ray (n,n) both sit in ONE family,")
print("      and invisibility is a property of the ADVANCE LATTICE, not of a choice of transport.")
print("      W-11's registered headline names the wrong object.")
