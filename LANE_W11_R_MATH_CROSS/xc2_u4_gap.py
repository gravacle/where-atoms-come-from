# LANE_W11_R_MATH_CROSS — LEG 2.  CLOSE THE GAP THE LANE UNDER TEST DECLARED IN ITS OWN
# STRONGEST RESULT -- AND IT CLOSES AGAINST THE LANE.
#
# LANE_W11_R_MATH self_flag #6: "MY THEOREM U4 IS PROVED ON K1 AND B0b, NOT IN GENERAL.  Its last
#   step needs either |class 11| = 1 (K1) or gcd(L_F, L_C) = 1 (B0b).  A carrier with
#   |class 11| >= 2 and gcd > 1 is not covered, and I did not construct one.  That is the live
#   gap in the strongest result I report."
# Its U4 is what licenses the report's headline sentence, quoted verbatim:
#   "no transport that transports anything is invisible, and the only invisible tick is the
#    convention itself" ... "no unitary tick that moves any fibre value achieves that off the
#    lcm sublattice".
#
# I CONSTRUCT THE CARRIER.  SHARE2: two 4-cycles sharing one edge, plus a spectator vertex.
#   gamma_F = 0->1->2->3->0,  gamma_C = 0->1->4->5->0,  v6 in neither.
#   class 11 = {0,1} (size 2), class 10 = {2,3}, class 01 = {4,5}, class 00 = {6}.  ALL FOUR
#   CLASSES OCCUPIED (W-09's condition), gcd(L_F,L_C) = gcd(4,4) = 4 > 1.  Exactly the gap.
#
# THE COUNTEREXAMPLE, exhibited below and verified:
#   U_F = V on class 11 (+) mu on class 10,  U_C = c V on class 11 (+) nu on class 01,
#   with V any 4th root of W_F * I_2 (NON-DIAGONAL ones exist), mu^4 = W_F, nu^4 = W_C,
#   c^4 = W_C / W_F.  Then A_n = U_F^-n U_C^n = diag(1, mu^-n, nu^n, c^n) by CLASS -- a
#   class-function multiplication operator at EVERY n -- so Z_n is pi-determined at every n,
#   while U_F is non-diagonal and MOVES fibre value between the fibres at v0 and v1.
import numpy as np, xc0_lib as X
rng=np.random.default_rng(20260817)
np.set_printoptions(linewidth=200)

car=X.SHARE2(); NV=car["NV"]
cl,F,C=X.classes(car)
LF,LC=len(car["walkF"]),len(car["walkC"])
print("== X2a  THE CARRIER THE LANE SAID IT DID NOT CONSTRUCT ==")
print(f"   V(gamma_F) = {sorted(F)}   V(gamma_C) = {sorted(C)}   L_F = {LF}  L_C = {LC}"
      f"   gcd = {np.gcd(LF,LC)}   lcm = {np.lcm(LF,LC)}")
print(f"   class vector (0=00,1=10,2=01,3=11) = {cl}   |class 11| = {int((cl==3).sum())}"
      f"   all four classes occupied: {sorted(set(cl.tolist()))==[0,1,2,3]}")
a=rng.uniform(0,2*np.pi,car["NE"])
TF,TC=X.Top(car["walkF"],a,NV),X.Top(car["walkC"],a,NV)
MF,MC=X.Mop(car["walkF"],a,NV),X.Mop(car["walkC"],a,NV)
WF,WC=X.hol(car["walkF"],a),X.hol(car["walkC"],a)
print(f"   COR-F's own edge tick on this carrier:  ||T_F^4 - M_F|| = "
      f"{np.linalg.norm(np.linalg.matrix_power(TF,4)-MF):.2e}"
      f"   ||T_C^4 - M_C|| = {np.linalg.norm(np.linalg.matrix_power(TC,4)-MC):.2e}")

pi=np.array([0.10,0.30,0.25,0.35])
states=X.same_pi_states(cl,pi,np.random.default_rng(11),40)
dmin=min(np.linalg.norm(states[i]-states[j]) for i in range(10) for j in range(i+1,10))
pimax=max(np.max(np.abs(X.pi_of(states[0],cl)-X.pi_of(s,cl))) for s in states[1:])
print(f"   ARMS DIFF over the 40-state ensemble: min ||s_i - s_j|| = {dmin:.4f}"
      f" ; pi equal to {pimax:.1e}  -> NOT a zero-variable control")

print("\n== X2b  COR-F's TICK ON SHARE2: visible, exactly as the lane's C1/C2/C3 predict ==")
print(f"   {'n':>3}{'  EDGE spread':>16}{'  CIRCUIT spread':>18}")
for n in (1,2,3,4,5,6,8,12):
    print(f"   {n:>3}{X.spread(TF,TC,n,states):>16.3e}{X.spread(MF,MC,n,states):>18.3e}")
print("   -> lcm(4,4) = 4, and the edge tick is invisible exactly on 4Z.  Nothing new so far.")

print("\n== X2c  THE COUNTEREXAMPLE TO THEOREM U4.  A NON-DIAGONAL TICK THAT MOVES FIBRE VALUE ==")
print("        AND IS pi-INVISIBLE AT EVERY SINGLE TICK.")
def blockop(vs,B):
    U=np.eye(NV,dtype=complex)
    for i,u in enumerate(vs):
        for j,v in enumerate(vs): U[u,v]=B[i,j]
    return U
def haar(n,r):
    z=(r.normal(size=(n,n))+1j*r.normal(size=(n,n)))/np.sqrt(2)
    q,rr=np.linalg.qr(z); return q@np.diag(np.diag(rr)/np.abs(np.diag(rr)))
cls11=sorted(np.where(cl==3)[0]); cls10=sorted(np.where(cl==1)[0]); cls01=sorted(np.where(cl==2)[0])
print(f"   class 11 = {cls11}   class 10 = {cls10}   class 01 = {cls01}")
worst_sp=0.0; worst_off=0.0; nbuilt=0
for trial in range(200):
    r=np.random.default_rng(1000+trial)
    Q=haar(2,r); ka,kb=r.integers(0,4,2)
    while ka==kb: ka,kb=r.integers(0,4,2)                       # non-diagonal 4th root of W_F I_2
    q4=np.exp(1j*np.angle(WF)/4)
    V=q4*(Q@np.diag([1j**int(ka),1j**int(kb)])@Q.conj().T)
    mu=q4*1j**int(r.integers(0,4))                              # mu^4 = W_F
    nu=np.exp(1j*np.angle(WC)/4)*1j**int(r.integers(0,4))       # nu^4 = W_C
    cc=np.exp(1j*(np.angle(WC)-np.angle(WF))/4)*1j**int(r.integers(0,4))   # c^4 = W_C/W_F
    UF=np.eye(NV,dtype=complex); UC=np.eye(NV,dtype=complex)
    UF[np.ix_(cls11,cls11)]=V
    for v in cls10: UF[v,v]=mu
    UC[np.ix_(cls11,cls11)]=cc*V
    for v in cls01: UC[v,v]=nu
    okF=np.linalg.norm(np.linalg.matrix_power(UF,4)-MF); okC=np.linalg.norm(np.linalg.matrix_power(UC,4)-MC)
    off=np.linalg.norm(UF-np.diag(np.diag(UF)))
    if max(okF,okC)>1e-9 or off<1e-6: continue
    nbuilt+=1
    sp=max(X.spread(UF,UC,n,states) for n in range(1,13))
    worst_sp=max(worst_sp,sp); worst_off=max(worst_off,off)
    if nbuilt==1:
        print(f"   FIRST WITNESS.  ||U_F^4 - M_F|| = {okF:.2e}   ||U_C^4 - M_C|| = {okC:.2e}")
        print(f"     U_F is NOT diagonal:  ||U_F - diag(U_F)|| = {off:.4f}")
        print(f"     it MOVES fibre value: |<e_1, U_F e_0>| = {abs(UF[cls11[1],cls11[0]]):.6f}"
              f"   (the amplitude at v{cls11[0]} is transported to v{cls11[1]})")
        print(f"     U_F block on class 11 =\n{np.round(UF[np.ix_(cls11,cls11)],6)}")
        A1=np.linalg.inv(UF)@UC
        print(f"     A_1 = U_F^-1 U_C off-diagonal norm = "
              f"{np.linalg.norm(A1-np.diag(np.diag(A1))):.2e}   diag = {np.round(np.diag(A1),6)}")
        print(f"     max spread over 40 same-pi states, n = 1..12 = {sp:.3e}")
print(f"   {nbuilt} witnesses built; WORST spread over all of them, n = 1..12 = {worst_sp:.3e}"
      f" ; LARGEST off-diagonal mass {worst_off:.4f}")
print("   -> pi-INVISIBLE AT EVERY TICK, NON-DIAGONAL, AND IT TRANSPORTS.  THEOREM U4 IS FALSE")
print("      OUTSIDE ITS TWO CARRIERS, exactly in the case the lane flagged and did not build.")

print("\n== X2d  WHAT SURVIVES THE COUNTEREXAMPLE, AND WHAT DOES NOT ==")
A=[np.linalg.inv(np.linalg.matrix_power(TF,n))@np.linalg.matrix_power(TC,n) for n in (1,2,3,4)]
print("   (i)  THE LANE'S THEOREM C1 SURVIVES: my witness is invisible because its RELATIVE")
print("        operator A_n is a class-function multiplication operator, which is C1 exactly.")
print("   (ii) WHAT FALLS is U4's corollary 'therefore U_F and U_C are DIAGONAL', and with it the")
print("        report's sentence 'no transport that transports anything is invisible'.  The")
print("        correct statement is about A_n, never about the individual ticks.")
print("   (iii) IT DOES NOT RESTORE READING A ON EITHER CORPUS CARRIER: K1 has |class 11| = 1 and")
print("        B0b has gcd(4,3) = 1, so on both of them U4's conclusion still holds and the")
print("        registrar's verdict is untouched.  And the mechanism is honest about itself: the")
print("        witness is invisible because U_C|_11 = c U_F|_11 makes the two branches move IN")
print("        STEP on the shared class.  It is the stipulation again, wearing a 2x2 block.")
print("        WHAT IT IS NOT is a free numerical block: X2e builds it from the connection's own")
print("        edge transport, so it is gauge-covariant, which is the corpus's own admissibility")
print("        test and the one LEG A uses to admit COR-F's T.")
print("   (iv) BUT IT IS A FOURTH FAMILY, and the lane reported the family count as ONE")
print("        ('Exactly one family restores invisibility: U = zeta M^(1/L)').  That census is")
print("        complete only under U4's unstated hypothesis.")

print("\n== X2e  AND THE WITNESS CAN BE MADE GAUGE-COVARIANT, WHICH IS THE CORPUS'S OWN TEST ==")
print("   LEG A of the lane under test admits COR-F's T because it is unitary, gauge-COVARIANT and")
print("   T^L = M.  My witness passes the same three.  Construction: on the shared edge e0 = (v0,v1)")
print("   put  V = Lam V0 Lam*  with Lam = diag(1, U_e0) and V0 ANY non-diagonal 4th root of W_F I_2.")
print("   Then V_{10} = U_e0 (V0)_{10} and V_{01} = conj(U_e0) (V0)_{01}: every off-diagonal entry is")
print("   a multiple of the connection's own parallel transport along a real edge of the complex, so")
print("   V(a + dtheta) = G V(a) G*.  V^4 = Lam V0^4 Lam* = W_F I.  Verified over 2000 gauge maps:")
E=car["E"]
def gauge_defect(build,ntrial=2000):
    worst=0.0
    r=np.random.default_rng(20260817)
    for _ in range(ntrial):
        b=r.uniform(0,2*np.pi,car["NE"]); th=r.uniform(0,2*np.pi,NV)
        bg=np.array([b[j]+th[t]-th[s] for j,(s,t) in enumerate(E)])
        s=r.normal(size=NV)+1j*r.normal(size=NV); s/=np.linalg.norm(s)
        G=np.exp(1j*th)
        worst=max(worst,np.linalg.norm(build(bg)@(G*s)-G*(build(b)@s)))
    return worst
r0=np.random.default_rng(4242)
Q=haar(2,r0); ka,kb=0,2
V0=Q@np.diag([1j**ka,1j**kb])@Q.conj().T
kmu,knu,kc=1,3,2
def make(bb):
    WFb,WCb=X.hol(car["walkF"],bb),X.hol(car["walkC"],bb)
    Lam=np.diag([1.0,np.exp(1j*bb[0])])
    V=np.exp(1j*np.angle(WFb)/4)*(Lam@V0@Lam.conj().T)
    U=np.eye(NV,dtype=complex); U[np.ix_(cls11,cls11)]=V
    for v in cls10: U[v,v]=np.exp(1j*np.angle(WFb)/4)*1j**kmu
    return U
def makeC(bb):
    WFb,WCb=X.hol(car["walkF"],bb),X.hol(car["walkC"],bb)
    Lam=np.diag([1.0,np.exp(1j*bb[0])])
    V=np.exp(1j*np.angle(WFb)/4)*(Lam@V0@Lam.conj().T)
    cc=np.exp(1j*(np.angle(WCb)-np.angle(WFb))/4)*1j**kc
    U=np.eye(NV,dtype=complex); U[np.ix_(cls11,cls11)]=cc*V
    for v in cls01: U[v,v]=np.exp(1j*np.angle(WCb)/4)*1j**knu
    return U
gF=gauge_defect(make); gC=gauge_defect(makeC)
UF=make(a); UC=makeC(a)
print(f"     max || U_F(g.a) (g.s) - g.(U_F(a) s) || over 2000 gauge transforms = {gF:.2e}")
print(f"     max || U_C(g.a) (g.s) - g.(U_C(a) s) || over 2000 gauge transforms = {gC:.2e}")
print(f"     || U_F* U_F - I || = {np.linalg.norm(UF.conj().T@UF-np.eye(NV)):.2e}"
      f"   || U_F^4 - M_F || = {np.linalg.norm(np.linalg.matrix_power(UF,4)-MF):.2e}")
print(f"     || U_C* U_C - I || = {np.linalg.norm(UC.conj().T@UC-np.eye(NV)):.2e}"
      f"   || U_C^4 - M_C || = {np.linalg.norm(np.linalg.matrix_power(UC,4)-MC):.2e}")
print(f"     U_F non-diagonal: ||U_F - diag|| = {np.linalg.norm(UF-np.diag(np.diag(UF))):.4f}"
      f"   |<e_1,U_F e_0>| = {abs(UF[cls11[1],cls11[0]]):.6f}")
print(f"     max pi-spread over 40 same-pi states, n = 1..24 = "
      f"{max(X.spread(UF,UC,n,states) for n in range(1,25)):.3e}")
print("   -> UNITARY, GAUGE-COVARIANT, U^L = M_gamma, NON-DIAGONAL, MOVES FIBRE VALUE ALONG A REAL")
print("      EDGE -- and pi-INVISIBLE AT EVERY TICK.  It passes every admissibility test the lane")
print("      under test applies to COR-F's own T, and it is not diagonal.")
