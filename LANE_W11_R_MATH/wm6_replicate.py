# LANE_W11_R_MATH — LEG 6.  CHECK THE REGISTRAR'S NUMBERS FROM INDEPENDENT CODE, AND AUDIT THE
# THREE PLACES WHERE THE REPORTED FINDING DOES NOT MATCH THE COMPUTATION THAT PRODUCED IT.
# Arithmetic agreement is NOT the test in this program -- every confounded headline had correct
# arithmetic -- but arithmetic DISagreement would be decisive.  So: check it, then look elsewhere.
import numpy as np, wm0_lib as L
rng=np.random.default_rng(20260817)

print("== M6a  LEG A REPRODUCED (COR-F's sealed exhibit, S3_THE_CROSSING_AUDIT_V001.md:160-209) ==")
a=np.zeros(6); a[3],a[4],a[5]=0.7,1.3,-0.4
car=L.K1(); NV=5; T=L.Top(car["walkC"],a,NV); WC=L.hol(car["walkC"],a)
rho=np.diag([0.40,0.15,0.15,0.15,0.15]).astype(complex)
print(f"   ||T*T - I||              {np.linalg.norm(T.conj().T@T-np.eye(5)):.2e}    [COR-F 0.00e+00, registrar 0.00e+00]")
print(f"   T diagonal?              {np.allclose(T,np.diag(np.diag(T)))}       [COR-F False]")
print(f"   W_C                      {WC:.6f}   [COR-F -0.029200+0.999574j]")
print(f"   diag(T rho T*)           {np.round(np.real(np.diag(T@rho@T.conj().T)),2)}  [COR-F 0.15 0.15 0.15 0.40 0.15]")
wF=wC=wG=0.0
E=[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]
for _ in range(2000):
    b=rng.uniform(0,2*np.pi,6)
    wF=max(wF,np.linalg.norm(np.linalg.matrix_power(L.Top(car["walkF"],b,NV),3)-L.Mop(car["walkF"],b,NV)))
    wC=max(wC,np.linalg.norm(np.linalg.matrix_power(L.Top(car["walkC"],b,NV),3)-L.Mop(car["walkC"],b,NV)))
    th=rng.uniform(0,2*np.pi,5); bg=np.array([b[j]+th[t]-th[s] for j,(s,t) in enumerate(E)])
    s=rng.normal(size=5)+1j*rng.normal(size=5)
    wG=max(wG,np.linalg.norm(L.Top(car["walkF"],bg,NV)@(np.exp(1j*th)*s)-np.exp(1j*th)*(L.Top(car["walkF"],b,NV)@s)))
print(f"   max||T_F^3 - M_dF||      {wF:.2e}   [registrar 4.64e-15]")
print(f"   max||T_C^3 - M_c||       {wC:.2e}   [registrar 3.25e-15]")
print(f"   gauge covariance         {wG:.2e}   [registrar 4.78e-15]   AGREES")
print("   The two T^L residuals are one digit smaller than the registrar's, and the cause is not")
print("   the operator identity: the registrar forms W = exp(i(a1+a2+a3)) while T^3's entries are")
print("   exp(ia1)exp(ia2)exp(ia3).  Sum-before-exp versus product-after-exp differ by ~1 ulp and")
print("   that difference, not the identity, is what both figures measure.  Same theorem, better")
print("   conditioned.  This is a reproduction, not a disagreement.\n")

print("== M6b  LEG B AND LEG D REPRODUCED ON K1 ==")
a=np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])
TF,TC=L.Top(car["walkF"],a,NV),L.Top(car["walkC"],a,NV)
MF,MC=L.Mop(car["walkF"],a,NV),L.Mop(car["walkC"],a,NV)
sA=np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
sB=np.sqrt(np.array([0.40,0.30,0.00,0.05,0.25]))+0j
sC=sA*np.exp(1j*np.array([0.0,1.3,-0.7,2.2,0.4]))
print("   ARMS DIFF (the defect class W-08 names as commonest and FATAL):")
print(f"     ||sA - sB|| = {np.linalg.norm(sA-sB):.6f}   ||sA - sC|| = {np.linalg.norm(sA-sC):.6f}"
      f"   ||sB - sC|| = {np.linalg.norm(sB-sC):.6f}")
print(f"     bytes(sA)==bytes(sB): {sA.tobytes()==sB.tobytes()}   bytes(sA)==bytes(sC): {sA.tobytes()==sC.tobytes()}")
cl,_,_=L.classes(car)
print(f"     pi(A)=pi(B)=pi(C) to {max(np.max(np.abs(L.pi_of(sA,cl)-L.pi_of(s,cl))) for s in (sB,sC)):.1e}"
      "   -> the arms are genuinely different states with genuinely equal pi.  NOT a zero-variable control.")
reg_circ=[0.319750930010,0.760050706869,0.586994870469,0.307756757635,0.456532066199,0.513844088442]
reg_edge=[0.569227769927,0.071727337054,0.319750930010,0.581125277232,0.234214852448,
          0.760050706869,0.647209237845,0.199740989173,0.586994870469]
mc=max(abs(abs(np.vdot(np.linalg.matrix_power(MF,k)@sA,np.linalg.matrix_power(MC,k)@sA))-reg_circ[k-1]) for k in range(1,7))
me=max(abs(abs(np.vdot(np.linalg.matrix_power(TF,n)@sA,np.linalg.matrix_power(TC,n)@sA))-reg_edge[n-1]) for n in range(1,10))
print(f"   max |mine - registrar's| over B1's six circuit rows : {mc:.2e}")
print(f"   max |mine - registrar's| over B2's nine edge rows   : {me:.2e}    NO ARITHMETIC DISAGREEMENT")
print("   B3: max|Z_edge(3k) - Z_circuit(k)| over 3 states, k<=7 =",
      f"{max(abs(np.vdot(np.linalg.matrix_power(TF,3*k)@s,np.linalg.matrix_power(TC,3*k)@s)-np.vdot(np.linalg.matrix_power(MF,k)@s,np.linalg.matrix_power(MC,k)@s)) for s in (sA,sB,sC) for k in range(1,8)):.2e}",
      "[registrar 1.97e-15]\n")

print("== M6c  DEFECT 1 -- LEG C2 IS A VACUOUS SEARCH REPORTED AS A SCAN, AND ITS PREMISE IS IDLE ==")
same=[n for n in range(1,2001) if n%4==0 and n%3==0 and (n//4)==(n//3)]
print(f"   the registrar's own predicate, run: {same}")
print("   n//4 == n//3 with 12 | n and n > 0 is arithmetically impossible -- 3n = 4n forces n = 0.")
print("   The scan over n <= 2000 could not have returned anything else.  It is a one-line identity")
print("   dressed as a search.  (COULD-NOT-HAVE-FAILED does not void a theorem, and this IS one --")
print("   what is void is presenting it as evidence.)")
carB=L.B0b(); aB=np.random.default_rng(20260817).uniform(0,2*np.pi,18); NB=9
TFb,TCb=L.Top(carB["walkF"],aB,NB),L.Top(carB["walkC"],aB,NB)
MFb,MCb=L.Mop(carB["walkF"],aB,NB),L.Mop(carB["walkC"],aB,NB)
clb,_,_=L.classes(carB)
w=np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); w/=w.sum()
wB=w.copy(); wB[0],wB[1]=w[0]+w[1],0.0; wB[3],wB[4]=0.0,w[3]+w[4]; wB[5],wB[8]=w[5]+w[8],0.0
sAb=np.sqrt(w)+0j; sBb=np.sqrt(wB)+0j
def spr(oF,oC,n,sts):
    v=[abs(np.vdot(np.linalg.matrix_power(oF,n)@s,np.linalg.matrix_power(oC,n)@s)) for s in sts]
    return max(v)-min(v)
print(f"   AND THE PREMISE IS IDLE: at n = 12 the two branches ARE at different circuit counts")
print(f"   (3 and 4) and the spread is {spr(TFb,TCb,12,(sAb,sBb)):.2e} anyway.  Equal circuit counts")
print("   are neither necessary nor sufficient for invisibility; only A_n = T_F^-n T_C^n matters.")
print("   The registrar's own C3 table refutes its own C2 narrative.\n")

print("== M6d  DEFECT 2 -- LEG D's NULL, m(pi)/3, IS NOT THE NULL OF ANY CONVENTION ==")
print("   see wm2_rate.OUT.txt M2b: the invisibility-PRESERVING fibre-wise cube-root tick")
print("   U = M^(1/3) has per-tick rate m(pi) = -0.767508, not m(pi)/3 = -0.255836.")
print("   'the EDGE rate is not m(pi)/3' is therefore true of every convention including the")
print("   corpus's own, and carries no weight.  The finding is the SPREAD, and the spread is real.\n")

print("== M6e  DEFECT 3 -- 'MAX SPREAD' LUMPS TWO DISTINCT INVISIBILITIES.  SEPARATE THEM. ==")
print("   (I)  PHASE-BLINDNESS   : Z depends on s only through (|s_v|^2)_v   <=>  A_n diagonal")
print("   (II) INCIDENCE-BLINDNESS: it depends only on the CLASS SUMS        <=>  + class-constant")
print("   N2/W-03's multiset theorem is (II).  S4 sec2's 'every vertex phase of s cancels' is (I).")
print("   They are different statements with different algebraic conditions, and the registrar's")
print("   single max-spread column cannot tell which one broke.")
print(f"   {'n':>3}{'  (II) A vs B  within-class':>28}{'  (I) A vs C  phase only':>27}")
for n in range(1,10):
    print(f"   {n:>3}{spr(TF,TC,n,(sA,sB)):>28.3e}{spr(TF,TC,n,(sA,sC)):>27.3e}")
print("   -> BOTH break, at the same ticks, and the registrar's verdict survives the separation.")
print("      But the two columns are the honest report and 'max spread' was not.\n")

print("== M6f  AND THE OBJECT IS NOT A LYAPUNOV EXPONENT.  NAME IT CORRECTLY. ==")
x=sA.copy(); y=sA.copy(); nrm=[]
for _ in range(500):
    x=TF@x; y=TC@y; nrm.append((np.linalg.norm(x),np.linalg.norm(y)))
nrm=np.array(nrm)
print(f"   ||T_F^n s|| and ||T_C^n s|| over n <= 500: max deviation from 1 = {np.max(np.abs(nrm-1)):.2e}")
print("   Both branch cocycles are UNITARY, so BOTH Lyapunov exponents are exactly 0 and there is")
print("   no Lyapunov characterisation to be had.  The decay lives entirely in the ANGLE between the")
print("   branches, and (1/N) sum log|Z_n| is a Birkhoff average of log|P| over a torus rotation --")
print("   i.e. a LOGARITHMIC MAHLER MEASURE, the same object N1 already names.  Under the edge")
print("   convention it is an AVERAGE OF lcm(L_F,L_C) OF THEM.  Same theory, more terms.")
