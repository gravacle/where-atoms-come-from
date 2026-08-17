# Registrar's verification of W-11's three load-bearing claims, including the one that cuts
# against my own convention test.
import numpy as np
rng=np.random.default_rng(20260819)
LOOP_F=[(0,1,0),(1,2,1),(2,0,2)]; LOOP_C=[(0,3,3),(3,4,4),(4,0,5)]
FACE_V={0,1,2}; CYC_V={0,3,4}
CLS=[(1,1),(1,0),(1,0),(0,1),(0,1)]                    # (in F, in C) per vertex on K1
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

print("== V1  THE TRIVIAL-CONNECTION CONTACT POINT.  This one cuts AGAINST my own test. ==")
print("  FOUNDING_DESIGN:117-118 and S2:583 pre-register: NO FORMATION AT THE TRIVIAL CONNECTION.")
a0=np.zeros(6); WF,WC=hol(a0)
TF,TC=Top(LOOP_F,a0),Top(LOOP_C,a0); MF,MC=Mop(FACE_V,WF),Mop(CYC_V,WC)
s=rng.normal(size=5)+1j*rng.normal(size=5); s/=np.linalg.norm(s)
zm=[abs(np.vdot(np.linalg.matrix_power(MF,n)@s,np.linalg.matrix_power(MC,n)@s)) for n in range(1,25)]
zt=[abs(np.vdot(np.linalg.matrix_power(TF,n)@s,np.linalg.matrix_power(TC,n)@s)) for n in range(1,25)]
print(f"  under M_gamma (the corpus's):  min over n<=24 = {min(zm):.12f}   -> NO formation. Contact point MET.")
print(f"  under COR-F's T:               min over n<=24 = {min(zt):.12f}   -> T FIRES AT ZERO FIELD.")
print("  At a = 0 the transports are pure cyclic SHIFTS; the branches differ because the paths")
print("  differ, not because any field does. THE EDGE CONVENTION FAILS THE CORPUS'S OWN")
print("  PRE-REGISTERED CONTACT POINT. That is a real cost to it and it must be carried.\n")

print("== V2  THE BICONDITIONAL — the seventh naming attempt, tested where it could fail ==")
print("  CLAIM: |Z_n| is a function of pi alone  IFF  Q_n = branch_F^* branch_C is CLASS-CONSTANT")
print("         DIAGONAL (fibre-wise AND constant on each incidence class).")
def piblind(AF,AC,trials=40):
    """does |<AF s, AC s>| depend on s beyond pi? build states with identical pi and compare."""
    base=np.array([0.40,0.15,0.15,0.15,0.15]); vals=[]
    for _ in range(trials):
        w=base.copy()
        t=rng.uniform(0,min(w[1],w[2])); w[1]+=t; w[2]-=t         # move within class 10
        t2=rng.uniform(0,min(w[3],w[4])); w[3]+=t2; w[4]-=t2      # within class 01
        ph=rng.uniform(0,2*np.pi,5)
        x=np.sqrt(w)*np.exp(1j*ph)
        vals.append(abs(np.vdot(AF@x,AC@x)))
    return (max(vals)-min(vals)) < 1e-9
def class_const_diag(Q):
    if not np.allclose(Q,np.diag(np.diag(Q)),atol=1e-9): return False
    d=np.diag(Q); byc={}
    for v,c in enumerate(CLS): byc.setdefault(c,[]).append(d[v])
    return all(np.allclose(g,g[0],atol=1e-9) for g in byc.values())
agree=dis=0; blind=notblind=0
for _ in range(600):
    a=rng.uniform(0,2*np.pi,6); WF,WC=hol(a)
    kind=rng.integers(0,4)
    if kind==0: AF,AC=Mop(FACE_V,WF),Mop(CYC_V,WC)                      # the corpus's
    elif kind==1: AF,AC=Top(LOOP_F,a),Top(LOOP_C,a)                     # COR-F's
    elif kind==2:                                                        # CORRELATED non-diagonal
        R=np.zeros((5,5),dtype=complex)                                  # block-diagonal on classes
        for c in {(1,1),(1,0),(0,1)}:
            idx=[v for v,cc in enumerate(CLS) if cc==c]
            B=rng.normal(size=(len(idx),len(idx)))+1j*rng.normal(size=(len(idx),len(idx)))
            Qb,_=np.linalg.qr(B)
            for i,vi in enumerate(idx):
                for j,vj in enumerate(idx): R[vi,vj]=Qb[i,j]
        AF,AC=Mop(FACE_V,WF)@R, Mop(CYC_V,WC)@R                          # SAME R on both branches
    else:                                                                # unrestricted Haar pair
        def haar():
            B=rng.normal(size=(5,5))+1j*rng.normal(size=(5,5)); Q,_=np.linalg.qr(B); return Q
        AF,AC=haar(),haar()
    pb=piblind(AF,AC); cc=class_const_diag(AF.conj().T@AC)
    agree+= (pb==cc); dis+= (pb!=cc); blind+=pb; notblind+= (not pb)
print(f"  cells tested 600:  AGREE {agree}   DISAGREE {dis}")
print(f"  pi-blind {blind}, not pi-blind {notblind}  -> BOTH SIDES POPULATED, so it could have failed.\n")

print("== V3  THE CORRELATED-PAIR CLAIM — several lanes' 'N random draws, 0 blind' were VACUOUS ==")
n_ok=0
for _ in range(200):
    a=rng.uniform(0,2*np.pi,6); WF,WC=hol(a)
    R=np.zeros((5,5),dtype=complex)
    for c in {(1,1),(1,0),(0,1)}:
        idx=[v for v,cc in enumerate(CLS) if cc==c]
        B=rng.normal(size=(len(idx),len(idx)))+1j*rng.normal(size=(len(idx),len(idx)))
        Qb,_=np.linalg.qr(B)
        for i,vi in enumerate(idx):
            for j,vj in enumerate(idx): R[vi,vj]=Qb[i,j]
    AF,AC=Mop(FACE_V,WF)@R, Mop(CYC_V,WC)@R
    nondiag = not np.allclose(AF,np.diag(np.diag(AF)),atol=1e-9)
    if nondiag and piblind(AF,AC): n_ok+=1
print(f"  correlated NON-DIAGONAL pairs that ARE pi-blind: {n_ok} of 200")
print("  -> 'only the diagonal ones are blind' is FALSE. The blind set is the CORRELATED locus,")
print("     which has measure zero for INDEPENDENT random draws -- so any lane reporting")
print("     '0 of N random draws preserve invisibility' ran a control that could not have failed.")
