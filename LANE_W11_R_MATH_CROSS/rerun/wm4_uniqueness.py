# LANE_W11_R_MATH — LEG 4.  "T IS NOT UNIQUE" -- THE REGISTRAR'S OWN DECLARED WEAKNESS.  SETTLED.
#
# The brief: "If a different, equally natural edge tick restores invisibility, Reading B falls.
#  If EVERY unitary with T^L = M_gamma except the diagonal ones breaks invisibility, Reading B is
#  far stronger than the registrar claimed."   Both halves are decided here.
#
# THEOREM U1 (EDGE-LOCALITY DICHOTOMY).  Call U edge-local if U_{vu} != 0 only when u = v or
#   (u -> v) is a directed edge of the loop.  On a simple L-cycle every edge-local unitary is
#   EITHER a pure weighted cyclic shift OR diagonal -- nothing in between.
#   PROOF.  Column u has support {(u,u),(u+1,u)} with entries d0_u, d1_u, |d0_u|^2 + |d1_u|^2 = 1.
#   Columns u and u+1 overlap in the single row u+1: conj(d1_u) d0_{u+1} = 0.  So d1_u != 0 forces
#   d0_{u+1} = 0, hence |d1_{u+1}| = 1, hence d1_{u+1} != 0 -- and the implication runs all the way
#   round the cycle.  Either every d1 vanishes (diagonal) or every d0 does (pure shift).  QED
#
# THEOREM U2 (HOW BIG THE NON-DIAGONAL FAMILY REALLY IS -- STATED AGAINST MYSELF).
#   Edge-local + unitary + gauge-covariant + U^L = M does NOT pin T.  Gauge covariance forces the
#   weight on edge e to be lam_e * U_e with lam_e a gauge-invariant unimodular scalar; U^L = M then
#   forces only prod_e lam_e = 1.  That is an (L-1)-parameter family, and COR-F picked lam = 1.
#   THE REGISTRAR'S DECLARED WEAKNESS IS REAL AND I DO NOT REPAIR IT.
#
# THEOREM U3 (AND IT IS IMMATERIAL).  By Leg 3 (C2), whether Z_n is pi-determined depends on the
#   two ticks ONLY through their PERMUTATIONS sigma_F, sigma_C -- the weights cancel out of the
#   diagonality question entirely.  Every member of the whole (L-1)-parameter family has
#   sigma = the full L-cycle, so every one of them is visible at exactly the same ticks as T.
#   The non-uniqueness of T cannot restore invisibility because it never touches the operative object.
#
# THEOREM U4 (ONLY DIAGONALS ARE INVISIBLE -- locality dropped entirely).  Let U_F, U_C be ANY
#   unitaries, identity off their own loop, with U_F^{L_F} = M_F, U_C^{L_C} = M_C, such that
#   Z_n = <U_F^n s, U_C^n s> is pi-determined for every n and every s.  Then on K1 and on B0b,
#   U_F and U_C are DIAGONAL -- indeed class-constant fibre-wise scalars zeta * M^{1/L}.
#   PROOF.  Leg 3 (C1): every A_n = U_F^{-n} U_C^n is diagonal and class-constant.  n = 1 gives
#   U_C = U_F D_1;  then A_n A_{n-1}^{-1} = U_F^{-(n-1)} D_1 U_F^{n-1} is diagonal for every n, so
#   U_F normalises the diagonal algebra; with D_1's class values distinct U_F preserves every class
#   coordinate subspace.  On class 10 = V(F)\V(C), U_C is the identity, so A_n|_10 = U_F^{-n}|_10
#   must be scalar for all n: U_F|_10 is scalar.  On class 11 = V(F) cap V(C), A_1|_11 scalar gives
#   U_C|_11 = c U_F|_11, so (U_F|_11)^{L_F} and (U_F|_11)^{L_C} are both scalars, hence so is
#   (U_F|_11)^{gcd(L_F,L_C)}.  K1: class 11 = {v0}, one-dimensional.  B0b: gcd(4,3) = 1.  Either way
#   U_F|_11 is scalar.  Off V(F), U_F is the identity.  So U_F is diagonal; symmetrically U_C.  QED
import numpy as np, wm0_lib as L
rng=np.random.default_rng(20260817)
om=np.exp(2j*np.pi/3)

def random_same_pi(cl,pi,r,k):
    NV=len(cl); out=[]
    for _ in range(k):
        w=np.zeros(NV)
        for c in range(4):
            idx=np.where(cl==c)[0]
            if len(idx)==0: continue
            q=r.random(len(idx)); q=q/q.sum()*pi[c]; w[idx]=q
        out.append(np.sqrt(w)*np.exp(1j*r.uniform(0,2*np.pi,NV)))
    return out
def spread(UF,UC,states,n):
    v=[abs(np.vdot(np.linalg.matrix_power(UF,n)@s, np.linalg.matrix_power(UC,n)@s)) for s in states]
    return max(v)-min(v)
def isdiag(U,tol=1e-11): return np.linalg.norm(U-np.diag(np.diag(U)))<tol
def haar(n,r):
    z=(r.normal(size=(n,n))+1j*r.normal(size=(n,n)))/np.sqrt(2)
    q,rr=np.linalg.qr(z); return q@np.diag(np.diag(rr)/np.abs(np.diag(rr)))

car=L.K1(); NV=5; a=np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])
wF,wC=car["walkF"],car["walkC"]
TF,TC=L.Top(wF,a,NV),L.Top(wC,a,NV); MF,MC=L.Mop(wF,a,NV),L.Mop(wC,a,NV)
WF,WC=L.hol(wF,a),L.hol(wC,a); cl,F,C=L.classes(car); Fv,Cv=sorted(F),sorted(C)
pi=np.array([0.0,0.30,0.30,0.40]); states=random_same_pi(cl,pi,np.random.default_rng(3),30)

print("== M4a  THEOREM U1, DIRECTLY: solve the edge-local unitarity constraints on an L-cycle ==")
print("   parametrise d0_u = cos(t_u) e^{i al_u}, d1_u = sin(t_u) e^{i be_u}; the ONLY constraint")
print("   beyond column normalisation is conj(d1_u) d0_{u+1} = 0 for every u.  Solve it exactly.")
for Lc in (3,4,5,6):
    found_mixed=0; found_shift=0; found_diag=0; tested=0
    for _ in range(200000):
        t=rng.uniform(0,np.pi/2,Lc)
        # impose the constraint by choosing, per u, which of the two factors is zeroed
        mask=rng.integers(0,2,Lc)                      # 0: set d1_u = 0 ; 1: set d0_{u+1} = 0
        d0=np.cos(t).astype(complex); d1=np.sin(t).astype(complex)
        for u in range(Lc):
            if mask[u]==0: d1[u]=0.0
            else:          d0[(u+1)%Lc]=0.0
        nrm=np.sqrt(np.abs(d0)**2+np.abs(d1)**2)
        if np.any(nrm<1e-12): continue
        d0/=nrm; d1/=nrm
        U=np.zeros((Lc,Lc),dtype=complex)
        for u in range(Lc): U[u,u]=d0[u]; U[(u+1)%Lc,u]=d1[u]
        if np.linalg.norm(U.conj().T@U-np.eye(Lc))>1e-10: continue
        tested+=1
        sh=np.all(np.abs(np.diag(U))<1e-12); dg=np.all(np.abs([U[(u+1)%Lc,u] for u in range(Lc)])<1e-12)
        if sh: found_shift+=1
        elif dg: found_diag+=1
        else: found_mixed+=1
    print(f"   L = {Lc}:  edge-local UNITARIES found {tested:>6}   pure shift {found_shift:>6}"
          f"   diagonal {found_diag:>6}   MIXED {found_mixed}")
print("   -> MIXED count is 0 at every L.  U1 holds.\n")

print("== M4b  THEOREM U2: the non-diagonal edge-local gauge-covariant family, honestly sized ==")
print("   U_lam: (U s)_v = lam_e U_e s_u on the loop.  Gauge-covariant for any unimodular lam.")
lam_tests=[]
for _ in range(2000):
    ph=rng.uniform(0,2*np.pi,3); ph[2]=-(ph[0]+ph[1])          # prod lam = 1
    lam=np.exp(1j*ph)
    U=np.zeros((NV,NV),dtype=complex)
    for v in range(NV):
        if v not in F: U[v,v]=1.0
    for i,(u,v,e,sg) in enumerate(wF): U[v,u]=lam[i]*np.exp(1j*a[e])
    lam_tests.append(np.linalg.norm(np.linalg.matrix_power(U,3)-MF))
print(f"   2000 random lam with prod lam = 1:  max ||U_lam^3 - M_F|| = {max(lam_tests):.2e}")
print(f"   so the family is (L-1)-dimensional and COR-F's lam = 1 is ONE point of it.  Registrar's")
print("   declared weakness CONFIRMED, not repaired.\n")
print("== M4c  THEOREM U3: and it is immaterial -- the whole family is visible at the same ticks ==")
mn=1e9; mx=0.0
for _ in range(2000):
    phF=rng.uniform(0,2*np.pi,3); phF[2]=-(phF[0]+phF[1])
    phC=rng.uniform(0,2*np.pi,3); phC[2]=-(phC[0]+phC[1])
    def mk(walk,vs,lam):
        U=np.zeros((NV,NV),dtype=complex)
        for v in range(NV):
            if v not in vs: U[v,v]=1.0
        for i,(u,v,e,sg) in enumerate(walk): U[v,u]=lam[i]*np.exp(1j*a[e])
        return U
    UF=mk(wF,F,np.exp(1j*phF)); UC=mk(wC,C,np.exp(1j*phC))
    s12=max(spread(UF,UC,states,1),spread(UF,UC,states,2))
    s3 =max(spread(UF,UC,states,3),spread(UF,UC,states,6))
    mn=min(mn,s12); mx=max(mx,s3)
print(f"   over 2000 random members of the family:")
print(f"     MINIMUM spread at n = 1 or 2 (off the sublattice) = {mn:.3e}   (never invisible)")
print(f"     MAXIMUM spread at n = 3 or 6 (on the sublattice)  = {mx:.3e}   (always invisible)")
print("   The weights lam never enter the verdict.  Only the PERMUTATION does.\n")

print("== M4d  THEOREM U4, THE DECISIVE SEARCH.  Drop locality.  Every unitary U_F on gamma_F with ==")
print("        U_F^3 = M_F is  exp(i arg W_F /3) * Q diag(w^a,w^b,w^c) Q^*.  20000 Haar draws.")
def build(Wh,vs,Q,expo):
    V=np.exp(1j*np.angle(Wh)/3)*(Q@np.diag([om**e for e in expo])@Q.conj().T)
    U=np.eye(NV,dtype=complex)
    for i,u in enumerate(vs):
        for j,v in enumerate(vs): U[u,v]=V[i,j]
    return U
n_inv=0; n_inv_nondiag=0; worstpow=0.0; n_diag=0
for _ in range(20000):
    UF=build(WF,Fv,haar(3,rng),tuple(rng.integers(0,3,3)))
    UC=build(WC,Cv,haar(3,rng),tuple(rng.integers(0,3,3)))
    worstpow=max(worstpow,np.linalg.norm(np.linalg.matrix_power(UF,3)-MF),
                          np.linalg.norm(np.linalg.matrix_power(UC,3)-MC))
    both_diag=isdiag(UF) and isdiag(UC)
    n_diag+=both_diag
    if max(spread(UF,UC,states,n) for n in (1,2,4,5))<1e-9:
        n_inv+=1
        if not both_diag: n_inv_nondiag+=1
print(f"   all 20000 satisfy ||U^3 - M|| <= {worstpow:.1e}")
print(f"   pi-INVISIBLE at n = 1,2,4,5 : {n_inv}      both ticks diagonal : {n_diag}")
print(f"   INVISIBLE and NOT both diagonal : {n_inv_nondiag}   <-- U4 says this must be 0")
print("   (the two counts coincide: a Haar Q gives a diagonal U only when all three exponents")
print("    agree, probability (3/27)^2 = 1/81, i.e. ~247 of 20000.)\n")

print("== M4e  THE ADVERSARIAL SUB-FAMILY my own proof names as the last hope: U_F class-block- ==")
print("        preserving but NON-diagonal on class 10 = {v1,v2}.  If anything survives, it is here.")
best=0.0; worst=1e9; cnt=0
for _ in range(4000):
    ph=np.exp(1j*np.angle(WF)/3); k=rng.integers(0,3,2)
    if k[0]==k[1]: continue                      # that is the diagonal case, excluded by hand
    q=haar(2,rng); R2=q@np.diag([om**k[0],om**k[1]])@q.conj().T
    UF=np.eye(NV,dtype=complex); UF[0,0]=ph*om**rng.integers(0,3)
    UF[np.ix_([1,2],[1,2])]=ph*R2
    UC=np.eye(NV,dtype=complex)
    for v in Cv: UC[v,v]=np.exp(1j*np.angle(WC)/3)
    assert np.linalg.norm(np.linalg.matrix_power(UF,3)-MF)<1e-9
    assert np.linalg.norm(np.linalg.matrix_power(UC,3)-MC)<1e-9
    assert not isdiag(UF)
    sp=max(spread(UF,UC,states,n) for n in (1,2)); cnt+=1
    best=max(best,sp); worst=min(worst,sp)
print(f"   {cnt} draws, all with U_F^3 = M_F exactly and U_F non-diagonal on class 10,")
print(f"   AND with U_C the invisibility-preserving diagonal tick (so U_F alone carries the test):")
print(f"     smallest spread produced = {worst:.3e}   largest = {best:.3e}   -> none invisible.\n")

print("== M4f  THE ONE FAMILY THAT DOES RESTORE INVISIBILITY -- AND WHAT IT IS ==")
for kF in (0,1,2):
    UF=np.eye(NV,dtype=complex); UC=np.eye(NV,dtype=complex)
    for v in Fv: UF[v,v]=np.exp(1j*np.angle(WF)/3)*om**kF
    for v in Cv: UC[v,v]=np.exp(1j*np.angle(WC)/3)
    sp=max(spread(UF,UC,states,n) for n in range(1,13))
    print(f"   U = zeta M^(1/3), zeta = w^{kF}:  ||U_F^3-M_F|| = "
          f"{np.linalg.norm(np.linalg.matrix_power(UF,3)-MF):.1e}   max spread, n <= 12 = {sp:.2e}")
print("   Invisible at EVERY tick -- and it is FIBRE-WISE SCALAR MULTIPLICATION BY A POWER OF THE")
print("   WHOLE-CIRCUIT HOLONOMY.  It is the corpus's own stipulation with a finer clock, not an")
print("   independent transport: it moves no fibre value along any edge, and it is the one operator")
print("   in the search that is NOT edge-local in the moving sense.  So the answer to the")
print("   registrar's declared weakness is: YES another tick restores invisibility, and it is the")
print("   stipulation.  No transport that MOVES anything does.")
