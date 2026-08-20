"""ADVERSARIAL VERIFY 1.
Rebuild the toric code L=2 FROM SCRATCH with the BRIEF's Hamiltonian
    H = -sum_v A_v - sum_p B_p        (ALL 4+4 terms, not the rref'd 6 generators the lane used)
and re-decide, independently of the lane's code:
  * clauses (i)-(iv) on the logical Z's
  * the joint dimension table
  * G_W by EXPLICIT construction+verification of every flip unitary
  * #orbits by Burnside AND dim(invariants) by an EXPLICIT NULLSPACE of the averaging map
  * the code distance by EXHAUSTIVE 4^8 Pauli scan
"""
import numpy as np, itertools
np.set_printoptions(linewidth=200)

I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex)
Z=np.array([[1,0],[0,-1]],dtype=complex); Y=1j*X@Z
def kron(ms):
    M=np.array([[1]],dtype=complex)
    for m in ms: M=np.kron(M,m)
    return M
def pau(n,spec):
    d={'X':X,'Y':Y,'Z':Z,'I':I2}
    return kron([d[spec[i]] for i in range(n)])

L=2; n=2*L*L
def h(i,j): return 2*((i%L)*L+(j%L))
def v(i,j): return 2*((i%L)*L+(j%L))+1

# --- vertex stars (X) and plaquettes (Z), ALL of them (brief's H)
stars=[]; plaqs=[]
for i in range(L):
    for j in range(L):
        s=['I']*n
        for e in (h(i,j),h(i-1,j),v(i,j),v(i,j-1)): s[e]='X'
        stars.append(''.join(s))
        p=['I']*n
        for e in (h(i,j),h(i,j+1),v(i,j),v(i+1,j)): p[e]='Z'
        plaqs.append(''.join(p))
S_all = stars+plaqs
H = -sum(pau(n,s) for s in S_all)
w,V=np.linalg.eigh(H)
# eigenspaces
def eigspaces(H,tol=1e-8):
    w,V=np.linalg.eigh(H); out=[]; i=0
    while i<len(w):
        j=i
        while j+1<len(w) and abs(w[j+1]-w[i])<tol: j+=1
        out.append((w[i], V[:,i:j+1])); i=j+1
    return out
ES=eigspaces(H)
print("BRIEF's H = -sum_v A_v - sum_p B_p  (all %d terms).  eigenvalues/dims:"%len(S_all))
for e,C in ES: print("   E=%7.3f  dim=%d"%(e,C.shape[1]))
print("   #shells =",len(ES),"   (the LANE reported 7 shells at E=-6..6, i.e. a DIFFERENT H:",
      "it used only the 6 rref'd independent generators)")

# --- logical operators, DERIVED (not nominated): centraliser of the stabiliser group mod stabiliser
def sym(a,b,n): return sum(a[i]*b[n+i]+a[n+i]*b[i] for i in range(n))%2
def vec(sstr):
    r=[0]*(2*n)
    for i,c in enumerate(sstr):
        if c in 'XY': r[i]=1
        if c in 'ZY': r[n+i]=1
    return r
Sv=[vec(s) for s in S_all]
def rref(rows,nn):
    rows=[r[:] for r in rows]; piv=[]; r=0
    for c in range(2*nn):
        p=next((i for i in range(r,len(rows)) if rows[i][c]),None)
        if p is None: continue
        rows[r],rows[p]=rows[p],rows[r]
        for i in range(len(rows)):
            if i!=r and rows[i][c]: rows[i]=[(x+y)%2 for x,y in zip(rows[i],rows[r])]
        piv.append(c); r+=1
    return rows[:r],piv
Sr,_=rref(Sv,n)
k = n-len(Sr)
print("\nindependent stabiliser generators =",len(Sr),"  => k =",k,"logical qubits")

# centraliser: all v in F_2^{2n} with sym(v,s)=0 for all s in Sr  -- enumerate exhaustively (4^8)
def tostr(vv):
    o=[]
    for i in range(n):
        a,b=vv[i],vv[n+i]
        o.append('I' if (a,b)==(0,0) else 'X' if (a,b)==(1,0) else 'Z' if (a,b)==(0,1) else 'Y')
    return ''.join(o)
cent=[]
for m in range(4**n):
    vv=[0]*(2*n); t=m
    for i in range(n):
        q=t%4; t//=4
        if q in (1,3): vv[i]=1
        if q in (2,3): vv[n+i]=1
    if all(sym(vv,s,n)==0 for s in Sr): cent.append(vv)
print("centraliser size (exhaustive over all %d Paulis) = %d   expected 2^(n+k)=%d"%(4**n,len(cent),2**(n+k)))
# logical = centraliser minus stabiliser (as cosets)
Sspan=set()
for msk in range(2**len(Sr)):
    a=[0]*(2*n)
    for i in range(len(Sr)):
        if (msk>>i)&1: a=[(x+y)%2 for x,y in zip(a,Sr[i])]
    Sspan.add(tuple(a))
logacting=[vv for vv in cent if tuple(vv) not in Sspan]
print("admissible Paulis acting non-trivially on the code space:",len(logacting))
wt=lambda vv: sum(1 for i in range(n) if vv[i] or vv[n+i])
print("MINIMUM WEIGHT of a logical (=code distance d) =",min(wt(vv) for vv in logacting),
      "   lane claimed d=2=L")

# ---------------------------------------------------------------- records: pick a symplectic
# basis of logicals by SEARCH (never nominate): find 2 commuting Z-type logicals with partners
reps=[vv for vv in logacting]
def commutes(a,b): return sym(a,b,n)==0
Zc=None
import random
# search for a pair Z1,Z2 commuting, and X1,X2 with sp(Xi,Zj)=delta
for a in reps:
    for b in reps:
        if not commutes(a,b): continue
        # need partners
        pa=[c for c in reps if sym(c,a,n)==1 and sym(c,b,n)==0]
        pb=[c for c in reps if sym(c,b,n)==1 and sym(c,a,n)==0]
        if pa and pb:
            # ensure a,b independent mod stabiliser
            ab=[(x+y)%2 for x,y in zip(a,b)]
            if tuple(ab) in Sspan: continue
            Zc=(a,b,pa[0],pb[0]); break
    if Zc: break
Z1,Z2,X1,X2 = Zc
print("\nRECORDS FOUND BY SEARCH (not nominated):")
print("   R_1 =",tostr(Z1),"  R_2 =",tostr(Z2))
print("   partners X_1 =",tostr(X1),"  X_2 =",tostr(X2),
      "  sp(X1,Z1)=%d sp(X1,Z2)=%d sp(X2,Z2)=%d sp(X2,Z1)=%d"%(sym(X1,Z1,n),sym(X1,Z2,n),sym(X2,Z2,n),sym(X2,Z1,n)))

def mat(vv):
    M=pau(n,tostr(vv))
    # fix global phase so it is Hermitian with square I
    if np.linalg.norm(M-M.conj().T)>1e-9: M=1j*M
    return M
R1,R2=mat(Z1),mat(Z2)
fam=[R1,R2]
print("\nCLAUSE AUDIT (D-18), independently recomputed on the BRIEF's H:")
for nm,R in (("R_1",R1),("R_2",R2),("R_1R_2",R1@R2)):
    ci = np.linalg.norm(R-R.conj().T)<1e-9 and np.linalg.norm(R@R-np.eye(2**n))<1e-9
    cii= np.linalg.norm(R@H-H@R)<1e-8
    ciii=False; mx=0.0
    for e,C in ES:
        P_=C@C.conj().T; m=C.shape[1]
        M=P_@R@P_
        if np.linalg.norm(M-(np.trace(M)/m)*P_)>1e-8: ciii=True
        mx=max(mx,abs(np.trace(P_@R)))
    print("   %-8s (i)=%s (ii)=%s (iii)=%s  max|Tr(P_E R)|=%.6f  (iv)=%s"%(nm,ci,cii,ciii,mx,mx<1e-7))

# ---------------------------------------------------------------- joint dimension table
def joint(H,fam):
    out={}
    for ei,(e,C) in enumerate(eigspaces(H)):
        groups={():C}
        for R in fam:
            ng={}
            for lab,Cc in groups.items():
                Rs=Cc.conj().T@R@Cc; ws,Vs=np.linalg.eigh(Rs)
                for s in (1,-1):
                    idx=[i for i in range(len(ws)) if (ws[i]>0)==(s>0)]
                    if idx: ng[lab+(s,)]=Cc@Vs[:,idx]
            groups=ng
        for lab,Cc in groups.items(): out[(ei,lab)]=Cc
    return out
B=joint(H,fam)
cfg=list(itertools.product((1,-1),repeat=2))
nE=len(ES)
d={(ei,s):0 for ei in range(nE) for s in cfg}
for (ei,lab),C in B.items(): d[(ei,lab)]=C.shape[1]
print("\nJOINT DIMENSION TABLE on the BRIEF's H:")
print("   shell    E        "+"".join("%9s"%str(s) for s in cfg)+"   uniform?")
for ei,(e,C) in enumerate(ES):
    row=[d[(ei,s)] for s in cfg]
    print("   %5d %7.3f     "%(ei,e)+"".join("%9d"%x for x in row)+"   %s"%(len(set(row))==1))

# ---------------------------------------------------------------- G_W by EXPLICIT unitaries
G=[]
for m_ in range(4):
    eps=tuple((m_>>i)&1 for i in range(2))
    ok=all(d[(ei,s)]==d[(ei,tuple(-x if e_ else x for e_,x in zip(eps,s)))] for ei in range(nE) for s in cfg)
    if not ok: continue
    # build U and VERIFY it
    U=np.zeros((2**n,2**n),dtype=complex); seen=set()
    for ei in range(nE):
        for s in cfg:
            if (ei,s) in seen: continue
            t=tuple(-x if e_ else x for e_,x in zip(eps,s))
            Cs=B.get((ei,s)); Ct=B.get((ei,t))
            if s==t: U+=Cs@Cs.conj().T
            else: U+=Cs@Ct.conj().T+Ct@Cs.conj().T
            seen.add((ei,s)); seen.add((ei,t))
    e1=np.linalg.norm(U.conj().T@U-np.eye(2**n)); e2=np.linalg.norm(U@H-H@U)
    e3=max(np.linalg.norm(U.conj().T@fam[i]@U-((-1)**eps[i])*fam[i]) for i in range(2))
    print("   eps=%s  |U*U-I|=%.2e  |[U,H]|=%.2e  |U*R_iU-eps_iR_i|=%.2e"%(eps,e1,e2,e3))
    assert e1<1e-8 and e2<1e-8 and e3<1e-8
    G.append(eps)
print("   |G_W| =",len(G),"  G_W =",G)

# ---------------------------------------------------------------- orbits: BURNSIDE vs EXPLICIT BASIS
fix=0
for g in G:
    fix+=sum(1 for s in cfg if tuple(-x if e_ else x for e_,x in zip(g,s))==s)
burn=fix/len(G)
# explicit basis: nullspace of (average - I) on R^{cfg}
Mavg=np.zeros((4,4))
idx={s:i for i,s in enumerate(cfg)}
for g in G:
    for s in cfg: Mavg[idx[tuple(-x if e_ else x for e_,x in zip(g,s))],idx[s]]+=1.0/len(G)
ev=np.linalg.eigvals(Mavg-np.eye(4))
dim_inv=int(sum(1 for z in ev if abs(z)<1e-9))
print("\n   #orbits by BURNSIDE = %g ;  dim(invariants) by EXPLICIT NULLSPACE = %d ;  AGREE = %s"
      %(burn,dim_inv,abs(burn-dim_inv)<1e-9))
