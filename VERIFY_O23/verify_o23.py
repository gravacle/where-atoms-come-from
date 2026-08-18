"""ADVERSARIAL VERIFICATION of LANE_O23_SPLITTING.  Independent re-implementation."""
import itertools, numpy as np
rng = np.random.default_rng(7)
I2=np.eye(2,dtype=complex); SX=np.array([[0,1],[1,0]],dtype=complex)
SY=np.array([[0,-1j],[1j,0]],dtype=complex); SZ=np.array([[1,0],[0,-1]],dtype=complex)
P={'I':I2,'X':SX,'Y':SY,'Z':SZ}
def pm(s):
    M=np.array([[1.0+0j]])
    for c in s: M=np.kron(M,P[c])
    return M
def wt(s): return sum(1 for c in s if c!='I')
def sy(a,b):
    c=0
    for x,y in zip(a,b):
        if x!='I' and y!='I' and x!=y: c^=1
    return c
def pmul(a,b):
    o=[]
    for x,y in zip(a,b):
        if x=='I': o.append(y)
        elif y=='I': o.append(x)
        elif x==y: o.append('I')
        else: o.append({'XY':'Z','YX':'Z','XZ':'Y','ZX':'Y','YZ':'X','ZY':'X'}[x+y])
    return ''.join(o)
def stabgrp(g,n):
    S={'I'*n}; fr=['I'*n]
    while fr:
        nx=[]
        for a in fr:
            for x in g:
                b=pmul(a,x)
                if b not in S: S.add(b); nx.append(b)
        fr=nx
    return S
def codespace(gens,n):
    N=2**n; H=np.zeros((N,N),complex)
    for g in gens: H-=pm(g)
    H=(H+H.conj().T)/2
    ev,U=np.linalg.eigh(H)
    sel=np.abs(ev-ev[0])<1e-8*max(1.0,abs(ev).max())
    Q=U[:,sel]
    return H,Q,Q@Q.conj().T
def phi(V,Q):
    B=Q.conj().T@V@Q; m=B.shape[0]
    return B-np.trace(B)/m*np.eye(m)
def paulis_w(n,w):
    out=[]
    for sup in itertools.combinations(range(n),w):
        for L in itertools.product('XYZ',repeat=w):
            s=['I']*n
            for k,c in zip(sup,L): s[k]=c
            out.append(''.join(s))
    return out

print("="*100); print("PROBE 1.  TORIC 2x2: does the ON-SITE FINITE symmetry X^(x)8 act NON-SCALARLY on E0?")
print("  If yes: a FINITE on-site 0-form symmetry with a NON-SCALAR action on the ground multiplet,")
print("  and NO weight-1 operator splits it -> the finite-G step is not merely open.")
print("="*100)
h=lambda i,j:(j%2)*2+(i%2); v=lambda i,j:4+(j%2)*2+(i%2)
gens=[]
for j in range(2):
    for i in range(2):
        s=['I']*8
        for L in [h(i,j),h(i-1,j),v(i,j),v(i,j-1)]: s[L]='X' if s[L]=='I' else 'I'
        gens.append(''.join(s))
        s=['I']*8
        for L in [h(i,j),v(i+1,j),h(i,j+1),v(i,j)]: s[L]='Z' if s[L]=='I' else 'I'
        gens.append(''.join(s))
def indep(gs,n):
    rows,keep=[],[]
    for g in gs:
        vec=0
        for k,c in enumerate(g):
            if c in 'XY': vec|=1<<(2*k)
            if c in 'ZY': vec|=1<<(2*k+1)
        cur=vec
        for r in rows:
            p=r.bit_length()-1
            if cur>>p&1: cur^=r
        if cur: rows.append(cur); rows.sort(reverse=True); keep.append(g)
    return keep
gens=indep(gens,8); n=8; S=stabgrp(gens,n)
H0,Q,Pj=codespace(gens,n)
print(f"  code space dim = {Q.shape[1]}, independent stabilisers = {len(gens)}")
GX='X'*8
print(f"  G = <X^(x)8>.  commutes with every stabiliser: {all(sy(GX,g)==0 for g in gens)}")
print(f"  X^(x)8 in the stabiliser group (i.e. acts as identity on E0): {GX in S}")
UG=pm(GX)
print(f"  ||[H0, X^(x)8]|| = {np.linalg.norm(H0@UG-UG@H0):.2e}   (a genuine ON-SITE symmetry: product of single-site X)")
B=Q.conj().T@UG@Q
print(f"  action on E0, eigenvalues of P X^(x)8 P = {np.round(np.linalg.eigvalsh((B+B.conj().T)/2),6)}")
print(f"  ||Phi(X^(x)8)|| on E0 = {np.linalg.norm(phi(UG,Q)):.4f}   (0 = SCALAR: the candidate fails)" )
best=max(float(np.linalg.norm(phi(pm(s),Q))) for s in paulis_w(8,1))
print(f"  D1 = max over ALL single-site Hermitian Paulis of ||Phi|| = {best:.3e}  -> NO local order parameter")
print("  RESULT: X^(x)8 IS in the stabiliser group here, so it acts as the IDENTITY on E0.")
print("          This candidate finite-G counterexample FAILS.  (See verify_finiteG.py for one")
print("          that works: the Steane code's TRANSVERSAL logical X^(x)7.)")

print()
print("="*100); print("PROBE 2.  CASE C / CASE D: is clause (iii) satisfied ON THE RECORD'S OWN BLOCK,")
print("  or only via untouched EXCITED eigenspaces?  (the lane's (iii) score is a max over ALL eigenspaces)")
print("="*100)
def blocks_from_E0(He,Q,tol=1e-9):
    """eigen-decompose He restricted to the E0 subspace spanned by Q."""
    Hs=Q.conj().T@He@Q; Hs=(Hs+Hs.conj().T)/2
    ev,U=np.linalg.eigh(Hs)
    out=[];i=0
    while i<len(ev):
        j=i
        while j+1<len(ev) and abs(ev[j+1]-ev[i])<1e-9*max(1.0,abs(ev[i])): j+=1
        out.append((ev[i],Q@U[:,i:j+1])); i=j+1
    return out
logic=[s for s in paulis_w(8,2) if all(sy(s,g)==0 for g in gens) and s not in S]
Z1=X1=Z2=X2=None
for a in logic:
    for b in logic:
        if sy(a,b)!=1: continue
        for c in logic:
            if sy(c,a) or sy(c,b): continue
            for dd in logic:
                if sy(dd,a) or sy(dd,b) or sy(dd,c)!=1: continue
                Z1,X1,Z2,X2=a,b,c,dd; break
            if X2: break
        if X2: break
    if X2: break
print(f"  toric logicals: Z1={Z1} X1={X1} Z2={Z2} X2={X2}")
R=pm(Z2)
for eps,lab in [(1e-3,"CASE C: H = H0 + 1e-3 Zbar1")]:
    He=H0+eps*pm(Z1)
    print(f"  {lab}")
    for e,Qb in blocks_from_E0(He,Q):
        b=Qb.conj().T@R@Qb
        dev=np.linalg.norm(b-np.trace(b)/b.shape[0]*np.eye(b.shape[0]))
        print(f"    E0-descended block at E={e:+.6f}, dim {b.shape[0]}: R deviation from scalar = {dev:.4f}"
              f"  {'(iii) HOLDS on this block' if dev>1e-9 else '(iii) FAILS on this block'}")
gensD=["ZZIIIIII","IZZIIIII","IIIXZZXI","IIIIXZZX","IIIXIXZZ","IIIZXIXZ"]
SD=stabgrp(gensD,8); HD,QD,PD=codespace(gensD,8)
ZA="ZIIIIIII"; ZB="IIIZZZZZ"; XB="IIIXXXXX"
print(f"  composite code space dim = {QD.shape[1]}")
HeD=HD+1e-3*pm(ZA)
for e,Qb in blocks_from_E0(HeD,QD):
    b=Qb.conj().T@pm(ZB)@Qb
    dev=np.linalg.norm(b-np.trace(b)/b.shape[0]*np.eye(b.shape[0]))
    print(f"    CASE D block at E={e:+.6f}, dim {b.shape[0]}: R=Zbar_B deviation from scalar = {dev:.4f}"
          f"  {'(iii) HOLDS' if dev>1e-9 else '(iii) FAILS'}")

print()
print("="*100); print("PROBE 3.  the 1.4e-30 entry, and d_R by brute force (independent).")
print("="*100)
def dR(gens,n,S,Rstr,wmax):
    for w in range(1,wmax+1):
        for s in paulis_w(n,w):
            if all(sy(s,g)==0 for g in gens) and s not in S and sy(s,Rstr)==1: return w,s
    return None,None
for lab,gg,SS,QQ,Rs in [("toric R=Zbar2",gens,S,Q,Z2),("composite R=Zbar_B",gensD,SD,QD,ZB),
                        ("composite R=Zbar_A",gensD,SD,QD,ZA)]:
    w,s=dR(gg,8,SS,Rs,4)
    R0=QQ.conj().T@pm(Rs)@QQ
    row=[]
    for ww in range(1,4):
        worst=max(float(np.linalg.norm(phi(pm(p),QQ)@R0-R0@phi(pm(p),QQ))) for p in paulis_w(8,ww))
        row.append(f"w={ww}: {worst:.2e}")
    print(f"  {lab:<22s} d_R = {w} (witness {s})   max||[Phi,R]||  " + "  ".join(row))

print()
print("="*100); print("PROBE 4.  SLOPE LAW ROBUSTNESS: 5 fresh random weight-1 perturbations per code,")
print("  and a SECOND, DISJOINT eps window.  Does slope = d survive, or is it window-tuned?")
print("="*100)
CODES={"[[3,1,1]] d=1":(["ZZI","IZZ"],1),"[[4,2,2]] d=2":(["XXXX","ZZZZ"],2),
       "toric2x2 d=2":(gens,2),"[[5,1,3]] d=3":(["XZZXI","IXZZX","XIXZZ","ZXIXZ"],3),
       "[[7,1,3]] d=3":(["IIIXXXX","IXXIIXX","XIXIXIX","IIIZZZZ","IZZIIZZ","ZIZIZIZ"],3)}
def splitting(H,P0,m):
    ev,U=np.linalg.eigh(H)
    ov=np.einsum('ij,jk,ki->i',U.conj().T,P0,U).real
    idx=np.argsort(-ov)[:m]
    e=np.sort(ev[idx]); return e[-1]-e[0]
for name,(gg,d) in CODES.items():
    nn=len(gg[0]); Hc,Qc,Pc=codespace(gg,nn); mm=Qc.shape[1]
    for wintag,eps_l in [("windowA",np.array([1e-4,3e-4,1e-3,3e-3])),
                         ("windowB",np.array([1e-2,2e-2,4e-2,8e-2]))]:
        sls=[]
        for t in range(5):
            V=sum(rng.normal()*pm(s) for s in paulis_w(nn,1)); V=(V+V.conj().T)/2
            V=V/np.linalg.norm(V)*np.linalg.norm(Hc)
            xs,ys=[],[]
            for e in eps_l:
                sp=splitting(Hc+e*V,Pc,mm)
                if sp>1e-12: xs.append(np.log(e)); ys.append(np.log(sp))
            sls.append(np.polyfit(xs,ys,1)[0] if len(xs)>1 else float('nan'))
        print(f"  {name:<16s} d={d}  {wintag}: slopes over 5 fresh draws = {[round(float(x),3) for x in sls]}")

print()
print("="*100); print("PROBE 5.  can the O-2 'rank Phi = 0 below d' claim be broken by a NON-PAULI local")
print("  perturbation?  sweep random single-site HERMITIAN operators (not just X,Y,Z) on [[5,1,3]].")
print("="*100)
gg=["XZZXI","IXZZX","XIXZZ","ZXIXZ"]; Hc,Qc,Pc=codespace(gg,5)
worst=0.0
for _ in range(500):
    k=rng.integers(0,5); A=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)); A=(A+A.conj().T)/2
    ms=[I2]*5; ms[k]=A
    M=np.array([[1.0+0j]])
    for x in ms: M=np.kron(M,x)
    worst=max(worst,float(np.linalg.norm(phi(M,Qc))/max(np.linalg.norm(A),1e-300)))
print(f"  500 random single-site Hermitian perturbations on [[5,1,3]]: max ||Phi||/||A|| = {worst:.3e}")
w2=0.0
for _ in range(500):
    ks=rng.choice(5,size=3,replace=False)
    ms=[I2]*5; nrm=1.0
    for k in ks:
        A=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)); A=(A+A.conj().T)/2
        ms[k]=A; nrm*=np.linalg.norm(A)
    M=np.array([[1.0+0j]])
    for x in ms: M=np.kron(M,x)
    w2=max(w2,float(np.linalg.norm(phi(M,Qc))/max(nrm,1e-300)))
print(f"  POSITIVE CONTROL, 500 random WEIGHT-3 product Hermitians:  max ||Phi||/||A|| = {w2:.3e}")
