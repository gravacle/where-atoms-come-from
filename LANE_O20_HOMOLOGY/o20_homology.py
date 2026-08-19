"""O-20: is chi a function of the coupling's HOMOLOGY CLASS alone?

O-18's dilemma was mis-stated. The coupling that failed (plaquettes) is a BOUNDARY, in B_1.
The coupling that succeeded (Zbar) is NOT a boundary, in Z_1 \\ B_1. Both are gauge-invariant.
They differ in homology class -- the quotient Z_1/B_1 the program calls Gamma.

PREDICTION: chi = 0 for every coupling in B_1; chi > 0 for every coupling in Z_1 \\ B_1;
and chi does not depend on weight or locality within a class."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py').read().split('say("="*104); say("0.')[0])
N=2**L; AV=[op({l:X for l in s},L) for s in STAR]
E0,V0=np.linalg.eigh(H0); gs=int(np.sum(np.abs(E0-E0[0])<1e-9)); Pg=V0[:,:gs]@V0[:,:gs].conj().T
# --- GF(2) cycle and boundary spaces on the link set ---
d1=np.zeros((nx*ny,L),dtype=np.int8)
for k_,(a,b) in enumerate(EDGES): d1[a,k_]^=1; d1[b,k_]^=1
d2=np.zeros((L,len(PLAQ)),dtype=np.int8)
for k_,p in enumerate(PLAQ):
    for e in p: d2[e,k_]^=1
def toint(v): return int(''.join(map(str,v)),2)
def nullspace(M):
    M=M.copy()%2; rows,cols=M.shape; pc=[]; r=0
    for c in range(cols):
        p=next((i for i in range(r,rows) if M[i,c]),None)
        if p is None: continue
        M[[r,p]]=M[[p,r]]
        for i in range(rows):
            if i!=r and M[i,c]: M[i]^=M[r]
        pc.append(c); r+=1
    out=[]
    for fc in [c for c in range(cols) if c not in pc]:
        v=np.zeros(cols,dtype=np.int8); v[fc]=1
        for i,p_ in enumerate(pc): v[p_]=M[i,fc]
        out.append(v)
    return out
Z1=[int(toint(v)) for v in nullspace(d1)]
B1raw=[int(toint(d2[:,c])) for c in range(d2.shape[1])]
piv={}
def insert(v):
    while v:
        h=v.bit_length()-1
        if h in piv: v^=piv[h]
        else: piv[h]=v; return True
    return False
for b in B1raw: insert(b)
def in_B1(v):
    for h in sorted(piv,reverse=True):
        if (v>>h)&1: v^=piv[h]
    return v==0
def zop(c): return op({l:Z for l in range(L) if (c>>(L-1-l))&1},L)
# --- bath and chi ---
nq=3; nB=2**nq; beta=2.0
I2b=np.eye(2); Xb2=np.array([[0,1],[1,0]],dtype=complex); Zb2=np.array([[1,0],[0,-1]],dtype=complex)
def bop(j,P):
    M=np.array([[1]],dtype=complex)
    for k in range(nq): M=np.kron(M,P if k==j else I2b)
    return M
HB=sum(w*bop(j,Zb2) for j,w in enumerate([1.0,1.4,0.7])); CB=sum(bop(j,Xb2) for j in range(nq))
def vN(r):
    e=np.linalg.eigvalsh(r); e=e[e>1e-13]; return float(-(e*np.log2(e)).sum())
def chi_of(A,lam=0.8,t=4.0):
    Ht=np.kron(H0,np.eye(nB))+np.kron(np.eye(N),HB)+lam*np.kron(A,CB)
    w,U=np.linalg.eigh(Ht); wB,VB=np.linalg.eigh(HB); pB=np.exp(-beta*wB); pB/=pB.sum()
    r0=np.kron(Pg/gs,(VB*pB)@VB.conj().T); Uc=U.conj().T@r0@U; ph=np.exp(-1j*w*t)
    r=U@(ph[:,None]*Uc*ph.conj()[None,:])@U.conj().T
    out=[]
    for s in (+1,-1):
        P=np.kron((np.eye(N)+s*Zbar)/2,np.eye(nB)); blk=P@r@P; p=np.real(np.trace(blk))
        if p<1e-12: continue
        out.append((p,(blk/p).reshape(N,nB,N,nB).trace(axis1=0,axis2=2)))
    if len(out)<2: return 0.0
    av=sum(p*rb for p,rb in out)
    return max(vN(av)-sum(p*vN(rb) for p,rb in out),0.0)

say("="*104); say("O-20  IS chi A FUNCTION OF THE COUPLING'S HOMOLOGY CLASS ALONE?"); say("="*104)
say(f"  toric 2x2, dim {N}; dim Z_1 = {len(Z1)}, dim B_1 = {len(piv)}, dim H_1 = {len(Z1)-len(piv)}")
say(f"  every coupling below is a Z-type CYCLE, hence GAUGE-INVARIANT by construction")
say("")
say(f"  {'coupling (link set)':<26}{'weight':>7}{'gauge-inv':>11}{'in B_1?':>9}{'homology class':>16}{'chi bits':>13}")
seen={}
allcyc=set()
for m in range(1,1<<len(Z1)):
    v=0
    for i in range(len(Z1)):
        if (m>>i)&1: v^=Z1[i]
    if v: allcyc.add(v)
for v in sorted(allcyc, key=lambda x:(bin(x).count('1'),x)):
    A=zop(v); gi=max(np.linalg.norm(A@G-G@A) for G in AV)
    cls='B_1 (boundary)' if in_B1(v) else 'NOT a boundary'
    c=chi_of(A)
    seen.setdefault(cls,[]).append(c)
    if len(seen[cls])<=4 or bin(v).count('1')<=2:
        links=[l for l in range(L) if (v>>(L-1-l))&1]
        say(f"  {str(links):<26}{len(links):>7}{gi:>11.1e}{str(in_B1(v)):>9}{cls:>16}{c:>13.8f}")
say("")
for cls,vals in seen.items():
    a=np.array(vals)
    say(f"  {cls:<18} n = {len(a):>3}   chi range [{a.min():.8f}, {a.max():.8f}]   spread {a.max()-a.min():.2e}")
say("")
b=np.array(seen.get('B_1 (boundary)',[0])); nb=np.array(seen.get('NOT a boundary',[0]))
say(f"  PREDICTION: chi = 0 on B_1, chi > 0 off it.")
say(f"    all boundaries give exactly 0 : {bool(b.max()<1e-10)}")
say(f"    all non-boundaries give > 0   : {bool(nb.min()>1e-10)}")
say(f"    -> {'CONFIRMED: chi is a function of the HOMOLOGY CLASS' if b.max()<1e-10 and nb.min()>1e-10 else 'NOT CONFIRMED'}")
