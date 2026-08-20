"""ADVERSARIAL VERIFY 3 -- AXIS 6: does a STANDARD fact reproduce every number?
Three candidate ordinary explanations, each tested."""
import numpy as np, itertools

print("="*100)
print("ORDINARY EXPLANATION 1.  On any stabiliser code, LEMMA B / 'the invariant space is")
print("  1-dimensional' reduces to: a NON-IDENTITY Pauli is TRACELESS.  Nothing more.")
print("="*100)
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
def kron(ms):
    M=np.array([[1]],dtype=complex)
    for m in ms: M=np.kron(M,m)
    return M
def pau(n,s):
    d={'X':X,'Z':Z,'I':I2}; return kron([d[c] for c in s])
n=8; L=2
def h(i,j): return 2*((i%L)*L+(j%L))
def v(i,j): return 2*((i%L)*L+(j%L))+1
St=[];Pl=[]
for i in range(L):
    for j in range(L):
        s=['I']*n
        for e in (h(i,j),h(i-1,j),v(i,j),v(i,j-1)): s[e]='X'
        St.append(''.join(s))
        p=['I']*n
        for e in (h(i,j),h(i,j+1),v(i,j),v(i+1,j)): p[e]='Z'
        Pl.append(''.join(p))
H=-sum(pau(n,s) for s in St+Pl)
w,V=np.linalg.eigh(H); shells=[];i=0
while i<len(w):
    j=i
    while j+1<len(w) and abs(w[j+1]-w[i])<1e-8: j+=1
    shells.append((w[i],V[:,i:j+1])); i=j+1
R1=pau(n,'XIXIIIII'); R2=pau(n,'IZIZIIII')
print("   S            R_S is a non-identity Pauli?   max_E |Tr(P_E R_S)|   plain Tr(R_S)")
for nm,R in (("{1}",R1),("{2}",R2),("{1,2}",R1@R2)):
    mx=max(abs(np.trace((C@C.conj().T)@R)) for _,C in shells)
    print("   %-8s %-30s %20.6f %15.6f"%(nm,"YES",mx,abs(np.trace(R))))
print("   READ: every R_S here is a non-identity Pauli, hence traceless, and P_E is a polynomial")
print("         in stabiliser-group elements, so Tr(P_E R_S) is a sum of traces of non-identity")
print("         Paulis = 0.  LEMMA B on a stabiliser code is that one line.  No record theory needed.")

print()
print("="*100)
print("ORDINARY EXPLANATION 2.  'THEOREM 1' on any stabiliser code is the textbook statement that")
print("  the LOGICAL PAULI GROUP acts regularly on the logical Z-basis.  Check: are the writers")
print("  the lane 'searched for' anything other than the logical X operators?")
print("="*100)
def sym(a,b,nn): return sum(a[t]*b[nn+t]+a[nn+t]*b[t] for t in range(nn))%2
def vecf(s):
    r=[0]*(2*n)
    for t,c in enumerate(s):
        if c in 'XY': r[t]=1
        if c in 'ZY': r[n+t]=1
    return r
Sv=[vecf(s) for s in St+Pl]
Z1=vecf('XIXIIIII'); Z2=vecf('IZIZIIII')
found={}
for m_ in range(4**n):
    vv=[0]*(2*n); t=m_
    for t2 in range(n):
        q=t%4; t//=4
        if q in (1,3): vv[t2]=1
        if q in (2,3): vv[n+t2]=1
    if not all(sym(vv,s,n)==0 for s in Sv): continue
    eps=(sym(vv,Z1,n),sym(vv,Z2,n))
    found.setdefault(eps,0); found[eps]+=1
print("   flip pattern eps  ->  #admissible Paulis realising it")
for eps in sorted(found): print("      %s : %d"%(str(eps),found[eps]))
print("   |G_W| from exhaustive Pauli search =",len(found),"  (matches the lane's 4)")
print("   READ: G_W = (Z_2)^2 because the code HAS logical X operators.  That is the definition")
print("         of a [[n,k,d]] stabiliser code.  Nothing is discovered by searching for them.")

print()
print("="*100)
print("ORDINARY EXPLANATION 3.  The 'TRADE LAW' dim inv = 2^(k - rank G_W) is LAGRANGE.")
print("  The 'THEOREM D' bound k_indep <= min_E v_2(dim E) is 'a free (Z_2)^j action has orbits")
print("  of size 2^j'.  Both tested against brute force on random subgroups.")
print("="*100)
rng=np.random.default_rng(1)
bad=0; tested=0
for k in range(1,8):
    cfg=list(itertools.product((0,1),repeat=k))
    for trial in range(300):
        r=rng.integers(0,k+1)
        gens=[tuple(int(x) for x in rng.integers(0,2,k)) for _ in range(r)]
        G=set()
        for msk in range(2**len(gens)):
            a=[0]*k
            for t in range(len(gens)):
                if (msk>>t)&1: a=[x^y for x,y in zip(a,gens[t])]
            G.add(tuple(a))
        # rank
        basis=[]
        for g in G:
            vv=list(g)
            for b in basis:
                hd=next((t for t in range(k) if b[t]),None)
                if hd is not None and vv[hd]: vv=[x^y for x,y in zip(vv,b)]
            if any(vv): basis.append(vv)
        rank=len(basis)
        # orbits by brute force
        seen=set(); orb=0
        for s in cfg:
            if s in seen: continue
            orb+=1
            for g in G: seen.add(tuple(x^y for x,y in zip(g,s)))
        tested+=1
        if orb != 2**(k-rank): bad+=1
print("   random subgroups of (Z_2)^k, k=1..7, %d cases:  #orbits != 2^(k-rank) in %d cases"%(tested,bad))
print("   READ: the 'TRADE LAW verified on every family' is the orbit-counting theorem for a")
print("         subgroup acting on itself by translation -- Lagrange.  It cannot fail.")
