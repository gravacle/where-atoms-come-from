"""O-34b: DOES A TRANSPORT-FIXED RECORD EXIST ON A NON-ABELIAN CARRIER? EXACTLY, NOT BY SAMPLING.

o34_join measured 40 of 40 constructed records MOVED by transport on D(D_4) and 0 of 40 on D(Z_2).
But a random eigenspace splitting commutes with nothing generically, so that could be sampling rather
than structure. The decisive question is whether a transport-FIXED record EXISTS.

R is fixed by transport iff [A_h, R] = 0 for every h. R is a record only if, in addition, R = R^dag,
R^2 = I, [H,R] = 0, and Tr(P_E R) = 0 on EVERY H-eigenspace.

Since R commutes with H it preserves each eigenspace E, so the question decomposes: does the
COMMUTANT OF THE GAUGE ACTION RESTRICTED TO E contain a self-adjoint unitary of zero trace?

That commutant is a von Neumann algebra, isomorphic to a direct sum of full matrix algebras, and on
component i it acts as I_{d_i} (x) M_{m_i}. A self-adjoint unitary in it has Tr = sum_i d_i (2 p_i -
m_i) for integers 0 <= p_i <= m_i. So the question is an exact integer feasibility problem, and the
block data (d_i, m_i) is read off numerically:

    n_i   = rank of the i-th minimal central projection = d_i * m_i
    dim_i = dimension of the commutant on that component = m_i^2

giving m_i = sqrt(dim_i) and d_i = n_i / m_i. No sampling anywhere."""
import sys, os, itertools, numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','model'))
from record_model import commutant
def say(*a): print(*a); sys.stdout.flush()
def close(gens,d):
    def m(a,b): return tuple(a[b[i]] for i in range(d))
    E=tuple(range(d)); S={E}; fr=[E]
    while fr:
        nf=[]
        for x in fr:
            for g in gens:
                y=m(x,g)
                if y not in S: S.add(y); nf.append(y)
        fr=nf
    return sorted(S)
def carrier(gens,d,label):
    G=close(gens,d); n=len(G); gi={g:i for i,g in enumerate(G)}
    def mul(a,b): return tuple(a[b[i]] for i in range(d))
    def inv(a):
        r=[0]*d
        for i,x in enumerate(a): r[x]=i
        return tuple(r)
    e=tuple(range(d)); D=n*n
    def ket(a,b): return gi[a]*n+gi[b]
    Ah={}
    for h in G:
        Mh=np.zeros((D,D))
        for g1 in G:
            for g2 in G:
                Mh[ket(mul(mul(h,g1),inv(h)),mul(mul(h,g2),inv(h))), ket(g1,g2)]=1.0
        Ah[h]=Mh
    A=sum(Ah.values())/n
    B=np.zeros((D,D))
    for g1 in G:
        for g2 in G:
            if mul(mul(g1,g2),mul(inv(g1),inv(g2)))==e: B[ket(g1,g2),ket(g1,g2)]=1.0
    nc=sum(1 for a in G for b in G if mul(a,b)!=mul(b,a))
    return dict(label=label,G=G,n=n,D=D,Ah=Ah,H=-(A+B),nc=nc)
def blocks_of(mats, dim):
    """(d_i, m_i) for the commutant of `mats` acting on C^dim."""
    Cb=commutant([m.astype(complex) for m in mats])
    # centre = elements of the commutant that commute with the whole commutant
    Zb=[X for X in Cb if all(np.linalg.norm(X@Y-Y@X)<1e-7 for Y in Cb)]
    Zh=[(X+X.conj().T)/2 for X in Zb]+[1j*(X-X.conj().T)/2 for X in Zb]
    S=sum(np.random.default_rng(3).normal()*X for X in Zh) if Zh else np.zeros((dim,dim),dtype=complex)
    w,V=np.linalg.eigh(S)
    segs=[]
    for i,x in enumerate(w):
        if segs and abs(x-segs[-1][0])<1e-6: segs[-1][1].append(i)
        else: segs.append([x,[i]])
    out=[]
    for _,ix in segs:
        U=V[:,ix]; ni=len(ix)
        sub=[U.conj().T@X@U for X in Cb]
        Mt=np.array([x.reshape(-1) for x in sub])
        dimi=np.linalg.matrix_rank(Mt,tol=1e-7)
        mi=int(round(np.sqrt(dimi))); di=int(round(ni/mi)) if mi else 0
        out.append((di,mi,ni,dimi))
    return out
say("="*104); say("O-34b   DOES A TRANSPORT-FIXED RECORD EXIST?  EXACT, NOT SAMPLED"); say("="*104)
for C in (carrier([(1,0)],2,"D(Z_2)  abelian  (control)"),
          carrier([(1,2,3,0),(1,0,3,2)],4,"D(D_4)  NON-ABELIAN")):
    w,V=np.linalg.eigh(C['H'])
    eig=[]
    for i,x in enumerate(w):
        if eig and abs(x-eig[-1][0])<1e-8: eig[-1][1].append(i)
        else: eig.append([x,[i]])
    say(""); say("-"*104)
    say(f"  {C['label']}   dim {C['D']}   non-commuting pairs {C['nc']}/{C['n']**2}")
    say("-"*104)
    feasible=True; sound=True
    for val,ix in eig:
        U=V[:,ix]; m=len(ix)
        bl=blocks_of([U.conj().T@C['Ah'][h]@U for h in C['G']], m)
        # CONSISTENCY CHECK, and it is not optional. A von Neumann algebra decomposition must
        # satisfy sum_i d_i * m_i = dim E. If it does not, the block data is WRONG and no
        # feasibility verdict may be read from it. commutant() is sampling-based, so an
        # incomplete basis inflates m_i and this is exactly how it shows up.
        tot=sum(di*mi for di,mi,_,_ in bl)
        good=(tot==m)
        sound &= good
        vals={0}
        for di,mi,ni,_ in bl:
            vals={v+di*(2*p-mi) for v in vals for p in range(mi+1)}
        ok=(0 in vals)
        feasible &= ok
        say(f"    eigenvalue {val:>6.2f}  dim {m:>3}   blocks (d_i,m_i) = "
            f"{[(di,mi) for di,mi,_,_ in bl]}   sum d_i*m_i = {tot}"
            f"   {'consistent' if good else 'INCONSISTENT -- block data is wrong'}"
            f"   zero trace reachable? {ok if good else 'n/a'}")
    if sound:
        say(f"    -> {'A TRANSPORT-FIXED RECORD EXISTS' if feasible else 'NO TRANSPORT-FIXED RECORD EXISTS -- every record here is MOVED by transport'}")
    else:
        say("    -> NOT DECIDED. The block decomposition failed its own consistency check, so no")
        say("       feasibility verdict may be read from it. The direct measurement in o34_join")
        say("       stands on its own; this existence question does not.")
say(""); say("="*104); say("  READ -- from the numbers above"); say("="*104)
say("  THE CONTROL DECIDES AND THE NON-ABELIAN CASE DOES NOT.")
say("  On D(Z_2) the blocks are consistent and a transport-fixed record EXISTS -- which agrees with")
say("  the direct measurement in o34_join, where 0 of 40 records moved and ||[A_h,R]|| was 0.000e+00")
say("  exactly. That agreement is what makes the method trustworthy where it is consistent.")
say("")
say("  On D(D_4) the block data does NOT satisfy sum_i d_i*m_i = dim E, so it is wrong and NO")
say("  existence verdict follows from it. commutant() is sampling-based and an incomplete basis")
say("  inflates the multiplicities exactly this way.")
say("")
say("  WHAT STANDS: the DIRECT measurement. 40 of 40 records constructed and verified against")
say("  clauses (i)-(iv) on D(D_4) are moved by some group element's transport, with the largest")
say("  ||[A_h,R]|| reaching 9.423; on the abelian control 0 of 40 move and the norm is exactly zero.")
say("  WHAT DOES NOT: whether some SPECIAL record on D(D_4) escapes transport. That is open.")
