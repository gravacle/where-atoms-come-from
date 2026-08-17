# ARE THE PAIRS OURS OR THE OBJECT'S?  Z_2 forces every charge to +-1, every sector to two, every
# operator to be its own inverse. If the structure this program has been finding is genuinely
# two-valued rather than an artefact of the smallest group, it must survive at Z_3.
# ONE VARIABLE: N. Same graph, same cut, same region, same code path.
import numpy as np, itertools
def clock(N): return np.diag([np.exp(2j*np.pi*k/N) for k in range(N)])
def shift(N):
    M=np.zeros((N,N),dtype=complex)
    for k in range(N): M[(k+1)%N,k]=1
    return M
E=[(0,1),(1,2),(2,0),(0,3)]        # triangle in the bulk + one dangling link to boundary vertex 3
BULK=[0,1,2]; L=len(E)

def build(N):
    d=N**L
    def emb(i,A):
        M=np.array([[1]],dtype=complex)
        for j in range(L): M=np.kron(M,A if j==i else np.eye(N))
        return M
    Z=[emb(i,clock(N)) for i in range(L)]; X=[emb(i,shift(N)) for i in range(L)]
    def G(v):
        M=np.eye(d,dtype=complex)
        for i,(a,b) in enumerate(E):
            if a==v: M=M@Z[i]
            if b==v: M=M@Z[i].conj().T
        return M
    P=np.eye(d,dtype=complex)
    for v in BULK:
        Pv=np.zeros((d,d),dtype=complex); Gv=G(v); acc=np.eye(d,dtype=complex)
        for k in range(N): Pv+=acc; acc=acc@Gv
        P=P@(Pv/N)
    return d,Z,X,G,P

def S(r,base):
    w=np.linalg.eigvalsh(r); w=w[w>1e-12]
    return float(-np.sum(w*np.log(w))/np.log(base))
def red(psi,keep,N):
    v=psi.reshape([N]*L)
    tr=[i for i in range(L) if i not in keep]
    r=np.tensordot(v,v.conj(),axes=(tr,tr))
    m=N**len(keep); return r.reshape(m,m)

rng=np.random.default_rng(20260831)
print(f"  {'N':>3}{'full dim':>10}{'phys dim':>10}{'sectors at bdry':>17}{'Gauss identity':>17}{'I(S:F)/H(S)':>26}")
for N in (2,3):
    d,Z,X,G,P = build(N)
    w_,vec=np.linalg.eigh(P); B=vec[:,np.abs(w_)>0.5]; dp=B.shape[1]
    Q=G(3)
    ev=np.round(np.linalg.eigvals(B.conj().T@Q@B),6)
    nsec=len({complex(np.round(x.real,4),np.round(x.imag,4)) for x in ev})
    # Gauss identity at vertex 0: link 0 outgoing, link 2 incoming, link 3 outgoing
    lhs=Z[0]@Z[2].conj().T@Z[3]
    gid=np.linalg.norm((lhs-np.eye(d))@P)
    # I(S:F)/H(S) with S = link 0, F = the rest of the cut at vertex 0
    cutF=[2,3]; vals=[]
    for _ in range(3):
        c=rng.normal(size=dp)+1j*rng.normal(size=dp); psi=B@c; psi/=np.linalg.norm(psi)
        hs=S(red(psi,[0],N),N); hf=S(red(psi,cutF,N),N); hj=S(red(psi,[0]+cutF,N),N)
        vals.append((hs+hf-hj)/hs if hs>1e-12 else float('nan'))
    print(f"  {N:>3}{d:>10}{dp:>10}{nsec:>17}{gid:>17.3e}{'  '.join(f'{v:.6f}' for v in vals):>26}")

print("\n== THE OPERATOR CLASSIFICATION (W-22's three-way split) AT EACH N ==")
for N in (2,3):
    d,Z,X,G,P = build(N)
    Gb=[G(v) for v in BULK]; Q=G(3)
    cats={}
    for i in range(L):
        for nm,A in (("Z",Z[i]),("X",X[i])):
            inv = max(np.linalg.norm(A@g-g@A) for g in Gb) < 1e-9
            # does it change the boundary sector? Q A Q^-1 != A
            moves = np.linalg.norm(Q@A@Q.conj().T - A) > 1e-9
            cat = "invariant" if inv and not moves else ("sector-changing" if moves and not inv else
                  ("invariant+moves" if inv and moves else "junk"))
            cats[cat]=cats.get(cat,0)+1
    print(f"  N={N}: {cats}")
print("\n  If the counts and the ratios keep the SAME SHAPE at N=3, the pairs are OURS.")
print("  If the structure changes shape, everything measured so far is a Z_2 statement.")
