# W-25d — DOES THE BOUNDARY EFFECT GROW WITH THE BOUNDARY?
# W-25c: the rim carries a real chirality at Z_3 (exactly 0 at Z_2, structurally) and the
# ENHANCEMENT |<W_rim>| / |prod_k <T_k>| is ~2.5 -- the triangles are positively correlated, so the
# boundary carries more than its pieces independently. Does that survive, or grow, with rim size?
# ONE VARIABLE: the number of rim links n. Same group, same coupling, same theta, same estimator.
import numpy as np, itertools
def wheel(n):
    """hub 0, rim 1..n. spokes 0..n-1 = (0,k+1); rim links n..2n-1 = (k+1, k+2 mod)"""
    E=[(0,k+1) for k in range(n)]+[(k+1,(k+1)%n+1) for k in range(n)]
    return E,n+1
def physical(E,V,N,L):
    keep=[]
    for s in itertools.product(range(N),repeat=L):
        ok=True
        for v in range(V):
            t=sum(s[i] for i,(a,b) in enumerate(E) if a==v)-sum(s[i] for i,(a,b) in enumerate(E) if b==v)
            if t%N: ok=False; break
        if ok: keep.append(s)
    return keep
def loop_op(st,idx,moves,N,L):
    D=len(st); M=np.zeros((D,D),dtype=complex)
    for j,s in enumerate(st):
        t=list(s)
        for l,sg in moves: t[l]=(t[l]+sg)%N
        k=idx.get(tuple(t))
        if k is None: return None
        M[k,j]=1.0
    return M
print(f"  {'N':>3}{'rim n':>7}{'links':>7}{'phys dim':>10}{'Im<W_rim>':>14}{'|<W_rim>|':>12}"
      f"{'|prod<T_k>|':>14}{'enhancement':>13}")
for N in (2,3):
    for n in (3,4,5,6):
        E,V=wheel(n); L=len(E)
        if N**L > 3_000_000: print(f"  {N:>3}{n:>7}{L:>7}   (skipped: {N**L} basis states)"); continue
        st=physical(E,V,N,L); idx={s:j for j,s in enumerate(st)}; D=len(st)
        TRI=lambda k:[(k,+1),(n+k,+1),((k+1)%n,-1)]
        T=[loop_op(st,idx,TRI(k),N,L) for k in range(n)]
        if any(x is None for x in T): print(f"  {N:>3}{n:>7}  triangle left the physical sector"); continue
        Wrim=loop_op(st,idx,[(n+k,+1) for k in range(n)],N,L)
        w=np.exp(2j*np.pi/N)
        Ee=np.diag([sum(w**s[l] for l in range(L)) for s in st])
        th=1.0
        Hm=sum(np.exp(1j*th)*x+np.exp(-1j*th)*x.conj().T for x in T)
        H=-(Hm)-1.0*(Ee+Ee.conj().T); H=(H+H.conj().T)/2
        ev,U=np.linalg.eigh(H); psi=U[:,0]
        f=lambda O: complex(psi.conj()@(O@psi))
        tv=[f(x) for x in T]; wr=f(Wrim)
        pr=np.prod(tv)
        enh = abs(wr)/abs(pr) if abs(pr)>1e-18 else float('nan')
        print(f"  {N:>3}{n:>7}{L:>7}{D:>10}{wr.imag:>14.6f}{abs(wr):>12.6f}{abs(pr):>14.3e}{enh:>13.2f}")
print()
print("  theta = 1.0, g^2 = 1.0 throughout. Enhancement = |<W_rim>| / |prod_k <T_k>|:")
print("  1.0 would mean the rim is exactly what independent triangles give.")
print("  Growth with n would mean the boundary becomes MORE than its pieces as it gets longer.")
