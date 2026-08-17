# W-25c — THE TRIANGLE MEASUREMENT. Is the small rim chirality a SCREENING effect, or did the
# registrar compare incommensurable objects?
# CLAIM UNDER TEST (the registrar's, from the previous turn): "Im<W_rim> is 20-30x smaller than
# Im<W_tri> despite Stokes saying the rim IS their product, so the four must nearly cancel."
# THE OBVIOUS HOLE, CHECKED HERE: Stokes is an OPERATOR identity, W_rim = T0 T1 T2 T3. It says
# NOTHING about <W_rim> versus <T_k>. Comparing a four-operator product to a single operator may be
# comparing <x^4> to <x>.
import numpy as np, itertools
E=[(0,1),(0,2),(0,3),(0,4),(1,2),(2,3),(3,4),(4,1)]
SPOKE=[0,1,2,3]; RIM=[4,5,6,7]; L=len(E); V=5
def physical(N):
    keep=[]
    for s in itertools.product(range(N),repeat=L):
        ok=True
        for v in range(V):
            t=sum(s[i] for i,(a,b) in enumerate(E) if a==v)-sum(s[i] for i,(a,b) in enumerate(E) if b==v)
            if t%N: ok=False; break
        if ok: keep.append(s)
    return keep
def loop_op(st,idx,moves,N):
    D=len(st); M=np.zeros((D,D),dtype=complex)
    for j,s in enumerate(st):
        t=list(s)
        for l,sg in moves: t[l]=(t[l]+sg)%N
        k=idx.get(tuple(t))
        if k is None: return None
        M[k,j]=1.0
    return M
def elec(st,links,N):
    w=np.exp(2j*np.pi/N)
    return np.diag([sum(w**s[l] for l in links) for s in st])
TRI=lambda k:[(k,+1),(4+k,+1),((k+1)%4,-1)]

for N in (2,3):
    st=physical(N); idx={s:j for j,s in enumerate(st)}; D=len(st)
    T=[loop_op(st,idx,TRI(k),N) for k in range(4)]
    Wrim=loop_op(st,idx,[(4,+1),(5,+1),(6,+1),(7,+1)],N)
    prod=T[0]
    for x in T[1:]: prod=prod@x
    Ee=elec(st,list(range(L)),N)
    print(f"\n  ===== N={N}, physical dim {D} =====")
    print(f"  operator identity  || T0 T1 T2 T3 - W_rim || = {np.linalg.norm(prod-Wrim):.2e}   (Stokes)")
    print(f"  {'theta':>7} | {'<T_0>':>20}{'<T_1>':>20}{'<T_2>':>20}{'<T_3>':>20}")
    print(f"  {'':>7} | {'<W_rim>':>20}{'prod of <T_k>':>20}{'|<W_rim>|':>20}{'|prod <T_k>|':>20}")
    for th in (0.4,1.0,2.0):
        Hm=sum(np.exp(1j*th)*x+np.exp(-1j*th)*x.conj().T for x in T)
        H=-(Hm)-1.0*(Ee+Ee.conj().T)
        H=(H+H.conj().T)/2
        ev,U=np.linalg.eigh(H); psi=U[:,0]
        f=lambda O: complex(psi.conj()@(O@psi))
        tv=[f(x) for x in T]; wr=f(Wrim)
        pr=tv[0]*tv[1]*tv[2]*tv[3]
        fmt=lambda z: f"{z.real:+.5f}{z.imag:+.5f}j"
        print(f"  {th:>7.2f} | " + "".join(f"{fmt(z):>20}" for z in tv))
        print(f"  {'':>7} | {fmt(wr):>20}{fmt(pr):>20}{abs(wr):>20.6f}{abs(pr):>20.6f}")
print()
print("  IF the four <T_k> are EQUAL (the wheel has 4-fold symmetry) they CANNOT cancel by sign.")
print("  IF |<W_rim>| exceeds |prod of <T_k>| then the rim is not suppressed at all -- a product of")
print("  four operators simply has a smaller expectation than one, and the registrar's comparison")
print("  was between incommensurable objects.")
