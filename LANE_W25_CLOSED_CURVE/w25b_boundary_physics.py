# W-25b — THE QUESTIONS THAT WERE UNASKABLE, ASKED ON A REAL BOUNDARY.
# Carrier: wheel W_4. Hub 0 = interior, rim 1-4 = a CLOSED CURVE = the boundary.
# All four rim links lie on cycles, the rim loop is nontrivial, and rim = product of the four
# interior triangles exactly (discrete Stokes, verified at N=2,3,4).
# Z_2 CONTROL THROUGHOUT: where a quantity is structurally zero at N=2, the arm must show it.
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
    Erim=elec(st,RIM,N); Espk=elec(st,SPOKE,N)
    def H(th,g2=1.0):
        Hm=sum(np.exp(1j*th)*x+np.exp(-1j*th)*x.conj().T for x in T)
        return -(Hm) - g2*(Erim+Erim.conj().T+Espk+Espk.conj().T)/1.0
    def gs(th,g2=1.0):
        h=H(th,g2); h=(h+h.conj().T)/2
        ev,U=np.linalg.eigh(h); return ev[0],U[:,0]
    print(f"\n  ===== N={N}  (physical dim {D}) =====")
    print(f"  {'theta':>7}{'Im<W_rim>':>13}{'Im<W_tri>':>13}{'J_rim = -dE/dth':>18}{'<E_rim>':>12}")
    for th in (0.0,0.4,1.0,2.0):
        e0,psi=gs(th)
        f=lambda O: complex(psi.conj()@(O@psi))
        h=1e-4
        J=-(gs(th+h)[0]-gs(th-h)[0])/(2*h)
        print(f"  {th:>7.3f}{f(Wrim).imag:>13.6f}{f(T[0]).imag:>13.6f}{J:>18.6f}{f(Erim).real:>12.6f}")
print()
print("  Im<W_rim> is the SENSE OF ROTATION ON THE BOUNDARY -- the thing that read exactly 0.000000")
print("  on every fringe carrier, because a fringe touches no plaquette. The rim touches four.")
