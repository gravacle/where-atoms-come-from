# DOES THE CHIRALITY REACH THE BOUNDARY, OR IS IT A BULK PLAQUETTE EFFECT?
# theta lives on the PLAQUETTE, which is entirely in the bulk (links 0,1,2). The dangling links
# (3,4) touch no plaquette at all. So there is no reason of construction for theta to be visible
# at the boundary -- and a reason to think it should not be. That makes this test able to fail.
# ONE VARIABLE: theta. Carrier, coupling, cut, region, code path all fixed.
import numpy as np
def clock(N): return np.diag([np.exp(2j*np.pi*k/N) for k in range(N)])
def shift(N):
    M=np.zeros((N,N),dtype=complex)
    for k in range(N): M[(k+1)%N,k]=1
    return M
E=[(0,1),(1,2),(2,0),(0,3),(1,4)]; BULK=[0,1,2]; L=len(E); PLAQ=[0,1,2]
BDRY_LINKS=[3,4]
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
        acc=np.eye(d,dtype=complex); Pv=np.zeros((d,d),dtype=complex); Gv=G(v)
        for k in range(N): Pv+=acc; acc=acc@Gv
        P=P@(Pv/N)
    W=np.eye(d,dtype=complex)
    for i in PLAQ: W=W@X[i]
    return d,Z,X,G,P,W

for N in (2,3):
    d,Z,X,G,P,W=build(N)
    w_,vec=np.linalg.eigh(P); B=vec[:,np.abs(w_)>0.5]
    Wp=B.conj().T@W@B; Ep=B.conj().T@sum(Z[i]+Z[i].conj().T for i in range(L))@B
    # boundary observables: the electric field on each dangling link, and their PRODUCT
    Zb=[B.conj().T@Z[i]@B for i in BDRY_LINKS]
    Zprod=B.conj().T@(Z[3]@Z[4].conj().T)@B          # a relative, oriented boundary quantity
    print(f"\n  N={N}   theta lives on the bulk plaquette; the dangling links touch NO plaquette")
    print(f"  {'theta':>8}{'Im<W> (bulk)':>15}{'Im<Z_3> (bdry)':>17}{'Im<Z_4> (bdry)':>17}{'Im<Z_3 Z_4-dag>':>18}")
    for th in (0.0,0.4,1.0,2.0):
        Hp=-(np.exp(1j*th)*Wp+np.exp(-1j*th)*Wp.conj().T)-1.0*Ep
        Hp=(Hp+Hp.conj().T)/2
        ev,U=np.linalg.eigh(Hp); psi=U[:,0]
        f=lambda O: complex(psi.conj()@(O@psi))
        print(f"  {th:>8.4f}{f(Wp).imag:>15.6f}{f(Zb[0]).imag:>17.6f}{f(Zb[1]).imag:>17.6f}{f(Zprod).imag:>18.6f}")
print()
print("  Im on a BOUNDARY observable would mean the sense of rotation reaches the surface.")
print("  Zero there with nonzero in the bulk would mean the chirality is confined to the plaquette.")
