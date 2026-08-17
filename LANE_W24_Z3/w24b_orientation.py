# DOES ORIENTATION DO ANYTHING?  In Z_2 every Wilson loop is real (+-1) and the two senses of a
# circuit are the same object. In Z_3 the loop is complex and the senses differ. The parameter that
# makes the difference physical is a THETA TERM -- dimensionless, periodic, carrying no scale.
# ONE VARIABLE: theta. Same carrier, same coupling, same everything else. N=2 is the control:
# there the term cannot do anything, because W = W^dagger identically.
import numpy as np
def clock(N): return np.diag([np.exp(2j*np.pi*k/N) for k in range(N)])
def shift(N):
    M=np.zeros((N,N),dtype=complex)
    for k in range(N): M[(k+1)%N,k]=1
    return M
E=[(0,1),(1,2),(2,0),(0,3),(1,4)]; BULK=[0,1,2]; L=len(E)
PLAQ=[0,1,2]                                    # the triangle, traversed 0 -> 1 -> 2 -> 0

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
    for i in PLAQ: W=W@X[i]                      # Wilson loop around the triangle
    return d,Z,X,G,P,W

print("  ground-state Wilson loop as theta is turned on. g^2 = 1.0 throughout.")
print(f"  {'N':>3}{'theta':>9}{'Re<W>':>12}{'Im<W>':>12}{'|<W>|':>10}   note")
for N in (2,3):
    d,Z,X,G,P,W = build(N)
    w_,vec=np.linalg.eigh(P); B=vec[:,np.abs(w_)>0.5]
    Wp=B.conj().T@W@B; Ep=B.conj().T@sum(Z[i]+Z[i].conj().T for i in range(L))@B
    for th in (0.0, 0.4, 1.0, 2.0, np.pi/2):
        Hp = -(np.exp(1j*th)*Wp + np.exp(-1j*th)*Wp.conj().T) - 1.0*Ep
        Hp = (Hp+Hp.conj().T)/2
        ev,U=np.linalg.eigh(Hp); psi=U[:,0]
        wv=complex(psi.conj()@(Wp@psi))
        note = "Z_2: W = W-dagger, Im is identically 0" if N==2 else ("chiral" if abs(wv.imag)>1e-9 else "")
        print(f"  {N:>3}{th:>9.4f}{wv.real:>12.6f}{wv.imag:>12.6f}{abs(wv):>10.6f}   {note}")
    print()
print("  Im<W> != 0 is a SENSE OF ROTATION. It cannot exist at N=2 for a structural reason:")
print("  the Wilson loop is its own adjoint there, so the two senses of the circuit are one object.")
