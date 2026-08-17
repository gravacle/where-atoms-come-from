# A CURRENT NEEDS A LOOP.  <psi|i[H,Q]|psi> = 0 in any eigenstate, so that measurement is forced.
# The real object is a PERSISTENT CURRENT: J = -dE0/dPhi, the ground state's response to a flux
# threaded through a CYCLE. And a flux can only thread a cycle.
# THE STRUCTURAL POINT: our "boundary" is a set of DANGLING LINKS. Dangling links lie on no cycle.
# So no flux parameter exists there, no current can circulate there, and NOTHING CAN FLOW ON IT.
import numpy as np
def clock(N): return np.diag([np.exp(2j*np.pi*k/N) for k in range(N)])
def shift(N):
    M=np.zeros((N,N),dtype=complex)
    for k in range(N): M[(k+1)%N,k]=1
    return M
def cycles_through(edges,i):
    a,b=edges[i]; adj={}
    for j,(u,v) in enumerate(edges):
        if j==i: continue
        adj.setdefault(u,[]).append(v); adj.setdefault(v,[]).append(u)
    seen={a}; st=[a]
    while st:
        x=st.pop()
        for y in adj.get(x,[]):
            if y not in seen: seen.add(y); st.append(y)
    return b in seen

E=[(0,1),(1,2),(2,0),(0,3),(1,4)]; BULK=[0,1,2]; L=len(E); PLAQ=[0,1,2]
print("  which links lie on a cycle, hence can be threaded by a flux?")
for i in range(L):
    print(f"    link {i} = {E[i]}   on a cycle: {cycles_through(E,i)}"
          + ("   <- BOUNDARY (dangling)" if i in (3,4) else ""))
print("  -> the two boundary links lie on NO cycle. There is no flux parameter for them,")
print("     therefore no persistent current, therefore nothing to flow. BY CONSTRUCTION.\n")

def build(N,th):
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
    w_,vec=np.linalg.eigh(P); B=vec[:,np.abs(w_)>0.5]
    Wp=B.conj().T@W@B; Ep=B.conj().T@sum(Z[i]+Z[i].conj().T for i in range(L))@B
    Hp=-(np.exp(1j*th)*Wp+np.exp(-1j*th)*Wp.conj().T)-1.0*Ep
    return float(np.linalg.eigvalsh((Hp+Hp.conj().T)/2)[0])

print("  THE BULK PERSISTENT CURRENT, J = -dE0/dtheta, by central difference (h=1e-4):")
print(f"  {'theta':>8}{'N=2':>14}{'N=3':>14}")
h=1e-4
for th in (0.0,0.4,1.0,2.0):
    row=[]
    for N in (2,3):
        row.append(-(build(N,th+h)-build(N,th-h))/(2*h))
    print(f"  {th:>8.4f}{row[0]:>14.6f}{row[1]:>14.6f}")
print()
print("  A nonzero J at N=3 is a genuine persistent current circulating the bulk plaquette.")
print("  At N=2 it must vanish where the two senses coincide.")
print("  AND THERE IS NO BOUNDARY COLUMN TO PRINT, because the boundary has no cycle to thread.")
