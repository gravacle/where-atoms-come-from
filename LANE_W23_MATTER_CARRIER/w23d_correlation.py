# CORRELATION, NOT PERSISTENCE. The value need not sit still if what survives is the RELATIONSHIP.
# THE TEST ATTEMPT ONE FAILED AT 0.0 (W-06: "carrier and record are never correlated"):
#   prepare two histories differing ONLY in what gets written; evolve BOTH under the SAME H;
#   ask whether a BOUNDARY measurement can still tell them apart later.
# WITH A CONTROL THAT MUST FAIL: two histories differing in something that is NOT written should
# become indistinguishable at the boundary. Without that arm, "distinguishable" means nothing.
import numpy as np
I2=np.eye(2); Xp=np.array([[0,1],[1,0]],dtype=complex); Zp=np.diag([1,-1]).astype(complex)
NL,NV=5,5; n=NL+NV
E=[(0,1),(1,2),(2,0),(0,3),(1,4)]; BULK=[0,1,2]
def q(i,P):
    M=np.array([[1]],dtype=complex)
    for j in range(n): M=np.kron(M,P if j==i else I2)
    return M
Zl=[q(i,Zp) for i in range(NL)]; Xl=[q(i,Xp) for i in range(NL)]
tau=[q(NL+v,Zp) for v in range(NV)]; mu=[q(NL+v,Xp) for v in range(NV)]
def G(v):
    M=tau[v].copy()
    for i,(a,b) in enumerate(E):
        if a==v or b==v: M=M@Zl[i]
    return M
def hop(i):
    a,b=E[i]; return mu[a]@Xl[i]@mu[b]
P=np.eye(2**n,dtype=complex)
for v in BULK: P=P@((np.eye(2**n)+G(v))/2)
w_,B0=np.linalg.eigh(P); B=B0[:,w_>0.5]; d=B.shape[1]
H=-(1.0)*(Xl[0]@Xl[1]@Xl[2])-1.0*sum(Zl)-0.8*sum(hop(i) for i in range(NL))-0.5*sum(tau)
Hp=B.conj().T@H@B; ev,U=np.linalg.eigh(Hp)
rng=np.random.default_rng(20260830)

# the BOUNDARY: link 3 and link 4 plus the two boundary matter sites -> qubits 3,4,8,9
BDRY=[3,4,NL+3,NL+4]
def red_full(psi,keep):
    v=(B@psi).reshape([2]*n)
    r=np.tensordot(v,v.conj(),axes=([i for i in range(n) if i not in keep],)*2)
    m=2**len(keep); return r.reshape(m,m)
def td(a,b):
    return 0.5*np.sum(np.abs(np.linalg.eigvalsh(a-b)))

def run(op_diff,tag):
    c=rng.normal(size=d)+1j*rng.normal(size=d); c/=np.linalg.norm(c)
    c2=B.conj().T@(op_diff@(B@c))
    nrm=np.linalg.norm(c2)
    if nrm<1e-9: return None
    c2/=nrm
    Uc,Uc2=U.conj().T@c,U.conj().T@c2
    out=[]
    for t in (0.0,1.0,5.0,20.0,100.0,400.0):
        a=U@(np.exp(-1j*ev*t)*Uc); b=U@(np.exp(-1j*ev*t)*Uc2)
        out.append(td(red_full(a,BDRY),red_full(b,BDRY)))
    print(f"  {tag:<46}" + "".join(f"{x:>10.5f}" for x in out))
    return out

print("  trace distance of the BOUNDARY reduced states, two histories differing only as labelled")
print(f"  {'':<46}" + "".join(f"{t:>10}" for t in ('t=0','t=1','t=5','t=20','t=100','t=400')))
run(mu[0], "WRITTEN: flip matter at bulk vertex 0")
run(mu[1], "WRITTEN: flip matter at bulk vertex 1")
run(Zl[1], "CONTROL: phase on link 1 (interior, not written)")
run(Zl[0]@Zl[1]@Zl[2], "CONTROL: phase on the whole plaquette")
print()
print("  A boundary that RECORDS keeps the written arms apart while the values themselves wander.")
print("  A boundary that records NOTHING sends every arm to 0. Attempt one measured 0.0 (W-06).")
print("  If the CONTROL arms stay high too, the boundary is not recording -- it is just not thermalising.")
