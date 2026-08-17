# THE PRINCIPAL'S REFRAME: free oscillation is the ALLOW signature, not a failure.
# A RECORD is a fourth modality neither allow, require nor forbid: ALLOWED BUT HELD.
# The value CAN change and DOESN'T. Measure whether any regime of this carrier exhibits it.
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
w_,vec=np.linalg.eigh(P); B=vec[:,w_>0.5]
Qm=B.conj().T@(tau[0]@tau[3])@B                    # matter charge in region {0,3}
def Hp(g2,J,m):
    plaq=Xl[0]@Xl[1]@Xl[2]
    H=-(1.0/g2)*plaq-g2*sum(Zl)-J*sum(hop(i) for i in range(NL))-m*sum(tau)
    return B.conj().T@H@B
rng=np.random.default_rng(20260829)

def probe(g2,J,m,T=400.0,K=800):
    """prepare the state in a definite charge sector, then evolve. Does it STAY?"""
    h=Hp(g2,J,m); ev,U=np.linalg.eigh(h)
    # project onto Q_matter = +1 and take a random state there: the value is WRITTEN at t=0
    wq,vq=np.linalg.eigh(Qm); sub=vq[:,wq>0.5]
    c=sub@(rng.normal(size=sub.shape[1])+1j*rng.normal(size=sub.shape[1])); c/=np.linalg.norm(c)
    ts=np.linspace(0,T,K); vals=[]
    Uc=U.conj().T@c
    for t in ts:
        ct=U@(np.exp(-1j*ev*t)*Uc); vals.append(float(np.real(ct.conj()@(Qm@ct))))
    v=np.array(vals)
    return v[0], v.mean(), v.min(), v.std()

print("  The value is WRITTEN to +1 at t=0, then left alone. Does it stay written?")
print(f"  {'g^2':>6}{'J':>6}{'m':>6}   {'<Q>(0)':>9}{'time-avg':>10}{'min over T':>12}{'std':>9}   reading")
for (g2,J,m) in [(1.0,0.8,0.5),(1.0,0.8,3.0),(1.0,0.2,3.0),(1.0,0.05,5.0),
                 (0.3,0.05,5.0),(3.0,0.05,5.0),(1.0,0.01,8.0),(1.0,0.0,5.0)]:
    q0,av,mn,sd=probe(g2,J,m)
    read = "FROZEN (J=0, forced)" if J==0 else ("HELD" if mn>0.5 else ("partly held" if av>0.3 else "free"))
    print(f"  {g2:>6.2f}{J:>6.2f}{m:>6.2f}   {q0:>9.4f}{av:>10.4f}{mn:>12.4f}{sd:>9.4f}   {read}")
print()
print("  FREE  : the value wanders the whole allowed range -> ALLOW-side, nothing holds it.")
print("  HELD  : the value CAN move (J>0, the commutator is nonzero) and DOESN'T -> the fourth")
print("          modality, and the only one that is a record.")
print("  FROZEN: J=0, the commutator vanishes -- pinned by the constraint, a LABEL not a record.")
