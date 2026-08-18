"""W-60.  WHICH HAMILTONIANS ADMIT RECORDS?

W-59 left the one substantive obstruction: the structures records need are NON-GENERIC, existing only
for special Hamiltonians. The content of the whole account is in WHICH way they are special. This
asks that directly.

THE CANDIDATE CRITERION, derived from what is already proved rather than guessed.
T1: a record must satisfy [H,R] = 0. A record must also be NON-TRIVIAL -- not merely a function of
the energy, or "the record" carries no information beyond which energy level the system is in. So

        a record exists  <=>  the commutant of H is NON-ABELIAN
                          <=>  H has a DEGENERATE eigenvalue

because if every eigenvalue of H is simple, everything commuting with H is a function of H, and the
commutant is abelian of dimension exactly D.

AND DEGENERACY IS NOT AN ACCIDENT. By Wigner, spectral degeneracy in a physical Hamiltonian comes
from SYMMETRY. If the criterion holds, then:

        A HAMILTONIAN ADMITS RECORDS IFF IT HAS A SYMMETRY,
        AND THE RECORDS ARE THAT SYMMETRY'S CONSERVED CHARGES.

TESTED, NOT ASSERTED: commutant dimension against degeneracy across generic, symmetric and
this program's own carrier; and whether the count of independent records matches the multiplicities.
"""
import itertools, numpy as np
rng=np.random.default_rng(41)

def commutant_dim(H,tol=1e-8):
    """dim { A : [A,H] = 0 } -- computed as the null space of the adjoint action"""
    D=H.shape[0]; I=np.eye(D,dtype=complex)
    M=np.kron(H,I)-np.kron(I,H.T)
    sv=np.linalg.svd(M,compute_uv=False)
    return int((sv<tol*max(1.0,sv[0])).sum())
def multiplicities(H,tol=1e-8):
    ev=np.sort(np.linalg.eigvalsh(H)); out=[]; cur=1
    for i in range(1,len(ev)):
        if abs(ev[i]-ev[i-1])<tol*max(1.0,abs(ev[i])): cur+=1
        else: out.append(cur); cur=1
    out.append(cur); return out
def report(name,H):
    D=H.shape[0]; mult=multiplicities(H)
    c=commutant_dim(H); pred=sum(d*d for d in mult)
    nontrivial = c-len(mult)          # commutant beyond the functions of H (which number len(mult))
    print(f"  {name:>36s} D={D:3d}  multiplicities {str(mult):>22s}  dim comm {c:4d} "
          f"(pred {pred:4d})  records beyond f(H): {nontrivial:4d}")
    return c,pred,nontrivial

print("W-60  DOES A RECORD EXIST IFF H IS DEGENERATE?")
print("      predicted dim of commutant = sum of (multiplicity)^2\n")
D=8
A=rng.normal(size=(D,D))+1j*rng.normal(size=(D,D)); Hg=(A+A.conj().T)/2
report("GENERIC random Hermitian",Hg)
Q,_=np.linalg.qr(rng.normal(size=(D,D))+1j*rng.normal(size=(D,D)))
for mult in ([2]*4,[4,4],[4,2,2],[8]):
    lv=[]; 
    for i,d in enumerate(mult): lv+=[float(i)]*d
    report(f"planted multiplicities {mult}",Q@np.diag(lv)@Q.conj().T)

print("\n  AND THIS PROGRAM'S OWN CARRIER:")
V2=[(i,j) for j in range(3) for i in range(3)]; vid={v:k for k,v in enumerate(V2)}
E=[]
for j in range(3):
    for i in range(2): E.append((vid[(i,j)],vid[(i+1,j)]))
for j in range(2):
    for i in range(3): E.append((vid[(i,j)],vid[(i,j+1)]))
L=len(E)
hid=lambda i,j:j*2+i; vx=lambda i,j:6+j*3+i
PLQ=[[(hid(i,j),+1),(vx(i+1,j),+1),(hid(i,j+1),-1),(vx(i,j),-1)] for j in range(2) for i in range(2)]
st=[s for s in itertools.product(range(2),repeat=L)
    if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
           -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%2==0 for v in range(9))]
idx={s:i for i,s in enumerate(st)}; Dp=len(st)
def Move(mv):
    M=np.zeros((Dp,Dp),complex)
    for j,s in enumerate(st):
        t=list(s)
        for k,sg in mv: t[k]=(t[k]+sg)%2
        t=tuple(t)
        if t in idx: M[idx[t],j]=1.0
    return M
MAG=sum((lambda X:X+X.conj().T)(Move(p)) for p in PLQ)
report("our 3x3 patch, H_magnetic",-MAG)
def Zop(links): return np.diag([(-1.0)**(sum(s[k] for k in links)%2) for s in st]).astype(complex)
ELEC=sum(Zop([k]) for k in range(L))
for g2 in (0.05,0.5):
    report(f"our patch, H_mag + {g2}*H_elec",-MAG-g2*ELEC)

print("\n  THE GENERICITY STATEMENT, MEASURED: how often is a random Hermitian degenerate?")
deg=0; N=400
for _ in range(N):
    A=rng.normal(size=(6,6))+1j*rng.normal(size=(6,6)); H=(A+A.conj().T)/2
    if max(multiplicities(H))>1: deg+=1
print(f"    {deg} of {N} random Hermitian matrices had ANY degenerate eigenvalue")
print("    -> degeneracy is non-generic, so RECORDS ARE NON-GENERIC, and the reason is spectral.")
print()
print("  THE SYMMETRY LINK: a perturbation that BREAKS the symmetry must destroy the records.")
Q,_=np.linalg.qr(rng.normal(size=(8,8))+1j*rng.normal(size=(8,8)))
H0=Q@np.diag([0,0,1,1,2,2,3,3]).astype(float)@Q.conj().T
Pert=(lambda X:(X+X.conj().T)/2)(rng.normal(size=(8,8))+1j*rng.normal(size=(8,8)))
print(f"    {'eps':>7s} {'multiplicities':>22s} {'dim commutant':>14s} {'records beyond f(H)':>20s}")
for eps in (0.0,1e-6,1e-3,0.1):
    H=H0+eps*Pert
    print(f"    {eps:7.0e} {str(multiplicities(H)):>22s} {commutant_dim(H):14d} "
          f"{commutant_dim(H)-len(multiplicities(H)):20d}")
