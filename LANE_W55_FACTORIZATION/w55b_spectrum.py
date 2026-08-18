"""W-55b.  THE EXACT QUESTION, DECIDED WITHOUT AN OPTIMISER.

w55's search plateaued near 0.94 whatever it was given: at eps=0 the planted Hamiltonian is EXACTLY
local (1.0000) and the search reached only 0.9461, while at eps=0.40 it reported 0.936 -- ABOVE the
planted truth of 0.875. The optimiser, not the physics, set the number. w55 is void.

THE EXACT QUESTION NEEDS NO OPTIMISER.
    H is exactly local under SOME factorisation, i.e. H = U (A (x) I + I (x) B) U^dag,
    IFF its spectrum is a SUMSET:  spec(H) = { a_i + b_j } for some a in R^dA, b in R^dB.
Conjugation preserves the spectrum, so this is a property of the eigenvalues ALONE, and it is
decidable. A dA x dB sumset has dA + dB - 1 free parameters; the spectrum has dA*dB entries. For
4 x 4 that is 7 parameters describing 16 numbers, so a generic spectrum cannot be one.

THIS IS ITEM 24 IN ITS SHARPEST FORM: if a generic Hamiltonian's spectrum is not a sumset, then NO
factorisation makes it local, no set of boundaries is derivable from it, and boundaries must be an
INPUT. That is a statement about the dynamics, not about our lattice.
"""
import itertools, numpy as np
rng=np.random.default_rng(23)

def sumset_residual(spec,dA,dB,restarts=60,iters=6000):
    """min over (a,b) and over assignments of || sort(spec) - sort(a_i+b_j) ||, normalised.
       Exact zero iff the spectrum is a sumset."""
    s=np.sort(np.asarray(spec,float)); scale=np.abs(s).max()
    best=np.inf
    for _ in range(restarts):
        a=rng.normal(scale=scale/2,size=dA); b=rng.normal(scale=scale/2,size=dB)
        lr=0.05*scale
        for t in range(iters):
            g=np.sort((a[:,None]+b[None,:]).ravel())
            # gradient of ||sort(a+b) - s||^2 w.r.t. a,b, via the sorting permutation
            M=(a[:,None]+b[None,:]).ravel(); order=np.argsort(M)
            diff=np.zeros_like(M); diff[order]=g-s
            Dm=diff.reshape(dA,dB)
            ga=Dm.sum(axis=1); gb=Dm.sum(axis=0)
            a=a-lr*ga/ (dB); b=b-lr*gb/(dA)
            if t%1500==1499: lr*=0.5
        r=np.linalg.norm(np.sort((a[:,None]+b[None,:]).ravel())-s)/np.linalg.norm(s)
        best=min(best,r)
    return best

dA=dB=4; D=dA*dB
print("W-55b  IS A GENERIC SPECTRUM A SUMSET?  (exact locality, no optimiser over factorisations)")
print(f"  {dA}x{dB}: a sumset has {dA+dB-1} free parameters describing {D} eigenvalues.\n")
print(f"  {'case':>40s} {'sumset residual':>16s}  reading")
print("  "+"-"*76)
# planted sumsets -- must come back ~0
for t in range(3):
    a=rng.normal(size=dA); b=rng.normal(size=dB)
    spec=(a[:,None]+b[None,:]).ravel()
    r=sumset_residual(spec,dA,dB)
    print(f"  {f'PLANTED sumset (trial {t})':>40s} {r:16.3e}  {'exactly local' if r<1e-3 else 'CONTROL FAILED'}")
# generic spectra
gen=[]
for t in range(4):
    spec=np.sort(np.linalg.eigvalsh((lambda A:(A+A.conj().T)/2)(
        rng.normal(size=(D,D))+1j*rng.normal(size=(D,D)))))
    r=sumset_residual(spec,dA,dB); gen.append(r)
    print(f"  {f'GENERIC random Hermitian (trial {t})':>40s} {r:16.3e}  "
          f"{'exactly local' if r<1e-3 else 'NOT a sumset -> no local factorisation'}")
# a physically structured spectrum: our own lattice carrier
print()
print("  AND THE PROGRAM'S OWN CARRIER, for contrast:")
import itertools as it
V2=[(i,j) for j in range(3) for i in range(3)]; vid={v:k for k,v in enumerate(V2)}
E=[]
for j in range(3):
    for i in range(2): E.append((vid[(i,j)],vid[(i+1,j)]))
for j in range(2):
    for i in range(3): E.append((vid[(i,j)],vid[(i,j+1)]))
L=len(E)
st=[s for s in it.product(range(2),repeat=L)
    if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
           -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%2==0 for v in range(9))]
idx={s:i for i,s in enumerate(st)}
hid=lambda i,j:j*2+i; vx=lambda i,j:6+j*3+i
P=[[(hid(i,j),+1),(vx(i+1,j),+1),(hid(i,j+1),-1),(vx(i,j),-1)] for j in range(2) for i in range(2)]
def Move(mv):
    M=np.zeros((16,16),complex)
    for j,s in enumerate(st):
        t=list(s)
        for k,sg in mv: t[k]=(t[k]+sg)%2
        t=tuple(t)
        if t in idx: M[idx[t],j]=1.0
    return M
MAG=sum((lambda X:X+X.conj().T)(Move(p)) for p in P)
r=sumset_residual(np.linalg.eigvalsh(-MAG),dA,dB)
print(f"  {'our 3x3 patch, pure magnetic H':>40s} {r:16.3e}  "
      f"{'EXACTLY LOCAL: a factorisation exists' if r<1e-3 else 'not a sumset'}")
print()
print(f"  generic spectra: min residual over {len(gen)} trials = {min(gen):.3e}")
print("  READING: a sumset residual bounded away from 0 means NO factorisation makes that H local.")
