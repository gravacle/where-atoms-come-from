# FRAME CHALLENGE ON THE FRAMEWORK ITSELF: is GAUGE doing the work, or is CONSTRAINT?
# We chose gauge theory because the physical algebra does not factorize across a cut. But ANY
# constrained system has that. If a plain spin chain with a conserved total magnetization shows the
# same structure, then "gauge theory" is narration and "constraint" is the ingredient.
# ONE VARIABLE: whether the constraint is a LOCAL gauge law or a GLOBAL conservation law.
import numpy as np, itertools
I2=np.eye(2); Xp=np.array([[0,1],[1,0]],dtype=complex); Zp=np.diag([1,-1]).astype(complex)
def op(i,P,n):
    M=np.array([[1]],dtype=complex)
    for j in range(n): M=np.kron(M,P if j==i else I2)
    return M
def S(r):
    w=np.linalg.eigvalsh(r); return float(-sum(x*np.log2(x) for x in w if x>1e-12))
def red(rho,keep,n):
    d=[2]*n; r=rho.reshape(d+d); tr=[i for i in range(n) if i not in keep]
    for k,i in enumerate(sorted(tr,reverse=True)): r=np.trace(r,axis1=i,axis2=i+n-k)
    return r.reshape(2**len(keep),2**len(keep))
rng=np.random.default_rng(20260828)

def probe(tag,n,P,region,frag,extra=""):
    w_,vec=np.linalg.eigh(P); B=vec[:,w_>0.5]; d=B.shape[1]
    vals=[]
    for _ in range(3):
        c=rng.normal(size=d)+1j*rng.normal(size=d); psi=B@c; psi/=np.linalg.norm(psi)
        rho=np.outer(psi,psi.conj())
        hs=S(red(rho,region,n)); hf=S(red(rho,frag,n)); hj=S(red(rho,sorted(region+frag),n))
        vals.append((hs+hf-hj)/hs if hs>1e-12 else float('nan'))
    print(f"  {tag:<44} dim {d:>4}   I(S:F)/H(S) = {'  '.join(f'{v:.6f}' for v in vals)}   {extra}")

print("== ARM A — LOCAL GAUGE CONSTRAINT (Z_2 gauge theory, what we have been doing) ==")
EG=[(0,1),(1,2),(2,0),(2,3),(3,4),(4,2),(4,5),(5,0),(0,4)]
nA=len(EG); ZA=[op(i,Zp,nA) for i in range(nA)]
def gauss(v):
    M=np.eye(2**nA,dtype=complex)
    for i,(a,b) in enumerate(EG):
        if a==v or b==v: M=M@ZA[i]
    return M
PA=np.eye(2**nA,dtype=complex)
for v in range(6): PA=PA@((np.eye(2**nA)+gauss(v))/2)
cutA=[2,7,8]                                              # the links incident to vertex 0, minus link 0
probe("gauge: system = link 0, fragment = its Gauss cut",nA,PA,[0],cutA)

print("\n== ARM B — GLOBAL CONSERVATION CONSTRAINT (no gauge law, no links, no vertices) ==")
nB=9                                                       # nine spins, same dimension count
M=sum(op(i,Zp,nB) for i in range(nB))                      # total magnetization
for target in (1,3):
    PB=np.zeros((2**nB,2**nB),dtype=complex)
    w,v=np.linalg.eigh(M)
    for k in range(len(w)):
        if abs(w[k]-target)<1e-9: PB+=np.outer(v[:,k],v[:,k].conj())
    probe(f"spins: fixed total magnetization M = {target}",nB,PB,[0],[1,2,3,4,5,6,7,8],
          "(fragment = ALL the rest)")

print("\n== AND THE STRUCTURAL COMPARISON THAT DECIDES IT ==")
print("  gauge:  knowing the flux on a CUT determines the link it encloses  -> LOCAL, a surface")
print("  spins:  knowing the magnetization of the COMPLEMENT determines the region -> GLOBAL, no surface")
print("  Both give a non-factorizing algebra. Only one of them makes the determining set a SURFACE.")
print("  THAT is what gauge buys, and it is the only thing it buys here.")
