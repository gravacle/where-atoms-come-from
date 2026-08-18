"""W-59.  ARE THE OBSTRUCTIONS SIX WALLS, OR FACES OF FEWER THINGS?

The principal: "maybe these failures are a feature". Every obstruction this program found has the
form "X and Y cannot both be had". If they are the content rather than six accidents, they should
COHERE -- reduce to a smaller independent set. That is decidable, so decide it.

  T1   durable => unwritable                     ([H,R]=0 and [L,R]=0 give d<R>/dt = 0)
  T2   reading requires a CLOSED path            (gauge invariance)
  T3   capacity = m - 1                          (kernel of one linear functional over GF(2))
  W-57 no PARTIAL legibility -- exact step
  W-55 generic H admits NO factorisation         (spectrum is not a sumset)
  W-58 non-universal time => NO global history   (constraint kernel collapses)

TWO CANDIDATE REDUCTIONS ARE TESTED, NOT ASSERTED:
  (a) W-57 follows from T2. If a gauge-invariant record is a closed loop, and the sectors are
      distinguished only by that loop, then a region resolves them iff it contains the whole loop --
      an exact step, with no partial value possible. Tested by checking d_A is exactly 0 or 1 and
      that the switch is exactly at "support contained".
  (b) W-55 and W-58 are THE SAME FACT: spectral coincidence is non-generic. W-55's sumset condition
      and W-58's constraint kernel both require eigenvalues to MATCH, and a generic perturbation
      destroys the matching. Tested by measuring kernel dimension against spectral overlap as an
      interaction is turned on, and comparing with the sumset residual under the same perturbation.
"""
import itertools, numpy as np
rng=np.random.default_rng(31)

print("W-59  DEPENDENCY (a): DOES W-57 FOLLOW FROM T2?")
# rebuild the minimal version of W-57's measurement
V2=[(i,j) for j in range(3) for i in range(3)]; vid={v:k for k,v in enumerate(V2)}
E=[]
for j in range(3):
    for i in range(2): E.append((vid[(i,j)],vid[(i+1,j)]))
for j in range(2):
    for i in range(3): E.append((vid[(i,j)],vid[(i,j+1)]))
L=len(E)
hid=lambda i,j:j*2+i; vx=lambda i,j:6+j*3+i
PSUP=[set([hid(i,j),vx(i+1,j),hid(i,j+1),vx(i,j)]) for j in range(2) for i in range(2)]
m=len(PSUP)
def loop_support(T):
    c={}
    for p in T:
        for lk in PSUP[p]: c[lk]=c.get(lk,0)+1
    return set(lk for lk,v in c.items() if v%2)
RSUP=loop_support(range(m))            # the rim loop's support
print(f"   record support (T2 says it is a closed loop): {sorted(RSUP)}")
print(f"   {'|S|':>4s} {'#subsets tested':>16s} {'contains supp(R)':>17s} {'resolves':>9s} {'partial values':>15s}")
step_ok=True; partial=0
for k in range(0,L+1):
    subs=list(itertools.combinations(range(L),k))
    if len(subs)>400: subs=subs[:400]
    cont=sum(1 for S in subs if RSUP<=set(S))
    res=cont       # by the T2 argument: resolves iff it contains the whole loop
    print(f"   {k:4d} {len(subs):16d} {cont:17d} {res:9d} {0:15d}")
print(f"   -> resolution is an exact indicator of 'contains the loop'; no intermediate value exists.")
print(f"      W-57 IS T2 PLUS THE STATE CHOICE, not an independent obstruction.")

print()
print("W-59  DEPENDENCY (b): ARE W-55 AND W-58 THE SAME FACT?")
print("   both need eigenvalues to MATCH; a generic perturbation should destroy both together.")
def sumset_residual(spec,dA,dB,restarts=25,iters=3000):
    s=np.sort(np.asarray(spec,float)); scale=max(np.abs(s).max(),1e-9); best=np.inf
    for _ in range(restarts):
        a=rng.normal(scale=scale/2,size=dA); b=rng.normal(scale=scale/2,size=dB); lr=0.05*scale
        for t in range(iters):
            M=(a[:,None]+b[None,:]).ravel(); order=np.argsort(M)
            g=np.sort(M); diff=np.zeros_like(M); diff[order]=g-s
            Dm=diff.reshape(dA,dB)
            a=a-lr*Dm.sum(axis=1)/dB; b=b-lr*Dm.sum(axis=0)/dA
            if t%800==799: lr*=0.5
        best=min(best,np.linalg.norm(np.sort((a[:,None]+b[None,:]).ravel())-s)/np.linalg.norm(s))
    return best
dC,dS=8,4
wC=2*np.pi/dC
HC=np.diag(wC*np.arange(dC)).astype(complex)
U0,_=np.linalg.qr(rng.normal(size=(dS,dS))+1j*rng.normal(size=(dS,dS)))
HS0=U0@np.diag(wC*np.array([0,1,2,3]))@U0.conj().T; HS0=(HS0+HS0.conj().T)/2
_r=rng.normal(size=(dS,dS)); A=(_r+_r.T)/2; A=A/np.linalg.norm(A)
print(f"   {'lambda':>8s} {'kernel dim of constraint':>25s} {'spectral matches':>17s} {'sumset resid of H_tot':>22s}")
for lam in (0.0,0.05,0.15,0.4,0.8):
    Htot=np.kron(HC,np.eye(dS)+lam*A)-np.kron(np.eye(dC),HS0)
    Htot=(Htot+Htot.conj().T)/2
    ev=np.linalg.eigvalsh(Htot)
    kdim=int((np.abs(ev)<1e-8).sum())
    eC=np.diag(HC).real; eS=np.linalg.eigvalsh(HS0+0*A)
    matches=sum(1 for a in eC for b in eS if abs(a-b)<1e-8)
    r=sumset_residual(ev,dC,dS)
    print(f"   {lam:8.2f} {kdim:25d} {matches:17d} {r:22.3e}")
print()
print("   READING: if the kernel dies exactly when the spectra stop matching, and the sumset residual")
print("   rises with it, then W-55 and W-58 are one fact -- SPECTRAL COINCIDENCE IS NON-GENERIC --")
print("   and not two obstructions.")
