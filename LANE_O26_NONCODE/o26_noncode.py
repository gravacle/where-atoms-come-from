"""O-26 route 2: DO NON-CODE RECORDS HAVE CURVED CONNECTIONS?

C-22 proved records exist with no qubits, no stabiliser group, no cells and no geometry -- in
dimension 6, which is not a power of two. Every curvature test so far has been run INSIDE the code
framework, and C-24 showed stabiliser records are flat NECESSARILY. So the flat result may be a
fact about CODES rather than about RECORDS.

The obstacle: minimality selected the flat connection by WEIGHT, and weight needs a tensor
factorisation that a dim-6 system does not have. So the question becomes: without weight, DOES
ANYTHING SELECT THE CONNECTION? Two natural candidates are measured --
   ||U - I||        how far the writer is from doing nothing
   ||log U||        the ACTION of the operation, the physically motivated one
If neither discriminates, nothing selects flat and curvature is freely available."""
import sys, itertools, numpy as np
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel
def say(*a): print(*a); sys.stdout.flush()
rng=np.random.default_rng(20260819)
say("="*100); say("O-26 ROUTE 2   NON-CODE RECORDS: WHAT SELECTS THE CONNECTION?"); say("="*100)

for dim,lbl in ((8,"C^8, H = 0 (3 records)"), (6,"C^6, H with three degenerate PAIRS -- NOT a power of 2")):
    if dim==8: H=np.zeros((8,8),dtype=complex)
    else:      H=np.diag(np.repeat(rng.normal(size=3),2)).astype(complex)
    m=RecordModel(H,[]); fam,comm,wr=m.independence(m.records())
    say(f"\n  {lbl}:  {len(fam)} independent records, all commuting {bool(comm.all())}")
    if len(fam)<2: say("    fewer than two records -- skipped"); continue
    R=fam[:min(3,len(fam))]
    # joint record basis (powers of two separate the sign patterns -- 1,2,3 would be degenerate)
    Mx=sum((2**i)*Rr for i,Rr in enumerate(R))
    w,V=np.linalg.eigh(Mx)
    labels=[tuple(int(round(np.real(V[:,k].conj()@Rr@V[:,k]))) for Rr in R) for k in range(dim)]
    if len(set(labels))<2**len(R):
        say(f"    records do not resolve the space ({len(set(labels))} labels) -- skipped"); continue
    def writer(idx, ph=None):
        U=np.zeros((dim,dim),dtype=complex)
        for k in range(dim):
            t=list(labels[k]); t[idx]=-t[idx]; j=labels.index(tuple(t))
            U += (ph[k] if ph is not None else 1.0)*np.outer(V[:,j],V[:,k].conj())
        return U
    def ok(U,i):
        if np.linalg.norm(U.conj().T@U-np.eye(dim))>1e-9: return False
        return all(np.linalg.norm(U.conj().T@R[j]@U-((-1 if j==i else 1)*R[j]))<1e-9 for j in range(len(R)))
    def hol(Us):
        return max(float(np.linalg.norm(Us[i]@Us[j]@np.linalg.inv(Us[i])@np.linalg.inv(Us[j])-np.eye(dim)))
                   for i,j in itertools.combinations(range(len(Us)),2))
    def act(U):
        """||log U|| for a unitary, from its eigenphases -- no scipy needed. U = W diag(e^{i.th}) W-dag,
           so log U = W diag(i.th) W-dag and the Frobenius norm is sqrt(sum th^2), with th in (-pi,pi]."""
        ev=np.linalg.eigvals(U)
        th=np.angle(ev)
        return float(np.sqrt((th**2).sum()))
    say(f"    {'writer choice':<26}{'||U - I||':>12}{'||log U||':>12}{'holonomy':>12}{'verdict':>10}")
    canon=[writer(i) for i in range(len(R))]
    if not all(ok(canon[i],i) for i in range(len(R))):
        say("    canonical writers failed verification -- not reporting"); continue
    say(f"    {'canonical (no phases)':<26}{np.linalg.norm(canon[0]-np.eye(dim)):>12.4f}"
        f"{act(canon[0]):>12.4f}{hol(canon):>12.3e}{'FLAT':>10}")
    best=None
    for _ in range(60):
        ph=[np.exp(2j*np.pi*rng.random(dim)) for _ in range(len(R))]
        Us=[writer(i,ph[i]) for i in range(len(R))]
        if not all(ok(Us[i],i) for i in range(len(R))): continue
        h=hol(Us)
        if h>1e-6 and (best is None or h>best[0]): best=(h,Us)
    if best is None:
        say("    no curved admissible alternative found"); continue
    h,Us=best
    say(f"    {'phase-modified (curved)':<26}{np.linalg.norm(Us[0]-np.eye(dim)):>12.4f}"
        f"{act(Us[0]):>12.4f}{h:>12.3e}{'CURVED':>10}")
    dn=abs(np.linalg.norm(Us[0]-np.eye(dim))-np.linalg.norm(canon[0]-np.eye(dim)))
    da=abs(act(Us[0])-act(canon[0]))
    say(f"    -> ||U-I|| discriminates: {'YES' if dn>1e-6 else 'NO  (identical, so it selects nothing)'}")
    say(f"    -> ||log U|| discriminates: {'YES, curved costs more' if da>1e-6 and act(Us[0])>act(canon[0]) else ('YES, curved costs LESS' if da>1e-6 else 'NO')}")
