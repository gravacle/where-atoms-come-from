"""O-26b: IS 'LEAST ACTION SELECTS CURVATURE' SYSTEMATIC, OR ONE LUCKY DRAW?

o26_noncode found the canonical FLAT writer at ||log U|| = 6.2832 = 2*pi and a curved one at
5.0921 -- curvature CHEAPER by action, the opposite of what weight-minimality suggested. That was
a single draw, chosen for MAXIMUM holonomy, so it proves nothing about the relationship.

Here: sweep many admissible writers, record (action, holonomy) for each, and ask whether they are
related. Also run a carrier of dimension 24 -- NOT a power of two, so not a qubit system -- with
three records, to check this is not an artefact of qubits."""
import sys, itertools, numpy as np
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel
def say(*a): print(*a); sys.stdout.flush()
rng=np.random.default_rng(11)
def act(U):
    return float(np.sqrt((np.angle(np.linalg.eigvals(U))**2).sum()))
def run(H,lbl,ntrial=400):
    dim=H.shape[0]
    m=RecordModel(H,[]); fam,comm,wr=m.independence(m.records())
    say(f"\n  {lbl}   dim {dim}   records {len(fam)}")
    if len(fam)<3: say("    fewer than three records -- cannot compare paths"); return
    R=fam[:3]
    Mx=sum((2**i)*Rr for i,Rr in enumerate(R))
    w,V=np.linalg.eigh(Mx)
    labels=[tuple(int(round(np.real(V[:,k].conj()@Rr@V[:,k]))) for Rr in R) for k in range(dim)]
    if len(set(labels))<8: say(f"    records resolve only {len(set(labels))} classes -- skipped"); return
    def writer(idx,ph=None):
        U=np.zeros((dim,dim),dtype=complex)
        for k in range(dim):
            t=list(labels[k]); t[idx]=-t[idx]; j=labels.index(tuple(t))
            U += (ph[k] if ph is not None else 1.0)*np.outer(V[:,j],V[:,k].conj())
        return U
    def ok(U,i):
        return (np.linalg.norm(U.conj().T@U-np.eye(dim))<1e-9 and
                all(np.linalg.norm(U.conj().T@R[j]@U-((-1 if j==i else 1)*R[j]))<1e-9 for j in range(3)))
    def hol(Us):
        return max(float(np.linalg.norm(Us[i]@Us[j]@np.linalg.inv(Us[i])@np.linalg.inv(Us[j])-np.eye(dim)))
                   for i,j in itertools.combinations(range(3),2))
    canon=[writer(i) for i in range(3)]
    if not all(ok(canon[i],i) for i in range(3)): say("    canonical writers failed -- not reporting"); return
    A0,H0v=act(canon[0]),hol(canon)
    say(f"    canonical (flat):  action {A0:.4f}   holonomy {H0v:.2e}")
    pts=[]
    for _ in range(ntrial):
        ph=[np.exp(2j*np.pi*rng.random(dim)) for _ in range(3)]
        Us=[writer(i,ph[i]) for i in range(3)]
        if not all(ok(Us[i],i) for i in range(3)): continue
        pts.append((sum(act(U) for U in Us)/3.0, hol(Us)))
    if not pts: say("    no admissible alternatives"); return
    A=np.array([p[0] for p in pts]); Hh=np.array([p[1] for p in pts])
    say(f"    {len(pts)} admissible alternatives sampled")
    say(f"      action   range [{A.min():.4f}, {A.max():.4f}]   canonical {A0:.4f}")
    say(f"      holonomy range [{Hh.min():.3e}, {Hh.max():.3e}]")
    lo=Hh[A<np.percentile(A,25)]; hi=Hh[A>np.percentile(A,75)]
    say(f"      mean holonomy of the LOWEST-action quartile : {lo.mean():.4f}")
    say(f"      mean holonomy of the HIGHEST-action quartile: {hi.mean():.4f}")
    r=float(np.corrcoef(A,Hh)[0,1])
    say(f"      correlation(action, holonomy) = {r:+.4f}")
    k=int((A<A0).sum())
    say(f"      admissible writers with action BELOW the flat one: {k} of {len(pts)}"
        f"   ({'flat is NOT least-action' if k>0 else 'flat IS least-action'})")
say("="*100); say("O-26b   IS LEAST ACTION SELECTING CURVATURE?"); say("="*100)
run(np.zeros((8,8),dtype=complex), "C^8, H = 0 (qubits)")
# dim 24, NOT a power of two: three eigenvalues each of multiplicity 8 -> k = min v2(m_E) = 3
d=np.repeat(rng.normal(size=3),8)
run(np.diag(d).astype(complex), "C^24, three 8-fold shells -- NOT a power of two")
