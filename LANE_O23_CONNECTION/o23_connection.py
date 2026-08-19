"""O-23, prior question: IS THE CONNECTION DETERMINED BY THE RECORDS?

T-20 measured holonomy using the minimal-weight Pauli writer and got exactly flat. But a writer is
only fixed up to a unitary that preserves EVERY record label: if U_A flips R_A and fixes the rest,
so does U_A * V for any V diagonal in the joint record basis. If different admissible choices give
different holonomy, then 'is the record geometry curved' is NOT WELL-POSED on the records alone --
the connection is extra data, exactly as it is in gauge theory.

Carrier: H = 0 on C^8, three independent records by the count law -- the most permissive case there
is, so if the writer is unique anywhere it is unique here."""
import sys, itertools, numpy as np
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel
def say(*a): print(*a); sys.stdout.flush()
rng=np.random.default_rng(20260819)
n=8
H=np.zeros((n,n),dtype=complex)
m=RecordModel(H,[])
fam,comm,writable=m.independence(m.records())
say("="*100); say("O-23   IS THE CONNECTION DETERMINED BY THE RECORDS?"); say("="*100)
say(f"  carrier H = 0 on C^{n}: {len(fam)} independent records, all commuting = {bool(comm.all())},")
say(f"  independently writable = {len(writable)} of {len(fam)}")
R=fam[:3]
# joint record basis: simultaneous eigenvectors of the three records
# coefficients must SEPARATE all 2^3 sign patterns. 1,2,3 does not -- 1+2-3 = -1-2+3 = 0 -- so
# the combination is degenerate and eigh does not return a simultaneous eigenbasis. Powers of two do.
Mx=sum((2**i)*Rr for i,Rr in enumerate(R))
w,V=np.linalg.eigh(Mx)
labels=[tuple(int(round(np.real(V[:,k].conj()@Rr@V[:,k]))) for Rr in R) for k in range(n)]
say(f"  joint record labels (all distinct = the records resolve the space): "
    f"{len(set(labels))} of {n} distinct")
def writer_for(idx, extra=None):
    """a unitary flipping record idx and fixing the others; `extra` is any label-preserving phase"""
    U=np.zeros((n,n),dtype=complex)
    for k in range(n):
        tgt=list(labels[k]); tgt[idx]=-tgt[idx]; tgt=tuple(tgt)
        j=labels.index(tgt)
        ph = extra[k] if extra is not None else 1.0
        U += ph * np.outer(V[:,j], V[:,k].conj())
    return U
def check(U, idx):
    oks=[np.linalg.norm(U.conj().T@R[i]@U - (-1 if i==idx else 1)*R[i]) for i in range(3)]
    return max(oks) < 1e-9 and np.linalg.norm(U.conj().T@U-np.eye(n)) < 1e-9
say("")
say("1.  THE CANONICAL CHOICE -- a plain permutation of the joint record basis")
U=[writer_for(i) for i in range(3)]
verified=all(check(U[i],i) for i in range(3))
say(f"    all three verified as writers: {verified}")
if not verified:
    say("    NOT REPORTING a holonomy on operators that are not writers."); sys.exit(1)
hol=[float(np.linalg.norm(U[i]@U[j]@np.linalg.inv(U[i])@np.linalg.inv(U[j])-np.eye(n)))
     for i,j in itertools.combinations(range(3),2)]
say(f"    closed-loop holonomies ||H - I|| = {[f'{h:.3e}' for h in hol]}   "
    f"-> {'FLAT' if max(hol)<1e-9 else 'CURVED'}")
say("")
say("2.  ANOTHER ADMISSIBLE CHOICE -- the same writers times a label-preserving PHASE")
say("    (every one still flips its own record and fixes the others: verified)")
best=0.0; found=None; nadm=0
for trial in range(200):
    ex=[np.exp(2j*np.pi*rng.random(n)) for _ in range(3)]
    U2=[writer_for(i, ex[i]) for i in range(3)]
    if not all(check(U2[i],i) for i in range(3)): continue
    nadm+=1
    h=[float(np.linalg.norm(U2[i]@U2[j]@np.linalg.inv(U2[i])@np.linalg.inv(U2[j])-np.eye(n)))
       for i,j in itertools.combinations(range(3),2)]
    if max(h)>best: best=max(h); found=h
say(f"    admissible choices FOUND: {nadm} of 200")
if nadm==0:
    say("    NO admissible alternative was found -- the loop never executed a valid case, so no")
    say("    conclusion can be drawn from it (D-8). Not reporting."); sys.exit(1)
say(f"    over {nadm} admissible phase choices, largest closed-loop holonomy = {best:.4f}")
if found: say(f"    that choice's three loops: {[f'{x:.4f}' for x in found]}")
say("")
say("3.  READ")
if best>1e-6:
    say("    THE CONNECTION IS NOT DETERMINED BY THE RECORDS. The same record family, with writers")
    say("    that every one of the five clauses accepts, gives holonomy 0 for one choice and")
    say(f"    {best:.3f} for another. 'Is the record geometry curved' is NOT WELL-POSED on the")
    say("    records alone: a connection is EXTRA DATA, exactly as it is in gauge theory.")
    say("    So O-23 is the wrong question. The question is WHAT DETERMINES THE CONNECTION.")
else:
    say("    Every admissible choice gives the same holonomy: the connection IS determined by the")
    say("    records, and T-20's flat result is a property of the record structure itself.")
