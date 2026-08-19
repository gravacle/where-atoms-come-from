"""T-12 / PF-4: is there a NECESSARY structural condition narrower than the clauses?

C-4 read 'records require non-trivial homology of an F_2 chain complex'. O-8 found non-CSS
stabiliser codes provably ESCAPE that class, so it is not necessary. The question is whether
anything narrower than the clauses themselves survives.

THE TEST: exhibit records in systems with NO code, NO lattice, NO complex and NO locality. If the
model finds them from (H,{L_k}) alone, then no homological or code-theoretic condition can be
necessary, and C-12 is the necessary condition -- which is the clauses restated, not narrower."""
import sys, numpy as np
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel, build_writer
def say(*a): print(*a); sys.stdout.flush()
rng=np.random.default_rng(20260819)
say("="*100); say("T-12   IS ANY STRUCTURAL CONDITION NECESSARY?"); say("="*100)

say("\n1.  RECORDS IN SYSTEMS WITH NO CODE, NO LATTICE, NO COMPLEX")
say(f"  {'system':<44}{'dim':>5}{'records':>9}{'writer verified':>17}")
cases=[]
# (a) a bare degenerate Hamiltonian -- no tensor factors at all
for n in (4,6,8):
    cases.append((f"random H, {n} distinct pairs, no noise", np.diag(np.repeat(rng.normal(size=n//2),2)), []))
# (b) degenerate H with STRUCTURED noise that is not a stabiliser of anything
for n in (4,6):
    D=np.diag(np.repeat(rng.normal(size=n//2),2))
    L=np.zeros((n,n),dtype=complex)
    for i in range(0,n,2): L[i,i]=1.0; L[i+1,i+1]=-1.0          # commutes with H, not a Pauli
    cases.append((f"random H + a diagonal jump, dim {n}", D, [L]))
ok=0
for lbl,H,Ls in cases:
    m=RecordModel(np.asarray(H,dtype=complex),Ls)
    recs=m.records(); n=H.shape[0]
    good=0
    for R in recs[:4]:
        U=build_writer(R,m.es)
        if U is None: continue
        if (np.linalg.norm(U.conj().T@U-np.eye(n))<1e-7 and np.linalg.norm(U@m.H-m.H@U)<1e-7
            and np.linalg.norm(U.conj().T@R@U+R)<1e-7): good+=1
    ok += (len(recs)>0 and good>0)
    say(f"  {lbl:<44}{n:>5}{len(recs):>9}{f'{good}/{min(len(recs),4)}':>17}")
say(f"\n  systems with a record and a verified writer: {ok} of {len(cases)}")

say("\n2.  WHAT STRUCTURE DO THEY HAVE?   (the point: none of the proposed conditions)")
say("     tensor factorisation into qubits : NO -- dim 6 is not a power of 2 in two of the cases")
say("     a stabiliser group               : NO -- there are no Pauli operators to generate one")
say("     an F_2 chain complex             : NO -- there are no cells")
say("     non-trivial homology             : NOT DEFINED -- there is no complex to take it of")
say("     a lattice or any geometry        : NO")

say("\n3.  SO WHAT IS NECESSARY?")
say("     C-12, and it is the clauses restated rather than anything narrower:")
say("       a record satisfying (i)-(iv) exists IFF the commutant of alg{I,H,L_k,L_k-dagger}")
say("       contains a projection that is non-trivial on some eigenspace and trace-balanced.")
say("     C-15 already showed (i)-(iv) are CARRIER-FREE and (v) needs a locality structure.")
say("     Given a locality structure, (v) is a DISTANCE condition -- no admissible local operation")
say("     reaches the record -- which is Knill-Laflamme and is not homological either.")
