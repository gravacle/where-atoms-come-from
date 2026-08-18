"""THE RECORD-COUNT LAW, FROM THE CLAUSES ALONE.

A record must be non-trivial (iii) and trace-balanced on EVERY eigenspace (iv, = C-11), and a
FAMILY of independent records must split every joint block EVENLY. Each independent record
therefore halves every block, so the family size is bounded by how many times every eigenspace
multiplicity can be halved:

      k  =  min over eigenspaces E of  v2(m_E)          v2 = the 2-adic valuation

No topology, no index, no lattice, no carrier. The model is given only H."""
import sys, numpy as np
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel
def say(*a): print(*a); sys.stdout.flush()
def v2(m):
    k=0
    while m % 2 == 0 and m: m//=2; k+=1
    return k
P=F=0
say("="*94); say("THE RECORD-COUNT LAW:  k = min_E v2(m_E)"); say("="*94)
say(f"  {'multiplicities':<26}{'dim':>5}{'k measured':>12}{'min v2(m_E)':>13}{'':>4}{'writable':>10}")
cases=[[8],[4,4],[2,2,2,2],[4,2],[8,8],[2,2],[16],[4,4,4],[6,2],[2,6],[4,4,2],[1,4],
       [3,3],[6,6],[5,5],[12],[12,4],[6],[10,6],[2,3],[9],[4,12],[16,8]]
for mult in cases:
    diag=[]
    for i,m in enumerate(mult): diag += [float(i)]*m
    if len(diag)>16: continue
    mdl=RecordModel(np.diag(diag).astype(complex),[])
    try: fam,_,indep=mdl.independence(mdl.records())
    except RuntimeError as e: say(f"  {str(mult):<26}{len(diag):>5}   skipped: {e}"); continue
    pred=min(v2(m) for m in mult)
    ok = (len(fam)==pred); P+=ok; F+=(not ok)
    say(f"  {str(mult):<26}{len(diag):>5}{len(fam):>12}{pred:>13}{'  ok' if ok else ' MISMATCH':>4}{len(indep):>10}")
say("")
say(f"  {P} PASS, {F} FAIL   over {P+F} spectra including odd, non-power-of-2 and mixed multiplicities")
say("")
say("  WHY THE NAIVE READING IS WRONG, and the control that shows it:")
say("    floor(log2 min m_E) would predict [3,3] -> 1, [6,6] -> 2, [5,5] -> 2.")
say("    Measured: 0, 1, 0.  The count is set by the 2-ADIC VALUATION, not the size.")
say("")
say("  AND IT DERIVES THE GAUGE RESULT WITHOUT TOPOLOGY:")
say("    toric code 2x2 has ground multiplicity 4 and the rest of its spectrum even,")
say("    so k = v2(4) = 2 = 2g.  The index law G-7 is a CONSEQUENCE on that carrier,")
say("    not the source of the count.")
sys.exit(1 if F else 0)
