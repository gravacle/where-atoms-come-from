"""Documents WHY the previous attempt failed: record_model.commutant() is SAMPLING-BASED and
   under-spans, which is exactly what inflates multiplicities. Exact prediction beside it."""
import sys
import numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O36_EXACT")
from record_model import commutant, eigenspaces
from o36_exact import (group_dihedral, group_Zn, conjugacy_classes, character_table, build_DG)

out = ["=" * 96,
       "O-36 DIAGNOSTIC: sampled commutant dimension vs the EXACT dimension sum_i m_i^2 (d_i=1) / sum_i m_i^2",
       "=" * 96,
       "  carrier      | eigenspace dims | EXACT dim of joint commutant | SAMPLED (record_model.commutant) | short by",
       "-" * 96]
for gs in (group_dihedral(4), group_Zn(2)):
    els, mul, _, name = gs
    classes, sizes, dims, table, cls_of, e = character_table(els, mul)
    _, _, inv, idx = conjugacy_classes(els, mul)
    H, A, B, As = build_DG(els, mul, idx, inv)
    N = H.shape[0]; es = eigenspaces(H)
    di = np.round(dims.real).astype(int)
    exact = 0
    edims = []
    for val, P, m in es:
        chiE = np.array([np.trace(P @ g) for g in As])
        # .real BEFORE float(): float(complex) emitted a ComplexWarning on stderr every run, and
        # reproduce.sh captures stderr (2>&1) -- the cast already used the real part, value bit-identical.
        mult = [int(round(float((((np.conj(table[r, cls_of[i]]) * chiE[i]).sum()
                                  if False else sum(np.conj(table[r, cls_of[i]]) * chiE[i]
                                                    for i in range(len(els)))) / len(els)).real)))
                for r in range(len(di))]
        exact += sum(x * x for x in mult)
        edims.append(int(m))
    gens = [np.eye(N, dtype=complex), H] + list(As)
    sampled = len(commutant(gens))
    out.append("  %-12s | %-15s | %28d | %32d | %d"
               % ("D(%s)" % name, str(edims), exact, sampled, exact - sampled))
out.append("-" * 96)
out.append("  A sampled basis that is SHORT does not merely lose dimensions -- fitting a von Neumann block")
out.append("  decomposition to it returns multiplicities that do not satisfy sum_i d_i*m_i = dim E.")
out.append("  That is the 20-vs-18 and 26-vs-22 seen in the previous attempt. Character theory has no such mode.")
out.append("=" * 96)
txt = "\n".join(out); print(txt)
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O36_EXACT/o36_whyfail.txt", "w").write(txt + "\n")
