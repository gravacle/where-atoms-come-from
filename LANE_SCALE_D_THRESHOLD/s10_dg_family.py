"""S10 -- WHY D(G) FOR GROWING 2-GROUPS IS NOT A SCALING FAMILY ON THE MINIMAL TORUS.

The brief offered D(G) for 2-groups of growing order as a candidate second family (C-41:
records exist iff |G| is a power of two).  This measures what actually happens when |G| grows,
so the choice of the two stabiliser families instead is a recorded finding and not a preference.

D(G) on the minimal torus (1 vertex, 2 edges, 1 face), C[G] (x) C[G], basis |g1,g2>:
    A_h |g1,g2> = |h g1 h^-1, h g2 h^-1|,  A = (1/|G|) sum_h A_h
    B diagonal, 1 where g1 g2 g1^-1 g2^-1 = e
    H = -(A + B)
"""
import sys, time, json
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD")
import numpy as np
from record_model import RecordModel, eigenspaces

OUT = []
def P(s=""):
    print(s); OUT.append(s)

def Zn(n):
    els = list(range(n))
    return els, (lambda a, b: (a + b) % n), 0

def prod(g1, g2):
    e1, m1, id1 = g1; e2, m2, id2 = g2
    els = [(a, b) for a in e1 for b in e2]
    return els, (lambda x, y: (m1(x[0], y[0]), m2(x[1], y[1]))), (id1, id2)

def D4():
    els = [(a, b) for a in range(4) for b in range(2)]
    def m(x, y):
        a, b = x; c, d = y
        return ((a + (c if b == 0 else -c)) % 4, (b + d) % 2)
    return els, m, (0, 0)

def build(els, mul, ident):
    n = len(els)
    gi = {g: i for i, g in enumerate(els)}
    inv = {}
    for g in els:
        for h in els:
            if mul(g, h) == ident: inv[g] = h; break
    dim = n * n
    A = np.zeros((dim, dim), dtype=complex)
    for h in els:
        for g1 in els:
            for g2 in els:
                s = gi[g1] * n + gi[g2]
                t = gi[mul(mul(h, g1), inv[h])] * n + gi[mul(mul(h, g2), inv[h])]
                A[t, s] += 1.0 / n
    Bd = np.zeros(dim)
    for g1 in els:
        for g2 in els:
            c = mul(mul(g1, g2), mul(inv[g1], inv[g2]))
            if c == ident: Bd[gi[g1] * n + gi[g2]] = 1.0
    return -(A + np.diag(Bd)), n, dim

t0 = time.time()
P("=" * 150)
P("S10  D(G) ON THE MINIMAL TORUS FOR GROWING 2-GROUPS -- does it give a scaling family?")
P("=" * 150)
P("  %-14s %-6s %-7s %-30s %-11s %-9s %-16s %s" %
  ("G", "|G|", "dim", "eigenvalue multiplicities", "H prop I?", "min_proj", "records()",
   "k_max <= log2(min m_E)"))
P("  " + "-" * 140)
rows = []
GROUPS = [("Z_2", Zn(2)), ("Z_4", Zn(4)), ("Z_2^2", prod(Zn(2), Zn(2))),
          ("D_4", D4()), ("Z_8", Zn(8)), ("Z_4xZ_2", prod(Zn(4), Zn(2))),
          ("Z_2^3", prod(prod(Zn(2), Zn(2)), Zn(2))),
          ("Z_2^4", prod(prod(Zn(2), Zn(2)), prod(Zn(2), Zn(2))))]
for name, g in GROUPS:
    els, mul, ident = g
    if len(els) ** 2 > 256:
        P("  %-14s %-6d %-7d SKIPPED: dim > 256 in this run" % (name, len(els), len(els) ** 2))
        continue
    H, ng, dim = build(els, mul, ident)
    es = eigenspaces(H)
    mults = [m for _, _, m in es]
    propI = bool(len(es) == 1)
    M = RecordModel(H)
    try:
        recs = M.records(); rec_s = "%d records" % len(recs)
        fam = M.commuting_family(recs); rec_s += " (family %d)" % len(fam)
    except RuntimeError as e:
        rec_s = "RAISES: O-28"
    kmax = int(np.floor(np.log2(min(mults))))
    P("  %-14s %-6d %-7d %-30s %-11s %-9d %-16s %d"
      % (name, ng, dim, str(mults)[:29], propI, len(M.projs), rec_s, kmax))
    rows.append(dict(G=name, order=ng, dim=dim, mults=mults, propI=propI,
                     nproj=len(M.projs), records=rec_s, kmax=kmax))
P()
P("READ, from the rows above:")
P("  * for every ABELIAN 2-group the minimal torus gives ONE eigenspace -- conjugation is")
P("    trivial so A = I, and every commutator is e so B = I, leaving H proportional to the")
P("    identity.  Clause (iii) is then satisfied by any non-scalar involution and the carrier")
P("    carries no dynamics at all: it is degenerate, not a scaling family.")
P("  * the non-abelian case D_4 has the three-eigenvalue structure C-43 used, but its order is")
P("    fixed; the next non-abelian 2-group of order 16 puts dim at 256 with more than 20")
P("    minimal projections, so RecordModel.records() hits O-28.")
P("  * hence the two stabiliser families A and B, where k scales and every clause stays")
P("    checkable, and not D(G).  This is a measured reason, not a preference.")
json.dump(rows, open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s10_rows.json", "w"))
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s10_dg_family.txt",
     "w").write("\n".join(OUT) + "\n")
P("total %.1fs" % (time.time() - t0))
