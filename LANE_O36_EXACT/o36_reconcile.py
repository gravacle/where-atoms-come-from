"""
O-36 ADDENDUM.

(1) RECONCILE WITH C-43 AT MATCHED SAMPLE SIZE.  C-43 found 40/40 SAMPLED records on D(D_4) moved
    by transport.  Here: 40 transport-BLIND records (random splittings of the commutant of H alone)
    beside 40 transport-FIXED records (built inside the commutant of {H} u {A_h} by exact group
    averaging), in the SAME table.  If 40/40 blind move and 0/40 fixed move, C-43 measured the
    sampling measure, not an obstruction.

(2) WHY.  The gauge action A_h on C[G] (x) C[G] is by simultaneous conjugation, so every CENTRAL
    element acts as the identity: the representation factors through G/Z(G).  For D_4 and Q_8,
    G/Z(G) = Z_2 x Z_2 is ABELIAN, so only 1-dimensional irreps can occur, every commutant block
    has d_i = 1, and any trace of the right parity -- including 0 -- is reachable.  Checked below,
    and checked NOT to be the whole story at larger scale (D-17): for D_8, G/Z(G) = D_4 is
    non-abelian, 2-dimensional blocks DO appear, and 0 is still reachable.
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import eigenspaces
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O36_EXACT")
from o36_exact import (group_dihedral, group_Zn, group_Q8, group_ZmxZn, conjugacy_classes,
                       character_table, build_DG, gauge_average, block_project,
                       minimal_projections_in, eigenbasis_of, signs_to_zero, verify_record)

def records_sample(res_els, mul, name, n_draw=40):
    els = res_els
    classes, sizes, dims, table, cls_of, e = character_table(els, mul)
    _, _, inv, idx = conjugacy_classes(els, mul)
    H, A, B, As = build_DG(els, mul, idx, inv)
    N = H.shape[0]
    es = eigenspaces(H)
    I = np.eye(N)

    def draw(seed, fixed):
        rng = np.random.default_rng(seed)
        Y = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
        Y = (Y + Y.conj().T) / 2
        X = block_project(Y, es)
        if fixed: X = gauge_average(X, As)
        R = np.zeros((N, N), dtype=complex)
        for (val, P, m) in es:
            projs = minimal_projections_in(X, eigenbasis_of(P))
            ranks = [int(round(float(np.trace(q).real))) for q in projs]
            if sum(ranks) != m: return None
            sg = signs_to_zero(ranks)
            if sg is None: return None
            R = R + sum(s * q for s, q in zip(sg, projs))
        return R

    rows = []
    for fixed in (False, True):
        moved = 0; built = 0; mx = 0.0; worst_clause = 0.0
        for s in range(n_draw):
            R = draw(1000 + s if not fixed else 5000 + s, fixed)
            if R is None: continue
            built += 1
            v = verify_record(R, H, es, As)
            worst_clause = max(worst_clause, v["herm"], v["invol"], v["commH"], v["trace"])
            if v["nonconst"] < 1e-8: continue          # would violate clause (iii)
            mx = max(mx, v["commA"])
            if v["commA"] > 1e-6: moved += 1
        rows.append((("BLIND (control)" if not fixed else "TRANSPORT-FIXED"), built, moved,
                     mx, worst_clause))
    return name, N, rows, (classes, dims, table, cls_of, els, As, es, idx, mul, inv)

def central_check(els, mul, As, idx, out, name):
    _, e, inv, ix = conjugacy_classes(els, mul)
    cen = [g for g in els if all(mul(g, h) == mul(h, g) for h in els)]
    I = np.eye(As[0].shape[0])
    err = max(np.linalg.norm(As[idx[z]] - I) for z in cen)
    out.append("  %s: |Z(G)| = %d   max ||A_z - I|| over z in Z(G) = %.2e   -> the gauge rep factors "
               "through G/Z(G) (order %d)" % (name, len(cen), err, len(els) // len(cen)))

def main():
    out = []
    out.append("=" * 112)
    out.append("O-36 ADDENDUM: reconciling the exact verdict with C-43's 40/40, at MATCHED SAMPLE SIZE")
    out.append("=" * 112)
    out.append("  carrier      | family            | built/40 | MOVED by some A_h | max ||[A_h,R]|| | worst (i)-(iv) residual")
    out.append("-" * 112)
    keep = {}
    for gs in (group_dihedral(4), group_Q8(), group_dihedral(8), group_Zn(2), group_ZmxZn(2, 2)):
        els, mul, _, name = gs
        nm, N, rows, extra = records_sample(els, mul, name)
        keep[name] = (els, mul, extra)
        for i, (fam, built, moved, mx, wc) in enumerate(rows):
            out.append("  %-12s | %-17s | %5d/40 | %8d/40      | %-15.3e | %.1e"
                       % (("D(%s) dim %d" % (nm, N)) if i == 0 else "", fam, built, moved, mx, wc))
        out.append("-" * 112)
    out.append("")
    out.append("=" * 112)
    out.append("WHY: the gauge action is by CONJUGATION, so Z(G) acts trivially and the rep factors through G/Z(G)")
    out.append("-" * 112)
    for name, (els, mul, extra) in keep.items():
        classes, dims, table, cls_of, els2, As, es, idx, mul2, inv = extra
        central_check(els, mul, As, idx, out, "D(%s)" % name)
    out.append("")
    out.append("  D_4/Z = Z_2xZ_2 and Q_8/Z = Z_2xZ_2 are ABELIAN  -> every commutant block has d_i = 1")
    out.append("     -> the achievable traces on E are ALL integers of the parity of dim E, so 0 is reachable")
    out.append("        exactly when dim E is even, which is C-41's condition and nothing more.")
    out.append("  D_8/Z = D_4 is NON-ABELIAN -> 2-dimensional blocks DO appear (see step 2: 2x12 and 2x36),")
    out.append("     so the D_4 result is NOT an artifact of every block being 1-dimensional (D-17).")
    out.append("=" * 112)
    txt = "\n".join(out)
    print(txt)
    with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O36_EXACT/o36_reconcile.txt", "w") as f:
        f.write(txt + "\n")

if __name__ == "__main__":
    main()
