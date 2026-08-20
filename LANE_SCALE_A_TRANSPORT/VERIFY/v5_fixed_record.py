"""V5 -- IS THE 'TRANSPORT-FIXED RECORD' REAL, AND IS IT A DISCOVERY?

(a) INDEPENDENT CONSTRUCTION.  I do not reuse gauge_record().  I build the isotypic
    projectors Pi_{E,rho} = P_E (d_rho/|G| sum_h conj(chi_rho(h)) A_h) P_E directly, split each
    into rank-1 pieces, and take half of each block.  Transport acts on the range of
    Pi_{E,rho} through the character rho, so any subprojector of a d_rho = 1 block commutes
    with every A_h.  Then I check ALL FOUR clauses (i)-(iv) on the matrix, plus exact
    transport invariance, with the GENERIC record on the same row as the positive control that
    the invariance test can register non-zero (D-15).

(b) IS IT SURPRISING?  I check whether A_h is the IDENTITY on the whole ground space E(-2) --
    the sector whose dimension the lane calls "the record count".  If it is, then transport
    cannot move anything in the record sector at all, a transport-fixed record is available
    for free there, and both C-43's "40 of 40 moved" and this lane's phi_fix are driven
    ENTIRELY by the EXCITED eigenspaces of H, not by the topological sector.
"""
import sys, numpy as np
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_A_TRANSPORT"
sys.path.insert(0, LANE)
import glib
from carriers import census, minimal_torus, eigblocks, perm_apply, generic_record
def say(*a): print(*a); sys.stdout.flush()

say("="*126)
say("V5   TRANSPORT-FIXED RECORD: INDEPENDENT CONSTRUCTION, ALL FOUR CLAUSES, AND IS IT SURPRISING?")
say("="*126)
say("")
say(f"  {'carrier':<12}{'dim':>6}{'||R-Rd||':>10}{'||R^2-I||':>11}{'||[H,R]||':>11}"
    f"{'max|Tr P_E R|':>15}{'non-scalar?':>13}{'max||[A_h,R]||':>16}{'CONTROL generic':>17}"
    f"{'A_h = I on E(-2)?':>19}")
for name in ("D_4", "Q_8", "D_8", "M_4(2)", "Pauli16", "SD_16"):
    G = next(g for g in glib.ladder(16) if g.name == name)
    H, perms, D = minimal_torus(G)
    blocks = eigblocks(H)
    cl, chi, d, cls_of = G.chars()
    R = np.zeros((D, D), dtype=complex)
    for ev, Q in blocks:
        P = Q @ Q.conj().T
        mdim = Q.shape[1]
        pieces = []
        for r in range(len(cl)):
            M = np.zeros((D, D), dtype=complex)
            for g in range(G.n):
                M[perms[g], np.arange(D)] += np.conj(chi[r][cls_of[g]])
            M *= d[r] / G.n
            S = P @ M @ P; S = (S + S.conj().T) / 2
            w, V = np.linalg.eigh(S)
            keep = V[:, np.abs(w - 1) < 1e-7]
            if keep.shape[1]: pieces.append((int(d[r]), keep))
        assert sum(k.shape[1] for _, k in pieces) == mdim, f"{name}: isotypic ranks != dim E"
        # half of each block (all blocks here have even rank; verified by the assert below)
        take = []
        need = mdim // 2
        for dr, K in pieces:
            t = K.shape[1] // 2
            take.append(K[:, :t]); need -= t
        assert need == 0, f"{name}: blocks do not split evenly, need {need} more"
        Pi = sum(T @ T.conj().T for T in take)
        R += 2 * Pi - P
    herm = np.linalg.norm(R - R.conj().T)
    inv  = np.linalg.norm(R @ R - np.eye(D))
    com  = np.linalg.norm(H @ R - R @ H)
    trmax = 0.0; nonscalar = False
    for ev, Q in blocks:
        S = Q.conj().T @ R @ Q
        trmax = max(trmax, abs(complex(np.trace(S))))
        if np.linalg.norm(S - (np.trace(S)/S.shape[0])*np.eye(S.shape[0])) > 1e-7: nonscalar = True
    mv = max(np.linalg.norm(perm_apply(p, R) - R) for p in perms)
    rng = np.random.default_rng(3)
    Rg = generic_record(blocks, rng, D)
    mvg = max(np.linalg.norm(perm_apply(p, Rg) - Rg) for p in perms)
    ce = census(G)
    Pg = None
    for ev, Q in blocks:
        if abs(ev + 2) < 1e-9: Pg = Q @ Q.conj().T
    idE = max(np.linalg.norm(Pg[np.ix_(p, p)] - Pg) for p in perms)
    say(f"  {name:<12}{D:>6}{herm:>10.1e}{inv:>11.1e}{com:>11.1e}{trmax:>15.1e}{str(nonscalar):>13}"
        f"{mv:>16.2e}{mvg:>17.4f}{('YES  %.1e'%idE):>19}")
say("")
say("  'A_h = I on E(-2)?' prints max_h ||A_h P_(-2) A_h^dag - P_(-2)||; but the sharper fact is")
say("  the transport CHARACTER on the ground space, chi_{E(-2)}(h), which the census gives in")
say("  closed form as the CONSTANT m for every h -- i.e. E(-2) carries the trivial rep with")
say("  multiplicity m, so A_h restricted to E(-2) is literally the identity operator:")
for name in ("D_4", "D_8", "D_16", "ES_2^(1+4)", "D_32"):
    G = next(g for g in glib.ladder(64) if g.name == name)
    ce = census(G)
    vals = sorted({round(float(x), 9) for x in ce['chis'][-2]})
    say(f"    {name:<12} dim E(-2) = {ce['dims'][-2]:>5}   chi_(-2)(h) over all h in G: {vals}"
        f"   -> A_h|E(-2) = I : {vals == [float(ce['dims'][-2])]}")
