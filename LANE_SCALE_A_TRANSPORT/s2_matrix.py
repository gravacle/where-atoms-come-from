"""S2 -- BUILD THE CARRIERS AS ACTUAL MATRICES AND VERIFY EVERYTHING S1 PREDICTED.

Three jobs.
  (A) CONFIRM the exact census numerically: eigenvalue multiplicities of H, and the transport
      character Tr(A_h P_E) against the closed form.  A disagreement kills S1.
  (B) CONSTRUCT records and check clauses (i)-(iv) on the matrices, for TWO families in the
      SAME TABLE (D-15):
        GENERIC  -- a random trace-balanced involution per eigenspace.  This is exactly the
                    protocol that produced C-43's "40 of 40 moved".
        GAUGE    -- built inside the commutant of the whole transport group.  If C-43's
                    "40 of 40" were structural rather than sampling, this family would be
                    EMPTY on a non-abelian carrier.  It is the positive control that a
                    NON-ZERO answer is reachable at all.
      Every carrier carries an ABELIAN control of the SAME ORDER in the same table.
  (C) SCALE CONTROL (D-17): repeat on the 1x2 torus (2 vertices, 4 edges, 2 faces), the
      venue's own next size up, and see whether the reading survives.
"""
import sys, time, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_A_TRANSPORT")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
import glib
from carriers import (census, isotypic, subset_sums, fixed_record_exists, phi,
                      minimal_torus, torus_1x2, eigblocks, perm_apply, check_clauses,
                      generic_record, gauge_record, moved)
def say(*a): print(*a); sys.stdout.flush()
TOL = 1e-8

say("=" * 122)
say("S2   MATRIX VERIFICATION OF THE EXACT CENSUS, AND TWO RECORD FAMILIES IN THE SAME TABLE")
say("=" * 122)

NREC = 40                      # matches LANE_O34's sample size exactly, so the numbers compare
NREC_BIG = 12                  # dim 1024 carriers: fewer draws, same protocol
CARRIERS = [G for G in glib.ladder(32)]
def nrec(G): return NREC if G.n <= 16 else NREC_BIG

say("")
say("  (A) DOES THE MATRIX AGREE WITH THE EXACT CENSUS?")
say(f"  {'carrier':<12}{'dim':>6}{'mult (numeric)':>28}{'mult (exact)':>28}{'match':>7}{'chi_E(h) max err':>18}")
agree = True
store = {}
for G in CARRIERS:
    t0 = time.time()
    H, perms, D = minimal_torus(G)
    blocks = eigblocks(H)
    mult_num = sorted([Q.shape[1] for _, Q in blocks], reverse=True)
    ce = census(G)
    mult_ex = sorted([ce['dims'][v] for v in (-2, -1, 0) if ce['dims'][v] > 0], reverse=True)
    ok = mult_num == mult_ex
    # transport character check
    err = 0.0
    vals = sorted(ce['dims'].items(), key=lambda kv: kv[0])       # -2, -1, 0
    ev_sorted = sorted(blocks, key=lambda b: b[0])                # ascending eigenvalue
    for (v, dv), (ev, Q) in zip(vals, [b for b in ev_sorted if b[1].shape[1] > 0]):
        P = Q @ Q.conj().T
        for h in range(G.n):
            num = float(np.sum(P[perms[h], np.arange(D)]))
            err = max(err, abs(num - ce['chis'][v][h]))
    agree &= ok and err < 1e-6
    say(f"  {G.name:<12}{D:>6}{str(mult_num):>28}{str(mult_ex):>28}{str(ok):>7}{err:>18.3e}")
    store[G.name] = (G, H, perms, D, blocks, ce)
say(f"  ALL CARRIERS AGREE: {agree}")
if not agree:
    say("  SELF-CHECK FAILED -- NO CONCLUSION IS DRAWN FROM THIS LANE."); sys.exit(1)

say("")
say("  (B) TWO RECORD FAMILIES, SAME TABLE.  'moved' = some A_h with [A_h,R] != 0.")
say(f"  {'carrier':<12}{'|G|':>4}{'abel':>6}{'records':>9}{'GENERIC moved':>15}{'max||[A,R]||':>14}"
    f"{'GAUGE moved':>13}{'max||[A,R]||':>14}{'clauses(i)(iii)(iv)':>21}")
res = {}
for G in CARRIERS:
    _, H, perms, D, blocks, ce = store[G.name]
    rng = np.random.default_rng(12345)
    gm = 0; gmax = 0.0; am = 0; amax = 0.0; cl_ok = True; nbuilt = 0; nga = 0
    for s in range(nrec(G)):
        R = generic_record(blocks, rng, D)
        if R is None: continue
        nbuilt += 1
        i1, i3, i4 = check_clauses(R, blocks)
        cl_ok &= (i1 and i3 and i4)
        mv, mg = moved(R, perms); gm += mv; gmax = max(gmax, mg)
        Rg = gauge_record(blocks, perms, rng, D)
        if Rg is not None:
            nga += 1
            j1, j3, j4 = check_clauses(Rg, blocks)
            cl_ok &= (j1 and j3 and j4)
            mv2, mg2 = moved(Rg, perms); am += mv2; amax = max(amax, mg2)
    say(f"  {G.name:<12}{G.n:>4}{str(G.abelian):>6}{nbuilt:>9}{f'{gm}/{nbuilt}':>15}{gmax:>14.4f}"
        f"{f'{am}/{nga}':>13}{amax:>14.3e}{str(cl_ok):>21}")
    res[G.name] = dict(n=G.n, abelian=G.abelian, m=ce['dims'][-2], nbuilt=nbuilt,
                       gen_moved=gm, gen_max=gmax, gauge_moved=am, gauge_n=nga, gauge_max=amax,
                       clauses=cl_ok, phi=phi(isotypic(G, ce), ce['dims'])[2])
say("")
say("  GENERIC is C-43's protocol.  GAUGE is the same carrier, same clauses, records built inside")
say("  the transport commutant.  Both columns are on the SAME ROW so neither can be read alone.")

say("")
say("  (C) SCALE CONTROL (D-17) -- THE VENUE'S OWN NEXT SIZE UP: the 1x2 torus, V=2, E=4, F=2.")
say(f"  {'carrier':<12}{'dim':>7}{'multiplicities (1x2 torus)':>34}{'all even':>10}{'phi_fix(1x2)':>14}{'phi_fix(min)':>14}{'secs':>7}")
scale_rows = []
for G in CARRIERS:
    if G.n > 8: continue
    if G.n ** 4 > 5000: continue
    t0 = time.time()
    H2, perms2, D2 = torus_1x2(G)
    blocks2 = eigblocks(H2)
    mult2 = [Q.shape[1] for _, Q in blocks2]
    # phi at this carrier: decompose each eigenspace under the GLOBAL conjugation action
    cl, chi, d, cls_of = G.chars()
    fixd = 0; alld = 0
    okdec = True
    for ev, Q in blocks2:
        P = Q @ Q.conj().T
        cf = np.array([float(np.sum(P[perms2[c[0]], np.arange(D2)])) for c in cl])
        mult, dd = G.decompose(cf)
        if np.max(np.abs(mult - np.round(mult))) > 1e-5: okdec = False
        mult = np.round(mult).astype(int)
        if int(np.sum(dd * mult)) != Q.shape[1]: okdec = False
        fixd += int(np.sum(mult ** 2)); alld += Q.shape[1] ** 2
    ce = census(G); ph_min = phi(isotypic(G, ce), ce['dims'])[2]
    say(f"  {G.name:<12}{D2:>7}{str(mult2):>34}{str(all(m%2==0 for m in mult2)):>10}"
        f"{fixd/alld:>14.4f}{ph_min:>14.4f}{time.time()-t0:>7.1f}"
        + ("" if okdec else "   <-- DECOMPOSITION SELF-CHECK FAILED"))
    scale_rows.append((G.name, G.abelian, D2, fixd / alld, ph_min, okdec))
say("")
np.save("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_A_TRANSPORT/s2_res.npy",
        np.array([res, scale_rows], dtype=object), allow_pickle=True)
say("=" * 122); say("  READ -- from the numbers above"); say("=" * 122)
na = [r for r in res.values() if not r['abelian']]; ab = [r for r in res.values() if r['abelian']]
say(f"  GENERIC records moved: abelian carriers {sum(r['gen_moved'] for r in ab)}/{sum(r['nbuilt'] for r in ab)}"
    f"   non-abelian carriers {sum(r['gen_moved'] for r in na)}/{sum(r['nbuilt'] for r in na)}")
say(f"  GAUGE   records moved: abelian carriers {sum(r['gauge_moved'] for r in ab)}/{sum(r['gauge_n'] for r in ab)}"
    f"   non-abelian carriers {sum(r['gauge_moved'] for r in na)}/{sum(r['gauge_n'] for r in na)}")
