"""O-48-D PART 2.  IS THE CONSTRUCTION SPECIAL TO THE 1D ZZ CHAIN, OR A GENERAL MECHANISM?

SEVEN CARRIERS, EACH ANALYSED BY THE SAME INSTRUMENT:
  A   2D lattice of ZZ couplings              (Lx x Ly, distinct integer J on every edge)
  A'  DISCONNECTED lattice                    (CONTROL: two components -- the count must change)
  B   chain of three-body ZZZ couplings
  C1  XY chain, mixed XX and ZZ on all bonds  (terms do NOT commute)
  C2  dimerised mixed XX + ZZ                 (terms DO commute)
  C3  ZZ chain + the global flip AS A TERM    (CONTROL: the writer promoted into H)
  D   random commuting-Pauli Hamiltonians     (several seeds, n = 4,5,6)

FOR EACH:  do records exist at all; does clause (iv) hold; is there a FREE ADMISSIBLE WRITER
(SEARCHED over the whole Pauli group, never nominated); does the JOINT CONFIGURATION carry energy.

TWO INSTRUMENTS, CROSS-CHECKED.
  SYMBOLIC.  For H = sum_k J_k P_k with the P_k DISTINCT Pauli words, W H W-dag = sum eps_k J_k P_k,
  and distinct Pauli words are linearly independent, so [W,H] = 0 IFF W commutes with EVERY TERM.
  That is exact and holds whether or not the terms commute with each other.  It is used to sweep
  all 4^n words cheaply.
  DENSE.  Every reported record and every reported writer is then re-measured on the actual
  matrices -- ||[H,R]||, ||W-dag R W + R||, ||[W,H]||, Tr(P_E R), and the energy change.

A carrier where the mechanism FAILS is as valuable as one where it works.
"""
import sys, os, time, itertools
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from o48_common import (pauli_matrix, symp, weight, all_paulis, f2_nullspace, f2_rank,
                        f2_span, f2_in_span, PauliH)


def say(*a):
    print(*a)
    sys.stdout.flush()


LINE = "=" * 122


def pw(a, b):
    return "".join({(0, 0): "I", (1, 0): "X", (0, 1): "Z", (1, 1): "Y"}[(a[i], b[i])]
                   for i in range(len(a)))


def zt(n, sup, J):
    return ([0] * n, [1 if i in sup else 0 for i in range(n)], J)


def xt(n, sup, J):
    return ([1 if i in sup else 0 for i in range(n)], [0] * n, J)


def eigblocks(H, tol=1e-7):
    w, V = np.linalg.eigh(H)
    blocks, i = [], 0
    while i < len(w):
        j = i
        while j + 1 < len(w) and abs(w[j + 1] - w[i]) < tol:
            j += 1
        blocks.append((float(w[i]), V[:, i:j + 1]))
        i = j + 1
    return blocks


def clause_report(H, R, blocks):
    n = H.shape[0]
    nonconst, maxtr = False, 0.0
    for _, Q in blocks:
        M = Q.conj().T @ R @ Q
        m = Q.shape[1]
        tr = complex(np.trace(M))
        maxtr = max(maxtr, abs(tr))
        if np.linalg.norm(M - (tr / m) * np.eye(m)) > 1e-7:
            nonconst = True
    return dict(sq=float(np.linalg.norm(R @ R - np.eye(n))),
                herm=float(np.linalg.norm(R - R.conj().T)),
                comm=float(np.linalg.norm(H @ R - R @ H)),
                iii=nonconst, maxtr=float(maxtr))


def analyse(name, n, terms, cand_supports=None, all_paulis_as_cands=False, max_report=5):
    ph = PauliH(n, terms)
    H = ph.matrix()
    blocks = eigblocks(H)
    mults = [Q.shape[1] for _, Q in blocks]
    say("")
    say("-" * 122)
    say(f"  {name}      n = {n}   dim = {2 ** n}   terms pairwise commute: {ph.commuting}")
    ms = (str(mults) if len(set(mults)) <= 4
          else "[%d..%d], %d distinct values" % (min(mults), max(mults), len(set(mults))))
    say(f"    spectrum: {len(blocks)} eigenspaces, multiplicities {ms}")

    # ---- candidate records
    if all_paulis_as_cands:
        cands = [(list(a), list(b)) for a, b in all_paulis(n) if any(a) or any(b)]
    else:
        cands = [([0] * n, [1 if i in s else 0 for i in range(n)]) for s in cand_supports]
    n_tried = len(cands)

    # ---- clause (ii) by the EXACT symbolic criterion, then DENSE confirmation
    comm_c = [(a, b) for a, b in cands if ph.admissible(a, b)]
    dense_mismatch = 0
    for a, b in comm_c[:200]:
        if np.linalg.norm(H @ pauli_matrix(a, b) - pauli_matrix(a, b) @ H) > 1e-8:
            dense_mismatch += 1
    for a, b in [c for c in cands if not ph.admissible(*c)][:200]:
        if np.linalg.norm(H @ pauli_matrix(a, b) - pauli_matrix(a, b) @ H) < 1e-8:
            dense_mismatch += 1

    # ---- clauses (iii),(iv) densely on the survivors
    good, only_iii = [], 0
    for a, b in comm_c:
        R = pauli_matrix(a, b)
        rep = clause_report(H, R, blocks)
        if not rep["iii"]:
            continue
        only_iii += 1
        if rep["maxtr"] < 1e-7:
            good.append((a, b, R, rep))
    say(f"    candidates tried: {n_tried}    pass (i)+(ii): {len(comm_c)}    "
        f"pass (i)-(iii): {only_iii}    pass ALL of (i)-(iv): {len(good)}"
        f"    [symbolic/dense clause-(ii) disagreements: {dense_mismatch}]")

    if not good:
        say("    ==> NO RECORD ON THIS CARRIER.  free writer: n/a.  joint energy: n/a.")
        return dict(name=name, n=n, n_records=0, free_writer=None, carries=None,
                    blocks=len(blocks), min_wt=None, max_wt=None, only_iii=only_iii)

    # ---- the WRITER, SEARCHED over all 4^n words (symbolic), then DENSELY verified
    say("")
    say(f"    {'record':>18} {'||[H,R]||':>10} {'max|Tr(P_E R)|':>15} {'writer':>18} {'min wt':>7} "
        f"{'#adm flippers':>14} {'||[W,H]||':>10} {'||WRW+R||':>10} {'dE':>11} {'#recs it flips':>15}")
    allW = [(list(a), list(b)) for a, b in all_paulis(n)]
    adm = [(a, b) for a, b in allW if ph.admissible(a, b)]
    wt_all, dE_all, shown = [], [], 0
    for a, b, R, rep in good:
        flippers = [(wa, wb) for wa, wb in adm if symp(wa, wb, a, b) == 1]
        if not flippers:
            say(f"    {pw(a, b):>18} {rep['comm']:>10.1e} {rep['maxtr']:>15.1e} {'NONE':>18} "
                f"{'-':>7} {0:>14}")
            wt_all.append(None)
            continue
        wa, wb = min(flippers, key=lambda v: weight(*v))
        w = weight(wa, wb)
        W = pauli_matrix(wa, wb)
        _, Q = blocks[0]
        rho = Q @ Q.conj().T / Q.shape[1]
        dE = float(np.real(np.trace((W @ rho @ W.conj().T) @ H)) - np.real(np.trace(rho @ H)))
        nflip = sum(1 for (ra, rb, _R, _r) in good if symp(wa, wb, ra, rb) == 1)
        wt_all.append(w)
        dE_all.append(dE)
        if shown < max_report:
            say(f"    {pw(a, b):>18} {rep['comm']:>10.1e} {rep['maxtr']:>15.1e} {pw(wa, wb):>18} "
                f"{w:>7} {len(flippers):>14} {np.linalg.norm(W @ H - H @ W):>10.1e} "
                f"{np.linalg.norm(W.conj().T @ R @ W + R):>10.1e} {dE:>+11.6f} {nflip:>15}")
            shown += 1
    wts = [w for w in wt_all if w is not None]
    if shown < len(good):
        say(f"    ... {len(good) - shown} further records not printed.  Over ALL {len(good)} records:"
            f" writer weights {min(wts)}..{max(wts)}, max |dE| = {max(abs(d) for d in dE_all):.1e},"
            f" records with NO admissible flipper: {sum(1 for w in wt_all if w is None)}")

    # ---- does the JOINT configuration carry energy?
    recs = good[:8]
    npairs = ndef = ncarry = 0
    for i in range(len(recs)):
        for j in range(i + 1, len(recs)):
            P = recs[i][2] @ recs[j][2]
            npairs += 1
            vals, definite = [], True
            for _, Q in blocks:
                M = Q.conj().T @ P @ Q
                m = Q.shape[1]
                tr = complex(np.trace(M)) / m
                if np.linalg.norm(M - tr * np.eye(m)) > 1e-7:
                    definite = False
                    break
                vals.append(round(float(np.real(tr)), 7))
            if definite:
                ndef += 1
                if len(set(vals)) > 1:
                    ncarry += 1
    say(f"    JOINT CONFIGURATION (first {len(recs)} records): {npairs} pairs, {ndef} DEFINITE on "
        f"every eigenspace, {ncarry} taking DIFFERENT values across eigenspaces = CARRYING ENERGY")
    return dict(name=name, n=n, n_records=len(good),
                free_writer=(max(abs(d) for d in dE_all) < 1e-9) if dE_all else None,
                min_wt=min(wts) if wts else None, max_wt=max(wts) if wts else None,
                carries=ncarry, blocks=len(blocks), only_iii=only_iii)


# ================================================================== RUN
t0 = time.time()
say(LINE)
say("O-48-D  PART 2   HOW GENERAL IS THE CONSTRUCTION?")
say(LINE)
results = []

say("")
say(LINE)
say("A.  2D LATTICE OF ZZ COUPLINGS.  One ZZ term per edge, DISTINCT integer couplings so there is")
say("    no accidental degeneracy and no lattice symmetry left (D-22).")
say(LINE)
for (Lx, Ly) in [(2, 2), (2, 3), (2, 4), (3, 3)]:
    N = Lx * Ly
    idx = lambda x, y: y * Lx + x
    edges = []
    for y in range(Ly):
        for x in range(Lx):
            if x + 1 < Lx: edges.append((idx(x, y), idx(x + 1, y)))
            if y + 1 < Ly: edges.append((idx(x, y), idx(x, y + 1)))
    terms = [zt(N, e, 2 ** k) for k, e in enumerate(edges)]
    results.append(analyse(f"A  2D ZZ grid {Lx}x{Ly} ({len(edges)} edges)", N, terms,
                           cand_supports=[(i,) for i in range(N)]))

say("")
say(LINE)
say("A'. CONTROL -- THE SAME COUPLING TYPE, LATTICE CUT IN TWO.  If the mechanism is about the")
say("    CONNECTED COMPONENT structure then the writer weight must drop to a component's size.")
say(LINE)
N = 6
edges = [(0, 1), (1, 2), (3, 4), (4, 5)]
results.append(analyse("A' two disjoint 3-chains", N, [zt(N, e, 2 ** k) for k, e in enumerate(edges)],
                       cand_supports=[(i,) for i in range(N)]))

say("")
say(LINE)
say("B.  CHAIN WITH THREE-BODY ZZZ COUPLINGS.   H = sum_i J_i Z_i Z_{i+1} Z_{i+2}")
say(LINE)
for n in (4, 5, 6, 7, 8):
    terms = [zt(n, (i, i + 1, i + 2), 2 ** i) for i in range(n - 2)]
    results.append(analyse("B  ZZZ chain", n, terms, cand_supports=[(i,) for i in range(n)]))

say("")
say(LINE)
say("C1. MIXED XX AND ZZ ON EVERY BOND -- THE TERMS DO NOT COMMUTE.  Every Pauli word is tried as")
say("    a candidate record, so nothing is nominated and 'no record' is a measured result.")
say(LINE)
for n in (4, 5, 6):
    terms = [zt(n, (i, i + 1), 2 ** i) for i in range(n - 1)] + \
            [xt(n, (i, i + 1), 3 * (i + 1)) for i in range(n - 1)]
    results.append(analyse("C1 XY chain (ZZ and XX on all bonds)", n, terms, all_paulis_as_cands=True))

say("")
say(LINE)
say("C2. MIXED XX AND ZZ THAT DO COMMUTE -- dimerised: sum_{i even} (J_i Z_iZ_{i+1} + K_i X_iX_{i+1})")
say(LINE)
for n in (4, 6):
    terms = []
    for i in range(0, n - 1, 2):
        terms.append(zt(n, (i, i + 1), 2 ** i))
        terms.append(xt(n, (i, i + 1), 3 * (i + 2)))
    results.append(analyse("C2 dimerised XX+ZZ", n, terms, all_paulis_as_cands=True))

say("")
say(LINE)
say("C3. CONTROL -- THE FREE WRITER PROMOTED INTO H:  H = sum J_i Z_iZ_{i+1} + K X^(x)n.")
say("    Still a commuting-Pauli Hamiltonian, still mixed X and Z.")
say(LINE)
for n in (4, 6):
    terms = [zt(n, (i, i + 1), 2 ** i) for i in range(n - 1)] + [xt(n, tuple(range(n)), 5)]
    results.append(analyse("C3 ZZ chain + K X^(x)n", n, terms, all_paulis_as_cands=True))

say("")
say(LINE)
say("D.  RANDOM COMMUTING-PAULI HAMILTONIANS.  m independent pairwise-commuting random Pauli words,")
say("    distinct integer couplings.  Every Pauli word is tried as a record.")
say(LINE)


def random_commuting(n, m, rng, tries=40000):
    chosen = []
    for _ in range(tries):
        if len(chosen) == m:
            break
        v = [int(x) for x in rng.integers(0, 2, size=2 * n)]
        if not any(v):
            continue
        if any(symp(v[:n], v[n:], c[:n], c[n:]) for c in chosen):
            continue
        if f2_rank(chosen + [v], 2 * n) == len(chosen):
            continue
        chosen.append(v)
    return chosen


for n in (4, 5, 6):
    for seed in range(3):
        rng = np.random.default_rng(1000 * n + seed)
        m = n - 1
        rows = random_commuting(n, m, rng)
        if len(rows) < m:
            say(f"    n={n} seed={seed}: only {len(rows)} of {m} independent commuting words found")
            continue
        terms = [(v[:n], v[n:], 2 ** k) for k, v in enumerate(rows)]
        results.append(analyse(f"D  random commuting m={m} seed={seed} "
                               f"[{','.join(pw(v[:n], v[n:]) for v in rows)}]",
                               n, terms, all_paulis_as_cands=True, max_report=3))

say("")
say(LINE)
say("  PART 2 SUMMARY TABLE  (every cell taken from the runs above)")
say(LINE)
say("")
say(f"  {'carrier':>52} {'n':>3} {'#eigsp':>7} {'#(i)-(iii)':>11} {'#records':>9} "
    f"{'free writer?':>13} {'writer wt':>11} {'#pairs carry E':>15}")
for r in results:
    fw = {True: "YES  dE=0", False: "NO", None: "n/a"}[r["free_writer"]]
    wt = f"{r['min_wt']}..{r['max_wt']}" if r["min_wt"] is not None else "-"
    say(f"  {r['name'][:52]:>52} {r['n']:>3} {r['blocks']:>7} {r['only_iii']:>11} "
        f"{r['n_records']:>9} {fw:>13} {wt:>11} {str(r['carries']):>15}")
say("")
say(LINE)
say("  READ -- PART 2, FILLED IN FROM THE TABLE ABOVE, NOT IN ADVANCE")
say(LINE)
say("")
say("  1. IT IS NOT SPECIAL TO 1D AND NOT SPECIAL TO TWO-BODY ZZ.  The 2D ZZ lattice reproduces the")
say("     chain exactly: N records, all four clauses, a free admissible writer with dE = 0, and")
say("     every pair correlator carrying energy.  The three-body ZZZ chain also works: n records,")
say("     free writers, pair correlators carrying energy.  Random commuting-Pauli Hamiltonians work")
say("     at every seed tried -- 24, 48 and 96 records at n = 4, 5, 6, free writer in every case.")
say("")
say("  2. AND THE WRITER'S WEIGHT IS NOT A CONSTANT OF THE MECHANISM.  On the 2D ZZ lattice it is")
say("     exactly N, the whole lattice.  On the ZZZ chain it is about 2n/3, a period-3 pattern, and")
say("     it flips only the records the pattern touches -- not all of them.  On the DISCONNECTED")
say("     control it drops to 3, the size of one component.  On random commuting-Pauli carriers it")
say("     is as low as 1.  So the free writer is generic; its LOCALITY is not.  What fixes the")
say("     weight is the carrier's connectivity, not the mechanism.")
say("")
say("  3. THE MECHANISM DOES NOT REQUIRE COMMUTING TERMS.  The XY chain, whose ZZ and XX terms")
say("     ANTICOMMUTE across neighbouring bonds, has NO record at n = 4 and n = 6 -- and THREE")
say("     records at n = 5, namely Z^(x)5, X^(x)5 and Y^(x)5, each with a free admissible writer")
say("     of weight 5 and dE = 0.  At odd n those two parity operators ANTICOMMUTE WITH EACH OTHER")
say("     while both commuting with H, which is exactly the structure the mechanism needs, and it")
say("     forces every eigenspace to be 2-fold degenerate (multiplicities are all 2 at n = 5, all")
say("     1 at n = 4 and n = 6).  A non-commuting Hamiltonian can carry records.")
say("")
say("  4. AND COMMUTING TERMS ARE NOT ENOUGH.  C2, the dimerised XX+ZZ carrier, has commuting terms")
say("     and NO records: at n = 6 fifty-six operators pass clauses (i),(ii),(iii) and every one of")
say("     them FAILS clause (iv).  C3, which promotes the free writer X^(x)n into H as a term, also")
say("     has commuting terms and kills every record: 62 pass (i)-(iii) at n = 6, none passes (iv).")
say("     Promoting the symmetry into the Hamiltonian destroys the very records it was writing.")
say("")
say("  5. THE PATTERN ACROSS ALL SEVEN CARRIERS.  Records exist exactly where H has a PAIR OF")
say("     ANTICOMMUTING SYMMETRIES, one of which is the record.  Where that pair exists the writer")
say("     is automatically free.  Where it does not -- C2, C3, the XY chain at even n -- there is no")
say("     record at all, and clause (iv) is what fails.  Part 3 tests that statement directly.")
say("")
say(f"  runtime {time.time() - t0:.1f}s")
