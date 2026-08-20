"""O-48-D PART 4.  CLAUSE (v).  IS A RECORD ON THE CHAIN *PROTECTED*?

Clause (v): NO ADMISSIBLE OPERATION ON A SINGLE CONTRACTIBLE REGION FLIPS THE RECORD.
Carrier data, as the model insists: a region is a CONTIGUOUS BLOCK of w sites (an arc, on a ring).

THE TEST IS EXACT AND COVERS EVERY UNITARY, NOT ONLY PAULIS.  Let U = U_B (x) I_out.
  * [U,H] = 0  <=>  <s|U|s'> = 0 unless E(s) = E(s').  Since U_B is one fixed matrix acting the
    same way for EVERY outside configuration o, (U_B)_{t,t'} != 0 requires E(t,o) = E(t',o) FOR ALL
    o.  Writing E(t,o) - E(t',o) as a function of o it is AFFINE in the outside spins, so it
    vanishes for all o iff
        (alpha) for every outside site, sum over its bonds into B of J * sigma_p is the same for
                t and t'      [when each outside site has one bond into B this just PINS that
                               boundary spin]
        (beta)  the INSIDE bond energy sum_{bonds inside B} J sigma sigma is the same for t and t'
  * U-dag Z_i U = -Z_i  <=>  (U_B)_{t,t'} != 0 requires sigma_i(t) = -sigma_i(t').
  * U UNITARY  =>  det U_B != 0  =>  by Frobenius-Konig its support CONTAINS A PERMUTATION.
    Conversely a permutation matrix obeying the two support conditions IS such a U.
  SO:  an admissible U on B flipping Z_i EXISTS  <=>  the allowed-support graph has a PERFECT
  MATCHING.  And the allowed graph is a disjoint union of COMPLETE BIPARTITE pieces -- t ~ t' iff
  they share the key ((alpha),(beta)) and disagree at i -- so a perfect matching exists IFF EVERY
  KEY CLASS CONTAINS EQUALLY MANY sigma_i = +1 AND sigma_i = -1 CONFIGURATIONS.
  That is an exact, O(2^w) decision procedure.  No search over unitaries, no tolerance.

CONTROLS IN THE SAME TABLE (D-15).  A carrier where a SMALL region MUST be able to flip:
the same chain with one coupling set to ZERO, which cuts it in two.  If the instrument reports
protection there too, it is broken.
"""
import sys, os, time, itertools
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from o48_common import pauli_matrix, symp, weight, all_paulis, PauliH


def say(*a):
    print(*a)
    sys.stdout.flush()


LINE = "=" * 118


def chain_bonds(n, J, ring=False):
    b = [(i, i + 1, J[i]) for i in range(n - 1)]
    if ring:
        b.append((n - 1, 0, J[n - 1]))
    return b


def region_flips(n, bonds, B, i):
    """EXACT: does some admissible unitary supported on the site set B flip Z_i?
       Returns (exists, n_key_classes, n_unbalanced_classes)."""
    if i not in B:
        return False, 0, 0
    Bl = sorted(B)
    pos = {s: k for k, s in enumerate(Bl)}
    w = len(Bl)
    inside = [(pos[a], pos[c], J) for a, c, J in bonds if a in pos and c in pos]
    # boundary bonds grouped by their OUTSIDE site
    outmap = {}
    for a, c, J in bonds:
        if (a in pos) != (c in pos):
            p, o = (a, c) if a in pos else (c, a)
            outmap.setdefault(o, []).append((pos[p], J))
    idx = np.arange(1 << w, dtype=np.int64)
    bits = ((idx[:, None] >> np.arange(w)[None, :]) & 1).astype(np.int64)
    sig = 1 - 2 * bits                                              # (2^w, w) of +-1
    keys = []
    for (p, q, J) in inside:                                        # (beta) each inside bond's
        keys.append(J * sig[:, p] * sig[:, q])                      # contribution, summed below
    Ein = sum(keys) if keys else np.zeros(1 << w, dtype=np.int64)
    cols = [Ein]
    for o, lst in sorted(outmap.items()):                           # (alpha) one column per
        cols.append(sum(J * sig[:, p] for p, J in lst))             # outside site
    K = np.stack(cols, axis=1)
    # group by key, then test balance in sigma_i
    order = np.lexsort(tuple(K[:, c] for c in range(K.shape[1])))
    Ks = K[order]
    si = sig[order, pos[i]]
    newk = np.ones(len(Ks), dtype=bool)
    newk[1:] = np.any(Ks[1:] != Ks[:-1], axis=1)
    lab = np.cumsum(newk) - 1
    nclass = int(lab[-1]) + 1
    plus = np.bincount(lab, weights=(si > 0).astype(np.float64), minlength=nclass)
    size = np.bincount(lab, minlength=nclass)
    unbal = int(np.sum(2 * plus != size))
    return (unbal == 0), nclass, unbal


def pauli_region_flips(n, terms, B, i):
    """CROSS-CHECK, Paulis only: is there an admissible Pauli supported in B flipping Z_i?"""
    ph = PauliH(n, terms)
    Ra = [0] * n
    Rb = [1 if k == i else 0 for k in range(n)]
    Bs = set(B)
    for a, b in all_paulis(n):
        if any((a[k] or b[k]) for k in range(n) if k not in Bs):
            continue
        if symp(a, b, Ra, Rb) != 1:
            continue
        if ph.admissible(a, b):
            return True, (list(a), list(b))
    return False, None


def blocks_of(n, w, ring):
    if ring:
        return [[(l + k) % n for k in range(w)] for l in range(n)] if w < n else [list(range(n))]
    return [list(range(l, l + w)) for l in range(0, n - w + 1)]


t0 = time.time()
say(LINE)
say("O-48-D  PART 4   CLAUSE (v):  IS THE RECORD PROTECTED ON THE CHAIN?")
say(LINE)

# ------------------------------------------------------------------ 4.0 instrument cross-check
say("")
say("4.0  THE EXACT CRITERION AGAINST AN EXHAUSTIVE PAULI SEARCH, on every block of every size,")
say("     n = 4..8.  The Pauli search can only ever find a SUBSET (Paulis are special unitaries),")
say("     so the exact test must be TRUE wherever the Pauli search succeeds, and it may be true in")
say("     places the Pauli search misses.  Any case where the Pauli search wins and the exact test")
say("     says NO would be a bug.")
say("")
say(f"  {'n':>3} {'J':>5} {'cases':>7} {'exact YES':>10} {'Pauli YES':>10} {'Pauli YES & exact NO':>22}")
for n in range(4, 9):
    for kind in ("UNI", "GEN"):
        J = [1] * (n - 1) if kind == "UNI" else [2 ** k for k in range(n - 1)]
        bonds = chain_bonds(n, J)
        terms = [([0] * n, [1 if s in (a, c) else 0 for s in range(n)], JJ) for a, c, JJ in bonds]
        cases = ex = pa = bad = 0
        for w in range(1, n + 1):
            for B in blocks_of(n, w, False):
                for i in B:
                    cases += 1
                    e, _, _ = region_flips(n, bonds, B, i)
                    p, _ = pauli_region_flips(n, terms, B, i)
                    ex += e
                    pa += p
                    if p and not e:
                        bad += 1
        say(f"  {n:>3} {kind:>5} {cases:>7} {ex:>10} {pa:>10} {bad:>22}")

# ------------------------------------------------------------------ 4.1 the sweep
say("")
say(LINE)
say("4.1  SWEEP OVER REGION SIZE w.  For each w, over EVERY contiguous block of that size and")
say("     EVERY record inside it: does an admissible operation on that one region flip the record?")
say(LINE)
say("")
say(f"  {'geometry':>14} {'J':>5} {'n':>3} " + " ".join(f"{('w=%d' % w):>5}" for w in range(1, 15))
    + f"  {'smallest w that flips':>22}")
summary = {}
for ring in (False, True):
    for kind in ("UNI", "GEN"):
        for n in range(3, 15):
            J = ([1] * n if kind == "UNI" else [2 ** k for k in range(n)])
            bonds = chain_bonds(n, J, ring=ring)
            cells, smallest = [], None
            for w in range(1, 15):
                if w > n:
                    cells.append("  .  ")
                    continue
                if ring and w == n:
                    cells.append(" WRAP")           # the whole ring is not a proper arc
                    continue
                hit = 0
                for B in blocks_of(n, w, ring):
                    for i in B:
                        e, _, _ = region_flips(n, bonds, B, i)
                        hit += int(e)
                cells.append(f"{hit:>5}")
                if hit and smallest is None:
                    smallest = w
            geom = "RING" if ring else "OPEN CHAIN"
            say(f"  {geom:>14} {kind:>5} {n:>3} " + " ".join(cells)
                + f"  {(str(smallest) if smallest else 'NONE up to w=n'):>22}")
            summary[(geom, kind, n)] = smallest
        say("")

say("  Cells count (block, record) pairs at that w for which an admissible operation on that")
say("  single region flips the record.  ' WRAP' marks w = n on a ring, where the region is the")
say("  whole ring and is NOT a proper arc -- it is not a contractible region at all.")

# ------------------------------------------------------------------ 4.2 controls
say("")
say(LINE)
say("4.2  CONTROLS (D-15).  Carriers where a SMALL region MUST be able to flip a record.")
say("     If the same instrument reports protection here, it is broken.")
say(LINE)
say("")
say(f"  {'control':>44} {'n':>3} {'smallest w that flips':>22} {'#(block,record) at that w':>26}")
CTRL = []
for n in (8, 10, 12):
    J = [2 ** k for k in range(n - 1)]
    Jc = list(J)
    Jc[n // 2] = 0                                              # CUT the chain in half
    CTRL.append((f"chain with J_{n // 2} = 0 (cut in two)", n, chain_bonds(n, Jc)))
for n in (8, 10):
    CTRL.append((f"NO Hamiltonian at all, H = 0", n, []))
for n in (8, 10):
    J = [2 ** k for k in range(n - 1)]
    b = chain_bonds(n, J)[: n // 2 - 1]                          # only the left half is coupled
    CTRL.append((f"only the left half coupled", n, b))
for tag, n, bonds in CTRL:
    smallest, cnt = None, 0
    for w in range(1, n + 1):
        hit = 0
        for B in blocks_of(n, w, False):
            for i in B:
                e, _, _ = region_flips(n, bonds, B, i)
                hit += int(e)
        if hit:
            smallest, cnt = w, hit
            break
    say(f"  {tag:>44} {n:>3} {(str(smallest) if smallest else 'NONE'):>22} {cnt:>26}")

# ------------------------------------------------------------------ 4.3 what the w=n operator is
say("")
say(LINE)
say("4.3  AT w = n ON THE OPEN CHAIN, WHAT IS THE OPERATOR?  Built explicitly and measured.")
say(LINE)
say("")
say(f"  {'n':>3} {'block':>10} {'operator':>18} {'||[U,H]||':>11} {'||U-dag Z_0 U + Z_0||':>22} "
    f"{'dE':>11} {'#records flipped':>18}")
for n in (4, 6, 8):
    J = [2 ** k for k in range(n - 1)]
    terms = [([0] * n, [1 if s in (i, i + 1) else 0 for s in range(n)], J[i]) for i in range(n - 1)]
    ph = PauliH(n, terms)
    H = ph.matrix()
    ok, W = pauli_region_flips(n, terms, list(range(n)), 0)
    U = pauli_matrix(*W)
    w_, V = np.linalg.eigh(H)
    k = int(np.sum(np.abs(w_ - w_[0]) < 1e-7))
    Q = V[:, :k]
    rho = Q @ Q.conj().T / k
    dE = float(np.real(np.trace((U @ rho @ U.conj().T) @ H)) - np.real(np.trace(rho @ H)))
    Z0 = pauli_matrix([0] * n, [1] + [0] * (n - 1))
    nf = sum(1 for s in range(n)
             if np.linalg.norm(U.conj().T @ pauli_matrix([0] * n, [1 if t == s else 0 for t in range(n)]) @ U
                               + pauli_matrix([0] * n, [1 if t == s else 0 for t in range(n)])) < 1e-8)
    name = "".join({(0, 0): "I", (1, 0): "X", (0, 1): "Z", (1, 1): "Y"}[(W[0][s], W[1][s])] for s in range(n))
    say(f"  {n:>3} {'all n sites':>10} {name:>18} {np.linalg.norm(U @ H - H @ U):>11.1e} "
        f"{np.linalg.norm(U.conj().T @ Z0 @ U + Z0):>22.1e} {dE:>+11.6f} {nf:>18}")

say("")
say(LINE)
say("  READ -- PART 4, FILLED IN FROM THE NUMBERS ABOVE, NOT IN ADVANCE")
say(LINE)
say("")
say("  1. THE INSTRUMENT AGREES WITH AN EXHAUSTIVE PAULI SEARCH.  Section 4.0: across 4..8 sites,")
say("     both coupling sets, every block of every size and every record inside it -- 315 cases per coupling set, 630 in all --")
say("     the exact criterion and the exhaustive Pauli search return the SAME count, and the")
say("     'Pauli says yes, exact says no' column is 0 everywhere.  Note the counts are also equal,")
say("     so on this carrier the extra freedom of non-Pauli unitaries buys nothing.")
say("")
say("  2. ON THE OPEN CHAIN THE RECORD IS PROTECTED AGAINST EVERY REGION SMALLER THAN THE WHOLE")
say("     CARRIER, AND AGAINST NOTHING ELSE.  Section 4.1: for n = 3..14 and both coupling sets,")
say("     every cell at w < n is 0 and the cell at w = n is n.  The smallest region size whose")
say("     admissible operations can flip a record is EXACTLY n, at every n.  This covers ALL")
say("     unitaries supported on the region, not only Paulis.")
say("")
say("  3. SO CLAUSE (v) FAILS ON THE OPEN CHAIN.  In one dimension a contiguous block of sites is")
say("     contractible, and the whole chain IS a contiguous block.  Section 4.3 exhibits the")
say("     operator: X^(x)n, ||[U,H]|| = 0, U-dag Z_0 U = -Z_0 exactly, dE = +0.000000, and it flips")
say("     all n records at once.  A single contractible region -- the entire carrier -- writes the")
say("     record for free.  THEREFORE THE Z_i ON THE OPEN CHAIN ARE NOT RECORDS IN THE PROGRAM'S")
say("     FULL FIVE-CLAUSE SENSE, and O-47, and Parts 1-3 above, are scoped to clauses (i)-(iv).")
say("")
say("  4. ON THE RING, UNDER THE ARC CONVENTION, CLAUSE (v) HOLDS.  If a contractible region is a")
say("     PROPER ARC -- the reading that matches T-11's 'a forest, i.e. no cycle' -- then no")
say("     contractible region flips a record at any w up to n-1, for n = 3..14 and both coupling")
say("     sets, and the flipping operator at w = n is the whole ring, which wraps and is not a")
say("     contractible region at all.  The carrier's answer to clause (v) therefore turns entirely")
say("     on a convention about what 'contractible' means on a 1D lattice.  That is CARRIER DATA,")
say("     exactly as the model warns, and it is a decision, not a measurement.")
say("")
say("  5. THE CONTROLS REGISTER.  Section 4.2: cutting the chain by setting one coupling to zero")
say("     drops the smallest flipping region from n to the size of the SMALLER COMPONENT -- 3 at")
say("     n = 8, 4 at n = 10, 5 at n = 12.  With H = 0 it drops to 1.  With only the left half")
say("     coupled it is 1, on the uncoupled sites.  So the zeros in 4.1 are measurements.")
say("")
say("  6. WHAT PROTECTION ACTUALLY COSTS HERE.  The protection distance equals the writer weight,")
say("     n, which Part 1 measured independently.  They are the same number for the same reason:")
say("     the only admissible operators that flip the record have full support.  Protection and")
say("     non-locality of the writer are not two properties -- they are one.")
say("")
say(f"  runtime {time.time() - t0:.1f}s")
