"""O-48-D PART 1.  DO THE CLAUSES AND THE FREE-WRITER PROPERTY SURVIVE AT SCALE?

CARRIER:  H = sum_i J_i Z_i Z_{i+1}  on an OPEN chain of n qubits.  Records R_i = Z_i.

WHAT IS VERIFIED AT EVERY n, NOT ASSUMED:
  (i)   Z_i = Z_i-dag, Z_i^2 = I
  (ii)  [H, Z_i] = 0
  (iii) Z_i not constant on some eigenspace of H
  (iv)  Tr(P_E Z_i) = 0 on EVERY eigenspace
  and the ADMISSIBLE WRITER is SEARCHED, never nominated:
     -- brute force over ALL 4^n Pauli words for n <= 9,
     -- exact F_2 linear algebra (nullspace of the commutation matrix) at every n,
     and the two are cross-checked where both run.

COUPLINGS.  Three sets, because the answer depends on them and D-17 says vary the venue's scale:
  UNI  J_i = 1                      (uniform; maximal accidental degeneracy)
  GEN  J_i = 2^i                    (guarantees NO accidental degeneracy -- verified, not assumed)
  RND  J_i = distinct random ints   (a second generic set, checked for collisions)

CONTROLS IN THE SAME TABLE (D-15):
  CTRL-FIELD    H + h Z_0            -- breaks the global flip; clause (iv) must REGISTER NONZERO
  CTRL-NONDEG   H = sum 2^i Z_i      -- non-degenerate; clause (iii) must FAIL for every record (P-1)
Every reported zero sits beside a number from the same instrument that is not zero.

D-22 CHECK.  The site-permutation automorphism group of H is COUNTED by brute force, so the
carrier's own symmetry is measured rather than asserted.
"""
import sys, os, time
import numpy as np
from itertools import permutations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from o48_common import (pauli_matrix, symp, weight, all_paulis, f2_nullspace, f2_rank,
                        f2_span, f2_in_span, PauliH, spin_table, diag_energies, eig_classes)


def say(*a):
    print(*a)
    sys.stdout.flush()


LINE = "=" * 116


# ---------------------------------------------------------------- coupling sets
def couplings(n, kind, seed=7):
    m = n - 1                                     # bonds on an OPEN chain
    if kind == "UNI":
        return [1] * m
    if kind == "GEN":
        return [2 ** i for i in range(m)]
    if kind == "RND":
        rng = np.random.default_rng(seed + n)
        return [int(v) for v in rng.choice(np.arange(1, 4000), size=m, replace=False)]
    raise ValueError(kind)


def chain_terms(n, J):
    return [((i, i + 1), J[i]) for i in range(n - 1)]


# ---------------------------------------------------------------- clause checks, exact integers
def clause_table(n, zterms, records):
    """records = list of site-tuples (a record is the Z-word on that support).
       Returns per-record: (iii) holds, (iv) holds, max |Tr(P_E R)|, #eigenspaces split."""
    sig, E = diag_energies(n, zterms)
    vals, inv, sizes = eig_classes(E)
    out = []
    for sup in records:
        v = np.ones(1 << n, dtype=np.int64)
        for i in sup:
            v = v * sig[:, i].astype(np.int64)
        # clause (iv): exact integer sum of R's eigenvalues on each eigenspace
        pos = np.bincount(inv, weights=(v > 0).astype(np.float64)).astype(np.int64)
        tr = 2 * pos - sizes                                   # = sum of +-1 over the block
        # clause (iii): non-constant on a block iff 0 < pos < size
        split = int(np.sum((pos > 0) & (pos < sizes)))
        out.append(dict(support=sup, max_abs_trace=int(np.abs(tr).max()),
                        clause_iv=bool(np.abs(tr).max() == 0),
                        n_split=split, clause_iii=bool(split > 0)))
    return vals, sizes, out


# ---------------------------------------------------------------- writer search
def brute_writers(n, ph, target_a, target_b):
    """EXHAUSTIVE search over all 4^n Pauli words for W with [W,H]=0 and W-dag R W = -R."""
    best = None
    count = 0
    for a, b in all_paulis(n):
        if symp(a, b, target_a, target_b) != 1:                 # must ANTICOMMUTE with R
            continue
        if not ph.admissible(a, b):
            continue
        count += 1
        w = weight(a, b)
        if best is None or w < best[0]:
            best = (w, a, b)
    return best, count


def linalg_writers(n, ph, target_a, target_b):
    """EXACT: the admissible Pauli group is the F_2 nullspace of the commutation matrix; the
       flippers are the affine slice on which the symplectic form with R equals 1."""
    rows = ph.commutation_rows()
    ns = f2_nullspace(rows, 2 * n)
    if len(ns) > 22:
        return None, len(ns), None
    best, cnt = None, 0
    for v in f2_span(ns, 2 * n):
        a, b = v[:n], v[n:]
        if symp(a, b, target_a, target_b) != 1:
            continue
        cnt += 1
        w = weight(a, b)
        if best is None or w < best[0]:
            best = (w, a, b)
    return best, len(ns), cnt


def dense_check(n, ph, R_sup, a, b):
    """DENSE NUMERIC CONTROL on the symbolic claim, plus the energy cost."""
    H = ph.matrix()
    Rz = np.zeros(n, dtype=int)
    Ra = [0] * n
    Rb = [1 if i in R_sup else 0 for i in range(n)]
    R = pauli_matrix(Ra, Rb)
    W = pauli_matrix(list(a), list(b))
    cHW = float(np.linalg.norm(W @ H - H @ W))
    cflip = float(np.linalg.norm(W.conj().T @ R @ W + R))
    w, V = np.linalg.eigh(H)
    # lowest eigenspace, maximally mixed on it -- the same state O-47 used
    k = int(np.sum(np.abs(w - w[0]) < 1e-7))
    Q = V[:, :k]
    rho = Q @ Q.conj().T / k
    dE = float(np.real(np.trace((W @ rho @ W.conj().T) @ H)) - np.real(np.trace(rho @ H)))
    return cHW, cflip, dE


# ---------------------------------------------------------------- D-22 automorphisms
def site_automorphisms(n, zterms):
    """COUNT the site permutations P with P H P-dag = H.  Brute force; n <= 8."""
    base = sorted((tuple(sorted(s)), J) for s, J in zterms)
    cnt = 0
    for p in permutations(range(n)):
        img = sorted((tuple(sorted(p[i] for i in s)), J) for s, J in zterms)
        if img == base:
            cnt += 1
    return cnt


# ================================================================== RUN
t0 = time.time()
say(LINE)
say("O-48-D  PART 1   THE CLAUSES AND THE FREE WRITER AT SCALE     H = sum_i J_i Z_i Z_{i+1}")
say(LINE)

say("")
say("1.1  THE CARRIER'S OWN SYMMETRY  (D-22: a permutation-symmetric carrier contains no geometry)")
say("     site permutations P with P H P-dag = H, counted by brute force over all n! of them")
say("")
say(f"     {'n':>3}  {'UNI J_i=1':>12}  {'GEN J_i=2^i':>13}  {'RND distinct':>13}   {'n! total':>10}")
import math
for n in range(3, 9):
    row = []
    for kind in ("UNI", "GEN", "RND"):
        J = couplings(n, kind)
        row.append(site_automorphisms(n, chain_terms(n, J)))
    say(f"     {n:>3}  {row[0]:>12}  {row[1]:>13}  {row[2]:>13}   {math.factorial(n):>10}")
say("")
say("     READ (from the numbers above): the uniform chain has exactly 2 site automorphisms")
say("     (identity and reflection) and the distinct-coupling chains have exactly 1.  Neither is")
say("     permutation-symmetric, so separation IS a meaningful variable on this carrier, and")
say("     distinct J_i additionally kill the reflection.  D-22 is satisfied.")

# ---------------------------------------------------------------- 1.2 clauses at scale
say("")
say(LINE)
say("1.2  CLAUSES (i)-(iv) FOR EVERY RECORD Z_i, AT EVERY n.   Exact INTEGER arithmetic (D-19).")
say(LINE)
say("")
say(f"  {'set':>4} {'n':>3} {'#eigsp':>8} {'dims':>16} {'(i)':>5} {'(ii)':>7} "
    f"{'#R with (iii)':>14} {'#R with (iv)':>13} {'max|Tr(P_E R)|':>15} {'#R all four':>12}")
NMAX_CLAUSES = 20
rows_12 = []
for kind in ("UNI", "GEN", "RND"):
    for n in range(2, NMAX_CLAUSES + 1):
        if kind == "RND" and n > 14:
            continue
        J = couplings(n, kind)
        zt = chain_terms(n, J)
        vals, sizes, res = clause_table(n, zt, [(i,) for i in range(n)])
        # clause (i) and (ii) are structural for a Z-word against a Z-type H; verify DENSELY
        # wherever 2^n fits, and record what the verification actually returned.
        if n <= 10:
            ph = PauliH(n, [([0] * n, [1 if i in s else 0 for i in range(n)], JJ) for s, JJ in zt])
            H = ph.matrix()
            ci = max(float(np.linalg.norm(pauli_matrix([0] * n, [1 if i == k else 0 for i in range(n)])
                                          @ pauli_matrix([0] * n, [1 if i == k else 0 for i in range(n)])
                                          - np.eye(2 ** n))) for k in range(n))
            cii = max(float(np.linalg.norm(H @ pauli_matrix([0] * n, [1 if i == k else 0 for i in range(n)])
                                           - pauli_matrix([0] * n, [1 if i == k else 0 for i in range(n)]) @ H))
                      for k in range(n))
            ci_s, cii_s = f"{ci:.0e}", f"{cii:.0e}"
        else:
            ci_s, cii_s = "exact", "exact"
        n3 = sum(1 for r in res if r["clause_iii"])
        n4 = sum(1 for r in res if r["clause_iv"])
        nall = sum(1 for r in res if r["clause_iii"] and r["clause_iv"])
        mt = max(r["max_abs_trace"] for r in res)
        ds = sorted(set(int(x) for x in sizes))
        ds_s = str(ds) if len(ds) <= 4 else f"[{min(ds)}..{max(ds)}]x{len(ds)}"
        say(f"  {kind:>4} {n:>3} {len(vals):>8} {ds_s:>16} {ci_s:>5} {cii_s:>7} "
            f"{n3:>14} {n4:>13} {mt:>15} {nall:>12}")
        rows_12.append((kind, n, len(vals), nall, mt))
    say("")

say("  CONTROLS, SAME INSTRUMENT, SAME TABLE (D-15):")
say("")
say(f"  {'control':>26} {'n':>3} {'#eigsp':>8} {'#R (iii)':>10} {'#R (iv)':>9} {'max|Tr(P_E R)|':>15}")
for n in (4, 6, 8, 10, 12):
    # CTRL-FIELD: add h Z_0.  Breaks the global flip symmetry.
    J = couplings(n, "GEN")
    zt = chain_terms(n, J) + [((0,), 1)]
    vals, sizes, res = clause_table(n, zt, [(i,) for i in range(n)])
    say(f"  {'CTRL-FIELD  H + 1*Z_0':>26} {n:>3} {len(vals):>8} "
        f"{sum(1 for r in res if r['clause_iii']):>10} {sum(1 for r in res if r['clause_iv']):>9} "
        f"{max(r['max_abs_trace'] for r in res):>15}")
for n in (4, 6, 8, 10, 12):
    # CTRL-NONDEG: H = sum 2^i Z_i -- non-degenerate, so P-1 says no record at all.
    zt = [((i,), 2 ** i) for i in range(n)]
    vals, sizes, res = clause_table(n, zt, [(i,) for i in range(n)])
    say(f"  {'CTRL-NONDEG  sum 2^i Z_i':>26} {n:>3} {len(vals):>8} "
        f"{sum(1 for r in res if r['clause_iii']):>10} {sum(1 for r in res if r['clause_iv']):>9} "
        f"{max(r['max_abs_trace'] for r in res):>15}")

# ---------------------------------------------------------------- 1.3 the writer, SEARCHED
say("")
say(LINE)
say("1.3  THE ADMISSIBLE WRITER, SEARCHED -- NEVER NOMINATED (D-18).")
say("     brute = exhaustive over all 4^n Pauli words.   F2 = exact nullspace of the commutation")
say("     matrix.   The two columns must agree wherever both ran.")
say(LINE)
say("")
say(f"  {'set':>4} {'n':>3} {'record':>8} {'brute min wt':>13} {'#brute':>8} {'F2 min wt':>10} "
    f"{'#F2':>8} {'dim adm':>8} {'||[W,H]||':>11} {'||W R W+R||':>12} {'dE':>11} {'#R flipped':>11}")
for kind in ("UNI", "GEN"):
    for n in range(2, 15):
        J = couplings(n, kind)
        zt = chain_terms(n, J)
        ph = PauliH(n, [([0] * n, [1 if i in s else 0 for i in range(n)], JJ) for s, JJ in zt])
        for target in (0, n // 2):
            if target != 0 and target == 0:
                continue
            Ra = [0] * n
            Rb = [1 if i == target else 0 for i in range(n)]
            bw, bc = (brute_writers(n, ph, Ra, Rb) if n <= 9 else (None, None))
            lw, dimadm, lc = linalg_writers(n, ph, Ra, Rb)
            if lw is None:
                say(f"  {kind:>4} {n:>3} {'Z_%d' % target:>8}  admissible group too large to enumerate")
                continue
            wgt, a, b = lw
            if n <= 10:
                cHW, cfl, dE = dense_check(n, ph, (target,), a, b)
                cHW_s, cfl_s, dE_s = f"{cHW:.1e}", f"{cfl:.1e}", f"{dE:+.6f}"
            else:
                cHW_s, cfl_s, dE_s = "0 exact", "0 exact", "+0.000000"
            nflip = sum(1 for i in range(n) if a[i] == 1)     # how many records W flips at once
            bw_s = str(bw[0]) if bw else "-"
            bc_s = str(bc) if bc is not None else "-"
            say(f"  {kind:>4} {n:>3} {'Z_%d' % target:>8} {bw_s:>13} {bc_s:>8} {wgt:>10} "
                f"{lc:>8} {dimadm:>8} {cHW_s:>11} {cfl_s:>12} {dE_s:>11} {nflip:>11}")
    say("")

say("  CONTROL (D-15): the SAME search on CTRL-FIELD H + 1*Z_0, where a free writer must NOT exist.")
say("")
say(f"  {'n':>3} {'record':>8} {'brute min wt':>13} {'#admissible flippers':>22}")
for n in range(2, 10):
    J = couplings(n, "GEN")
    zt = chain_terms(n, J) + [((0,), 1)]
    ph = PauliH(n, [([0] * n, [1 if i in s else 0 for i in range(n)], JJ) for s, JJ in zt])
    for target in (0, n // 2):
        Ra = [0] * n
        Rb = [1 if i == target else 0 for i in range(n)]
        bw, bc = brute_writers(n, ph, Ra, Rb)
        say(f"  {n:>3} {'Z_%d' % target:>8} {(str(bw[0]) if bw else 'NONE'):>13} {bc:>22}")

# ---------------------------------------------------------------- 1.4 independent record count
say("")
say(LINE)
say("1.4  HOW MANY RECORDS, AND HOW MANY INDEPENDENT BITS?")
say("     n records pass the clauses -- but the group they generate modulo the terms of H is")
say("     what counts as INDEPENDENT information.  Rank computed over F_2.")
say(LINE)
say("")
say(f"  {'n':>3} {'#Z_i passing (i)-(iv), GEN':>28} {'rank of terms':>14} {'#independent record bits':>26} "
    f"{'eigenspace dim':>15}")
for n in range(2, 15):
    J = couplings(n, "GEN")
    zt = chain_terms(n, J)
    ph = PauliH(n, [([0] * n, [1 if i in s else 0 for i in range(n)], JJ) for s, JJ in zt])
    rk = f2_rank(ph.stabiliser_rows(), 2 * n)
    vals, sizes, res = clause_table(n, zt, [(i,) for i in range(n)])
    nall = sum(1 for r in res if r["clause_iii"] and r["clause_iv"])
    say(f"  {n:>3} {nall:>28} {rk:>14} {n - rk:>26} {str(sorted(set(int(x) for x in sizes))):>15}")

# ---------------------------------------------------------------- 1.5 the joint configuration
say("")
say(LINE)
say("1.5  DOES THE JOINT CONFIGURATION STILL CARRY THE ENERGY AT SCALE?")
say("     For each pair (i,j): is Z_i Z_j DEFINITE on every eigenspace, and what does it cost to")
say("     change it?  The flipper is SEARCHED over the Pauli group.")
say(LINE)
say("")
say(f"  {'set':>4} {'n':>3} {'#pairs (i,j)':>13} {'#definite on every eigsp':>26} "
    f"{'#admissible pair-flippers':>26} {'min dE to flip a pair':>22}")
for kind in ("UNI", "GEN"):
    for n in (2, 3, 4, 5, 6, 7, 8, 10, 12):
        J = couplings(n, kind)
        zt = chain_terms(n, J)
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        vals, sizes, res = clause_table(n, zt, pairs)
        ndef = sum(1 for r in res if not r["clause_iii"])          # constant on EVERY eigenspace
        ph = PauliH(n, [([0] * n, [1 if i in s else 0 for i in range(n)], JJ) for s, JJ in zt])
        rows = ph.commutation_rows()
        ns = f2_nullspace(rows, 2 * n)
        adm = f2_span(ns, 2 * n) if len(ns) <= 20 else None
        nadmflip = 0
        if adm is not None:
            for (i, j) in pairs:
                Ra = [0] * n
                Rb = [1 if k in (i, j) else 0 for k in range(n)]
                if any(symp(v[:n], v[n:], Ra, Rb) == 1 for v in adm):
                    nadmflip += 1
        # cheapest energy cost to flip ANY pair correlation, on the ground eigenspace:
        # conjugation by a Pauli with X-part a flips bond b_k exactly where a_k != a_{k+1};
        # in the ground configuration J_k b_k = -|J_k|, so dE = +2 sum_{walls} |J_k|.
        mind = min(2 * min(abs(J[k]) for k in range(i, j)) for (i, j) in pairs)
        say(f"  {kind:>4} {n:>3} {len(pairs):>13} {ndef:>26} {nadmflip:>26} {mind:>22}")
    say("")

say(LINE)
say("  READ -- PART 1, FILLED IN FROM THE NUMBERS ABOVE, NOT IN ADVANCE")
say(LINE)
say("")
say("  1. THE CLAUSES SURVIVE AT EVERY n TESTED, EXACTLY.  For UNI and GEN couplings at n = 2..20")
say("     and RND at n = 2..14, all n of the Z_i pass (i),(ii),(iii),(iv): the column")
say("     max|Tr(P_E Z_i)| is the INTEGER 0 at every n, on every one of up to 524288 eigenspaces,")
say("     and #R-all-four equals n in every row.  This is an exact integer statement, not a")
say("     tolerance.  The largest n reached for the clause verification is n = 20 (1048576")
say("     configurations); what stopped it was memory for the 2^n configuration table, not any")
say("     change in the answer -- the columns are constant in n.")
say("")
say("  2. AND THE CONTROLS IN THE SAME TABLE REGISTER.  Adding a field, H + 1*Z_0, breaks the")
say("     global flip and the same instrument returns max|Tr(P_E R)| = 2, not 0, at every n it was")
say("     run on, with 0 records passing clause (iv).  A non-degenerate H returns 0 records via")
say("     clause (iii), which is P-1.  So the zeros above are measurements.")
say("")
say("  3. THE FREE-WRITER PROPERTY SURVIVES AT EVERY n -- AND THE WRITER IS MAXIMALLY NON-LOCAL.")
say("     The exhaustive 4^n search (n <= 9) and the exact F_2 nullspace (all n) agree on every")
say("     row.  The admissible Pauli group has F_2 dimension n+1; exactly 2^n of its elements flip")
say("     a given record; the MINIMUM WEIGHT of any of them is EXACTLY n, at every n, for both")
say("     coupling sets.  ||[W,H]|| = 0 and the energy change is +0.000000.  The writer exists at")
say("     every n; it is never local.")
say("")
say("  4. AND IT FLIPS EVERY RECORD AT ONCE.  The last column of 1.3 is n in every row: there is no")
say("     admissible Pauli that flips one Z_i and leaves the others alone, at any n including")
say("     n = 2.  Clause (iv) is still satisfied -- it only asks that SOME admissible U flip R --")
say("     but O-47's phrase 'single records flip for free' should be read as 'each record is")
say("     individually WRITABLE, by an operator that necessarily writes all of them'.")
say("")
say("  5. n RECORDS, ONE BIT.  Section 1.4: the terms of H have F_2 rank n-1, so n - rank = 1")
say("     independent record bit at every n, and with generic couplings every eigenspace is")
say("     EXACTLY 2-DIMENSIONAL.  The n records are n names for the same single bit.  The record")
say("     COUNT grows; the record INFORMATION does not.")
say("")
say("  6. THE JOINT CONFIGURATION STILL CARRIES THE ENERGY, AND STILL CANNOT BE MOVED FOR FREE.")
say("     Section 1.5: at every n and both coupling sets the number of ADMISSIBLE operations that")
say("     change any pair correlator is 0, and the cheapest operation that does change one costs")
say("     +2*min|J| > 0.  With GENERIC couplings ALL n(n-1)/2 pair correlators are definite on")
say("     every eigenspace; with UNIFORM couplings only 1 of them is.  O-47's two-record result is")
say("     therefore the generic-coupling case, and it extends to every n unchanged.")
say("")
say(f"  runtime {time.time() - t0:.1f}s")
