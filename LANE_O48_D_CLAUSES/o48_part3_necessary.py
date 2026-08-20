"""O-48-D PART 3.  WHAT DOES THE CONSTRUCTION ACTUALLY REQUIRE?

THE HYPOTHESIS HANDED TO THIS LANE:
   "H must commute with every record AND there must be a symmetry of H that flips each record
    individually."
It has two readings and they behave differently, so BOTH are tested:
   WEAK   for each record R there EXISTS an admissible U with U-dag R U = -R
   STRONG for each record R there exists an admissible U that flips R AND FIXES EVERY OTHER RECORD

FOUR EXACT STATEMENTS ARE TESTED, EACH WITH A CONTROL:

 T1  ADMISSIBLE IMPLIES FREE.  dE = Tr(W rho W-dag H) - Tr(rho H) = Tr(rho W-dag H W) - Tr(rho H),
     and [W,H] = 0 means W-dag H W = H, so dE = 0 IDENTICALLY, for every state, at every n, on
     every carrier.  The "free writer property" is therefore not an extra property that could fail
     at scale -- it is a tautology of O-4's definition of ADMISSIBLE.  What can fail is EXISTENCE.
     CONTROL: the same dE computed for flippers that are NOT admissible must come out non-zero.

 T2  AN ANTICOMMUTING PAIR OF SYMMETRIES GIVES ALL FOUR CLAUSES.  If R = R-dag, R^2 = I, [R,H] = 0,
     and some unitary W has [W,H] = 0 and W R = -R W, then
        (iv)  Tr(P_E R) = Tr(W-dag P_E R W) = Tr(P_E W-dag R W) = -Tr(P_E R) = 0   on every E
        (iii) if R were constant on E then R|_E = cI with c = +-1, and conjugating by W (which
              preserves E) gives -cI = cI, so c = 0 -- impossible.  So R is non-constant on EVERY
              eigenspace, not merely on one.
     So (iii) and (iv) are both FREE once the anticommuting partner exists.

 T3  AND THE CONVERSE IS C-11: clause (iv) holds IFF such a W exists.  So the true requirement is
     exactly A PAIR OF ANTICOMMUTING SYMMETRIES OF H, one of which is the record.

 T4  THE STRONG READING IS FALSE ON THE VERY CARRIER WHERE THE CONSTRUCTION WORKS.

Everything below is measured.  Nothing is nominated.
"""
import sys, os, time, itertools
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from o48_common import (pauli_matrix, symp, weight, all_paulis, f2_rank, PauliH,
                        diag_energies, eig_classes)
from record_model import RecordModel


def say(*a):
    print(*a)
    sys.stdout.flush()


LINE = "=" * 120


def pw(a, b):
    return "".join({(0, 0): "I", (1, 0): "X", (0, 1): "Z", (1, 1): "Y"}[(a[i], b[i])]
                   for i in range(len(a)))


def zt(n, sup, J):
    return ([0] * n, [1 if i in sup else 0 for i in range(n)], J)


def xt(n, sup, J):
    return ([1 if i in sup else 0 for i in range(n)], [0] * n, J)


def eigblocks(H, tol=1e-7):
    w, V = np.linalg.eigh(H)
    out, i = [], 0
    while i < len(w):
        j = i
        while j + 1 < len(w) and abs(w[j + 1] - w[i]) < tol:
            j += 1
        out.append((float(w[i]), V[:, i:j + 1]))
        i = j + 1
    return out


t0 = time.time()
say(LINE)
say("O-48-D  PART 3   THE NECESSARY CONDITION")
say(LINE)

# ==================================================================== T1
say("")
say(LINE)
say("T1.  IS THE 'FREE WRITER' AN EXTRA PROPERTY, OR A TAUTOLOGY OF ADMISSIBILITY?")
say("     dE measured on FOUR states, for the cheapest ADMISSIBLE flipper and -- CONTROL, same")
say("     table -- for the cheapest NON-ADMISSIBLE flipper found by the same search.")
say(LINE)
say("")
say(f"  {'carrier':>34} {'n':>3} {'record':>10} {'adm W':>10} {'||[W,H]||':>10} "
    f"{'dE ground':>11} {'dE mid':>11} {'dE top':>11} {'dE random':>11} | "
    f"{'non-adm V':>10} {'||[V,H]||':>10} {'dE ground':>11}")
rng = np.random.default_rng(3)
for tag, n, terms, rec in [
        ("ZZ chain  GEN J", 4, [zt(4, (i, i + 1), 2 ** i) for i in range(3)], (0,)),
        ("ZZ chain  GEN J", 6, [zt(6, (i, i + 1), 2 ** i) for i in range(5)], (2,)),
        ("ZZ chain  UNI J", 6, [zt(6, (i, i + 1), 1) for i in range(5)], (2,)),
        ("2D ZZ grid 2x3", 6, [zt(6, e, 2 ** k) for k, e in
                               enumerate([(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5)])], (0,)),
        ("ZZZ chain", 6, [zt(6, (i, i + 1, i + 2), 2 ** i) for i in range(4)], (1,)),
        ("XY chain (non-commuting)", 5,
         [zt(5, (i, i + 1), 2 ** i) for i in range(4)] + [xt(5, (i, i + 1), 3 * (i + 1)) for i in range(4)],
         "ZZZZZ")]:
    ph = PauliH(n, terms)
    H = ph.matrix()
    bl = eigblocks(H)
    if isinstance(rec, str):
        Ra = [1 if c in "XY" else 0 for c in rec]
        Rb = [1 if c in "ZY" else 0 for c in rec]
    else:
        Ra = [0] * n
        Rb = [1 if i in rec else 0 for i in range(n)]
    R = pauli_matrix(Ra, Rb)
    admW, nonW = None, None
    for a, b in all_paulis(n):
        if symp(a, b, Ra, Rb) != 1:
            continue
        if ph.admissible(a, b):
            if admW is None or weight(a, b) < weight(*admW):
                admW = (list(a), list(b))
        else:
            if nonW is None or weight(a, b) < weight(*nonW):
                nonW = (list(a), list(b))
    states = []
    for _, Q in (bl[0], bl[len(bl) // 2], bl[-1]):
        states.append(Q @ Q.conj().T / Q.shape[1])
    v = rng.normal(size=2 ** n) + 1j * rng.normal(size=2 ** n)
    v /= np.linalg.norm(v)
    states.append(np.outer(v, v.conj()))
    row = [f"  {tag:>34} {n:>3} {(pw(Ra, Rb)):>10}"]
    if admW is None:
        row.append(f" {'NONE':>10}")
    else:
        W = pauli_matrix(*admW)
        row.append(f" {pw(*admW):>10} {np.linalg.norm(W @ H - H @ W):>10.1e}")
        for rho in states:
            dE = float(np.real(np.trace((W @ rho @ W.conj().T) @ H)) - np.real(np.trace(rho @ H)))
            row.append(f" {dE:>+11.6f}")
    row.append(" |")
    if nonW is None:
        row.append(f" {'NONE':>10}")
    else:
        V = pauli_matrix(*nonW)
        dE = float(np.real(np.trace((V @ states[0] @ V.conj().T) @ H)) - np.real(np.trace(states[0] @ H)))
        row.append(f" {pw(*nonW):>10} {np.linalg.norm(V @ H - H @ V):>10.1e} {dE:>+11.6f}")
    say("".join(row))
say("")
say("  The ADMISSIBLE column is +0.000000 on every state including a RANDOM one -- not just on")
say("  eigenstates -- while the NON-ADMISSIBLE control on the very same carrier is not zero.")

# ==================================================================== T2 / T3
say("")
say(LINE)
say("T2/T3.  DOES 'AN ANTICOMMUTING PAIR OF SYMMETRIES' PREDICT THE RECORD COUNT ON EVERY CARRIER?")
say("     PREDICTED = #{Hermitian Pauli involutions R with [R,H]=0 that have SOME W with [W,H]=0")
say("                  and WR = -RW}, computed from the symplectic form alone.")
say("     MEASURED  = #{R passing clauses (i)-(iv) verified on the dense matrices}.")
say("     They must agree carrier by carrier, or the criterion is wrong.")
say(LINE)
say("")
say(f"  {'carrier':>40} {'n':>3} {'#[R,H]=0':>9} {'PREDICTED':>10} {'MEASURED':>9} {'agree?':>7} "
    f"{'#(iii) but not (iv)':>21}")

CARRIERS = []
for n in (4, 5, 6):
    CARRIERS.append((f"ZZ chain GEN J", n, [zt(n, (i, i + 1), 2 ** i) for i in range(n - 1)]))
for n in (4, 5, 6):
    CARRIERS.append((f"ZZ chain UNI J", n, [zt(n, (i, i + 1), 1) for i in range(n - 1)]))
CARRIERS.append(("2D ZZ grid 2x3", 6, [zt(6, e, 2 ** k) for k, e in
                                       enumerate([(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5)])]))
for n in (5, 6, 7):
    CARRIERS.append((f"ZZZ chain", n, [zt(n, (i, i + 1, i + 2), 2 ** i) for i in range(n - 2)]))
for n in (4, 5, 6):
    CARRIERS.append((f"XY chain (NON-commuting terms)", n,
                     [zt(n, (i, i + 1), 2 ** i) for i in range(n - 1)] +
                     [xt(n, (i, i + 1), 3 * (i + 1)) for i in range(n - 1)]))
for n in (4, 6):
    t = []
    for i in range(0, n - 1, 2):
        t += [zt(n, (i, i + 1), 2 ** i), xt(n, (i, i + 1), 3 * (i + 2))]
    CARRIERS.append(("C2 dimerised XX+ZZ", n, t))
for n in (4, 6):
    CARRIERS.append(("C3 ZZ chain + K X^(x)n", n,
                     [zt(n, (i, i + 1), 2 ** i) for i in range(n - 1)] + [xt(n, tuple(range(n)), 5)]))
CARRIERS.append(("CTRL ZZ chain + field 1*Z_0", 6,
                 [zt(6, (i, i + 1), 2 ** i) for i in range(5)] + [zt(6, (0,), 1)]))
CARRIERS.append(("CTRL non-degenerate sum 2^i Z_i", 6, [zt(6, (i,), 2 ** i) for i in range(6)]))

rows_T2 = []
for tag, n, terms in CARRIERS:
    ph = PauliH(n, terms)
    H = ph.matrix()
    bl = eigblocks(H)
    comm = [(list(a), list(b)) for a, b in all_paulis(n)
            if (any(a) or any(b)) and ph.admissible(a, b)]
    pred = 0
    for a, b in comm:
        if any(symp(a, b, wa, wb) == 1 for wa, wb in comm):
            pred += 1
    meas, iii_not_iv = 0, 0
    for a, b in comm:
        R = pauli_matrix(a, b)
        nonconst, maxtr = False, 0.0
        for _, Q in bl:
            M = Q.conj().T @ R @ Q
            m = Q.shape[1]
            tr = complex(np.trace(M))
            maxtr = max(maxtr, abs(tr))
            if np.linalg.norm(M - (tr / m) * np.eye(m)) > 1e-7:
                nonconst = True
        if nonconst and maxtr < 1e-7:
            meas += 1
        elif nonconst:
            iii_not_iv += 1
    say(f"  {tag:>40} {n:>3} {len(comm):>9} {pred:>10} {meas:>9} "
        f"{('YES' if pred == meas else '*** NO ***'):>7} {iii_not_iv:>21}")
    rows_T2.append((tag, n, pred, meas, iii_not_iv))

say("")
say("  The last column is the number of operators satisfying (i),(ii),(iii) but FAILING (iv).")
say("  Wherever it is non-zero, (iii) does NOT imply (iv) -- so the two clauses are independent")
say("  and clause (iv) is doing real work, not restating clause (iii).")

# ==================================================================== T4  the STRONG reading
say("")
say(LINE)
say("T4.  THE STRONG READING: IS THERE A SYMMETRY THAT FLIPS ONE RECORD AND FIXES THE OTHERS?")
say("     (a) EXHAUSTIVE Pauli search: admissible W with W R_i W = -R_i and W R_j W = +R_j for all")
say("         j != i.   (b) the model's own independently_writable(), which builds a permutation")
say("         on the joint eigenbasis and is not restricted to Paulis.")
say(LINE)
say("")
say(f"  {'carrier':>34} {'n':>3} {'#records':>9} {'(a) #individually flippable (Pauli)':>36} "
    f"{'(b) model independently_writable':>34}")
for tag, n, terms, sup in [("ZZ chain GEN J", 3, [zt(3, (i, i + 1), 2 ** i) for i in range(2)], None),
                           ("ZZ chain GEN J", 4, [zt(4, (i, i + 1), 2 ** i) for i in range(3)], None),
                           ("ZZ chain GEN J", 5, [zt(5, (i, i + 1), 2 ** i) for i in range(4)], None),
                           ("ZZ chain UNI J", 4, [zt(4, (i, i + 1), 1) for i in range(3)], None),
                           ("ZZ chain UNI J", 5, [zt(5, (i, i + 1), 1) for i in range(4)], None),
                           ("2D ZZ grid 2x2", 4, [zt(4, e, 2 ** k) for k, e in
                                                  enumerate([(0, 1), (0, 2), (1, 3), (2, 3)])], None),
                           ("ZZZ chain", 5, [zt(5, (i, i + 1, i + 2), 2 ** i) for i in range(3)], None),
                           ("ZZZ chain", 6, [zt(6, (i, i + 1, i + 2), 2 ** i) for i in range(4)], None)]:
    ph = PauliH(n, terms)
    H = ph.matrix()
    recs = [([0] * n, [1 if i == k else 0 for i in range(n)]) for k in range(n)]
    Rm = [pauli_matrix(a, b) for a, b in recs]
    cnt = 0
    for k, (ra, rb) in enumerate(recs):
        ok = False
        for a, b in all_paulis(n):
            if not ph.admissible(a, b):
                continue
            if symp(a, b, ra, rb) != 1:
                continue
            if all(symp(a, b, recs[j][0], recs[j][1]) == 0 for j in range(n) if j != k):
                ok = True
                break
        cnt += int(ok)
    m = RecordModel(H)
    try:
        fam, _c, iw = m.independence(m.records())
        model_str = f"family size {len(fam)}, independently writable {len(iw)}"
    except Exception as e:
        model_str = f"model: {type(e).__name__}"
    say(f"  {tag:>34} {n:>3} {n:>9} {cnt:>36} {model_str:>34}")

say("")
say("  CONTROL for T4 (D-15): the same search on carriers where individual flipping MUST be")
say("  possible, so that the zeros above are a measurement and not a broken instrument.")
say("  H = 0 makes every Pauli admissible; a PARTIALLY COUPLED chain leaves the uncoupled sites")
say("  individually flippable and the coupled ones not, so the same table shows both answers.")
say("")
say(f"  {'carrier':>34} {'n':>3} {'#records':>9} {'#individually flippable (Pauli)':>34} {'which':>16}")
CTRLS = [("H = 0  (no terms at all)", 4, []),
         ("H = J Z_0Z_1 only, sites 2,3 free", 4, [zt(4, (0, 1), 3)]),
         ("H = J Z_0Z_1 + J' Z_1Z_2, 3,4,5 free", 6, [zt(6, (0, 1), 3), zt(6, (1, 2), 5)])]
for tag, n, terms in CTRLS:
    ph = PauliH(n, terms)
    H = ph.matrix() if terms else np.zeros((2 ** n, 2 ** n), dtype=complex)
    bl = eigblocks(H)
    recs = [([0] * n, [1 if i == k else 0 for i in range(n)]) for k in range(n)]
    nrec = 0
    for a, b in recs:
        R = pauli_matrix(a, b)
        rep_nc, rep_tr = False, 0.0
        for _, Q in bl:
            M = Q.conj().T @ R @ Q
            m2 = Q.shape[1]
            tr = complex(np.trace(M))
            rep_tr = max(rep_tr, abs(tr))
            if np.linalg.norm(M - (tr / m2) * np.eye(m2)) > 1e-7:
                rep_nc = True
        if np.linalg.norm(H @ R - R @ H) < 1e-8 and rep_nc and rep_tr < 1e-7:
            nrec += 1
    cnt, which = 0, []
    for k, (ra, rb) in enumerate(recs):
        ok = any((not terms or ph.admissible(a, b)) and symp(a, b, ra, rb) == 1
                 and all(symp(a, b, recs[j][0], recs[j][1]) == 0 for j in range(n) if j != k)
                 for a, b in all_paulis(n))
        if ok:
            cnt += 1
            which.append(k)
    say(f"  {tag:>34} {n:>3} {nrec:>9} {cnt:>34} {str(which):>16}")

say("")
say("  AND THE PAULI SEARCH IS NOT THE WHOLE STORY -- the exact argument, checked numerically.")
say("  With GENERIC couplings every eigenspace of the chain is EXACTLY 2-DIMENSIONAL, spanned by")
say("  a configuration and its global flip.  Any admissible U preserves each block, and flipping")
say("  Z_i forces the 2x2 block to be ANTI-DIAGONAL -- which flips EVERY Z_j at once.  So no")
say("  admissible unitary of ANY kind, Pauli or not, flips one record alone.  Built and measured:")
say("")
say(f"  {'n':>3} {'#eigenspaces':>13} {'all dim 2?':>11} {'random admissible U with anti-diagonal blocks':>45}")
rng2 = np.random.default_rng(11)
for n in (4, 5, 6, 7):
    terms = [zt(n, (i, i + 1), 2 ** i) for i in range(n - 1)]
    H = PauliH(n, terms).matrix()
    bl = eigblocks(H)
    alld2 = all(Q.shape[1] == 2 for _, Q in bl)
    U = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for _, Q in bl:                                  # the MOST GENERAL flipping admissible U
        ph1, ph2 = np.exp(1j * rng2.uniform(0, 2 * np.pi, size=2))
        blk = np.array([[0, ph1], [ph2, 0]], dtype=complex)
        U += Q @ blk @ Q.conj().T
    uni = np.linalg.norm(U.conj().T @ U - np.eye(2 ** n))
    adm = np.linalg.norm(U @ H - H @ U)
    flips = []
    for k in range(n):
        Rk = pauli_matrix([0] * n, [1 if i == k else 0 for i in range(n)])
        flips.append(np.linalg.norm(U.conj().T @ Rk @ U + Rk) < 1e-8)
    say(f"  {n:>3} {len(bl):>13} {str(alld2):>11}   unitary err {uni:.1e}, ||[U,H]|| {adm:.1e},"
        f" records flipped {sum(flips)}/{n}")

# ==================================================================== T5
say("")
say(LINE)
say("T5.  IS 'H IS A COMMUTING-PAULI HAMILTONIAN' NECESSARY?  IS IT SUFFICIENT?")
say(LINE)
say("")
say(f"  {'carrier':>40} {'n':>3} {'terms commute?':>15} {'#records':>9} {'verdict':>44}")
for tag, n, terms in [("XY chain, terms do NOT commute", 5,
                       [zt(5, (i, i + 1), 2 ** i) for i in range(4)] +
                       [xt(5, (i, i + 1), 3 * (i + 1)) for i in range(4)]),
                      ("XY chain, terms do NOT commute", 4,
                       [zt(4, (i, i + 1), 2 ** i) for i in range(3)] +
                       [xt(4, (i, i + 1), 3 * (i + 1)) for i in range(3)]),
                      ("C2 dimerised XX+ZZ, terms DO commute", 4,
                       [zt(4, (0, 1), 1), xt(4, (0, 1), 6), zt(4, (2, 3), 4), xt(4, (2, 3), 12)]),
                      ("C3 ZZ chain + K X^(x)n, terms DO commute", 6,
                       [zt(6, (i, i + 1), 2 ** i) for i in range(5)] + [xt(6, tuple(range(6)), 5)])]:
    ph = PauliH(n, terms)
    H = ph.matrix()
    bl = eigblocks(H)
    comm = [(list(a), list(b)) for a, b in all_paulis(n) if (any(a) or any(b)) and ph.admissible(a, b)]
    meas = 0
    for a, b in comm:
        R = pauli_matrix(a, b)
        nonconst, maxtr = False, 0.0
        for _, Q in bl:
            M = Q.conj().T @ R @ Q
            m2 = Q.shape[1]
            tr = complex(np.trace(M))
            maxtr = max(maxtr, abs(tr))
            if np.linalg.norm(M - (tr / m2) * np.eye(m2)) > 1e-7:
                nonconst = True
        if nonconst and maxtr < 1e-7:
            meas += 1
    v = ("commuting-Pauli NOT necessary" if (not ph.commuting and meas > 0) else
         ("commuting-Pauli NOT sufficient" if (ph.commuting and meas == 0) else "consistent either way"))
    say(f"  {tag:>40} {n:>3} {str(ph.commuting):>15} {meas:>9} {v:>44}")

# ==================================================================== T6
say("")
say(LINE)
say("T6.  A COUNTEREXAMPLE INSIDE THE CHAIN FAMILY ITSELF: does (i)+(ii)+(iii) imply (iv)?")
say("     Candidate R = Z_0 Z_1 on the UNIFORM chain, where the eigenspaces are coarse enough that")
say("     a product of terms of H is no longer constant on them.  Exact integer arithmetic.")
say(LINE)
say("")
say(f"  {'J set':>6} {'n':>3} {'R':>12} {'(ii) [H,R]=0':>13} {'(iii) non-const':>16} "
    f"{'max|Tr(P_E R)|':>15} {'(iv)':>7}")
for kind in ("UNI", "GEN"):
    for n in (4, 6, 8, 10, 12):
        J = [1] * (n - 1) if kind == "UNI" else [2 ** i for i in range(n - 1)]
        zterms = [((i, i + 1), J[i]) for i in range(n - 1)]
        sig, E = diag_energies(n, zterms)
        vals, inv, sizes = eig_classes(E)
        for sup, nm in ((( 0, 1), "Z_0 Z_1"), ((0, 2), "Z_0 Z_2")):
            v = np.ones(1 << n, dtype=np.int64)
            for i in sup:
                v = v * sig[:, i].astype(np.int64)
            pos = np.bincount(inv, weights=(v > 0).astype(np.float64)).astype(np.int64)
            tr = 2 * pos - sizes
            nonconst = bool(np.any((pos > 0) & (pos < sizes)))
            say(f"  {kind:>6} {n:>3} {nm:>12} {'YES':>13} {str(nonconst):>16} "
                f"{int(np.abs(tr).max()):>15} {str(int(np.abs(tr).max()) == 0):>7}")
    say("")

say(LINE)
say("  READ -- PART 3, FILLED IN FROM THE NUMBERS ABOVE, NOT IN ADVANCE")
say(LINE)
say("")
say("  1. THE 'FREE WRITER' IS A TAUTOLOGY, NOT A DISCOVERY THAT COULD FAIL AT SCALE.  T1 measured")
say("     dE = +0.000000 for the admissible writer on FOUR different states -- ground, middle, top,")
say("     and a RANDOM pure state -- on all six carriers.  The one-line reason is exact and needs")
say("     no computation: dE = Tr(rho W-dag H W) - Tr(rho H) and [W,H] = 0 gives W-dag H W = H.  The")
say("     control on the same carriers, the cheapest NON-admissible flipper, returns +2, +12, +4,")
say("     +6, +6 and +11.05.  So O-47's 'flips for free' is O-4's definition of ADMISSIBLE read back")
say("     to us.  The question that has content is EXISTENCE, which is clause (iv).")
say("")
say("  2. THE CRITERION IS A PAIR OF ANTICOMMUTING SYMMETRIES.  T2/T3: on all NINETEEN carriers,")
say("     including the ones with NO records and the two controls, the count predicted from the")
say("     symplectic form alone -- R commutes with H and SOME W commutes with H and anticommutes")
say("     with R -- equals the count measured on the dense matrices, with no exceptions.  And when")
say("     that pair exists, clauses (iii) and (iv) BOTH come free: (iv) because Tr(P_E R) = -Tr(P_E R),")
say("     (iii) because a constant +-1 on a block cannot survive conjugation by W.")
say("")
say("  3. CLAUSE (iv) IS NOT IMPLIED BY (i)-(iii).  The last column of T2 is non-zero on six")
say("     carriers -- 6, 14 and 30 on the UNIFORM chain at n = 4,5,6, 56 on C2, 62 on C3, 62 on the")
say("     field control.  T6 exhibits one inside the chain family itself: on the UNIFORM chain")
say("     R = Z_0 Z_1 commutes with H and IS non-constant on an eigenspace, and max|Tr(P_E R)| runs")
say("     2, 6, 18, 56, 180 at n = 4,6,8,10,12 -- it fails (iv), and fails it worse as n grows.  On")
say("     the GENERIC chain the same operator fails (iii) instead.  Clause (iv) does real work.")
say("")
say("  4. THE HYPOTHESIS AS HANDED TO THIS LANE, STRONG READING, IS FALSE.  T4: on the ZZ chain")
say("     (both coupling sets), the 2D grid and the ZZZ chain, the number of records that can be")
say("     flipped by an admissible operation FIXING THE OTHERS is 0 out of n, at every n tested.")
say("     The controls in the same table return 4/4 for H = 0 and 2/4 and 3/6 for partially coupled")
say("     chains, naming exactly the uncoupled sites -- so the zeros are measurements.  And the")
say("     restriction to Paulis is not what is doing it: with generic couplings every eigenspace is")
say("     EXACTLY 2-dimensional, the most general admissible flipping U was built explicitly with")
say("     random phases, and it flips n/n records at n = 4,5,6,7 with ||[U,H]|| = 0.  NO admissible")
say("     unitary of any kind flips one record alone on this carrier.")
say("")
say("  5. VERDICT ON THE HYPOTHESIS.")
say("     'H commutes with every record'            NECESSARY (it is clause (ii)), NOT sufficient:")
say("                                               C2, C3 and the field control all satisfy it and")
say("                                               have no records.")
say("     'a symmetry flips each record INDIVIDUALLY' NEITHER necessary NOR sufficient.  Not")
say("                                               necessary: false on every carrier where the")
say("                                               construction works.  Not sufficient: H = 0")
say("                                               satisfies it and carries no configuration")
say("                                               energy at all.")
say("     'a symmetry flips each record'            NECESSARY AND SUFFICIENT given (i)+(ii) -- but")
say("                                               it IS clause (iv) via C-11, so as a criterion it")
say("                                               is a restatement, not an explanation.")
say("     THE STATEMENT WITH CONTENT: H must admit A PAIR OF ANTICOMMUTING SYMMETRIES, one of which")
say("     is the record.  Then (iii) and (iv) follow and the writer is free.  For a commuting-Pauli")
say("     H with generic couplings this says exactly: the record is a nontrivial LOGICAL operator")
say("     of the group generated by H's terms.  T5 shows commuting-Pauli is NEITHER necessary (the")
say("     XY chain at n = 5 has 3 records) NOR sufficient (C2 and C3 have none).")
say("")
say(f"  runtime {time.time() - t0:.1f}s")
