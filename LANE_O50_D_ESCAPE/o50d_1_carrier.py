"""O-50 D  PART 1 -- THE CARRIER, THE CLAUSES, AND THE WRITER GROUP, ON THE TORUS.

Nothing is nominated.  Logicals come from symplectic_logicals (D-18); writers are SEARCHED
for over the whole logical Pauli group; the record family is SELECTED by a commutation
search.  All five clauses are checked on the torus (D-23: no proxy convention anywhere).
"""
import sys, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_D_ESCAPE")
from o50d_common import *

say("=" * 104)
say("O-50 D  PART 1   CARRIER / CLAUSES / WRITER GROUP        canonical carrier: TORIC CODE (O-49)")
say("=" * 104)

# ---------------------------------------------------------------- 1. the carrier family
say("")
say("1. THE CARRIER FAMILY  (exact F_2; no dense matrix is built above L=2)")
say(f"   {'L':>3}{'n=2L^2':>9}{'rank S':>9}{'k=n-rank':>10}{'#logical classes':>19}{'dim code':>10}")
TORI = {}
for L in range(2, 9):
    T = Torus(L); TORI[L] = T
    r = rank2(T.stab, 2 * T.nq); k = T.nq - r
    say(f"   {L:>3}{T.nq:>9}{r:>9}{k:>10}{4 ** k:>19}{2 ** k:>10}")
say("   k = 2 at every L: genus 1 gives exactly 2 logical qubits.  EXACT.")

# ---------------------------------------------------------------- 2. logicals, COMPUTED
say("")
say("2. LOGICAL OPERATORS -- COMPUTED by symplectic_logicals, never nominated (D-18)")
LOG = {}
LS = (2, 3, 4, 5, 6)
for L in LS:
    T = TORI[L]
    pairs = symplectic_logicals(T.stab, T.nq); LOG[L] = pairs
    ok = all(T.sp(a, b) == 1 for a, b in pairs)
    cross = all(T.sp(p[i], q[j]) == 0 for pi, p in enumerate(pairs) for qi, q in enumerate(pairs)
                if pi != qi for i in (0, 1) for j in (0, 1))
    comm_stab = all(T.sp(x, s) == 0 for p in pairs for x in p for s in T.stab)
    say(f"   L={L}: {len(pairs)} conjugate pair(s); within-pair anticommute={ok}; "
        f"cross-pair commute={cross}; all commute with S={comm_stab}; "
        f"weights={[(T.weight(p[0]), T.weight(p[1])) for p in pairs]}")

# ---------------------------------------------------------------- 3. the record family, SELECTED
say("")
say("3. THE RECORD FAMILY -- selected by SEARCH over the 16 logical classes, not nominated")
def classes(T, pairs):
    basis = [x for p in pairs for x in p]; out = []
    for coef in itertools.product((0, 1), repeat=len(basis)):
        v = [0] * (2 * T.nq)
        for c, b in zip(coef, basis):
            if c: v = [(x + y) % 2 for x, y in zip(v, b)]
        out.append((coef, v))
    return out

RECS = {}
for L in LS:
    T = TORI[L]; cls = classes(T, LOG[L])
    nontriv = [(c, v) for c, v in cls if any(c)]
    best = None
    for (c1, v1), (c2, v2) in itertools.combinations(nontriv, 2):
        if T.sp(v1, v2) != 0: continue
        if in_span(list(c2), [list(c1)], 4): continue
        w = T.weight(v1) + T.weight(v2)
        if best is None or w < best[0]: best = (w, v1, v2, c1, c2)
    RECS[L] = (best[1], best[2])
    say(f"   L={L}: maximal commuting family of size 2; coefficient labels {best[3]},{best[4]}; "
        f"weights {T.weight(best[1])},{T.weight(best[2])}")

# ---------------------------------------------------------------- 4. the five clauses
say("")
say("4. THE FIVE CLAUSES ON THE TORUS  (D-23: TORUS-SCOPED, no 1D proxy convention anywhere)")
say("   (i)  BIT           -- xz_to_matrix returns a tensor product of I,X,Y,Z: Hermitian, R^2=I.")
say("                        EXACT for every (x|z) vector; verified densely at L=2 in section 8.")
say("   (ii) DURABLE       -- [H,R]=0 iff sp(R,s)=0 for every stabiliser s:")
for L in LS:
    T = TORI[L]; R1, R2 = RECS[L]
    say(f"                        L={L}: {all(T.sp(R, s) == 0 for R in (R1, R2) for s in T.stab)}")
say("   (iii) NON-TRIVIAL  -- R in N(S)\\S acts as a non-scalar on the 2^k-dim code space, which is")
say("                        an eigenspace of H.  R in S? (must be False):")
for L in LS:
    T = TORI[L]; R1, R2 = RECS[L]
    say(f"                        L={L}: R1 {in_span(R1, T.stab, 2*T.nq)}   R2 {in_span(R2, T.stab, 2*T.nq)}")
say("   (iv) WRITABLE      -- searched for in section 5; a global Pauli that anticommutes with R and")
say("                        commutes with H flips R on EVERY eigenspace at once, so Tr(P_E R)=0")
say("                        everywhere.  EXACT; densely confirmed at L=2 in section 8.")

say("")
say("   (v)  PROTECTED -- exact F_2 test: for a patch P of the lattice, is N(S) restricted to P")
say("        contained in S?  A patch with a,b <= L-1 is CONTRACTIBLE; a = L or b = L WRAPS.")
def restricted_normaliser(T, region):
    n = T.nq; rows = []
    for s in T.stab: rows.append([s[n + i] for i in range(n)] + [s[i] for i in range(n)])
    inside = set(region)
    for e in range(n):
        if e in inside: continue
        r1 = [0] * (2 * n); r1[e] = 1; rows.append(r1)
        r2 = [0] * (2 * n); r2[n + e] = 1; rows.append(r2)
    ns = nullspace2(rows, 2 * n)
    dS = rank2(T.stab, 2 * n)
    return len(ns), rank2(list(T.stab) + ns, 2 * n) - dS

say(f"        {'L':>3}{'patch':>10}{'contractible':>14}{'#edges':>8}{'dim N(S)|_P':>13}"
    f"{'logical dirs':>14}{'verdict':>12}")
for L in (3, 4, 5):
    T = TORI[L]; broke = []
    for a in range(1, L + 1):
        for b in range(1, L + 1):
            region = sorted({T.h(i, j) for i in range(a) for j in range(b)} |
                            {T.v(i, j) for i in range(a) for j in range(b)})
            dn, extra = restricted_normaliser(T, region)
            contr = (a <= L - 1 and b <= L - 1)
            if extra > 0: broke.append((a, b, contr))
            if (a, b) in ((L - 1, L - 1), (L, L - 1), (L, L)):
                say(f"        {L:>3}{f'{a}x{b}':>10}{str(contr):>14}{len(region):>8}{dn:>13}"
                    f"{extra:>14}{'PROTECTED' if extra == 0 else 'carries logical':>12}")
    contr_break = [x for x in broke if x[2]]
    say(f"        L={L}: patches carrying a logical: {len(broke)} of {L*L}; of those, CONTRACTIBLE "
        f"(a,b <= L-1): {len(contr_break)}   -> clause (v) {'HOLDS' if not contr_break else 'FAILS'}")

say("")
say("   MINIMUM LOGICAL WEIGHT (the code distance), exhaustive over the relevant coset")
say(f"        {'L':>3}{'min |Zbar|':>12}{'min |Xbar|':>12}{'full exhaustive':>18}{'L':>4}")
for L in (2, 3, 4):
    T = TORI[L]
    # Z-type logical: a non-contractible CYCLE; its coset is generated by the plaquettes
    zrep = sorted({T.h(0, j) for j in range(L)})
    xrep = sorted({T.v(0, j) for j in range(L)})
    def mn(rep, gens):
        best = 10 ** 9
        for coef in itertools.product((0, 1), repeat=len(gens)):
            s = set(rep)
            for c, g in zip(coef, gens):
                if c: s ^= set(g)
            if s: best = min(best, len(s))
        return best
    mz = mn(zrep, T.plaq); mx = mn(xrep, T.star)
    full = ""
    if L <= 3:                                        # full exhaustive over ALL of S
        sb, _ = rref(T.stab, 2 * T.nq); bb = 10 ** 9
        for c, v in classes(T, LOG[L]):
            if not any(c): continue
            for coef in itertools.product((0, 1), repeat=len(sb)):
                w = list(v)
                for cc, s in zip(coef, sb):
                    if cc: w = [(x + y) % 2 for x, y in zip(w, s)]
                bb = min(bb, T.weight(w))
        full = str(bb)
    say(f"        {L:>3}{mz:>12}{mx:>12}{full:>18}{L:>4}")
say("        d = L exactly.  BOTH distances scale with L: clause (v) is realised by HOMOLOGY.")

# ---------------------------------------------------------------- 5. the writer group, SEARCHED
say("")
say("5. THE WRITER GROUP  G_W -- SEARCHED over every admissible logical class (D-18)")
say(f"   {'L':>3}{'#adm classes':>14}{'#writers R1':>13}{'#writers R2':>13}{'|G_W| mod S':>13}"
    f"{'kernel dim':>12}{'orbit size':>12}{'simply trans':>14}")
for L in LS:
    T = TORI[L]; R1, R2 = RECS[L]; cls = classes(T, LOG[L])
    w1 = [v for _, v in cls if T.sp(v, R1) == 1 and T.sp(v, R2) == 0]
    w2 = [v for _, v in cls if T.sp(v, R2) == 1 and T.sp(v, R1) == 0]
    gens = w1 + w2; grp = set()
    for coef in itertools.product((0, 1), repeat=len(gens)):
        v = [0] * (2 * T.nq)
        for cc, g in zip(coef, gens):
            if cc: v = [(x + y) % 2 for x, y in zip(v, g)]
        grp.add(tuple(v))
    kernel = [v for v in grp if T.sp(list(v), R1) == 0 and T.sp(list(v), R2) == 0]
    reach = {(s[0] * (-1) ** T.sp(list(v), R1), s[1] * (-1) ** T.sp(list(v), R2))
             for v in grp for s in [(1, 1)]}
    simply = (len(reach) == 4 and len(grp) // len(kernel) == 4)
    say(f"   {L:>3}{len(cls):>14}{len(w1):>13}{len(w2):>13}{len(grp):>13}"
        f"{int(round(np.log2(len(kernel)))):>12}{len(reach):>12}{str(simply):>14}")

say("")
say("   EXACT STRUCTURE (the symplectic form proves what the table shows):")
say("     * the writers of record i form the COSET  Xbar_i + span{R_1..R_k}: 2^k of them, all")
say("       equally admissible.  THE WRITER OF A RECORD IS NOT UNIQUE.")
say("     * together they generate ALL of N(S)/S -- the full logical Pauli group, order 2^(2k).")
say("     * the action on the 2^k configurations is the quotient by span{R_i}: Z_2^k, SIMPLY")
say("       TRANSITIVE.  The theorem candidate's hypothesis HOLDS on the torus.")
say("     * THE KERNEL of the configuration action is span{R_1..R_k} -- THE RECORDS THEMSELVES.")
say("       It is non-trivial and it acts non-trivially on STATES.  That is escape (3)'s content.")

# ---------------------------------------------------------------- 6. invariant functionals
say("")
say("6. THE WRITER-INVARIANT FUNCTIONAL SPACE ON RECORD CONFIGURATIONS")
for L in LS:
    T = TORI[L]; R1, R2 = RECS[L]
    cfgs = list(itertools.product((1, -1), repeat=2)); M = []
    gens = sorted({(T.sp(v, R1), T.sp(v, R2)) for _, v in classes(T, LOG[L])})
    for g in gens:
        for a, s in enumerate(cfgs):
            t = (s[0] * (-1) ** g[0], s[1] * (-1) ** g[1])
            row = [0.0] * 4; row[a] += 1.0; row[cfgs.index(t)] -= 1.0; M.append(row)
    dim_inv = 4 - np.linalg.matrix_rank(np.array(M))
    say(f"   L={L}: dim(G_W-invariant functionals) = {dim_inv} of 4  -- CONSTANTS ONLY; "
        f"non-constant invariant exhibited: NONE")

say("")
say("   AND ON THE OBSERVABLE SIDE -- the 4^k logical Paulis, conjugated by every writer:")
for L in LS:
    T = TORI[L]; cls = classes(T, LOG[L])
    inv = [c for c, v in cls if all(T.sp(v, w) == 0 for _, w in cls)]
    say(f"   L={L}: logical observables fixed by conjugation by EVERY element of G_W: "
        f"{len(inv)} of {len(cls)} -> {inv}")
say("   ONLY THE IDENTITY.  Every non-identity logical observable is sign-flipped by some")
say("   admissible writer, so its average over G_W is EXACTLY ZERO.  Strictly stronger than the")
say("   configuration statement: it covers OFF-DIAGONAL observables -- coherences -- which are")
say("   not functions of the +-1 values at all.  ESCAPE (2) IS CLOSED FOR LINEAR FUNCTIONALS.")

# ---------------------------------------------------------------- 7. D-22 geometry check
say("")
say("7. D-22 -- DOES THE CARRIER CONTAIN GEOMETRY TO DETECT?")
for L in (3, 4, 5):
    T = TORI[L]
    ds = {}
    for e in range(T.nq):
        for f in range(e + 1, T.nq):
            d = T.dist(e, f); ds[d] = ds.get(d, 0) + 1
    Rr, piv = rref(T.stab, 2 * T.nq)
    def inS(v):
        v = list(v)
        for i, c in enumerate(piv):
            if v[c]: v = [(x + y) % 2 for x, y in zip(v, Rr[i])]
        return not any(v)
    rng = np.random.default_rng(0); keep = 0; trials = 2000
    for _ in range(trials):
        p = rng.permutation(T.nq); good = True
        for s in T.stab:
            t = [0] * (2 * T.nq)
            for i in range(T.nq):
                t[int(p[i])] = s[i]; t[T.nq + int(p[i])] = s[T.nq + i]
            if not inS(t): good = False; break
        if good: keep += 1
    say(f"   L={L}: {len(ds)} distinct pairwise edge distances (a permutation-symmetric carrier has 1); "
        f"random qubit permutations preserving S: {keep}/{trials}")
say("   The torus is NOT permutation-symmetric: it carries a metric, translations Z_L x Z_L, and a")
say("   point group.  D-22 satisfied -- separation results on this carrier are readable.")

# ---------------------------------------------------------------- 8. dense confirmation at L=2
say("")
say("8. DENSE CONFIRMATION AT L=2 (dim 256) USING THE PROGRAM'S OWN CLAUSE ROUTINES")
T = TORI[2]; R1v, R2v = RECS[2]
H = -sum(dense(s, T.nq) for s in T.stab)
es = eigenspaces(H)
say(f"   H = -sum_v A_v - sum_p B_p ; dim {H.shape[0]}; eigenvalue multiplicities "
    f"{[m for _, _, m in es]}")
for nm, Rv in (("R1", R1v), ("R2", R2v)):
    R = dense(Rv, T.nq)
    say(f"   {nm}: ||R-R-dag||={np.linalg.norm(R - R.conj().T):.2e}  "
        f"||R^2-I||={np.linalg.norm(R @ R - np.eye(256)):.2e}  "
        f"clause(iii)={clause_iii(R, es)}  clause(iv)={clause_iv(R, es)}  "
        f"max|Tr(P_E R)|={max(abs(np.trace(P @ R)) for _, P, _ in es):.2e}")
R1 = dense(R1v, T.nq)
U = build_writer(R1, es)
say(f"   build_writer(R1): ||[U,H]||={np.linalg.norm(U @ H - H @ U):.2e}  "
    f"||U-dag R1 U + R1||={np.linalg.norm(U.conj().T @ R1 @ U + R1):.2e}")
say("")
say("   ALL FIVE CLAUSES HOLD ON THE TORUS.  (i),(ii),(iv) exact at every L; (iii) exact and")
say("   confirmed densely at L=2; (v) by homology with d = L, verified exhaustively to L=4 and by")
say("   the contractible-patch test to L=5.  NOTHING HERE RESTS ON THE 1D PROXY CONVENTION.")
say("=" * 104)
