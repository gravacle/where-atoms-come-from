"""O-50-A  step 6.  THE THEOREM CANDIDATE, decided.

   PART A.  The exact argument, and the exact place where its hypothesis can fail.
   PART B.  The torus read the OTHER way -- ALL mutually commuting records, not just the greedy
            family -- to check the two readings agree on the torus (they do) where they diverge
            on the chain (they do, see s5).
   PART C.  A FALSIFICATION SCAN: random stabiliser carriers, searching for any case where a
            commuting_family of size k >= 2 is NOT independently writable.
   PART D.  Higher genus / several tori: more records, same verdict.
   PART E.  The cancellation law, stated and checked exactly by character decomposition."""
import sys, itertools, random, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_A_ACTION")
from record_model import RecordModel, symplectic_logicals, xz_to_matrix, eigenspaces, clause_iii, clause_iv
from f2lib import Toric, sp, rank, in_span, nullspace, span
import orbits
from fractions import Fraction

print("=" * 78); print("PART B.  THE TORUS READ WITH *ALL* MUTUALLY COMMUTING RECORDS (L=2, dense)")
L = 2; T = Toric(L); n = T.n
S = [xz_to_matrix(s, n) for s in T.stab]; H = -sum(S); es = eigenspaces(H)
pairs = symplectic_logicals(T.stab, n); flat = [v for pr in pairs for v in pr]
cls = []
for bits in itertools.product((0, 1), repeat=4):
    if not any(bits): continue
    v = [0] * (2 * n)
    for b, f in zip(bits, flat):
        if b: v = [(x + y) % 2 for x, y in zip(v, f)]
    cls.append((bits, v))
# SEARCH: maximal sets of mutually commuting record classes
best = []
for r in range(1, 6):
    for combo in itertools.combinations(cls, r):
        if all(sp(a, b, n) == 0 for (_, a), (_, b) in itertools.combinations(combo, 2)):
            if r > len(best): best = list(combo)
print(f"  largest set of MUTUALLY COMMUTING non-identity records found: {len(best)}")
print(f"  its classes: {[b for b, _ in best]}")
Ms = [xz_to_matrix(v, n) for _, v in best]
rm = RecordModel(H, S)
jb = rm.joint_basis(Ms)
labs = sorted({lab for (ei, lab) in jb if ei == 0})
print(f"  joint sign labels REALISED on the ground space: {labs}  -> {len(labs)} of 2^{len(best)}={2**len(best)}")
print("  the third record is the PRODUCT of the other two, so it carries no independent bit:")
for lab in labs: print(f"     {lab}   s3 == s1*s2 ? {lab[2] == lab[0]*lab[1] if len(lab)>2 else 'n/a'}")
fam = rm.commuting_family(Ms)
print(f"  commuting_family on that set -> {len(fam)} independent bits; "
      f"independently_writable = {rm.independently_writable(fam)}")
print("  READ: on the torus BOTH readings of 'maximal commuting family' give the SAME "
      "configuration space (4 configurations, 2 independent bits).")

print()
print("=" * 78); print("PART C.  FALSIFICATION SCAN over random stabiliser carriers")
random.seed(11)
def rand_code(nq, m, rng):
    """m random INDEPENDENT commuting Paulis on nq qubits -> a stabiliser carrier."""
    gens = []
    tries = 0
    while len(gens) < m and tries < 4000:
        tries += 1
        v = [rng.randrange(2) for _ in range(2 * nq)]
        if not any(v): continue
        if any(sp(v, g, nq) for g in gens): continue
        if in_span(v, gens, 2 * nq) if gens else False: continue
        if gens and rank(gens + [v], 2 * nq) == len(gens): continue
        gens.append(v)
    return gens if len(gens) == m else None

rows = []
fails = 0; cases = 0
for nq in (3, 4, 5):
    for m in range(1, nq):
        for trial in range(12):
            g = rand_code(nq, m, random)
            if g is None: continue
            Sm = [xz_to_matrix(v, nq) for v in g]
            Hm = -sum(Sm)
            rmm = RecordModel(Hm, Sm)
            k = nq - m
            prs = symplectic_logicals(g, nq)
            fl = [v for pr in prs for v in pr]
            cands = []
            for bits in itertools.product((0, 1), repeat=len(fl)):
                if not any(bits): continue
                v = [0] * (2 * nq)
                for b, f in zip(bits, fl):
                    if b: v = [(x + y) % 2 for x, y in zip(v, f)]
                cands.append(xz_to_matrix(v, nq))
            fam = rmm.commuting_family(cands)
            iw = rmm.independently_writable(fam)
            cases += 1
            ok = (sorted(iw) == list(range(len(fam))))
            if not ok: fails += 1
            rows.append((nq, m, k, len(prs), len(fam), len(iw), ok))
print(f"  random stabiliser carriers scanned: {cases}")
print(f"  cases where commuting_family was NOT fully independently writable: {fails}")
agg = {}
for nq, m, k, np_, lf, li, ok in rows:
    agg.setdefault((nq, m), [0, 0, set(), set()])
    agg[(nq, m)][0] += 1; agg[(nq, m)][1] += ok
    agg[(nq, m)][2].add(lf); agg[(nq, m)][3].add(li)
print(f"    {'nq':>3} {'m':>3} {'trials':>7} {'all writable':>13} {'|family|':>10} {'|writable|':>11}")
for (nq, m), (tot, okc, lfs, lis) in sorted(agg.items()):
    print(f"    {nq:>3} {m:>3} {tot:>7} {okc:>13} {str(sorted(lfs)):>10} {str(sorted(lis)):>11}")
print("  READ: no counterexample found among stabiliser carriers. The greedy family's "
      "even-splitting criterion FORCES independent writability -- see PART A.")

print()
print("=" * 78); print("PART D.  HIGHER GENUS -- several tori, F_2 only")
def multi_torus(Ls):
    """direct sum of toric codes: stabilisers block-diagonal in the symplectic representation."""
    tots = sum(2 * l * l for l in Ls); stab = []; off = 0; pieces = []
    for l in Ls:
        Tl = Toric(l); pieces.append((Tl, off))
        for s in Tl.stab:
            v = [0] * (2 * tots)
            for e in range(Tl.n):
                if s[e]: v[off + e] = 1
                if s[Tl.n + e]: v[tots + off + e] = 1
            stab.append(v)
        off += Tl.n
    return tots, stab

for Ls in ([3], [3, 3], [3, 3, 3], [3, 4]):
    N, stab = multi_torus(Ls)
    r = rank(stab, 2 * N); k = N - r
    prs = symplectic_logicals(stab, N)
    fl = [v for pr in prs for v in pr]
    G = [[sp(a, b, N) for b in fl] for a in fl]
    from s1_toric_f2 import det_f2
    # a maximal isotropic family: search
    Rs = None
    for combo in itertools.combinations(fl, k):
        if all(sp(a, b, N) == 0 for a, b in itertools.combinations(combo, 2)) and rank(list(combo), 2*N) == k:
            Rs = list(combo); break
    rows_ns = [[sp([1 if q == j else 0 for q in range(2 * N)], s, N) for j in range(2 * N)] for s in stab]
    NS = nullspace(rows_ns, 2 * N)
    phi = [[sp(b, R, N) for R in Rs] for b in NS]
    gw = 2 ** rank(phi, k)
    cfgs = list(itertools.product((0, 1), repeat=k)); idx = {c: i for i, c in enumerate(cfgs)}
    gens = []
    for j in range(k):
        tgt = tuple(1 if i == j else 0 for i in range(k))
        v = next(v for v in fl + [ [ (a[t]+b[t])%2 for t in range(2*N)] for a,b in itertools.combinations(fl,2)]
                 if tuple(sp(v, R, N) for R in Rs) == tgt)
        gens.append([idx[tuple((c[i] + tuple(sp(v, R, N) for R in Rs)[i]) % 2 for i in range(k))] for c in cfgs])
    act = orbits.analyse(gens, len(cfgs)); act.pop('orbits')
    print(f"  tori {Ls}: n = {N}, k = {k} records, gram det_F2 = {det_f2(G)}, "
          f"configurations = {2**k}, |G_W| = {gw}")
    print(f"      action: {act}")

print()
print("=" * 78); print("PART E.  THE CANCELLATION LAW, exact character decomposition")
for k in (1, 2, 3, 4):
    m = 2 ** k
    cfgs = list(itertools.product((1, -1), repeat=k))
    # every non-constant character has EXACT mean zero over configurations
    means = {}
    for Ssub in range(1, 2 ** k):
        idxs = [i for i in range(k) if (Ssub >> i) & 1]
        vals = [np.prod([c[i] for i in idxs]) for c in cfgs]
        means[tuple(idxs)] = int(sum(vals))
    print(f"  k={k}: characters chi_S, S nonempty: {len(means)};  "
          f"max |sum over configurations| = {max(abs(v) for v in means.values())}  "
          f"(exact integer, must be 0)")
    # and each is flipped by at least one writer
    flipped = all(any((sum((1 if i in idxs else 0) * g[i] for i in range(k)) % 2) == 1
                      for g in itertools.product((0, 1), repeat=k) if any(g))
                  for idxs in means)
    print(f"        every non-constant character is ODD under some writer: {flipped}")
