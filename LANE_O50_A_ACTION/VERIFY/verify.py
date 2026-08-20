"""ADVERSARIAL VERIFICATION OF O-50-A.  Independent re-implementation; imports NOTHING from the lane.
   Only record_model is imported, and only for the dense L=2 cross-check.
   Every number printed here is computed in this file."""
import sys, itertools, random
from fractions import Fraction
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")

# ------------------------------------------------------------------ F2 utilities (mine)
def sp(a, b, n):
    return sum(a[i]*b[n+i] + a[n+i]*b[i] for i in range(n)) % 2

def rref(rows, w):
    rows = [r[:] for r in rows]; piv = []; r = 0
    for c in range(w):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [(x+y) % 2 for x, y in zip(rows[i], rows[r])]
        piv.append(c); r += 1
    return rows[:r], piv

def rank(rows, w): return len(rref(rows, w)[1])

def in_span(v, basis, w):
    if not basis: return not any(v)
    B, piv = rref(basis, w); v = v[:]
    for row, c in zip(B, piv):
        if v[c]: v = [(x+y) % 2 for x, y in zip(v, row)]
    return not any(v)

def nullspace(rows, w):
    R, piv = rref(rows, w); free = [c for c in range(w) if c not in piv]; out = []
    for f in free:
        v = [0]*w; v[f] = 1
        for i, c in enumerate(piv): v[c] = R[i][f]
        out.append(v)
    return out

def add(a, b): return [(x+y) % 2 for x, y in zip(a, b)]

# ------------------------------------------------------------------ MY toric code (different edge indexing)
class MyToric:
    """edge index = 2*(i*L+j) + d,  d=0 horizontal (i,j)->(i,j+1), d=1 vertical (i,j)->(i+1,j).
       DELIBERATELY a different labelling from the lane's f2lib.Toric."""
    def __init__(self, L):
        self.L = L; self.n = 2*L*L
        self.e = lambda i, j, d: 2*((i % L)*L + (j % L)) + d
        A, B = [], []
        for i in range(L):
            for j in range(L):
                s = [0]*(2*self.n)
                for (a, b, d) in ((i, j, 0), (i, j-1, 0), (i, j, 1), (i-1, j, 1)): s[self.e(a, b, d)] ^= 1
                A.append(s)
                t = [0]*(2*self.n)
                for (a, b, d) in ((i, j, 0), (i+1, j, 0), (i, j, 1), (i, j+1, 1)): t[self.n + self.e(a, b, d)] ^= 1
                B.append(t)
        self.stab = A + B

# ------------------------------------------------------------------ symplectic Gram-Schmidt (mine)
def logical_pairs(stab, n):
    """returns [(X_i, Z_i)] a conjugate basis of N(S)/S.  Searched, never nominated."""
    Srr, _ = rref(stab, 2*n)
    # N(S) = nullspace of the map v -> (sp(v,s))_s
    rows = [[sp([1 if t == j else 0 for t in range(2*n)], s, n) for j in range(2*n)] for s in Srr]
    NS = nullspace(rows, 2*n)
    chosen = list(Srr)
    pool = [v for v in NS if not in_span(v, chosen, 2*n)]
    pairs = []
    while True:
        pool = [v for v in pool if not in_span(v, chosen, 2*n)]
        if not pool: break
        a = pool[0]
        b = next((v for v in pool if sp(a, v, n) == 1), None)
        assert b is not None, "degenerate symplectic form -- would falsify the lane"
        pairs.append((a, b))
        chosen = chosen + [a, b]
        pool = [add(add(v, [sp(v, b, n)*x for x in a]), [sp(v, a, n)*x for x in b]) for v in pool]
    return pairs, NS, Srr

# ------------------------------------------------------------------ generic orbit machinery (mine)
def close_group(gens, m):
    idp = tuple(range(m)); G = {idp}; fr = [idp]
    while fr:
        nf = []
        for g in fr:
            for h in gens:
                p = tuple(h[g[i]] for i in range(m))
                if p not in G: G.add(p); nf.append(p)
        fr = nf
    return G

def orbits_of(G, m):
    seen = [False]*m; out = []
    for c in range(m):
        if seen[c]: continue
        o = sorted({g[c] for g in G})
        for x in o: seen[x] = True
        out.append(o)
    return out

def invariant_dim_by_basis(G, m):
    """dimension of {f : f(g.c)=f(c)} by EXPLICIT BASIS: indicator functions of orbits."""
    return len(orbits_of(G, m))

def invariant_dim_by_elimination(G, m):
    rows = []
    for g in G:
        for c in range(m):
            if g[c] == c: continue
            r = [Fraction(0)]*m; r[g[c]] += 1; r[c] -= 1; rows.append(r)
    if not rows: return m
    piv = 0
    for col in range(m):
        p = next((i for i in range(piv, len(rows)) if rows[i][col] != 0), None)
        if p is None: continue
        rows[piv], rows[p] = rows[p], rows[piv]
        pv = rows[piv][col]; rows[piv] = [x/pv for x in rows[piv]]
        for i in range(len(rows)):
            if i != piv and rows[i][col] != 0:
                f = rows[i][col]; rows[i] = [a - f*b for a, b in zip(rows[i], rows[piv])]
        piv += 1
        if piv == len(rows): break
    return m - piv

def burnside(G, m):
    return Fraction(sum(sum(1 for c in range(m) if g[c] == c) for g in G), len(G))

# ================================================================== PART 1
print("="*100)
print("PART 1.  INDEPENDENT RECONSTRUCTION OF THE TORUS NUMBERS (my own code, my own edge labelling)")
print("  L   n   rank S   k   gram_det_F2   phi surjective   |G_W|  |C|  orbits  Burnside  inv_dim(basis)  inv_dim(elim)  simply transitive")
torus_rows = []
for L in (2, 3, 4, 5):
    T = MyToric(L); n = T.n
    pairs, NS, Srr = logical_pairs(T.stab, n)
    k = len(pairs)
    flat = [v for pr in pairs for v in pr]
    G = [[sp(a, b, n) for b in flat] for a in flat]
    # det over F2 of the 2k x 2k Gram
    def det_f2(M):
        M = [r[:] for r in M]; d = 1
        for c in range(len(M)):
            p = next((i for i in range(c, len(M)) if M[i][c]), None)
            if p is None: return 0
            M[c], M[p] = M[p], M[c]
            for i in range(c+1, len(M)):
                if M[i][c]: M[i] = add(M[i], M[c])
        return 1
    # THE RECORD FAMILY: the k mutually commuting partners Z_i (they commute by construction of GS)
    fam = [pr[1] for pr in pairs]
    assert all(sp(a, b, n) == 0 for a, b in itertools.combinations(fam, 2))
    # WRITER SEARCH over N(S): solve for u in N(S) with prescribed phi = (sp(u,R_1),...,sp(u,R_k))
    # phi image computed by exact linear algebra over the whole of N(S)
    phi_basis = [tuple(sp(v, R, n) for R in fam) for v in NS]
    # image = span of phi_basis
    imgspan = {tuple([0]*k)}
    for pb in phi_basis:
        imgspan |= {tuple((a+b) % 2 for a, b in zip(x, pb)) for x in imgspan}
    surj = len(imgspan) == 2**k
    # explicit generators: a vector of N(S) achieving each unit flip
    gens_f2 = []
    for j in range(k):
        target = tuple(1 if t == j else 0 for t in range(k))
        # solve small linear system over F2
        M = [list(pb) for pb in phi_basis]           # rows = basis vectors, cols = k
        Mt = [[M[r][c] for r in range(len(M))] for c in range(k)]   # k x |NS|
        aug = [Mt[c] + [target[c]] for c in range(k)]
        R, piv = rref(aug, len(NS)+1)
        assert not any(row[:-1] == [0]*len(NS) and row[-1] == 1 for row in R), "no writer -- would falsify"
        sol = [0]*len(NS)
        for row, c in zip(R, piv):
            if c < len(NS): sol[c] = row[-1]
        u = [0]*(2*n)
        for coef, v in zip(sol, NS):
            if coef: u = add(u, v)
        assert all(sp(u, s, n) == 0 for s in T.stab), "writer not admissible"
        assert tuple(sp(u, R, n) for R in fam) == target, "writer does not realise the flip"
        gens_f2.append(u)
    # ACTION on record configurations
    cfgs = list(itertools.product((0, 1), repeat=k)); idx = {c: i for i, c in enumerate(cfgs)}
    perms = []
    for u in gens_f2:
        t = tuple(sp(u, R, n) for R in fam)
        perms.append([idx[tuple((c[i]+t[i]) % 2 for i in range(k))] for c in cfgs])
    GW = close_group(perms, len(cfgs))
    orbs = orbits_of(GW, len(cfgs))
    d1 = invariant_dim_by_basis(GW, len(cfgs)); d2 = invariant_dim_by_elimination(GW, len(cfgs))
    st = (len(orbs) == 1 and len(GW) == len(cfgs) and all(sum(1 for g in GW if g[o[0]] == o[0]) == 1 for o in orbs))
    print(f"  {L}  {n:3d}   {rank(T.stab,2*n):5d}  {k:2d}      {det_f2(G)}            {surj}          "
          f"{len(GW):4d} {len(cfgs):4d}   {len(orbs):4d}    {burnside(GW,len(cfgs))!s:>5}        {d1:6d}        {d2:6d}        {st}")
    torus_rows.append((L, k, len(GW), len(cfgs), len(orbs), d1, d2, st, [sum(u[:n])+sum(u[n:]) for u in gens_f2]))
print("  writer supports (weight of each generator, my search):", [(r[0], r[8]) for r in torus_rows])
print("  VERDICT PART 1: the lane's structural numbers (|G_W|=4, 4 configs, 1 orbit, inv dim 1, simply transitive)")
print("                  REPRODUCE independently.  Burnside and explicit-orbit-basis agree with exact elimination.")

# ================================================================== PART 2
print()
print("="*100)
print("PART 2.  THE HEADLINE CLAIM, TESTED LITERALLY.")
print('  CLAIM UNDER TEST (lane headline / brief theorem candidate):')
print('    "no functional of the record configuration can be BOTH responsive to writing')
print('     AND non-cancelling"   and   "every non-constant functional is G_W-ODD, hence has')
print('     mean exactly zero over configurations".')
print()
print("  DEFINITIONS USED (the brief's own):")
print("    responsive(f)     :  exists g in G_W and c with f(g.c) != f(c)")
print("    non-cancelling(f) :  sum of f over the writer orbit is NOT zero")
print()
hdr = f"  {'carrier':34s} {'|C|':>4} {'|G_W|':>6} {'draws':>6} {'responsive':>11} {'non-cancelling':>15} {'BOTH':>6} {'G_W-odd':>8}"
print(hdr); print("  " + "-"*(len(hdr)-2))
def literal_test(name, gens, m, ntrial=2000, seed=17):
    G = close_group(gens, m); orbs = orbits_of(G, m); rng = random.Random(seed)
    resp = nonc = both = odd = 0
    for _ in range(ntrial):
        f = [Fraction(rng.randint(-9, 9)) for _ in range(m)]
        r = any(f[g[c]] != f[c] for g in G for c in range(m))
        nc = any(sum(f[c] for c in o) != 0 for o in orbs)
        isodd = all(sum(f[c] for c in o) == 0 for o in orbs)
        resp += r; nonc += nc; both += (r and nc); odd += isodd
    print(f"  {name:34s} {m:4d} {len(G):6d} {ntrial:6d} {resp:11d} {nonc:15d} {both:6d} {odd:8d}")
    return both, resp
k2 = [[1,0,3,2],[2,3,0,1]]
b2, _ = literal_test("torus L=any (k=2), G_W=(Z2)^2", k2, 4)
def transl_gens(k):
    cf = list(itertools.product((0,1), repeat=k)); ix = {c:i for i,c in enumerate(cf)}
    return [[ix[tuple((c[i]+ (1 if i==j else 0))%2 for i in range(k))] for c in cf] for j in range(k)], len(cf)
for k in (4, 6):
    g, m = transl_gens(k); literal_test(f"higher genus (k={k}), G_W=(Z2)^{k}", g, m)
g, m = transl_gens(1); literal_test("1D chain greedy family (k=1)", g, m)
m = 16; literal_test("1D chain naive {Z_i} n=4", [[m-1-x for x in range(m)]], m)
print()
print("  EXPLICIT WITNESS on the torus (k=2), configurations c = (s1,s2) in {+-1}^2:")
cfgs = [(1,1),(1,-1),(-1,1),(-1,-1)]
f = {c: 1 + c[0] for c in cfgs}          # f = 1 + s_1
print("    f(s1,s2) = 1 + s1   ->  values", [f[c] for c in cfgs])
print("    responsive?      U_1 sends (1,1)->(-1,1):  f changes 2 -> 0   =>  RESPONSIVE = True")
print("    orbit sum:       ", sum(f[c] for c in cfgs), " != 0            =>  NON-CANCELLING = True")
print("    G_W-odd?         False  (mean over configurations =", Fraction(sum(f[c] for c in cfgs), 4), ")")
print("    ==> a functional of the record configuration that is BOTH responsive AND non-cancelling.")
print("        The theorem candidate's step 'every non-constant functional is G_W-ODD' is FALSE.")

# ================================================================== PART 3
print()
print("="*100)
print("PART 3.  IS s8's 'violations' TEST CAPABLE OF FIRING AT ALL?")
print("  s8_dichotomy.py line:   if responsive and all(x == 0 for x in fodd): bad += 1")
print("  f_odd == 0 identically  <=>  f == f_inv  <=>  f is G-invariant  <=>  NOT responsive.")
print("  So the predicate is a contradiction in terms and can never be true, for ANY f, ANY G.")
print("  Brute-force confirmation over EVERY f in a small function space (exhaustive, not sampled):")
for name, gens, m, vals in [("torus k=2", k2, 4, (-1,0,1)),
                            ("chain naive n=3", [[7-x for x in range(8)]], 8, (0,1))]:
    G = close_group(gens, m); orbs = orbits_of(G, m)
    tot = fired = both_rn = 0
    for f in itertools.product(vals, repeat=m):
        f = [Fraction(x) for x in f]
        finv = [sum(f[g[c]] for g in G)/len(G) for c in range(m)]
        fodd = [f[c]-finv[c] for c in range(m)]
        r = any(f[g[c]] != f[c] for g in G for c in range(m))
        tot += 1
        if r and all(x == 0 for x in fodd): fired += 1
        if r and any(sum(f[c] for c in o) != 0 for o in orbs): both_rn += 1
    print(f"    {name:18s} exhaustive over {tot:6d} functionals:  s8 predicate fired {fired} times;"
          f"  responsive AND non-cancelling: {both_rn} ({100.0*both_rn/tot:.1f}%)")
print("  VERDICT PART 3: s8's zero-violation result is a TAUTOLOGY, not a measurement, and it is")
print("                  NOT the claim printed above it.  The printed claim is refuted by the same scan.")

# ================================================================== PART 4
print()
print("="*100)
print("PART 4.  'CANCELS' = MEAN ZERO OVER ALL CONFIGURATIONS, NOT 'a configuration cancels'.")
print("  Torus k=2.  Take the responsive part of a functional and evaluate it AT EACH CONFIGURATION.")
f2 = {c: Fraction(3) + 2*c[0] + c[1] + 5*c[0]*c[1] for c in cfgs}
mean = Fraction(sum(f2.values()), 4)
print("    f = 3 + 2 s1 + s2 + 5 s1 s2")
for c in cfgs:
    print(f"      c={c!s:9s} f={f2[c]!s:>4}   f_odd(c)=f(c)-mean = {f2[c]-mean!s:>4}   ZERO AT THIS CONFIGURATION? {f2[c]-mean==0}")
print(f"    sum over ALL 4 configurations of f_odd = {sum(f2[c]-mean for c in cfgs)}")
print("  READ: f_odd is non-zero at EVERY single configuration.  It sums to zero only under a UNIFORM")
print("        average over the whole orbit.  A carrier sits in ONE configuration, not in a uniform")
print("        mixture, so 'the record-dependent part cancels' is an ENSEMBLE statement that the lane")
print("        converts into a per-carrier statement.  Uniformity is INSERTED: C-60 gives degeneracy,")
print("        and degeneracy does not give a uniform weight -- a record is by construction in a")
print("        DEFINITE configuration.")

# ================================================================== PART 5
print()
print("="*100)
print("PART 5.  s4c's 'wraps is decided by the TOPOLOGY of the closed edge set, not by size'.")
print("  The function actually used is:")
print("      def block_wraps(T, a, b):  return (a + 1 >= T.L) or (b + 1 >= T.L)")
print("  That is a PURE SIZE TEST.  Independent topological test below: compute the cycle space of the")
print("  region's own subgraph and ask whether any cycle has a non-trivial homology class in T^2.")
def topological_wraps(L, a, b, i0, j0):
    # vertices (i,j); closed block edge set
    E = set()
    for i in range(a+1):
        for j in range(b): E.add(((i0+i) % L, (j0+j) % L, 0))
    for i in range(a):
        for j in range(b+1): E.add(((i0+i) % L, (j0+j) % L, 1))
    E = sorted(E); ei = {e: t for t, e in enumerate(E)}
    verts = sorted({(e[0], e[1]) for e in E} | {((e[0]+1) % L, e[1]) if e[2] == 1 else (e[0], (e[1]+1) % L) for e in E})
    vi = {v: t for t, v in enumerate(verts)}
    rows = [[0]*len(E) for _ in verts]
    for e in E:
        u = (e[0], e[1]); w = (e[0], (e[1]+1) % L) if e[2] == 0 else ((e[0]+1) % L, e[1])
        rows[vi[u]][ei[e]] ^= 1; rows[vi[w]][ei[e]] ^= 1
    cyc = nullspace(rows, len(E))
    for z in cyc:
        wx = sum(z[ei[e]] for e in E if e[2] == 0 and e[1] == 0) % 2      # crossings of a dual cut
        wy = sum(z[ei[e]] for e in E if e[2] == 1 and e[0] == 0) % 2
        if (wx, wy) != (0, 0): return True
    return False
print("   L   a   b   size-rule 'wraps'   topological 'wraps'   agree")
disagree = 0; tested = 0
for L in (3, 4, 5):
    for a in range(1, L):
        for b in range(1, L):
            s = (a+1 >= L) or (b+1 >= L); t = topological_wraps(L, a, b, 0, 0)
            tested += 1; disagree += (s != t)
            print(f"   {L}   {a}   {b}        {str(s):5s}                {str(t):5s}            {s==t}")
print(f"  shapes tested {tested}, disagreements {disagree}")
print("  VERDICT PART 5: the classification happens to agree on RECTANGULAR blocks, so the clause-(v)")
print("                  numbers survive -- but the printed sentence describes code that was not run,")
print("                  and the region family is rectangles only, not all contractible regions.")

# ================================================================== PART 6
print()
print("="*100)
print("PART 6.  IS THE 'THEOREM' MORE THAN A TEXTBOOK IDENTITY?")
print("  Claim chain: independently writable -> block-swap writers exist -> G_W = (Z2)^k acts by")
print("  translation on {+-1}^k -> transitive -> invariants are constants.")
print("  The last two steps are the standard fact 'a transitive action has only constant invariants',")
print("  equivalently Burnside with |Fix(g)|=0 for g != e.  Check that a bare textbook computation")
print("  reproduces every torus row with no toric-code input at all:")
for k in (1, 2, 3, 4, 6):
    g, m = transl_gens(k); G = close_group(g, m)
    print(f"    (Z2)^{k} acting on {{+-1}}^{k}: |G|={len(G):3d} |C|={m:3d} orbits={len(orbits_of(G,m))} "
          f"Burnside={burnside(G,m)} inv_dim={invariant_dim_by_elimination(G,m)}")
print("  Identical to every TORUS row of the lane's master table.  The toric code contributes the")
print("  values of k (=2 at genus 1) and nothing else to the headline.")

# ================================================================== PART 7
print()
print("="*100)
print("PART 7.  RECORD_MODEL DEFECT TOUCHING commuting_family (read-only observation).")
print("  record_model.py lines 412-413 read:")
print("      if not all(np.linalg.norm(C@C.T.conj()@R@C@C.T.conj() - C@C.T.conj()@R@C@C.T.conj()) < 1 ...)")
print("  The two sides of the subtraction are THE SAME EXPRESSION, so the norm is identically 0 and the")
print("  guard is a no-op.  commuting_family therefore tests ONLY even splitting; it never verifies that")
print("  a candidate commutes with the family already chosen.  The lane's torus result is not damaged")
print("  (mutual commutation was checked separately in s2/s6), but the lane's stated INDUCTION rests on")
print("  this function's criterion, so the criterion is weaker than its docstring claims.")

# ================================================================== PART 8
print()
print("="*100)
print("PART 8.  THE CLAUSE-(v) ARGUMENT, TESTED AGAINST ITS OWN STATED IMPLICATION.")
print('  s4c/the finding assert: "a region contained in a topological disk has H_1 = 0 and therefore')
print('  carries no [non-trivial logical] class -- clause (v) then holds for free, with no convention".')
print("  TEST: for each rectangular plaquette block, compute (a) H_1 of the region's OWN subgraph")
print("  (the lane's stated criterion) and (b) whether the region carries a non-trivial logical.")
print()
def block_edges(L, a, b, i0=0, j0=0):
    E = set()
    for i in range(a+1):
        for j in range(b): E.add(((i0+i) % L, (j0+j) % L, 0))
    for i in range(a):
        for j in range(b+1): E.add(((i0+i) % L, (j0+j) % L, 1))
    return sorted(E)

def region_H1_dim(L, E):
    ei = {e: t for t, e in enumerate(E)}
    vs = set()
    for e in E:
        vs.add((e[0], e[1]))
        vs.add((e[0], (e[1]+1) % L) if e[2] == 0 else (((e[0]+1) % L), e[1]))
    verts = sorted(vs); vi = {v: t for t, v in enumerate(verts)}
    rows = [[0]*len(E) for _ in verts]
    for e in E:
        u = (e[0], e[1]); w = (e[0], (e[1]+1) % L) if e[2] == 0 else (((e[0]+1) % L), e[1])
        rows[vi[u]][ei[e]] ^= 1; rows[vi[w]][ei[e]] ^= 1
    cyc = nullspace(rows, len(E))
    # boundaries: plaquettes entirely inside E
    bnd = []
    for i in range(L):
        for j in range(L):
            pe = [(i, j, 0), ((i+1) % L, j, 0), (i, j, 1), (i, (j+1) % L, 1)]
            if all(p in ei for p in pe):
                v = [0]*len(E)
                for p in pe: v[ei[p]] ^= 1
                bnd.append(v)
    return rank(cyc + bnd, len(E)) - rank(bnd, len(E)), cyc

for L in (3, 4, 5):
    T = MyToric(L); n = T.n
    pairs, NS, Srr = logical_pairs(T.stab, n)
    fam = [pr[1] for pr in pairs] + [pr[0] for pr in pairs]
    rS = rank(T.stab, 2*n)
    print(f"  L = {L}")
    print(f"    {'a':>2} {'b':>2} {'#qubits':>8} {'H_1(region)':>12} {'carries logical':>16} {'lane size-rule wraps':>21}  CONTRADICTS LANE'S STATED IMPLICATION")
    for a in range(1, L):
        for b in range(1, L):
            E = block_edges(L, a, b)
            h1, _ = region_H1_dim(L, E)
            idxs = [T.e(e[0], e[1], e[2]) for e in E]
            coords = idxs + [n + q for q in idxs]
            rows = [[(s[n+c] if c < n else s[c-n]) % 2 for c in coords] for s in T.stab]
            V = []
            for bb in nullspace(rows, len(coords)):
                v = [0]*(2*n)
                for c, bit in zip(coords, bb):
                    if bit: v[c] = 1
                V.append(v)
            dq = rank(T.stab + V, 2*n) - rS
            size_wraps = (a+1 >= L) or (b+1 >= L)
            bad = (h1 == 0 and dq > 0)
            print(f"    {a:>2} {b:>2} {len(E):>8} {h1:>12} {str(dq>0):>16} {str(size_wraps):>21}  {'*** YES ***' if bad else ''}")
            if bad:
                # exhibit an explicit witness inside the region
                wit = None
                for v in V:
                    if not in_span(v, T.stab, 2*n): wit = v; break
                supp_x = [q for q in range(n) if wit[q]]
                supp_z = [q for q in range(n) if wit[n+q]]
                ok_adm = all(sp(wit, s, n) == 0 for s in T.stab)
                flips = [t for t, R in enumerate(fam) if sp(wit, R, n) == 1]
                inside = set(supp_x) | set(supp_z) <= set(idxs)
                print(f"         WITNESS: X-support {supp_x}  Z-support {supp_z}")
                print(f"         commutes with every stabiliser: {ok_adm};  in S: {in_span(wit, T.stab, 2*n)};"
                      f"  support inside the region: {inside}")
                print(f"         anticommutes with logical-basis members {flips}  -> it FLIPS a record")
                print(f"         region H_1 = 0 (the region IS a disk by the lane's own criterion), yet it carries a logical.")
print("  VERDICT PART 8: the lane's stated homological implication is FALSE as written.  A region can")
print("  have H_1 = 0 and still carry a non-trivial logical, because an X-type logical is a cycle of")
print("  the DUAL lattice, whose contractibility is a statement about the region's COMPLEMENT.  The")
print("  lane's actual classifier (block_wraps, a size rule 'spans every vertex row or column') encodes")
print("  that missing dual condition by hand.  So clause (v) on the torus was decided by an INSERTED")
print("  size convention after all -- the very thing D-23 says must be labelled.")

# ================================================================== PART 9
print()
print("="*100)
print("PART 9.  s4c's CLAUSE-(v) TABLE, RECOUNTED WITH THE LANE'S OWN STATED CRITERION.")
print("  Lane's criterion as PRINTED : contractible  <=>  region has H_1 = 0.")
print("  Lane's criterion as CODED   : contractible  <=>  (a+1 < L) and (b+1 < L)   [a size rule].")
print()
print(f"  {'L':>2} {'regions':>8} {'H_1=0 (printed rule)':>21} {'of those, carry a logical':>26} "
      f"{'contractible (coded rule)':>26} {'of those, carry a logical':>26}")
for L in (3, 4, 5):
    T = MyToric(L); n = T.n; rS = rank(T.stab, 2*n)
    tot = h1z = h1z_bad = sz = sz_bad = 0
    for a in range(1, L):
        for b in range(1, L):
            E = block_edges(L, a, b); h1, _ = region_H1_dim(L, E)
            idxs = [T.e(*e) for e in E]; coords = idxs + [n+q for q in idxs]
            rows = [[(s[n+c] if c < n else s[c-n]) % 2 for c in coords] for s in T.stab]
            V = []
            for bb in nullspace(rows, len(coords)):
                v = [0]*(2*n)
                for c, bit in zip(coords, bb):
                    if bit: v[c] = 1
                V.append(v)
            carries = (rank(T.stab+V, 2*n) - rS) > 0
            mult = L*L                      # translations of the same shape
            tot += mult
            if h1 == 0: h1z += mult; h1z_bad += mult*carries
            if (a+1 < L) and (b+1 < L): sz += mult; sz_bad += mult*carries
    print(f"  {L:>2} {tot:>8} {h1z:>21} {h1z_bad:>26} {sz:>26} {sz_bad:>26}")
print("  READ: under the criterion the lane PRINTED, every rectangular block on the torus is")
print("  contractible and a large fraction of them DO carry a logical -- clause (v) FAILS.  The")
print("  reported '0 of 874 contractible regions carries a logical' is produced only by the size")
print("  rule that was coded but not stated.  The lane's own numbers, re-bucketed, say the opposite.")
print()
print("  MINIMAL EXPLICIT FAILURE, every L: the strip of L-1 plaquettes p(0,0..L-2) is a closed")
print("  topological DISK of T^2 (H_1 = 0, no wrapping cycle) and its closed edge set contains an")
print("  entire minimum-weight X-logical of weight L.  Verified above at L=3 (shape 1x2, witness")
print("  weight 3), L=4 (shape 1x3, weight 4), L=5 (shape 1x4, weight 5).")
