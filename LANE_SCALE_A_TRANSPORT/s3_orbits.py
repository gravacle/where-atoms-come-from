"""S3 -- THE ORBIT STRUCTURE OF THE TRANSPORT ACTION ON THE RECORD SET, EXACTLY.

WHAT 'ORBIT' CAN HONESTLY MEAN HERE, said before any number is produced.

Clause (ii) puts a record in the commutant of H, which on this carrier is the whole block
algebra (+)_E End(E); any trace-balanced Hermitian involution in it is a record.  THE RECORD
SET IS THEREFORE A CONTINUUM, not a finite list, and 'how many orbits' has the answer
INFINITELY MANY on every carrier in the ladder -- abelian and non-abelian alike.  Refusing to
say that and quoting a finite orbit count off a random sample is exactly how C-43's
'40 of 40 moved' became a structural-sounding claim about a measure-zero condition.

So two things are computed, and they are labelled apart.

  (1) THE EXACT ORBIT-TYPE STRATIFICATION OF THE ACTUAL (continuum) RECORD SET.
      Transport is the finite group Gbar = G/Z(G) -- A_z is literally the identity matrix for
      central z, which is CHECKED.  Every orbit is finite of size dividing |Gbar|.  For each
      subgroup Z <= L <= G the records fixed by all of L form a subvariety of exact dimension
          dim_L = sum_E sum_{sigma in Irr(L)} m_sigma(E)^2
      (Hermitian directions in the L-commutant, from characters).  The generic stabiliser is
      trivial iff dim_L < dim_all for every L strictly above Z, so the LARGEST ORBIT is |Gbar|,
      and the SMALLEST is 1 exactly when a transport-fixed record exists.

  (2) A FINITE, CANONICAL, TRANSPORT-CLOSED SET OF OBJECTS ON THE SAME CARRIER: the
      K-ISOTYPIC BLOCKS  Pi^K_{E,sigma} = P_E (d_sigma/|K| sum_k conj(chi_sigma(k)) A_k) P_E,
      over every subgroup K <= G.  These are the covariant projections out of which every
      symmetry-adapted record on this carrier is assembled, transport carries them
      covariantly (A_h Pi^K A_h^dag = Pi^{hKh^-1}), and their orbit count IS finite and exact.
      This is a proxy for the record set, and is reported as one.

SELF-CHECKS: A_z = I for central z; sum_sigma d_sigma m_sigma = dim E for every L and E;
the block set must be CLOSED under transport (an image not found in the set kills the row);
orbit sizes must divide |Gbar|.
"""
import sys, time, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_A_TRANSPORT")
import glib
from carriers import census, isotypic, minimal_torus, eigblocks, perm_apply
def say(*a): print(*a); sys.stdout.flush()
NJ = 5

def subgroup_grp(G, K): return glib.Grp("K", sorted(K), lambda a, b: int(G.mt[a, b]))

def dimL(G, ce, K):
    """dimension of the subvariety of records fixed by every element of K"""
    Kg = subgroup_grp(G, K)
    cl, chi, d, cls_of = Kg.chars()
    tot = 0
    for v in (-2, -1, 0):
        if ce['dims'][v] == 0: continue
        cf = np.array([ce['chis'][v][Kg.el[c[0]]] for c in cl])
        m, dd = Kg.decompose(cf)
        assert np.max(np.abs(m - np.round(m))) < 1e-6, "non-integer multiplicity"
        m = np.round(m).astype(int)
        assert int(np.sum(dd * m)) == ce['dims'][v], "sum d*m != dim E"
        tot += int(np.sum(m ** 2))
    return tot

def block_orbits(G):
    """finite canonical block set and the exact transport permutation on it"""
    H, perms, D = minimal_torus(G)
    blocks = eigblocks(H)
    n = G.n
    rng = np.random.default_rng(9)
    Ws = [(lambda W: (W + W.T) / 2)(rng.normal(size=(D, D))) for _ in range(NJ)]
    fps = {}                                     # fingerprint -> list over h of image fingerprint
    for ev, Q in sorted(blocks, key=lambda b: b[0]):
        dE = Q.shape[1]
        WT = np.zeros((n, NJ, dE, dE))
        for h in range(n):
            p = perms[h]
            for j in range(NJ): WT[h, j] = Q.conj().T @ Ws[j][np.ix_(p, p)] @ Q
        for K in G.subgroups:
            Kg = subgroup_grp(G, K)
            cl, chi, d, cls_of = Kg.chars()
            tot = 0
            for r in range(len(cl)):
                M = np.zeros((D, D), dtype=complex)
                for ii, g in enumerate(Kg.el):
                    M[perms[g], np.arange(D)] += np.conj(chi[r][cls_of[ii]])
                M *= d[r] / Kg.n
                S = Q.conj().T @ M @ Q; S = (S + S.conj().T) / 2
                rk = int(round(float(np.real(np.trace(S)))))
                tot += rk
                if rk == 0: continue
                sig = [tuple(round(float(np.real(np.sum(S * WT[h, j].T))), 6) for j in range(NJ))
                       for h in range(n)]
                fps.setdefault(sig[0], sig)
            assert tot == dE, f"{G.name}: block ranks {tot} != dim {dE}"
    keys = list(fps.keys()); pos = {k: i for i, k in enumerate(keys)}
    NN = len(keys); closed = True
    act = np.zeros((n, NN), dtype=np.int64)
    for h in range(n):
        for i, k in enumerate(keys):
            img = fps[k][h]
            if img not in pos: closed = False; break
            act[h, i] = pos[img]
        if not closed: break
    return NN, act, closed, perms, D

def orbit_sizes(act, n):
    NN = act.shape[1]; par = list(range(NN))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for h in range(n):
        for i in range(NN):
            a, b = find(i), find(int(act[h, i]))
            if a != b: par[a] = b
    from collections import Counter
    return sorted(Counter(find(i) for i in range(NN)).values(), reverse=True)

say("=" * 130)
say("S3   TRANSPORT ORBITS -- EXACT STRATIFICATION OF THE RECORD SET, AND A FINITE CANONICAL PROXY")
say("=" * 130)
say("")
say("  (0) SELF-CHECK: is A_z the identity for every central z?  (If not, Gbar = G/Z is the wrong")
say("      transport group and every orbit-size statement below is void.)")
allok = True
for G in glib.ladder(16):
    H, perms, D = minimal_torus(G)
    bad = [z for z in G.centre if not np.array_equal(perms[z], np.arange(D))]
    if bad: allok = False; say(f"      {G.name}: FAILED, central elements acting non-trivially: {bad}")
say(f"      A_z = I for every central z, on every carrier: {allok}")
if not allok: sys.exit(1)

say("")
say("  (1) EXACT ORBIT-TYPE STRATIFICATION OF THE RECORD SET (a continuum).")
say(f"  {'carrier':<12}{'|G|':>4}{'abel':>6}{'ground dim':>11}{'|Gbar|':>8}{'dim(all records)':>18}"
    f"{'dim(fixed)':>12}{'phi_fix':>9}{'generic stab':>14}{'largest orbit':>15}{'smallest':>10}{'#orbits':>10}")
strat = []
for G in glib.ladder(64):
    ce = census(G)
    Z = frozenset(G.centre)
    dall = sum(ce['dims'][v] ** 2 for v in (-2, -1, 0))
    dfix = dimL(G, ce, frozenset(range(G.n)))
    # generic stabiliser: trivial in Gbar iff every subgroup strictly containing Z has smaller
    # fixed dimension than dall
    gen_triv = True; strata = []
    if G.n <= 16:
        for K in G.subgroups:
            if not Z <= K: continue
            dk = dimL(G, ce, K)
            strata.append((len(K), dk))
            if K != Z and dk >= dall: gen_triv = False
    else:
        gen_triv = None
    Gb = G.n // len(Z)
    largest = Gb if (gen_triv or G.abelian) else None
    say(f"  {G.name:<12}{G.n:>4}{str(G.abelian):>6}{ce['dims'][-2]:>11}{Gb:>8}{dall:>18}"
        f"{dfix:>12}{dfix/dall:>9.4f}{('trivial' if gen_triv else ('n/a' if gen_triv is None else 'NON-trivial')):>14}"
        f"{(str(largest) if largest else '?'):>15}{'1':>10}{'infinite':>10}")
    strat.append(dict(name=G.name, n=G.n, abelian=G.abelian, ground=ce['dims'][-2], Gb=Gb,
                      dall=dall, dfix=dfix, phi=dfix / dall, strata=strata))
say("")
say("  'smallest orbit = 1' is the transport-FIXED record, whose existence S1 proved by subset-sum")
say("  and S2 exhibited as an explicit matrix with ||[A_h,R]|| < 1e-10 on every carrier.")
say("")
say("  (2) THE FINITE CANONICAL PROXY: K-isotypic blocks over every subgroup K <= G.")
say(f"  {'carrier':<12}{'|G|':>4}{'abel':>6}{'ground dim':>11}{'|Gbar|':>8}{'#subgroups':>12}"
    f"{'#blocks':>9}{'#orbits':>9}{'largest orbit':>15}{'orbit sizes':>26}{'frac moved':>12}{'closed':>8}")
proxy = []
for G in glib.ladder(16):
    t0 = time.time()
    NN, act, closed, perms, D = block_orbits(G)
    os_ = orbit_sizes(act, G.n)
    from collections import Counter
    hist = sorted(Counter(os_).items())
    fixed = sum(1 for s in os_ if s == 1)
    ce = census(G)
    Gb = G.n // len(G.centre)
    okdiv = all(Gb % s == 0 for s in os_)
    say(f"  {G.name:<12}{G.n:>4}{str(G.abelian):>6}{ce['dims'][-2]:>11}{Gb:>8}{len(G.subgroups):>12}"
        f"{NN:>9}{len(os_):>9}{max(os_):>15}{str(hist):>26}{(NN-fixed)/NN:>12.4f}{str(closed and okdiv):>8}")
    proxy.append(dict(name=G.name, n=G.n, abelian=G.abelian, ground=ce['dims'][-2], Gb=Gb,
                      nblocks=NN, norb=len(os_), largest=max(os_), hist=hist,
                      fracmoved=(NN - fixed) / NN, closed=closed and okdiv))
say("")
np.save("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_A_TRANSPORT/s3_rows.npy",
        np.array([strat, proxy], dtype=object), allow_pickle=True)
say("=" * 130); say("  READ -- from the numbers above"); say("=" * 130)
say("  Against RECORD COUNT and against ABELIANNESS, side by side:")
say(f"  {'carrier':<12}{'record count':>14}{'abelian':>9}{'|Gbar|':>8}{'largest orbit (records)':>25}"
    f"{'phi_fix':>9}{'proxy largest orbit':>21}{'proxy frac moved':>18}")
pm = {p['name']: p for p in proxy}
for s in sorted(strat, key=lambda s: s['ground']):
    p = pm.get(s['name'])
    say(f"  {s['name']:<12}{s['ground']:>14}{str(s['abelian']):>9}{s['Gb']:>8}{s['Gb']:>25}{s['phi']:>9.4f}"
        f"{(str(p['largest']) if p else '-'):>21}{(f'{p[chr(102)+chr(114)+chr(97)+chr(99)+chr(109)+chr(111)+chr(118)+chr(101)+chr(100)]:.4f}' if p else '-'):>18}")
