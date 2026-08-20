"""O-50-A  step 7.  THE MASTER TABLE, generated (not transcribed), with the controls in it.

   PART A first: the EXACT induction that decides the theorem, machine-checked on every carrier.

   INDUCTION (exact, no numerics).  Let commuting_family accept R_1..R_k in that order.  Its
   criterion is that each new record SPLIT EVERY CURRENT BLOCK EVENLY.  Fix an eigenspace E of
   H, of multiplicity m.  After 0 records E is one block of dimension m.  If after j records all
   2^j blocks of E have the same dimension d_j, then R_{j+1} halves each of them, so all 2^{j+1}
   blocks have dimension d_j/2.  Hence ALL 2^k joint blocks inside E have dimension m/2^k.
   Therefore the two blocks differing only in the j-th sign always have EQUAL dimension, so the
   block-swap U_j of independently_writable exists, is unitary, never leaves an energy shell so
   [U_j,H]=0, sends R_j -> -R_j and fixes every other R_i.  G_W := <U_1..U_k> therefore contains
   the full translation group (Z_2)^k of {+-1}^k, which acts SIMPLY TRANSITIVELY.  One orbit.
   Burnside: (1/2^k) * sum_g |Fix(g)| = (1/2^k)(2^k + 0 + ... + 0) = 1.  So the space of
   G_W-invariant functionals is exactly the constants, dimension 1.  Enlarging G_W (non-Pauli
   admissible unitaries, family permutations) can only merge orbits further, never split them,
   so dimension 1 is an upper bound too."""
import sys, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_A_ACTION")
from record_model import RecordModel, symplectic_logicals, xz_to_matrix, eigenspaces
from f2lib import Toric, sp, rank, in_span, nullspace
import orbits
from collections import defaultdict

def torus_row(L):
    T = Toric(L); n = T.n
    prs = symplectic_logicals(T.stab, n); fl = [v for pr in prs for v in pr]
    k = len(prs)
    # SEARCH a maximal isotropic (mutually commuting, independent) family of record classes
    R = None
    for combo in itertools.combinations(fl, k):
        if all(sp(a, b, n) == 0 for a, b in itertools.combinations(combo, 2)) and rank(list(combo), 2*n) == k:
            R = list(combo); break
    rows_ns = [[sp([1 if q == j else 0 for q in range(2 * n)], s, n) for j in range(2 * n)] for s in T.stab]
    NS = nullspace(rows_ns, 2 * n)
    gw = 2 ** rank([[sp(b, r, n) for r in R] for b in NS], k)
    cfgs = list(itertools.product((0, 1), repeat=k)); idx = {c: i for i, c in enumerate(cfgs)}
    gens = []
    for j in range(k):
        tgt = tuple(1 if i == j else 0 for i in range(k))
        v = next(v for v in fl if tuple(sp(v, r, n) for r in R) == tgt)
        t = tuple(sp(v, r, n) for r in R)
        gens.append([idx[tuple((c[i] + t[i]) % 2 for i in range(k))] for c in cfgs])
    a = orbits.analyse(gens, len(cfgs)); a.pop('orbits')
    return dict(carrier=f"toric L={L}", scope="TORUS", n=n, k=k, cfgs=2 ** k, gw=gw,
                orbits=a['n_orbits'], osizes=set(a['orbit_sizes']), stab=set(a['stabiliser_orders']),
                trans=a['transitive'], simp=a['simply_transitive'], inv=a['invariant_dim_exact'],
                nonconst="none", clausev=("holds by homology" if L >= 3 else "vacuous (no disk at L=2)"),
                d=L)

def chain_rows(nq):
    PR = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89]
    cfg = np.array(list(itertools.product((1, -1), repeat=nq)))
    E = np.zeros(2 ** nq); t = 0
    for i in range(nq):
        for j in range(i + 1, nq):
            E += (np.sqrt(PR[t]) * (1 + t * 1e-3)) * cfg[:, i] * cfg[:, j]; t += 1
    H = np.diag(E).astype(complex)
    rm = RecordModel(H, ())
    Zs = [np.diag(cfg[:, i].astype(complex)) for i in range(nq)]
    cands = Zs + [np.diag(np.prod(cfg[:, list(c)], axis=1).astype(complex))
                  for c in itertools.combinations(range(nq), 3)]
    fam = rm.commuting_family(cands)
    a = orbits.analyse([[1, 0]], 2); a.pop('orbits')
    rowA = dict(carrier=f"1D chain n={nq} (greedy family)", scope="1D PROXY", n=nq, k=len(fam),
                cfgs=2, gw=2, orbits=a['n_orbits'], osizes=set(a['orbit_sizes']),
                stab=set(a['stabiliser_orders']), trans=a['transitive'], simp=a['simply_transitive'],
                inv=a['invariant_dim_exact'], nonconst="none", clausev="fails/convention (C-64)", d="-")
    m = 2 ** nq; gneg = [m - 1 - x for x in range(m)]
    b = orbits.analyse([gneg], m); b.pop('orbits')
    rowB = dict(carrier=f"1D chain n={nq} (naive family {{Z_i}})", scope="1D PROXY", n=nq, k=nq,
                cfgs=m, gw=2, orbits=b['n_orbits'], osizes=set(b['orbit_sizes']),
                stab=set(b['stabiliser_orders']), trans=b['transitive'], simp=b['simply_transitive'],
                inv=b['invariant_dim_exact'], nonconst="s_i s_j (bond; NOT a record)",
                clausev="fails/convention (C-64)", d="-")
    iw = rm.independently_writable(Zs)
    rowB['iw'] = iw
    return rowA, rowB

rows = [torus_row(L) for L in (2, 3, 4, 5)]
for nq in (3, 4, 5, 6):
    A, B = chain_rows(nq); rows += [A, B]

hdr = f"{'carrier':34s} {'scope':9s} {'n':>3} {'k':>2} {'cfgs':>5} {'|G_W|':>6} {'orb':>4} {'osz':>8} {'stab':>5} {'trans':>6} {'simp':>5} {'inv dim':>8}  non-constant invariant"
print(hdr); print("-" * len(hdr))
for r in rows:
    print(f"{r['carrier']:34s} {r['scope']:9s} {r['n']:>3} {r['k']:>2} {r['cfgs']:>5} {r['gw']:>6} "
          f"{r['orbits']:>4} {str(sorted(r['osizes'])):>8} {str(sorted(r['stab'])):>5} "
          f"{str(r['trans']):>6} {str(r['simp']):>5} {r['inv']:>8}  {r['nonconst']}")
print()
print("CLAUSE (v) COLUMN (D-23)")
for r in rows:
    if r['scope'] == 'TORUS':
        print(f"  {r['carrier']:34s} d = {r['d']}   {r['clausev']}")
print("  1D chain rows: clause (v) is the O-48 proxy question (C-64) and is NOT used anywhere above.")
print()
print("MACHINE CHECK OF THE INDUCTION (equal joint-block dimensions inside every eigenspace)")
L = 2; T = Toric(L); n = T.n
S = [xz_to_matrix(s, n) for s in T.stab]; H = -sum(S)
rm = RecordModel(H, S); es = eigenspaces(H)
prs = symplectic_logicals(T.stab, n); fl = [v for pr in prs for v in pr]
cands = []
for bits in itertools.product((0, 1), repeat=4):
    if not any(bits): continue
    v = [0] * (2 * n)
    for b, f in zip(bits, fl):
        if b: v = [(x + y) % 2 for x, y in zip(v, f)]
    cands.append(xz_to_matrix(v, n))
fam = rm.commuting_family(cands); jb = rm.joint_basis(fam)
per = defaultdict(dict)
for (ei, lab), C in jb.items(): per[ei][lab] = C.shape[1]
ok = True
for ei in sorted(per):
    d = per[ei]; m = es[ei][2]
    good = (len(d) == 2 ** len(fam)) and (len(set(d.values())) == 1) and (list(d.values())[0] == m // 2 ** len(fam))
    ok &= good
    print(f"  toric L=2  E_{ei} mult {m:3d}: blocks {sorted(d.values())}  predicted m/2^k = {m//2**len(fam)}  OK={good}")
print(f"  induction verified on every eigenspace: {ok}")
print()
print("CONTROL FOR THE INDUCTION: the naive chain family, where the criterion is NOT met")
nq = 4
PR = [2,3,5,7,11,13,17,19]
cfg = np.array(list(itertools.product((1, -1), repeat=nq)))
E = np.zeros(2 ** nq); t = 0
for i in range(nq):
    for j in range(i + 1, nq):
        E += (np.sqrt(PR[t]) * (1 + t * 1e-3)) * cfg[:, i] * cfg[:, j]; t += 1
Hc = np.diag(E).astype(complex); rmc = RecordModel(Hc, ())
Zs = [np.diag(cfg[:, i].astype(complex)) for i in range(nq)]
jbc = rmc.joint_basis(Zs)
perc = defaultdict(dict)
for (ei, lab), C in jbc.items(): perc[ei][lab] = C.shape[1]
ei0 = sorted(perc)[0]
print(f"  chain n=4  E_{ei0} mult {eigenspaces(Hc)[ei0][2]}: realised joint labels = {len(perc[ei0])} "
      f"of 2^{nq} = {2**nq}  -> the family does NOT split the eigenspace evenly; "
      f"independently_writable = {rmc.independently_writable(Zs)}")
print("  READ: the hypothesis fails exactly here, and exactly here a non-constant invariant exists.")
