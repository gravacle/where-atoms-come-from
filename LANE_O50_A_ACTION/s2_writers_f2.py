"""O-50-A  step 2.  SEARCH for the admissible writer group on the torus, and characterise its
   action on record configurations.  Exact F_2 / integer arithmetic only.

   ADMISSIBILITY, exactly:  for a Pauli P, P H P^dag = -sum_g eps_g g with eps_g = +-1, so
   [P,H] = 0  iff  sum_g (1-eps_g) g = 0  iff  eps_g = +1 for every stabiliser generator g,
   because DISTINCT PAULI MATRICES ARE LINEARLY INDEPENDENT.  (Distinctness is verified below.)
   Hence {admissible Paulis} = N(S) EXACTLY -- a computed set, not a nominated one.

   ACTION:  P^dag R P = (-1)^{sp(P,R)} R, so P translates the record configuration by the
   vector phi(P) = (sp(P,R_1),...,sp(P,R_k)) in F_2^k.  phi is computed, never assumed."""
import sys, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_A_ACTION")
from record_model import symplectic_logicals
from f2lib import Toric, sp, rref, rank, in_span, nullspace, span
import orbits

def logical_classes(T, pairs):
    """all 2^{2k} logical classes as F_2 combinations of the computed conjugate pairs."""
    n = T.n; flat = [v for pr in pairs for v in pr]
    out = []
    for bits in itertools.product((0, 1), repeat=len(flat)):
        v = [0] * (2 * n)
        for b, f in zip(bits, flat):
            if b: v = [(x + y) % 2 for x, y in zip(v, f)]
        out.append((bits, v))
    return out

def maximal_commuting_families(T, pairs):
    """SEARCH: every maximal isotropic subspace of the logical symplectic space, i.e. every
       maximal set of mutually commuting, S-independent logical classes."""
    n = T.n; k = len(pairs)
    cls = logical_classes(T, pairs)
    nz = [(b, v) for b, v in cls if any(b)]
    fams = []
    for combo in itertools.combinations(nz, k):
        bits = [b for b, _ in combo]
        if rank([list(b) for b in bits], 2 * k) != k: continue       # independent classes
        if all(sp(a, b, n) == 0 for (_, a), (_, b) in itertools.combinations(combo, 2)):
            fams.append(combo)
    return fams

def analyse(L, exhaustive_pauli_scan=False):
    T = Toric(L); n = T.n
    pairs = symplectic_logicals(T.stab, n)
    k = len(pairs)
    rec = dict(L=L, n=n, k=k)

    # -- distinctness of the stabiliser generators as Paulis (needed for the admissibility proof)
    rec['stabiliser_generators_distinct'] = (len({tuple(s) for s in T.stab}) == len(T.stab))

    # -- SEARCH for maximal commuting families of records
    fams = maximal_commuting_families(T, pairs)
    rec['n_maximal_commuting_families_of_logicals'] = len(fams)
    rec['max_commuting_family_size'] = k
    fam = fams[0]
    R = [v for _, v in fam]
    rec['chosen_family_class_bits'] = [list(b) for b, _ in fam]
    rec['chosen_family_weights'] = [sum(1 for i in range(n) if v[i] or v[n + i]) for v in R]
    rec['chosen_family_mutually_commute'] = all(sp(a, b, n) == 0 for a, b in itertools.combinations(R, 2))

    # -- admissible Paulis = N(S): computed as the symplectic complement (nullspace)
    rows = [[sp([1 if q == j else 0 for q in range(2 * n)], s, n) for j in range(2 * n)] for s in T.stab]
    NS = nullspace(rows, 2 * n)
    rec['dim_N(S)'] = len(NS)
    rec['n_admissible_pauli_classes'] = 2 ** len(NS)
    rec['n_all_pauli_classes'] = 2 ** (2 * n)

    # -- exhaustive brute-force scan over the WHOLE Pauli group (only feasible at L=2)
    if exhaustive_pauli_scan:
        adm = 0; hits = {}
        for bits in itertools.product((0, 1), repeat=2 * n):
            v = list(bits)
            if any(sp(v, s, n) for s in T.stab): continue
            adm += 1
            f = tuple(sp(v, r, n) for r in R)
            hits.setdefault(f, []).append(v)
        rec['brute_force_admissible_count'] = adm
        rec['brute_force_matches_N(S)'] = (adm == 2 ** len(NS))
        rec['brute_force_phi_image'] = sorted(hits.keys())
        rec['brute_force_phi_fibre_sizes'] = {str(f): len(vs) for f, vs in sorted(hits.items())}
        wgen = {}
        for f, vs in hits.items():
            if any(f): wgen[f] = min(vs, key=lambda v: sum(1 for i in range(n) if v[i] or v[n + i]))
        rec['brute_force_min_weight_writer_per_flip'] = {
            str(f): dict(weight=sum(1 for i in range(n) if v[i] or v[n + i]),
                         support=[i for i in range(n) if v[i] or v[n + i]])
            for f, v in sorted(wgen.items())}
    # -- image of phi, by exact linear algebra over F_2 (all L)
    phi_rows = [[sp(b, r, n) for r in R] for b in NS]
    rec['phi_image_dim'] = rank(phi_rows, k)
    rec['phi_image_order'] = 2 ** rank(phi_rows, k)
    rec['phi_kernel_order'] = 2 ** (len(NS) - rank(phi_rows, k))
    rec['writer_group_G_W_order'] = 2 ** rank(phi_rows, k)

    # -- explicit generators: SEARCH N(S) for minimum-weight elements realising each basis flip
    gens = {}
    if len(NS) <= 12:
        pool = span(NS, 2 * n)
    else:
        pool = None
    for j in range(k):
        target = tuple(1 if i == j else 0 for i in range(k))
        best = None
        if pool is not None:
            for v in pool:
                if tuple(sp(v, r, n) for r in R) == target:
                    w = sum(1 for i in range(n) if v[i] or v[n + i])
                    if best is None or w < best[0]: best = (w, v)
        else:                                        # L>=4: search the logical coset directly
            for bits, v in logical_classes(T, pairs):
                if tuple(sp(v, r, n) for r in R) == target:
                    w = sum(1 for i in range(n) if v[i] or v[n + i])
                    if best is None or w < best[0]: best = (w, v)
        gens[j] = best
    rec['generators'] = {}
    for j, (w, v) in gens.items():
        rec['generators'][f'U_{j+1}'] = dict(
            weight=w,
            support=[i for i in range(n) if v[i] or v[n + i]],
            commutes_with_every_stabiliser=all(sp(v, s, n) == 0 for s in T.stab),
            in_S=in_span(v, T.stab, 2 * n),
            flips=[i + 1 for i in range(k) if sp(v, R[i], n) == 1],
            fixes=[i + 1 for i in range(k) if sp(v, R[i], n) == 0],
            x_class=T.x_class(v[:n]) if any(v[:n]) else None,
            z_class=T.z_class(v[n:]) if any(v[n:]) else None)
    rec['generators_flip_exactly_one_record_each'] = all(
        len(d['flips']) == 1 for d in rec['generators'].values())

    # -- the ACTION on record configurations, computed as explicit permutations
    cfgs = list(itertools.product((0, 1), repeat=k))
    idx = {c: i for i, c in enumerate(cfgs)}
    perms = []
    for j, (w, v) in gens.items():
        t = tuple(sp(v, r, n) for r in R)
        perms.append([idx[tuple((c[i] + t[i]) % 2 for i in range(k))] for c in cfgs])
    rec['action'] = orbits.analyse(perms, len(cfgs))
    rec['action'].pop('orbits')
    return rec

if __name__ == "__main__":
    for L in (2, 3, 4, 5):
        r = analyse(L, exhaustive_pauli_scan=(L == 2))
        print("=" * 78)
        for kk, vv in r.items():
            if kk == 'generators':
                print("  generators:")
                for g, d in vv.items(): print(f"      {g}: {d}")
            else:
                print(f"  {kk:44s} {vv}")
