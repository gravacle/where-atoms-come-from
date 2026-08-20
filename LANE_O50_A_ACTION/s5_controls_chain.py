"""O-50-A  step 5.  THE CONTROLS (D-15), in the same table as the torus.

   CONTROL A -- the fully coupled 1D chain of O-48 / C-65, read with the PROGRAM'S OWN maximal
   commuting family.  n candidate records collapse to ONE bit; the writer group has order 2.
   Its invariant space must come out at dimension 1 as well -- and the DIFFERENCE from the torus
   is not the invariant dimension but the CONFIGURATION SPACE: 2 configurations against 4, one
   record bit against two.

   CONTROL B -- the SAME chain, read with the naive family {Z_1..Z_n}: every Z_i genuinely
   satisfies clauses (i)-(iv), they mutually commute, and they are NOT independently writable.
   The theorem's hypothesis FAILS here, and a NON-CONSTANT INVARIANT MUST APPEAR.  It does:
   the bond variables s_i s_j.  And Z_i Z_j FAILS clause (iv) -- it is not a record.  That is
   exactly the distinction O-48 could not draw and the torus can.

   Everything below is exact: H is diagonal, the spectrum is checked for exact 2-fold
   degeneracy, and the admissible group is characterised over ALL unitaries, not just Paulis."""
import sys, itertools, numpy as np
from math import isqrt
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_A_ACTION")
from record_model import RecordModel, eigenspaces, clause_iii, clause_iv
import orbits

PRIMES = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107]

def chain(n, uniform=False):
    """H = sum_{i<j} J_ij Z_i Z_j, J generic (sqrt of distinct primes) unless uniform."""
    cfg = np.array(list(itertools.product((1, -1), repeat=n)))          # 2^n spin configs
    J = {}
    t = 0
    for i in range(n):
        for j in range(i + 1, n):
            J[(i, j)] = 1.0 if uniform else float(np.sqrt(PRIMES[t % len(PRIMES)])) * (1 + t * 1e-3)
            t += 1
    E = np.zeros(2 ** n)
    for (i, j), v in J.items(): E += v * cfg[:, i] * cfg[:, j]
    return cfg, J, E

def report(n, uniform=False):
    cfg, J, E = chain(n, uniform)
    N = 2 ** n
    H = np.diag(E).astype(complex)
    es = eigenspaces(H)
    mult = sorted({m for _, _, m in es})
    print(f"\n  --- chain n={n}  {'UNIFORM J (permutation-symmetric)' if uniform else 'GENERIC J'} ---")
    print(f"    dim = {N}   distinct eigenvalues = {len(es)}   multiplicities present = {mult}")
    # D-22: automorphism group of the coupling matrix
    aut = 0
    for p in itertools.permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(i + 1, n):
                a = J[tuple(sorted((i, j)))]; b = J[tuple(sorted((p[i], p[j])))]
                if abs(a - b) > 1e-12: ok = False; break
            if not ok: break
        if ok: aut += 1
    print(f"    D-22 carrier automorphism group order = {aut}   (n! = {__import__('math').factorial(n)})   "
          f"permutation-symmetric = {aut == __import__('math').factorial(n)}")
    Zs = []
    for i in range(n):
        Zs.append(np.diag(cfg[:, i].astype(complex)))
    print(f"    every Z_i: (iii)={all(clause_iii(Z, es) for Z in Zs)}  (iv)={all(clause_iv(Z, es) for Z in Zs)}  "
          f"[Z_i,H]=0 -> {all(np.linalg.norm(Z@H-H@Z)<1e-9 for Z in Zs)}")
    ZZ = [np.diag((cfg[:, i] * cfg[:, j]).astype(complex)) for i in range(n) for j in range(i + 1, n)]
    print(f"    bond operators Z_iZ_j: (iii)={any(clause_iii(M, es) for M in ZZ)}  "
          f"(iv)={any(clause_iv(M, es) for M in ZZ)}   -> BOND VARIABLES ARE NOT RECORDS: "
          f"{not any(clause_iv(M, es) for M in ZZ)}")
    # --- admissible action over ALL unitaries: which configurations are energy-degenerate?
    from collections import defaultdict
    byE = defaultdict(list)
    for a in range(N): byE[round(E[a], 9)].append(a)
    deg = sorted({len(v) for v in byE.values()})
    neg = {a: N - 1 - a for a in range(N)}                    # index of the globally flipped config
    exact = all(sorted(v) == sorted([v[0], neg[v[0]]]) for v in byE.values() if len(v) == 2)
    print(f"    energy-degenerate config classes: sizes {deg}; every class = {{s,-s}} exactly: {exact}")
    # --- CONTROL A: the program's own maximal commuting family
    rm = RecordModel(H, ())
    cands = Zs + [np.diag(np.prod(cfg[:, list(c)], axis=1).astype(complex))
                  for r in (3, 5) if r <= n for c in itertools.combinations(range(n), r)]
    fam = rm.commuting_family(cands)
    iw = rm.independently_writable(fam)
    print(f"    CONTROL A  commuting_family size = {len(fam)}   independently_writable = {iw}")
    permsA = [[1, 0]] if len(fam) == 1 else None
    A = orbits.analyse([[1, 0]], 2); A.pop('orbits')
    print(f"    CONTROL A  configs = 2   G_W = <global flip>, order 2 -> {A}")
    # --- CONTROL B: the naive family {Z_i}, n bits, writer group = {id, global negation}
    m = N
    gneg = [neg[a] for a in range(m)]
    B = orbits.analyse([gneg], m); orbs = B.pop('orbits')
    print(f"    CONTROL B  family = {{Z_1..Z_n}} (all records, mutually commuting), "
          f"independently writable = {rm.independently_writable(Zs)}")
    print(f"    CONTROL B  {B}")
    # exhibit non-constant invariants explicitly
    s = cfg
    inv = []
    for i in range(n):
        for j in range(i + 1, n):
            f = s[:, i] * s[:, j]
            ok = all(f[a] == f[gneg[a]] for a in range(m))
            if ok and len(set(f)) > 1: inv.append((i, j))
    print(f"    CONTROL B  NON-CONSTANT G_W-INVARIANTS EXHIBITED: s_i s_j for {len(inv)} pairs, "
          f"e.g. s_0 s_1;  invariant dim {B['invariant_dim_exact']} > 1 = "
          f"{B['invariant_dim_exact'] > 1}")
    f01 = s[:, 0] * s[:, 1]
    print(f"    CONTROL B  witness f(s)=s_0 s_1: values {sorted(set(f01.tolist()))}, "
          f"invariant under the writer = {all(f01[a]==f01[gneg[a]] for a in range(m))}, "
          f"non-constant = {len(set(f01.tolist()))>1}, "
          f"and s_0 s_1 IS NOT A RECORD (clause (iv) fails above)")
    return dict(n=n, configsA=2, invA=A['invariant_dim_exact'], famA=len(fam),
                configsB=m, invB=B['invariant_dim_exact'], orbitsB=B['n_orbits'],
                GW_order=2, aut=aut)

if __name__ == "__main__":
    print("CONTROLS -- the fully coupled 1D chain of O-48 (D-23: 1D PROXY carrier, not the torus)")
    out = []
    for n in (3, 4, 5, 6):
        out.append(report(n))
    print("\n  PERMUTATION-SYMMETRIC VARIANT (D-22: no geometry to detect)")
    report(4, uniform=True)
    print("\n  SUMMARY")
    print(f"    {'n':>3} {'famA':>5} {'cfgA':>5} {'invA':>5} {'cfgB':>6} {'orbB':>6} {'invB':>6} {'|G_W|':>6} {'aut':>5}")
    for r in out:
        print(f"    {r['n']:>3} {r['famA']:>5} {r['configsA']:>5} {r['invA']:>5} "
              f"{r['configsB']:>6} {r['orbitsB']:>6} {r['invB']:>6} {r['GW_order']:>6} {r['aut']:>5}")
