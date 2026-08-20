"""Carrier families for the k-threshold sweep.  SHARED MODULE -- no output of its own.

Two structurally different families that both scale the number of records k:

  A  [[n, n-2, 2]]   n even, stabilisers X^(x)n and Z^(x)n, k = n-2 records, dim 2^n.
                     ONE code block: every record lives on the same block.
  B  m blocks of [[4,2,2]] tensored, k = 2m records, dim 4^(2m)... = 2^(4m).
                     m DISJOINT blocks: records in different blocks share no qubit.

Both are stabiliser carriers so every relation among records is exact over F_2 -- no sampling.
Records are NEVER nominated: they come from symplectic_logicals(stab_xz, n) (a list of
CONJUGATE PAIRS) and the commuting family is one member of each pair.
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import symplectic_logicals, xz_to_matrix   # noqa: F401

# ------------------------------------------------------------------ F_2 Pauli helpers
def sp(a, b, n):
    """symplectic form: 1 iff the two Paulis ANTICOMMUTE"""
    return sum(a[i] * b[n + i] + a[n + i] * b[i] for i in range(n)) % 2

def support(v, n):
    return frozenset(i for i in range(n) if v[i] or v[n + i])

def weight(v, n):
    return len(support(v, n))

# ------------------------------------------------------------------ family A
def family_A(n):
    """[[n, n-2, 2]]: stabilisers X^(x)n and Z^(x)n.  n must be EVEN so they commute."""
    assert n % 2 == 0
    S = [[1] * n + [0] * n, [0] * n + [1] * n]
    return dict(name="A", label="[[%d,%d,2]]" % (n, n - 2), n=n, dim=2 ** n, stabs=S)

# ------------------------------------------------------------------ family B
def family_B(m):
    """m disjoint [[4,2,2]] blocks."""
    n = 4 * m
    S = []
    for b in range(m):
        x = [0] * n + [0] * n
        z = [0] * n + [0] * n
        for j in range(4 * b, 4 * b + 4):
            x[j] = 1
            z[n + j] = 1
        S.append(x); S.append(z)
    return dict(name="B", label="[[4,2,2]]^%d" % m, n=n, dim=2 ** n, stabs=S)

# ------------------------------------------------------------------ records from the carrier
def records_of(car, verbose=False):
    """The COMMUTING FAMILY of records: one member of each conjugate pair from
       symplectic_logicals.  Returns (family, pairs) and runs the mandatory self-checks."""
    n = car["n"]
    pairs = symplectic_logicals(car["stabs"], n)
    k = len(pairs)
    fam = [a for a, b in pairs]
    part = [b for a, b in pairs]
    checks = {}
    # SELF-CHECK 1: the symplectic pairing matrix must be the IDENTITY (non-degenerate).
    M = [[sp(fam[i], part[j], n) for j in range(k)] for i in range(k)]
    checks["pairing_is_identity"] = all(M[i][j] == (1 if i == j else 0)
                                        for i in range(k) for j in range(k))
    # SELF-CHECK 2: the chosen family must be MUTUALLY COMMUTING.
    checks["family_commutes"] = all(sp(fam[i], fam[j], n) == 0
                                    for i in range(k) for j in range(i + 1, k))
    # SELF-CHECK 3: every record must COMMUTE WITH EVERY STABILISER (clause ii, symbolic).
    checks["in_normaliser"] = all(sp(r, s, n) == 0 for r in fam + part for s in car["stabs"])
    # SELF-CHECK 4: no record may BE a stabiliser (clause iii would fail).  Tested by RANK,
    # not by enumerating the stabiliser group -- 2^(rank S) is 2^32 at the top of family B.
    from f2 import rank as _rank
    r0 = _rank(list(car["stabs"]), 2 * n)
    checks["not_a_stabiliser"] = all(_rank(list(car["stabs"]) + [r], 2 * n) > r0 for r in fam)
    # SELF-CHECK 5: k must equal n - rank(S)
    checks["k_equals_n_minus_rank"] = (k == n - len(car["stabs"]))
    return fam, part, checks

def all_checks_pass(checks):
    return all(checks.values())

# ------------------------------------------------------------------ representative gauge
def _stab_groups(car):
    """Partition the stabiliser generators into groups with disjoint joint support, so a
       coset minimisation factorises.  For A there is ONE group; for B one per block."""
    n = car["n"]
    gens = car["stabs"]
    groups, sups = [], []
    for s in gens:
        su = support(s, n)
        hit = [i for i, t in enumerate(sups) if t & su]
        if not hit:
            groups.append([s]); sups.append(set(su))
        else:
            g = [s]; t = set(su)
            for i in reversed(hit):
                g += groups.pop(i); t |= sups.pop(i)
            groups.append(g); sups.append(t)
    return groups

def min_weight_rep(v, car):
    """Least-weight representative of v modulo the stabiliser group (exact, factorised)."""
    n = car["n"]
    out = v[:]
    for g in _stab_groups(car):
        best, bw = out, weight(out, n)
        for bits in range(1, 2 ** len(g)):
            w = out[:]
            for i, s in enumerate(g):
                if (bits >> i) & 1:
                    w = [(x + y) % 2 for x, y in zip(w, s)]
            ww = weight(w, n)
            if ww < bw:
                best, bw = w, ww
        out = best
    return out

def random_rep(v, car, rng):
    n = car["n"]
    w = v[:]
    for s in car["stabs"]:
        if rng.integers(0, 2):
            w = [(x + y) % 2 for x, y in zip(w, s)]
    return w

# ------------------------------------------------------------------ family C
def family_C(m):
    """m disjoint [[5,1,3]] perfect-code blocks.  k = m records, n = 5m, DISTANCE 3.

       Structurally different from A and B in the two ways that matter to a threshold claim:
         * record density k/n = 1/5   (A -> 1, B -> 1/2)
         * protection radius 3        (A and B -> 2)
       Stabilisers XZZXI, IXZZX, XIXZZ, ZXIXZ per block."""
    n = 5 * m
    pat = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]
    S = []
    for b in range(m):
        for p in pat:
            v = [0] * (2 * n)
            for j, ch in enumerate(p):
                q = 5 * b + j
                if ch in "XY": v[q] = 1
                if ch in "ZY": v[n + q] = 1
            S.append(v)
    return dict(name="C", label="[[5,1,3]]^%d" % m, n=n, dim=2 ** n, stabs=S)
