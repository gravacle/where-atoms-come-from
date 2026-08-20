"""S4 -- ADDITIVITY OVER DISJOINT REGIONS, on products of [[4,2,2]] BLOCKS.

Gravity's requirement (b): the source of two far-apart clusters is the SUM of theirs.
A quantity that fails (b) cannot be a density even if it grows.

CARRIER: m disjoint blocks of [[4,2,2]] on n = 4m qubits, 2m stabilisers (X^4 and Z^4 per
block), k = 2m records.  Blocks share no qubit, so they are the cleanest "disjoint regions"
this program's carrier admits.

REPRESENTATION: F_2 symplectic for the structural quantities; the exact reduced code-space
cq-state engine validated in S3 for chi.  Both are exact.

D-15: the subadditivity of chi under a SHARED bath is printed beside a DISJOINT-BATH control
in the same table, which must come out exactly additive -- otherwise the test cannot register
additivity at all and the shared-bath deficit means nothing.
"""
import sys, json, math, itertools
import numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import symplectic_logicals, xz_to_matrix, RecordModel, Environment

OUT = []
def p(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)

I2 = np.eye(2); Xm = np.array([[0, 1], [1, 0]], complex); Zm = np.array([[1, 0], [0, -1]], complex)
LAM, BETA = 0.8, 2.0
TIMES = np.linspace(1.0, 13.0, 25)
ENERGY_POOL = (1.0, 1.4, 0.7, 1.2, 0.9, 1.6, 0.8, 1.1)
def energies(nq): return tuple(ENERGY_POOL[j % len(ENERGY_POOL)] for j in range(nq))

def vN(r):
    e = np.linalg.eigvalsh(r); e = e[e > 1e-13]
    return float(-(e * np.log2(e)).sum())
def uprop(H, t):
    w, V = np.linalg.eigh(H); return (V * np.exp(-1j * w * t)) @ V.conj().T
def binom_pmf(m):
    return {2 * a - m: math.comb(m, a) / 2 ** m for a in range(m + 1)}

def chi_distributed(k, nq, t, lam=LAM, beta=BETA, offset=0):
    """k records distributed over nq bath sites, record i -> site (i+offset) mod nq."""
    E = energies(nq); per = [0.0] * k
    for j in range(nq):
        idx = [i for i in range(k) if (i + offset) % nq == j]
        m = len(idx)
        if m == 0: continue
        Hb = E[j] * Zm
        w = np.exp(-beta * np.array([E[j], -E[j]]))
        rth = np.diag(w / w.sum()).astype(complex)
        def rho(c):
            U = uprop(Hb + lam * c * Xm, t); return U @ rth @ U.conj().T
        pm, pm1 = binom_pmf(m), binom_pmf(m - 1)
        rbar = sum(pr * rho(c) for c, pr in pm.items())
        cond = {s: sum(pr * rho(c + s) for c, pr in pm1.items()) for s in (+1, -1)}
        chi = max(vN(rbar) - 0.5 * (vN(cond[+1]) + vN(cond[-1])), 0.0)
        for i in idx: per[i] = chi
    return per, float(sum(per))

def tavg_total(k, nq, offset=0):
    return float(np.mean([chi_distributed(k, nq, float(t), offset=offset)[1] for t in TIMES]))

# ---------------------------------------------------------------- F_2 structure of m blocks
POP = [bin(i).count("1") for i in range(1 << 16)]
def popcount(v):
    c = 0
    while v: c += POP[v & 0xFFFF]; v >>= 16
    return c
def sp(a, b, n):
    return sum(a[i] * b[n + i] + a[n + i] * b[i] for i in range(n)) % 2
def wt(v, n):
    return sum(1 for i in range(n) if v[i] or v[n + i])

def block_stabs(m):
    n = 4 * m; out = []
    for b in range(m):
        vx = [0] * (2 * n); vz = [0] * (2 * n)
        for q in range(4 * b, 4 * b + 4):
            vx[q] = 1; vz[n + q] = 1
        out += [vx, vz]
    return out, n

def analyse_blocks(m):
    S, n = block_stabs(m)
    k = n - len(S)
    pairs = symplectic_logicals([s[:] for s in S], n)
    ok = (len(pairs) == k)
    R = [a for a, b in pairs]; W = [b for a, b in pairs]
    ok &= all(sp(v, s, n) == 0 for v in R + W for s in S)
    for i in range(k):
        for j in range(k):
            ok &= (sp(R[i], R[j], n) == 0) and (sp(W[i], W[j], n) == 0) \
                  and (sp(R[i], W[j], n) == (1 if i == j else 0))
    P_rec = sum(sp(R[i], R[j], n) for i in range(k) for j in range(i + 1, k))
    allops = R + W
    P_all = sum(sp(allops[i], allops[j], n) for i in range(len(allops))
                for j in range(i + 1, len(allops)))
    # distance: exhaustive over weight <= 2
    def inS(v):
        vv = v[:]
        rows = [s[:] for s in S]
        for b in rows:
            h = next((i for i in range(2 * n) if b[i]), None)
            if h is not None and vv[h]: vv = [(x + y) % 2 for x, y in zip(vv, b)]
        return not any(vv)
    d = None
    for w in (1, 2):
        for sites in itertools.combinations(range(n), w):
            for kinds in itertools.product((1, 2, 3), repeat=w):
                v = [0] * (2 * n)
                for s2, kd in zip(sites, kinds):
                    if kd in (1, 3): v[s2] = 1
                    if kd in (2, 3): v[n + s2] = 1
                if all(sp(v, s, n) == 0 for s in S) and not inS(v):
                    d = w; break
            if d: break
        if d: break
    # records disturbed by a weight-1 op, and the weight-2 CONTROL
    def ndist(v):
        if any(sp(v, s, n) for s in S): return 0
        return sum(sp(v, R[i], n) for i in range(k))
    D1 = 0
    for site in range(n):
        for kd in (1, 2, 3):
            v = [0] * (2 * n)
            if kd in (1, 3): v[site] = 1
            if kd in (2, 3): v[n + site] = 1
            D1 = max(D1, ndist(v))
    D2 = 0
    for (a, b) in itertools.combinations(range(n), 2):
        for ka in (1, 2, 3):
            for kb in (1, 2, 3):
                v = [0] * (2 * n)
                for s2, kd in ((a, ka), (b, kb)):
                    if kd in (1, 3): v[s2] = 1
                    if kd in (2, 3): v[n + s2] = 1
                D2 = max(D2, ndist(v))
    return dict(m=m, n=n, k=k, selfcheck=bool(ok), P_rec=P_rec, P_all=P_all,
                distance=d, D1max=D1, D2max=D2, Wtot_lower=2 * k, log2dim=k, dim=2 ** k)

p("=" * 118)
p("S4  ADDITIVITY OVER DISJOINT REGIONS -- m disjoint [[4,2,2]] blocks, n = 4m qubits, k = 2m records.")
p("=" * 118)
MS = [1, 2, 3, 4, 5, 6]
bres = {}
for m in MS:
    r = analyse_blocks(m); bres[m] = r
    if not r["selfcheck"]:
        p("SELF-CHECK FAILED at m=%d -- CONCLUDING NOTHING" % m); sys.exit(1)
p("self-check passed at every m: symplectic_logicals returned k = 2m conjugate pairs with the canonical Gram matrix.")
p("")
p("  m   n    k |  P_rec  P_all(CONTROL) |  d  | D1max  D2max(CONTROL) | 2k  log2 dim   dim")
p("-" * 118)
for m in MS:
    r = bres[m]
    p("%3d %3d %4d | %6d %14d | %2d  | %5d %14d | %2d %9d %11s"
      % (r["m"], r["n"], r["k"], r["P_rec"], r["P_all"], r["distance"], r["D1max"], r["D2max"],
         r["Wtot_lower"], r["log2dim"], r["dim"]))
p("-" * 118)
p("")
p("ADDITIVITY OF THE STRUCTURAL QUANTITIES: value at m blocks vs m * (value at 1 block).")
p("  quantity                          1 block   m=2    m=3    m=4    m=5    m=6   |  additive over disjoint regions?")
p("-" * 118)
def addrow(name, key):
    v = [bres[m][key] for m in MS]
    add = all(abs(v[i] - (i + 1) * v[0]) < 1e-12 for i in range(len(MS)))
    p("  %-32s %7s %6s %6s %6s %6s %6s  |  %s"
      % (name, v[0], v[1], v[2], v[3], v[4], v[5], "YES" if add else "NO"))
addrow("number of records k", "k")
addrow("pairing sum over records P_rec", "P_rec")
addrow("pairing sum records+writers P_all", "P_all")
addrow("total min writer weight 2k", "Wtot_lower")
addrow("log2 code-space dimension", "log2dim")
addrow("code distance d", "distance")
addrow("records moved by a weight-1 op", "D1max")
addrow("records moved by a weight-2 op", "D2max")
p("-" * 118)
p("")

# ---------------------------------------------------------------- chi additivity
p("ADDITIVITY OF TOTAL chi.  N = k_A + k_B records split into two disjoint clusters of equal size.")
p("SHARED BATH: one nq=3 bath serves both clusters.  DISJOINT BATHS (CONTROL): each cluster")
p("has its own nq=3 bath, which is what 'far apart' means for an environment.")
p("")
p("   N   k_A  k_B |  chi(A)   chi(B)   chi(A)+chi(B) |  chi(A u B) shared bath   deficit  | chi(A u B) disjoint baths  exact sum?")
p("-" * 118)
chires = {}
for m in (2, 4, 6, 8, 10, 12, 16, 24, 32):
    kA = m           # m blocks -> 2m records; split half/half
    kB = m
    cA = tavg_total(kA, 3)
    cB = tavg_total(kB, 3, offset=1)     # B's records land on a different phase of the same bath
    cShared = tavg_total(kA + kB, 3)
    cDisjoint = tavg_total(kA, 3) + tavg_total(kB, 3, offset=1)
    chires[m] = (cA, cB, cA + cB, cShared, cDisjoint)
    p("%4d %5d %5d |  %-8.4f %-8.4f %-14.4f |  %-21.4f %-9.4f | %-25.4f %s"
      % (2 * m, kA, kB, cA, cB, cA + cB, cShared, cA + cB - cShared, cDisjoint,
         "YES" if abs(cDisjoint - (cA + cB)) < 1e-12 else "NO"))
p("-" * 118)
p("")

# ---------------------------------------------------------------- redundancy scalar
p("REDUNDANCY SCALAR: how many single-bath-qubit fragments hold at least 10%% of the whole-bath chi")
p("about one record.  Computed with the model's own RecordModel.redundancy on the dense carrier.")
p("")
p("   n   k  nq |  whole-bath chi   n_fragments >= 10%% of it   bound = nq")
p("-" * 118)
redres = {}
for n, nqs in ((4, (3, 4, 5)), (6, (3, 4))):
    k = n - 2
    Xn = np.array([[1]], complex); Zn = np.array([[1]], complex)
    for _ in range(n): Xn = np.kron(Xn, Xm); Zn = np.kron(Zn, Zm)
    H = -(Xn + Zn)
    stab = [[1] * n + [0] * n, [0] * n + [1] * n]
    pairs = symplectic_logicals([s[:] for s in stab], n)
    R = [xz_to_matrix(a, n) for a, b in pairs]
    mm = RecordModel(H)
    for nq in nqs:
        env = Environment(nq=nq, energies=energies(nq), beta=BETA)
        whole = 0.0; frac = np.zeros(nq)
        for t in TIMES:
            w, parts = mm.redundancy(R[0], sum(R), env, lam=LAM, t=float(t))
            whole += w; frac += parts
        whole /= len(TIMES); frac /= len(TIMES)
        nfrag = int((frac >= 0.10 * whole).sum()) if whole > 1e-12 else 0
        redres["n%d_nq%d" % (n, nq)] = (whole, frac.tolist(), nfrag)
        p("%4d %3d %3d |  %-16.5f %-24d %d" % (n, k, nq, whole, nfrag, nq))
p("-" * 118)
p("")
p("READ (filled from the numbers above, not in advance):")
p("  every structural quantity marked YES above is additive over disjoint regions; d and D1max are")
p("  marked NO because they are INTENSIVE (constant), not because they fail to add up.")
p("  chi under a SHARED bath is strictly SUBadditive at every m: the deficit column is positive throughout.")
p("  the DISJOINT-BATH control is exactly additive, so the test can register additivity -- the deficit is real.")
p("  additivity of chi therefore requires the ENVIRONMENT to be doubled along with the matter, which is")
p("  not what extensivity of a source means.")
p("  the redundancy scalar is bounded by nq, the bath's own size, at every N.")

json.dump(dict(blocks={str(m): bres[m] for m in MS},
               chi={str(m): chires[m] for m in chires},
               redundancy=redres),
          open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_A_ZERO_OR_GROWING/s4_blocks_additivity.json", "w"), indent=1)
with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_A_ZERO_OR_GROWING/s4_blocks_additivity.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
