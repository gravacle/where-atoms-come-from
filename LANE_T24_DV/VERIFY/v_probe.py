"""ADVERSARIAL VERIFICATION, part 2.

(A) EXACT: the full diagonal sector of the admissible algebra.  A diagonal D commutes
    with H iff D is constant on gauge orbits (all off-diagonal H entries are gauge-perm
    entries of one sign -- no cancellation).  So the COMPLETE set of diagonal record
    candidates is the +-1 orbit-constant functions; a diagonal record exists iff some
    +-1 orbit assignment has Tr(P_E D) = 0 on every eigenspace.  This SUPERSEDES the
    lane's 16383-member family search on the diagonal side: it covers ALL diagonals.
    Search is exact (meet-in-the-middle over per-orbit eigen-trace vectors, x64 integer
    numerators).

(B) NUMERIC (labeled; demonstration only): dense eigh of H, eigenvalue multiplicities;
    an R_moved witness (random balanced eigenspace split): clause residuals, transport
    movement ||[A_v(k), R]||, and ||{T, R}|| != 0 as the exact obstruction requires.
"""
import itertools, json, time
from fractions import Fraction
import numpy as np
from v_core_defs import group_D4, Car   # re-use my own fresh definitions

t0 = time.time()
say = lambda s="": print(s, flush=True)
REPORT = {}

G4 = group_D4()
n = G4["n"]; MUL, INV = G4["MUL"], G4["INV"]
car = Car(G4)
N = car.N
idx = np.arange(N)

# ---------------------------------------------------------------- (A) gauge orbits
gens = [car.A0(k) for k in range(n)] + [car.A1(k) for k in range(n)]
# no-cancellation check for the orbit argument: all off-diagonal H entries share a sign
# (they are -(1/8) times counts of gauge perms; B terms are diagonal) -- structural, noted.
parent = np.arange(N)
def find(x):
    root = x
    while parent[root] != root: root = parent[root]
    while parent[x] != root:
        parent[x], x = root, parent[x]
    return root
for g in gens:
    for c in range(N):
        a, b = find(c), find(int(g[c]))
        if a != b: parent[a] = b
roots = np.array([find(c) for c in range(N)])
orbit_ids = {r: i for i, r in enumerate(sorted(set(int(r) for r in roots)))}
orb_of = np.array([orbit_ids[int(r)] for r in roots])
n_orb = len(orbit_ids)
say("(A) gauge orbits on the 4096 configs: %d" % n_orb)
REPORT["n_gauge_orbits"] = int(n_orb)

# per-config eigenspace weights: 64 * diag(P_E) as exact integers
# diag of sector (a0,a1,b0,b1) at c = sum over (x0,x1) sign (1/8)^{x0+x1} sum_{k,k'} [A0A1(c)=c] d(c)
db0, db1 = car.dB0(), car.dB1()
diag64 = {}   # eig -> np.int64 vector = 64 * diag(P_E)
for s in itertools.product((0, 1), repeat=4):
    a0, a1, b0, b1 = s
    d = (db0 if b0 else 1 - db0) * (db1 if b1 else 1 - db1)
    acc = np.zeros(N, dtype=np.int64)
    terms0 = [(1, True)] if a0 else [(1, False), (-1, True)]
    terms1 = [(1, True)] if a1 else [(1, False), (-1, True)]
    for s0, useA0 in terms0:
        for s1, useA1 in terms1:
            scale = (1 if useA0 else 8) * (1 if useA1 else 8)   # (1/8)^{x} -> x64 numerator
            sub = np.zeros(N, dtype=np.int64)
            for k0 in (range(n) if useA0 else [None]):
                p0 = car.A0(k0) if k0 is not None else idx
                for k1 in (range(n) if useA1 else [None]):
                    p01 = p0[car.A1(k1)] if k1 is not None else p0
                    sub += (p01 == idx).astype(np.int64)
            acc += s0 * s1 * scale * sub * d
    k = -(a0 + a1 + b0 + b1)
    diag64[k] = diag64.get(k, np.zeros(N, dtype=np.int64)) + acc
# sanity: sums = 64 * dims
dims_chk = {k: int(v.sum()) for k, v in diag64.items()}
assert dims_chk == {0: 64*2686, -1: 64*864, -2: 64*476, -3: 64*48, -4: 64*22}, dims_chk
# per-orbit eigen-trace numerators (x64)
Es = [0, -1, -2, -3, -4]
orb_vec = np.zeros((n_orb, 5), dtype=np.int64)
for ei, Ek in enumerate(Es):
    np.add.at(orb_vec[:, ei], orb_of, diag64[Ek])
tot = orb_vec.sum(axis=0)
say("    per-orbit eigen-trace numerators computed; totals/64 = %s" % (tot // 64).tolist())
# balanced +-1 assignment: orbits with EQUAL eigen-trace vectors are interchangeable, so
# group by distinct vector u (multiplicity m_u) and enumerate a_u = (#plus - #minus) in
# {-m_u, -m_u+2, ..., m_u}; a diagonal record exists iff sum_u a_u * u = 0 has a solution.
# (all +1 / all -1 give sums = +-64*dims != 0, so triviality is excluded automatically)
from collections import Counter
cnt = Counter(tuple(int(x) for x in v) for v in orb_vec)
types = sorted(cnt.items())
say("    distinct orbit vectors: %d : %s" % (len(types), types))
assert len(types) <= 7, "unexpectedly many types"
U = np.array([t[0] for t in types], dtype=np.int64)          # (T, 5)
Ms = [t[1] for t in types]
ranges = [np.arange(-m, m + 1, 2, dtype=np.int64) for m in Ms]
found = 0
n_combos = int(np.prod([len(r) for r in ranges]))
say("    exact enumeration over %s combinations" % n_combos)
# vectorize over the last two types
last2 = np.stack(np.meshgrid(ranges[-2], ranges[-1], indexing="ij"), axis=-1).reshape(-1, 2)
tail = last2[:, :1] * U[-2] + last2[:, 1:] * U[-1]           # (K, 5)
import itertools as _it
for head in _it.product(*ranges[:-2]):
    base = sum(a * u for a, u in zip(head, U[:-2]))
    hit = np.all(tail == -base, axis=1)
    found += int(hit.sum())
say("    balanced +-1 orbit assignments (diagonal records) found: %d" % found)
say("    => %s" % ("NO diagonal record exists on this carrier -- the lane's 'necessarily"
                   " non-diagonal' claim is TRUE (now proved over ALL diagonals, not just its family)"
                   if found == 0 else
                   "DIAGONAL RECORD EXISTS -- lane's 'necessarily non-diagonal' claim is FALSE"))
REPORT["diagonal_records_found"] = int(found)

# ---------------------------------------------------------------- (B) numeric witnesses
say("")
say("(B) NUMERIC (floats; demonstration only)")
H = np.zeros((N, N))
for k in range(n):
    for p in (car.A0(k), car.A1(k)):
        H[p, idx] += -1.0 / n
H[idx, idx] += -(db0 + db1)
w, V = np.linalg.eigh(H)
wr = np.round(w).astype(int)
assert np.max(np.abs(w - wr)) < 1e-9
mult = {int(k): int((wr == k).sum()) for k in sorted(set(wr.tolist()))}
say("    eigh multiplicities: %s (exact: {-4:22,-3:48,-2:476,-1:864,0:2686})" % mult)
REPORT["eigh_mult_match"] = (mult == {-4: 22, -3: 48, -2: 476, -1: 864, 0: 2686})

rng = np.random.default_rng(20260820)
R = np.zeros((N, N))
for k in sorted(mult):
    cols = V[:, wr == k]
    m = cols.shape[1]
    signs = np.array([1.0] * (m // 2) + [-1.0] * (m // 2))
    rng.shuffle(signs)
    Q = cols @ np.diag(signs) @ cols.T
    R += Q
resid_inv = np.linalg.norm(R @ R - np.eye(N))
resid_sym = np.linalg.norm(R - R.T)
resid_H = np.linalg.norm(R @ H - H @ R)
bal = max(abs(np.trace(V[:, wr == k].T @ R @ V[:, wr == k])) for k in sorted(mult))
# transport movement
move = 0.0
for k in range(n):
    P = np.zeros((N, N)); P[car.A0(k), idx] = 1.0
    move = max(move, np.linalg.norm(P @ R - R @ P))
# {T, R} on h0
Tp = car.with_comp("h0", MUL[car.h0, G4["r2"]])
TP = np.zeros((N, N)); TP[Tp, idx] = 1.0
anti = np.linalg.norm(TP @ R + R @ TP)
say("    R_moved witness: ||R^2-I||=%.2e ||R-R^T||=%.2e ||[R,H]||=%.2e max|Tr(P_E R)|=%.2e" %
    (resid_inv, resid_sym, resid_H, bal))
say("    max_k ||[A_v0(k),R]|| = %.3e (MOVED) ; ||{T,R}||_F = %.3e (NONZERO as the exact obstruction demands)" %
    (move, anti))
REPORT["R_moved_witness"] = dict(inv=float(resid_inv), sym=float(resid_sym), commH=float(resid_H),
                                 bal=float(bal), moved=float(move), antiT=float(anti))
ok_witness = resid_inv < 1e-8 and resid_H < 1e-8 and bal < 1e-8 and move > 1.0 and anti > 1.0
say("    witness consistent with lane's section 9: %s" % ok_witness)
REPORT["witness_ok"] = bool(ok_witness)

REPORT["elapsed_s"] = round(time.time() - t0, 1)
with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_T24_DV/VERIFY/v_probe.json", "w") as f:
    json.dump(REPORT, f, indent=1)
say("")
say("REPORT: %s" % json.dumps(REPORT))
