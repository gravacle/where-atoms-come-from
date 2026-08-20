"""O-50-A  step 3.  DENSE verification at L=2 (dim 256 -- the largest torus a dense matrix can
   reach; L=3 would be 262144 and is refused by the brief).  Everything here goes through the
   PROGRAM'S OWN commuting_family / independently_writable / joint_basis -- nothing asserted."""
import sys, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_A_ACTION")
from record_model import (RecordModel, symplectic_logicals, xz_to_matrix,
                          eigenspaces, clause_iii, clause_iv)
from f2lib import Toric, sp, rank, in_span, nullspace, span

np.set_printoptions(linewidth=200)
L = 2; T = Toric(L); n = T.n
S = [xz_to_matrix(s, n) for s in T.stab]
H = -sum(S)
print(f"L={L}  n={n}  dim={H.shape[0]}")
print("  H Hermitian                       ", np.linalg.norm(H - H.conj().T) < 1e-9)
print("  every stabiliser commutes with H  ", all(np.linalg.norm(g @ H - H @ g) < 1e-9 for g in S))
es = eigenspaces(H)
print("  eigenvalues / multiplicities      ", [(round(float(w), 6), m) for w, _, m in es])

pairs = symplectic_logicals(T.stab, n)
flat = [v for pr in pairs for v in pr]
print(f"  symplectic_logicals returned {len(pairs)} CONJUGATE PAIRS (D-18)")

# all 15 non-identity logical classes, as dense Hermitian involutions -- CANDIDATE records
cands, cand_bits = [], []
for bits in itertools.product((0, 1), repeat=len(flat)):
    if not any(bits): continue
    v = [0] * (2 * n)
    for b, f in zip(bits, flat):
        if b: v = [(x + y) % 2 for x, y in zip(v, f)]
    cands.append((bits, v, xz_to_matrix(v, n)))
print(f"  candidate record operators (non-identity logical classes): {len(cands)}")

I = np.eye(2 ** n)
rows = []
for bits, v, M in cands:
    c1 = np.linalg.norm(M - M.conj().T) < 1e-9 and np.linalg.norm(M @ M - I) < 1e-9
    c2 = np.linalg.norm(M @ H - H @ M) < 1e-9 and all(np.linalg.norm(M @ g - g @ M) < 1e-9 for g in S)
    c3 = clause_iii(M, es); c4 = clause_iv(M, es)
    rows.append((bits, c1, c2, c3, c4))
print("\n  CLAUSE TABLE over all 15 logical classes  (i) bit  (ii) durable  (iii) non-trivial  (iv) writable")
for bits, c1, c2, c3, c4 in rows:
    print(f"    class {bits}   (i)={c1}  (ii)={c2}  (iii)={c3}  (iv)={c4}")
print("  ALL 15 satisfy (i)-(iv):", all(all(r[1:]) for r in rows))

# CONTROL (D-15): operators that should FAIL a clause
ctrl = []
X0 = xz_to_matrix([1 if i == 0 else 0 for i in range(2 * n)], n)         # single-qubit X, not in N(S)
ctrl.append(("single X on edge 0 (not in N(S))", X0))
ctrl.append(("a stabiliser A_v (in S: trivial on code space)", xz_to_matrix(T.stab[0], n)))
ctrl.append(("identity", I.astype(complex)))
print("\n  CONTROLS")
for name, M in ctrl:
    print(f"    {name:46s} (i)={np.linalg.norm(M@M-I)<1e-9}  "
          f"(ii)={np.linalg.norm(M@H-H@M)<1e-9}  (iii)={clause_iii(M,es)}  (iv)={clause_iv(M,es)}")

# ---- the program's own multi-record machinery
rm = RecordModel(H, S)
print(f"\n  RecordModel: dim={rm.n}  minimal projections in A' = {len(rm.projs)}")
try:
    rm.records(); print("  .records() returned")
except RuntimeError as e:
    print("  .records() RAISED as expected (O-28):", e)

fam = rm.commuting_family([M for _, _, M in cands])
print(f"\n  commuting_family -> MAXIMAL COMMUTING FAMILY SIZE = {len(fam)}")
fam_bits = []
for R in fam:
    for bits, v, M in cands:
        if np.linalg.norm(M - R) < 1e-9: fam_bits.append(bits)
print("  family classes:", fam_bits)
print("  family members mutually commute:",
      all(np.linalg.norm(a @ b - b @ a) < 1e-9 for a, b in itertools.combinations(fam, 2)))

iw = rm.independently_writable(fam)
print(f"  independently_writable -> {iw}   (all of range({len(fam)})? {sorted(iw)==list(range(len(fam)))})")

jb = rm.joint_basis(fam)
from collections import defaultdict
per_E = defaultdict(dict)
for (ei, lab), C in jb.items(): per_E[ei][lab] = C.shape[1]
print("\n  JOINT BLOCK DIMENSIONS per eigenspace (theorem part (a): all equal within an eigenspace)")
allequal = True
for ei in sorted(per_E):
    d = per_E[ei]
    eq = (len(set(d.values())) == 1 and len(d) == 2 ** len(fam))
    allequal &= eq
    print(f"    E_{ei} (mult {es[ei][2]:3d}): {dict(sorted(d.items()))}   all 2^k blocks present & equal = {eq}")
print("  every eigenspace has 2^k equal joint blocks:", allequal)

# ---- verify the F_2-found writers act as claimed, densely
print("\n  DENSE VERIFICATION OF THE SEARCHED WRITERS (from s2, recomputed here by search)")
rows_ns = [[sp([1 if q == j else 0 for q in range(2 * n)], s, n) for j in range(2 * n)] for s in T.stab]
NS = nullspace(rows_ns, 2 * n)
Rv = []
for R in fam:
    for bits, v, M in cands:
        if np.linalg.norm(M - R) < 1e-9: Rv.append(v)
for j in range(len(fam)):
    target = tuple(1 if i == j else 0 for i in range(len(fam)))
    best = min((v for v in span(NS, 2 * n)
                if tuple(sp(v, r, n) for r in Rv) == target),
               key=lambda v: sum(1 for i in range(n) if v[i] or v[n + i]))
    U = xz_to_matrix(best, n)
    ok_H = np.linalg.norm(U @ H - H @ U) < 1e-9
    ok_u = np.linalg.norm(U.conj().T @ U - I) < 1e-9
    flips = [i for i in range(len(fam)) if np.linalg.norm(U.conj().T @ fam[i] @ U + fam[i]) < 1e-9]
    fixes = [i for i in range(len(fam)) if np.linalg.norm(U.conj().T @ fam[i] @ U - fam[i]) < 1e-9]
    print(f"    U_{j+1}: support={[i for i in range(n) if best[i] or best[n+i]]}  unitary={ok_u}  "
          f"[U,H]=0 -> {ok_H}  flips={flips}  fixes={fixes}")

# ---- CONTROL: a NON-admissible operation that does move a record
Ubad = xz_to_matrix([1 if i == 0 else 0 for i in range(2 * n)], n)
print(f"    CONTROL single-X: [U,H]=0 -> {np.linalg.norm(Ubad@H-H@Ubad)<1e-9}  "
      f"||[U,H]|| = {np.linalg.norm(Ubad@H-H@Ubad):.6f}")

# ---- the action, densely: apply the writers to a joint eigenvector and read the signs
print("\n  ACTION ON RECORD CONFIGURATIONS, read off dense states")
labs = sorted({lab for (ei, lab) in jb if ei == 0})
print("    ground eigenspace joint labels:", labs)
gvecs = {lab: jb[(0, lab)] for lab in labs}
Us = []
for j in range(len(fam)):
    target = tuple(1 if i == j else 0 for i in range(len(fam)))
    best = min((v for v in span(NS, 2 * n) if tuple(sp(v, r, n) for r in Rv) == target),
               key=lambda v: sum(1 for i in range(n) if v[i] or v[n + i]))
    Us.append(xz_to_matrix(best, n))
for lab, C in gvecs.items():
    psi = C[:, 0]
    for j, U in enumerate(Us):
        phi = U @ psi
        sig = tuple(int(np.sign(np.real(np.vdot(phi, R @ phi)))) for R in fam)
        print(f"    config {lab}  --U_{j+1}-->  {sig}")
