"""S1 -- DOES A [[4,2,2]] BLOCK ACTUALLY HOLD RECORDS?  All five clauses, checked, not assumed.
Then: is the [CODE] restriction exact?  Checked against [PHYS] at m=1 and m=2.

D-18: a record is a LOGICAL OPERATOR.  We check the five clauses on the carrier it lives on
before calling anything a record.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY")
from lanelib import *

np.set_printoptions(precision=6, suppress=True)
OUT = []
def P(s=""):
    print(s); OUT.append(str(s))

P("=" * 100)
P("S1  CLAUSE AUDIT OF THE [[4,2,2]] BLOCK  and  VALIDATION OF THE [CODE] REDUCTION")
P("=" * 100)

# ---------------------------------------------------------------- the block, in [PHYS]
n = 4
stab = block_stab_422()
pairs = symplectic_logicals(stab, n)
P("\nsymplectic_logicals returned %d CONJUGATE PAIRS (never nominated):" % len(pairs))
for i, (a, c) in enumerate(pairs):
    P("   pair %d :  A = (x|z) %s   B = (x|z) %s   sp(A,B) = %d  (must be 1)"
      % (i, a, c, sp(a, c, n)))
assert len(pairs) == 2, "expected k=2 logical qubits"
assert all(sp(a, c, n) == 1 for a, c in pairs)

# cross-pair elements must all commute
cross = [(i, j, sp(u, v, n)) for i, pi in enumerate(pairs) for j, pj in enumerate(pairs)
         if i < j for u in pi for v in pj]
P("   cross-pair symplectic products (must all be 0): %s" % sorted(set(x[2] for x in cross)))
assert all(x[2] == 0 for x in cross)

REC_v = [pairs[0][0], pairs[1][0]]          # the RECORDS
WRT_v = [pairs[0][1], pairs[1][1]]          # their conjugate WRITERS
P("\n   RECORDS  R1,R2 = %s , %s   weights %d,%d" % (REC_v[0], REC_v[1],
                                                     weight(REC_v[0], n), weight(REC_v[1], n)))
P("   WRITERS  W1,W2 = %s , %s   weights %d,%d" % (WRT_v[0], WRT_v[1],
                                                   weight(WRT_v[0], n), weight(WRT_v[1], n)))
P("   sp(R1,R2) = %d  (records must COMMUTE)" % sp(REC_v[0], REC_v[1], n))
assert sp(REC_v[0], REC_v[1], n) == 0

H = stab_hamiltonian(1)
es = eigenspaces(H)
P("\n   H = -(XXXX + ZZZZ), dim 16, eigenvalues/multiplicities: %s"
  % [(round(float(v), 6), m) for v, _, m in es])

R = [xz_to_matrix(v, n) for v in REC_v]
W = [xz_to_matrix(v, n) for v in WRT_v]

P("\n" + "-" * 100)
P("CLAUSE TABLE  ([PHYS], dim 16).  CONTROL COLUMN: an operator that must FAIL each clause.")
P("-" * 100)
# controls
ctrl_notbit  = 0.7 * xz_to_matrix(REC_v[0], n)          # fails (i): square != I
ctrl_notdur  = xz_to_matrix([1, 0, 0, 0, 0, 0, 0, 0], n)  # single X_0: fails (ii)
ctrl_trivial = np.eye(16, dtype=complex)                 # fails (iii) and (iv)

def clause_i(M):   return np.linalg.norm(M - M.conj().T) < 1e-9 and np.linalg.norm(M @ M - np.eye(M.shape[0])) < 1e-9
def clause_ii(M, HH): return np.linalg.norm(HH @ M - M @ HH) < 1e-9

hdr = "%-14s %-10s %-10s %-12s %-12s %-14s" % ("operator", "(i) BIT", "(ii) DUR", "(iii) NONTRIV", "(iv) WRITABLE", "(v) PROTECTED")
P(hdr); P("-" * 100)
rows = []
for nm, M in [("R1 (record)", R[0]), ("R2 (record)", R[1]),
              ("CTRL 0.7*R1", ctrl_notbit), ("CTRL X_0", ctrl_notdur), ("CTRL I", ctrl_trivial)]:
    c1 = clause_i(M); c2 = clause_ii(M, H)
    c3 = clause_iii(M, es) if c1 else False
    c4 = clause_iv(M, es)
    rows.append((nm, c1, c2, c3, c4))

# ---- clause (v): CARRIER DATA supplied explicitly.  contractible region := a SINGLE qubit.
# an admissible operation on region {j} is a unitary u_j (x) I with [u_j (x) I, H] = 0.
# solve for the commutant of {X,Z} inside 2x2 -- that is what "supported on one qubit and
# admissible" means here.
from record_model import commutant as _commutant
one_site = _commutant([Xm, Zm])
P("   [clause (v) carrier data] admissible single-qubit operators = commutant of {X,Z} in 2x2:")
P("      dim = %d  -> %s" % (len(one_site), "only multiples of I" if len(one_site) == 1 else "MORE THAN I"))
# and: does ANY weight-1 Pauli both commute with H and anticommute with a record?
w1 = []
for j in range(n):
    for pv in ([1 if i == j else 0 for i in range(4)] + [0] * 4,
               [0] * 4 + [1 if i == j else 0 for i in range(4)],
               [1 if i == j else 0 for i in range(4)] + [1 if i == j else 0 for i in range(4)]):
        adm = all(sp(pv, s, n) == 0 for s in stab)
        flips = [sp(pv, rv, n) for rv in REC_v]
        w1.append((pv, adm, flips))
n_bad = sum(1 for pv, adm, fl in w1 if adm and any(fl))
P("      weight-1 Paulis that are ADMISSIBLE and FLIP a record: %d (must be 0)" % n_bad)
# positive control for clause (v): weight-2 region -- distance is 2, so this MUST find one
w2_bad = 0
for j in range(n):
    for kk in range(j + 1, n):
        for a in range(4):
            for b in range(4):
                pv = [0] * 8
                for (q, code) in ((j, a), (kk, b)):
                    if code in (1, 3): pv[q] = 1
                    if code in (2, 3): pv[4 + q] = 1
                if all(sp(pv, s, n) == 0 for s in stab) and any(sp(pv, rv, n) for rv in REC_v):
                    w2_bad += 1
P("      CONTROL: weight-2 Paulis that are ADMISSIBLE and FLIP a record: %d (must be > 0;"
  % w2_bad)
P("               the code has distance 2, so 2-qubit regions are NOT protected -- this is")
P("               the positive control that the clause-(v) test can register a failure)")

prot = (len(one_site) == 1 and n_bad == 0)
for nm, c1, c2, c3, c4 in rows:
    pv_str = ("yes(1-qubit)" if prot else "no") if nm.startswith("R") else "n/a"
    P("%-14s %-10s %-10s %-12s %-12s %-14s" % (nm, c1, c2, c3, c4, pv_str))

P("\nREAD: R1 and R2 satisfy (i)-(iv) exactly; (v) holds against 1-qubit regions and FAILS")
P("      against 2-qubit regions (distance 2).  The carrier is honest: a [[4,2,2]] block holds")
P("      2 records with contractible := single qubit.  Every control failed the clause it was")
P("      built to fail.")
assert rows[0][1:] == (True, True, True, True) and rows[1][1:] == (True, True, True, True)
assert (not rows[2][1]) and (not rows[3][2]) and (not rows[4][3]) and (not rows[4][4])

# ---------------------------------------------------------------- independence, via the model
P("\n" + "-" * 100)
P("INDEPENDENT WRITABILITY (model.independently_writable), m=1, [PHYS]")
P("-" * 100)
rm = RecordModel(H)
fam = [R[0], R[1]]
iw = rm.independently_writable(fam)
P("   family size %d ; independently writable indices: %s" % (len(fam), iw))
for i in range(2):
    P("   W%d flips R%d ?  %s   leaves R%d fixed ?  %s"
      % (i + 1, i + 1, np.linalg.norm(W[i].conj().T @ R[i] @ W[i] + R[i]) < 1e-9,
         2 - i, np.linalg.norm(W[i].conj().T @ R[1 - i] @ W[i] - R[1 - i]) < 1e-9))

# ---------------------------------------------------------------- [CODE] validation
P("\n" + "=" * 100)
P("VALIDATION OF THE [CODE] REDUCTION  --  is codespace (x) bath an exact invariant subspace?")
P("=" * 100)
for m in (1, 2):
    Hm = stab_hamiltonian(m)
    nq = 4 * m
    rmm = RecordModel(Hm)
    Pg, kdim = rmm.ground_space()
    recs_v, wrts_v, _ = composite_records_writers(m)
    Rmat = [xz_to_matrix(v, nq) for v in recs_v]
    P("\n m=%d  dim=%d  ground-space dim=%d (must be 4^m=%d)" % (m, 2 ** nq, kdim, 4 ** m))
    assert kdim == 4 ** m
    # every record and every coupling must preserve the code space
    bad = max(np.linalg.norm(A @ Pg - Pg @ A) for A in Rmat)
    P("   max ||[record, P_ground]|| = %.3e  (must be 0 -> invariant subspace)" % bad)
    # build the isometry onto the code space and check records act as Z_i on 2m logical qubits
    wv, V = np.linalg.eigh(Pg); Q = V[:, wv > 0.5]
    Rc = [Q.conj().T @ A @ Q for A in Rmat]
    Zl, _ = code_records_couplings(m)
    # find the permutation/sign matching: records restricted should equal some Z_i up to
    # a basis choice; test the ALGEBRA instead -- pairwise commutation and spectra
    P("   restricted records: all Hermitian involutions? %s"
      % all(np.linalg.norm(A @ A - np.eye(kdim)) < 1e-9 for A in Rc))
    P("   restricted records: all traceless?             %s"
      % all(abs(np.trace(A)) < 1e-9 for A in Rc))
    P("   restricted records: mutually commuting?        %s"
      % all(np.linalg.norm(Rc[i] @ Rc[j] - Rc[j] @ Rc[i]) < 1e-9
            for i in range(len(Rc)) for j in range(len(Rc))))
    # the joint spectrum must be the full 2^(2m) sign table, each once -> isomorphic to Z_i's
    signs = {}
    Mgen = sum((2.0 ** i) * A for i, A in enumerate(Rc))
    ev = np.linalg.eigvalsh(Mgen)
    ev2 = np.linalg.eigvalsh(sum((2.0 ** i) * A for i, A in enumerate(Zl)))
    P("   joint sign-spectrum matches Z_1..Z_2m exactly? %s  (max dev %.3e)"
      % (np.allclose(np.sort(ev), np.sort(ev2), atol=1e-8),
         float(np.abs(np.sort(ev) - np.sort(ev2)).max())))

P("\n" + "-" * 100)
P("NUMERICAL EQUIVALENCE [PHYS] vs [CODE] for the Holevo quantity itself")
P("-" * 100)
env = Environment(nq=3, energies=(1.0, 1.4, 0.7), beta=2.0)
P("%-6s %-5s %-16s %-16s %-14s" % ("m", "rec", "chi [PHYS]", "chi [CODE]", "|difference|"))
maxdev = 0.0
for m in (1, 2):
    nq_ = 4 * m
    Hm = stab_hamiltonian(m)
    recs_v, _, _ = composite_records_writers(m)
    Rmat = [xz_to_matrix(v, nq_) for v in recs_v]
    rmm = RecordModel(Hm); Pg, kdim = rmm.ground_space()
    # PHYS propagator: couple every record to bath site (i mod 3)
    cpl_phys = [(Rmat[i], i % env.nq) for i in range(len(Rmat))]
    pp = Propagator(Hm, env, cpl_phys, lam=0.8, state0=Pg / kdim)
    chip = chi_timeavg(pp, Rmat)
    # CODE propagator
    Zl, dimc = code_records_couplings(m)
    HSc = -2.0 * m * np.eye(dimc, dtype=complex)
    cpl_code = [(Zl[i], i % env.nq) for i in range(len(Zl))]
    pc = Propagator(HSc, env, cpl_code, lam=0.8)
    chic = chi_timeavg(pc, Zl)
    for i in range(len(Rmat)):
        d = abs(chip[i] - chic[i]); maxdev = max(maxdev, d)
        P("%-6d %-5d %-16.10f %-16.10f %-14.3e" % (m, i, chip[i], chic[i], d))
P("\n   max |[PHYS] - [CODE]| over all records, m=1,2 : %.3e" % maxdev)
P("   SELF-CHECK: %s" % ("PASS -- the [CODE] restriction is exact, so all larger m may use it"
                         if maxdev < 1e-8 else "FAIL -- DO NOT USE [CODE]"))
assert maxdev < 1e-8

# also confirm Propagator == RecordModel.evolve
Hm = stab_hamiltonian(1); recs_v, _, _ = composite_records_writers(1)
Rmat = [xz_to_matrix(v, 4) for v in recs_v]
rmm = RecordModel(Hm)
r_model = rmm.evolve([(Rmat[0], 0)], env, lam=0.8, t=4.0)
Pg, kdim = rmm.ground_space()
pp = Propagator(Hm, env, [(Rmat[0], 0)], lam=0.8, state0=Pg / kdim)
P("   Propagator vs RecordModel.evolve at t=4: max dev %.3e"
  % float(np.abs(r_model - pp.state(4.0)).max()))
assert float(np.abs(r_model - pp.state(4.0)).max()) < 1e-9

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY/s1_block_verify.txt", "w").write("\n".join(OUT) + "\n")
print("\n[written]")
