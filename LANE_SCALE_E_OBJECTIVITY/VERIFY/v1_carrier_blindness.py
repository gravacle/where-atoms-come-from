"""V1 -- IS ANY REPORTED NUMBER A FUNCTION OF RECORD-HOOD AT ALL?  (D-15 / D-18)

THE ATTACK.  Every quantity the lane reports (chi_site, chi_whole, E, R_delta, Phi, EXCESS,
depth*, N_better, SYN, C3) is computed by common.Broadcast, whose constructor signature is
    Broadcast(k, nq, W, lam, times, energies, beta)
-- the CARRIER NEVER ENTERS.  s1 verifies the [[n,n-2,2]] clauses and s2 verifies the fast
path, but nothing downstream of s2 consumes either.  So the numbers cannot distinguish
"k protected records" from "k arbitrary commuting +-1 bits".

THIS SCRIPT RUNS THE MISSING CONTROL through record_model itself:

  CARRIER A  [[4,2,2]] -- H = -(XXXX + ZZZZ), records = the Z_i of symplectic_logicals.
             Satisfies (i)-(v).  PROTECTED against every single-qubit region (s1's claim).
  CARRIER B  2 FREE QUBITS -- H = 0, "records" = Z_0, Z_1.
             Satisfies (i)-(iv) but clause (v) FAILS: X_0 is a single-qubit (contractible)
             admissible operation that flips Z_0.  THESE ARE NOT RECORDS.

Same bath, same W, same lam, same times, same coupling operator sum_ij W[i,j] R_i (x) X_j.
If every chi agrees to machine precision, then every number in the lane is blind to clause (v)
-- i.e. blind to the only clause that makes a bit a RECORD -- and the finding is a statement
about k classical coins in a shared bath, not about records.

POSITIVE CONTROL (D-15): the SAME comparison run with a coupling that does NOT commute with
the observables must give DIFFERENT numbers on A and B, otherwise this test could not tell any
two carriers apart.
"""
import sys, itertools
import numpy as np

sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY")
from record_model import (RecordModel, Environment, symplectic_logicals, xz_to_matrix,
                          eigenspaces, clause_iii, clause_iv)
from common import weights, Broadcast, TIMES, ENERGIES, BETA, code_stabilisers, carrier

I2 = np.eye(2, dtype=complex)
XM = np.array([[0, 1], [1, 0]], dtype=complex)
ZM = np.array([[1, 0], [0, -1]], dtype=complex)
YM = 1j * XM @ ZM


def op(n, j, P):
    M = np.array([[1]], dtype=complex)
    for q in range(n):
        M = np.kron(M, P if q == j else I2)
    return M


def clause_v_single_qubit(n, H, recs):
    """Enumerate every single-qubit Pauli; admissible = commutes with H; flips = anticommutes
       with some record.  Returns (#paulis, #admissible, #flips)."""
    tot = adm = flip = 0
    for j in range(n):
        for P in (XM, YM, ZM):
            A = op(n, j, P)
            tot += 1
            if np.linalg.norm(A @ H - H @ A) < 1e-9:
                adm += 1
                if any(np.linalg.norm(A @ R + R @ A) < 1e-9 for R in recs):
                    flip += 1
    return tot, adm, flip


def build_A():
    n = 4
    car = carrier(n)
    Sx = xz_to_matrix(code_stabilisers(n)[0], n)
    Sz = xz_to_matrix(code_stabilisers(n)[1], n)
    H = -(Sx + Sz)
    recs = [xz_to_matrix(v, n) for v in car['recs_xz']]
    return n, H, recs


def build_B():
    n = 2
    H = np.zeros((4, 4), dtype=complex)
    recs = [op(2, 0, ZM), op(2, 1, ZM)]
    return n, H, recs


def clause_report(name, n, H, recs):
    es = eigenspaces(H)
    ci = all(np.linalg.norm(R @ R - np.eye(R.shape[0])) < 1e-9 and
             np.linalg.norm(R - R.conj().T) < 1e-9 for R in recs)
    cii = all(np.linalg.norm(H @ R - R @ H) < 1e-9 for R in recs)
    ciii = all(clause_iii(R, es) for R in recs)
    civ = all(clause_iv(R, es) for R in recs)
    tot, adm, flip = clause_v_single_qubit(n, H, recs)
    cv = (flip == 0)
    print(f"  {name:28s} (i)={ci} (ii)={cii} (iii)={ciii} (iv)={civ} | "
          f"single-qubit Paulis {tot:3d} admissible {adm:3d} flips {flip:3d} -> "
          f"(v)={'PROTECTED' if cv else 'FAILS -- NOT A RECORD'}")
    return dict(i=ci, ii=cii, iii=ciii, iv=civ, v=cv, flips=flip)


def chi_table(n, H, recs, W, lam, nq, times, commuting=True):
    """chi(R_i : fragment) time-averaged, plus the joint 4-outcome pair chi, through the FULL
       record_model (RecordModel.evolve + Environment.holevo).  If commuting=False the coupling
       uses X-type system operators that do NOT commute with the observables (positive control)."""
    env = Environment(nq=nq, energies=ENERGIES[:nq], beta=BETA)
    nS = H.shape[0]
    k = len(recs)
    if commuting:
        Cops = recs
    else:
        Cops = [op(n, i % n, XM) for i in range(k)]
    HINT = sum(W[i, j] * np.kron(Cops[i], env.site[j]) for i in range(k) for j in range(nq))
    M = RecordModel(H)
    frags = [[0], [1], list(range(nq))]
    acc = np.zeros((k, len(frags)))
    accpair = np.zeros(len(frags))
    for t in times:
        r = M.evolve(HINT, env, lam=lam, t=float(t))
        for fi, f in enumerate(frags):
            for i in range(k):
                acc[i, fi] += env.holevo(r, recs[i], nS, fragment=f)
            # 4-outcome joint Holevo on (R_0, R_1), assembled from the model's own partial traces
            outs = []
            for s0 in (+1, -1):
                for s1 in (+1, -1):
                    P = np.kron((np.eye(nS) + s0 * recs[0]) / 2 @ (np.eye(nS) + s1 * recs[1]) / 2,
                                np.eye(env.dim))
                    blk = P @ r @ P
                    p = float(np.real(np.trace(blk)))
                    if p < 1e-12:
                        continue
                    rB = env._trace_system(blk / p, nS)
                    outs.append((p, env._fragment(rB, f)))
            av = sum(p * rb for p, rb in outs)

            def vn(m):
                e = np.linalg.eigvalsh(m); e = e[e > 1e-13]
                return float(-(e * np.log2(e)).sum())
            accpair[fi] += max(vn(av) - sum(p * vn(rb) for p, rb in outs), 0.0)
    return acc / len(times), accpair / len(times)


print("=" * 110)
print("V1  CARRIER BLINDNESS -- do the lane's numbers know whether the bits are RECORDS?")
print("=" * 110)
print()
print("STEP 0  clause audit of the two carriers (D-18: check the clauses on the carrier used)")
nA, HA, recA = build_A()
nB, HB, recB = build_B()
cA = clause_report("A: [[4,2,2]] code", nA, HA, recA)
cB = clause_report("B: 2 free qubits (H=0)", nB, HB, recB)
print()
print("  SELF-CHECK: the two carriers must DIFFER on clause (v), or this control is dead.")
print(f"     A (v) = {cA['v']}   B (v) = {cB['v']}   differ = {cA['v'] != cB['v']}")
if cA['v'] == cB['v']:
    print("     SELF-CHECK FAILED -- draw no conclusion.")
    sys.exit(0)
print()

NQ = 4
LAM = 0.8
TS = np.linspace(1.0, 13.0, 9)          # 9 of the lane's 25 times; time-averaged either way
for kind in ('crowded', 'sym', 'separate'):
    W = weights(kind, 2, NQ, seed=7)
    a, ap = chi_table(nA, HA, recA, W, LAM, NQ, TS, commuting=True)
    b, bp = chi_table(nB, HB, recB, W, LAM, NQ, TS, commuting=True)
    d = float(np.abs(a - b).max())
    dp = float(np.abs(ap - bp).max())
    print(f"STEP 1  geometry {kind:9s}  k=2, nq={NQ}, lam={LAM}, 9 times in [1,13]")
    print(f"   fragment          |   site0     site1     whole  |  (RECORDS on [[4,2,2]])")
    for i in range(2):
        print(f"   chi(R_{i}) carrier A |  {a[i,0]:.6f}  {a[i,1]:.6f}  {a[i,2]:.6f}")
    print(f"   chi(R0,R1) carr. A |  {ap[0]:.6f}  {ap[1]:.6f}  {ap[2]:.6f}")
    print(f"   fragment          |   site0     site1     whole  |  (NOT RECORDS: clause (v) fails)")
    for i in range(2):
        print(f"   chi(R_{i}) carrier B |  {b[i,0]:.6f}  {b[i,1]:.6f}  {b[i,2]:.6f}")
    print(f"   chi(R0,R1) carr. B |  {bp[0]:.6f}  {bp[1]:.6f}  {bp[2]:.6f}")
    print(f"   MAX |A - B| single-record = {d:.3e}   pair = {dp:.3e}")
    print()

print("STEP 2  POSITIVE CONTROL (D-15) -- a THIRD carrier that is also a genuine record carrier")
print("        (all five clauses PASS) but whose ground space has a different SECTOR structure:")
print("        H = -Z_0 Z_1 on 2 qubits, ground space span{|00>,|11>} (dim 2, not 4), so the two")
print("        records are perfectly correlated in the initial state -- 2 joint sectors, not 4.")
print("        Same coupling, same bath.  If the instrument is live this MUST differ from A.")
nC, HC, recC = 2, -np.kron(ZM, ZM), [op(2, 0, ZM), op(2, 1, ZM)]
cC = clause_report("C: H = -Z0 Z1", nC, HC, recC)
W = weights('crowded', 2, NQ, seed=7)
a2, ap2 = chi_table(nA, HA, recA, W, LAM, NQ, TS, commuting=True)
c2, cp2 = chi_table(nC, HC, recC, W, LAM, NQ, TS, commuting=True)
print(f"   carrier A (4 sectors) chi(R_0) [site0 site1 whole] = {np.round(a2[0],6)}")
print(f"   carrier C (2 sectors) chi(R_0) [site0 site1 whole] = {np.round(c2[0],6)}")
print(f"   MAX |A - C| = {float(np.abs(a2-c2).max()):.6f}   (pair: {float(np.abs(ap2-cp2).max()):.6f})")
print()

print("STEP 3  and the lane's own fast path reproduces BOTH, since it only ever sees (k,nq,W,lam):")
W = weights('crowded', 2, NQ, seed=7)
B = Broadcast(2, NQ, W, LAM, times=TS)
fp = np.array([[np.mean([B.chi_single([j], i, ti) for ti in range(len(TS))])
                for j in (0, 1)] + [np.mean([B.chi_single(list(range(NQ)), i, ti)
                                             for ti in range(len(TS))])] for i in range(2)])
a, _ = chi_table(nA, HA, recA, W, LAM, NQ, TS, commuting=True)
b, _ = chi_table(nB, HB, recB, W, LAM, NQ, TS, commuting=True)
print(f"   fast path        chi(R_0) = {fp[0]}")
print(f"   |fast - A| = {float(np.abs(fp-a).max()):.3e}     |fast - B| = {float(np.abs(fp-b).max()):.3e}")
print()
print("READ (filled from the numbers above):")
same = float(np.abs(a - b).max()) < 1e-10
diff = float(np.abs(a2 - c2).max())
print(f"   commuting coupling: carriers agree to {float(np.abs(a-b).max()):.1e}  -> "
      f"{'IDENTICAL' if same else 'DIFFERENT'}")
print(f"   sector-structure positive control: carriers differ by {diff:.6f} -> "
      f"{'instrument is live' if diff > 1e-3 else 'INSTRUMENT DEAD, no conclusion'}")
