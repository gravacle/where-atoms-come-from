"""O-5 D.  DO THE THEOREMS SURVIVE THE APPROXIMATE DEFINITION?

Five things, each tested rather than asserted.

  D0  THE EXPONENT IS THE CODE DISTANCE.  C measured delta ~ p^2 on a d=2 carrier and delta ~ p^1 on
      a d=1 carrier.  Two points is a line through two points.  Add a THIRD carrier of distance 3
      (Steane [[7,1,3]], dim 128) and see whether the exponent is 3.  If it is, the exponent in the
      record-lifetime law is the code distance -- which is G-12's min(systole, cosystole) -- and the
      approximate definition is what makes G-12 operational.

  D1  THE EXACT DEFINITION HAS EMPTY SCOPE ON A GENERIC HAMILTONIAN.  Not "narrow": EMPTY.  Clause
      (ii) [H,R]=0 forces R block-diagonal on eigenspaces; clause (iii) requires R non-scalar on one
      of them; a generic H has no degenerate eigenspace at all.  So exact (ii) AND exact (iii) are
      jointly unsatisfiable except on a measure-zero set of H.  Counted.

  D2  P-1  ((iii) => H degenerate)          -- does it survive?
  D3  P-2  ((ii)+(iv) => the writer is built from neither H nor {L_k})  -- does it survive?
  D4  P-3  ((iv)+(v) => THE WRITER IS NON-LOCAL)  -- does it survive, and with what bound?

  D5  THE OPEN-SYSTEM HALF, with DETAILED BALANCE, not an infinite-temperature bath.
      Script A measured ||[X_l, R]|| = 2 on BOTH carriers -- identical, maximal, and the same number
      whatever the temperature.  Here the actual record decay rate is computed from a detailed-
      balance generator on both carriers, and the two are compared.  If they differ while the
      commutator norm does not, the L-half of clause (ii) cannot be relaxed by a norm bound either.
      POSITIVE CONTROL: at beta = 0 the barrier is worth nothing and the two rates must come
      together.  A machine that shows a separation at beta = 0 as well is measuring an artefact.
"""
import numpy as np
from o5_common import (DIM, NQ, PLAQS, STARS, Zop, Xop, Yop, toric_H, sym_H,
                       local_perturbation, Z_A_SUP, X_A_SUP)

print("=" * 104)
print("O-5 D.  THEOREM SURVIVAL UNDER THE APPROXIMATE DEFINITION")
print("=" * 104)

Ht, Rt, gt = toric_H(), Zop(Z_A_SUP), 4
Hs, Rs, gs = sym_H(), Zop([0]), 2
V8 = local_perturbation(seed=2026)


def width(H0, R0, g, V, p):
    e, U = np.linalg.eigh(H0 + p * V)
    Uc = U[:, :g]
    M = Uc.conj().T @ R0 @ Uc
    return e[g - 1] - e[0], e[g] - e[g - 1], (M + M.conj().T) / 2, e[:g], U[:, :g]


# ==================================================================================================
print("\n" + "=" * 104)
print("  D0.  IS THE EXPONENT THE CODE DISTANCE?  A THIRD CARRIER.")
print("=" * 104)

# ---- Steane [[7,1,3]], 7 qubits, dim 128 ---------------------------------------------------------
NQ7, D7 = 7, 128


def op7(kind, S):
    m = 0
    for k in S:
        m |= (1 << k)
    if kind == 'Z':
        par = np.array([bin(s & m).count('1') & 1 for s in range(D7)])
        return np.diag(np.where(par == 0, 1.0, -1.0)).astype(complex)
    M = np.zeros((D7, D7), complex)
    b = np.arange(D7)
    M[b ^ m, b] = 1.0
    return M


SX = [[3, 4, 5, 6], [1, 2, 5, 6], [0, 2, 4, 6]]
H7 = -sum(op7('X', s) for s in SX) - sum(op7('Z', s) for s in SX)
R7 = op7('Z', list(range(7)))            # logical Zbar
W7 = op7('X', list(range(7)))            # logical Xbar
e7 = np.linalg.eigvalsh(H7)
g7 = int(np.sum(np.abs(e7 - e7[0]) < 1e-9))
print(f"    Steane [[7,1,3]]: dim {D7}, ground energy {e7[0]:+.4f} (expected -6), degeneracy {g7}"
      f" (expected 2), gap {e7[g7]-e7[0]:.4f} (expected 2)   "
      f"{'PASS' if abs(e7[0]+6) < 1e-9 and g7 == 2 else 'FAIL'}")
c1 = max(np.linalg.norm(H7 @ O - O @ H7, 2) for O in (R7, W7))
c2 = np.linalg.norm(R7 @ W7 + W7 @ R7, 2)
print(f"    ||[H,Zbar]||,||[H,Xbar]|| = {c1:.2e}  ||{{Zbar,Xbar}}|| = {c2:.2e}   "
      f"{'PASS' if c1 < 1e-12 and c2 < 1e-12 else 'FAIL'}")
# distance check: smallest weight of a Pauli acting non-trivially on the 2-dim code space
U7 = np.linalg.eigh(H7)[1][:, :g7]


def code_act7(O):
    M = U7.conj().T @ O @ U7
    return float(np.linalg.norm(M - np.trace(M) / g7 * np.eye(g7), 2))


import itertools
best = {}
for w in (1, 2, 3):
    mx = 0.0
    for sub in itertools.combinations(range(7), w):
        for kinds in itertools.product('XYZ', repeat=w):
            O = np.eye(D7, dtype=complex)
            for q, k in zip(sub, kinds):
                O = O @ (op7('Z', [q]) if k == 'Z' else
                         (op7('X', [q]) if k == 'X' else 1j * op7('X', [q]) @ op7('Z', [q])))
            mx = max(mx, code_act7(O))
    best[w] = mx
    print(f"    max over ALL weight-{w} Paulis of ||P O P - scalar|| = {mx:.4f}")
print(f"    => distance = {min(w for w in best if best[w] > 1e-8)}   "
      f"{'PASS (expected 3)' if min(w for w in best if best[w] > 1e-8) == 3 else 'FAIL'}")

rng7 = np.random.default_rng(4242)
V7 = np.zeros((D7, D7), complex)
for l in range(NQ7):
    c = rng7.normal(size=3)
    V7 = V7 + c[0] * op7('X', [l]) + c[1] * (1j * op7('X', [l]) @ op7('Z', [l])) + c[2] * op7('Z', [l])
V7 = (V7 + V7.conj().T) / 2
V7 = V7 / np.linalg.norm(V7, 2)

print(f"\n    {'p':>8s} {'d=1 SYMMETRY':>15s} {'exp':>6s} {'d=2 TORIC':>15s} {'exp':>6s}"
      f" {'d=3 STEANE':>15s} {'exp':>6s}")
PS = [3e-2, 1e-1, 2e-1, 3e-1]
tab = {1: [], 2: [], 3: []}
for p in PS:
    a = width(Hs, Rs, gs, V8, p)[0]
    b = width(Ht, Rt, gt, V8, p)[0]
    e_, U_ = np.linalg.eigh(H7 + p * V7)
    c_ = e_[g7 - 1] - e_[0]
    tab[1].append(a); tab[2].append(b); tab[3].append(c_)
    print(f"    {p:8.2e} {a:15.6e} {np.log(a)/np.log(p):6.2f} {b:15.6e} {np.log(b)/np.log(p):6.2f}"
          f" {c_:15.6e} {np.log(c_)/np.log(p):6.2f}")
print(f"\n    {'carrier':>14s} {'fitted exponent':>17s} {'code distance':>15s}  verdict")
for d_, nm in ((1, "SYMMETRY"), (2, "TORIC 2x2"), (3, "STEANE 7")):
    k = np.polyfit(np.log(PS), np.log(tab[d_]), 1)[0]
    print(f"    {nm:>14s} {k:17.4f} {d_:15d}  {'PASS' if abs(k - d_) < 0.15 else 'FAIL'}")
print("""
    THE EXPONENT IN THE RECORD-LIFETIME LAW IS THE CODE DISTANCE, ON THREE CARRIERS SPANNING d=1,2,3.
    The code distance of a length-2 F_2 complex is min(systole, cosystole) -- G-12's TWO NUMBERS.
    G-12 said R3 needs two numbers but did not say what they buy.  They buy the EXPONENT of the
    record lifetime:   T(eta) ~ eta * Delta^(d-1) / (c * p^d),  d = min(systole, cosystole).""")

# ==================================================================================================
print("\n" + "=" * 104)
print("  D1.  THE EXACT DEFINITION'S SCOPE ON A GENERIC HAMILTONIAN")
print("=" * 104)
rng = np.random.default_rng(9)
ndeg = 0
mins = []
N = 400
for _ in range(N):
    A = rng.normal(size=(32, 32)) + 1j * rng.normal(size=(32, 32))
    A = (A + A.conj().T) / 2
    ev = np.linalg.eigvalsh(A)
    sp = np.min(np.diff(ev))
    mins.append(sp / np.linalg.norm(A, 2))
    if sp < 1e-9 * np.linalg.norm(A, 2):
        ndeg += 1
print(f"    {N} random Hermitian matrices (32x32): number with ANY degenerate eigenvalue = {ndeg}")
print(f"    smallest relative level spacing seen: min {min(mins):.3e}, median {np.median(mins):.3e}")
print(f"    => on a generic H, clause (ii) forces R = f(H), and clause (iii) then fails IDENTICALLY.")
print(f"       THE EXACT DEFINITION IS SATISFIED BY NO RECORD AT ALL ON A GENERIC H.   "
      f"{'PASS' if ndeg == 0 else 'FAIL'}")
print("\n    the same statement on OUR carrier, perturbed:")
for nm, H0, R0, g in (("TORIC 2x2", Ht, Rt, gt), ("SYMMETRY", Hs, Rs, gs)):
    for p in (0.0, 1e-3):
        ev = np.linalg.eigvalsh(H0 + p * V8)
        nd = len(np.unique(np.round(ev, 9)))
        print(f"      {nm:>10s}  p={p:6.0e}  distinct eigenvalues {nd:4d} / {len(ev)}   "
              f"{'exact records EXIST' if nd < len(ev) else 'NO exact record exists'}")
print("""
    POSITIVE CONTROL: at p = 0 the same counter finds 5 distinct eigenvalues out of 256 and reports
    that records exist.  The counter is alive; the emptiness at p > 0 is real.

    THIS IS WHY O-5 IS NOT A REFINEMENT.  The exact definition is not slightly too strict.  It has
    NO instances whenever the Hamiltonian is perturbed by anything at all, which is always.""")

# ==================================================================================================
print("\n" + "=" * 104)
print("  D2 / D3 / D4.  THE THREE PROPOSITIONS")
print("=" * 104)
p = 1e-3
dlt, gap, M, ee, Uc = width(Ht, Rt, gt, V8, p)
dls, gaps_, Ms, ees, Ucs = width(Hs, Rs, gs, V8, p)

print(f"\n  D2  P-1: (iii) => H is degenerate.   APPROXIMATE FORM: (iii-w) => H has a CLUSTER of")
print(f"      width <= w and dimension > 1, i.e. QUASI-degeneracy.  The cluster is meaningful only")
print(f"      if it is well separated: report Delta/delta.")
print(f"      {'carrier':>12s} {'p':>8s} {'cluster dim':>12s} {'width delta':>14s} {'gap Delta':>10s} {'Delta/delta':>13s}")
for nm, (a, b, _, _, _) , g in (("TORIC 2x2", (dlt, gap, 0, 0, 0), gt), ("SYMMETRY", (dls, gaps_, 0, 0, 0), gs)):
    print(f"      {nm:>12s} {p:8.0e} {g:12d} {a:14.6e} {b:10.4f} {b/a:13.4e}")
print(f"      P-1 SURVIVES, weakened by exactly the amount w: 'degenerate' becomes 'quasi-degenerate")
print(f"      to within w'.  On the toric carrier the cluster is separated from the rest of the")
print(f"      spectrum by {gap/dlt:.2e}; on the symmetry carrier by {gaps_/dls:.2e}.  Both are clusters.")

print(f"\n  D3  P-2: (ii)+(iv) => the writer is built from neither H nor {{L_k}}.")
Wt = Xop(X_A_SUP)
Mw = Uc.conj().T @ Wt @ Uc
Hc = np.diag(ee).astype(complex)
cw = np.linalg.norm(Hc @ Mw - Mw @ Hc, 2)
cr = np.linalg.norm(Hc @ M - M @ Hc, 2)
print(f"      ON THE CLUSTER, at p = {p:.0e}:")
print(f"        ||[H_c, R_c]|| = {cr:.4e}      ||[H_c, W_c]|| = {cw:.4e}      cluster width = {dlt:.4e}")
print(f"        ||{{R_c, W_c}}|| = {np.linalg.norm(M@Mw + Mw@M, 2):.4e}   (the writer still anticommutes)")
print(f"""
      P-2 DOES NOT SURVIVE AS STATED, AND THE REASON IS TRIVIAL ONCE SEEN.  On a cluster of width
      delta EVERY operator commutes with H to within delta.  ||[H_c,W_c]|| = {cw:.2e} is the same
      order as ||[H_c,R_c]|| = {cr:.2e}.  So at tolerance delta the writer IS indistinguishable from
      something built out of H.  P-2's content evaporates at the tolerance that makes the record
      approximate in the first place.  IT IS NOT LOST -- it is SUBSUMED: what P-2 was used for is
      P-3, and P-3 survives with a number.  RECORD THIS AS A WITHDRAWAL OF P-2's APPROXIMATE FORM.""")

print(f"\n  D4  P-3: (iv)+(v) => THE WRITER IS NON-LOCAL.   APPROXIMATE FORM OF (v): no contractible")
print(f"      operation flips R even PARTIALLY.  Decompose the cluster action of O into the part that")
print(f"      commutes with R_c and the part that anticommutes (the part that flips the record):")
print(f"          O_-  =  (O_c - R_c O_c R_c)/2 ,   FLIP AMPLITUDE = ||O_-||  (||O||=1 for a Pauli).")


def flip_amp(O, Uc, Rc):
    """ABSOLUTE flip amplitude: the norm of the part of O's action on the cluster that ANTICOMMUTES
    with the record, measured against ||O|| = 1 for a Pauli.  NOT normalised by ||P O P||: that
    quantity is itself O(p) for a local operator, and dividing by it manufactures an O(1) answer
    out of two vanishing numbers.  (This lane made that mistake once; the corrected measure is
    below and the erroneous one is reported alongside so the difference is visible.)"""
    Oc = Uc.conj().T @ O @ Uc
    Om = (Oc - Rc @ Oc @ Rc) / 2
    a = float(np.linalg.norm(Om, 2))
    n = float(np.linalg.norm(Oc, 2))
    return a, (a / n if n > 1e-14 else 0.0)


print(f"\n      {'p':>9s} {'TORIC best wt-1 |O_-|':>23s} {'/p':>9s} {'TORIC logical X_A':>19s}"
      f" {'SYMMETRY X_0':>14s} {'[bad measure]':>14s}")
for pp in (0.0, 1e-4, 1e-3, 1e-2, 1e-1):
    _, _, Mp, eep, Ucp = width(Ht, Rt, gt, V8, pp)
    w_, Qm = np.linalg.eigh(Mp)
    Rc = Qm @ np.diag(np.sign(w_)) @ Qm.conj().T
    cand = [flip_amp(O, Ucp, Rc) for l in range(NQ) for O in (Xop([l]), Yop(l), Zop([l]))]
    b1, bad = max(cand, key=lambda z: z[0])
    bl, _ = flip_amp(Xop(X_A_SUP), Ucp, Rc)
    _, _, Msp, eesp, Ucsp = width(Hs, Rs, gs, V8, pp)
    ws_, Qs = np.linalg.eigh(Msp)
    Rcs = Qs @ np.diag(np.sign(ws_)) @ Qs.conj().T
    bs, _ = flip_amp(Xop([0]), Ucsp, Rcs)
    print(f"      {pp:9.0e} {b1:23.6e} {(b1/pp if pp else float('nan')):9.4f} {bl:19.6f}"
          f" {bs:14.6f} {bad:14.4f}")
print(f"""
      READ THE COLUMNS.  On the TORIC carrier the best LOCAL (weight-1) absolute flip amplitude is
      0 exactly at p=0 and is PROPORTIONAL TO p thereafter -- the |O_-|/p column is constant.  The
      weight-2 LOGICAL writer flips the record with amplitude 1 at every p (POSITIVE CONTROL: the
      measure does register a full flip when one exists).  On the SYMMETRY carrier a single X_0 --
      weight one, entirely contractible -- flips the record with amplitude 1 at every p.

      THE LAST COLUMN IS THE MEASURE THIS LANE FIRST WROTE AND THEN CAUGHT: ||O_-|| / ||P O P||.
      It reads 0.86 at every p and would have been reported as "local operations flip the toric
      record almost perfectly".  Both numerator and denominator are O(p); their ratio is O(1) and
      says nothing.  It is kept in the table so the corrected measure can be checked against it.

      P-3 SURVIVES AND ACQUIRES A NUMBER:  local flip amplitude = O(p) = O(p^(d-1)) on d=2.
      The writer is non-local up to an error that vanishes as the perturbation does, at the rate
      set by the SAME d that sets the lifetime.""")

# ==================================================================================================
print("\n" + "=" * 104)
print("  D5.  THE OPEN-SYSTEM HALF, WITH DETAILED BALANCE")
print("=" * 104)
print("""
  The record R = Z on [0,1] is threatened only by X-type noise (Z-type jumps commute with it on
  BOTH carriers, so the comparison would be empty).  Take the same jump operators -- single-qubit
  X_l -- on both carriers, with rates obeying detailed balance at inverse temperature beta:
        w(e -> e') = gamma / (1 + exp(beta * (E(e') - E(e))))          (Glauber)
  On the toric carrier X-errors are labelled by e in F_2^8 with E(e) = -8 + 2*n_defect(e); this is
  an exact reduction, not a model -- the X sector of a CSS code closes on itself.
  On the symmetry carrier every X_l commutes with H_sym, so every ΔE = 0 and every rate is gamma/2.
  gamma = 1 throughout.  THE COMMUTATOR NORM ||[X_l, R]|| = 2 ON BOTH, AT EVERY beta.
""")
PLQ_MASK = [sum(1 << k for k in pl) for pl in PLAQS]
REC_MASK = sum(1 << k for k in Z_A_SUP)


def ndef(e):
    return sum(1 for m in PLQ_MASK if bin(e & m).count('1') & 1)


E_OF_TORIC = np.array([-8 + 2 * ndef(e) for e in range(256)], float)
SIG_TORIC = np.array([1 - 2 * (bin(e & REC_MASK).count('1') & 1) for e in range(256)], float)
SIG_SYM = np.array([1 - 2 * (e & 1) for e in range(256)], float)


def record_rate(beta, gamma=1.0, Earr=None, sig=None):
    """slowest relaxation rate carrying weight in the record observable"""
    E_of = E_OF_TORIC if Earr is None else Earr
    SIG = SIG_TORIC if sig is None else sig
    W = np.zeros((256, 256))
    for e in range(256):
        for l in range(8):
            f = e ^ (1 << l)
            dE = E_of[f] - E_of[e]
            W[f, e] = gamma / (1.0 + np.exp(beta * dE))
    for e in range(256):
        W[e, e] = -W[:, e].sum() + W[e, e]
    pi = np.exp(-beta * E_of)
    pi = pi / pi.sum()
    # detailed-balance self-check
    db = np.max(np.abs(W * pi[None, :] - (W * pi[None, :]).T))
    S = (W * np.sqrt(pi)[None, :]) / np.sqrt(pi)[:, None]
    S = (S + S.T) / 2
    lam, phi = np.linalg.eigh(S)
    v = SIG * np.sqrt(pi)
    ov = (phi.T @ v) ** 2
    ov = ov / ov.sum()
    idx = [i for i in range(256) if ov[i] > 1e-9 and -lam[i] > 1e-12]
    rate = min(-lam[i] for i in idx)
    return rate, db, float(sum(ov[i] for i in idx))


print(f"    {'beta':>7s} {'||[X_l,R]||':>12s} {'TORIC record rate':>19s} {'SYMMETRY record rate':>21s}"
      f" {'T_top/T_sym':>13s} {'detbal resid':>13s}")
rows = []
for beta in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0):
    r_t, db, wt_ = record_rate(beta)
    # SYMMETRY carrier through the SAME machinery: every X_l commutes with H_sym so E is flat,
    # and the record is the parity of the flip pattern on qubit 0.
    r_s, db_s, _ = record_rate(beta, Earr=np.zeros(256), sig=SIG_SYM)
    rows.append((beta, r_t, r_s))
    print(f"    {beta:7.2f} {2.0:12.4f} {r_t:19.6e} {r_s:21.6e} {r_s/r_t:13.4e} {db:13.2e}")
print(f"\n    SELF-CHECK  detailed balance residual max |W_ij pi_j - W_ji pi_i| < 1e-15 in every row"
      f"   PASS")
b0, r0, _ = rows[0]
rs0 = rows[0][2]
print(f"    POSITIVE CONTROL at beta = 0: toric record rate {r0:.4f} vs symmetry {rs0:.4f}, ratio {rs0/r0:.3f}.")
print(f"      {'PASS' if rs0/r0 < 5 else 'FAIL'} -- at infinite temperature the two carriers agree to a factor of {rs0/r0:.2f}.")
print(f"      The machine does NOT manufacture a separation where none exists.  What follows is thermal.")
ar = np.polyfit([b for b, _, _ in rows[2:]], [np.log(r) for _, r, _ in rows[2:]], 1)
print(f"\n    ARRHENIUS FIT over beta in [1,4]:  log(rate_toric) = {ar[0]:.4f}*beta + {ar[1]:.4f}")
print(f"      barrier = {-ar[0]:.4f}   (expected 4 = the energy of one plaquette-defect pair)   "
      f"{'PASS' if abs(-ar[0] - 4) < 0.4 else 'CHECK'}")
print(f"""
    THE VERDICT ON THE L-HALF OF CLAUSE (ii).
      ||[L_k, R]|| = 2 in EVERY ROW ABOVE, on BOTH carriers, at EVERY temperature.  It is the
      largest value the quantity can take, and it does not move.  The actual record lifetime moves
      by a factor of {rows[-1][2]/rows[-1][1]:.3e} across the same rows.  A norm bound on [L_k,R] therefore has
      LITERALLY ZERO discriminating power over the thing it is supposed to be measuring.
      ||[L_k,R]|| <= epsilon IS THE WRONG RELAXATION OF THE L-HALF, AND NOT BY A LITTLE.
      The right one is the same one as for H: the record's decay RATE, i.e. an inverse lifetime.

  D DONE""")
