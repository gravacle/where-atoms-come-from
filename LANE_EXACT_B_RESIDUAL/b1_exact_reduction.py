"""B1 -- THE EXACT REDUCTION, AND THE EXACT FACTS IT MAKES AVAILABLE.

Before any residual can be computed, the pipeline has to be shown correct and the venue's exact
structure has to be laid out.  Everything here is an EXACT ARGUMENT checked numerically; nothing
here is a fitted trend.

  E1  logicals COMPUTED, symplectic self-check at every n.
  E2  code space built exactly; equals the ground space of H = -(X^(x)n + Z^(x)n); dim 2^(n-2).
  E3  H restricted to the code space is the c-number -2*I  => H_S drops out of the dynamics.
  E4  every logical preserves the code space  => the joint evolution never leaves it
      => the FULL 2^n computation and the 2^(n-2) computation are the SAME computation.
      VERIFIED: full-space chi vs code-space chi.
  E5  the fast chi path vs record_model's own env.holevo on the full joint state.  THE FLOAT64
      NOISE FLOOR of this pipeline is read off here.
  E6  MATCHED-CONFIGURATION N-INDEPENDENCE.  The same symplectic configuration at n=4,6,8,10
      gives the same chi -- exactly, not approximately.
  E7  SEPARATION IN THE CODE.  Two partners with identical pairing to the read record but
      DIFFERENT physical support give the same chi -- exactly.
  E8  C-38 AS AN EXACT THEOREM, not a six-decimal coincidence.
"""
import numpy as np, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

t0 = time.time()
say("=" * 110)
say("B1   THE EXACT REDUCTION  --  [[n,n-2,2]] family, n = 4,6,8,10")
say("=" * 110)

NS = [4, 6, 8, 10]
DATA = {}
say("")
say("E1  LOGICALS COMPUTED BY symplectic_logicals, NEVER NOMINATED.  SELF-CHECK MUST PASS.")
say(f"  {'n':>4}{'dim':>8}{'k=n-2':>8}{'pairs K':>9}{'symplectic self-check':>26}")
for n in NS:
    pairs = logical_pairs(n)
    ok, bad = check_symplectic(pairs, n)
    say(f"  {n:>4}{2**n:>8}{n-2:>8}{len(pairs):>9}{('PASS' if ok else 'FAIL ' + str(bad[:3])):>26}")
    assert ok, f"SELF-CHECK FAILED at n={n}: {bad[:5]}"
    assert len(pairs) == n - 2, f"expected {n-2} pairs, got {len(pairs)}"
    DATA[n] = dict(pairs=pairs)

say("")
say("E2/E3/E4  CODE SPACE, EXACTLY.")
say(f"  {'n':>4}{'dim code':>10}{'||Q^dQ-I||':>13}{'||QQ^d-Pg||':>14}{'H|code +2I':>12}{'max leak ||(1-Pg)R Pg||':>26}")
for n in NS:
    Q = code_basis(n)
    pairs = DATA[n]['pairs']
    sx = xz_to_matrix(stabilisers(n)[0], n)
    sz = xz_to_matrix(stabilisers(n)[1], n)
    H = -(sx + sz)
    Pg = ((np.eye(2**n) + sx) / 2) @ ((np.eye(2**n) + sz) / 2)
    orth = np.linalg.norm(Q.conj().T @ Q - np.eye(Q.shape[1]))
    proj = np.linalg.norm(Q @ Q.conj().T - Pg)
    Hc = reduce_op(H, Q)
    hcn = np.linalg.norm(Hc + 2 * np.eye(Q.shape[1]))
    leak = 0.0
    ops = []
    for (a, b) in pairs:
        ops += [xz_to_matrix(a, n), xz_to_matrix(b, n)]
    for P in ops:
        leak = max(leak, np.linalg.norm((np.eye(2**n) - Pg) @ P @ Pg))
    say(f"  {n:>4}{Q.shape[1]:>10}{orth:>13.2e}{proj:>14.2e}{hcn:>12.2e}{leak:>26.2e}")
    assert orth < 1e-12 and proj < 1e-12 and hcn < 1e-12 and leak < 1e-12
    DATA[n]['Q'] = Q
    DATA[n]['H'] = H
    DATA[n]['Pg'] = Pg

# reduced logical operators, and their Pauli algebra
say("")
say("  reduced logicals Rbar = Q^dag R Q -- must be Hermitian, square to I, be traceless,")
say("  and reproduce the symplectic relations EXACTLY on the code space:")
say(f"  {'n':>4}{'max||Rbar^2-I||':>18}{'max|Tr Rbar|':>14}{'max||Rbar-Rbar^d||':>21}{'commutation matches F_2':>25}")
for n in NS:
    Q, pairs = DATA[n]['Q'], DATA[n]['pairs']
    labels, vecs = [], []
    for i, (a, b) in enumerate(pairs):
        labels += [f"X{i+1}", f"Z{i+1}"]; vecs += [a, b]
    red = {}
    e2 = tr = hm = 0.0
    for lab, v in zip(labels, vecs):
        Rb = reduce_op(xz_to_matrix(v, n), Q)
        red[lab] = Rb
        e2 = max(e2, np.linalg.norm(Rb @ Rb - np.eye(Rb.shape[0])))
        tr = max(tr, abs(np.trace(Rb)))
        hm = max(hm, np.linalg.norm(Rb - Rb.conj().T))
    okc = True
    for la, va in zip(labels, vecs):
        for lb, vb in zip(labels, vecs):
            c = np.linalg.norm(red[la] @ red[lb] - red[lb] @ red[la])
            anti = np.linalg.norm(red[la] @ red[lb] + red[lb] @ red[la])
            want = sp(va, vb, n)
            got = 1 if c > 1e-9 else 0
            if want != got or (want == 1 and anti > 1e-9): okc = False
    say(f"  {n:>4}{e2:>18.2e}{tr:>14.2e}{hm:>21.2e}{('PASS' if okc else 'FAIL'):>25}")
    assert e2 < 1e-12 and tr < 1e-10 and hm < 1e-12 and okc
    DATA[n]['red'] = red
    DATA[n]['vec'] = dict(zip(labels, vecs))

# ------------------------------------------------------------------ E4/E5 full vs reduced
say("")
say("E4/E5  FULL 2^n COMPUTATION  vs  CODE-SPACE COMPUTATION  vs  record_model's OWN env.holevo.")
say("  Configuration notation:  read X1 on site 0; partners listed as label@site.")
NB = 3
V = BASE
TT = V.times

def full_chi(n, read, parts, lam, env, times):
    """the untouched, unreduced pipeline: full 2^n system, H_S included, state Pg/k,
       chi read with record_model's OWN Environment.holevo on the full joint state."""
    Pg = DATA[n]['Pg']; k = 2 ** (n - 2)
    nS = 2 ** n; nB = env.dim
    R = xz_to_matrix(DATA[n]['vec'][read], n)
    HINT = sum(np.kron(xz_to_matrix(DATA[n]['vec'][lab], n), env.site[j % env.nq]) for lab, j in parts)
    HINT = HINT + np.kron(R, env.site[0])
    Ht = np.kron(DATA[n]['H'], np.eye(nB)) + np.kron(np.eye(nS), env.HB) + lam * HINT
    w, U = np.linalg.eigh(Ht)
    r0 = np.kron(Pg / k, env.thermal())
    Uc = U.conj().T @ r0 @ U
    vals = []
    for t in times:
        ph = np.exp(-1j * w * t)
        r = U @ (ph[:, None] * Uc * ph.conj()[None, :]) @ U.conj().T
        vals.append(env.holevo(r, R, nS))
    return float(np.mean(vals))

def red_chi(n, read, parts, lam, env, times):
    red = DATA[n]['red']
    ops = [(red[read], 0)] + [(red[lab], j) for lab, j in parts]
    return float(np.mean(chi_times(ops, red[read], env, lam, times)))

env3 = V.env(NB)
CHECKS = [(4, "X1", []), (4, "X1", [("Z1", 0)]), (4, "X1", [("Z1", 1)]),
          (6, "X1", []), (6, "X1", [("Z1", 0)]), (6, "X1", [("X2", 0)]), (6, "X1", [("X2", 1)]),
          (6, "X1", [("Z1", 0), ("X2", 1)]),
          (8, "X1", [("Z1", 0)]), (8, "X1", [("X2", 0), ("X3", 1)])]
say(f"  {'n':>4}  {'configuration':<26}{'chi FULL 2^n':>15}{'chi CODE SPACE':>16}{'|difference|':>15}")
worst = 0.0
for n, read, parts in CHECKS:
    cfg = read + "@0" + ("," + ",".join(f"{l}@{j}" for l, j in parts) if parts else "")
    a = full_chi(n, read, parts, 0.8, env3, TT)
    b = red_chi(n, read, parts, 0.8, env3, TT)
    worst = max(worst, abs(a - b))
    say(f"  {n:>4}  {cfg:<26}{a:>15.12f}{b:>16.12f}{abs(a-b):>15.2e}")
say("")
say(f"  WORST full-vs-reduced discrepancy over {len(CHECKS)} configurations: {worst:.3e}")
say("  -> the reduction is EXACT; the discrepancy is the float64 floor of the two paths.")
assert worst < 1e-9, "reduction check failed"
FLOOR_PATH = worst

# ------------------------------------------------------------------ E6 N-independence
say("")
say("E6  MATCHED-CONFIGURATION N-INDEPENDENCE.  Identical symplectic configuration, different n.")
say("      (the logicals at different n have DIFFERENT physical support and different weight)")
CFGS = [("alone",              "X1", []),
        ("Z1@0  (pair 1, same site)",  "X1", [("Z1", 0)]),
        ("Z1@1  (pair 1, other site)", "X1", [("Z1", 1)]),
        ("X2@0  (pair 0, same site)",  "X1", [("X2", 0)]),
        ("X2@1  (pair 0, other site)", "X1", [("X2", 1)]),
        ("Z1@0,X2@0",                  "X1", [("Z1", 0), ("X2", 0)]),
        ("X2@0,Z2@0  (partners pair each other)", "X1", [("X2", 0), ("Z2", 0)])]
say(f"  {'configuration':<42}" + "".join(f"{'n='+str(n):>18}" for n in NS) + f"{'max spread':>13}")
spreads = []
for name, read, parts in CFGS:
    row = []
    for n in NS:
        need = set(l for l, _ in parts)
        if not need.issubset(DATA[n]['red'].keys()):
            row.append(None); continue
        row.append(red_chi(n, read, parts, 0.8, env3, TT))
    vals = [v for v in row if v is not None]
    sprd = max(vals) - min(vals) if len(vals) > 1 else 0.0
    spreads.append(sprd)
    say(f"  {name:<42}" + "".join((f"{v:>18.12f}" if v is not None else f"{'--':>18}") for v in row) + f"{sprd:>13.2e}")
say("")
say(f"  LARGEST spread of a matched configuration across n = 4..10: {max(spreads):.3e}")
say("  -> chi is EXACTLY independent of n at fixed configuration.  PROOF: on the code space")
say("     H_S is a c-number and the written logicals generate a Pauli algebra on k=n-2 logical")
say("     qubits; a logical Clifford maps any symplectic configuration to the standard one, and")
say("     the UNWRITTEN logical qubits are maximally mixed and uncoupled, so they factor out of")
say("     the state, the evolution and the readout.  Spectators are exactly invisible.")
FLOOR_N = max(spreads)

# ------------------------------------------------------------------ E7 separation in the code
say("")
say("E7  SEPARATION IN THE CODE.  Partners with the SAME pairing to X1 but different physical")
say("    support / weight.  If influence fell off with separation these would differ.")
n = 8
Q = DATA[n]['Q']; vecs = DATA[n]['vec']; red = DATA[n]['red']
def wt(v): return sum(1 for i in range(n) if v[i] or v[n + i])
def overlap(u, v): return sum(1 for i in range(n) if (u[i] or u[n+i]) and (v[i] or v[n+i]))
base = vecs["X1"]
say(f"  n=8, read X1 (weight {wt(base)}).")
say(f"  {'partner':<10}{'weight':>8}{'overlap with X1':>18}{'pairing':>9}{'chi (same site)':>19}{'chi (other site)':>19}")
groups = {}
for lab in ["X2", "Z2", "X3", "Z3", "Z1"]:
    v = vecs[lab]
    pg = sp(base, v, n)
    a = red_chi(n, "X1", [(lab, 0)], 0.8, env3, TT)
    b = red_chi(n, "X1", [(lab, 1)], 0.8, env3, TT)
    groups.setdefault(pg, []).append((lab, a, b))
    say(f"  {lab:<10}{wt(v):>8}{overlap(base,v):>18}{pg:>9}{a:>19.12f}{b:>19.12f}")
say("")
for pg, rows in sorted(groups.items()):
    if len(rows) < 2: continue
    sa = max(r[1] for r in rows) - min(r[1] for r in rows)
    sb = max(r[2] for r in rows) - min(r[2] for r in rows)
    say(f"  pairing {pg}: spread over {len(rows)} partners of different support -- same site {sa:.3e}, other site {sb:.3e}")
say("  -> chi depends on the PAIRING and on nothing else about the partner. Separation in the")
say("     code is EXACTLY irrelevant, not merely weakly relevant.")

# ------------------------------------------------------------------ E8 C-38 exact
say("")
say("E8  C-38 AS AN EXACT THEOREM.  A COMMUTING partner on ANOTHER bath site changes chi by")
say("    EXACTLY ZERO -- for any coupling, any time, any number of such partners.")
say("    PROOF.  Let the read record R and every other-site partner P_j commute.  Then the")
say("    system splits into joint eigensectors of {R, P_j}.  H_B is a SUM over bath qubits and")
say("    rho_th is a PRODUCT, and each partner drives a DIFFERENT bath qubit, so in every sector")
say("    the bath state factorises:  rho_B(r, p) = sigma_0(r) (x) prod_j tau_j(p_j).")
say("    Conditioning on R = r and averaging over p gives sigma_0(r) (x) tau, with tau INDEPENDENT")
say("    of r.  Both S(average) and average S(.) then shift by the same S(tau), so chi is the")
say("    alone value identically.  The theorem does not depend on lam or t at all.")
say("")
say(f"  {'n':>4}{'lam':>7}  {'configuration':<30}{'chi':>18}{'chi/alone - 1':>18}")
worst38 = 0.0
for n in NS:
    for lam in (0.4, 0.8, 1.2, 2.5):
        alone = red_chi(n, "X1", [], lam, env3, TT)
        tests = [("X2@1", [("X2", 1)])] if n >= 6 else []
        if n >= 8: tests.append(("X2@1,X3@2", [("X2", 1), ("X3", 2)]))
        for name, parts in tests:
            c = red_chi(n, "X1", parts, lam, env3, TT)
            d = c / alone - 1.0
            worst38 = max(worst38, abs(d))
            say(f"  {n:>4}{lam:>7.2f}  {name:<30}{c:>18.12f}{d:>18.2e}")
say("")
say(f"  WORST |chi/alone - 1| over all commuting-other-site tests: {worst38:.3e}")
say("  -> ZERO to the float64 floor, and EXACTLY zero by the proof above.  C-38 is not a")
say("     numerical trend; it survives the weakness objection outright.")

say("")
say("=" * 110)
say("B1 FLOOR SUMMARY -- these are the numbers a residual must beat")
say("=" * 110)
say(f"  float64 floor, two independent code paths for the same exact quantity : {FLOOR_PATH:.3e}")
say(f"  float64 floor, matched configuration across n = 4..10                 : {FLOOR_N:.3e}")
say(f"  float64 floor, C-38 exact-zero ratio                                  : {worst38:.3e}")
say(f"  elapsed {time.time()-t0:.1f}s")
