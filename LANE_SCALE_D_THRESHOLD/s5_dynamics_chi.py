"""S5 -- THE DYNAMICAL ARM: TIME-AVERAGED HOLEVO chi vs THE NUMBER OF RECORDS, EXACTLY.

WHY AN ANALYTIC SOLVER IS ALLOWED HERE.  The couplings used in this program's density work are
the records themselves: H_int = lam * sum_i R_i (x) X_(site i).  The R_i all commute with H and
with each other, so from the maximally mixed code state the joint state stays DIAGONAL in the
joint record basis: it is a classical mixture over sign strings a in {+-1}^k of
    |a><a| (x) rho_B(a),      rho_B(a) = product over bath sites j of rho_j(c_j(a)),
    c_j(a) = sum of a_i over the records assigned to site j.
Each bath site therefore evolves as a SINGLE QUBIT under  e_j Z + lam c_j X.  chi about record 1
involves only the site record 1 sits on, because every other site's state is independent of a_1
and factors out of both entropies.  The whole dynamical question collapses to an EXACT
one-qubit computation at ANY k -- so this arm is not limited by Hilbert-space dimension.

THAT CLAIM IS NOT ASSUMED.  Part 1 below validates the analytic solver against the model's own
RecordModel.evolve + Environment.holevo on the real Hilbert space at dim 16, 64 and 256, on
records verified by S1.  If the validation fails, the sweep is not reported.
"""
import sys, time, json
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD")
import numpy as np
from record_model import RecordModel, Environment, xz_to_matrix
import carriers as C

OUT = []
def P(s=""):
    print(s); OUT.append(s)

Z = np.array([[1, 0], [0, -1]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)

def vN(r):
    e = np.linalg.eigvalsh(r); e = e[e > 1e-13]
    return float(-(e * np.log2(e)).sum())

def site_rho(e, lam, c, beta, t):
    """one bath qubit: thermal in e*Z, then evolved under e*Z + lam*c*X"""
    p = np.exp(-beta * np.array([e, -e])); p = p / p.sum()
    r0 = np.diag(p).astype(complex)
    Hh = e * Z + lam * c * X
    w, V = np.linalg.eigh(Hh)
    U = (V * np.exp(-1j * w * t)) @ V.conj().T
    return U @ r0 @ U.conj().T

def chi_site(q, e, lam, beta, t, about="one"):
    """EXACT chi at a site carrying q records, from the maximally mixed code state.
       about='one'      -> chi(R_1 : bath)
       about='register' -> chi(whole q-record register : bath)

       Only c = sum of the q signs enters the bath Hamiltonian, so the 2^q sign strings are
       collapsed to BINOMIAL MULTIPLICITIES -- exact, and it is what makes k = 40 reachable.
       (Enumerating the strings is what capped the first version of this sweep at k ~ 22.)"""
    from math import comb
    cs = {}
    for j in range(q + 1):
        cs[q - 2 * j] = cs.get(q - 2 * j, 0) + comb(q, j)
    states = {c: site_rho(e, lam, c, beta, t) for c in cs}
    tot = float(2 ** q)
    if about == "register":
        av = sum((m / tot) * states[c] for c, m in cs.items())
        cond = sum((m / tot) * vN(states[c]) for c, m in cs.items())
        return max(vN(av) - cond, 0.0)
    halves = []
    for s in (1, -1):
        cnt = {}
        for j in range(q):
            c = s + (q - 1 - 2 * j)
            cnt[c] = cnt.get(c, 0) + comb(q - 1, j)
        h = float(2 ** (q - 1))
        for c in cnt:
            if c not in states: states[c] = site_rho(e, lam, c, beta, t)
        halves.append(sum((m / h) * states[c] for c, m in cnt.items()))
    av = 0.5 * (halves[0] + halves[1])
    return max(vN(av) - 0.5 * (vN(halves[0]) + vN(halves[1])), 0.0)

TIMES = np.linspace(1.0, 13.0, 25)
def chi_avg(q, e, lam, beta, about="one"):
    return float(np.mean([chi_site(q, e, lam, beta, t, about=about) for t in TIMES]))

# ---------------------------------------------------------------- PART 1: VALIDATION
t0 = time.time()
P("=" * 150)
P("S5  TIME-AVERAGED HOLEVO chi vs NUMBER OF RECORDS")
P("=" * 150)
P()
P("PART 1  VALIDATION OF THE ANALYTIC SOLVER AGAINST RecordModel.evolve + Environment.holevo")
P("  full Hilbert space, records verified in S1, 25 times averaged in [1,13].  If any row's")
P("  |analytic - model| is not ~1e-12 the sweep below is NOT reported.")
P("  %-14s %-4s %-4s %-9s %-6s %-14s %-14s %-12s" %
  ("carrier", "k", "nq", "layout", "lam", "chi_model", "chi_analytic", "|diff|"))
P("  " + "-" * 96)
val_rows = []
worst = 0.0
for n, nq, layout in ((4, 3, "crowded"), (4, 3, "spread"), (6, 3, "crowded"),
                      (6, 3, "spread"), (8, 3, "crowded")):
    car = C.family_A(n)
    fam, part, ch = C.records_of(car)
    assert C.all_checks_pass(ch)
    k = len(fam)
    Hm = -sum(xz_to_matrix(s, n) for s in car["stabs"])
    M = RecordModel(Hm)
    Rs = [xz_to_matrix(v, n) for v in fam]
    energies = tuple(1.0 + 0.0 * j for j in range(nq))    # all sites identical -> analytic e
    env = Environment(nq=nq, energies=energies, beta=2.0)
    e0, beta, lam = energies[0], 2.0, 0.8
    if layout == "crowded":
        assign = [0] * k
    else:
        assign = [i % nq for i in range(k)]
    coupling = [(Rs[i], assign[i]) for i in range(k)]
    # ONE eigendecomposition, all 25 times (the model redoes eigh per call; this does not)
    nS, nB = M.n, env.dim
    HINT = sum(np.kron(A, env.site[j]) for A, j in coupling)
    Ht = np.kron(Hm, np.eye(nB)) + np.kron(np.eye(nS), env.HB) + lam * HINT
    w, U = np.linalg.eigh(Ht)
    Pg, kg = M.ground_space()
    r0 = np.kron(Pg / kg, env.thermal())
    Uc = U.conj().T @ r0 @ U
    chis = []
    for t in TIMES:
        ph = np.exp(-1j * w * t)
        r = U @ (ph[:, None] * Uc * ph.conj()[None, :]) @ U.conj().T
        chis.append(env.holevo(r, Rs[0], nS))
    cm = float(np.mean(chis))
    q = sum(1 for a in assign if a == assign[0])
    ca = chi_avg(q, e0, lam, beta)
    P("  %-14s %-4d %-4d %-9s %-6.2f %-14.10f %-14.10f %-12.2e"
      % (car["label"], k, nq, layout, lam, cm, ca, abs(cm - ca)))
    worst = max(worst, abs(cm - ca))
    val_rows.append((car["label"], k, nq, layout, cm, ca))
    P("     (t=%.0fs)" % (time.time() - t0))
P("  worst |analytic - model| = %.3e  ->  %s" % (worst, "VALIDATED" if worst < 1e-9 else "FAILED"))
VALID = worst < 1e-9

# also a SANITY CONTROL: a coupling that opens NO channel must give chi = 0
car = C.family_A(4); fam, part, ch = C.records_of(car)
Hm = -sum(xz_to_matrix(s, 4) for s in car["stabs"])
M = RecordModel(Hm); Rs = [xz_to_matrix(v, 4) for v in fam]
env = Environment(nq=3, energies=(1.0, 1.0, 1.0), beta=2.0)
ident = np.eye(16, dtype=complex)
r = M.evolve([(ident, 0)], env, lam=0.8, t=4.0)
P("  NEGATIVE CONTROL: coupling = identity (opens no channel) -> chi = %.3e   [must be 0]"
  % env.holevo(r, Rs[0], 16))
r = M.evolve([(Rs[0], 0)], env, lam=0.8, t=4.0)
P("  POSITIVE CONTROL: coupling = R_1                        -> chi = %.6f   [must be > 0]"
  % env.holevo(r, Rs[0], 16))

# ---------------------------------------------------------------- PART 2: THE SWEEP
if not VALID:
    P(); P("VALIDATION FAILED -- NO SWEEP REPORTED, NO CONCLUSION DRAWN.")
else:
    P(); P("PART 2  chi vs k, EXACT, at three venue scales (D-17) and with the SPREAD control")
    P("  in the same table (D-16: suppression is measured against SPREAD, never against ALONE).")
    P("  CROWDED = all k records on ONE bath site.  SPREAD = one record per site (q=1).")
    P("  PAIRED  = family B's structure: 2 records per site, whatever k is.")
    rows = []
    for lam in (0.4, 0.8, 1.2):
        for e0 in (1.0, 0.5):
            for beta in (2.0, 0.5):
                sp = chi_avg(1, e0, lam, beta)
                pr = chi_avg(2, e0, lam, beta)
                for k in list(range(1, 25)) + [28, 32, 40, 48, 64]:
                    cr = chi_avg(k, e0, lam, beta)
                    reg = chi_avg(k, e0, lam, beta, about="register")
                    rows.append(dict(lam=lam, e=e0, beta=beta, k=k, crowded=cr, spread=sp,
                                     paired=pr, ratio=(cr / sp if sp > 0 else float("nan")),
                                     register=reg))
    P()
    P("  %-6s %-5s %-6s %-4s %-13s %-13s %-13s %-13s %-13s" %
      ("lam", "e", "beta", "k", "chi_CROWDED", "chi_SPREAD", "chi_PAIRED", "crowded/spread",
       "chi_register"))
    P("  " + "-" * 104)
    for r in rows:
        P("  %-6.2f %-5.2f %-6.2f %-4d %-13.8f %-13.8f %-13.8f %-13.8f %-13.8f" %
          (r["lam"], r["e"], r["beta"], r["k"], r["crowded"], r["spread"], r["paired"],
           r["ratio"], r["register"]))
    json.dump(rows, open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s5_rows.json", "w"))

    P()
    P("PART 3  MARGINAL COST OF THE k-TH RECORD, and its SIGN -- a sign change would be a")
    P("  qualitative threshold.  d1(k) = chi_crowded(k) - chi_crowded(k-1);")
    P("  dR(k) = chi_register(k) - chi_register(k-1).")
    P("  %-6s %-5s %-6s %-4s %-14s %-14s %-9s %-9s" %
      ("lam", "e", "beta", "k", "d1(k)", "dR(k)", "sgn d1", "sgn dR"))
    P("  " + "-" * 76)
    for lam in (0.4, 0.8, 1.2):
        for e0 in (1.0,):
            for beta in (2.0,):
                prev1 = prevR = None
                for k in range(1, 33):
                    c1 = chi_avg(k, e0, lam, beta)
                    cR = chi_avg(k, e0, lam, beta, about="register")
                    if prev1 is not None:
                        d1 = c1 - prev1; dR = cR - prevR
                        P("  %-6.2f %-5.2f %-6.2f %-4d %-14.8f %-14.8f %-9s %-9s" %
                          (lam, e0, beta, k, d1, dR,
                           "+" if d1 > 1e-12 else ("-" if d1 < -1e-12 else "0"),
                           "+" if dR > 1e-12 else ("-" if dR < -1e-12 else "0")))
                    prev1, prevR = c1, cR

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s5_dynamics_chi.txt",
     "w").write("\n".join(OUT) + "\n")
P("total %.1fs" % (time.time() - t0))
