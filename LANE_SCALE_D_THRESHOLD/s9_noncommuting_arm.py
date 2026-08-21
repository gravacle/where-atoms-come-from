"""S9 -- THE ARM THE ANALYTIC SOLVER CANNOT REACH, AND HOW FAR IT ACTUALLY GOES.

S5's exact solver works because the couplings COMMUTE with the records, which is what makes the
joint state stay diagonal.  Extend C-39 instead: keep record R_1 coupled to bath site 0 and put
the OTHER k-1 couplings on the same site, choosing them two ways --

   COMMUTING CROWD  the other records R_2..R_k          (they commute with R_1)
   PAIRING CROWD    the conjugate partners P_1..P_(k-1)  (P_1 ANTICOMMUTES with R_1)

and measure chi(R_1 : bath) time-averaged over 25 times in [1,13].  C-39 measured exactly this
contrast at k = 2..4 and found the pairing crowd cheaper.  Here it is pushed as far as the full
Hilbert space allows.  Nothing is collapsed analytically: this is RecordModel + Environment on
the real joint space, so the reach is set by eigh at dim (2^n * 2^nq) and is reported as such.

D-15: the ALONE column and the NO-CHANNEL negative control are in the same table.
D-16: suppression is reported against SPREAD, not against ALONE.
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

TIMES = np.linspace(1.0, 13.0, 25)

def chi_series(M, Hm, env, coupling, R, lam=0.8):
    """ONE eigendecomposition, all 25 times -- and the record's value displacement."""
    nS, nB = M.n, env.dim
    HINT = sum(np.kron(A, env.site[j]) for A, j in coupling)
    Ht = np.kron(Hm, np.eye(nB)) + np.kron(np.eye(nS), env.HB) + lam * HINT
    w, U = np.linalg.eigh(Ht)
    Pg, kg = M.ground_space(); rho0 = Pg / kg
    Uc = U.conj().T @ np.kron(rho0, env.thermal()) @ U
    chis = []; moved = []
    base = float(np.real(np.trace(rho0 @ R)))
    for t in TIMES:
        ph = np.exp(-1j * w * t)
        r = U @ (ph[:, None] * Uc * ph.conj()[None, :]) @ U.conj().T
        chis.append(env.holevo(r, R, nS))
        rS = r.reshape(nS, nB, nS, nB).trace(axis1=1, axis2=3)
        moved.append(abs(float(np.real(np.trace(rS @ R))) - base))
    return float(np.mean(chis)), float(np.mean(moved))

t0 = time.time()
P("=" * 150)
P("S9  C-39 EXTENDED: chi(R_1) WITH A COMMUTING CROWD vs A PAIRING CROWD ON THE SAME BATH SITE")
P("=" * 150)
P("  full Hilbert space, no analytic collapse.  lam = 0.8, bath nq = 3, beta = 2.0,")
P("  25 times averaged in [1,13].  Records verified in S1.")
P()
P("  %-12s %-4s %-8s %-6s %-13s %-13s %-13s %-13s %-13s %-10s" %
  ("carrier", "k", "dim_tot", "lam", "chi_ALONE", "chi_SPREAD", "chi_COMMUTING",
   "chi_PAIRING", "commut/spread", "pair/spread"))
P("  " + "-" * 126)
rows = []
for n in (4, 6, 8):
    car = C.family_A(n)
    fam, part, ch = C.records_of(car)
    if not C.all_checks_pass(ch):
        P("  SELF-CHECK FAILED on %s -- no conclusion" % car["label"]); continue
    k = len(fam)
    Hm = -sum(xz_to_matrix(s, n) for s in car["stabs"])
    M = RecordModel(Hm)
    Rs = [xz_to_matrix(v, n) for v in fam]
    Ps = [xz_to_matrix(v, n) for v in part]
    env = Environment(nq=3, energies=(1.0, 1.0, 1.0), beta=2.0)
    lam = 0.8
    alone, _ = chi_series(M, Hm, env, [(Rs[0], 0)], Rs[0], lam)
    spread, _ = chi_series(M, Hm, env, [(Rs[i], i % 3) for i in range(k)], Rs[0], lam)
    commut, mv_c = chi_series(M, Hm, env, [(Rs[i], 0) for i in range(k)], Rs[0], lam)
    # C-39's contrast: the crowd must include R_1's OWN conjugate partner P_1, which is the
    # only operator here that does NOT commute with R_1.  Using P_2..P_k instead would make the
    # pairing crowd commute with R_1 and the two columns would be the same measurement twice.
    pairing, mv_p = chi_series(M, Hm, env, [(Rs[0], 0)] + [(Ps[i], 0) for i in range(k - 1)],
                               Rs[0], lam)
    P("  %-12s %-4d %-8d %-6.2f %-13.8f %-13.8f %-13.8f %-13.8f %-13.8f %-10.6f"
      % (car["label"], k, M.n * env.dim, lam, alone, spread, commut, pairing,
         commut / spread, pairing / spread))
    rows.append(dict(car=car["label"], k=k, dim=M.n * env.dim, alone=alone, spread=spread,
                     commuting=commut, pairing=pairing, moved_commuting=mv_c,
                     moved_pairing=mv_p))
    P("     record value displacement |d<R_1>|: commuting crowd %.3e   pairing crowd %.3e"
      % (mv_c, mv_p))
    # T-35: "(t=Ns)" reformatted to "t=Ns" -- the parenthesised form escapes reproduce.sh's norm() and differed every run
    P("     t=%.0fs" % (time.time() - t0))

# D-15 negative control, in the same run
car = C.family_A(4); fam, part, ch = C.records_of(car)
Hm = -sum(xz_to_matrix(s, 4) for s in car["stabs"]); M = RecordModel(Hm)
Rs = [xz_to_matrix(v, 4) for v in fam]
env = Environment(nq=3, energies=(1.0, 1.0, 1.0), beta=2.0)
nc, _ = chi_series(M, Hm, env, [(np.eye(16, dtype=complex), 0)], Rs[0])
pc, _ = chi_series(M, Hm, env, [(Rs[0], 0)], Rs[0])
P()
P("  NEGATIVE CONTROL  coupling opens no channel (identity)  -> chi = %.3e  [must be 0]" % nc)
P("  POSITIVE CONTROL  coupling = R_1                        -> chi = %.6f  [must be > 0]" % pc)
P()
P("REACH OF THIS ARM: the largest carrier here is dim %d joint (system 2^8 x bath 2^3).  The"
  % (256 * 8))
P("  next step up, n = 10, is a joint dimension of 8192; eigh at 4096 already costs minutes on")
P("  this machine and 8192 needs ~1 GB per dense copy.  So the NON-COMMUTING arm stops at k = 6.")
P("  That is the honest ceiling on any dynamical statement in this lane that is not analytic.")

json.dump(rows, open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s9_rows.json", "w"))
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s9_noncommuting_arm.txt",
     "w").write("\n".join(OUT) + "\n")
# T-35: post-write "total Ns" stdout line removed -- it was captured by reproduce.sh but never part of the sealed output, so it differed every run
