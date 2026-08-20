"""S3 -- HOLEVO chi ON THE [[n, n-2, 2]] FAMILY, EXACTLY, TO LARGE k.

REPRESENTATION.  The code space (dimension 2^k) only -- the 2^n dense space is never built.
Inside it the k records Xbar_i are commuting +-1 observables, so the code space is the
classical k-bit register {+-1}^k, and the joint state after evolution is EXACTLY

    rho(t) = 2^-k SUM_s |s><s| (x) U_s(t) rho_B U_s(t)-dag ,
    U_s(t) = exp(-i t (H_B + lam SUM_i s_i b_i)) ,   b_i = the bath site record i couples to.

S1/SC-3 validated this against RecordModel.formation()'s dense 2^n x 2^nq computation.

THE COMPRESSION.  With equal coupling strengths the bath sees only the nq integers
c_j = SUM_{i : b_i = j} s_i, so the 2^k sum collapses to a product of binomial distributions
over O((k/nq)^nq) counts.  S1/SC-4 validated this to 1e-12 against the 2^k enumeration.
This is what makes k = 40 reachable where 2^k enumeration dies at k = 20.

TWO VENUES, and the contrast is the point:
  FIXED BATH   nq = 3 for every N; record i couples to site i mod 3.  Records SHARE sites.
  GROWN BATH   nq = N; record i couples to its OWN site.  The bath grows with the records.

D-15: the GROWN BATH is the positive control that would register linear growth if it were there.
D-17: lam, nq, beta and the time window are all varied before any effect is called new.
TIME-AVERAGE: chi is averaged over 25 times in [1,13]; a fixed-t snapshot recurs.
"""
import sys, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE")
from chi_lib import (vN, energies_for, bath_states, group_dist, chi_suite,
                     total_chi_fixed, total_chi_grown, TIMES, Environment)

OUT = []
def say(s=""):
    print(s); OUT.append(s)

# =====================================================================================
say("="*112)
say("S3   HOLEVO chi ON THE [[n, n-2, 2]] FAMILY.  Code-space representation; counts compression.")
say("     Every number is time-averaged over 25 times in [1,13]; +- is the standard error of that mean.")
say("="*112)
say()

# ---- SC-5: chi_i depends only on i mod nq -- verified, not assumed
say("SC-5  chi_i depends only on the site i couples to.  Verified against a per-record computation.")
say("      (if this failed, the group compression used below would be invalid)")
say()
say("        k   per-group chi (g=0,1,2)                        max spread within a group")
sc5_ok = True
for k in (4, 6, 9):
    per, cj, m = chi_suite(k, 3, 0.8, 2.0, 4.0)
    # independent recomputation record-by-record via full 2^k sectors, small k only
    env = Environment(nq=3, energies=energies_for(3), beta=2.0)
    rB = env.thermal(); acc = {}
    for s in itertools.product((1, -1), repeat=k):
        HB = env.HB + 0.8*sum(s[i]*env.site[i % 3] for i in range(k))
        w, U = np.linalg.eigh(HB); ph = np.exp(-1j*w*4.0)
        Uc = U.conj().T @ rB @ U
        acc[s] = U @ (ph[:, None]*Uc*ph.conj()[None, :]) @ U.conj().T
    chis = []
    for i in range(k):
        P = sum(acc[s] for s in acc if s[i] > 0)/2**(k-1)
        M = sum(acc[s] for s in acc if s[i] < 0)/2**(k-1)
        chis.append(max(vN(0.5*(P+M)) - 0.5*(vN(P)+vN(M)), 0.0))
    spread = max(max(abs(chis[i]-per[i % 3]) for i in range(k)), 0.0)
    sc5_ok &= spread < 1e-10
    say("     %4d   %-45s %.2e" % (k, "  ".join("%.9f" % x for x in per), spread))
say()
say("   SC-5 %s" % ("PASS" if sc5_ok else "FAIL"))
say()

KS = list(range(2, 41, 2))
say("-"*112)
say("TABLE 2.  chi vs N, FIXED 3-qubit bath (saturating venue) beside GROWN bath (linear control).")
say("          lam = 0.8, beta = 2.0, energies (1.0,1.4,0.7) repeating.")
say()
say("      n     N=k   |  SUM_i chi_i (fixed bath)   sd(t)   |  chi_joint (fixed)  |"
    "  SUM_i chi_i (grown bath)  |  SUM/joint")
rows = []
for k in KS:
    tm, tse, tsd, jm, jse = total_chi_fixed(k)
    gm, gse = total_chi_grown(k)
    rows.append((k, tm, tse, tsd, jm, jse, gm, gse))
    say("   %4d %6d   | %14.6f +- %.4f  %6.3f  | %10.6f +- %.4f |  %14.6f +- %.4f  | %8.3f"
        % (k+2, k, tm, tse, tsd, jm, jse, gm, gse, tm/jm if jm > 0 else float('nan')))
say()
say("  EXACT BOUND ON THE FIXED-BATH COLUMN.  The records are independent classical bits, so")
say("      SUM_i chi(R_i : B)  <=  chi(s_1..s_N : B)  =  S(rho_B_bar) - <S(rho_B(s))>  <=  S(rho_B_bar)  <=  nq bits.")
say("  With nq = 3 the cap is 3.000 bits AT EVERY N, by an exact argument, not by extrapolation.")
say("  Observed max of the fixed-bath column over N = 2..40: %.6f bits.  Cap respected: %s"
    % (max(r[1] for r in rows), max(r[1] for r in rows) <= 3.0 + 1e-9))
say("  Observed max of chi_joint: %.6f bits (<= 3): %s"
    % (max(r[4] for r in rows), max(r[4] for r in rows) <= 3.0 + 1e-9))
say("  SUM <= joint holds at every N: %s" % all(r[1] <= r[4] + 1e-9 for r in rows))
say("-"*112)
say()

# ------------------------------------------------------------------ D-17 venue scale sweep
say("-"*112)
say("D-17  VARY THE VENUE'S OWN SCALE BEFORE CALLING ANYTHING NEW.")
say("      Fixed-bath SUM_i chi_i at N = 2, 8, 20, 40 across coupling, bath size, temperature,")
say("      and time window.  If the shape were an artifact of one setting it would move here.")
say()
say("     lam    nq   beta   window      N=2       N=8      N=20      N=40    ratio N=40/N=8   cap(nq bits)")
sweep = []
for (lam, nq, beta, win) in [(0.4, 3, 2.0, (1, 13)), (0.8, 3, 2.0, (1, 13)), (1.6, 3, 2.0, (1, 13)),
                             (0.8, 2, 2.0, (1, 13)), (0.8, 4, 2.0, (1, 13)),
                             (0.8, 3, 0.5, (1, 13)), (0.8, 3, 5.0, (1, 13)),
                             (0.8, 3, 2.0, (5, 40))]:
    tms = np.linspace(win[0], win[1], 25)
    vals = []
    for k in (2, 8, 20, 40):
        tm, *_ = total_chi_fixed(k, nq=nq, lam=lam, beta=beta, times=tms)
        vals.append(tm)
    sweep.append((lam, nq, beta, win, vals))
    say("   %5.2f %5d %6.2f  [%2d,%2d] %9.4f %9.4f %9.4f %9.4f %14.3f %14d"
        % (lam, nq, beta, win[0], win[1], vals[0], vals[1], vals[2], vals[3],
           vals[3]/vals[1] if vals[1] > 0 else float('nan'), nq))
say()
say("  READ: gravity's requirement (a) is S(2N)/S(N) -> 2.  The N=40/N=8 column is a FIVEFOLD")
say("  increase in N; linear extensivity would put that ratio at 5.  Observed range: %.3f to %.3f."
    % (min(v[4][3]/v[4][1] for v in sweep), max(v[4][3]/v[4][1] for v in sweep)))
say("-"*112)
say()

np.save("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE/s3_data.npy", np.array(rows))
say("largest N reached in S3: %d  (n = %d).  What stopped it: nothing in this venue -- the counts" % (max(KS), max(KS)+2))
say("compression is polynomial.  The dense 2^n route stops at n = 8 (2048-dim joint state); the 2^k")
say("sector route stops near k = 20; the compression is what carried it to k = 40 and would carry further.")
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE/s3_chi.txt", "w").write("\n".join(OUT)+"\n")
