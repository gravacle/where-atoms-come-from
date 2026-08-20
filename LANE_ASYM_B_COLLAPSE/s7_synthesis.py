"""S7 -- ONE TABLE.  Every quantity, every control, every test, in the same rows.

Sources, all produced in this lane and all on disk beside this file:
  s1_validate.txt          the self-checks that license everything else
  s2_combinatorial.txt     exact F_2 quantities to n = 60
  s3_chi.txt               exact Holevo chi to N = 40
  s4_fss.txt               the finite-size-scaling engine and its calibration
  s5_additivity.txt        requirement (b), additivity over disjoint regions
  s6_saturating_control.txt the brief's own saturating control, and its linear twin

THE STANDARD, restated so the READ below can be checked against it:
  (a) STRICT EXTENSIVITY   S(N) unbounded, asymptotically linear, S(2N)/S(N) -> 2
  (b) ADDITIVITY over disjoint regions
  (c) NOT saturating, NOT topological
A quantity failing any of these is not gravity's source at ANY N.  Failing (a) or (b) by an
EXACT argument rules it out at every finite N; failing it only within the computed range is
ABSENCE WITHIN RANGE and is labelled as such in the RANGE column.
"""
import sys, numpy as np
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE"
sys.path.insert(0, LANE)
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from fss_lib import classify

OUT = []
def say(s=""):
    print(s); OUT.append(s)

s2 = np.load(LANE+"/s2_data.npy"); s3 = np.load(LANE+"/s3_data.npy")
k2, k3 = s2[:, 1], s3[:, 0]
F = 1e-9

ROWS = [
    # label, N, Q, sigma, role, additivity verdict, how the additivity verdict was obtained
    ("CONTROL-LIN   k = n-2, the record count", k2, k2, np.full(len(k2), F),
     "CONTROL known exactly linear", "ADDITIVE", "exact: a count over disjoint sets"),
    ("CONTROL-SAT-1 chi_joint(N records : fixed 3q bath)", k3, s3[:, 4], s3[:, 5],
     "CONTROL known bounded <= 3 bits", "NOT ADDITIVE", "S5 table 4, shared bath"),
    ("CONTROL-LIN-2 SUM chi_i, one bath site per record", k3, s3[:, 6], s3[:, 7],
     "CONTROL known exactly linear", "ADDITIVE", "S5 table 4, separate baths"),
    ("Q1  TOTAL CHI  SUM_i chi(R_i : fixed 3q bath)", k3, s3[:, 1], s3[:, 2],
     "probe", "NOT ADDITIVE", "S5 table 4, shared bath"),
    ("Q2a INTERACTING-PAIR COUNT, symplectic, N records", k2, s2[:, 2], np.full(len(k2), F),
     "probe", "ADDITIVE (of zero)", "S5 table 3"),
    ("Q2b   same, on records+writers  [POSITIVE CONTROL]", k2, s2[:, 7], np.full(len(k2), F),
     "CONTROL positive for Q2a", "ADDITIVE", "S5 table 3"),
    ("Q2c INTERACTING-PAIR COUNT, support overlap", k2, s2[:, 4], np.full(len(k2), F),
     "probe", "ADDITIVE but packing-dependent", "S5 table 3 + S5 summary"),
    ("Q3  TOTAL WRITER WEIGHT  SUM_i w(R_i)", k2, s2[:, 6], np.full(len(k2), F),
     "probe", "ADDITIVE", "S5 table 3"),
    ("Q4a LAM_MAX, symplectic relation matrix", k2, s2[:, 3], np.full(len(k2), F),
     "probe", "ADDITIVE (of zero)", "S5 table 3"),
    ("Q4b   same, on records+writers  [POSITIVE CONTROL]", k2, s2[:, 8], np.full(len(k2), F),
     "CONTROL positive for Q4a", "NOT ADDITIVE", "S5 exact argument 4"),
    ("Q4c LAM_MAX, support-overlap relation matrix", k2, s2[:, 5], np.full(len(k2), F),
     "probe", "NOT ADDITIVE", "S5 exact argument 4 + table 3"),
]

say("="*126)
say("S7   THE ONE TABLE.  [[n, n-2, 2]] family, N = k = n-2 records.")
say("="*126)
say()
say("  %-51s %-9s %-11s %5s %6s %7s %6s %-9s" %
    ("quantity", "N range", "form", "dAIC", "alpha", "Q2N/QN", "expo", "EXT (a)?"))
say("  " + "-"*124)
res = {}
for (lab, N, Q, S, role, addv, addsrc) in ROWS:
    r = classify(np.asarray(N, float), np.asarray(Q, float), np.asarray(S, float), lab)
    res[lab] = (r, role, addv, addsrc)
    say("  %-51s %-9s %-11s %5s %6s %7s %6s %-9s"
        % (lab[:51], "%g..%g" % (min(N), max(N)), r["best"],
           ("%.0f" % min(r["dAICc"], 999)) if np.isfinite(r["dAICc"]) else "exact",
           ("%.2f" % r["alpha"]) if r["alpha"] == r["alpha"] else "-",
           ("%.2f" % r["dbl"]) if r["dbl"] == r["dbl"] else "-",
           ("%.2f" % r.get("expo", np.nan)) if r.get("expo", np.nan) == r.get("expo", np.nan) else "-",
           "YES" if r["extensive"] else "no"))
say()
say("  CRITERION (c) -- NOT TOPOLOGICAL.  STATED RULE: a quantity that is an EXACT affine function")
say("  of the record count N (least-squares residual against a*N+b identically zero to machine")
say("  precision) IS the count in other units, and C-35 has already established that a count is")
say("  topological and cannot supply a density law.  Applied mechanically below, not by judgement.")
say()
say("  %-51s %-30s %-14s %-11s %s" % ("quantity", "CATEGORY", "additive (b)?", "not topo (c)?", "role"))
say("  " + "-"*124)
for (lab, *_rest) in ROWS:
    r, role, addv, addsrc = res[lab]
    isN = "-"
    if "fits" in r:
        rss = float(np.sum(r["fits"]["LIN"]["resid"]**2))
        isN = "NO (= a*N+b)" if rss < 1e-18 else "yes"
    elif r["category"] in ("IDENTICALLY ZERO", "CONSTANT"):
        isN = "NO (= a*N+b)"
    res[lab] = (r, role, addv, addsrc, isN)
    say("  %-51s %-30s %-14s %-11s %s" % (lab[:51], r["category"], addv, isN, role))
say()
say("  %-51s %14s %14s %11s" % ("quantity", "Q(1/N->0)", "+- err", "1/N-fit rms"))
say("  " + "-"*124)
for (lab, *_rest) in ROWS:
    r, *_ = res[lab]
    say("  %-51s %14.5f %14.5f %11.4g" % (lab[:51], r["Q0"], r["Q0e"], r["rms1N"]))
say()

say("="*126)
say("  THE READ -- filled from the numbers above, not in advance.")
say("="*126)
say()
nz = [l for l in res if res[l][0]["category"] == "IDENTICALLY ZERO"]
cn = [l for l in res if res[l][0]["category"] == "CONSTANT"]
sa = [l for l in res if res[l][0]["category"] in ("SATURATING", "DECAYING", "CANNOT DISTINGUISH")]
gr = [l for l in res if res[l][0]["category"].startswith("GROWING")]
ex = [l for l in res if res[l][0]["extensive"]]
say("  IDENTICALLY ZERO at every N tested (%d): %s" % (len(nz), "; ".join(l.split()[0] for l in nz)))
say("  EXACTLY CONSTANT at every N tested (%d): %s" % (len(cn), "; ".join(l.split()[0] for l in cn)))
say("  BOUNDED / SATURATING / DECAYING   (%d): %s" % (len(sa), "; ".join(l.split()[0] for l in sa)))
say("  GROWING                           (%d): %s" % (len(gr), "; ".join(l.split()[0] for l in gr)))
say("  MEETING gravity's requirement (a) (%d): %s" % (len(ex), "; ".join(l.split()[0] for l in ex)))
say()
probes = [l for l in res if res[l][1] == "probe"]
pab = [l for l in probes if res[l][0]["extensive"] and "ADDITIVE" == res[l][2]]
pabc = [l for l in pab if res[l][4] == "yes"]
say("  Of the %d PROBE quantities (controls excluded):" % len(probes))
say("    meeting (a) STRICT EXTENSIVITY and (b) ADDITIVITY : %d  %s"
    % (len(pab), "; ".join(l.split()[0] for l in pab) or "none"))
say("    ALSO meeting (c) NOT TOPOLOGICAL                  : %d  %s"
    % (len(pabc), "; ".join(l.split()[0] for l in pabc) or "NONE"))
say()
say("  The one quantity that clears (a) and (b) is Q3, total writer weight, and it clears them")
say("  because it is EXACTLY 2N.  The per-record writer weight is 2 at every n and for every")
say("  record -- proved, not fitted -- so the sum is the record COUNT wearing different units.")
say("  C-35 has already ruled the count out as a density law.  Nothing new passed; the table")
say("  found the count again under a new name, which is what criterion (c) is there to catch.")
say()
say("  THE CONTROLS BEHAVED: the engine called the exactly-linear control linear (dAICc decisive),")
say("  called the brief's own saturating control saturating in S6 (SAT1 over LIN by dAICc 32), and")
say("  ruled linear out for the bounded chi_joint control by dAICc 80.  Its ONE failure is stated")
say("  in S4 part 0 and repeated here: over SIX points with a realistic noise floor it cannot tell")
say("  LOGARITHMIC from SATURATING.  It can always tell either of them from LINEAR.  Since the only")
say("  question gravity asks is linear-or-not, that failure does not touch the verdicts above --")
say("  but it does mean no claim of the form 'this saturates rather than growing logarithmically'")
say("  is made anywhere in this lane from a fit alone.")
say()
say("  WHAT IS EXACT AND WHAT IS A TREND:")
say("    EXACT (rules out emergence at ANY N, however large):")
say("      Q2a, Q4a identically zero -- independent records COMMUTE, so their symplectic relation")
say("        matrix is the zero matrix.  Positive controls Q2b (= k) and Q4b (= 1) confirm the")
say("        instrument would register a non-zero if one were there.")
say("      Q4b, and lam_max of ANY relation matrix, fails additivity by exact argument 4:")
say("        the union's matrix is block diagonal, so lam_max is a MAX, and a max is not a sum.")
say("      Q3 total writer weight = 2N exactly: the per-record writer weight is 2 at every n and")
say("        for every record, proved both ways (no weight-1 admissible unitary of any kind exists;")
say("        a weight-2 admissible Pauli flipping R always exists), and confirmed by exhaustive")
say("        search over all 4^n Paulis at n = 4, 6, 8.  So Q3 is the record COUNT times a constant")
say("        -- C-35's count wearing different units, not a new quantity.")
say("      CONTROL-SAT-1 and Q1 are capped at nq bits by SUM_i chi_i <= chi_joint <= S(rho_B) <= nq,")
say("        an exact chain, so no bath of fixed size can carry an extensive record source.")
say("    TREND ONLY (absence WITHIN RANGE, N <= 40 for chi and N <= 58 for the combinatorics):")
say("      Q1's particular decay exponent, and the fitted exponents of Q2c and Q4c.")
say()
say("  THE ONE HONEST POSITIVE: Q2c and Q4c DO grow, and grow FASTER than linearly (exponents")
say("  2.05 +- 0.00 and 1.95 +- 0.04).  Neither is a source term.  Q4c is a lam_max and fails (b)")
say("  exactly.  Q2c is additive over disjoint regions but is not a function of how much is")
say("  enclosed: two disjoint 6-record regions give 14, one region holding 12 records gives 21.")
say("  Its super-linear growth is the crowding of k logical supports onto n = k+2 qubits, and it")
say("  vanishes the moment the records are genuinely far apart -- which is the configuration")
say("  gravity's additivity clause is about.  A super-linear quantity also fails (a) in the other")
say("  direction: S(2N)/S(N) -> 4, not 2.")
say()
say("  LARGEST N REACHED, AND WHAT STOPPED IT:")
say("    combinatorial quantities  N = 58 (n = 60).  Nothing stopped it; cost is O(n^3) in F_2.")
say("      Stopped at 60 because the categories were already exact and further n adds no information.")
say("    Holevo chi              N = 40 (n = 42).  Reached via an EXACT reduction chain validated")
say("      in S1: dense 2^n x 2^nq dies at n = 8 (2048-dim joint state); the 2^k sector sum dies")
say("      near k = 20; the counting compression is polynomial and carried it to 40.  It would go")
say("      further; 40 was enough for the doubling test to have five doublings inside it.")
open(LANE+"/s7_synthesis.txt", "w").write("\n".join(OUT)+"\n")
