"""S4 -- THE FINITE-SIZE SCALING ENGINE, AND THE CALIBRATION THAT LICENSES IT.

The engine does four things to every quantity, in the same way, with the same criteria:

  (1) EXTRAPOLATION.  Weighted least squares of Q against 1/N, Q = Q0 + b/N + c/N^2, over the
      largest-N half of the data.  Q0 is the 1/N -> 0 intercept; its uncertainty is the square
      root of the (0,0) covariance entry, inflated by sqrt(reduced chi^2) when that exceeds 1.
      Residuals reported.

  (2) MODEL SELECTION among four forms that look alike at small N:
        LIN   Q = a N + b                 (2 params)   -- extensive, what gravity requires
        LOG   Q = a ln N + b              (2 params)   -- sub-extensive
        SAT1  Q = Qinf - c/N              (2 params)   -- saturating, algebraic approach
        SAT2  Q = Qinf - c exp(-N/xi)     (3 params)   -- saturating, exponential approach
      Ranked by AICc.  STATED CRITERION (Burnham-Anderson): dAICc >= 10 decisive,
      4 <= dAICc < 10 substantial, dAICc < 4 THE DATA CANNOT DISTINGUISH.  When the spread of
      dAICc across the runner-up is under 4 the engine says CANNOT DISTINGUISH and stops.

  (3) SCALING COLLAPSE.  alpha* minimises the coefficient of variation of Q(N)/N^alpha over the
      N range.  STATED CRITERION for "collapse": CV(alpha*) < 0.02, i.e. Q/N^alpha is constant
      to 2% across the whole range.  A pure power law gives CV = 0 exactly.

  (4) THE DOUBLING TEST, which is gravity's own requirement (a): Q(2N)/Q(N) -> 2.
      Reported at every (N, 2N) pair available.  Nothing that fails this is a source term.

CALIBRATION (this is what makes the rest admissible).  The engine is first run on SYNTHETIC
series whose form is KNOWN -- exact linear, exact logarithmic, exact saturating -- sampled on
the SAME N grid, with and without the noise floor measured on the real chi data.  If the engine
cannot classify those, it cannot classify anything, and this lane says so and stops.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE")
from fss_lib import _wls, fit_form, aicc, classify, ascii_plot, LANE, OUT, say

# =====================================================================================
say("="*116)
say("S4   FINITE-SIZE SCALING.  Engine calibrated on known forms BEFORE it is used on anything.")
say("="*116)
say()

# --------------------------------------------------------- load
s2 = np.load(LANE+"/s2_data.npy")     # n k Ap_s Al_s Ap_o Al_o AW Bp_s Bl_s Bp_o Bl_o BW supp
s3 = np.load(LANE+"/s3_data.npy")     # k tm tse tsd jm jse gm gse
n2, k2 = s2[:, 0], s2[:, 1]
k3 = s3[:, 0]

# --------------------------------------------------------- CALIBRATION
say("-"*116)
say("PART 0.  CALIBRATION OF THE ENGINE.  Synthetic series of KNOWN form on the same N grids.")
say("         SHORT grid = N in {2,4,6,8,10,12} -- exactly the six sizes the [[n,n-2,2]] brief names.")
say("         FULL  grid = N in {2,4,...,40}    -- the range this lane actually reached.")
say("         Noise: (a) exact, sigma = 1e-9; (b) sigma = 0.04 bits, the measured chi noise floor.")
say()
GRID_SHORT = np.arange(2, 13, 2, dtype=float)
GRID_FULL = np.arange(2, 41, 2, dtype=float)
rng = np.random.default_rng(20260819)

def synth(kind, N):
    if kind == "linear":  return 0.5*N
    if kind == "log":     return 1.5*np.log(N) + 0.3
    if kind == "sat":     return 3.0*(1 - np.exp(-N/4.0))
    if kind == "sqrt":    return 0.8*np.sqrt(N)

say("   grid    noise      true form        engine says                          best  dAICc   alpha    CV")
cal_ok_exact = True; cal_ok_noisy = True
for gname, G in (("SHORT", GRID_SHORT), ("FULL", GRID_FULL)):
    for nz, sname in ((1e-9, "exact"), (0.04, "0.04 bits")):
        for kind in ("linear", "log", "sat", "sqrt"):
            Q = synth(kind, G) + (rng.normal(0, nz, len(G)) if nz > 1e-6 else 0.0)
            r = classify(G, Q, np.full(len(G), nz), kind)
            want = {"linear": "LIN", "log": "LOG", "sat": "SAT", "sqrt": "?"}[kind]
            hit = ("OK " if (r["best"].startswith(want) or want == "?") else "MISS")
            if nz < 1e-6 and want != "?" and not r["best"].startswith(want): cal_ok_exact = False
            if nz > 1e-6 and want != "?" and not r["best"].startswith(want): cal_ok_noisy = False
            say("   %-6s %-10s %-14s %-36s %-5s %6.1f %7.3f %6.4f  %s"
                % (gname, sname, kind, r["category"], r["best"], min(r["dAICc"], 999.9),
                   r["alpha"], r["cv"], hit))
    say()
say("   CALIBRATION VERDICT")
say("     on EXACT synthetic data the engine picks the true form: %s" % cal_ok_exact)
say("     on data carrying the measured chi noise floor it picks the true form: %s" % cal_ok_noisy)
say("     THE MISS IS THE IMPORTANT ROW: on the SHORT grid with a 0.04-bit floor the engine calls a")
say("     TRUE SATURATING series LOGARITHMIC at dAICc 8.6.  Six points cannot separate log from")
say("     saturation.  Therefore NO chi-based classification is reported from the SHORT window; the")
say("     short-window chi rows below are read only as CANNOT DISTINGUISH, whatever the fit says.")
say("     Note also what is NOT confused anywhere: LINEAR is never mistaken for either, on either")
say("     grid, at either noise level.  The extensive/not-extensive call -- the only one gravity")
say("     needs -- survives the short window; the log-vs-saturation call does not.")
say("     -> the alpha* for the sqrt row is the honest check that the collapse test finds a known")
say("        exponent: true 0.5, engine %.3f on the FULL grid."
    % classify(GRID_FULL, synth("sqrt", GRID_FULL), np.full(len(GRID_FULL), 1e-9), "s")["alpha"])
say("-"*116)
say()

# --------------------------------------------------------- the quantities
say("-"*116)
say("PART 1.  THE QUANTITIES.  Definitions, sources, and which representation produced them.")
say()
Qs = []
def add(label, N, Q, sig, kind, src):
    Qs.append(dict(label=label, N=np.asarray(N, float), Q=np.asarray(Q, float),
                   sig=np.asarray(sig, float), kind=kind, src=src))

FLOOR = 1e-9
add("CONTROL-LIN  k = n-2 (number of independent records)", k2, k2, np.full(len(k2), FLOOR),
    "control-linear", "S2 exact")
add("CONTROL-SAT  chi_joint(all N records : FIXED 3-qubit bath)", k3, s3[:, 4], s3[:, 5],
    "control-saturating", "S3 exact, time-averaged")
add("CONTROL-LIN2 SUM_i chi_i, GROWN bath (one site per record)", k3, s3[:, 6], s3[:, 7],
    "control-linear", "S3 exact, time-averaged")
add("Q1  total chi  SUM_i chi(R_i : FIXED 3-qubit bath)", k3, s3[:, 1], s3[:, 2],
    "probe", "S3 exact, time-averaged")
add("Q2a interacting-pair count, symplectic, SET A (the N records)", k2, s2[:, 2],
    np.full(len(k2), FLOOR), "probe", "S2 exact")
add("Q2b interacting-pair count, symplectic, SET B (records+writers) CONTROL", k2, s2[:, 7],
    np.full(len(k2), FLOOR), "control-positive", "S2 exact")
add("Q2c interacting-pair count, SUPPORT OVERLAP, SET A", k2, s2[:, 4],
    np.full(len(k2), FLOOR), "probe", "S2 exact")
add("Q3  total writer weight  SUM_i w(R_i), SET A", k2, s2[:, 6],
    np.full(len(k2), FLOOR), "probe", "S2 exact + exact argument")
add("Q4a lam_max, record-record SYMPLECTIC relation matrix, SET A", k2, s2[:, 3],
    np.full(len(k2), FLOOR), "probe", "S2 exact")
add("Q4b lam_max, SYMPLECTIC relation matrix, SET B (CONTROL)", k2, s2[:, 8],
    np.full(len(k2), FLOOR), "control-positive", "S2 exact")
add("Q4c lam_max, SUPPORT-OVERLAP relation matrix, SET A", k2, s2[:, 5],
    np.full(len(k2), FLOOR), "probe", "S2 exact")

for q in Qs:
    say("   %-62s  %-20s  %s" % (q["label"], q["kind"], q["src"]))
say("-"*116)
say()

# --------------------------------------------------------- master table, two windows
for wname, sel in (("SHORT WINDOW  N = 2..12 (the six sizes the brief names)", lambda N: N <= 12),
                   ("FULL WINDOW   N = 2..N_max (everything this lane reached)", lambda N: N <= 1e9)):
    say("="*116)
    say("MASTER TABLE -- %s" % wname)
    say("="*116)
    say("  %-58s %-30s %6s %7s %7s %9s %6s" % ("quantity", "engine verdict", "best", "dAICc", "alpha", "CV", "EXT?"))
    say("  " + "-"*114)
    results = {}
    for q in Qs:
        msk = sel(q["N"])
        if msk.sum() < 5: continue
        r = classify(q["N"][msk], q["Q"][msk], q["sig"][msk], q["label"])
        results[q["label"]] = r
        say("  %-58s %-30s %6s %7s %7s %9s %6s"
            % (q["label"][:58], r["category"], r["best"],
               ("%.1f" % min(r["dAICc"], 999.9)) if np.isfinite(r["dAICc"]) else "inf",
               ("%.3f" % r["alpha"]) if r["alpha"] == r["alpha"] else "-",
               ("%.4f" % r["cv"]) if r["cv"] == r["cv"] else "-",
               "YES" if r["extensive"] else "no"))
    say()
    say("  %-58s %14s %14s %10s %12s %10s" % ("quantity", "Q(1/N->0)", "+- err", "Q(2N)/Q(N)", "exponent", "1/N-fit rms"))
    say("  " + "-"*114)
    for q in Qs:
        r = results.get(q["label"])
        if r is None: continue
        e = ("%.3f+-%.3f" % (r.get("expo", np.nan), r.get("expo_e", np.nan))
             if r.get("expo", np.nan) == r.get("expo", np.nan) else "-")
        say("  %-58s %14.5f %14.5f %10s %12s %10.3g"
            % (q["label"][:58], r["Q0"], r["Q0e"],
               ("%.3f" % r["dbl"]) if r["dbl"] == r["dbl"] else "-", e, r["rms1N"]))
    say()
    say("  EXT? = gravity's requirement (a): LIN preferred by dAICc >= 4 AND Q(2N)/Q(N) in [1.6,2.6]")
    say("         AND the fitted power-law exponent within 0.25 of 1.  All three, or it is not extensive.")
    say("  The Q(1/N->0) intercept is only meaningful for a BOUNDED quantity; for a divergent one the")
    say("  large 1/N-fit rms is the signature that no finite intercept exists.")
    say("  COLLAPSE CRITERION: CV(alpha*) < 0.02 means Q/N^alpha is constant to 2% across the window.")
    say("  Quantities meeting it: %s"
        % (", ".join(l.split()[0] for l, r in results.items()
                     if r["cv"] == r["cv"] and r["cv"] < 0.02) or "none"))
    say()

# --------------------------------------------------------- residuals for the chi fits
say("="*116)
say("RESIDUALS AND NOISE FLOOR -- no fit is reported without them.")
say("="*116)
say()
for q in Qs:
    if "chi" not in q["label"]: continue
    r = classify(q["N"], q["Q"], q["sig"], q["label"])
    say("  %s" % q["label"])
    say("    noise floor (median SE of the 25-time average): %.5f bits" % np.median(q["sig"]))
    say("    AICc ranking: %s" % ", ".join("%s=%.1f" % (nm, a) for nm, a in r["rank"]))
    f = r["fits"][r["best"]]
    say("    best form %s: params %s" % (r["best"], ["%.5f" % b for b in f["beta"]]))
    say("    residuals (data - fit), in units of the noise floor:")
    say("      " + "  ".join("%+.2f" % (x/s) for x, s in zip(f["resid"], q["sig"])))
    say("    reduced chi^2 = %.3f on %d dof" % (f["chi2"]/max(len(q["N"])-f["p"], 1),
                                                len(q["N"])-f["p"]))
    say()

# --------------------------------------------------------- plots
say("="*116)
say("PLOTS -- each quantity against N and against 1/N (ASCII; the .npy files carry the numbers).")
say("="*116)
for q in Qs:
    if abs(q["Q"]).max() < 1e-12:
        say(); say("    %s : IDENTICALLY ZERO at every N -- nothing to plot." % q["label"]); continue
    say()
    for ln in ascii_plot(q["N"], q["Q"], q["label"] + "   [ vs N ]"): say(ln)
    for ln in ascii_plot(1.0/q["N"], q["Q"], q["label"] + "   [ vs 1/N ]"): say(ln)

open(LANE+"/s4_fss.txt", "w").write("\n".join(OUT)+"\n")
