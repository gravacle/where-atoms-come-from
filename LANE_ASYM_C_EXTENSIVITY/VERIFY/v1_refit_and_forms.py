"""V1 -- ADVERSARIAL RE-FIT.  Do the lane's 5-point windows actually separate the asymptotic
forms they are said to separate?  And do the fitted exponents agree with EXACT asymptotics
derived independently here?

Competing forms fitted to each series over the lane's own window and over the full range:
   POW      y = A N^p
   POWLOG   y = A N^p log N          (a power law with a log correction)
   SAT      y = S_inf - B N^-q       (saturating)
   LOG      y = A + B log N
Reported: max |log-space residual| for each.  If two forms fit within a factor of a few of
each other on the window, the window does NOT separate them and any claim that it does is
refuted.  Exact asymptotics are derived below and compared to the fits.
"""
import sys, math, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY")
from qcore import (chi_of_n, chi_register_site, J2_same_site, phi_moments,
                   loglog_fit, LAM, BETA, ENERGIES)
# no scipy in this environment: every competing form below is fitted by an EXACT linear
# least-squares in its linear parameters, with the one non-linear parameter swept on a fine grid.

OUT = []
def P(s=""):
    print(s); OUT.append(s)

P("=" * 110)
P("V1  ADVERSARIAL RE-FIT OF THE LANE'S GROWTH LAWS")
P("=" * 110)

NN  = [256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]
BIG = [4096, 8192, 16384, 32768, 65536]           # the lane's own window: 5 points, factor 16

series = {
 "chi_total SHARED 1 site": [N * chi_of_n(N) for N in NN],
 "chi_reg SHARED nq=3":     [sum(chi_register_site(N, e)[0] for e in ENERGIES[:3]) for N in NN],
 "sum_pairs|J_2| 1 site":   [(N*(N-1)/2)*abs(J2_same_site(N)) for N in NN],
 "|J_2(N)| 1 site":         [abs(J2_same_site(N)) for N in NN],
 "Var(Phi) 1 shared site":  [phi_moments(N)[1]**2 for N in NN],
 "spread(Phi) 1 shared site":[phi_moments(N)[2] for N in NN],
 "std(Phi) 1 shared site":  [phi_moments(N)[1] for N in NN],
}

def resid_pow(x, y):
    lx, ly = np.log(x), np.log(y)
    A = np.vstack([lx, np.ones_like(lx)]).T
    c, *_ = np.linalg.lstsq(A, ly, rcond=None)
    return float(np.abs(ly - A @ c).max()), c[0]

def _linfit_logresid(cols, x, y):
    """least squares y ~ cols (design matrix), residual measured in LOG space like the others"""
    A = np.asarray(cols, float).T
    yv = np.asarray(y, float)
    c, *_ = np.linalg.lstsq(A, yv, rcond=None)
    pred = A @ c
    if np.any(pred <= 0): return np.inf, c
    return float(np.abs(np.log(yv) - np.log(pred)).max()), c

def res_powlog(x, y):
    """y = A N^p log N   ->  log y - log log N = log A + p log N   (exact linear)"""
    lx = np.log(np.asarray(x, float)); ly = np.log(np.asarray(y, float)) - np.log(lx)
    A = np.vstack([lx, np.ones_like(lx)]).T
    c, *_ = np.linalg.lstsq(A, ly, rcond=None)
    return float(np.abs(ly - A @ c).max()), c

def res_sat(x, y):
    """y = S - B N^-q ; q swept on a grid, (S,B) exact by linear LS"""
    x = np.asarray(x, float); best = (np.inf, None)
    for q in np.concatenate([np.linspace(0.005, 3.0, 1200), np.linspace(3.0, 12.0, 200)]):
        r, c = _linfit_logresid([np.ones_like(x), -x**(-q)], x, y)
        if r < best[0]: best = (r, (c[0], c[1], q))
    return best

def res_log(x, y):
    """y = A + B log N (exact linear)"""
    x = np.asarray(x, float)
    return _linfit_logresid([np.ones_like(x), np.log(x)], x, y)

P("\n--- (1a) COMPETING FORMS on the LANE'S OWN 5-POINT WINDOW N in [4096, 65536] ---")
P("%-27s %-9s %-11s %-11s %-11s %-11s %s" % ("series", "p(POW)", "res POW", "res POWLOG",
                                             "res SAT", "res LOG", "separated?"))
P("-" * 110)
for nm, ys in series.items():
    yb = [ys[NN.index(N)] for N in BIG]
    rp, p = resid_pow(BIG, yb)
    rpl, _ = res_powlog(BIG, yb)
    rs,  _ = res_sat(BIG, yb)
    rl,  _ = res_log(BIG, yb)
    others = [r for r in (rpl, rs, rl) if not np.isnan(r)]
    best_other = min(others) if others else np.nan
    sep = "YES (POW wins x%.0f)" % (best_other/max(rp,1e-18)) if best_other > 10*rp else \
          "NO -- forms indistinguishable"
    P("%-27s %-9.4f %-11.2e %-11.2e %-11.2e %-11.2e %s" % (nm, p, rp, rpl, rs, rl, sep))

P("\n--- (1b) SAME FORMS on the FULL range N in [256, 65536] (9 points, factor 256) ---")
P("%-27s %-9s %-11s %-11s %-11s %-11s %s" % ("series", "p(POW)", "res POW", "res POWLOG",
                                             "res SAT", "res LOG", "verdict"))
P("-" * 110)
for nm, ys in series.items():
    rp, p = resid_pow(NN, ys)
    rpl, _ = res_powlog(NN, ys)
    rs,  _ = res_sat(NN, ys)
    rl,  _ = res_log(NN, ys)
    others = [r for r in (rpl, rs, rl) if not np.isnan(r)]
    best_other = min(others) if others else np.nan
    verdict = "POW clearly best" if best_other > 10*rp else "NOT separated"
    P("%-27s %-9.4f %-11.2e %-11.2e %-11.2e %-11.2e %s" % (nm, p, rp, rpl, rs, rl, verdict))

# ------------------------------------------------------------------ exact asymptotics
P("\n" + "-" * 110)
P("(1c) INDEPENDENT EXACT ASYMPTOTICS.  Derived here, not taken from the lane.")
P("-" * 110)
P("  f(c) = -(1/b) ln 2cosh(b E),  E = sqrt(e^2 + lam^2 c^2)  ->  f(c) ~ -lam|c| for large |c|.")
P("  c = sum of N iid +-1  ->  c/sqrt(N) -> Normal(0,1),  density at 0  p(0) = 1/sqrt(2 pi N).")
P("  J_2(N) = 1/2 E[f(c''+2) - f(c'')] with c'' a sum of N-2 signs")
P("         ~ -(lam/2) E[|c+2|-|c|] = -(lam/2)(4 p(0)) = -2 lam / sqrt(2 pi N)   -> EXACTLY N^-1/2")
P("  sum_pairs |J_2| ~ (N^2/2)(2 lam/sqrt(2 pi N)) = lam N^{3/2}/sqrt(2 pi)      -> EXACTLY N^+3/2")
P("  std(Phi)  ~ lam * std(|c|) = lam sqrt(N) sqrt(1 - 2/pi)                     -> EXACTLY N^+1/2")
P("  spread(Phi) ~ lam N + const                                                 -> EXACTLY N^+1")
P("")
P("%-24s %-22s %-22s %-11s %s" % ("quantity", "value at N=65536", "exact asymptotic", "rel err", "exact p"))
P("-" * 110)
Nb = 65536
pred = {
 "|J_2(N)|":        (abs(J2_same_site(Nb)), 2*LAM/math.sqrt(2*math.pi*Nb), -0.5),
 "sum_pairs|J_2|":  ((Nb*(Nb-1)/2)*abs(J2_same_site(Nb)), LAM*Nb**1.5/math.sqrt(2*math.pi), 1.5),
 "std(Phi)":        (phi_moments(Nb)[1], LAM*math.sqrt(Nb)*math.sqrt(1-2/math.pi), 0.5),
 "spread(Phi)":     (phi_moments(Nb)[2], LAM*Nb, 1.0),
}
for nm, (v, a, pex) in pred.items():
    P("%-24s %-22.10g %-22.10g %-11.2e %+.1f" % (nm, v, a, abs(v-a)/abs(a), pex))

P("\n  => the lane's fitted exponents (-0.4998, 1.5003, 0.5019, 1.0001) reproduce these EXACT")
P("     values.  The 5-point window is NOT the evidence: the closed forms are.  Their")
P("     sigma(p) column is however NOT a statistical uncertainty -- the series are exact")
P("     deterministic evaluations, so sigma(p) measures CURVATURE, not noise.  Reported as")
P("     'uncertainty on the exponent' it is a category error, though a harmless one here.")

# ------------------------------------------------------------------ chi_total shared: which power?
P("\n" + "-" * 110)
P("(1d) chi_total SHARED: the lane's own two statements disagree.")
P("-" * 110)
ys = series["chi_total SHARED 1 site"]
for a, b in zip(NN[:-1], NN[1:]):
    ya, yb = ys[NN.index(a)], ys[NN.index(b)]
    P("   N %6d -> %6d   local slope %+.4f   ratio %.6f" % (a, b, math.log(yb/ya)/math.log(2.0), yb/ya))
P("   s5 TABLE 17 prints exponent -1.011201 for this series.")
P("   s5 TABLE 18 prints 'grows as N^-0.50' for the SAME series, and READ OF S5 item 4 repeats")
P("   'chi_total on a shared bath DECAYS as N^-0.50'.  The local slopes above are ~ -1.0.")
P("   => a CONCLUSION CONTRADICTS THE TABLE ABOVE IT.  (The sign of the verdict is unaffected:")
P("      the series decays either way, but the printed exponent in TABLE 18/READ is wrong.)")

# ------------------------------------------------------------------ chi_reg bound is ordinary
P("\n" + "-" * 110)
P("(1e) IS THE 'SATURATION' AN ORDINARY BOUND?  chi(register:bath) <= S(bath) <= nq bits.")
P("-" * 110)
def Stau(e):
    z = math.tanh(BETA*e); p = (1+z)/2
    return -(p*math.log2(p) + (1-p)*math.log2(1-p))
tot = sum(1-Stau(e) for e in ENERGIES[:3])
P("   lane's 'exact bound' nq(1-S(tau))  = %.10f" % tot)
P("   textbook Holevo ceiling log2(dim bath) = %.10f bits (3 qubits)" % 3.0)
P("   per-site  1 - S(tau_j) for e = 1.0, 1.4, 0.7 : %s" %
      ", ".join("%.6f" % (1-Stau(e)) for e in ENERGIES[:3]))
P("   => the ceiling is S(rho_bar) - S(tau) summed over bath qubits: this is the ORDINARY")
P("      Holevo/entropy bound on a fixed-size quantum channel, tightened only by the bath's")
P("      thermal purity.  It is a statement about the BATH's capacity, not about records.")
P("      The lane says as much; but 'SATURATING' is therefore not a discovery about records.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY/VERIFY/v1_refit_and_forms.txt",
     "w").write("\n".join(OUT) + "\n")
