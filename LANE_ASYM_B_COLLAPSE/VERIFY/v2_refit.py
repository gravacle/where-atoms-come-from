"""V2 -- ADVERSARIAL RE-FIT OF THE LANE'S OWN chi DATA.

The lane's sigma is the STANDARD ERROR OF A 25-POINT TIME AVERAGE.  chi(t) is a deterministic
oscillating function, not a noisy measurement, so that "sigma" is the amplitude of an
oscillation, not a noise floor.  AICc is linear in chi^2 = sum (r/sigma)^2, so the whole model
selection is proportional to 1/sigma^2 and can be moved at will by that choice.

THIS SCRIPT: refit CONTROL-SAT (chi_joint) and Q1 (total chi) and CONTROL-LIN2 under
  LIN   a N + b
  LOG   a ln N + b
  SAT1  Qinf - c/N
  SAT2  Qinf - c exp(-N/xi)
  POW   a N^b          <-- a form the lane's engine does NOT carry, and the one a decaying
                           series most plausibly follows
with FOUR sigma choices:
  (s1) the lane's own SE of the time average
  (s2) the STD over the 25 times (25x larger)
  (s3) a flat sigma = median SE
  (s4) a flat sigma = 0.005 bits (near-exact: chi itself is computed exactly)
and report whether the lane's classification survives.

Also: the SHORT-WINDOW test.  Refit on N = 2..12 only (the six sizes the brief names) and see
what the engine would have said there.
"""
import sys, numpy as np
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE"
sys.path.insert(0, LANE)
from fss_lib import fit_form, aicc, _wls

OUT = []
def say(s=""):
    print(s); OUT.append(s)

s3 = np.load(LANE+"/s3_data.npy")     # k tm tse tsd jm jse gm gse
N = s3[:, 0]
SER = {
    "Q1  SUM_i chi_i (FIXED 3q bath)":  (s3[:, 1], s3[:, 2], s3[:, 3]),
    "CONTROL-SAT chi_joint (FIXED)":    (s3[:, 4], s3[:, 5], s3[:, 5]*np.sqrt(25)),
    "CONTROL-LIN2 SUM chi_i (GROWN)":   (s3[:, 6], s3[:, 7], s3[:, 7]*np.sqrt(25)),
}

def fit_pow(N, Q, sig):
    """a N^b by weighted LS in log space, then chi^2 in LINEAR space so AICc is comparable."""
    X = np.c_[np.log(N), np.ones_like(N)]
    beta, _, _, _ = _wls(X, np.log(Q), np.maximum(sig/Q, 1e-12))
    pred = np.exp(beta[1])*N**beta[0]
    chi2 = float(np.sum(((Q-pred)/sig)**2))
    return dict(name="POW", p=2, beta=list(beta), chi2=chi2, resid=Q-pred)

def rank(N, Q, sig):
    fits = {nm: fit_form(nm, N, Q, sig) for nm in ("LIN", "LOG", "SAT1", "SAT2")}
    fits["POW"] = fit_pow(N, Q, sig)
    m = len(N)
    sc = sorted(((aicc(f["chi2"], m, f["p"]), nm) for nm, f in fits.items()))
    return sc, fits

say("="*118)
say("V2  RE-FIT OF THE LANE'S chi DATA UNDER FOUR SIGMA CHOICES AND ONE EXTRA MODEL (POW = a N^b)")
say("    AICc is linear in 1/sigma^2, so the sigma choice IS the model-selection knob.  Here it is turned.")
say("="*118)
say()

for label, (Q, se, sd) in SER.items():
    say("-"*118)
    say(label)
    say()
    sigmas = [("lane's SE of time-avg", se),
              ("STD over the 25 times", sd),
              ("flat = median SE", np.full(len(N), float(np.median(se)))),
              ("flat = 0.005 bits", np.full(len(N), 0.005))]
    say("    sigma choice              best   dAICc(best-2nd)   full AICc ranking")
    for sname, sg in sigmas:
        sc, fits = rank(N, Q, sg)
        d = sc[1][0]-sc[0][0]
        say("    %-24s  %-5s  %14.2f   %s" % (sname, sc[0][1], d,
            "  ".join("%s=%.1f" % (nm, a) for a, nm in sc)))
    say()
    # SHORT window
    msk = N <= 12
    say("    SHORT WINDOW N=2..12 (the six sizes the brief names), lane's own sigma:")
    sc, fits = rank(N[msk], Q[msk], se[msk])
    say("      best %s   dAICc %.2f   ranking %s" % (sc[0][1], sc[1][0]-sc[0][0],
        "  ".join("%s=%.1f" % (nm, a) for a, nm in sc)))
    say()

say("-"*118)
say("RESIDUALS OF THE LANE'S PREFERRED FORM FOR Q1 (SAT2) AGAINST THE POWER LAW, lane's sigma")
Q, se, sd = SER["Q1  SUM_i chi_i (FIXED 3q bath)"]
sc, fits = rank(N, Q, se)
for nm in ("SAT2", "POW", "LOG", "SAT1", "LIN"):
    r = fits[nm]["resid"]
    say("   %-5s chi2=%9.2f  p=%d  AICc=%9.2f  max|r|/sigma=%6.2f  rms r = %.5f bits"
        % (nm, fits[nm]["chi2"], fits[nm]["p"], aicc(fits[nm]["chi2"], len(N), fits[nm]["p"]),
           float(np.max(np.abs(r/se))), float(np.sqrt(np.mean(r**2)))))
say()
say("   The N=4 point carries residual %+.2f sigma under SAT2 -- one point drives the chi^2."
    % float((Q-(Q-fits['SAT2']['resid']))[1]/se[1]))
say()

# -------- how far does an N<=40 window let you extrapolate?
say("-"*118)
say("EXTRAPOLATION HONESTY CHECK.  Fit Q1 on N<=20 only, predict N=40, compare to truth.")
m2 = N <= 20
for nm in ("SAT2", "POW", "LOG", "SAT1", "LIN"):
    if nm == "POW":
        f = fit_pow(N[m2], Q[m2], se[m2]); b = f["beta"]
        pred = np.exp(b[1])*40.0**b[0]
    else:
        f = fit_form(nm, N[m2], Q[m2], se[m2]); b = f["beta"]
        if nm == "LIN":    pred = b[0]*40+b[1]
        elif nm == "LOG":  pred = b[0]*np.log(40)+b[1]
        elif nm == "SAT1": pred = b[0]-b[1]/40
        else:              pred = b[0]-b[1]*np.exp(-40/b[2])
    say("   %-5s predicts chi(N=40) = %9.5f   truth = %9.5f   miss = %+8.5f bits (%.1f sigma)"
        % (nm, pred, Q[-1], pred-Q[-1], (pred-Q[-1])/se[-1]))
say("-"*118)
open(LANE+"/VERIFY/v2_refit.txt", "w").write("\n".join(OUT)+"\n")
