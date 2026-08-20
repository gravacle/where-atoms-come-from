"""V2 -- ADVERSARIAL CHECK OF THE CAPACITY EXPONENT  chi ~ q^(-1.5438 +- 0.0206),
"26 sigma steeper than the p = -1 that equipartition of a fixed site capacity requires"
(lane B9 Part 3).

TWO OBJECTIONS, BOTH RUNNABLE.

(1) IS IT A POWER LAW AT ALL?  q = 2..8 is 7 points over less than one decade.  A local slope
    is not an exponent.  The lane's "+- 0.0206" is a REGRESSION standard error computed from
    the scatter of a DETERMINISTIC curve -- there is no statistical noise in this computation
    at all, so "26 sigma" measures how smoothly a non-power-law bends, not how well an
    exponent is determined.  Extend q far past 8 and see whether the local slope DRIFTS.

(2) IS THE ORDINARY EXPLANATION SUFFICIENT?  For q mutually commuting records on one bath
    site the coupling is lam*(R + S) (x) X_0 with S a sum of q-1 independent commuting
    +-1 Paulis, so the bath sees the read bit r shifted by a BINOMIAL random variable of
    variance q-1.  That is the textbook "signal 1 against noise sqrt(q)" channel, whose
    mutual information falls as 1/q ASYMPTOTICALLY but is nowhere near its asymptote at
    q = 8.  Compute the exact large-q behaviour and compare.

METHOD.  The spectral reduction validated in V1 (agreement with the lane's own pipeline to
1e-14) makes any q reachable: chi depends only on the binomial distribution of S.
POSITIVE CONTROL (D-15): the reduction is re-checked here against the lane's published
q = 2..8 numbers before it is used to extrapolate.
"""
import numpy as np, sys, os, time
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from common import *          # noqa

t0 = time.time()
TIMES = np.linspace(1.0, 13.0, 25)

def _vn(rho):
    e = np.linalg.eigvalsh(rho); e = e[e > 1e-13]
    return float(-(e * np.log2(e)).sum())

def chi_crowd(q, lam, env, times=TIMES):
    """q mutually commuting records on ONE bath site; read one of them.
       S = sum of the other q-1, binomially distributed."""
    m = q - 1
    dist = [(m - 2 * k, comb(m, k) / 2 ** m) for k in range(m + 1)] if m > 0 else [(0, 1.0)]
    rho0 = env.thermal(); X0 = env.site[0]
    fields = sorted({round(lam * (r + s), 12) for r in (1, -1) for s, _ in dist})
    prop = {f: np.linalg.eigh(env.HB + f * X0) for f in fields}
    acc = []
    for t in times:
        U = {}
        for f, (w, V) in prop.items():
            M = (V * np.exp(-1j * w * t)[None, :]) @ V.conj().T
            U[f] = M @ rho0 @ M.conj().T
        rr = {r: sum(p * U[round(lam * (r + s), 12)] for s, p in dist) for r in (1, -1)}
        av = 0.5 * (rr[1] + rr[-1])
        acc.append(max(_vn(av) - 0.5 * (_vn(rr[1]) + _vn(rr[-1])), 0.0))
    return float(np.mean(acc))

def fitp(qs, cs):
    x = np.log(np.array(qs, float)); y = np.log(np.array(cs))
    A = np.vstack([np.ones_like(x), x]).T
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ c
    dof = max(len(x) - 2, 1)
    s2 = float(res @ res) / dof
    cov = s2 * np.linalg.inv(A.T @ A)
    return float(c[1]), float(np.sqrt(cov[1, 1])), float(np.sqrt(np.mean(res ** 2)))

say("=" * 118)
say("V2   IS chi ~ q^-1.54 AN EXPONENT, OR THE LOCAL SLOPE OF A CROSSOVER?")
say("=" * 118)

env = BASE.env(3)
say("")
say("  POSITIVE CONTROL: reproduce the lane's own B9 Part-3 column (n = 10, lam = 0.8, V0).")
LANE = {1: 0.521527300760, 2: 0.136408688972, 3: 0.071783219844, 4: 0.045162574990,
        5: 0.033530458642, 6: 0.024307043226, 7: 0.020221343237, 8: 0.015623024871}
say(f"  {'q':>4}{'lane chi':>20}{'this reduction':>20}{'|diff|':>12}")
worst = 0.0
for q in sorted(LANE):
    c = chi_crowd(q, 0.8, env); worst = max(worst, abs(c - LANE[q]))
    say(f"  {q:>4}{LANE[q]:>20.12f}{c:>20.12f}{abs(c-LANE[q]):>12.2e}")
say(f"  worst disagreement with the lane's published numbers: {worst:.2e}   (their FLOOR-M = 3e-15)")

say("")
say("  THE CURVE EXTENDED.  Same venue, lam = 0.8, q = 2 .. 200.")
QS = [2,3,4,5,6,7,8,10,12,14,16,20,24,28,32,40,48,56,64,80,96,112,128,160,200]
CH = {q: chi_crowd(q, 0.8, env) for q in QS}
say(f"  {'q':>6}{'chi':>18}{'q*chi':>14}{'local slope d ln chi / d ln q':>32}")
for i, q in enumerate(QS):
    if i == 0: sl = ""
    else:
        qa, qb = QS[i - 1], q
        sl = f"{(np.log(CH[qb])-np.log(CH[qa]))/(np.log(qb)-np.log(qa)):>32.4f}"
    say(f"  {q:>6}{CH[q]:>18.12f}{q*CH[q]:>14.9f}{sl}")

say("")
say("  POWER-LAW FITS OVER SLIDING WINDOWS.  If -1.5438 were an exponent these would agree.")
say(f"  {'window in q':<18}{'fitted p':>12}{'quoted 1-sigma':>16}{'rms resid (log)':>18}")
for lo, hi in [(2, 8), (4, 16), (8, 32), (16, 64), (32, 128), (64, 200), (2, 200)]:
    qs = [q for q in QS if lo <= q <= hi]
    p, s, r = fitp(qs, [CH[q] for q in qs])
    say(f"  {f'{lo} .. {hi}':<18}{p:>12.4f}{s:>16.4f}{r:>18.2e}")
say("  the lane's window is 2 .. 8.")

say("")
say("  IS THE 'EXPONENT' EVEN STABLE ACROSS THE VENUE'S OWN SCALE (D-17)?  Refit q = 2..8 only.")
VENUES = [("V0 baseline, lam 0.8", BASE.env(3), 0.8),
          ("lam 0.4", BASE.env(3), 0.4),
          ("lam 1.2", BASE.env(3), 1.2),
          ("lam 2.5", BASE.env(3), 2.5),
          ("beta 1.2", Venue("b", BASE.energies, BASE.times, beta=1.2).env(3), 0.8),
          ("energies (0.4,0.9,1.7)", Venue("e", (0.4, 0.9, 1.7, 1.0), BASE.times).env(3), 0.8),
          ("times [3,15]", BASE.env(3), 0.8),
          ("bath 1 qubit", Venue("q1", BASE.energies, BASE.times).env(1), 0.8)]
say(f"  {'venue':<28}{'p (q=2..8)':>14}{'quoted 1-sigma':>16}{'p (q=32..200)':>16}")
for nm, e, lam in VENUES:
    tt = np.linspace(3.0, 15.0, 25) if "3,15" in nm else TIMES
    small = [chi_crowd(q, lam, e, tt) for q in range(2, 9)]
    p1, s1, _ = fitp(list(range(2, 9)), small)
    big_qs = [32, 40, 48, 56, 64, 80, 96, 112, 128, 160, 200]
    big = [chi_crowd(q, lam, e, tt) for q in big_qs]
    p2, _, _ = fitp(big_qs, big)
    say(f"  {nm:<28}{p1:>14.4f}{s1:>16.4f}{p2:>16.4f}")

say("")
say(f"  elapsed {time.time()-t0:.1f}s")
