"""V1 -- ADVERSARIAL CHECK: is the lane's "THREE-RECORD residual" a three-record quantity at all?

THE ORDINARY EXPLANATION UNDER TEST.  In common.chi_times the coupling is
        lam * sum_i A_i (x) X_{site(i)}
so every operator written on bath site 0 enters ONLY through the SINGLE SUM
        T = R + S,      S = sum of the partner operators on site 0.
When every one of those partners COMMUTES with the read record R (which is true of every
configuration in the lane's three-record channel (a) and of its entire pure-crowding family),
R and S can be simultaneously diagonalised, and the bath sees nothing but the classical
random variable  lam*(r + s).  chi is then a ZERO-PARAMETER functional of the eigenvalue
DISTRIBUTION of S -- no record count, no occupancy, no pairing count, and nothing three-body.

The distribution of S is fixed by elementary Pauli algebra:
   two COMMUTING partners  X2+X3   -> s in {+2, 0, 0, -2}
   two ANTICOMMUTING ones  X2+Z2   -> s in {+sqrt2, -sqrt2}          (since (X2+Z2)^2 = 2 I)
That is the WHOLE of "channel (a)".  If the zero-parameter spectral model reproduces the
lane's numbers to the float64 floor, the lane's "missing three-record quantity" is not a
three-body interaction between records; it is the spectrum of a sum of Paulis.

CONTROL (D-15): the same predictor is run on configurations it should NOT reproduce
(partners that ANTICOMMUTE with the read record, e.g. Z1@0), where R and S do not commute
and the reduction is invalid.  Those rows must MISS.
"""
import numpy as np, sys, os, time, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from common import *                      # noqa
from battery import build_ops, READ       # noqa

t0 = time.time()
N = 10
VEN = BASE
NB = 3
TIMES = VEN.times

say("=" * 118)
say("V1   IS THE 'THREE-RECORD RESIDUAL' A THREE-RECORD QUANTITY, OR THE SPECTRUM OF A SUM OF PAULIS?")
say("=" * 118)

ops, Q, _ = build_ops(N)
R = ops[READ][0]
env = VEN.env(NB)
d = R.shape[0]

# ---------------------------------------------------------------- the zero-parameter predictor
def s_distribution(labels):
    """eigenvalue distribution of S = sum(partner ops), resolved inside each R eigensector.
       Returns {r: [(s, weight), ...]}.  Requires [R, S] = 0 (checked by the caller)."""
    S = sum(ops[l][0] for l in labels) if labels else np.zeros_like(R)
    out = {}
    for r in (+1, -1):
        P = (np.eye(d) + r * R) / 2
        # basis of the sector
        w, V = np.linalg.eigh(P)
        B = V[:, w > 0.5]                      # d/2 columns
        Sr = B.conj().T @ S @ B
        ev = np.linalg.eigvalsh(Sr)
        # bin identical eigenvalues
        ev = np.sort(np.real(ev))
        binned = []
        for x in ev:
            if binned and abs(x - binned[-1][0]) < 1e-9:
                binned[-1][1] += 1
            else:
                binned.append([float(x), 1])
        tot = sum(c for _, c in binned)
        out[r] = [(x, c / tot) for x, c in binned]
    return out

def _vn(rho):
    e = np.linalg.eigvalsh(rho); e = e[e > 1e-13]
    return float(-(e * np.log2(e)).sum())

def chi_spectral(dist, lam, times=TIMES, environ=None):
    """chi from the eigenvalue distribution ALONE.  Zero fitted parameters."""
    environ = env if environ is None else environ
    rho0 = environ.thermal()
    X0 = environ.site[0]
    # cache one eigendecomposition per distinct field value
    fields = sorted({round(lam * (r + s), 12) for r in dist for s, _ in dist[r]})
    prop = {}
    for f in fields:
        w, V = np.linalg.eigh(environ.HB + f * X0)
        prop[f] = (w, V)
    def evolve(f, t):
        w, V = prop[f]
        U = (V * np.exp(-1j * w * t)[None, :]) @ V.conj().T
        return U @ rho0 @ U.conj().T
    acc = []
    for t in times:
        rr = {}
        for r in dist:
            rr[r] = sum(p * evolve(round(lam * (r + s), 12), t) for s, p in dist[r])
        av = 0.5 * (rr[+1] + rr[-1])
        acc.append(max(_vn(av) - 0.5 * (_vn(rr[+1]) + _vn(rr[-1])), 0.0))
    return float(np.mean(acc))

def chi_full(partners, lam):
    so = [(R, 0)] + [(ops[l][0], s) for l, s in partners]
    return float(np.mean(chi_times(so, R, env, lam, TIMES)))

# ---------------------------------------------------------------- the test set
COMMUTING = [
    ("alone",                        []),
    ("X2@0",                         ["X2"]),
    ("X2@0,X3@0",                    ["X2", "X3"]),
    ("X2@0,Z2@0",                    ["X2", "Z2"]),                     # channel (a), B
    ("X2@0,X3@0,X4@0",               ["X2", "X3", "X4"]),
    ("X2@0,Z2@0,X3@0",               ["X2", "Z2", "X3"]),
    ("X2@0,X3@0,X4@0,X5@0",          ["X2", "X3", "X4", "X5"]),
    ("X2@0,Z2@0,X3@0,Z3@0",          ["X2", "Z2", "X3", "Z3"]),         # the 'factor 4.28' row
    ("X2..X6@0",                     ["X2", "X3", "X4", "X5", "X6"]),
    ("X2..X7@0",                     ["X2", "X3", "X4", "X5", "X6", "X7"]),
    ("X2..X8@0",                     ["X2", "X3", "X4", "X5", "X6", "X7", "X8"]),
    ("X2,Z2,X3,Z3,X4,Z4@0",          ["X2", "Z2", "X3", "Z3", "X4", "Z4"]),
]
CONTROL = [   # partners that ANTICOMMUTE with R: the reduction is INVALID, these must MISS
    ("Z1@0            (CONTROL)",    ["Z1"]),
    ("Z1@0,X2@0       (CONTROL)",    ["Z1", "X2"]),
    ("Z1@0,Z1X2@0     (CONTROL)",    ["Z1", "Z1X2"]),
]

say("")
say("  ZERO-PARAMETER SPECTRAL PREDICTION vs THE LANE'S OWN PIPELINE, lam = 0.4 / 0.8 / 1.2.")
say("  The predictor uses ONLY the eigenvalue distribution of S = sum(same-site partner ops).")
say("")
hdr = f"  {'configuration':<28}{'[R,S]=0?':>10}" + "".join(
    f"{'chi full ' + str(l):>20}{'chi spectral':>16}{'|diff|':>12}" for l in (0.8,))
say(f"  {'configuration':<28}{'[R,S]':>7}{'#recs':>7}" +
    "".join(f"{'lam='+str(l)+' full':>17}{'spectral':>14}{'|diff|':>11}" for l in (0.4, 0.8, 1.2)))
worst_ok, worst_ctrl = 0.0, 0.0
for name, labs in COMMUTING + CONTROL:
    S = sum(ops[l][0] for l in labs) if labs else np.zeros_like(R)
    comm = float(np.linalg.norm(R @ S - S @ R))
    dist = s_distribution(labs)
    cells = ""
    for lam in (0.4, 0.8, 1.2):
        cf = chi_full([(l, 0) for l in labs], lam)
        cs = chi_spectral(dist, lam)
        diff = abs(cf - cs)
        cells += f"{cf:>17.12f}{cs:>14.12f}{diff:>11.2e}"
        if "CONTROL" in name: worst_ctrl = max(worst_ctrl, diff)
        else: worst_ok = max(worst_ok, diff)
    say(f"  {name:<28}{comm:>7.1f}{len(labs)+1:>7}" + cells)

say("")
say(f"  worst |full - spectral| over the configurations where [R,S] = 0 : {worst_ok:.3e}")
say(f"  worst |full - spectral| over the CONTROL rows ([R,S] != 0)      : {worst_ctrl:.3e}")
say("  FLOOR-M (the lane's own float64 floor)                          : 3.0e-15")

# ---------------------------------------------------------------- the exact algebra behind it
say("")
say("  THE SPECTRA, IN EXACT INTEGER ARITHMETIC (Fraction-free: all entries of S^2 are integers).")
say("  This is what actually differs between the lane's 'three-record' configurations.")
from fractions import Fraction
say(f"  {'partner set on site 0':<26}{'S^2 = ?':>34}{'eigenvalues of S (distribution)':>46}")
for name, labs in COMMUTING[1:]:
    S = sum(ops[l][0] for l in labs)
    S2 = S @ S
    # is S^2 a c-number?
    cnum = np.allclose(S2, S2[0, 0] * np.eye(d), atol=1e-10)
    tag = f"{np.real(S2[0,0]):.0f} * I  (EXACT)" if cnum else "not a c-number"
    dist = s_distribution(labs)[+1]
    ds = ", ".join(f"{s:+.4f}:{p:.3f}" for s, p in dist)
    say(f"  {name:<26}{tag:>34}{ds:>46}")

say("")
say(f"  elapsed {time.time()-t0:.1f}s")
