"""B3 -- THE NOISE FLOOR, ESTABLISHED BEFORE ANY RESIDUAL IS CALLED A SIGNAL.

Two floors, and they are different objects:

  FLOOR-M  THE MACHINE FLOOR.  Re-run the IDENTICAL mathematical quantity by a route float64
           cannot distinguish in principle but does in practice: relabel the bath qubits and
           permute their energies together (an exact symmetry of the whole construction).  The
           spread is pure float64.

  FLOOR-V  THE VENUE FLOOR.  Re-run the whole pipeline -- battery, fit, residual -- in venues
           that differ in bath energy assignment, energy values, inverse temperature, time
           window and time sampling.  chi itself moves a lot between venues; the question is
           whether the RESIDUAL of the explained model moves with it.  A residual that scatters
           across venues is venue noise.  A residual that is reproduced across venues is
           structure, and only then is it allowed to be called one.
"""
import numpy as np, sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from battery import *

t0 = time.time()
NB = 3
LAM = 0.8
N = 8
say("=" * 126)
say("B3   THE NOISE FLOOR")
say("=" * 126)
ops, _, _ = build_ops(N)

# ---------------------------------------------------------------- FLOOR-M
say("")
say("FLOOR-M   MACHINE FLOOR.  The same exact quantity computed twice, with the bath qubits")
say("          relabelled (0,1,2)->(1,2,0) and their energies permuted to match.  Exactly equal")
say("          by symmetry; any difference is float64.")
E = BASE.energies[:NB]
vA = Venue("A", E, BASE.times)
vB = Venue("B", (E[2], E[0], E[1]), BASE.times)     # energies rolled to follow the relabelling
envA, envB = vA.env(NB), vB.env(NB)
say(f"  {'configuration':<38}{'chi (labels 0,1,2)':>22}{'chi (labels 1,2,0)':>22}{'|difference|':>16}")
worstM = 0.0
for name, partners in CONFIGS:
    if cfg_maxq(partners) > N - 2: continue
    if any(s >= NB for _, s in partners): continue
    a = float(np.mean(chi_times([(ops[READ][0], 0)] + [(ops[l][0], s) for l, s in partners],
                                ops[READ][0], envA, LAM, vA.times)))
    b = float(np.mean(chi_times([(ops[READ][0], 1)] + [(ops[l][0], (s + 1) % NB) for l, s in partners],
                                ops[READ][0], envB, LAM, vB.times)))
    worstM = max(worstM, abs(a - b))
    if len(partners) <= 2 or abs(a - b) > 1e-14:
        say(f"  {name:<38}{a:>22.15f}{b:>22.15f}{abs(a-b):>16.2e}")
say("")
say(f"  FLOOR-M (machine floor, worst over the whole battery) = {worstM:.3e}")

# ---------------------------------------------------------------- FLOOR-V
VENUES = [
    Venue("V0 baseline      E=(1.0,1.4,0.7) t in [1,13] x25  b=2.0",  (1.0, 1.4, 0.7), np.linspace(1, 13, 25), 2.0),
    Venue("V1 energies permuted E=(0.7,1.0,1.4)",                     (0.7, 1.0, 1.4), np.linspace(1, 13, 25), 2.0),
    Venue("V2 energies changed  E=(0.85,1.25,1.55)",                  (0.85, 1.25, 1.55), np.linspace(1, 13, 25), 2.0),
    Venue("V3 time window       t in [3,15] x25",                     (1.0, 1.4, 0.7), np.linspace(3, 15, 25), 2.0),
    Venue("V4 time sampling     t in [1,13] x37",                     (1.0, 1.4, 0.7), np.linspace(1, 13, 37), 2.0),
    Venue("V5 temperature       beta = 1.2",                          (1.0, 1.4, 0.7), np.linspace(1, 13, 25), 1.2),
    Venue("V6 seed 12345 random E, random t window",
          tuple(np.random.default_rng(12345).uniform(0.6, 1.7, 4)),
          np.sort(np.random.default_rng(12345).uniform(1.0, 15.0, 25)), 2.0),
    Venue("V7 seed 987 random E, random t window",
          tuple(np.random.default_rng(987).uniform(0.6, 1.7, 4)),
          np.sort(np.random.default_rng(987).uniform(1.0, 15.0, 25)), 2.0),
]
say("")
say("FLOOR-V   VENUE FLOOR.  The whole pipeline re-run and REFIT in each venue.")
say("")
RES = {}
say(f"  {'venue':<52}{'chi alone':>12}{'a':>8}{'gamma':>8}{'beta':>8}{'delta':>8}{'c0':>10}{'rms resid':>12}{'max|resid|':>12}")
for V in VENUES:
    rows = run_battery(N, V, LAM, NB=NB, ops=ops)
    F = fit(rows)
    RES[V.name] = (rows, F)
    alone = [r['chi'] for r in rows if r['name'] == 'alone'][0]
    say(f"  {V.name:<52}{alone:>12.6f}" + "".join(f"{c:>8.4f}" for c in F['coef'][:4]) +
        f"{F['coef'][4]:>10.2e}{F['rms']:>12.3e}{F['maxabs']:>12.3e}")
say("")
say(f"  chi(alone) ranges over {min(r[0][0]['chi'] for r in [(RES[V.name][0],) for V in VENUES]):.4f} .. "
    f"{max(r[0][0]['chi'] for r in [(RES[V.name][0],) for V in VENUES]):.4f}"
    "   -- the venues really are different venues (D-17).")

# residual per configuration across venues
names = [r['name'] for r in RES[VENUES[0].name][0]]
say("")
say("  RESIDUAL PER CONFIGURATION, VENUE BY VENUE.  Read the SPREAD column: that is FLOOR-V for")
say("  that configuration.  Read the MEAN column against it.")
say(f"  {'configuration':<38}" + "".join(f"{'V'+str(i):>10}" for i in range(len(VENUES))) +
    f"{'mean':>11}{'spread':>10}{'|mean|/spread':>15}")
ratios = []
for j, nm in enumerate(names):
    vals = np.array([RES[V.name][1]['resid'][j] for V in VENUES])
    sprd = vals.max() - vals.min()
    r = abs(vals.mean()) / sprd if sprd > 0 else np.inf
    ratios.append((nm, vals.mean(), sprd, r))
    say(f"  {nm:<38}" + "".join(f"{v:>+10.4f}" for v in vals) +
        f"{vals.mean():>+11.5f}{sprd:>10.5f}{r:>15.2f}")
say("")
allsp = [x[2] for x in ratios]
say(f"  FLOOR-V, median configuration spread across venues : {np.median(allsp):.5f}")
say(f"  FLOOR-V, largest configuration spread             : {max(allsp):.5f}")
say(f"  rms residual, averaged over venues                : {np.mean([RES[V.name][1]['rms'] for V in VENUES]):.5f}")

# ---------------------------------------------------------------- reproducibility of the SIGN
say("")
say("  IS THE RESIDUAL REPRODUCED, OR DOES IT SCATTER?  A residual that is real must keep its")
say("  SIGN across venues that share nothing but the algebra.")
say(f"  {'configuration':<38}{'signs agree in':>18}{'of':>5}{'venues':>8}{'mean resid':>13}")
agree = 0
for j, nm in enumerate(names):
    vals = np.array([RES[V.name][1]['resid'][j] for V in VENUES])
    s = np.sign(vals)
    k = max((s > 0).sum(), (s < 0).sum())
    if k == len(VENUES): agree += 1
    say(f"  {nm:<38}{k:>18}{'/':>5}{len(VENUES):>8}{vals.mean():>+13.5f}")
say("")
say(f"  {agree} of {len(names)} configurations keep the SAME SIGN of residual in all {len(VENUES)} venues.")

# ---------------------------------------------------------------- three-body across venues
say("")
say("  THE MODEL-FREE THREE-BODY CONTRAST, ACROSS VENUES.  No fit involved.")
PAIRS3 = [("X2@0,X3@0", "X2@0,Z2@0"),
          ("Z1@0,Z1X2@0", "Z1X2@0,Z1Z2@0"),
          ("X2@0,X3@0,X4@0", "X2@0,Z2@0,X3@0"),
          ("Z1@0,X2@0,X3@0", "Z1X2@0,Z1Z2@0,X3@0"),
          ("X2@1,X3@2", "X2@1,X3@1"),
          ("Z1@1,Z1X2@2", "Z1@1,Z1X2@1")]
say(f"  {'contrast (A vs B)':<52}" + "".join(f"{'V'+str(i):>11}" for i in range(len(VENUES))) + f"{'sign stable':>13}")
for a, b in PAIRS3:
    ds = []
    for V in VENUES:
        rows = RES[V.name][0]
        ca = [r['chi'] for r in rows if r['name'] == a]
        cb = [r['chi'] for r in rows if r['name'] == b]
        ds.append(ca[0] - cb[0] if ca and cb else np.nan)
    ds = np.array(ds)
    st = "YES" if (np.all(ds > 0) or np.all(ds < 0)) else "no"
    say(f"  {a+'  vs  '+b:<52}" + "".join(f"{d:>+11.6f}" for d in ds) + f"{st:>13}")

json.dump({V.name: dict(names=names, resid=list(map(float, RES[V.name][1]['resid'])),
                        chi=[float(r['chi']) for r in RES[V.name][0]],
                        coef=list(map(float, RES[V.name][1]['coef'])))
           for V in VENUES},
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'b3_venues.json'), 'w'), indent=0)
say("")
say(f"  elapsed {time.time()-t0:.1f}s")
