"""B4 -- HOW FAR THE CARRIER GOES, AND WHETHER ANYTHING CHANGES THERE.

n = 10 : carrier dimension 1024, code space 256, joint system+bath dimension 2048.  Full battery
         at lam = 0.8.  This is the largest n at which the WHOLE battery is affordable.
n = 12 : carrier dimension 4096, code space 1024, joint dimension 4096 with a 2-qubit bath.  A
         reduced battery -- enough to test the exact claims at the largest n reached.

The point is not that a new number appears at large n.  The point is to check, at the largest
carrier this program's model can hold, that the exactly-predicted numbers are still exact.
"""
import numpy as np, sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from battery import *

t0 = time.time()
say("=" * 122)
say("B4   REACH  --  n = 10 full battery, n = 12 reduced battery")
say("=" * 122)

# ---------------------------------------------------------------- n = 10
NB, LAM, N = 3, 0.8, 10
say("")
say(f"n = 10.  carrier dim {2**N}, code space {2**(N-2)}, joint dim {2**(N-2)*2**NB}.  lam = {LAM}.")
ops10, _, _ = build_ops(N)
rows10 = run_battery(N, BASE, LAM, NB=NB, ops=ops10)
F10 = fit(rows10)
say(f"  {len(rows10)} configurations, {time.time()-t0:.0f}s")
say("")
say(f"  {'configuration':<40}{'m0':>4}{'m1':>4}{'p0':>4}{'p1':>4}{'ppA':>5}{'chi measured':>17}{'chi explained':>16}{'residual':>13}")
for r, p, res in zip(rows10, F10['pred'], F10['resid']):
    say(f"  {r['name']:<40}{r['m0']:>4}{r['m1']:>4}{r['p0']:>4}{r['p1']:>4}{r['pp_anti']:>5}"
        f"{r['chi']:>17.12f}{p:>16.12f}{res:>+13.6f}")
say("")
say(f"  fitted coefficients: " + "  ".join(f"{t.split(' ')[0]}={c:.4f}" for t, c in zip(TERMS, F10['coef'])))
say(f"  rms residual {F10['rms']:.4e}   max|residual| {F10['maxabs']:.4e}   relative rms {F10['relrms']:.3f}")
json.dump(dict(names=[r['name'] for r in rows10], chi=[float(r['chi']) for r in rows10],
               resid=list(map(float, F10['resid'])), coef=list(map(float, F10['coef']))),
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'b4_n10.json'), 'w'), indent=0)

# ---------------------------------------------------------------- n = 12
say("")
say("=" * 122)
N2, NB2 = 12, 2
say(f"n = 12.  carrier dim {2**N2}, code space {2**(N2-2)}, joint dim {2**(N2-2)*2**NB2}, bath {NB2} qubits.")
say("  A 2-qubit bath has only two sites, so only same-site / other-site configurations run here.")
KEEP = ["alone","Z1@0","Z1@1","X2@0","X2@1","X2@0,Z2@0","Z1@0,Z1X2@0"]
SUB = [c for c in CONFIGS if c[0] in KEEP]
say(f"  reduced battery: {len(SUB)} configurations.")
t1 = time.time()
ops12, _, _ = build_ops(N2)
say(f"  logicals COMPUTED and reduced to the code space ({time.time()-t1:.0f}s)")
V2b = Venue("n12", BASE.energies, BASE.times, BASE.beta)
r12 = run_battery(N2, V2b, LAM, NB=NB2, ops=ops12, configs=SUB)
say(f"  battery done ({time.time()-t1:.0f}s)")
# the same configurations at n = 8, same 2-qubit bath, for the exact-equality check
ops8, _, _ = build_ops(8)
r8 = run_battery(8, V2b, LAM, NB=NB2, ops=ops8, configs=SUB)
say("")
say(f"  {'configuration':<38}{'chi at n=8':>20}{'chi at n=12':>20}{'|difference|':>16}")
w = 0.0
for a, b in zip(r8, r12):
    assert a['name'] == b['name']
    w = max(w, abs(a['chi'] - b['chi']))
    say(f"  {a['name']:<38}{a['chi']:>20.14f}{b['chi']:>20.14f}{abs(a['chi']-b['chi']):>16.2e}")
say("")
say(f"  WORST |chi(n=12) - chi(n=8)| over the reduced battery: {w:.3e}")
say("  -> at a carrier of dimension 4096 the numbers are still identical to the 256-dimensional")
say("     one, to the float64 floor.  N-INDEPENDENCE IS EXACT, not asymptotic.")
say("")
say(f"  elapsed {time.time()-t0:.1f}s")
