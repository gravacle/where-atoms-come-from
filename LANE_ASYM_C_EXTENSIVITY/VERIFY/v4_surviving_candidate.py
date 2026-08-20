"""V4 -- WHAT SURVIVES.  The lane's own numbers contain a quantity that meets (a), (b) AND (c)
as the brief states them, which the lane discarded by calling it 'a count'.  V2 already showed
it is not a count.  Here it is put through the full three-requirement test, with controls.

CANDIDATE:  Var(Phi) -- the variance, over record configurations s, of the bath-induced
            effective potential Phi(s) on the record register.  Units of energy^2; the
            corresponding energy scale is std(Phi).

  (a) strict extensivity      S(2N)/S(N) -> 2
  (b) additivity over disjoint regions   defect = 0
  (c) not saturating, not topological    varies continuously with the couplings

CONTROLS in the same tables:
  * N itself (a genuine count): flat under every coupling change.
  * the SHARED-site configuration: additivity defect non-zero, so the defect metric is live.
  * 1 record per site: Var(Phi) identically zero -- a dead configuration, shown so that a
    non-zero elsewhere is not an artefact of the estimator.
"""
import sys, math, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY")
from qcore import phi_moments, chi_of_n, loglog_fit, LAM, BETA

OUT = []
def P(s=""):
    print(s); OUT.append(s)

P("=" * 112)
P("V4  THE CANDIDATE THE LANE DISCARDED:  Var(Phi) / std(Phi) on a heterogeneous body")
P("=" * 112)

def blk_var(e, lam, k=2):  return phi_moments(k, e=e, lam=lam)[1] ** 2

# a heterogeneous 'lump': 4 blocks with different local content, 2 records each
LUMP = [(0.5, 0.4), (1.0, 0.8), (2.0, 1.6), (1.4, 0.6)]

P("\n" + "-" * 112)
P("TABLE V4-1  --  (a) and (b) for a HETEROGENEOUS body: r copies of a fixed 4-block lump,")
P("                each block on its OWN bath site.  N = 8r records.")
P("-" * 112)
P("%-7s %-8s %-18s %-13s %-15s %-16s %-15s"
  % ("r", "N", "Var(Phi)", "S(2N)/S(N)", "additivity def", "CTRL shared-site", "CTRL N (count)"))
P("-" * 112)
prev = None; Ns = []; Vs = []
for r in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]:
    cfg = LUMP * r
    V = sum(blk_var(e, l) for e, l in cfg)
    ratio = (V / prev) if prev else float("nan"); prev = V
    h = len(cfg) // 2
    defect = V - sum(blk_var(e, l) for e, l in cfg[:h]) - sum(blk_var(e, l) for e, l in cfg[h:])
    N = 2 * len(cfg)
    # CONTROL: the same N records crowded onto ONE bath site -- additivity must FAIL there
    vs_full = phi_moments(N)[1] ** 2
    vs_half = 2 * (phi_moments(N // 2)[1] ** 2)
    P("%-7d %-8d %-18.10f %-13.6f %-15.3e %-16.6e %-15d"
      % (r, N, V, ratio, defect, vs_full - vs_half, N))
    Ns.append(N); Vs.append(V)
p, sg, rs = loglog_fit(Ns[-6:], Vs[-6:])
P("   log-log fit over the last 6 points: exponent p = %.9f, spread of fit sigma %.2e, max resid %.2e"
  % (p, sg, rs))
P("   READ: (a) ratio is exactly 2.000000 at every r; (b) defect is 0 to 1e-11 while the")
P("   shared-site CONTROL column is non-zero and growing, so the defect metric is live.")

P("\n" + "-" * 112)
P("TABLE V4-2  --  (c) NOT TOPOLOGICAL.  Same body, same N, couplings deformed continuously.")
P("-" * 112)
P("%-10s %-10s %-20s %-20s %-14s" % ("lam scale", "N (CTRL)", "Var(Phi)", "std(Phi)", "d ln V/d ln lam"))
P("-" * 112)
r = 64; prevV = None; prevs = None
for sc in [0.5, 0.6, 0.75, 1.0, 1.25, 1.5, 2.0]:
    cfg = [(e, l * sc) for e, l in LUMP] * r
    V = sum(blk_var(e, l) for e, l in cfg)
    slope = (math.log(V / prevV) / math.log(sc / prevs)) if prevV else float("nan")
    P("%-10.2f %-10d %-20.10f %-20.10f %-14.4f" % (sc, 2 * len(cfg), V, math.sqrt(V), slope))
    prevV, prevs = V, sc
P("   READ: N is exactly constant (that is what a count/topological invariant does); Var(Phi)")
P("   moves smoothly, roughly as lam^2.  It is not topological.  It passes (c).")

P("\n" + "-" * 112)
P("TABLE V4-3  --  THE DEAD-CONFIGURATION CONTROL.  1 record per bath site: Phi is EVEN in")
P("                each s_i, so Var(Phi) is identically 0.  Printed beside the live case.")
P("-" * 112)
P("%-10s %-24s %-24s" % ("N", "Var(Phi), 1 rec/site", "Var(Phi), 2 rec/site (live)"))
P("-" * 112)
for N in [8, 32, 128, 512, 2048]:
    dead = N * (phi_moments(1)[1] ** 2)
    live = (N / 2) * (phi_moments(2)[1] ** 2)
    P("%-10d %-24.3e %-24.10f" % (N, dead, live))
P("   READ: the estimator returns exactly 0 where the physics says 0, and a stable non-zero")
P("   where the physics says non-zero.  D-15 satisfied.")

P("\n" + "=" * 112)
P("READ OF V4")
P("=" * 112)
P(" * Var(Phi) on disjoint regions satisfies (a) exactly, (b) exactly, and (c) -- it is a")
P("   continuous function of the local couplings, not a topological invariant and not a count.")
P(" * The lane classified this SAME quantity as 'GROWING linear ... but a count'.  V2 TABLE")
P("   PART 1 shows it moves by a factor of 72.9 at FIXED N.  The classification is wrong.")
P(" * This does NOT make Var(Phi) gravity.  It has the wrong units (energy^2), it is defined")
P("   relative to a chosen coupling operator, and nothing here shows it sources any geometry.")
P("   What it does is remove the lane's stated ground for excluding it.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY/VERIFY/v4_surviving_candidate.txt",
     "w").write("\n".join(OUT) + "\n")
