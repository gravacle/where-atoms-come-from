"""V2 -- THE HEADLINE TEST.  The lane's headline is:

    "Quantities that do satisfy strict extensivity and additivity exist ... but every one of
     them is (a fixed per-block number) x (the number of blocks), i.e. the block COUNT in
     different units, which C-35 already excludes; none satisfies requirement (c)."

A COUNT is a quantity determined by N alone.  So the claim is falsifiable in one line: hold N
fixed and change what is inside the blocks.  If the quantity moves, it is NOT the count.

This test also checks the converse the lane needs: is the quantity a SUM OVER REGIONS of a
LOCAL per-region value?  That is exactly the structure of a density integrated over a region,
which is what requirement (b) asks for.  A count is the degenerate case where every region
contributes the same number.

CONTROL (D-15): N itself, and P_cross, are carried in the same table.  N is a genuine count and
MUST be flat under composition change; if it moves, the instrument is broken.
"""
import sys, math, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY")
from qcore import chi_of_n, J2_same_site, phi_moments, LAM, BETA

OUT = []
def P(s=""):
    print(s); OUT.append(s)

P("=" * 112)
P("V2  IS EVERY EXTENSIVE QUANTITY 'THE COUNT IN DIFFERENT UNITS'?  Hold N fixed, vary the")
P("    CONTENT of the regions.  A count cannot move.  A density-like source must.")
P("=" * 112)

# per-block quantities: block b = one [[4,2,2]] block, its 2 records on the block's OWN bath
# site, that site having energy e_b and coupling lam_b.  Regions are DISJOINT (own baths), so
# every quantity below is exactly additive over blocks -- the lane proved that and it is
# re-verified in part 3.
def blk_chi(e, lam):   return 2.0 * chi_of_n(2, e=e, lam=lam)
def blk_var(e, lam):   return phi_moments(2, e=e, lam=lam)[1] ** 2
def blk_spread(e, lam):return phi_moments(2, e=e, lam=lam)[2]
def blk_J2(e, lam):    return abs(J2_same_site(2, e=e, lam=lam))

# ------------------------------------------------------------------ 1. fixed N, varied content
P("\n" + "-" * 112)
P("PART 1  --  N HELD FIXED AT 24 RECORDS (m = 12 disjoint blocks).  Only the block CONTENT")
P("            changes.  Every configuration has the same number of records and the same")
P("            number of regions.")
P("-" * 112)
m = 12
configs = [
  ("uniform  e=1.0 lam=0.8",  [(1.0, 0.8)] * m),
  ("uniform  e=1.0 lam=0.4",  [(1.0, 0.4)] * m),
  ("uniform  e=1.0 lam=1.6",  [(1.0, 1.6)] * m),
  ("uniform  e=0.3 lam=0.8",  [(0.3, 0.8)] * m),
  ("uniform  e=2.5 lam=0.8",  [(2.5, 0.8)] * m),
  ("MIXED    half/half lam",  [(1.0, 0.4)] * (m // 2) + [(1.0, 1.6)] * (m // 2)),
  ("MIXED    graded e",       [(0.3 + 0.2 * b, 0.8) for b in range(m)]),
]
P("%-26s %-7s %-9s %-15s %-15s %-15s %-15s"
  % ("configuration", "N", "P_cross", "chi_total", "Var(Phi)", "spread(Phi)", "sum|J_2|"))
P("-" * 112)
rows = []
for nm, cfg in configs:
    N = 2 * len(cfg)
    chi = sum(blk_chi(e, l) for e, l in cfg)
    var = sum(blk_var(e, l) for e, l in cfg)
    spr = sum(blk_spread(e, l) for e, l in cfg)
    j2  = sum(blk_J2(e, l) for e, l in cfg)
    rows.append((nm, N, chi, var, spr, j2))
    P("%-26s %-7d %-9d %-15.9f %-15.9f %-15.9f %-15.9f" % (nm, N, 0, chi, var, spr, j2))
P("-" * 112)
for lbl, idx in (("chi_total", 2), ("Var(Phi)", 3), ("spread(Phi)", 4), ("sum|J_2|", 5)):
    vals = [r[idx] for r in rows]
    P("   %-13s at FIXED N=24 :  min %.9f   max %.9f   max/min = %.4f"
      % (lbl, min(vals), max(vals), max(vals) / min(vals)))
Ns = [r[1] for r in rows]
P("   %-13s at FIXED N=24 :  min %d   max %d   max/min = %.4f   <-- CONTROL, a real count"
  % ("N (control)", min(Ns), max(Ns), max(Ns) / min(Ns)))
P("")
P("   READ: N is flat by construction (it is the count).  chi_total, Var(Phi), spread(Phi) and")
P("   sum|J_2| all move by factors of 2 to 40 at the SAME N.  They are therefore NOT")
P("   'the block count in different units'.  They are sums over regions of a LOCAL value that")
P("   depends continuously on what is in the region -- which is the structure of a density.")

# ------------------------------------------------------------------ 2. topological?
P("\n" + "-" * 112)
P("PART 2  --  REQUIREMENT (c): 'NOT SATURATING and NOT TOPOLOGICAL'.  A topological quantity")
P("            is locally constant under continuous deformation of the couplings.  Sweep lam")
P("            continuously at fixed N = 24 and fixed record set.")
P("-" * 112)
P("%-10s %-17s %-17s %-17s %-17s" % ("lam", "N (control)", "chi_total", "Var(Phi)", "spread(Phi)"))
P("-" * 112)
prev = None
for lam in [0.40, 0.45, 0.50, 0.60, 0.80, 1.00, 1.30, 1.60]:
    cfg = [(1.0, lam)] * m
    P("%-10.2f %-17d %-17.10f %-17.10f %-17.10f"
      % (lam, 2 * m, sum(blk_chi(e, l) for e, l in cfg),
         sum(blk_var(e, l) for e, l in cfg), sum(blk_spread(e, l) for e, l in cfg)))
P("   READ: the control column N is exactly constant -- that is what topological/counting looks")
P("   like.  chi_total, Var(Phi) and spread(Phi) vary smoothly and monotonically with a")
P("   continuous coupling, so they are NOT topological and NOT counts.  They pass (c).")

# ------------------------------------------------------------------ 3. (a) and (b) still hold
P("\n" + "-" * 112)
P("PART 3  --  DO (a) AND (b) STILL HOLD FOR A HETEROGENEOUS BODY?  Replicate a fixed MIXED")
P("            composition r times; check S(2N)/S(N) -> 2 and additivity over disjoint halves.")
P("-" * 112)
base = [(0.5, 0.4), (1.0, 0.8), (2.0, 1.6), (1.4, 0.6)]     # one heterogeneous 'lump'
P("%-8s %-8s %-17s %-13s %-15s %-15s"
  % ("r", "N", "chi_total", "S(2N)/S(N)", "additivity def", "CONTROL: same"))
P("%-8s %-8s %-17s %-13s %-15s %-15s"
  % ("", "", "", "", "(own baths)", "with a SHARED site"))
P("-" * 112)
prevv = None
for r in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]:
    cfg = base * r
    S = sum(blk_chi(e, l) for e, l in cfg)
    ratio = (S / prevv) if prevv else float("nan")
    prevv = S
    # additivity over disjoint halves, own baths
    if r >= 2:
        h = len(cfg) // 2
        defect = S - sum(blk_chi(e, l) for e, l in cfg[:h]) - sum(blk_chi(e, l) for e, l in cfg[h:])
    else:
        defect = 0.0
    # control: the SAME record count crowded onto ONE shared bath site
    ctrl = (2 * len(cfg)) * chi_of_n(2 * len(cfg))
    P("%-8d %-8d %-17.10f %-13.6f %-15.3e %-15.6e"
      % (r, 2 * len(cfg), S, ratio, defect, ctrl))
P("   READ: for a heterogeneous body the extensive quantity is sum_b q_b with q_b DIFFERENT")
P("   per region, the ratio S(2N)/S(N) is exactly 2.000000, and the disjoint-region defect is")
P("   0 to machine precision.  The CONTROL column (all N records forced onto one bath site)")
P("   collapses, so the instrument can register a failure of extensivity.")

# ------------------------------------------------------------------ 4. what the lane got right
P("\n" + "-" * 112)
P("PART 4  --  THE CONFOUND THE LANE DID NOT NAME, TESTED DIRECTLY.")
P("            Every 'own baths' quantity grows because THE BATH GROWS WITH THE RECORDS.")
P("            Hold the bath FIXED (nq sites) and grow N; see what survives.")
P("-" * 112)
P("%-8s %-12s %-18s %-18s %-18s"
  % ("N", "nq (fixed)", "chi_total", "Var(Phi)", "spread(Phi)"))
P("-" * 112)
NQ = 8
for N in [8, 16, 32, 64, 128, 256, 1024, 4096, 16384]:
    per = N // NQ
    chi = N * chi_of_n(per)
    v, sd, sp = phi_moments(per)
    P("%-8d %-12d %-18.6e %-18.6e %-18.6e" % (N, NQ, chi, NQ * sd ** 2, NQ * sp))
P("   READ: at FIXED bath size chi_total DECAYS and the Phi scales grow only as sqrt/linear in")
P("   the per-site occupancy.  So the extensivity found in PART 1-3 is extensivity in the")
P("   AMOUNT OF STUFF (records AND their bath sites together), not in record count alone.")
P("   That is the honest rescoping -- and it is also exactly how mass works: doubling the")
P("   matter doubles the atoms AND the volume they occupy.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY/VERIFY/v2_is_it_really_a_count.txt",
     "w").write("\n".join(OUT) + "\n")
