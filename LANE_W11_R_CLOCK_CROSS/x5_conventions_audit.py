# X5 — CUSTODY: DOES LANE C's CODE USE THE CONNECTIONS ITS PUBLISHED_CONVENTIONS DECLARES?
# The program's rule: "Unreproducible numbers are treated as absent.  Publish conventions, seeds,
# grid, code."  Lane C's PUBLISHED_CONVENTIONS declares TWO connections and says where each is used.
import numpy as np, re, sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W11_R_CLOCK")
from w11c_lib import K1, B0b, ops, generic_conn

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_W11_R_CLOCK"
print("DECLARED (PUBLISHED_CONVENTIONS.txt):")
print("   Primary : generic_conn -> (f,c) = (1.0, sqrt2)  [the corpus's only generic point]")
print("   Secondary (LEG 3A ONLY): the registrar's a = (1.0,0.37,0.91,sqrt2,0.23,1.77)")
print("\nWHAT THE CODE ACTUALLY DOES, per file:")
for fn in ("w11c_1_lattice.py","w11c_2_clockrays.py","w11c_3_rate_decomp.py",
           "w11c_4_gauge.py","w11c_5_correspondences.py","w11c_7_closedform.py"):
    src = open(f"{LANE}/{fn}").read()
    hits = [l.strip() for l in src.split("\n")
            if re.search(r"^\s*a[KB]?\s*=|generic_conn\(", l) and "def " not in l]
    print(f"   {fn:<28} {hits}")

K, B = K1(), B0b()
aK_leg1 = np.array([1.0, 0.37, 0.91, 2**0.5, 0.23, 1.77])
aB_leg1 = np.random.default_rng(20260817).uniform(0, 2*np.pi, 18)
for nm, Kx, a in (("leg1 K1", K, aK_leg1), ("leg1 B0b", B, aB_leg1)):
    *_ , WF, WC = ops(Kx, a)
    print(f"   {nm:<10} (f,c) = ({np.angle(WF)%(2*np.pi):.9f}, {np.angle(WC)%(2*np.pi):.9f})")
for Kx in (K, B):
    *_ , WF, WC = ops(Kx, generic_conn(Kx, np.random.default_rng(7+Kx.nv)))
    print(f"   declared primary on {Kx.name:<4} (f,c) = ({np.angle(WF)%(2*np.pi):.9f}, "
          f"{np.angle(WC)%(2*np.pi):.9f})")
print("""
FINDING X5-1.  LEG 1 -- the leg that produces the lane's headline set-equality, its 81/676 and
63/676 counts, its 9.18e-01 / 9.04e-01 off-sublattice spreads, and its named operative variable --
USES NEITHER DECLARED CONNECTION.  On K1 it uses the SECONDARY connection, which the conventions
restrict to "leg 3A only"; on B0b it uses a THIRD connection, an undeclared uniform(0,2pi) draw
that appears nowhere in PUBLISHED_CONVENTIONS.  The lane's own self_flag says "ONE pi PER CARRIER
AND ONE GENERIC CONNECTION ... Leg 3A additionally uses the registrar's".  There are THREE.
This does not change leg 1's arithmetic (its content is a diagonality theorem, connection-free),
but the conventions block is false as written and the B0b connection is unledgered.

FINDING X5-2.  LEG 1's B0b STATE FAMILY IS NOT REPRODUCIBLE FROM ITS PUBLISHED SEED ALONE.
w11c_1_lattice.py draws both state families from ONE module-level rng = default_rng(20260817),
so the B0b family depends on the K1 block having consumed 64 dirichlet+uniform draws first.
Re-running the B0b block alone gives a different family (my x1 got only-F-mid-loop min 1.566e-01
against the lane's 1.667e-01; only-C 3.271e-01 against 3.669e-01).  The set-equality, the counts
and the orders of magnitude reproduce exactly; the quadrant minima do not, and those minima are
quoted as figures in the lane's what_survives.
""")
