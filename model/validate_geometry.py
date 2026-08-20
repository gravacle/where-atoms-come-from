"""EVERY SEALED GEOMETRY HEADLINE (the C-77 increments), REPRODUCED THROUGH THE ONE ENTRY POINT.

Each check recomputes a registered headline number through ProjectModel's GEOMETRY layer
(model/geometry.py, machinery ported from the sealed lanes) and gates it against the
SEALED value copied from the lane's own sealed OUT file.  Sealed sources per section:

  C-78   LANE_T42_A_DISTANCE/t42_a_toric.OUT.txt        (d_W class matrices, L = 2,3,4)
  C-79   LANE_T42_D_DERIVE/t42d_derive.OUT.txt          (IR2/capP sequences at L=12; the
                                                         corner form; world IFACE/CAP rows)
  C-80   LANE_O54_A_CORNER/o54a_two_region.OUT.txt      (separated zeros, seam law, winding)
         LANE_O54_C_ATTEMPT/o54c_attempts.OUT.txt       (w_min = d rows, venue (3,7))
  C-81   LANE_T43_A_CORNER/t43_corner.OUT.txt           (CERT = 8s-10 = IR2, STORED anchors,
                                                         WRITE0 = 0)
  C-83   LANE_T45_CLOCK/t45_clock.txt                   (CERT_W table; W=1 and W=inf limits)
  C-71/  LANE_T32_REALSOURCE/t32_source.txt             (orientation accumulation/screening)
  C-72   LANE_T34_NAND/t34_nand.txt                     (occupancy ratio 1; unwritten null)

Exit code: 0 iff every gate passes AND (unless --no-chain) validate_project.py still
returns its 11 PASS -- the full-model conjunction."""
import sys, os, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from project_model import ProjectModel

M = ProjectModel()
n_pass = 0
n_fail = 0


def check(name, cond, detail=""):
    global n_pass, n_fail
    if cond:
        n_pass += 1
        print(f"  PASS  {name}  {detail}")
    else:
        n_fail += 1
        print(f"  FAIL  {name}  {detail}")


print("VALIDATE THE GEOMETRY LAYER (C-77 increments) THROUGH ProjectModel")
print("=" * 78)

# ---------------------------------------------------------------- C-78: earned distance
# SEALED (t42_a_toric.OUT.txt): class minima {(1,0): L, (0,1): L, (1,1): 2L}, d_code = L,
# and the full 4x4 matrix in computed-label order 00,10,01,11.
for L in (2, 3, 4):
    dm = M.distance_matrix(L)
    sealed_cm = {(0, 0): 0, (1, 0): L, (0, 1): L, (1, 1): 2 * L}
    check(f"C-78 d_W class minima at L={L} == sealed {{{L}, {L}, {2*L}}}",
          dm["class_min"] == sealed_cm and dm["dcode"] == L,
          f"class_min={dm['class_min']} dcode={dm['dcode']}")
    order = dm["order"]
    sealed_mat = [[sealed_cm[(a[0] ^ b[0], a[1] ^ b[1])] for b in order] for a in order]
    got = [[M.distance(L, a, b) for b in order] for a in order]
    check(f"C-78 full 4x4 d_W matrix at L={L} == sealed matrix (via .distance)",
          got == sealed_mat, f"rows={got}")

# ---------------------------------------------------------------- C-79: corner tier
# SEALED (t42d_derive.OUT.txt, L=12): IR2(s) = [6,14,22,30,38,46,54] for s=2..8,
# capP(s) = [7,19,35,55,79,107,139]; IR2 = 8s-10; capP = 2s^2+2s-5; the corner form.
SEALED_IR2 = [6, 14, 22, 30, 38, 46, 54]
SEALED_CAPP = [7, 19, 35, 55, 79, 107, 139]
S_RANGE = list(range(2, 9))
rows = [M.interface_rank(12, s) for s in S_RANGE]
check("C-79 IR2(s), L=12, s=2..8 == sealed sequence",
      [r["IR2"] for r in rows] == SEALED_IR2, f"{[r['IR2'] for r in rows]}")
check("C-79 IR2 == 8s-10 on every s (the cut-rank boundary law)",
      all(r["IR2"] == 8 * s - 10 for s, r in zip(S_RANGE, rows)))
check("C-79 capP(s), L=12, s=2..8 == sealed sequence",
      [r["capP"] for r in rows] == SEALED_CAPP, f"{[r['capP'] for r in rows]}")
check("C-79 capP == 2s^2+2s-5 on every s",
      all(r["capP"] == 2 * s * s + 2 * s - 5 for s, r in zip(S_RANGE, rows)))
check("C-79 parameter-free corner form 32*capP == (IR2+10)^2 + 8*(IR2+10) - 160",
      all(32 * r["capP"] == (r["IR2"] + 10) ** 2 + 8 * (r["IR2"] + 10) - 160 for r in rows))
check("C-79 PER == 4s (the perimeter count behind LAW-2: IR2 = 2*PER - 10)",
      all(r["PER"] == 4 * s and r["IR2"] == 2 * r["PER"] - 10
          for s, r in zip(S_RANGE, rows)))

# ---------------------------------------------------------------- C-79: world tier
# SEALED (t42d_derive.OUT.txt, n=2..12; gated here on the asked range n=2..10):
SEALED_IFACE = [24, 54, 96, 150, 216, 294, 384, 486, 600]
SEALED_CAPD1 = [8, 26, 56, 98, 152, 218, 296, 386, 488]
N_RANGE = list(range(2, 11))
wrows = [M.world_interface(n) for n in N_RANGE]
check("C-79 world IFACE(n), n=2..10 == sealed sequence (= 6n^2)",
      [w["IFACE"] for w in wrows] == SEALED_IFACE and
      all(w["IFACE"] == 6 * n * n for n, w in zip(N_RANGE, wrows)),
      f"{[w['IFACE'] for w in wrows]}")
check("C-79 world CAP_d1(n) == sealed sequence == 6n^2-12n+8 == n^3-(n-2)^3",
      [w["CAP_d1"] for w in wrows] == SEALED_CAPD1 and
      all(w["CAP_d1"] == 6 * n * n - 12 * n + 8 == n ** 3 - (n - 2) ** 3
          for n, w in zip(N_RANGE, wrows)),
      f"{[w['CAP_d1'] for w in wrows]}")
check("C-79 world tier form IFACE^3 == 216*CAP_total^2 at every n (CAP_total == n^3)",
      all(w["IFACE"] ** 3 == 216 * w["CAP_total"] ** 2 and w["CAP_total"] == n ** 3
          for n, w in zip(N_RANGE, wrows)))

# ---------------------------------------------------------------- C-80: two-region taxonomy
# SEALED (o54a_two_region.OUT.txt): separated rows identically zero; single-contact seam
# I_IR = 2(s-2), I_MI = s-2; winding pair the exact constant w.
sep_cases = [(8, 2, 2), (8, 3, 1), (10, 4, 1)]
sep_rows = [M.mutual_interface(L, s, g) for (L, s, g) in sep_cases]
check("C-80 I(A:B) == 0 separated, three sealed cases (I_IR, I_MI, joint content all 0)",
      all(r["I_IR"] == 0 and r["I_MI"] == 0 and r["JC"] == 0 and
          r["JOINT_NEW"] == 0 and r["IDENT"] == 0 for r in sep_rows),
      f"cases={sep_cases}")
c3 = M.mutual_interface(8, 3, 0)
c4 = M.mutual_interface(10, 4, 0)
check("C-80 single-contact seam law I_IR == 2(s-2) for s=3,4 (sealed 2 and 4)",
      c3["I_IR"] == 2 and c3["I_MI"] == 1 and c4["I_IR"] == 4 and c4["I_MI"] == 2,
      f"s=3: I_IR={c3['I_IR']} I_MI={c3['I_MI']}; s=4: I_IR={c4['I_IR']} I_MI={c4['I_MI']}")
wnd = M.winding_interface(8, 2, 4)
check("C-80 winding constant: width-2 rings, L=8 -- I_IR == I_MI == IDENT == w == 2, "
      "JOINT_NEW == 0 (sealed: distance-independent constant)",
      wnd["I_IR"] == 2 and wnd["I_MI"] == 2 and wnd["IDENT"] == 2 and
      wnd["JOINT_NEW"] == 0 and wnd["capA"] == wnd["capB"] == wnd["capAB"] == 2,
      f"{ {k: wnd[k] for k in ('I_IR','I_MI','IDENT','JOINT_NEW')} }")
# SEALED (o54c_attempts.OUT.txt, venue (3,7), u=(0,0)): v=(0,1) d=1 w=1 N=1;
# v=(0,2) d=2 w=2 N=1; v=(0,3) d=3 w=3 N=1 -- exhaustive coset 2^22 per pair.
sealed_pairs = [((0, 1), 1, 1, 1), ((0, 2), 2, 2, 1), ((0, 3), 3, 3, 1)]
got_pairs = []
wd_ok = True
for v, sd, sw, sn in sealed_pairs:
    r = M.coupling_cost(3, 7, (0, 0), v)
    got_pairs.append((v, r["d_gen"], r["w_min"], r["N_min"]))
    wd_ok &= (r["d_gen"] == sd and r["w_min"] == sw and r["N_min"] == sn
              and r["w_min"] == r["d_gen"])
check("C-80 coupling-writer law w_min == d, three sealed pairs on the (3,7) venue "
      "(exhaustive coset scans)", wd_ok, f"(v, d, w, N)={got_pairs}")

# ---------------------------------------------------------------- C-81: certifiability
# SEALED (t43_corner.OUT.txt): squares -- CERT = 8s-10 (= the cut-rank identity),
# STORED anchor 19 at s=3; WRITE0 = 0 on every region (the locality theorem).
cert_ok = True
w0_ok = True
detail = []
for s in (3, 4, 5):
    c = M.certifiability(8, s)
    ir = M.interface_rank(8, s)
    cert_ok &= (c["CERT"] == 8 * s - 10 == ir["IR2"])
    w0_ok &= (c["WRITE0"] == 0)
    detail.append((s, c["CERT"], ir["IR2"], c["WRITE0"]))
check("C-81 CERT == 8s-10 == IR2 identity, three cases (s=3,4,5 at L=8)",
      cert_ok, f"(s, CERT, IR2, WRITE0)={detail}")
check("C-81 WRITE0 == 0 on all three regions (the locality theorem, computed)", w0_ok)
c33 = M.certifiability(8, 3)
check("C-81 sealed anchor: 3x3 STORED == 19, CERT == 14, HIDDEN == 5",
      c33["STORED"] == 19 and c33["CERT"] == 14 and c33["HIDDEN"] == 5,
      f"STORED={c33['STORED']} CERT={c33['CERT']}")

# ---------------------------------------------------------------- C-83: the window law
# SEALED (t45_clock.txt): the CERT_W table rows for W = 1, 10, 1000, inf at n = 10,100,1000.
SEALED_CW = {1: [488, 58808, 5988008],
             10: [1000, 588080, 59880080],
             1000: [1000, 1000000, 1000000000],
             float('inf'): [1000, 1000000, 1000000000]}
cw_ok = all([M.cert_window(n, W) for n in (10, 100, 1000)] == row
            for W, row in SEALED_CW.items())
check("C-83 CERT_W(n) = min(n^3, W*(6n^2-12n+8)): sealed table rows W=1,10,1000,inf",
      cw_ok, f"W=1 -> {[M.cert_window(n, 1) for n in (10, 100, 1000)]}")
check("C-83 W=1 limit reproduces C-82's per-epoch law L1(n) = 6n^2-12n+8 exactly",
      all(M.cert_window(n, 1) == 6 * n * n - 12 * n + 8 for n in (10, 100, 1000)))
check("C-83 W=inf limit reproduces the volume corner n^3 (DEF-A, t_m -> infinity)",
      all(M.cert_window(n, float('inf')) == n ** 3 for n in (10, 100, 1000)))

# ---------------------------------------------------------------- C-71/C-72: formation
# SEALED (t32_source.txt, seed 7): DC-SATURATED ratio 1.000000 at every N; random data
# sums +4 (N=10) and +270 (N=1e5); AC-erased -2 and +262; DC-free coded exactly 0.
import geometry as GE
ori = M.formation_orientation()
check("C-72 orientation: DC-SATURATED ratio == 1 at every sealed N (accumulates)",
      all(ori[("DC-SATURATED", N)]["ratio"] == 1.0
          for N in (10, 100, 1000, 10000, 100000)))
r10 = round(ori[("random data", 10)]["sum"] / GE.M_GRAIN)
r1e5 = round(ori[("random data", 100000)]["sum"] / GE.M_GRAIN)
a10 = round(ori[("AC-erased ctrl", 10)]["sum"] / GE.M_GRAIN)
a1e5 = round(ori[("AC-erased ctrl", 100000)]["sum"] / GE.M_GRAIN)
check("C-72 orientation: random data SCREENS -- sealed sums +4 (N=10), +270 (N=1e5) "
      "(ratios 0.400000, 0.002700)",
      r10 == 4 and r1e5 == 270 and
      abs(ori[("random data", 10)]["ratio"] - 0.4) < 1e-12 and
      abs(ori[("random data", 100000)]["ratio"] - 0.0027) < 1e-12,
      f"sums {r10}, {r1e5}")
check("C-72 orientation: AC-erased control screens -- sealed sums -2 (N=10), +262 (N=1e5)",
      a10 == -2 and a1e5 == 262, f"sums {a10}, {a1e5}")
check("C-72 orientation: DC-free coded pattern screens EXACTLY (sum 0, ratio 0)",
      ori[("DC-free coded", 10)]["sum"] == 0.0 and
      ori[("DC-free coded", 100000)]["sum"] == 0.0 and
      ori[("DC-free coded", 100000)]["ratio"] == 0.0)
# SEALED (t34_nand.txt, seed 11): written ratio exactly 1 for every pattern (min over
# 1000 draws == 1); unwritten sum -38 e vs tolerance 5000 e -> null within tolerance.
occ = M.formation_occupancy()
check("C-71 occupancy: written ratio == 1 for EVERY pattern (1000 sealed draws + "
      "all-programmed at N=1000) -- mechanism-constitutive one-signedness",
      min(occ["written_ratios"]) == 1.0 == max(occ["written_ratios"]) and
      occ["all_programmed"]["ratio"] == 1.0,
      f"min ratio {min(occ['written_ratios'])}")
e_sum = round(occ["unwritten"]["sum"] / GE.E_CHARGE)
check("C-71 occupancy: unwritten null -- sealed sum -38 e within declared tolerance "
      "+-5e/cell (5000 e)",
      e_sum == -38 and occ["null_within_tolerance"],
      f"sum {e_sum} e, tol {round(occ['tolerance'] / GE.E_CHARGE)} e")

# ---------------------------------------------------------------- summary + chain
print("=" * 78)
print(f"  GEOMETRY: {n_pass} PASS, {n_fail} FAIL")

chain_ok = True
if "--no-chain" not in sys.argv:
    print()
    print("CHAIN: validate_project.py (the existing 11 model gates must still pass)")
    print("-" * 78)
    sys.stdout.flush()   # keep the chained child's output in order in file captures
    here = os.path.dirname(os.path.abspath(__file__))
    r = subprocess.run([sys.executable, os.path.join(here, "validate_project.py")],
                       cwd=here)
    chain_ok = (r.returncode == 0)
    print("-" * 78)
    print(f"  CHAIN validate_project.py: {'PASS' if chain_ok else 'FAIL'}")

print("=" * 78)
overall = (n_fail == 0) and chain_ok
print(f"  OVERALL: {'PASS' if overall else 'FAIL'} "
      f"(geometry {n_pass}/{n_pass + n_fail}"
      + ("" if "--no-chain" in sys.argv else f", chain {'ok' if chain_ok else 'FAILED'}")
      + ")")
sys.exit(0 if overall else 1)
