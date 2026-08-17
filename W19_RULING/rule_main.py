# rule_main.py -- W-19 RULING.  Independent checks that decide the ruling.
import numpy as np, time
from rule_verify import *

t0 = time.time()
P("=" * 118)
P("W-19 RULING -- INDEPENDENT VERIFICATION.  Z_2 pure gauge, exact, numpy only, no lane code imported.")
P("=" * 118)

# ================================================================== 0  INSTRUMENT: reproduce the sealed T1 null
P("")
P("[0] INSTRUMENT VALIDATION -- reproduce the SEALED LANE_T1_NEW_PROGRAM null by a third route.")
P("    T1: theta graph, 3 links, G = X0 X1 X2, plaquettes W01 = Z0 Z1 and W12 = Z1 Z2,")
P("    H = -0.7 sum_l X_l - (W01 + W12), ground state in the physical sector.")
c = Carrier("theta_3(T1)", *theta(3))
c.plaq = [0b011, 0b110]; c.cyc = c.plaq
c.tvec = [sum((1 << i) for i, cc in enumerate(c.cyc) if cc >> l & 1) for l in range(c.L)]
D = c.dimP; H = np.zeros((D, D)); m = np.arange(D)
for i in range(c.C): H[m, m] += -(1.0 - 2.0 * ((m >> i) & 1))
for l in range(c.L): H[m ^ c.tvec[l], m] += -0.7
w, U = np.linalg.eigh(H); vec = c.lift(U[:, 0])
P("    dim_phys = %d   gauss residual on the lifted state = %.3e" % (c.dimP, gauss_residual(c, vec)))
def I2(vec, L, A, B):
    return S_ax(vec, L, AX(L, A)) + S_ax(vec, L, AX(L, B)) - S_ax(vec, L, AX(L, A + B))
P("    I(0:1) = %.6f    I(0:2) = %.6f    I(0:{1,2}) = %.6f" %
  (I2(vec, 3, [0], [1]), I2(vec, 3, [0], [2]), I2(vec, 3, [0], [1, 2])))
P("    SEALED T1 :  0.690763            0.384496            1.075259")

# ================================================================== 1  carrier table
P("")
P("[1] CARRIERS.  Built here from their definitions; d = dist_{G-l}(tail,head) for l = link 0.")
defs = [("theta_3", theta(3)), ("theta_6", theta(6)), ("theta_8", theta(8)),
        ("dbl_chain9", mg_chain(5)), ("tri_chain12", tri_chain12()),
        ("petersen", petersen()), ("heawood", heawood())]
cars = {}
P("    %-14s %3s %3s %3s %7s %7s %6s %4s   plaquette weights" %
  ("carrier", "V", "L", "C", "dim_ph", "mindeg", "girth", "d"))
for nm, (V, E) in defs:
    car = Carrier(nm, V, E); cars[nm] = car
    _, d = rule_A_fragments(V, E, 0)
    P("    %-14s %3d %3d %3d %7d %7d %6d %4d   %s" %
      (nm, car.V, car.L, car.C, car.dimP, car.mindeg, car.girth, d,
       sorted(bin(p).count("1") for p in car.plaq)))

# ================================================================== 2  Perron-Frobenius
P("")
P("[2] IS THE GAUSS-LAW PROJECTION INERT ON THE GROUND STATE?  (lane B refuter's Perron-Frobenius)")
P("    Full 2^L Hamiltonian, NO projection, vs the physical-sector Hamiltonian.  Tested here on")
P("    lane A's own min-degree-3 carriers, which the refuter did not test.")
P("    %-14s %5s %16s %16s %12s %10s" % ("carrier", "g2", "E0(full 2^L)", "E0(physical)", "|<free|proj>|", "min<G_v>"))
for nm in ["theta_6", "dbl_chain9", "tri_chain12"]:
    car = cars[nm]
    for g2 in (0.30, 1.00, 3.00):
        Hf = H_full(car, g2)
        wf, Uf = np.linalg.eigh(Hf); free = Uf[:, 0]
        wp, psi = car.ground(g2); proj = car.lift(psi)
        ov = abs(float(np.dot(free, proj)))
        idx = np.arange(len(free))
        mg = min(float(np.dot(free, free[idx ^ msk])) for msk in gauss_masks(car))
        P("    %-14s %5.2f %16.9f %16.9f %12.9f %10.6f" % (nm, g2, wf[0], wp[0], ov, mg))
P("    CONTROL THAT CAN FAIL: flip the sign of the electric term (H' = -(1/g2)sum W + g2 sum X).")
for nm in ["dbl_chain9", "tri_chain12"]:
    car = cars[nm]
    Dd = 1 << car.L; z = np.arange(Dd); Hf = np.zeros((Dd, Dd))
    for p in car.plaq:
        Hf[z, z] += -(1.0) * (1.0 - 2.0 * (np.bitwise_count(z & np.int64(p)) & 1))
    for l in range(car.L): Hf[z ^ (1 << l), z] += +1.0
    wf, Uf = np.linalg.eigh(Hf); free = Uf[:, 0]
    mg = min(float(np.dot(free, free[z ^ msk])) for msk in gauss_masks(car))
    P("    %-14s L=%2d  min<G_v> of the free ground state = %+.9f   physical? %s"
      % (nm, car.L, mg, mg > 0.999))

# ================================================================== 3  three channels on ground states
P("")
P("[3] THE THREE ALGEBRA CHANNELS ON GROUND STATES.  g2 = 0.50 (lane A's exhibit coupling).")
P("    Rule-A nested fragments.  ONE VARIABLE MOVES between EXT / CHI / CL: the system algebra.")
res = {}
for nm in ["theta_6", "theta_8", "dbl_chain9", "tri_chain12", "petersen", "heawood"]:
    car = cars[nm]
    wp, psi = car.ground(0.50); vec = car.lift(psi)
    gr = gauss_residual(car, vec)
    P("  -- %-12s L=%2d  E0=%.9f  gap=%.6f  gauss residual=%.2e" % (nm, car.L, wp[0], wp[1] - wp[0], gr))
    res[nm] = measure(car, vec, 0, nm + " ground g2=0.50")

# ================================================================== 4  small-g2 (magnetic) limit
P("")
P("[4] THE SAME CARRIERS DEEP IN THE MAGNETIC PHASE, g2 = 0.10.")
for nm in ["theta_6", "theta_8", "dbl_chain9", "tri_chain12"]:
    car = cars[nm]
    wp, psi = car.ground(0.10); vec = car.lift(psi)
    measure(car, vec, 0, nm + " ground g2=0.10")
P("    and the EXACT g2->0 state on theta_8: the magnetic GHZ (|0..0>+|1..1>)/sqrt(2).")
car = cars["theta_8"]; v = np.zeros(1 << 8); v[0] = v[(1 << 8) - 1] = 1 / np.sqrt(2)
P("      gauss residual = %.3e" % gauss_residual(car, v))
measure(car, v, 0, "theta_8 magnetic GHZ")

# ================================================================== 5  discrimination
P("")
P("[5] DOES THE CRITERION DISCRIMINATE?  Same carrier, same cut, same fragments, same estimator.")
P("    ONLY THE STATE MOVES: ground state vs three HAAR-RANDOM PHYSICAL states (complex Ginibre).")
for nm in ["tri_chain12", "theta_8"]:
    car = cars[nm]
    for s in (11, 22, 33):
        vec = car.lift(haar_physical(car, s))
        measure(car, vec, 0, "%s HAAR seed=%d" % (nm, s))
P("    and a SCRAMBLED-but-physical control on theta_8 that lane B used: Bell pair on links 0,1")
P("    lifted into the physical sector by gauge averaging.")
car = cars["theta_8"]; L = 8
raw = np.zeros([2] * L)
raw[0, 0] = 1 / np.sqrt(2); raw[1, 1] = 1 / np.sqrt(2)
raw = raw.reshape(-1)
idx = np.arange(1 << L); acc = np.zeros(1 << L)
for msk in [0] + gauss_masks(car):
    acc += raw[idx ^ msk]
acc /= np.linalg.norm(acc)
P("      gauss residual = %.3e" % gauss_residual(car, acc))
measure(car, acc, 0, "theta_8 Bell(0,1) gauge-averaged")

P("")
P("elapsed %.1f s" % (time.time() - t0))
open("OUT_rule_main.txt", "w").write("\n".join(LOG) + "\n")
