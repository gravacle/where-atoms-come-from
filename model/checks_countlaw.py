"""EVERY SEALED COUNT-LAW HEADLINE (C-86), GATED THROUGH model/countlaw.py.

Idiom: validate_geometry.py -- run_countlaw_checks(check) makes check(name, cond, detail="")
calls; the registrar chains it into a validator (see INTEGRATION_countlaw.md). Standalone:
`python3 checks_countlaw.py` runs the block with its own harness, exit 0 iff every gate passes.

Sealed sources, per section (each literal below is a SEALED ANCHOR, stated as such -- D-8:
computed-vs-computed comparisons run beside every anchor; no literal stands where a computed
comparison is possible):

  A  instrument fidelity            LANE_T47_A_WIDTH/t47_a_width.txt  section 0
  B  the two widths                 LANE_T47_A_WIDTH sections 1, 2, 4, 5;
                                    LANE_T47_B_STAIRCASE/t47b_staircase.txt section A
  C  the staircase / census         LANE_T47_A_WIDTH section 6; LANE_T47_B_STAIRCASE
                                    sections C, D, G
  D  the T-31 control carrier       LANE_T47_B_STAIRCASE section F;
                                    LANE_T31_ASYMMETRY/t31_asymmetry.txt (ratio class)
  E  the C-76 gate                  structural -- the chosen-width form is unreachable
  F  observation entry (D-25)       project_model.URM.surface -> countlaw.census
  G  API-fidelity probes BEYOND the gated range (definition-not-shortcut)

D-15: every zero below is gated with a positive control beside it, in the same or the
adjacent check (no-crossing beside crossing, k = 0 beside k > 0, v2 = 0 beside v2 = 2 = k,
refusal beside acceptance, departure 0 beside departure 5.4339, declined beside counted)."""
import sys
import os
import inspect
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import grounded as G
import countlaw as CL
from project_model import URM

EV, YEAR = CL.EV, CL.YEAR
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
SP = np.array([[0, 1], [0, 0]], dtype=complex)


def corner(name, B, dE, T, f0):
    """Declared corner-grid surface through the D-25 gate (T-55: no raw RecordSurface
       anywhere in this file). B from the lower well -> E_b = B - dE, the exact map."""
    return URM.surface(name, "thermal two-well (declared corner grid)", dE, B - dE, T, f0,
                       provenance="DEF-A", tier="corner")


def crossing_instr(B, T, f0, t_m):
    """Bisect the INSTRUMENT's record-mode rate against 1/t_m in dE, to fp spacing
       (LANE_T47_A_WIDTH crossing_num, ported as a probe). None where no crossing."""
    def rate(dE):
        return CL.record_rate(corner("probe", B, dE, T, f0))
    if rate(0.0) >= 1.0 / t_m:
        return None
    kT = G.KB * T
    naive = B - kT * np.log(f0 * t_m)
    lo, hi = 0.0, min(0.999 * B, max(naive, 0.0) + 40.0 * kT)
    if rate(hi) < 1.0 / t_m:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if rate(mid) <= 1.0 / t_m:
            lo = mid
        else:
            hi = mid
        if hi - lo <= np.spacing(hi):
            break
    return 0.5 * (lo + hi)


def rel(a, b):
    return abs(a - b) / abs(b)


def run_countlaw_checks(check):
    # ================================================================ A. INSTRUMENT
    # SEALED (t47_a_width.txt lines 16-19): rates at B = 1.2 eV (lower well), T = 300 K,
    # f0 = 1e9 Hz, printed to 12 digits.
    B0, T0, F0, TM0 = 1.2 * EV, 300.0, 1e9, 10.0 * YEAR
    kT0 = G.KB * T0
    SEALED_RATES = {0.0: 1.386491636889e-11, 0.05: 5.488926847605e-11,
                    0.158: 3.134238877190e-09, 0.30: 7.597477224034e-07}
    ok, worst = True, 0.0
    for dE_eV, sealed in SEALED_RATES.items():
        r = CL.record_rate(corner("a", B0, dE_eV * EV, T0, F0))
        worst = max(worst, rel(r, sealed))
        ok &= rel(r, sealed) < 1e-9
    check("CL-A1 record_rate == sealed instrument rates, 4 dE values (sealed anchors, "
          "t47_a_width.txt sec 0)", ok, f"worst rel {worst:.1e}")
    ok, worst = True, 0.0
    for dE_eV in SEALED_RATES:
        s = corner("a", B0, dE_eV * EV, T0, F0)
        gu = s.f0 * np.exp(-s.E_b / kT0)
        gl = s.f0 * np.exp(-(s.E_b + s.dE) / kT0)
        worst = max(worst, rel(CL.record_rate(s), gu + gl))
        ok &= rel(CL.record_rate(s), gu + gl) < 1e-12
    check("CL-A2 instrument rate == g_u + g_l closed form (computed comparison)",
          ok, f"worst rel {worst:.1e}")

    # ================================================================ B. THE TWO WIDTHS
    # SEALED (t47_a_width.txt line 39): exact threshold at (1.2 eV, 300 K, 1e9, 10 y).
    d10 = CL.delta_pop(B0, T0, F0, TM0)
    check("CL-B1 delta_pop(1.2 eV, 300 K, 1e9 Hz, 10 y) == sealed 2.535960987846706e-20 J "
          "(sealed anchor)", rel(d10, 2.535960987846706e-20) < 1e-12, f"{d10:.15e} J")
    naive = B0 - kT0 * np.log(F0 * TM0)
    corr = kT0 * np.log1p(np.exp(-d10 / kT0))
    check("CL-B2 exact width == naive width minus derived two-sided correction (computed "
          "identity; sealed naive 0.158338858318 eV beside it)",
          abs((naive - corr) - d10) < 1e-12 * naive and rel(naive / EV, 0.158338858318) < 1e-9,
          f"naive {naive/EV:.12f} eV, correction {corr/kT0:.6f} kT")
    xb = crossing_instr(B0, T0, F0, TM0)
    check("CL-B3 bisected INSTRUMENT crossing lands on delta_pop (definition-not-shortcut)",
          rel(xb, d10) < 1e-11, f"rel {rel(xb, d10):.1e}")
    # SEALED (t47_a_width.txt lines 46, 77, 70): grid anchors, printed to 9 decimals eV.
    g1 = CL.delta_pop(B0, 200.0, 1e9, 1.0)
    g2 = CL.delta_pop(B0, 300.0, 1e11, 3.156e8)
    g3 = CL.delta_pop(B0, 200.0, 1e13, 3.156e7)
    check("CL-B4 grid anchors: (1 s, 1e9, 200 K), (10 y, 1e11, 300 K), (1 y, 1e13, 200 K) "
          "== sealed 0.842841424 / 0.032902903 / 0.386506380 eV (sealed anchors)",
          abs(g1 / EV - 0.842841424) < 1e-9 and abs(g2 / EV - 0.032902903) < 1e-9
          and abs(g3 / EV - 0.386506380) < 1e-9,
          f"{g1/EV:.9f} {g2/EV:.9f} {g3/EV:.9f}")
    nc = CL.delta_pop(B0, 300.0, 1e13, 3.156e7)
    nc_inst = crossing_instr(B0, 300.0, 1e13, 3.156e7)
    check("CL-B5 no-crossing zero WITH positive control (D-15): delta_pop None at "
          "(1 y, 1e13, 300 K) -- sealed no-crossing row -- both routes, beside the 200 K "
          "crossing of CL-B4", nc is None and nc_inst is None and g3 is not None,
          f"None/None beside {g3/EV:.6f} eV")
    # SEALED (t47b_staircase.txt sec A): four widths at (1.2 eV, 350 K, 1e9), 12 decimals eV.
    SEALED_DELTA_350 = {1e2: 0.436077474379, 1e4: 0.297180902789,
                        1e6: 0.158128505564, 1e7: 0.087211155930}
    ok, worst = True, 0.0
    for t_m, sealed in SEALED_DELTA_350.items():
        d = CL.delta_pop(1.2 * EV, 350.0, 1e9, t_m) / EV
        worst = max(worst, abs(d - sealed))
        ok &= abs(d - sealed) < 1e-11
    check("CL-B6 delta_pop at (1.2 eV, 350 K, 1e9), four t_m == sealed 12-decimal values "
          "(t47b sec A, sealed anchors)", ok, f"worst abs {worst:.1e} eV")
    # SEALED (t47_a_width.txt line 110): delta_coh at t_m = 1e-6 s; identity beside it.
    dc = CL.delta_coh(1e-6)
    check("CL-B7 delta_coh = hbar/t_m: sealed 1.054572e-28 J at 1e-6 s AND the exact "
          "identity delta_coh*t_m == hbar", rel(dc, 1.054572e-28) < 1e-6
          and CL.delta_coh(1e-6) * 1e-6 == G.HBAR, f"{dc:.6e} J")
    # SEALED (t47_a_width.txt sec 5): one modulus, two corners; g_u=5, g_l=3, dE/hbar=6.
    w, _ = np.linalg.eig(G.liouvillian(-(6.0 / 2) * SZ,
                                       [np.sqrt(3.0) * SP.conj().T, np.sqrt(5.0) * SP]))
    pop = min(w, key=lambda z: abs(z - (-8.0)))
    coh = max(w, key=lambda z: abs(z.imag))
    m = CL.coh_modulus(8.0, 6.0 * G.HBAR)
    check("CL-B8 coh_modulus == the Liouvillian's own coherence |lambda| (computed) == "
          "sealed 7.211102550928; population face |lambda| == 8 (sealed anchors)",
          rel(abs(coh), m) < 1e-10 and rel(m, 7.211102550928) < 1e-10
          and abs(abs(pop) - 8.0) < 1e-10, f"|coh|={abs(coh):.12f} modulus={m:.12f}")
    # SEALED (t47_a_width.txt sec 4): coherence face -- Hamiltonian-only bisection lands
    # on hbar/t_m (register note: H in generator units on this face); slow-mode counts.
    tm_c = 1.0e-6
    lo, hi = 0.0, 2.0 * G.HBAR / tm_c
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if G.clause_ii(-(mid / 2) * SZ / G.HBAR, [], SP, tm_c)['rate'] <= 1.0 / tm_c:
            lo = mid
        else:
            hi = mid
        if hi - lo <= np.spacing(hi):
            break
    xc = 0.5 * (lo + hi)
    n_lo = len(G.slow_modes(-(0.25 * G.HBAR / tm_c) * SZ / G.HBAR, [], tm_c)[0])
    n_hi = len(G.slow_modes(-(1.0 * G.HBAR / tm_c) * SZ / G.HBAR, [], tm_c)[0])
    check("CL-B9 coherence-face crossing == hbar/t_m (sealed 1.054571817000000e-28 J); "
          "slow-mode counts 4 (durable) / 2 (dropped), sealed",
          rel(xc, G.HBAR / tm_c) < 1e-12 and n_lo == 4 and n_hi == 2,
          f"dE* = {xc:.15e} J, counts {n_lo}/{n_hi}")

    # ================================================================ C. THE CENSUS
    # SEALED (t47_a_width.txt sec 6): the six-record ensemble, its t*_i and its staircase.
    ens6 = [corner("r%d" % i, B0, dE * EV, T0, F0)
            for i, dE in enumerate((0.0, 0.05, 0.10, 0.15, 0.20, 0.25))]
    SEALED_TSTAR = [7.212449e+10, 1.821850e+10, 2.952602e+09,
                    4.344239e+08, 6.296091e+07, 9.104795e+06]
    c1y = CL.census(ens6, 3.156e7)
    ts = [r['t_star'] for r in c1y['schedule']]
    check("CL-C1 six-record drop times == sealed t*_i (6 anchors, rel <= 1e-6)",
          all(rel(a, b) < 1e-6 for a, b in zip(ts, SEALED_TSTAR)),
          f"worst rel {max(rel(a, b) for a, b in zip(ts, SEALED_TSTAR)):.1e}")
    check("CL-C2 t_star (instrument) == t_star_formula on every record (computed, the "
          "closed form is the check never the source)",
          all(rel(r['t_star'], r['t_star_formula']) < 1e-12 for r in c1y['schedule']))
    SEALED_K6 = [(3.156e4, 6), (3.156e6, 6), (3.156e7, 5), (3.156e8, 4),
                 (3.156e9, 2), (3.156e10, 1), (3.156e11, 0)]
    ks = [CL.census(ens6, t)['k'] for t, _ in SEALED_K6]
    check("CL-C3 the sealed staircase [6,6,5,4,2,1,0] over the seven sealed t_m -- the "
          "k = 0 zero gated beside k = 6 (D-15), monotone decreasing",
          ks == [k for _, k in SEALED_K6] and all(a >= b for a, b in zip(ks, ks[1:])),
          f"{ks}")
    check("CL-C4 both census routes agree at all seven sealed t_m (instrument vs width)",
          all(CL.census(ens6, t)['k'] == CL.census(ens6, t)['k_formula']
              for t, _ in SEALED_K6))
    check("CL-C5 drop_order is the dated schedule: first-to-die first == sealed order "
          "r5..r0 (largest dE dies first on the common-B carrier)",
          c1y['drop_order'] == [5, 4, 3, 2, 1, 0], f"{c1y['drop_order']}")

    # SEALED (t47b_staircase.txt sec C): UNIFORM N=12 and LOGNORMAL-LIKE N=10 ensembles.
    Bs_u = np.array([0.95 + 0.40 * i / 11.0 for i in range(12)]) * EV
    dEs_u = np.array([0.18 * ((7 * i + 3) % 12) / 11.0 for i in range(12)]) * EV
    uni = [corner("u%d" % i, Bs_u[i], dEs_u[i], 350.0, 1e9) for i in range(12)]
    SEALED_UNI = [7.845972e+03, 6.995896e+02, 3.315836e+04, 8.896269e+05, 1.302688e+05,
                  5.009554e+06, 4.979281e+05, 2.265904e+07, 1.884988e+06, 9.156432e+07,
                  3.026220e+09, 3.536011e+08]
    zs = np.array([-1.5, -1.1, -0.7, -0.35, 0.0, 0.0, 0.35, 0.7, 1.1, 1.5])
    zp = np.array([0.7, -0.35, 1.5, 0.0, -1.5, 1.1, -0.7, 0.35, -1.1, 0.0])
    logn = [corner("l%d" % i, 1.10 * np.exp(0.12 * zs[i]) * EV,
                   0.05 * np.exp(0.55 * zp[i]) * EV, 350.0, 1e9) for i in range(10)]
    SEALED_LOGN = [1.366468e+03, 1.542332e+04, 8.137589e+03, 2.466660e+05, 2.251496e+06,
                   3.165801e+05, 8.073097e+06, 1.994613e+07, 3.417463e+08, 1.470174e+09]
    tgrid = np.logspace(0, 12, 1201)
    for name, ens, sealed_t, N in (("UNIFORM N=12", uni, SEALED_UNI, 12),
                                   ("LOGNORMAL N=10", logn, SEALED_LOGN, 10)):
        c0 = CL.census(ens, 1.0)
        taus = np.array([r['t_star'] for r in c0['schedule']])
        check(f"CL-C6[{name}] drop times == sealed tau_i ({N} anchors, rel <= 1e-6)",
              all(rel(a, b) < 1e-6 for a, b in zip(taus, sealed_t)),
              f"worst rel {max(rel(a, b) for a, b in zip(taus, sealed_t)):.1e}")
        k1 = (taus[None, :] >= tgrid[:, None]).sum(axis=1)
        k2 = []
        for t in tgrid:
            n = 0
            for r in c0['schedule']:
                d = CL.delta_pop(r['B'], r['T'], r['f0'], t)
                n += (d is not None and r['dE'] <= d)
            k2.append(n)
        k2 = np.array(k2)
        check(f"CL-C7[{name}] two routes agree at ALL 1201 grid t_m (sealed: max "
              f"discrepancy 0); k(1 s) = {N} beside k(1e12 s) = 0 (D-15)",
              int(np.abs(k1 - k2).max()) == 0 and k1[0] == N and k1[-1] == 0
              and bool(np.all(np.diff(k1) <= 0)),
              f"max |k1-k2| = {int(np.abs(k1-k2).max())}")

    # SEALED (t47b_staircase.txt sec D): the symmetric D-15 control, factor-2 corner.
    sym = [corner("s%d" % i, 1.0 * EV, 0.0, 350.0, 1e9) for i in range(8)]
    csym = CL.census(sym, 1.0)
    taus_s = [r['t_star'] for r in csym['schedule']]
    bound = np.exp(1.0 * EV / (G.KB * 350.0)) / (2 * 1e9)
    ks_s = sorted({CL.census(sym, t)['k'] for t in np.logspace(0, 8, 801)}, reverse=True)
    check("CL-C8 symmetric control (D-15): all 8 tau equal sealed 1.254112e+05 s == "
          "e^(B/kT)/(2 f0) computed (the derived factor-2 corner), staircase flat {8, 0}",
          all(rel(t, 1.254112e+05) < 1e-6 for t in taus_s)
          and rel(taus_s[0], bound) < 1e-12 and ks_s == [8, 0],
          f"tau {taus_s[0]:.6e}, bound {bound:.6e}, values {ks_s}")

    # SEALED (t47b_staircase.txt sec G): the departure term on the 10-grain carrier.
    kTg = G.KB * 350.0
    xs = [0.2, 0.5, 0.8, 1.2, 1.8, 2.5, 3.2, 4.0, 5.0, 6.5]
    grains = [corner("g%d" % i, 1.1 * EV, x * kTg, 350.0, 1e9) for i, x in enumerate(xs)]
    cg = CL.census(grains, 1e6)
    dead = [r['index'] for r in cg['schedule'] if not r['alive']]
    check("CL-C9 departure carrier (sealed sec G): k = 4 at t = 1e6 s, dead grains "
          "exactly [4..9], tau_0/tau_9 == sealed 3.1094e+06 / 1.0369e+04 s",
          cg['k'] == 4 and dead == [4, 5, 6, 7, 8, 9]
          and rel(cg['schedule'][0]['t_star'], 3.1094e+06) < 1e-4
          and rel(cg['schedule'][9]['t_star'], 1.0369e+04) < 1e-4,
          f"k={cg['k']} dead={dead}")
    csym_dead = CL.census(sym, 1e8)          # symmetric carrier past its threshold
    check("CL-C10 departure == sealed 5.4339 (sum over dead of tanh(dE/2kT)) WITH the "
          "D-15 zero beside it: symmetric carrier, all dead, departure exactly 0",
          abs(cg['departure'] - 5.4339) < 5e-4
          and csym_dead['k'] == 0 and csym_dead['departure'] == 0.0,
          f"departure {cg['departure']:.4f}; symmetric k={csym_dead['k']} "
          f"departure={csym_dead['departure']}")

    # ================================================================ D. THE T-31 CONTROL
    # SEALED (t47b_staircase.txt sec F): the kill becomes a staircase.
    b0 = CL.t31_basin(0.0)
    check("CL-D1 T-31 corner (eps=0): basin taus both == sealed 5.5066e+03; exact-"
          "multiplicity corner reading v2 = 2 == the record count 2 (the corner law is "
          "model/count_law.py's, referenced not duplicated)",
          all(rel(t, 5.5066e+03) < 5e-4 for t in b0['taus'])
          and b0['v2_exact'] == 2 and b0['multiplicities'] == [4, 8, 4],
          f"taus {b0['taus'][0]:.4e}/{b0['taus'][1]:.4e}, mults {b0['multiplicities']}")
    SEALED_T31 = {0.05: (0.0698, 0.0641, 5.5008e+03, 5.5147e+03),
                  0.10: (0.2550, 0.2325, 5.2452e+03, 5.4219e+03),
                  0.16: (0.5684, 0.5124, 3.6284e+03, 4.1648e+03)}
    ok = True
    detail = []
    for eps, (d1, d2, t1, t2) in SEALED_T31.items():
        b = CL.t31_basin(eps)
        ok &= (abs(b['dE'][0] - d1) < 5e-4 and abs(b['dE'][1] - d2) < 5e-4
               and rel(b['taus'][0], t1) < 5e-4 and rel(b['taus'][1], t2) < 5e-4)
        detail.append((eps, round(b['dE'][0], 4), round(b['taus'][0], 1)))
    check("CL-D2 T-31 asymmetric grid: dE_i read off the carrier's own spectrum and "
          "basin taus == sealed rows at eps = 0.05, 0.10, 0.16 (sealed anchors)",
          ok, f"{detail}")
    gen_worst = max(CL.t31_basin(e)['gen_dev'] / max(CL.t31_basin(e)['slow_rates'][-1], 1e-300)
                    for e in (0.0, 0.05, 0.10, 0.16))
    check("CL-D3 basin Q is a CLASSICAL MARKOV GENERATOR at every eps (sealed worst "
          "3.73e-11; gate <= 1e-6 of the slow scale) -- measured, not assumed",
          gen_worst <= 1e-6, f"worst {gen_worst:.2e}")
    b16 = CL.t31_basin(0.16)
    tlo, thi = min(b16['taus']), max(b16['taus'])
    grid16 = np.logspace(np.log10(tlo) - 2, np.log10(thi) + 2, 401)
    vals16 = sorted(set(CL.t31_staircase(0.16, grid16)), reverse=True)
    check("CL-D4 the staircase at eps = 0.16: values [2, 1, 0] over the sealed grid; "
          "step order tracks the derived exponent (smaller B - dE dies first)",
          vals16 == [2, 1, 0]
          and ((b16['B_minus_dE'][0] < b16['B_minus_dE'][1])
               == (b16['taus'][0] < b16['taus'][1])),
          f"values {vals16}, taus {b16['taus'][0]:.4e}/{b16['taus'][1]:.4e}")
    check("CL-D5 THE C-76 KILL AS CONTROL (D-15 pair): exact-multiplicity v2 reads 0 at "
          "EVERY eps > 0 (T-31's kill, all multiplicities 1) WHILE the record-mode census "
          "reads 2 at small t_m (what survived it)",
          all(CL.t31_basin(e)['v2_exact'] == 0 for e in (0.05, 0.10, 0.16))
          and CL.t31_staircase(0.16, [grid16[0]]) == [2],
          f"v2 {[CL.t31_basin(e)['v2_exact'] for e in (0.05, 0.10, 0.16)]}, "
          f"k(small t) = {CL.t31_staircase(0.16, [grid16[0]])[0]}")
    # SEALED (t31_asymmetry.txt): splitting/delta ratio class 2.99e36..2.99e41 at 10 y --
    # the derived coherence width can never recover population clustering; the kill's own
    # arithmetic, anchored.
    d10y = CL.delta_coh(10.0 * YEAR)
    ratios = [split / d10y for split in (1e-6, 1e-3, 1e-1)]
    check("CL-D6 sealed T-31 ratio class: splitting/delta_coh(10 y) == 2.99e36/39/41 for "
          "the declared 1e-6/1e-3/1e-1 J splittings; delta == sealed 3.341e-43 J",
          all(rel(r, s) < 5e-3 for r, s in zip(ratios, (2.99e36, 2.99e39, 2.99e41)))
          and rel(d10y, 3.341e-43) < 5e-4,
          f"delta {d10y:.4e} J, ratios {['%.2e' % r for r in ratios]}")

    # ================================================================ E. THE C-76 GATE
    check("CL-E1 census signature is EXACTLY (surfaces, t_m) -- the C-76 lesson in the "
          "API's shape", list(inspect.signature(CL.census).parameters) == ['surfaces', 't_m'],
          f"{list(inspect.signature(CL.census).parameters)}")
    blocked = []
    for kw in ('width', 'tol', 'delta', 'cluster_width', 'margin'):
        try:
            CL.census(ens6, 1.0, **{kw: 1e-2})
            blocked.append((kw, False))
        except TypeError:
            blocked.append((kw, True))
    ctl = CL.census(ens6, 1.0)
    check("CL-E2 the chosen-width form is UNREACHABLE: census(..., width/tol/delta/"
          "cluster_width/margin=...) all TypeError, beside the same call without -> k "
          "returned (D-15 pair)",
          all(b for _, b in blocked) and ctl['k'] == 6,
          f"blocked {blocked}, control k = {ctl['k']}")
    offenders = []
    for nm in dir(CL):
        if nm.startswith('_'):
            continue
        f = getattr(CL, nm)
        if callable(f) and getattr(f, '__module__', None) == 'countlaw':
            for p in inspect.signature(f).parameters:
                if any(bad in p.lower() for bad in ('width', 'tol', 'cluster', 'margin')):
                    offenders.append((nm, p))
    check("CL-E3 no public countlaw callable accepts a width/tolerance/cluster/margin "
          "parameter (module-wide scan)", offenders == [], f"offenders {offenders}")

    # ================================================================ F. OBSERVATION ENTRY
    # How a real device census enters: URM.surface (D-25 provenance) -> census. The NAND
    # constants are check-side declared values of the registry's provenance class.
    nand = URM.surface("NAND floating gate", "trapped charge",
                       0.05 * EV, 1.0 * EV, 358.0, 1e9)
    cw_short = CL.census([nand], 1.0e3)
    cw_long = CL.census([nand], 1.0 * YEAR)
    check("CL-F1 world-tier entry: registry provenance carried into the census row; "
          "alive at 1e3 s (positive) beside dead at 1 y (zero, D-15); both routes agree "
          "both times",
          bool(cw_short['schedule'][0]['provenance'])
          and 'Weller' not in str(cw_short['schedule'][0]['provenance'])
          and cw_short['k'] == 1 and cw_long['k'] == 0
          and cw_short['k'] == cw_short['k_formula']
          and cw_long['k'] == cw_long['k_formula'],
          f"k(1e3 s)={cw_short['k']} k(1 y)={cw_long['k']}, provenance="
          f"{str(cw_short['schedule'][0]['provenance'])[:40]}...")
    refused = False
    try:
        URM.surface("unregistered mystery device", "unknown", 0.05 * EV, 1.0 * EV, 300.0, 1e9)
    except ValueError:
        refused = True
    check("CL-F2 D-25 refusal at the gate: an unregistered world surface without "
          "provenance raises, beside the accepted registered one (D-15 pair)",
          refused and cw_short['k'] == 1)
    nonth = URM.surface("NAND floating gate", "trapped charge",
                        0.05 * EV, 1.0 * EV, 358.0, 1e9, thermal=False)
    cmix = CL.census([nand, nonth], 1.0e3)
    check("CL-F3 decline is DECLARED, never silent: thermal=False surface lands in "
          "`declined` (len 1), k unchanged beside the thermal control that counts",
          len(cmix['declined']) == 1 and cmix['k'] == 1
          and len(cmix['schedule']) == 1 and cmix['declined'][0]['index'] == 1,
          f"declined {cmix['declined'][0]['name']}, k = {cmix['k']}")

    # ================================================================ G. BEYOND THE GATED RANGE
    t3 = c1y['schedule'][3]['t_star']
    kA = CL.census(ens6, 0.999 * t3)['k']
    kB = CL.census(ens6, 1.001 * t3)['k']
    check("CL-G1 probe BEYOND sealed grid: the step sits exactly at t*_3 -- k drops "
          "4 -> 3 across (0.999, 1.001) x t*_3 (the step is where the law says, not "
          "hard-coded)", kA == 4 and kB == 3, f"k = {kA} -> {kB} across t*_3 = {t3:.4e} s")
    # fresh N=5 ensemble, MIXED per-record f0 (no sealed lane ever mixed f0), T = 320 K:
    fresh = [corner("f%d" % i, (1.02 + 0.03 * i) * EV, 0.017 * ((3 * i + 1) % 5) * EV,
                    320.0, 10.0 ** (9 + i % 3)) for i in range(5)]
    cf = CL.census(fresh, 1.0)
    ok_form = all(rel(r['t_star'], r['t_star_formula']) < 1e-10 for r in cf['schedule'])
    gridf = np.logspace(0, 13, 301)
    kf1 = [CL.census(fresh, t)['k'] for t in gridf]
    kf2 = [CL.census(fresh, t)['k_formula'] for t in gridf]
    check("CL-G2 probe BEYOND every sealed case: fresh declared N=5 ensemble with MIXED "
          "per-record f0 -- routes agree at all 301 t_m, drop times land on the closed "
          "form <= 1e-10, staircase monotone",
          kf1 == kf2 and ok_form and all(a >= b for a, b in zip(kf1, kf1[1:]))
          and kf1[0] == 5 and kf1[-1] == 0, f"k(1 s) = {kf1[0]}, k(1e13 s) = {kf1[-1]}")
    b20 = CL.t31_basin(0.20)
    grid20 = np.logspace(np.log10(min(b20['taus'])) - 2, np.log10(max(b20['taus'])) + 2, 401)
    vals20 = sorted(set(CL.t31_staircase(0.20, grid20)), reverse=True)
    check("CL-G3 T-31 probe BEYOND the sealed eps grid (eps = 0.20): Q still classical, "
          "staircase still [2, 1, 0], order still tracks B - dE",
          b20['gen_dev'] / max(b20['slow_rates'][-1], 1e-300) <= 1e-6
          and vals20 == [2, 1, 0]
          and ((b20['B_minus_dE'][0] < b20['B_minus_dE'][1])
               == (b20['taus'][0] < b20['taus'][1])),
          f"values {vals20}, taus {b20['taus'][0]:.4e}/{b20['taus'][1]:.4e}")
    x500 = crossing_instr(1.2 * EV, 500.0, 1e9, 1.0)
    d500 = CL.delta_pop(1.2 * EV, 500.0, 1e9, 1.0)
    check("CL-G4 width probe BEYOND the sealed 200-400 K grid: bisected instrument "
          "crossing at 500 K lands on delta_pop (rel <= 1e-11)",
          rel(x500, d500) < 1e-11, f"rel {rel(x500, d500):.1e} at delta {d500/EV:.6f} eV")


# =====================================================================================
if __name__ == '__main__':
    import time
    t0 = time.time()
    counts = {'pass': 0, 'fail': 0}

    def check(name, cond, detail=""):
        counts['pass' if cond else 'fail'] += 1
        print(f"  {'PASS' if cond else 'FAIL'}  {name}  {detail}")

    print("VALIDATE THE COUNT-LAW LAYER (C-86) THROUGH model/countlaw.py")
    print("=" * 78)
    run_countlaw_checks(check)
    print("=" * 78)
    print(f"  COUNTLAW: {counts['pass']} PASS, {counts['fail']} FAIL"
          f"   ({time.time() - t0:.1f} s)")
    sys.exit(0 if counts['fail'] == 0 else 1)
