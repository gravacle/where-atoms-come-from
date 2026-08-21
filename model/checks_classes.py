"""EVERY SEALED REACHABLE-CLASS HEADLINE (C-87, C-90), REPRODUCED THROUGH model/classes.py.

Check-block for the classes family (T-54 fold-in), in the validate_geometry.py idiom:
run_classes_checks(check) issues check(name, cond, detail) calls; the registrar chains it
into the validator stack.  Standalone: `python3 checks_classes.py` runs the block with its
own counter and exits 0 iff every check passes.

Sealed sources gated here:
  C-87   LANE_T44_A_CORNER/t44a_corner.OUT.txt        (corner venue row sums; resolvent
                                                       singular at 1/4, solvable at 1/8 and
                                                       23/100; chain mu_c = 1/2; sandwich)
  C-87/  LANE_T44_B_WORLD/t44b_world.OUT.txt +        (world venue; mu_c = 1/6 three ways;
  C-90   t44b_world.RESULT.json                        subcritical rows; witnesses; the
                                                       M=1400 kernel brackets; G(0) bracket;
                                                       discriminator windows; coefficient at
                                                       M=1400)
  C-90   register row C-90 / LANE_T44_B_WORLD/        (the M=2800 deepened coefficient
         VERIFY/adv_verify.OUT.txt E7                  bracket [0.476369, 0.487321] --
                                                       d*G(d) converging onto 3/(2 pi))

Where a sealed number IS the anchor it is gated as the sealed anchor it is, stated as such
(D-8); every zero is gated with a positive control beside it (D-15); owner values appear
only in checks labeled COMPARISON.  API-fidelity probes BEYOND the gated range: the
resolvent identity at the unswept mu = 1/9; the untested kernel target (10,0,0); the
untested venue n = 3 entered through the provenance gate end-to-end; kernel_pass ==
crit_kernel_3d at a spot depth (one instrument, two drivers)."""
import os as _os
import sys as _sys
from fractions import Fraction
from math import comb, isqrt

_HERE = _os.path.dirname(_os.path.abspath(__file__))
if _HERE not in _sys.path:
    _sys.path.insert(0, _HERE)

import classes as CL
from classes import Fraction, ff  # noqa: F811


def run_classes_checks(check):
    # ================================================================ A: venue entry (D-25)
    w7, w7_cells, w7_idx = CL.world_venue(7)
    check("C-87 world venue n=7: row sums all exactly 6 == sealed (deg-6 venue, Perron "
          "landmark input)", all(r == 6 for r in w7.row_sums()),
          f"deg={w7.degree()}")
    src7 = w7_idx[(0, 0, 0)]
    dist7 = CL.bfs_venue(w7, src7)
    check("C-87 world venue n=7: BFS distance == earned separation (L1 with wraps) on "
          "every grain == sealed",
          all(dist7[i] == sum(min(abs(c), 7 - abs(c)) for c in w7_cells[i])
              for i in range(w7.n)))
    c46, _c46_cells, _c46_idx = CL.corner_venue(4, 6)
    c37, _, _ = CL.corner_venue(3, 7)
    check("C-87 corner venue from the model's own carrier: (4,6) and (3,7) row sums all "
          "exactly 4 == sealed (dual lattice of geometry.Torus plaquette supports)",
          all(r == 4 for r in c46.row_sums()) and all(r == 4 for r in c37.row_sums()),
          f"deg(4,6)={c46.degree()} deg(3,7)={c37.degree()}")
    ch12 = CL.chain_venue(12)
    check("C-87 chain venue L=12: row sums all exactly 2 == sealed (D=1 control)",
          all(r == 2 for r in ch12.row_sums()))
    refused = False
    try:
        CL.venue("anonymous geometry", ch12.adj)
    except ValueError:
        refused = True
    accepted = CL.venue("named geometry", ch12.adj,
                        provenance="declared test source (this check)")
    corner_refused = False
    try:
        CL.venue("idealisation", ch12.adj, tier="corner")
    except ValueError:
        corner_refused = True
    check("D-25 gate: venue entry REFUSES a world-tier graph without provenance and a "
          "corner graph without the DEF-A self-declaration; POSITIVE CONTROL beside: with "
          "provenance the same adjacency enters and carries it",
          refused and corner_refused and
          accepted.provenance == "declared test source (this check)" and
          c46.provenance == "DEF-A" and CL.WORLD_VENUE_PROVENANCE in (w7.provenance,))

    # ================================================================ B: mu_c three ways
    w4, _w4_cells, _w4_idx = CL.world_venue(4)
    loc4 = CL.mu_c_of(w4, certify="full")
    check("C-90 mu_c of the world venue (n=4) LOCATED == sealed 1/6 EXACT: Perron row-sum "
          "candidate 1/deg, (I - A/6) annihilates the constant vector (pole AT 1/6)",
          loc4.get("located") and loc4["mu_c"] == Fraction(1, 6)
          and loc4["pole_at_candidate"],
          f"mu_c={loc4.get('mu_c')} deg={loc4.get('degree')}")
    check("C-90 D-15 zero + control: exact resolvent SINGULAR at the computed mu_c (the "
          "zero) beside SOLVABLE at 19/120 and 7/40 (the positive controls -- the sealed "
          "probe points, = mu_c*(19/20) and mu_c*(21/20))",
          loc4.get("resolvent_singular_at_mu_c") and loc4.get("solvable_below")
          and loc4.get("solvable_above")
          and loc4["mu_c"] * Fraction(19, 20) == Fraction(19, 120)
          and loc4["mu_c"] * Fraction(21, 20) == Fraction(7, 40))
    sand_ok = all(2 * 36 ** m <= CL.n3_even_row(m, 0, 0, 0) * (2 * m + 1) * (m + 1) * (m + 2)
                  and CL.n3_even_row(m, 0, 0, 0) <= 36 ** m for m in range(1, 41))
    binom_ok = all(comb(2 * m, m) * (2 * m + 1) >= 4 ** m for m in range(1, 301))
    check("C-90 mu_c third way == sealed sector sandwich: 2*36^m <= N_2m(0)(2m+1)(m+1)(m+2) "
          "and N_2m(0) <= 36^m (m <= 40); C(2m,m)(2m+1) >= 4^m (m <= 300)",
          sand_ok and binom_ok)
    loc7 = CL.mu_c_of(w7, certify="rowsum")
    check("C-90 world venue n=7 (row-sum tier): its own mu_c == 1/6 == sealed (mu_c is "
          "the venue's own computed number)",
          loc7.get("located") and loc7["mu_c"] == Fraction(1, 6))
    locc = CL.mu_c_of(c46, certify="full")
    check("C-87 corner venue (4,6): its own mu_c LOCATED == sealed 1/4 EXACT (resolvent "
          "singular at 1/4, solvable beside)",
          locc.get("located") and locc["mu_c"] == Fraction(1, 4),
          f"mu_c={locc.get('mu_c')}")
    check("C-87 corner venue (4,6): resolvent nonsingular at the SEALED probe points 1/8 "
          "and 23/100 (t44a S2), singular at the computed 1/4",
          CL.resolvent_exact(c46.adj, Fraction(1, 8), 0) is not None and
          CL.resolvent_exact(c46.adj, Fraction(23, 100), 0) is not None and
          CL.resolvent_exact(c46.adj, Fraction(1, 4), 0) is None)
    locch = CL.mu_c_of(ch12, certify="full")
    check("C-87 chain venue: its own mu_c LOCATED == sealed 1/2 EXACT (the critical point "
          "is the venue's own, dimension-dependent, computed)",
          locch.get("located") and locch["mu_c"] == Fraction(1, 2),
          f"mu_c={locch.get('mu_c')}")

    # ================================================================ C: the verdict triple
    v_sub = CL.class_verdict(w4, Fraction(1, 8))
    v_cri = CL.class_verdict(w4, Fraction(1, 6))
    v_sup = CL.class_verdict(w4, Fraction(13, 72))
    check("C-87 class verdict triple on the world venue: mu=1/8 -> (EXPONENTIAL,-,-), "
          "mu=1/6 -> (-,CRITICAL,-), mu=13/72 -> (-,-,DIVERGENT); each exactly one True; "
          "each verdict's mu_c is the venue's own computed 1/6 == sealed taxonomy",
          (v_sub["exponential"], v_sub["critical"], v_sub["divergent"]) == (True, False, False)
          and (v_cri["exponential"], v_cri["critical"], v_cri["divergent"]) == (False, True, False)
          and (v_sup["exponential"], v_sup["critical"], v_sup["divergent"]) == (False, False, True)
          and v_sub["mu_c"] == v_cri["mu_c"] == v_sup["mu_c"] == Fraction(1, 6))
    v_c2 = CL.class_verdict(c46, Fraction(1, 6))
    v_c1 = CL.class_verdict(ch12, Fraction(1, 6))
    check("C-87 mu = 1/6 is CRITICAL on the D=3 venue and SUBCRITICAL on the D=2 and D=1 "
          "venues == sealed S5 (mu_c is the venue's own computed number)",
          v_cri["critical"] and v_c2["exponential"] and not v_c2["critical"]
          and v_c1["exponential"],
          f"D2 mu_c={v_c2['mu_c']} D1 mu_c={v_c1['mu_c']}")
    sweep = [Fraction(1, 12), Fraction(1, 8), Fraction(3, 20), Fraction(4, 25),
             Fraction(1, 6), Fraction(13, 72), Fraction(1, 5), Fraction(1, 4)]
    check("C-87 the triple is exhaustive and exclusive over the full sealed sweep "
          "(exactly one class boolean True at every declared mu)",
          all(sum([CL.class_verdict(w4, mu)[k] for k in
                   ("exponential", "critical", "divergent")]) == 1 for mu in sweep))

    # ================================================================ D: subcritical rows
    r12 = CL.subcritical_row(Fraction(1, 12), K=90, dmax=12)
    rr = r12["ratios"][11]
    check("C-87 mu=1/12 row == sealed: EXPONENTIAL by computed booleans; r(dmax) == "
          "[0.11478359, 0.11478359]; d(1-r) == 9.737; G(1) == 0.0935147685, "
          "G(4) == 0.0000807666 (RESULT.json, 8/10dp)",
          r12["cls"] == "EXPONENTIAL"
          and (ff(rr[0], 8), ff(rr[1], 8)) == ("0.11478359", "0.11478359")
          and ff(r12["power_exclusion"], 3) == "9.737"
          and ff(r12["G"][1][0], 10) == "0.0935147685"
          and ff(r12["G"][4][0], 10) == "0.0000807666",
          f"r=[{ff(rr[0], 8)},{ff(rr[1], 8)}] qpow={ff(r12['power_exclusion'], 3)}")
    check("C-87 mu=1/12 COMPARISON (owner: Ornstein-Zernike rate, labeled comparison "
          "only): within 1/25; owner lo == sealed 0.12701665",
          r12["comparison_oz_within_tol"]
          and ff(r12["comparison_owner_rate"][0], 8) == "0.12701665")
    r18 = CL.subcritical_row(Fraction(1, 8), K=160, dmax=14)
    rr = r18["ratios"][13]
    check("C-87 mu=1/8 row == sealed: EXPONENTIAL by computed booleans; r(dmax) lo == "
          "0.24792841; d(1-r) == 9.776; G(1,2,4,8,14) == 0.1699358732 / 0.0280775518 / "
          "0.0009845743 / 0.0000022863 / 0.0000000004 (RESULT.json, 10dp)",
          r18["cls"] == "EXPONENTIAL" and ff(rr[0], 8) == "0.24792841"
          and ff(r18["power_exclusion"], 3) == "9.776"
          and ff(r18["G"][1][0], 10) == "0.1699358732"
          and ff(r18["G"][2][0], 10) == "0.0280775518"
          and ff(r18["G"][4][0], 10) == "0.0009845743"
          and ff(r18["G"][8][0], 10) == "0.0000022863"
          and ff(r18["G"][14][0], 10) == "0.0000000004"
          and ff(r18["comparison_owner_rate"][0], 8) == "0.26794919"
          and r18["comparison_oz_within_tol"],
          f"r_lo={ff(rr[0], 8)} qpow={ff(r18['power_exclusion'], 3)}")
    mu, K = Fraction(1, 8), 160
    Sx, _ = CL.series_3d(mu, (2, 0, 0), K)
    Sn1, _ = CL.series_3d(mu, (1, 0, 0), K - 1)
    Sn3, _ = CL.series_3d(mu, (3, 0, 0), K - 1)
    Sn2, _ = CL.series_3d(mu, (2, 1, 0), K - 1)
    check("C-87 resolvent identity S_K(2,0,0) = mu [S_{K-1}(1,0,0) + S_{K-1}(3,0,0) + "
          "4 S_{K-1}(2,1,0)] EXACT at mu=1/8 == sealed S3 (the measured G is the venue "
          "equation's own solution)", Sx == mu * (Sn1 + Sn3 + 4 * Sn2))
    mu, K = Fraction(1, 9), 120
    Sx, _ = CL.series_3d(mu, (2, 0, 0), K)
    Sn1, _ = CL.series_3d(mu, (1, 0, 0), K - 1)
    Sn3, _ = CL.series_3d(mu, (3, 0, 0), K - 1)
    Sn2, _ = CL.series_3d(mu, (2, 1, 0), K - 1)
    check("API-FIDELITY PROBE (beyond the sealed sweep): the same identity EXACT at the "
          "UNSWEPT mu = 1/9 -- the series is the venue's own walk sum by definition, not "
          "a table of sealed rows", Sx == mu * (Sn1 + Sn3 + 4 * Sn2))

    # ================================================================ E: supercritical
    sealed_super = {Fraction(13, 72): ("3398653.92", "1.1619"),
                    Fraction(1, 5): ("71952995756302043835.71", "1.4257"),
                    Fraction(1, 4):
                    ("8512414893074179630698564026681656120355239071322.37", "2.2277")}
    sup_ok = True
    detail = []
    for mu, (st, sr) in sealed_super.items():
        wit = CL.divergence_witness(mu)
        sup_ok &= (wit["grows"] and ff(wit["term"], 2) == st
                   and ff(wit["term_ratio"], 4) == sr)
        detail.append((str(mu), ff(wit["term_ratio"], 4)))
    check("C-87 divergence witnesses at every sealed supercritical mu (13/72, 1/5, 1/4): "
          "terms GROW (computed booleans); witness term and term ratio == sealed strings",
          sup_ok, f"ratios={detail}")

    # ================================================================ F: leading stratum
    DP = CL.dp3_counts(10)
    check("C-87 axis-split reference N3_ref == brute-force DP, k <= 10, all |x|_inf <= 4 "
          "== sealed gate (trimmed range; the identity is the measurement route's ground)",
          all(CL.N3_ref(k, a, b, c) == DP.get((k, a, b, c), 0)
              for k in range(11) for a in range(-4, 5) for b in range(-4, 5)
              for c in range(-4, 5)))
    probe_targets = [(0, 0, 0), (2, 0, 0), (4, 0, 0), (8, 0, 0), (10, 0, 0), (16, 0, 0)]
    check("C-90 fast even-row == reference on every target in use here PLUS the probe "
          "target (10,0,0), m <= 5 == sealed gate idiom",
          all(CL.n3_even_row(m, *t) == CL.N3_ref(2 * m, *t)
              for m in range(6) for t in probe_targets))
    LEAD = [(1, 0, 0), (2, 0, 0), (3, 0, 0), (1, 1, 0), (2, 1, 0), (1, 1, 1), (2, 2, 2)]
    lead_ok = True
    for t in LEAD:
        d = sum(abs(c) for c in t)
        nmin = comb(d, abs(t[0])) * comb(d - abs(t[0]), abs(t[1]))
        below = all(DP.get((k,) + t, 0) == 0 for k in range(d))          # the D-15 zero
        atd = DP.get((d,) + t, 0) == nmin and nmin > 0                    # positive control
        S, _tail = CL.series_3d(Fraction(1, 1000), t, d + 6)
        dev = abs(S / Fraction(1, 1000) ** d - nmin)
        bnd = Fraction(6) ** d * (6 * Fraction(1, 1000)) / (1 - 6 * Fraction(1, 1000))
        lead_ok &= below and atd and dev <= bnd
    check("C-87 leading stratum == sealed: N_k(x) = 0 below d (D-15 zeros) WITH the "
          "positive control N_d(x) = d!/(a!b!c!) > 0 beside each; G/mu^d -> N_min within "
          "the proven bound at mu=1/1000 (w_min = d, the Gamma price; C-80/O-54 standing)",
          lead_ok)
    check("C-87 parity zeros (D-15): N_k((2,0,0)) == 0 at every odd k <= 9, positive "
          "even-k counts beside them",
          all(DP.get((k, 2, 0, 0), 0) == 0 for k in (1, 3, 5, 7, 9)) and
          all(DP.get((k, 2, 0, 0), 0) > 0 for k in (2, 4, 6, 8, 10)))
    check("C-87 octahedral symmetry of the venue counts, k <= 8 == sealed gate (trimmed)",
          all(DP.get((k, a, b, c), 0) == DP.get((k, abs(c), abs(a), abs(b)), 0)
              for k in range(9) for a in range(-3, 4) for b in range(-3, 4)
              for c in range(-3, 4)))

    # ================================================================ G: the critical kernel
    tc = CL.tail_constants()
    B5 = tc["B5"]
    check("C-90 tail constants == sealed displays: Q3 == 0.826995, B5 == 0.466824; "
          "p_2m(0)(m-2)^{3/2} <= B5 EXACT at m in [250, 600] step 25 == sealed S1",
          ff(tc["Q3"]) == "0.826995" and ff(B5) == "0.466824" and
          all(Fraction(CL.n3_even_row(m, 0, 0, 0), 36 ** m) * (m - 2) * (isqrt(m - 2) + 1)
              <= B5 for m in range(250, 601, 25)),
          f"Q3={ff(tc['Q3'])} B5={ff(B5)}")
    A2, A4, A8, A16 = (2, 0, 0), (4, 0, 0), (8, 0, 0), (16, 0, 0)
    ker60a = CL.crit_kernel_3d([A2, A4], 60)
    ker60b = CL.kernel_pass({A2: 60, A4: 60}, 60)
    check("API-FIDELITY: kernel_pass == crit_kernel_3d EXACTLY at spot depth M=60 (one "
          "instrument, two drivers -- the combined pass is no shortcut)",
          ker60b["ker"] == ker60a[0] and ker60b["s0"][60] == ker60a[1]
          and ker60b["p2m0"][60] == ker60a[2])
    KP = CL.kernel_pass({A2: 1400, A4: 1400, A8: 2800, A16: 2800}, 2800,
                        gate_range=(30, 240), gate_bound=CL.gate_bound_L3,
                        snapshots=(350, 700, 1400))
    check("C-90 assembled L3-3D difference bound holds EXACTLY, m in [30, 240], all four "
          "axis targets == sealed S4 gate", KP["asm_ok"])
    S0_1400 = KP["s0"][1400]
    TABS = CL.abs_tail_bound(1400, B5, KP["p2m0"][1400])
    G0 = (S0_1400, S0_1400 + TABS)
    inc1 = KP["s0"][700] - KP["s0"][350]
    inc2 = KP["s0"][1400] - KP["s0"][700]
    check("C-90 G(0) at mu_c CONVERGES == sealed: bracket [1.503919, 1.554391] (6dp), "
          "doubling increments SHRINK 0.007287 -> 0.005159, tail < 1/15",
          (ff(G0[0]), ff(G0[1])) == ("1.503919", "1.554391")
          and ff(inc1) == "0.007287" and ff(inc2) == "0.005159"
          and inc2 < inc1 and inc2 < Fraction(1, 50) and TABS < Fraction(1, 15),
          f"G0=[{ff(G0[0])},{ff(G0[1])}] inc {ff(inc1)}->{ff(inc2)}")
    check("C-90 COMPARISON (owner: Watson 1939, stated as the sealed anchor it is): "
          "G(0) = 1.5163860591... lies INSIDE the computed bracket",
          G0[0] <= CL.WATSON[0] and CL.WATSON[1] <= G0[1])
    KER1 = KP["ker_at"][1400]
    TL1 = {t: CL.diff_tail_bound(t, 1400, B5) for t in (A2, A4, A8, A16)}
    r24 = CL.doubling_ratio(KER1, TL1, (A2, A4), (A4, A8))
    r48 = CL.doubling_ratio(KER1, TL1, (A4, A8), (A8, A16))
    check("C-90 axis doubling pair 2->4 at M=1400 == sealed [0.453369, 0.458705], in the "
          "INV window [2/5, 3/5], outside LOG and LIN",
          (ff(r24[0]), ff(r24[1])) == ("0.453369", "0.458705")
          and CL.INV_LO <= r24[0] and r24[1] <= CL.INV_HI and r24[1] < CL.LOG_LO,
          f"[{ff(r24[0])},{ff(r24[1])}]")
    check("C-90 KEY == sealed: the critical class is a POWER LAW with EXPONENT 1 -- "
          "deepest axis pair 4->8 == [0.469869, 0.507905], bracket contains 1/2 within "
          "the INV window (2^-p, p=1: the exponent bracket contains 1)",
          (ff(r48[0]), ff(r48[1])) == ("0.469869", "0.507905")
          and r48[0] <= Fraction(1, 2) <= r48[1]
          and CL.INV_LO <= r48[0] and r48[1] <= CL.INV_HI,
          f"[{ff(r48[0])},{ff(r48[1])}]")
    H816_1 = CL.increment_interval(KER1, TL1, A8, A16)
    c_ax1 = (16 * H816_1[0], 16 * H816_1[1])
    check("C-90 coefficient 2*d*H at axis d=8, M=1400 == sealed [0.467230, 0.500095]; "
          "COMPARISON: owner 3/(2 pi) = 0.477465 within 1/25 (after the class)",
          (ff(c_ax1[0]), ff(c_ax1[1])) == ("0.467230", "0.500095")
          and c_ax1[0] - CL.COMP_TOL <= CL.C_3D[0] and CL.C_3D[1] <= c_ax1[1] + CL.COMP_TOL,
          f"[{ff(c_ax1[0])},{ff(c_ax1[1])}]")
    KER2 = KP["ker"]
    TL2 = {t: CL.diff_tail_bound(t, 2800, B5) for t in (A8, A16)}
    H816_2 = CL.increment_interval(KER2, TL2, A8, A16)
    c_ax2 = (16 * H816_2[0], 16 * H816_2[1])
    check("C-90 THE DEEPENED COEFFICIENT BRACKET (M=2800) == register row C-90: d*G(d) "
          "converging onto [0.476369, 0.487321]; COMPARISON: owner 3/(2 pi) INSIDE; "
          "strictly tighter than the M=1400 bracket",
          (ff(c_ax2[0]), ff(c_ax2[1])) == ("0.476369", "0.487321")
          and c_ax2[0] <= CL.C_3D[0] and CL.C_3D[1] <= c_ax2[1]
          and (c_ax2[1] - c_ax2[0]) < (c_ax1[1] - c_ax1[0]),
          f"[{ff(c_ax2[0])},{ff(c_ax2[1])}]")
    check("C-90 tail honesty (computed comparison, the verifier's deepening idiom): "
          "a_2800(t) - a_1400(t) in [0, published-tail(1400)] for the deep targets",
          all(0 <= KER2[t] - KER1[t] <= TL1[t] for t in (A8, A16)))
    kp350 = CL.kernel_pass({A8: 350, (10, 0, 0): 350, A16: 350}, 350)
    k350 = kp350["ker"]
    check("API-FIDELITY PROBE (beyond every sealed grid): UNTESTED target (10,0,0) at "
          "M=350 lies STRICTLY between a(8,0,0) and a(16,0,0) -- the kernel is the "
          "definition's own monotone sum, not a table of sealed targets",
          k350[A8] < k350[(10, 0, 0)] < k350[A16],
          f"a(8)={ff(k350[A8])} a(10)={ff(k350[(10, 0, 0)])} a(16)={ff(k350[A16])}")
    ev = CL.critical_evidence_3d(M=350)
    check("C-90 evidence tier (the class_verdict instrument at M=350): the 1/d-class "
          "signature emerges at reduced depth (increment ratio in INV, outside LOG and "
          "LIN, deepening-stabilized) -- computed booleans, no sealed number consumed",
          ev["inv_class_signature"],
          f"r24={ff(ev['ratio_24'])} stab={ff(ev['stabilization'])}")

    # ================================================================ H: cross-dimension
    disc = CL.discriminator(K2=6000, K1=80000)
    r2, r2b, r1 = disc["ratio_D2"], disc["ratio_D2_deeper"], disc["ratio_D1"]
    check("C-87 D=2 venue at ITS computed mu_c=1/4: increment ratios == sealed "
          "[0.967771, 0.986029], [0.953910, 1.028887], in the LOG window [4/5, 5/4]",
          (ff(r2[0]), ff(r2[1])) == ("0.967771", "0.986029")
          and (ff(r2b[0]), ff(r2b[1])) == ("0.953910", "1.028887")
          and disc["D2_in_log_window"],
          f"[{ff(r2[0])},{ff(r2[1])}]")
    check("C-87 D=1 venue at ITS computed mu_c=1/2: increment ratio == sealed "
          "[1.940862, 2.047476], in the LIN window [9/5, 11/5]",
          (ff(r1[0]), ff(r1[1])) == ("1.940862", "2.047476")
          and disc["D1_in_lin_window"],
          f"[{ff(r1[0])},{ff(r1[1])}]")
    check("C-87 the three windows are pairwise DISJOINT and each venue's critical class "
          "lands in its own: D=3 INV (pair 4->8 above), D=2 LOG, D=1 LIN == sealed S5 "
          "(one instrument, three venues; classes differ across earned dimension)",
          disc["windows_disjoint"] and r48[1] <= CL.INV_HI and r48[0] >= CL.INV_LO
          and disc["D2_in_log_window"] and disc["D1_in_lin_window"])
    S8, t8 = CL.series_target_2d(Fraction(1, 6), 8, 0, 260)
    S10, t10 = CL.series_target_2d(Fraction(1, 6), 10, 0, 260)
    rr16 = CL.ratio_interval((S10, t10), (S8, t8))
    check("C-87 mu=1/6 SUBCRITICAL on the D=2 venue == sealed [0.063870, 0.063870] "
          "(exponential there, critical on D=3: mu_c is the venue's own number)",
          (ff(rr16[0]), ff(rr16[1])) == ("0.063870", "0.063870")
          and rr16[1] <= 1 - CL.MARGIN_LT1,
          f"r=[{ff(rr16[0])},{ff(rr16[1])}]")

    # ================================================================ I: observation entry
    w3, w3_cells, w3_idx = CL.world_venue(3)
    src3 = w3_idx[(0, 0, 0)]
    CT3 = CL.walk_counts(w3, src3, 8)
    wrap_ok = True
    for k in (5, 8):
        for i, (x, y, z) in enumerate(w3_cells):
            free = 0
            for wx in (-3, -2, -1, 0, 1, 2, 3):
                for wy in (-3, -2, -1, 0, 1, 2, 3):
                    for wz in (-3, -2, -1, 0, 1, 2, 3):
                        free += DP.get((k, x + wx * 3, y + wy * 3, z + wz * 3), 0)
            if CT3[k][i] != free:
                wrap_ok = False
    check("OBSERVATION-ENTRY PROBE (a NEW venue, n=3, off every sealed grid): enters "
          "through the provenance gate; torus walk counts == wrap-summed Z^3 counts "
          "k <= 8 (the universal-cover identity re-earned on the new venue)",
          wrap_ok and w3.provenance == CL.WORLD_VENUE_PROVENANCE)
    dist3 = CL.bfs_venue(w3, src3)
    check("OBSERVATION-ENTRY PROBE: BFS distance == earned separation (L1 with wraps) on "
          "the new venue",
          all(dist3[i] == sum(min(abs(c), 3 - abs(c)) for c in w3_cells[i])
              for i in range(w3.n)))
    loc3 = CL.mu_c_of(w3, certify="full")
    check("OBSERVATION-ENTRY PROBE: the new venue's OWN mu_c located full-tier == 1/6 "
          "(computed on this venue, not imported from the sealed one)",
          loc3.get("located") and loc3["mu_c"] == Fraction(1, 6))
    mu = Fraction(1, 12)
    K = 40
    VS = CL.venue_series(w3, mu, src3, list(range(w3.n)), K)
    res = CL.resolvent_exact(w3.adj, mu, src3)
    check("OBSERVATION-ENTRY PROBE: on the new venue the EXACT resolvent lies inside "
          "[partial sum, partial sum + geometric tail] on every node at mu=1/12 (the "
          "series machinery verified against exact linear algebra -- the sealed t44a S2 "
          "gate re-earned on an untested venue)",
          all(VS[t][0] <= res[t] <= VS[t][0] + VS[t][1] for t in range(w3.n)))
    v_new = CL.class_verdict(w3, mu)
    check("OBSERVATION-ENTRY PROBE: the verdict triple on the new venue at mu=1/12 is "
          "(EXPONENTIAL, -, -) against the venue's own computed mu_c",
          (v_new["exponential"], v_new["critical"], v_new["divergent"])
          == (True, False, False) and v_new["mu_c"] == Fraction(1, 6))


if __name__ == "__main__":
    import time

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

    print("VALIDATE THE CLASSES LAYER (C-87/C-90) THROUGH model/classes.py")
    print("=" * 78)
    t0 = time.time()
    run_classes_checks(check)
    print("=" * 78)
    print(f"  CLASSES: {n_pass} PASS, {n_fail} FAIL   [{time.time() - t0:.1f} s]")
    _sys.exit(0 if n_fail == 0 else 1)
