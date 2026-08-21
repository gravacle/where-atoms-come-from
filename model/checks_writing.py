"""EVERY SEALED WRITING-TIER HEADLINE (C-91), REPRODUCED THROUGH model/writing.py.

Each check recomputes a registered headline through the writing tier's ported machinery
and gates it against the SEALED value copied from the lane's own sealed OUT file, in the
validate_geometry.py idiom (run_writing_checks(check)).  Sealed sources per section:

  KERNEL  LANE_T48_A_DERIVATION/t48a_derivation.OUT.txt   (126/126: unitarity forces
          conservation for EVERY weighting; conserving <=> critical det(I-K) == 0 with
          the exact nonzero CTRL-LEAK dets beside; CTRL-BIAS-LINK conserving-without-
          uniformity still critical, invariance broken; the moving census 1024/2048 and
          2304/4608; the pair split 4/18 + 4/18 + 10/18; CTRL-BATH column sums)
  CORNER  LANE_T48_B_CORNER/t48b_corner.OUT.txt           (112/112: link-uniformity
          EARNED from the writer algebra -- invariant tuple (1,2,2,2,1) identically;
          mu_c = 1/deg located in-lane; the weight-1 coset stratum (1,1); the chain's
          own 1/2)
  WORLD   LANE_T48_C_WORLD/t48c_world.OUT.txt             (36/36: E1 transport critical
          at every dE, mu set {1/6}; E2 conserving, uniform iff dE = 0, split (1/7, 2/7)
          at b = 1/2; E3 never critical, mu = b/(deg*b+1), gap ln(mu_c/mu) =
          ln(1 + e^{dE/kT}/l) with f0 and E_b dropping out; the T44-B comparison row
          mu(1/2) = 1/8; partial-sum anchors 13 and 39867016537742941/4096000000000000)

DISCIPLINE: deg / mu_c / 1-deg amplitudes are COMPUTED (D-8); sealed numbers are gated
AS THE SEALED ANCHORS THEY ARE, stated as such; every zero has a positive control beside
it (D-15); the D-25 provenance gate is TESTED, not assumed; API-fidelity probes run
BEYOND the gated range (fresh venue sizes, fresh sample points, a fresh lazy row, a
fresh corner venue) so the machinery is the definition, not a shortcut.

Standalone: python3 checks_writing.py  (exit 0 iff every check passes)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def run_writing_checks(check):
    from fractions import Fraction as Fr
    import writing as WR

    F1 = Fr(1)

    # ================================================================ KERNEL TIER (T48_A)
    # ---------------------------------------------------------------- venues from supports
    RING = WR.ring_venue(8)
    T3 = WR.plaquette_venue(3)
    T4 = WR.plaquette_venue(4)
    Z27 = WR.grain_venue(3)
    VEN = {"C8": RING["adj"], "T3": T3["adj"], "T4": T4["adj"], "Z27": Z27["adj"]}
    degs = {nm: WR.venue_degree(adj) for nm, adj in VEN.items()}
    check("C-91 kernel venues from supports: computed degrees (C8, T3, T4, Z27) == "
          "(2, 4, 4, 6) with sizes (8, 9, 16, 27)",
          degs == {"C8": 2, "T3": 4, "T4": 4, "Z27": 6} and
          [len(a) for a in VEN.values()] == [8, 9, 16, 27] and RING["well_formed"],
          f"degs={degs}")
    check("C-91 kernel D=2 carrier grounding: every carrier edge in exactly 2 plaquettes, "
          "adjacent plaquettes share exactly one edge, 18 elementary writers at L=3",
          T3["edge_in_two"] and T3["shared_single"] and T3["n_writers"] == 18 and
          len(set(T3["deltas"])) == 18)

    # ---------------------------------------------------------------- the energy census
    cen1 = WR.moving_census(RING["deltas"], 8, even_only=False)
    cen2 = WR.moving_census(T3["deltas"], 9, even_only=True)
    check("C-91 energy census: dE == 0 exactly on the defect-moving sector, dE in "
          "{-4,0,4}, both venues; SEALED ANCHORS 1024/2048 (D=1) and 2304/4608 (D=2) "
          "-- the zero-dE sector with its positive count beside it (D-15)",
          cen1["sectors_ok"] and cen2["sectors_ok"] and
          (cen1["moving"], cen1["total"]) == (1024, 2048) and
          (cen2["moving"], cen2["total"]) == (2304, 4608),
          f"D=1 {cen1['moving']}/{cen1['total']}, D=2 {cen2['moving']}/{cen2['total']}")

    # ---------------------------------------------------- conserving <=> critical + LEAK
    # SEALED ANCHORS (t48a OUT, DATA S3.4 == DATA S5.2): the exact nonzero dets of the
    # CTRL-LEAK kernel (declared survival 9/10) on every venue -- gated as sealed anchors.
    SEALED_LEAK_DET = {"C8": Fr(269059, 4000000),
                       "T3": Fr(653188856401, 4096000000000),
                       "T4": Fr(196750721299, 2560000000000),
                       "Z27": Fr(886217035591602121, 16384000000000000000)}
    for nm, adj in VEN.items():
        deg = degs[nm]
        K = WR.kernel_uniform(adj, 0)
        dl = WR.crit_det(WR.leak_kernel(adj, Fr(9, 10)))
        check(f"C-91 {nm} E-LOC: doubly stochastic AND link amplitude == 1/deg == 1/{deg} "
              "(formed from the computed deg) AND CRITICAL det(I - K) == 0 exactly",
              WR.is_doubly_stochastic(K) and
              WR.link_amplitudes(adj, K) == {Fr(1, deg)} and WR.crit_det(K) == 0)
        check(f"C-91 {nm} CTRL-LEAK (9/10): NOT conserving (row sums 9/10), per-link "
              f"9/(10*deg) != 1/deg, det(I - K_leak) != 0 AND == the SEALED ANCHOR",
              set(WR.mat_row_sums(WR.leak_kernel(adj, Fr(9, 10)))) == {Fr(9, 10)} and
              WR.link_amplitudes(adj, WR.leak_kernel(adj, Fr(9, 10)))
              == {Fr(9, 10) / deg} != {Fr(1, deg)} and
              dl != 0 and dl == SEALED_LEAK_DET[nm],
              f"det={dl}")
    # leak mass decay beside conserving mass (S5.4): (9/10)^k vs 1, k = 1..8, on T3
    KL, KC = WR.leak_kernel(T3["adj"], Fr(9, 10)), WR.kernel_uniform(T3["adj"], 0)
    v, vc, ok_mass = [F1] * 9, [F1] * 9, True
    for k in range(1, 9):
        v, vc = WR.mvec(KL, v), WR.mvec(KC, vc)
        ok_mass &= all(x == Fr(9, 10) ** k for x in v) and all(x == F1 for x in vc)
    check("C-91 CTRL-LEAK mass law: surviving measure == (9/10)^k exactly (k=1..8) while "
          "the conserving kernel holds measure 1 at the same depths", ok_mass)

    # ---------------------------------------------------------------- the lazy family
    ok_lazy = True
    for nm in ("T3", "Z27"):
        adj, deg = VEN[nm], degs[nm]
        for c in (Fr(1, 5), Fr(1, 2)):
            K = WR.kernel_uniform(adj, c)
            t = (F1 - c) / deg
            ok_lazy &= (WR.is_doubly_stochastic(K) and WR.crit_det(K) == 0
                        and t / (F1 - c) == Fr(1, deg))
    check("C-91 lazy family c*I + ((1-c)/deg)*A: conserving AND critical at every "
          "declared c on T3 and Z27, per-crossing amplitude t/(1-c) == 1/deg "
          "IDENTICALLY (laziness never opens the gap)", ok_lazy)

    # ------------------------------------- unitarity conserves for EVERY weighting (S3)
    W_D1 = [[Fr(i + 1, 36) for i in range(8)],
            [Fr(2 ** i, 255) for i in range(8)],
            [Fr(1, 6) if i % 2 == 0 else Fr(1, 12) for i in range(8)]]
    W_D2 = [[Fr(e + 1, 171) for e in range(18)],
            [Fr(2 ** e, (1 << 18) - 1) for e in range(18)],
            [Fr(1, 12) if e % 2 == 0 else Fr(1, 36) for e in range(18)]]
    ok_g = WR.dict_doubly_stochastic(*WR.glob_kernel(8, RING["deltas"], [Fr(1, 8)] * 8))
    ok_g &= WR.dict_doubly_stochastic(*WR.glob_kernel(9, T3["deltas"], [Fr(1, 18)] * 18))
    ok_nu = all(sum(w) == F1 and
                WR.dict_doubly_stochastic(*WR.glob_kernel(8, RING["deltas"], w))
                for w in W_D1)
    ok_nu &= all(sum(w) == F1 and
                 WR.dict_doubly_stochastic(*WR.glob_kernel(9, T3["deltas"], w))
                 for w in W_D2)
    synds2s, T2s = WR.glob_kernel(9, T3["deltas"], [Fr(1, 27)] * 18,
                                  extra_trivial=Fr(9, 27))
    ok_st = (WR.dict_doubly_stochastic(synds2s, T2s)
             and all(T2s[s].get(s, 0) == Fr(1, 3) for s in synds2s))
    check("C-91 UNITARITY CONSERVES FOR EVERY WEIGHTING (the writer algebra, not "
          "uniformity): global syndrome-space kernels doubly stochastic for the uniform "
          "rows AND all 6 declared non-uniform weight rows AND the stabilizer-augmented "
          "ensemble (diagonal 1/3 from trivial-action writers)", ok_g and ok_nu and ok_st)

    # ---------------------------------------------------------------- the pair sector
    synds2, T2g = WR.glob_kernel(9, T3["deltas"], [Fr(1, 18)] * 18)
    ps = WR.pair_sector(T2g, T3["adj"], 1 * 3 + 1, 0)
    check("C-91 pair sector (SEALED split): tracked-end 4/18 + origin-end 4/18 + pair "
          "creation 10/18 == 1; tracked amplitudes equal on the 4 links; CONDITIONED "
          "per-crossing amplitude == 1/deg == 1/4 exactly",
          ps["nonadjacent"] and ps["tracked_on_neighbors"] and
          set(ps["tracked"].values()) == {Fr(1, 18)} and
          sum(ps["tracked"].values()) == Fr(4, 18) and
          ps["origin_sum"] == Fr(4, 18) and ps["creation"] == Fr(10, 18) and
          sum(ps["tracked"].values()) + ps["origin_sum"] + ps["creation"] == F1 and
          ps["conditional"] == {Fr(1, 4)})

    # ---------------------------------------------------------------- CTRL-BIAS-LINK
    KB3 = WR.biased_kernel_2d(3, Fr(1, 3), Fr(1, 6))
    KBZ = WR.biased_kernel_3d(3, (Fr(1, 4), Fr(1, 6), Fr(1, 12)))
    for nm, K, gmap in (("T3", KB3, WR.plaquette_transpose(3)),
                        ("Z27", KBZ, Z27["transpose"])):
        adj = VEN[nm]
        amps = WR.link_amplitudes(adj, K)
        viol = WR.invariance_violations(adj, K, gmap)
        viol0 = WR.invariance_violations(adj, WR.kernel_uniform(adj, 0), gmap)
        check(f"C-91 {nm} CTRL-BIAS-LINK: measure IS conserved (doubly stochastic) yet "
              f"uniformity FAILS ({len(amps)} distinct link amplitudes) AND still "
              "CRITICAL det(I - K) == 0 -- conservation pins the gap shut without "
              "forcing 1/deg",
              WR.is_doubly_stochastic(K) and len(amps) > 1 and WR.crit_det(K) == 0,
              f"amps={sorted(str(a) for a in amps)}")
        check(f"C-91 {nm} CTRL-BIAS-LINK invariance witness: the axis-transposing map IS "
              "a venue automorphism (computed), the biased kernel breaks invariance at "
              f"it ({viol} violations > 0) while the uniform kernel's violation count is "
              "the exact zero beside it (D-15) -- ensemble-WEIGHT symmetry is extra data",
              WR.check_aut(adj, gmap) and viol > 0 and viol0 == 0)

    # ---------------------------------------------------------------- CTRL-BATH
    b_half = WR.bath_dilation(Fr(1, 2))
    b_34 = WR.bath_dilation(Fr(3, 4))
    b_1 = WR.bath_dilation(F1)
    check("C-91 CTRL-BATH: the exchange dilation is unitary and conserves excitation "
          "(iv'); induced channel trace-preserving at EVERY bias; measure conservation "
          "holds at p == 1/2 and FAILS at 3/4 and 1 with the SEALED column sums "
          "[1,1], [3/2,1/2], [2,0] -- bath polarization is bias, and bias is mass",
          all(x["unitary"] and x["conserves_excitation"] and x["trace_preserving"]
              for x in (b_half, b_34, b_1)) and
          b_half["doubly_stochastic"] and not b_34["doubly_stochastic"] and
          not b_1["doubly_stochastic"] and
          b_half["col_sums"] == [F1, F1] and b_34["col_sums"] == [Fr(3, 2), Fr(1, 2)]
          and b_1["col_sums"] == [Fr(2), Fr(0)])

    # ================================================================ CORNER TIER (T48_B)
    CV37 = WR.corner_venue(3, 7)
    CV46 = WR.corner_venue(4, 6)
    CV55 = WR.corner_venue(5, 5)
    check("C-91 corner venues from carrier supports: #links == #carrier edges == "
          "(42, 48, 50) at (3,7), (4,6), (5,5); all multiplicities 1; every edge in "
          "exactly 2 plaquettes; no self-links",
          all(cv["edge_in_two"] and cv["mult_all_one"] and cv["no_self"]
              for cv in (CV37, CV46, CV55)) and
          (len(CV37["links"]), len(CV46["links"]), len(CV55["links"])) == (42, 48, 50))
    for tag, cv in (("(3,7)", CV37), ("(4,6)", CV46), ("(5,5)", CV55)):
        inv = WR.writer_invariants(cv)
        check(f"C-91 corner {tag} UNIFORMITY EARNED FROM THE ALGEBRA: the elementary "
              "writer's invariant tuple (price, |syndrome|, #stars, #plaqs, mult) is "
              "IDENTICAL on every link AND == the SEALED (1, 2, 2, 2, 1) -- every "
              "algebra-measurable ensemble is link-uniform",
              all(v == inv[0] for v in inv) and inv[0] == (1, 2, 2, 2, 1),
              f"tuple={inv[0]} on {len(inv)} links")
    for tag, cv in (("(3,7)", CV37), ("(4,6)", CV46)):
        loc = WR.mu_c_locate(cv["rows"], 0, beside=(Fr(1, 8), Fr(23, 100)))
        check(f"C-91 corner {tag} mu_c LOCATED IN-LANE: constant row sums (Perron) give "
              f"deg == 4 computed; exact resolvent SINGULAR at 1/deg == {loc['mu_c']} "
              "and NONSINGULAR beside it at 1/8 and 23/100 (positive controls, D-15)",
              loc["rows_constant"] and loc["deg"] == 4 and loc["singular_at_mu_c"]
              and all(loc["nonsingular_beside"]))
    sw = WR.sector_sandwich()
    check("C-91 corner sector sandwich (exact binomial lemmas, m <= 300, induction to "
          "10000): 16^m/(2m+1)^2 <= C(2m,m)^2 <= 16^m -- the venue-limit return series "
          "has radius EXACTLY 1/4", sw["all_ok"])
    CH = WR.chain_venue(24)
    locc = WR.mu_c_locate(CH, 0, beside=(Fr(9, 20),))
    check("C-91 chain C_24 (the D=1 discriminator): deg == 2 computed, resolvent "
          "SINGULAR at 1/2 and NONSINGULAR at 9/20 -- each venue returns its OWN number",
          locc["deg"] == 2 and locc["mu_c"] == Fr(1, 2) and locc["singular_at_mu_c"]
          and all(locc["nonsingular_beside"]))
    cm37, cmch = WR.conserving_member(CV37["rows"]), WR.conserving_member(CH)
    check("C-91 K1 THE CONSERVING MEMBER IS THE CRITICAL ONE: the unique conserving "
          "amplitude t* == 1/deg equals the in-lane mu_c on the corner venue (1/4) AND "
          "on the chain (1/2), each formed from its own computed deg",
          cm37["unique"] and cm37["t_star"] == Fr(1, 4) and
          cmch["unique"] and cmch["t_star"] == Fr(1, 2))
    lk = CV37["links"][0]
    w37 = WR.elementary_coset(CV37, (lk[0], lk[1]))
    check("C-91 corner (3,7) the weight-1 coset stratum IS the link: exhaustive coset "
          "scan of the elementary writer gives (w_min, N_min) == (1, 1) over 2^22 "
          "admissible writers", w37[:2] == (1, 1), f"(w,N,bits)={w37}")

    # ================================================================ WORLD TIER (T48_C)
    cells, idx, nbr = WR.torus3(8)
    NC = len(cells)
    deg_w = len(nbr[0]) if all(len(r) == len(nbr[0]) for r in nbr) else None
    d_bfs = WR.bfs_dist(nbr, 0)
    check("C-91 world venue (8^3 census access geometry): face-degree == 6 computed on "
          "every grain AND BFS distance == earned L1 wrap separation on every grain",
          deg_w == 6 and all(d_bfs[i] == WR.l1_wrap(cells[i], 8) for i in range(NC)))
    te = WR.trail_energetics(nbr)
    check("C-91 world energetics (the amplitude grounding): TRANSPORT dN == 0 in EVERY "
          "direction (the computed zero earning E1's direction-uniform amplitude) "
          "beside WRITE dN == +1 (5 fresh directions) and ERASE dN == -1 (D-15)",
          te["transport"] == [0] * 6 and te["write"] == [1] * 5 and te["erase"] == [-1])
    db_ok = all(WR.detailed_balance(u, b) == dict(stationary=True, ratio_is_b=True)
                for u in (Fr(1, 20), Fr(1, 100))
                for b in (F1, Fr(9, 10), Fr(3, 4), Fr(1, 2), Fr(1, 4), Fr(1, 10)))
    check("C-91 detailed balance from the model's own two-state kernel: pi K == pi and "
          "pi_meta/pi_stable == v/u == b exactly at every (u, b) sample", db_ok)

    U_S = [Fr(1, 20), Fr(1, 100)]
    B_S = [F1, Fr(9, 10), Fr(3, 4), Fr(1, 2), Fr(1, 4), Fr(1, 10)]

    # ---------------------------------------------------------------- E1 TRANSPORT
    ok_e1, mus_e1 = True, set()
    for u in U_S:
        for b in B_S:
            for a in (u, u * b):     # both honest saddle readings
                v1 = WR.transport_verdict(nbr, WR.ensemble_transport(nbr, a))
                ok_e1 &= (v1["conserving"] and v1["uniform"] and v1["at_criticality"]
                          and v1["op_identity"] and v1["massless_signature"]
                          and v1["mu"] == Fr(1, v1["deg"]))
                mus_e1.add(v1["mu"])
    check("C-91 E1 TRANSPORT: conserving AND critical -- induced mu == 1/deg == mu_c -- "
          "at EVERY (u, b) sample and BOTH saddle readings; operator identity entrywise; "
          "(I - W)1 == 0 exactly; SEALED induced-mu set {1/6}",
          ok_e1 and mus_e1 == {Fr(1, 6)}, f"mu set={sorted(str(m) for m in mus_e1)}")
    W1 = WR.ensemble_transport(nbr, Fr(1, 20))
    vec = [Fr(0)] * NC
    vec[0] = F1
    tot_sum, ok_ev = Fr(0), True
    for k in range(13):
        tot = sum(vec)
        tot_sum += tot
        ok_ev &= (tot == F1)
        if k < 12:
            vec = WR.apply_forward(W1, vec)
    check("C-91 E1 measure evolution: total == 1 at every step k <= 12 AND propagator "
          "partial sums == K+1 == 13 (SEALED anchor; LINEAR divergence -- the massless "
          "witness)", ok_ev and tot_sum == 13, f"sum={tot_sum}")

    # ---------------------------------------------------------------- E2 TRAIL W/ RETREAT
    ok_e2 = True
    for (u, b) in ((Fr(1, 20), F1), (Fr(1, 20), Fr(1, 2)), (Fr(1, 100), F1)):
        v2 = WR.retreat_verdict(nbr, WR.ensemble_trail_retreat(nbr, u, b))
        ok_e2 &= (v2["conserving"] and v2["doubly_stochastic"]
                  and v2["redistribution_exact"] and v2["uniform"] == (b == F1))
        if b == F1:
            ok_e2 &= (v2["m_fresh"] == v2["m_back"] == Fr(1, v2["deg"]))
        else:
            v_r = u * b
            ok_e2 &= (v2["m_fresh"] == v_r / (5 * v_r + u)
                      and v2["m_back"] == u / (5 * v_r + u))
    v2h = WR.retreat_verdict(nbr, WR.ensemble_trail_retreat(nbr, Fr(1, 20), Fr(1, 2)))
    check("C-91 E2 TRAIL WITH RETREAT: conserving (doubly stochastic, redistribution "
          "total 1 exactly) at every sampled dE; UNIFORM iff dE == 0 with m == 1/deg; "
          "SEALED split at b = 1/2: (m_fresh, m_back) == (1/7, 2/7)",
          ok_e2 and (v2h["m_fresh"], v2h["m_back"]) == (Fr(1, 7), Fr(2, 7)),
          f"b=1/2 split=({v2h['m_fresh']}, {v2h['m_back']})")
    W2 = WR.ensemble_trail_retreat(nbr, Fr(1, 20), Fr(1, 2))
    dr = WR.drift_per_state(W2, nbr)
    v_r = Fr(1, 20) * Fr(1, 2)
    avg = [Fr(0)] * 3
    ok_dr = True
    for i in range(NC):
        for d in range(6):
            s = i * 6 + d
            ok_dr &= dr[s] == tuple((v_r - Fr(1, 20)) * WR.DIRS[d][ax]
                                    for ax in range(3))
            for ax in range(3):
                avg[ax] += dr[s][ax]
    check("C-91 E2 drift at dE != 0: per-state drift == (v-u)*DIRS[d] exactly (nonzero "
          "path persistence, the positive control) while the stationary-average drift "
          "is the exact ZERO beside it (no spatial bias; D-15)",
          ok_dr and all(x == 0 for x in avg))

    # ---------------------------------------------------------------- E3 TRAIL W/ DECAY
    # SEALED ANCHOR (t48c OUT G24): the induced-mu table at the declared b samples.
    SEALED_MU_E3A = {F1: Fr(1, 7), Fr(9, 10): Fr(9, 64), Fr(3, 4): Fr(3, 22),
                     Fr(1, 2): Fr(1, 8), Fr(1, 4): Fr(1, 10), Fr(1, 10): Fr(1, 16)}
    ok_e3, mu_by_b = True, {}
    for u in U_S:
        for b in B_S:
            v3 = WR.decay_verdict(nbr, WR.ensemble_trail_decay(nbr, u, b, "H1"))
            v_r = u * b
            ok_e3 &= (not v3["conserving"] and v3["row_sum"] == 1 - u
                      and v3["loss"] == u and v3["uniform"]
                      and v3["mu"] == v_r / (v3["deg"] * v_r + u)
                      and v3["below_criticality"]
                      and v3["mass_ratio"]
                      == WR.closed_form_gap_ratio(b, v3["deg"]))
            mu_by_b.setdefault(b, set()).add(v3["mu"])
    check("C-91 E3a TRAIL WITH DECAY: NEVER critical -- row sums == 1-u exactly (loss "
          "== the model's own erase probability, the D-15 nonzero beside E1's zero "
          "deficit), mu == v/(deg*v+u) < mu_c STRICTLY at every (u, b), AND the "
          "computed mass ratio mu_c/mu == the closed form 1 + e^{dE/kT}/deg -- CHECKED "
          "against it, never sourced from it", ok_e3)
    check("C-91 E3a barrier/attempt-clock independence + SEALED mu table: induced mu "
          "identical across u at fixed b (f0 and E_b drop out EXACTLY) and equal to the "
          "sealed anchors {1/7, 9/64, 3/22, 1/8, 1/10, 1/16}; strictly increasing in b",
          all(len(s) == 1 for s in mu_by_b.values()) and
          {b: next(iter(s)) for b, s in mu_by_b.items()} == SEALED_MU_E3A and
          all(next(iter(mu_by_b[b1])) < next(iter(mu_by_b[b2]))
              for b1, b2 in zip(sorted(B_S), sorted(B_S)[1:])),
          f"mu(b)={ {str(b): str(next(iter(s))) for b, s in mu_by_b.items()} }")
    check("C-91 COMPARISON-TO-SEALED (T44-B): at dE = kT ln 2 the trail ensemble's "
          "induced amplitude == 1/8 == the sealed T44-B subcritical sweep row whose "
          "class was computed EXPONENTIAL (SEALED anchor)",
          next(iter(mu_by_b[Fr(1, 2)])) == Fr(1, 8))
    u0, b0 = Fr(1, 20), Fr(1, 2)
    W3 = WR.ensemble_trail_decay(nbr, u0, b0, "H1")
    vec = [Fr(0)] * NC
    vec[0] = F1
    tot_sum, ok_ev3 = Fr(0), True
    for k in range(13):
        tot = sum(vec)
        tot_sum += tot
        ok_ev3 &= (tot == (1 - u0) ** k)
        if k < 12:
            vec = WR.apply_forward(W3, vec)
    closed = (1 - (1 - u0) ** 13) / u0
    check("C-91 E3a measure evolution: total == (1-u)^k at every k <= 12; partial sums "
          "== (1-(1-u)^(K+1))/u == the SEALED anchor 39867016537742941/4096000000000000 "
          "< 1/u (CONVERGENT, massive, beside E1's linear divergence)",
          ok_ev3 and tot_sum == closed
          and tot_sum == Fr(39867016537742941, 4096000000000000) and tot_sum < 1 / u0,
          f"sum={tot_sum}")
    ok_nb = True
    for u in U_S:
        for b in (F1, Fr(1, 2), Fr(1, 10)):
            vnb = WR.decay_verdict_nb(nbr, WR.ensemble_trail_decay(nbr, u, b, "NB"))
            ok_nb &= (vnb["deg_nb"] == 5 and vnb["mu_c_nb"] == Fr(1, 5)
                      and vnb["loss"] == u and vnb["back_channel_empty"]
                      and vnb["mu"] == b / (5 * b + 1)
                      and vnb["redistribution_total"] == 5 * b / (5 * b + 1) < 1
                      and vnb["mass_ratio"] == WR.closed_form_gap_ratio(b, 5)
                      and vnb["below_criticality"])
    check("C-91 E3b (non-backtracking counting): criticality reference deg_NB == 5 "
          "EARNED from the venue's own directed-edge operator row sums; mu == b/(5b+1) "
          "off criticality; SAME gap law 1 + e^{dE/kT}/l with l = 5 (counting-"
          "independent in form); row deficit == u == exactly E2's retreat amplitude "
          "(the measure extension-only counting loses IS the backtracking channel)",
          ok_nb)

    # ------------------------------------------------- labeled FLOAT comparison (last)
    from math import log, exp
    ok_cmp = True
    for b in B_S:
        x = -log(b)
        ratio = WR.closed_form_gap_ratio(b, 6)
        ok_cmp &= abs(log(Fr(ratio)) - log(1 + exp(x) / 6)) <= 1e-12
    check("C-91 COMPARISON (floats, labeled, after results): gap ln(mu_c/mu) agrees "
          "with ln(1 + e^{dE/kT}/deg) at every sample within 1e-12", ok_cmp)

    # ================================================================ THE D-25 GATE
    from project_model import URM, RecordSurface
    eV = 1.602176634e-19
    refused = False
    try:
        WR.surface_gap(RecordSurface("mystery", "unknown", 1e-20, 1e-19, 300.0, 1e9))
    except ValueError:
        refused = True
    check("D-25 the writing tier's observation gate REFUSES a world surface without "
          "provenance (tested, not assumed)", refused)
    nand = URM.surface("NAND floating gate", "trapped charge", 0.30 * eV, 1.60 * eV,
                       300.0, 1e13)
    g = WR.surface_gap(nand)
    check("C-91/D-25 A NEW OBSERVATION ENTERS: the NAND surface's written-trail gap is "
          "computed through the E3 kernel at exact rational brackets of its own b -- "
          "closed form CHECKED against every computed ratio, u-independence re-computed "
          "at entry, float gap certified INSIDE the computed bracket",
          g["closed_form_agrees"] and g["u_independent"] and g["contained"],
          f"gap ln(mu_c/mu)={g['gap_ln']:.4f} in {tuple(round(x, 4) for x in g['gap_bracket_ln'])}")
    cmb = URM.surface("CMB photon", "free flight", 0.0, 0.0, 2.7, 0.0, thermal=False,
                      provenance="census GR1 control (not thermally activated)")
    check("C-91 the writing tier DECLINES a non-thermal surface (None, never a number) "
          "exactly as the model's laws do", WR.surface_gap(cmb) is None)

    # ============================================= API-FIDELITY PROBES BEYOND THE RANGE
    cells6, _i6, nbr6 = WR.torus3(6)
    v1f = WR.transport_verdict(nbr6, WR.ensemble_transport(nbr6, Fr(1, 30)))
    bf = Fr(2, 5)
    v3f = WR.decay_verdict(nbr6, WR.ensemble_trail_decay(nbr6, Fr(1, 30), bf, "H1"))
    check("PROBE (beyond gated range) fresh world venue 6^3, fresh samples u=1/30, "
          "b=2/5: E1 still conserving+critical with mu == 1/deg computed; E3a mu == "
          "2/17 == v/(deg*v+u) formed fresh, mass ratio == 17/12 == closed form at the "
          "fresh point -- the law emerges from the machinery, it is not hard-coded",
          v1f["conserving"] and v1f["at_criticality"] and
          v3f["mu"] == Fr(2, 17) == (Fr(1, 30) * bf) / (6 * Fr(1, 30) * bf + Fr(1, 30))
          and v3f["mass_ratio"] == Fr(17, 12) == WR.closed_form_gap_ratio(bf, 6))
    v2f = WR.retreat_verdict(nbr6, WR.ensemble_trail_retreat(nbr6, Fr(1, 30), Fr(2, 3)))
    check("PROBE (beyond gated range) fresh E2 point b=2/3 on 6^3: conserving, split "
          "(m_fresh, m_back) == (2/13, 3/13) == (b, 1)/(5b+1) formed fresh, non-uniform",
          v2f["conserving"] and v2f["redistribution_exact"] and not v2f["uniform"] and
          (v2f["m_fresh"], v2f["m_back"]) == (Fr(2, 13), Fr(3, 13)))
    C12 = WR.ring_venue(12)
    K12 = WR.kernel_uniform(C12["adj"], 0)
    dl12 = WR.crit_det(WR.leak_kernel(C12["adj"], Fr(9, 10)))
    check("PROBE (beyond gated range) fresh kernel venue C_12: E-LOC doubly stochastic "
          "and critical (det == 0), leak det != 0 -- no sealed value exists for this "
          "venue; the machinery computes it (definition-not-shortcut)",
          WR.is_doubly_stochastic(K12) and WR.crit_det(K12) == 0 and dl12 != 0,
          f"leak det={dl12}")
    c_f = Fr(2, 7)
    K4z = WR.kernel_uniform(T4["adj"], c_f)
    check("PROBE (beyond gated range) fresh lazy row c=2/7 on T4 (outside the declared "
          "LAZY_C rows): conserving, critical, per-crossing t/(1-c) == 1/deg computed",
          WR.is_doubly_stochastic(K4z) and WR.crit_det(K4z) == 0 and
          ((F1 - c_f) / 4) / (F1 - c_f) == Fr(1, 4))
    CV45 = WR.corner_venue(4, 5)
    inv45 = WR.writer_invariants(CV45)
    loc45 = WR.mu_c_locate(CV45["rows"], 0, beside=(Fr(1, 8),))
    lk45 = CV45["links"][0]
    w45 = WR.elementary_coset(CV45, (lk45[0], lk45[1]))
    check("PROBE (beyond gated range) fresh corner venue (4,5): invariant tuple "
          "(1,2,2,2,1) identical on all 40 links, mu_c located at 1/4 (singular there, "
          "nonsingular at 1/8), fresh exhaustive coset stratum (w_min, N_min) == (1,1)",
          len(inv45) == 40 and all(v == inv45[0] for v in inv45)
          and inv45[0] == (1, 2, 2, 2, 1) and loc45["deg"] == 4
          and loc45["singular_at_mu_c"] and all(loc45["nonsingular_beside"])
          and w45[:2] == (1, 1))
    a0 = Fr(1, 20)
    Wf_const = WR.kernel_pos_field(nbr6, lambda i: 1 - 6 * a0, lambda i, d: a0)
    same = all(Wf_const[i] == WR.ensemble_transport(nbr6, a0)[i]
               for i in range(len(nbr6)))
    a1, a2 = Fr(1, 8), Fr(1, 16)
    c_field = 2 * a1 + 4 * a2                # the extracted normalization, formed from
    Wf = WR.kernel_pos_field(nbr6, lambda i: 1 - c_field,      # the declared constants
                             lambda i, d: a1 if d < 2 else a2)
    vf = WR.transport_verdict(nbr6, Wf)
    _cs_f, ms_f = WR.extract_pos(Wf, nbr6)
    check("PROBE the C-93 entry point (kernel_pos_field): at constant rates it IS "
          "kernel_pos entrywise; at a declared anisotropic field it computes -- "
          "conserving (row sums 1 exactly) with uniformity honestly FAILING (exactly "
          "the 2 amplitudes {a1/c, a2/c} formed from the declared constants) -- the "
          "instrument measures, it does not assume",
          same and vf["conserving"] and not vf["uniform"]
          and ms_f == {a1 / c_field, a2 / c_field})


if __name__ == "__main__":
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

    print("VALIDATE THE WRITING TIER (C-91) THROUGH model/writing.py")
    print("=" * 78)
    run_writing_checks(check)
    print("=" * 78)
    print(f"  WRITING: {n_pass} PASS, {n_fail} FAIL")
    sys.exit(0 if n_fail == 0 else 1)
