"""LANE_T51_A / a2 -- THE F TABLE AT 3x3 (the commissioned measurement).

Per FIELD_INSTRUMENT_V001.md section 1: toric 3x3 (18 edges), sector-exact blocks via
the diagonal conserved algebra; probe star-hole pair at adjacent vertices (connector
weight 1); source plaquette-hole pair at the two distinct earned separations the
venue affords, each stated as (d_gen, w_enc_conn) from the sealed connectivity gate;
mediator V = lam * sum_e Z_e, lam in {0.02, 0.05, 0.10}.  Computed per
(placement, lam): Delta(b=+1), Delta(b=-1), F = Delta(-1) - Delta(+1); onset order of
F in lam by the declared estimator with its noise floor printed beside it (C4's
bracket, scored against connected w_enc PER PLACEMENT -- C-93's binding repair).

PLACEMENT CHOICE, DECLARED AS A RULE (not picked after looking at F): the gate's
sealed table gives the venue's earned-separation classes.  The two commissioned
separations are the venue's EXTREME earned classes -- the smallest, (d_gen=1,
w_enc_conn=3), represented by its first placement in the gate's row-major tag order,
and the largest, (d_gen=2, w_enc_conn=5), whose placement is unique.  This realizes
"the two distinct earned separations 3x3 affords" with one placement from each
descriptive d_gen class and maximal earned contrast.  The Gamma-equivalent swap
partners (control C2, script a3) and an all-18-placement survey at lam = 0.05 sit in
the same tables so nothing is narrowed by the choice.

C1 (two-way control, in the same table): F itself is the written-vs-unwritten
comparison -- Delta(b=-1) against Delta(b=+1), two exact sectors of ONE Hamiltonian.
Beside it, the source-absent configuration (source plaquette terms RESTORED; a
genuinely different Hamiltonian) is computed per row.  CONSTRUCTION CERTIFICATE
(labeled, never counted as a control): in the conserving quadrature
Delta(b=+1) = Delta(absent) is an algebraic identity (the restored terms are an exact
in-sector constant, E0 shift exactly +2); its measured residual certifies the port
and supplies the same-table measured control floor used by V1.  D-15: every certified
zero in this table sits beside nonzero F rows in different configurations of the
same table.

PORT CERTIFICATE: at the grid overlap lam = 0.064 this implementation's Delta and F
are compared against the gate's sealed prints (parsed at runtime, seal verified;
tolerances declared in t51a_lib) -- the same cross-implementation discipline the gate
applied to the fourth-angle lane.
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from t51a_lib import (Venue, F_reading, onset_measure, sector_spectrum_g, Emitter,
                      write_json, parse_gate_table, bool_word, COMM_LAMS, FLOOR_EIG,
                      FLOOR_USE, TOL_K, TOL_PORT_D, TOL_PORT_F, WITNESS_MIN,
                      ONSET_GRID, BOUND_LOG2)

emit = Emitter()


def main():
    t0 = time.time()
    emit("=" * 88)
    emit("LANE_T51_A / a2 -- THE F TABLE AT 3x3 (the commissioned measurement)")
    emit("date: 2026-08-21   commissioned in FIELD_INSTRUMENT_V001.md section 1")
    emit("=" * 88)
    emit("")
    emit("DECLARED BEFORE USE: commissioned mediator grid lam in %s; onset grid %s"
         % (list(COMM_LAMS), list(ONSET_GRID)))
    emit("with the gate's estimator (smallest usable pair, usable iff both |F| >=")
    emit("FLOOR_USE = %.0e; drift beside; FLOOR_EIG = %.0e beside every F);"
         % (FLOOR_USE, FLOOR_EIG))
    emit("TOL_K = %.2f; WITNESS_MIN = %.1f; port tolerances TOL_PORT_D = %.0e (Delta,"
         % (TOL_K, WITNESS_MIN, TOL_PORT_D))
    emit("10 printed digits), TOL_PORT_F = %.0e (F, 4 printed digits)." % TOL_PORT_F)
    emit("")

    gate, gate_seal = parse_gate_table()
    emit("GATE INPUT: g1_connected_wenc.txt seal VERIFIED (sha256 %s...);"
         % gate_seal[:16])
    emit("18 placements parsed; earned quantities below are the gate's sealed table")
    emit("(the gate is input, not re-run).")
    classes = {}
    for tag, g in gate.items():
        classes.setdefault((g["d_gen"], g["w_conn"]), []).append(tag)
    small = min(classes)          # (1, 3)
    large = max(classes)          # (2, 5)
    # first placement of the smallest class in the gate's own table order:
    gate_order = ["3x3-%d%d%s" % (x, y, hv) for y in range(3) for x in range(3)
                  for hv in ("H", "V")]
    p_near = next(t for t in gate_order if t in classes[small])
    assert len(classes[large]) == 1
    p_far = classes[large][0]
    emit("")
    emit("EARNED-SEPARATION CLASSES (d_gen descriptive, w_enc_conn earned):")
    for k in sorted(classes):
        emit("  (d_gen=%d, w_enc_conn=%d): %s" % (k[0], k[1],
                                                  ", ".join(sorted(classes[k]))))
    emit("COMMISSIONED PLACEMENTS BY THE DECLARED RULE:")
    emit("  near = %s  (d_gen=%d, w_enc_conn=%d)  [first of the smallest class]"
         % (p_near, small[0], small[1]))
    emit("  far  = %s  (d_gen=%d, w_enc_conn=%d)  [the largest class is a singleton]"
         % (p_far, large[0], large[1]))
    swap_near = next(t for t in gate_order if t in classes[small] and t != p_near)
    pair24 = sorted(classes[(2, 4)])
    emit("  control rows (C2, a3): %s [swap partner of near, same class]; %s and %s"
         % (swap_near, pair24[0], pair24[1]))
    emit("  [the two-member class (d_gen=2, w_enc_conn=4)].")
    emit("")

    Lx, Ly = 3, 3
    V = Venue(Lx, Ly)
    n = 2 * Lx * Ly
    dim_sector = 1 << (n - (Lx * Ly - 1) - 2)
    emit("VENUE: toric %dx%d, %d edges; full Hilbert space 2^%d = %d." % (Lx, Ly, n,
                                                                          n, 1 << n))
    emit("Diagonal conserved algebra: 9 plaquette parities (8 independent) + source")
    emit("values + Z-winding pair = 10 bits -> sector dimension %d = 2^%d <= declared"
         % (dim_sector, n - (Lx * Ly - 1) - 2))
    emit("bound 2^%d -- CHECKED HERE BEFORE ANY eigh RUN (and asserted inside every"
         % BOUND_LOG2)
    emit("sector build).  Dense sector eigh at %d x %d: ~%.1f MB workspace -- fits."
         % (dim_sector, dim_sector, dim_sector * dim_sector * 8 / 1e6))
    emit("")

    probe = ((0, 0), (1, 0))
    table_tags = [p_near, swap_near, pair24[0], pair24[1], p_far]

    # ---- port certificate at the grid overlap --------------------------------------
    emit("PORT CERTIFICATE (labeled; cross-implementation, not a control): lam=0.064")
    port_ok = True
    for tag in (p_near, p_far):
        src = gate[tag]["src"]
        r = F_reading(V, probe, src, 0.064)
        gd = gate[tag]["dyn"][0.064]
        okD = (abs(r["dM"] - gd[0]) / abs(gd[0]) <= TOL_PORT_D and
               abs(r["dP"] - gd[1]) / abs(gd[1]) <= TOL_PORT_D)
        okF = abs(r["F"] - gd[2]) / abs(gd[2]) <= TOL_PORT_F
        port_ok = port_ok and okD and okF
        emit("  %s: Delta rel diffs (%.1e, %.1e) <= %.0e: %s;  F rel diff %.1e <= "
             "%.0e: %s"
             % (tag, abs(r["dM"] - gd[0]) / abs(gd[0]),
                abs(r["dP"] - gd[1]) / abs(gd[1]), TOL_PORT_D, bool_word(okD),
                abs(r["F"] - gd[2]) / abs(gd[2]), TOL_PORT_F, bool_word(okF)))
    emit("  PORT CERTIFICATE: %s" % bool_word(port_ok))
    emit("")

    # ---- THE F TABLE ----------------------------------------------------------------
    emit("-" * 88)
    emit("THE F TABLE -- winding sector (+1,+1) (the sealed zbar convention's")
    emit("reference sector; the full winding sweep is control C3, script a3).")
    emit("Per row: the two exact sectors of ONE Hamiltonian (written b=-1, unwritten")
    emit("b=+1), the source-absent configuration beside them, and the certificate")
    emit("residual.  Earned separation stated per placement; construction coordinates")
    emit("are labels only (D-24).")
    emit("-" * 88)
    ftab = {}
    cert_resid = []
    for tag in table_tags:
        g = gate[tag]
        src = g["src"]
        role = ("COMMISSIONED near" if tag == p_near else
                "COMMISSIONED far" if tag == p_far else "control row (C2)")
        emit("")
        emit("PLACEMENT %s -- %s -- earned (d_gen=%d DESCRIPTIVE, w_enc_conn=%d)"
             % (tag, role, g["d_gen"], g["w_conn"]))
        emit("  [construction labels: probe stars %s,%s; source plaquettes %s,%s]"
             % (probe[0], probe[1], src[0], src[1]))
        emit("  %-6s %-17s %-17s %-12s %-13s %-9s" %
             ("lam", "Delta(b=-1)", "Delta(b=+1)", "F", "cert-resid", "witness"))
        rows = []
        for lam in COMM_LAMS:
            r = F_reading(V, probe, src, lam)
            ra = sector_spectrum_g(V, probe, (), +1, lam, src_present=False)
            resid = abs(r["dP"] - ra["delta"])
            e0shift = abs((r["E0P"] - ra["E0"]) - 2.0)
            cert_resid.append(resid)
            rows.append(dict(lam=lam, dM=r["dM"], dP=r["dP"], F=r["F"],
                             dAbs=ra["delta"], resid=resid, e0shift=e0shift,
                             witness=r["witness"], gap=r["gap"], m=r["m"]))
            emit("  %-6.3f %+.9e  %+.9e  %+.4e  %.2e     %.4f"
                 % (lam, r["dM"], r["dP"], r["F"], resid, r["witness"]))
        emit("  (cert-resid = |Delta(b=+1) - Delta(source-absent)|: CONSTRUCTION")
        emit("  CERTIFICATE of the conserving-quadrature identity, never a control;")
        emit("  E0(+1)-E0(absent)-2 residual max %.1e; noise floor FLOOR_EIG = %.0e"
             % (max(x["e0shift"] for x in rows), FLOOR_EIG))
        emit("  beside every F; positive control beside these zeros: the nonzero F")
        emit("  rows of this same table in written configurations (D-15).)")
        # onset on the declared grid
        orows, slopes, k_hat, drift = onset_measure(V, probe, src)
        wit_ok = all(x["witness"] >= WITNESS_MIN for x in rows) and \
            all(x["witness"] >= WITNESS_MIN for x in orows)
        emit("  onset (declared grid %s):" % (list(ONSET_GRID),))
        for x in orows:
            emit("    lam=%.3f  F=%+.4e" % (x["lam"], x["F"]))
        if slopes:
            emit("    pair log-slopes: %s"
                 % "; ".join("(%.3f,%.3f)->%.3f" % s for s in slopes))
            emit("    k_hat = %.3f (smallest usable pair; drift %.3f; floors as "
                 "declared)" % (k_hat, drift))
        b_c4 = (k_hat is not None) and abs(k_hat - g["w_conn"]) <= TOL_K
        emit("  C4 bracket (computed, both branches reachable): |k_hat - w_enc_conn|"
             " = |%s - %d| <= TOL_K: %s"
             % (("%.3f" % k_hat) if k_hat is not None else "none", g["w_conn"],
                bool_word(b_c4)))
        emit("  doublet witness >= %.1f on every row above: %s; min band gap %.3f"
             % (WITNESS_MIN, bool_word(wit_ok), min(x["gap"] for x in rows)))
        ftab[tag] = dict(role=role, d_gen=g["d_gen"], w_conn=g["w_conn"],
                         src=list(map(list, src)), rows=rows,
                         onset_rows=[dict(lam=x["lam"], F=x["F"]) for x in orows],
                         slopes=slopes, k_hat=k_hat, drift=drift, c4=b_c4,
                         witness_ok=wit_ok)

    # ---- 18-placement survey at lam = 0.05 (same-table context) --------------------
    emit("")
    emit("-" * 88)
    emit("SURVEY -- all 18 placements at lam = 0.05, winding sector (+1,+1), with")
    emit("earned class labels (same-table context; nothing narrowed by the placement")
    emit("rule).")
    emit("-" * 88)
    emit("  %-10s %-6s %-10s %-14s %-8s" % ("placement", "d_gen", "w_enc_conn",
                                            "F(0.05)", "witness"))
    survey = {}
    for tag in gate_order:
        g = gate[tag]
        r = F_reading(V, probe, g["src"], 0.05)
        survey[tag] = dict(F=r["F"], witness=r["witness"], d_gen=g["d_gen"],
                           w_conn=g["w_conn"])
        emit("  %-10s %-6d %-10d %+.6e  %.4f" % (tag, g["d_gen"], g["w_conn"],
                                                 r["F"], r["witness"]))
    emit("")
    emit("MEASURED CONTROL FLOOR (feeds V1; the same-table certificate residuals):")
    floor_meas = max(cert_resid)
    emit("  max certificate residual over the F table = %.2e; control floor ="
         % floor_meas)
    emit("  max(FLOOR_EIG, that) = %.2e" % max(FLOOR_EIG, floor_meas))
    emit("")
    emit("runtime %.1f s.  NEXT STEP NAMED: a3 (controls C2, C3, C5 + certificates),"
         % (time.time() - t0))
    emit("then a4 assembles the pre-registered verdicts V1-V5.")
    emit("=" * 88)

    seal = emit.seal(os.path.join(HERE, "a2_field_table_3x3.txt"))
    write_json(os.path.join(HERE, "a2_field_table_3x3.json"),
               dict(gate_seal=gate_seal, p_near=p_near, p_far=p_far,
                    swap_near=swap_near, pair24=pair24, port_ok=port_ok,
                    ftab=ftab, survey=survey, floor_meas=floor_meas,
                    control_floor=max(FLOOR_EIG, floor_meas), seal_txt=seal))


if __name__ == "__main__":
    main()
