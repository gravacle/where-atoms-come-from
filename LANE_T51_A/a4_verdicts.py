"""LANE_T51_A / a4 -- THE PRE-REGISTERED DECISION RULE, ASSEMBLED.

V1-V5 exactly as pre-registered in FIELD_INSTRUMENT_V001.md section 1, computed from
the sealed sidecars of a1 (calibration), a2 (the F table), a3 (controls); every
sidecar seal verified before use.  All verdicts are computed booleans with both
branches reachable (D-8); no literal expected value sits on any decision path; the
calibration verdict is a reproduction comparison and is reported beside, not inside,
V1-V5.  Per the principal's binding directive (2026-08-20, quoted in C-92), no
outcome here is a failure against an imported standard: every branch registers as
the surface's own law, and the branch taken names its next step.
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from t51a_lib import (Emitter, load_json_verified, bool_word, COMM_LAMS, FLOOR_EIG,
                      FLOOR_USE, TOL_K, TOL_SWAP, TOL_BACK, TOL_CAL_REL, TOL_PORT_D,
                      TOL_PORT_F, WITNESS_MIN, ONSET_GRID)

emit = Emitter()


def main():
    t0 = time.time()
    a1 = load_json_verified(os.path.join(HERE, "a1_calibration_3x2.json"))
    a2 = load_json_verified(os.path.join(HERE, "a2_field_table_3x3.json"))
    a3 = load_json_verified(os.path.join(HERE, "a3_controls_3x3.json"))
    assert a2["gate_seal"] == a3["gate_seal"]

    emit("=" * 88)
    emit("LANE_T51_A / a4 -- THE MEASUREMENT: PRE-REGISTERED VERDICTS V1-V5")
    emit("date: 2026-08-21   (FIELD_INSTRUMENT_V001.md section 1; C-92; gate seal")
    emit("%s... verified; sidecar seals of a1/a2/a3 verified)" % a2["gate_seal"][:16])
    emit("=" * 88)
    emit("")
    p_near, p_far = a2["p_near"], a2["p_far"]
    near = a2["ftab"][p_near]
    far = a2["ftab"][p_far]

    emit("CALIBRATION (a1, reproduction, stated first as commissioned): CAL = %s"
         % bool_word(a1["CAL"]))
    emit("  contact: quoted %+.2e, matched %s rel diff %.1e; far: quoted %+.2e,"
         % (a1["targets"]["F_contact"], a1["match_contact"]["tag"],
            a1["match_contact"]["rel"], a1["targets"]["F_far"]))
    emit("  matched %s rel diff %.1e (TOL_CAL_REL = %.0e).  PORT CERTIFICATE (a2,"
         % (a1["match_far"]["tag"], a1["match_far"]["rel"], TOL_CAL_REL))
    emit("  lam=0.064 vs the gate's sealed prints): %s" % bool_word(a2["port_ok"]))
    emit("")

    emit("THE COMMISSIONED READINGS (winding sector (+1,+1); full tables in a2/a3):")
    emit("  %-10s %-24s %-6s %-14s" % ("placement", "earned (d_gen, w_conn)", "lam",
                                       "F"))
    for tag, blk in ((p_near, near), (p_far, far)):
        for r in blk["rows"]:
            emit("  %-10s (d_gen=%d, w_enc_conn=%d)    %-6.3f %+.6e"
                 % (tag, blk["d_gen"], blk["w_conn"], r["lam"], r["F"]))
    emit("")

    # ---- V1 -------------------------------------------------------------------------
    cf = a2["control_floor"]
    minF_far = min(abs(r["F"]) for r in far["rows"])
    v1 = minF_far > cf
    emit("V1  |F| exceeds the same-table measured control floor beyond contact at")
    emit("    the larger separation: min over commissioned lams |F(%s)| = %.3e >"
         % (p_far, minF_far))
    emit("    control floor %.3e (= max(FLOOR_EIG, max certificate residual %.2e))"
         % (cf, a2["floor_meas"]))
    emit("    V1 = %s" % bool_word(v1))
    emit("")

    # ---- V2 -------------------------------------------------------------------------
    v2 = a3["c2_all"]
    emit("V2  the reading follows earned geometry under the placement swap (C2,")
    emit("    both Gamma-equivalent pairs; TOL_K = %.2f, TOL_SWAP = %.2f):"
         % (TOL_K, TOL_SWAP))
    for label, blk in a3["c2"].items():
        emit("    %s %s <-> %s: %s" % (label, blk["pair"][0], blk["pair"][1],
                                       bool_word(blk["verdict"])))
    emit("    [computed fact: the far class is a singleton -- no swap partner exists")
    emit("    at (d_gen=2, w_enc_conn=5); the far-side swap ran at (2, 4).]")
    emit("    V2 = %s" % bool_word(v2))
    emit("")

    # ---- V3 -------------------------------------------------------------------------
    v3 = near["c4"] and far["c4"]
    emit("V3  the onset-order bracket contains CONNECTED w_enc at both placements")
    emit("    (scored per placement against the sealed gate table, TOL_K = %.2f):"
         % TOL_K)
    for tag, blk in ((p_near, near), (p_far, far)):
        emit("    %s: k_hat = %s (drift %s) vs w_enc_conn = %d: %s"
             % (tag,
                ("%.3f" % blk["k_hat"]) if blk["k_hat"] is not None else "none",
                ("%.3f" % blk["drift"]) if blk["drift"] is not None else "n/a",
                blk["w_conn"], bool_word(blk["c4"])))
    emit("    V3 = %s" % bool_word(v3))
    emit("")

    # ---- V4 -------------------------------------------------------------------------
    v4 = a3["c5_all"]
    emit("V4  back-action below the declared tolerance (C5, TOL_BACK = %.0e):"
         % TOL_BACK)
    for tag in (p_near, p_far):
        mx = max(r["BA"] for r in a3["c5"][tag]["rows"])
        emit("    %s: max BA over commissioned lams = %.3e: %s"
             % (tag, mx, bool_word(a3["c5"][tag]["ok"])))
    emit("    V4 = %s" % bool_word(v4))
    emit("")

    # ---- V5 -------------------------------------------------------------------------
    v5 = a3["c3_all"]
    emit("V5  sign attribution licensed only if the winding-sector sweep separates")
    emit("    content sign from winding sign (C3, pre-registered character")
    emit("    criterion):")
    for tag in (p_near, p_far):
        blk = a3["c3"][tag]
        if blk["ok"]:
            scope = ("convention-free" if blk["c"] == [0, 0] else
                     "relative to the sealed zbar1/zbar2 convention ONLY")
            emit("    %s: SEPARATED -- winding class c=%s, content sign s0=%+d "
                 "(%s)" % (tag, tuple(blk["c"]), blk["s0"], scope))
        else:
            emit("    %s: NOT SEPARATED -- no stable character factorization"
                 % tag)
    emit("    V5 = %s" % bool_word(v5))
    emit("")

    # ---- assembly -------------------------------------------------------------------
    allpass = v1 and v2 and v3 and v4 and v5
    emit("=" * 88)
    emit("ALL PASS = V1 AND V2 AND V3 AND V4 AND V5 = %s" % bool_word(allpass))
    emit("=" * 88)
    if allpass:
        emit("")
        emit("BRANCH TAKEN (as pre-registered): the record surface has a field side")
        emit("in the only sense it has earned -- a test record's own law is modified")
        emit("at a place where the written content is not, by an amount set by the")
        emit("content's value and earned separation, through the declared, priced,")
        emit("content-blind mediator (C-80's computed division of labor).  Scope")
        emit("carried with the claim: sign attribution is licensed per placement")
        emit("exactly with the scope C3 printed beside it (convention-free only for")
        emit("the winding-even class; otherwise relative to the sealed zbar")
        emit("convention); the reading is Z_2-valued in the source and nothing here")
        emit("says how it composes -- MEDIATION, COMPOSITION, SHAPE, O-39's full")
        emit("ask, and empirical contact all remain exactly as carried in")
        emit("FIELD_INSTRUMENT_V001.md section 4.")
        emit("")
        emit("NEXT STEP NAMED (the commission's own order): (a) the two-source")
        emit("composition at 4x4 FIRST -- it measures what composition IS on this")
        emit("surface, and every outcome (superposition, mod-2 saturation,")
        emit("screening) registers as the surface's own law; then (b) the Tier-2")
        emit("signed kernel over the connected/dynamical string class; then (c) the")
        emit("Tier-3 grounded tau readout.  The surviving operator combination")
        emit("(source Zbar x probe writer amplitude) feeds O-53 its field-side")
        emit("candidate.")
    else:
        emit("")
        emit("BRANCH TAKEN (as pre-registered): one or more verdicts FALSE -- per")
        emit("the decision rule, sweep the priced mediator family next (quadratures,")
        emit("staggered signs, weight-2 local terms) before any negative closure is")
        emit("registered; the null, if it stands, reads TWO ways (no field present,")
        emit("or field-side concept differently shaped) and the register entry must")
        emit("carry both.  Every outcome registers as the surface's own law.")
    emit("")
    emit("TOLERANCE REGISTRY (every tolerance used in this lane, declared before")
    emit("use in the script that used it):")
    emit("  FLOOR_EIG   = %.0e  absolute eigenvalue-difference noise floor" % FLOOR_EIG)
    emit("  FLOOR_USE   = %.0e  onset-pair usability floor (= 100 x FLOOR_EIG)"
         % FLOOR_USE)
    emit("  TOL_K       = %.2f  onset-vs-integer bracket half-width" % TOL_K)
    emit("  TOL_SWAP    = %.2f  C2 onset agreement inside a swap pair" % TOL_SWAP)
    emit("  TOL_BACK    = %.0e  C5 back-action bound (energy units of the unit"
         % TOL_BACK)
    emit("                      stabilizer term)")
    emit("  TOL_CAL_REL = %.0e  a1 calibration reproduction (3-significant-figure"
         % TOL_CAL_REL)
    emit("                      targets)")
    emit("  TOL_PORT_D  = %.0e  port certificate, Delta vs gate prints" % TOL_PORT_D)
    emit("  TOL_PORT_F  = %.0e  port certificate, F vs gate prints" % TOL_PORT_F)
    emit("  WITNESS_MIN = %.1f  doublet-identity witness floor" % WITNESS_MIN)
    emit("  eigenvalue solver: dense sector-exact eigh (LAPACK), sector dims 32/256")
    emit("  (printed and checked against the declared 2^11 bound before every run);")
    emit("  onset grid %s; commissioned mediator grid %s." % (list(ONSET_GRID),
                                                              list(COMM_LAMS)))
    emit("")
    emit("runtime %.1f s (assembly only; all numbers from sealed sidecars)."
         % (time.time() - t0))
    emit("=" * 88)
    emit.seal(os.path.join(HERE, "a4_verdicts.txt"))


if __name__ == "__main__":
    main()
