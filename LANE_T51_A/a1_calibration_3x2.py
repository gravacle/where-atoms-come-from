"""LANE_T51_A / a1 -- THE 3x2 CALIBRATION (runs FIRST, before any new number).

The commission (FIELD_INSTRUMENT_V001.md section 1): reproduce the design
exploration's numbers -- F(contact) = -1.99e-3, F(far) = +1.94e-4 at lam = 0.05 --
as a calibration CHECK, stated as reproduction, before any new number.

METHOD.  The quoted values are PRIOR DATA: they are parsed at runtime from the sealed
design document LANE_T51_IDENT/T51_DESIGN_the-second-lump.json (nothing hardcoded),
labeled as reproduction targets, and never enter any physics verdict of this lane
(D-8: the calibration verdict is a reproduction comparison, the gate's prior-art
pattern).  The exploration did not record its placement coordinates, so the check is
exhaustive: with the probe stars at (0,0),(1,0) (the same probe every Second Lump
computation uses), ALL adjacent source plaquette pairs the 3x2 torus affords are
computed at lam = 0.05, and the reproduction verdict is: some placement matches the
quoted contact value and some placement matches the quoted far value, each within the
declared TOL_CAL_REL with matching sign.  Earned separations (connected w_enc; d_gen
descriptive) are computed per placement by the gate's own geometry_block machinery
and stated beside every row; the matched placements' onset orders on the declared
grid are reported beside the exploration's quoted onsets as labeled cross-references.

CAVEAT CARRIED (Second Lump critique, measurability): 3x2 is a degenerate proxy --
at Ly = 2 the vertical plaquette pairs share TWO edges and wrap classes are cheap.
The calibration certifies the port; 3x3 is the measurement venue.
"""
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from t51a_lib import (Venue, geometry_block, F_reading, onset_measure, Emitter,
                      write_json, bool_word, TOL_CAL_REL, FLOOR_EIG, FLOOR_USE,
                      TOL_K, WITNESS_MIN, ONSET_GRID, BOUND_LOG2)

REPO = os.path.dirname(HERE)
DESIGN_JSON = os.path.join(REPO, "LANE_T51_IDENT", "T51_DESIGN_the-second-lump.json")

emit = Emitter()


def parse_targets():
    """Parse the design exploration's quoted numbers at runtime (prior data)."""
    txt = open(DESIGN_JSON).read()
    m = re.search(r"giving F = (-[\d.]+e-\d+) at lam=0\.05 contact and "
                  r"\+([\d.]+e-\d+) at the far placement", txt)
    F_contact = float(m.group(1))
    F_far = float(m.group(2))
    m2 = re.search(r"computed ([\d.]+) and ([\d.]+) against w_enc = 3 and 4", txt)
    onsets = (float(m2.group(1)), float(m2.group(2)))
    return F_contact, F_far, onsets


def main():
    t0 = time.time()
    emit("=" * 88)
    emit("LANE_T51_A / a1 -- THE 3x2 CALIBRATION (reproduction, before any new number)")
    emit("date: 2026-08-21   commissioned in FIELD_INSTRUMENT_V001.md section 1")
    emit("=" * 88)
    emit("")
    emit("DECLARED BEFORE USE: TOL_CAL_REL = %.0e (quoted targets carry 3 significant"
         % TOL_CAL_REL)
    emit("figures; half-ulp quantization <= 2.5e-3 relative, 2x taken); onset grid =")
    emit("%s with FLOOR_USE = %.0e and FLOOR_EIG = %.0e (the gate's declared"
         % (list(ONSET_GRID), FLOOR_USE, FLOOR_EIG))
    emit("estimator, unchanged); WITNESS_MIN = %.1f; TOL_K = %.2f." % (WITNESS_MIN,
                                                                       TOL_K))
    emit("")

    F_contact_q, F_far_q, onsets_q = parse_targets()
    emit("REPRODUCTION TARGETS (parsed at runtime from the sealed design file")
    emit("%s --" % os.path.relpath(DESIGN_JSON, REPO))
    emit("prior data, labeled; never on any physics decision path of this lane):")
    emit("  F(contact, lam=0.05) quoted = %+.2e" % F_contact_q)
    emit("  F(far,     lam=0.05) quoted = %+.2e" % F_far_q)
    emit("  onset orders quoted         = %.1f and %.1f (against w_enc 3 and 4)"
         % onsets_q)
    emit("")

    Lx, Ly = 3, 2
    V = Venue(Lx, Ly)
    n = 2 * Lx * Ly
    emit("VENUE: toric %dx%d, %d edges; full Hilbert space 2^%d = %d; sector-exact"
         % (Lx, Ly, n, n, 1 << n))
    emit("blocks via the diagonal conserved algebra (plaquettes, source values,")
    emit("Z-winding pair): dimension %d = 2^%d <= declared bound 2^%d -- checked"
         % (1 << (n - (Lx * Ly - 1) - 2), n - (Lx * Ly - 1) - 2, BOUND_LOG2))
    emit("before any eigh call (asserted again inside every sector build).")
    emit("")

    probe = ((0, 0), (1, 0))
    pairs = []
    for y in range(Ly):
        for x in range(Lx):
            pairs.append(("%d%dH" % (x, y), ((x, y), ((x + 1) % Lx, y))))
    for x in range(Lx):
        pairs.append(("%d0V" % x, ((x, 0), (x, 1))))
    emit("ALL %d adjacent source plaquette pairs at 3x2 (probe stars %s,%s), each"
         % (len(pairs), probe[0], probe[1]))
    emit("with earned separation from the gate's geometry machinery (connected w_enc")
    emit("earned; d_gen descriptive -- the mixed-type demotion is binding); vertical")
    emit("pairs at Ly=2 are doubly adjacent (share two edges) -- the degenerate-proxy")
    emit("caveat above.")
    emit("")
    rows = []
    for tag, src in pairs:
        g = geometry_block(Lx, Ly, probe, src)
        r = F_reading(V, probe, src, 0.05)
        rows.append(dict(tag=tag, src=src, w_conn=g["w_conn"], w_old=g["w_old"],
                         d_gen=g["d_gen"], F=r["F"], witness=r["witness"],
                         gap=r["gap"], m=r["m"], dM=r["dM"], dP=r["dP"]))
    emit("  %-5s %-18s %-7s %-6s %-6s %-14s %-8s %-6s"
         % ("tag", "src (constr.)", "w_conn", "w_old", "d_gen", "F(0.05)", "wit",
            "dim"))
    for r in rows:
        emit("  %-5s %-18s %-7d %-6d %-6d %+.6e  %.4f  %d"
             % (r["tag"], str(r["src"]), r["w_conn"], r["w_old"], r["d_gen"],
                r["F"], r["witness"], r["m"]))
    emit("  (noise floor FLOOR_EIG = %.0e beside every F above; doublet witness"
         % FLOOR_EIG)
    emit("  >= %.1f on every row: %s)" % (WITNESS_MIN,
         bool_word(all(r["witness"] >= WITNESS_MIN for r in rows))))
    emit("")

    # ---- the reproduction comparison ------------------------------------------------
    def match(rowlist, target):
        best = min(rowlist, key=lambda r: abs(r["F"] - target))
        rel = abs(best["F"] - target) / abs(target)
        ok = rel <= TOL_CAL_REL and (best["F"] * target > 0)
        return best, rel, ok

    mc, rel_c, ok_c = match(rows, F_contact_q)
    mf, rel_f, ok_f = match(rows, F_far_q)
    emit("CALIBRATION REPRODUCTION (comparison against prior data, labeled):")
    emit("  contact target %+.2e : best match %s (w_enc_conn=%d) F=%+.6e,"
         % (F_contact_q, mc["tag"], mc["w_conn"], mc["F"]))
    emit("    rel diff %.1e <= TOL_CAL_REL, sign match: %s" % (rel_c,
                                                               bool_word(ok_c)))
    emit("  far target     %+.2e : best match %s (w_enc_conn=%d) F=%+.6e,"
         % (F_far_q, mf["tag"], mf["w_conn"], mf["F"]))
    emit("    rel diff %.1e <= TOL_CAL_REL, sign match: %s" % (rel_f,
                                                               bool_word(ok_f)))
    emit("")

    # onset orders at the matched placements, beside the quoted onsets (cross-ref)
    onset_out = {}
    for r, quoted in ((mc, onsets_q[0]), (mf, onsets_q[1])):
        orows, slopes, k_hat, drift = onset_measure(V, probe, dict(pairs)[r["tag"]])
        onset_out[r["tag"]] = dict(k_hat=k_hat, drift=drift)
        emit("  onset at matched %s placement %s: k_hat = %.3f (drift %.3f; grid %s;"
             % ("contact" if r is mc else "far", r["tag"], k_hat, drift,
                list(ONSET_GRID)))
        emit("    floors beside as declared); exploration quoted %.1f -- labeled"
             % quoted)
        emit("    cross-reference, |diff| = %.3f; my earned w_enc_conn here = %d,"
             % (abs(k_hat - quoted), r["w_conn"]))
        emit("    |k_hat - w_enc_conn| <= TOL_K: %s (computed, informative)"
             % bool_word(abs(k_hat - r["w_conn"]) <= TOL_K))
    emit("")

    cal = ok_c and ok_f
    emit("CALIBRATION VERDICT (reproduction; computed booleans, both branches")
    emit("reachable -- a port that drifted would fail either match):")
    emit("  CAL = (contact value reproduced) AND (far value reproduced) = %s"
         % bool_word(cal))
    emit("")
    emit("runtime %.1f s.  NEXT STEP NAMED: a2 (the commissioned F table at 3x3)"
         % (time.time() - t0))
    emit("runs only with this reproduction stated first, as commissioned.")
    emit("=" * 88)

    out_txt = os.path.join(HERE, "a1_calibration_3x2.txt")
    emit_seal = emit.seal(out_txt)
    payload = dict(rows=[{k: v for k, v in r.items() if k != "src"} |
                         {"src": list(map(list, r["src"]))} for r in rows],
                   targets=dict(F_contact=F_contact_q, F_far=F_far_q,
                                onsets=list(onsets_q)),
                   match_contact=dict(tag=mc["tag"], rel=rel_c, ok=ok_c),
                   match_far=dict(tag=mf["tag"], rel=rel_f, ok=ok_f),
                   onsets=onset_out, CAL=cal, seal_txt=emit_seal)
    write_json(os.path.join(HERE, "a1_calibration_3x2.json"), payload)


if __name__ == "__main__":
    main()
