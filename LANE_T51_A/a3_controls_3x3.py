"""LANE_T51_A / a3 -- THE CONTROL BATTERY (C2, C3, C5) AND THE CERTIFICATES.

All controls two-way (D-8), in the same tables as the readings they control; the
construction certificates are labeled certificates and never counted as controls
(D-15).  Placements and earned quantities come from a2's sealed sidecar and the
seal-verified gate table.

C2 -- GAMMA-EQUIVALENT PLACEMENT SWAP at equal earned separation.  A swap pair is two
placements in the SAME earned class (equal d_gen and equal connected w_enc).  Verdict
per pair (declared before the run): both members' onsets bracket the class's
w_enc_conn (|k_hat - w_enc_conn| <= TOL_K) AND the two onsets agree with each other
(|k1 - k2| <= TOL_SWAP).  F magnitudes are reported side by side as data (on a
3-column torus the magnitude carries winding furniture; the earned-geometry reading
is the onset).  COMPUTED FACT reported plainly: the largest earned class is a
singleton, so no swap partner exists at the far separation; the far-side swap runs at
the venue's largest multi-member class, (d_gen=2, w_enc_conn=4).

C3 -- THE WINDING-SECTOR SWEEP (mandatory before any sign of F is attributed to
content).  F is computed in EVERY winding sector (w1, w2) of the sealed zbar1/zbar2
labels, at both commissioned placements, at every commissioned lam.  PRE-REGISTERED
SEPARATION CRITERION (declared here, before results): the sweep SEPARATES content
sign from winding sign at a placement iff the sign pattern of F across the four
winding sectors is EXACTLY a winding character times a fixed residual sign --
sign(F(w1,w2)) = s0 * (-1)^(c1*w1 + c2*w2) for one class c = (c1,c2) and s0 in
{+1,-1} -- with the SAME (c, s0) at every commissioned lam at which all four sectors
are usable (|F| >= FLOOR_USE), and at least one such lam.  Then c is the winding
class carrying the reading and s0 is the content-attributed sign.  SCOPE, stated
with the verdict: for c = (0,0) the attribution is convention-free; for any other c
it is licensed ONLY relative to the sealed zbar1/zbar2 winding representatives (the
gate's computed exhibit: enclosure parity of winding-odd strings is
representative-dependent).  Both branches reachable: a non-character sign pattern, a
lam-dependent pattern, or unusable sectors all refuse the license.

C5 -- BACK-ACTION.  The source's written/unwritten sector-energy differential
W = E0(b=-1) - E0(b=+1), computed with the probe present (star holes) and absent
(full star set), same sectors otherwise.  BA = |W(present) - W(absent)| must sit
below the declared TOL_BACK at both commissioned placements at every commissioned
lam.  Two-way: a stronger probe or coupling fails it.

CERTIFICATES (labeled; never counted as controls):
  - conserving quadrature: every single-edge Z of the mediator commutes with both
    source plaquette operators (symplectic pairing identically 0) -- the source value
    is an exact quantum number; zero back-action on the WRITTEN VALUE by construction.
  - X-quadrature probe quiescence: every single-edge X commutes with both probe star
    operators -- the probe's record bit is conserved by the complementary quadrature
    (commutation identity; the critiques' relabeling repair is binding).
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from t51a_lib import (Venue, F_reading, onset_measure, sector_spectrum_g, Emitter,
                      write_json, load_json_verified, parse_gate_table, bool_word,
                      sp_pair, COMM_LAMS, FLOOR_EIG, FLOOR_USE, TOL_K, TOL_SWAP,
                      TOL_BACK, WITNESS_MIN, ONSET_GRID)

emit = Emitter()
CHARS = ((0, 0), (0, 1), (1, 0), (1, 1))


def char_fit(signs):
    """signs: dict (w1,w2) -> +1/-1.  Return (c, s0) if the pattern is exactly a
       winding character times a fixed sign, else None."""
    for c in CHARS:
        s0 = signs[(0, 0)]  # character value at (0,0) is +1
        if all(signs[w] == s0 * (-1) ** (c[0] * w[0] + c[1] * w[1]) for w in CHARS):
            return c, s0
    return None


def main():
    t0 = time.time()
    emit("=" * 88)
    emit("LANE_T51_A / a3 -- CONTROLS C2, C3, C5 AND THE CONSTRUCTION CERTIFICATES")
    emit("date: 2026-08-21   commissioned in FIELD_INSTRUMENT_V001.md section 1")
    emit("=" * 88)
    emit("")
    emit("DECLARED BEFORE USE: TOL_SWAP = %.2f (C2 onset agreement); TOL_BACK = %.0e"
         % (TOL_SWAP, TOL_BACK))
    emit("(C5, energy units of the unit stabilizer term; basis: three orders below")
    emit("the sector gap ~2 and one order above the largest commissioned far-side")
    emit("reading; scale surveyed in the scratchpad BEFORE this sealed run, logged in")
    emit("D24_AUDIT.txt); C3 separation criterion as pre-registered in this script's")
    emit("header; floors FLOOR_EIG = %.0e, FLOOR_USE = %.0e; WITNESS_MIN = %.1f."
         % (FLOOR_EIG, FLOOR_USE, WITNESS_MIN))
    emit("")

    gate, gate_seal = parse_gate_table()
    a2 = load_json_verified(os.path.join(HERE, "a2_field_table_3x3.json"))
    assert a2["gate_seal"] == gate_seal
    p_near, p_far = a2["p_near"], a2["p_far"]
    swap_near, pair24 = a2["swap_near"], a2["pair24"]
    emit("INPUTS: gate seal %s... verified; a2 sidecar seal verified." % gate_seal[:16])
    emit("Commissioned placements: near %s (d_gen=%d, w_enc_conn=%d), far %s"
         % (p_near, gate[p_near]["d_gen"], gate[p_near]["w_conn"], p_far))
    emit("(d_gen=%d, w_enc_conn=%d)." % (gate[p_far]["d_gen"], gate[p_far]["w_conn"]))
    emit("")

    V = Venue(3, 3)
    probe = ((0, 0), (1, 0))
    T = V.T
    n = T.n

    # ---- CERTIFICATES ---------------------------------------------------------------
    emit("-" * 88)
    emit("CONSTRUCTION CERTIFICATES (labeled; computed; never counted as controls)")
    emit("-" * 88)
    src_ops = {tag: [T.plaq(*p) for p in gate[tag]["src"]]
               for tag in (p_near, p_far)}
    probe_ops = [T.star(*v) for v in probe]
    cz = all(sp_pair(1 << (n + e), B, n) == 0
             for e in range(n) for tag in src_ops for B in src_ops[tag])
    cx = all(sp_pair(1 << e, A, n) == 0 for e in range(n) for A in probe_ops)
    emit("  [V, B_src] = 0: symplectic pairing of every single-edge Z (18 edges)")
    emit("  with both source plaquette operators, both commissioned placements,")
    emit("  identically zero: %s  -- the source value is an exact sector label;"
         % bool_word(cz))
    emit("  zero back-action on the written VALUE by construction.")
    emit("  X-quadrature probe quiescence: symplectic pairing of every single-edge X")
    emit("  with both probe star operators identically zero: %s  -- the probe bit"
         % bool_word(cx))
    emit("  is conserved by the complementary quadrature (commutation identity).")
    emit("")

    # ---- C2 -------------------------------------------------------------------------
    emit("-" * 88)
    emit("C2 -- GAMMA-EQUIVALENT PLACEMENT SWAP AT EQUAL EARNED SEPARATION (two-way)")
    emit("-" * 88)
    singleton = len([t for t, g in gate.items()
                     if (g["d_gen"], g["w_conn"]) ==
                     (gate[p_far]["d_gen"], gate[p_far]["w_conn"])]) == 1
    emit("  COMPUTED FACT: the far class (d_gen=%d, w_enc_conn=%d) has exactly one"
         % (gate[p_far]["d_gen"], gate[p_far]["w_conn"]))
    emit("  placement (singleton: %s) -- no swap partner exists at the far"
         % bool_word(singleton))
    emit("  separation; the far-side swap therefore runs at the largest multi-member")
    emit("  class, (d_gen=2, w_enc_conn=4): %s <-> %s." % tuple(pair24))
    swap_pairs = [("near-class swap", p_near, swap_near),
                  ("far-side swap", pair24[0], pair24[1])]
    c2 = {}
    for label, ta, tb in swap_pairs:
        ga, gb = gate[ta], gate[tb]
        assert (ga["d_gen"], ga["w_conn"]) == (gb["d_gen"], gb["w_conn"])
        res = {}
        for tag in (ta, tb):
            orows, slopes, k_hat, drift = onset_measure(V, probe, gate[tag]["src"])
            F05 = F_reading(V, probe, gate[tag]["src"], 0.05)
            res[tag] = dict(k_hat=k_hat, drift=drift, F05=F05["F"],
                            witness=F05["witness"])
        ka, kb = res[ta]["k_hat"], res[tb]["k_hat"]
        w = ga["w_conn"]
        ok_class = (ka is not None and kb is not None and
                    abs(ka - w) <= TOL_K and abs(kb - w) <= TOL_K)
        ok_agree = ka is not None and kb is not None and abs(ka - kb) <= TOL_SWAP
        verdict = ok_class and ok_agree
        emit("")
        emit("  %s: %s <-> %s  [class (d_gen=%d, w_enc_conn=%d)]"
             % (label, ta, tb, ga["d_gen"], w))
        for tag in (ta, tb):
            emit("    %s: k_hat = %s (drift %s)  F(0.05) = %+.6e  wit %.4f"
                 % (tag,
                    ("%.3f" % res[tag]["k_hat"])
                    if res[tag]["k_hat"] is not None else "none",
                    ("%.3f" % res[tag]["drift"])
                    if res[tag]["drift"] is not None else "n/a",
                    res[tag]["F05"], res[tag]["witness"]))
        dk = abs(ka - kb) if (ka is not None and kb is not None) else None
        emit("    both onsets bracket w_enc_conn = %d (TOL_K): %s;  |k_a - k_b| = "
             "%s <= TOL_SWAP: %s"
             % (w, bool_word(ok_class),
                ("%.3f" % dk) if dk is not None else "n/a", bool_word(ok_agree)))
        emit("    F-magnitude ratio (data, not a verdict -- magnitude carries")
        emit("    winding furniture on a 3-column torus): %.3f"
             % (abs(res[ta]["F05"]) / abs(res[tb]["F05"])))
        emit("    SWAP VERDICT (%s) = %s" % (label, bool_word(verdict)))
        c2[label] = dict(pair=(ta, tb), ok_class=ok_class, ok_agree=ok_agree,
                         verdict=verdict, res=res)
    c2_all = all(v["verdict"] for v in c2.values())
    emit("")
    emit("  C2 = both swap verdicts = %s" % bool_word(c2_all))
    emit("")

    # ---- C3 -------------------------------------------------------------------------
    emit("-" * 88)
    emit("C3 -- THE WINDING-SECTOR SWEEP (mandatory before any sign attribution)")
    emit("-" * 88)
    emit("  F in every winding sector (w1, w2) of the sealed zbar1/zbar2 labels")
    emit("  [parity p <-> eigenvalue (-1)^p], both commissioned placements, every")
    emit("  commissioned lam.  Any sector with |F| < FLOOR_USE would be excluded")
    emit("  from the pattern (indeterminate) -- exclusions are printed.")
    c3 = {}
    for tag in (p_near, p_far):
        src = gate[tag]["src"]
        emit("")
        emit("  PLACEMENT %s (d_gen=%d, w_enc_conn=%d):" % (tag, gate[tag]["d_gen"],
                                                            gate[tag]["w_conn"]))
        emit("    %-6s %-15s %-15s %-15s %-15s %-9s"
             % ("lam", "F(+ +)", "F(+ -)", "F(- +)", "F(- -)", "wit_min"))
        fits = []
        rows_out = []
        for lam in COMM_LAMS:
            Fs, wits = {}, []
            for w in CHARS:
                r = F_reading(V, probe, src, lam, w[0], w[1])
                Fs[w] = r["F"]
                wits.append(r["witness"])
            rows_out.append(dict(lam=lam, F=dict((str(k), v)
                                                 for k, v in Fs.items()),
                                 wit=min(wits)))
            emit("    %-6.3f %+.6e  %+.6e  %+.6e  %+.6e  %.4f"
                 % (lam, Fs[(0, 0)], Fs[(0, 1)], Fs[(1, 0)], Fs[(1, 1)],
                    min(wits)))
            usable = all(abs(Fs[w]) >= FLOOR_USE for w in CHARS)
            if not usable:
                below = [w for w in CHARS if abs(Fs[w]) < FLOOR_USE]
                emit("      lam=%.3f EXCLUDED from the pattern: |F| < FLOOR_USE in "
                     "sectors %s" % (lam, below))
                continue
            signs = {w: (1 if Fs[w] > 0 else -1) for w in CHARS}
            fits.append((lam, char_fit(signs)))
        for lam, fit in fits:
            emit("      lam=%.3f sign pattern fit: %s"
                 % (lam, ("character c=%s, residual sign s0=%+d" % fit)
                    if fit else "NO character fits"))
        ok = (len(fits) > 0 and all(f[1] is not None for f in fits)
              and len({f[1] for f in fits}) == 1)
        c_found = fits[0][1] if ok else None
        emit("    SEPARATION VERDICT at %s (pre-registered criterion): %s"
             % (tag, bool_word(ok)))
        if ok:
            scope = ("convention-free (winding-even class)" if c_found[0] == (0, 0)
                     else "licensed ONLY relative to the sealed zbar1/zbar2 "
                          "representatives (winding-odd class; the gate's exhibit "
                          "makes this scope binding)")
            emit("    -> the reading's winding class is c=%s, content-attributed "
                 "sign s0=%+d;" % (c_found[0], c_found[1]))
            emit("       sign attribution scope: %s." % scope)
            emit("    labeled cross-reference (data, not on the decision path): the")
            emit("    gate's connected-minimum winding classes at this placement =")
            emit("    %s; dynamical class matches the static minimum: %s"
                 % (gate[tag]["conn_windings"],
                    bool_word(list(c_found[0]) in
                              [list(x) for x in gate[tag]["conn_windings"]])))
        c3[tag] = dict(rows=rows_out, ok=ok,
                       c=list(c_found[0]) if ok else None,
                       s0=c_found[1] if ok else None)
    c3_all = all(c3[tag]["ok"] for tag in c3)
    emit("")
    emit("  C3 = sweep separates content sign from winding sign at both commissioned")
    emit("  placements = %s" % bool_word(c3_all))
    emit("")

    # ---- C5 -------------------------------------------------------------------------
    emit("-" * 88)
    emit("C5 -- BACK-ACTION (two-way, declared tolerance)")
    emit("-" * 88)
    emit("  W = E0(b=-1) - E0(b=+1) in the source's sectors (winding (+1,+1)),")
    emit("  probe present (star holes) vs absent (full star set); BA = |difference|.")
    c5 = {}
    for tag in (p_near, p_far):
        src = gate[tag]["src"]
        emit("")
        emit("  PLACEMENT %s:" % tag)
        emit("    %-6s %-16s %-16s %-12s %-10s"
             % ("lam", "W(probe present)", "W(probe absent)", "BA", "<=TOL_BACK"))
        rows_out = []
        for lam in COMM_LAMS:
            Wp = (sector_spectrum_g(V, probe, src, -1, lam)["E0"]
                  - sector_spectrum_g(V, probe, src, +1, lam)["E0"])
            Wa = (sector_spectrum_g(V, probe, src, -1, lam,
                                    probe_present=False)["E0"]
                  - sector_spectrum_g(V, probe, src, +1, lam,
                                      probe_present=False)["E0"])
            ba = abs(Wp - Wa)
            ok = ba <= TOL_BACK
            rows_out.append(dict(lam=lam, Wp=Wp, Wa=Wa, BA=ba, ok=ok))
            emit("    %-6.3f %+.6e    %+.6e    %.3e   %s"
                 % (lam, Wp, Wa, ba, bool_word(ok)))
        allok = all(r["ok"] for r in rows_out)
        emit("    BACK-ACTION VERDICT at %s = %s" % (tag, bool_word(allok)))
        c5[tag] = dict(rows=rows_out, ok=allok)
    c5_all = all(c5[tag]["ok"] for tag in c5)
    emit("")
    emit("  C5 = back-action below TOL_BACK = %.0e at both commissioned placements,"
         % TOL_BACK)
    emit("  every commissioned lam = %s" % bool_word(c5_all))
    emit("")
    emit("runtime %.1f s.  NEXT STEP NAMED: a4 assembles the pre-registered verdicts"
         % (time.time() - t0))
    emit("V1-V5 from the sealed sidecars of a1, a2, a3.")
    emit("=" * 88)

    seal = emit.seal(os.path.join(HERE, "a3_controls_3x3.txt"))
    write_json(os.path.join(HERE, "a3_controls_3x3.json"),
               dict(gate_seal=gate_seal, cert_cz=cz, cert_cx=cx,
                    singleton_far=singleton, c2=c2, c2_all=c2_all, c3=c3,
                    c3_all=c3_all, c5=c5, c5_all=c5_all, seal_txt=seal))


if __name__ == "__main__":
    main()
