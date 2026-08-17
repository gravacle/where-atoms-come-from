#!/usr/bin/env python3
"""
LANE W10-A — SCRIPT 1.  BUILD B0b AND B4 CONCRETELY AND VALIDATE AGAINST S4'S PUBLISHED ROW.

Everything below is computed from MY incidence matrices.  Nothing is asserted.
Ranks are EXACT over Q (Fraction elimination); class multisets are EXACT (integer counts).
S4's published rows are hard-coded here as the TARGET and every mismatch is printed loudly.
"""
import sys
from fractions import Fraction
from w10a_lib import (Carrier, K1, B1q, B0b, B4, exact_rank, in_column_span,
                      CLASSES, CLASS_NAME, L_S, hnf2)

LOG = []


def out(s=""):
    print(s)
    LOG.append(s)


out("=" * 100)
out("W10-A SCRIPT 1 — B0b AND B4 BUILT AND VALIDATED AGAINST S4:511-590")
out("=" * 100)
out("Exact integer/Fraction arithmetic throughout this script.  No float appears in any check.")
out()

# S4's PUBLISHED rows, transcribed from S4_THE_MEASUREMENT_V001.md:511-516, :538-542, :574-578.
# columns: V E F chi b0 b1 b2 gaugeparams invariants curvature flat
S4_ROW = {
    "B0b": dict(V=9, E=18, F=9, chi=0, b0=1, b1=2, b2=1, gauge=8, inv=10, curv=8, flat=2,
                multiset={(0, 0): 4, (0, 1): 1, (1, 0): 2, (1, 1): 2},
                gF_bounds=True, gC_bounds=False, independent=True),
    "B4":  dict(V=6, E=8, F=4, chi=2, b0=1, b1=1, b2=2, gauge=5, inv=3, curv=2, flat=1,
                multiset={(0, 0): 1, (0, 1): 1, (1, 0): 1, (1, 1): 3},
                gF_bounds=True, gC_bounds=False, independent=True),
    "B1":  dict(V=5, E=6, F=1, chi=0, b0=1, b1=1, b2=0, gauge=4, inv=2, curv=1, flat=1,
                multiset={(0, 0): 0, (0, 1): 2, (1, 0): 2, (1, 1): 1},
                gF_bounds=True, gC_bounds=False, independent=True),
    "B1q": dict(V=7, E=8, F=1, chi=0, b0=1, b1=1, b2=0, gauge=6, inv=2, curv=1, flat=1,
                multiset={(0, 0): 1, (0, 1): 3, (1, 0): 3, (1, 1): 0},
                gF_bounds=True, gC_bounds=False, independent=True),
}

CARRIERS = [("B1", K1()), ("B1q", B1q()), ("B0b", B0b()), ("B4", B4())]

FAILS = []


def check(tag, got, want, label):
    ok = (got == want)
    if not ok:
        FAILS.append(f"{tag}: {label} got {got} want {want}")
    return "OK " if ok else "**MISMATCH**"


for tag, K in CARRIERS:
    R = S4_ROW[tag]
    out("-" * 100)
    out(f"CARRIER {tag}  —  {K.name}")
    out("-" * 100)
    out("VERTICES: " + ", ".join(f"{i}:{n}" for i, n in enumerate(K.vnames)))
    out("EDGES (tail -> head):")
    out("   " + "  ".join(f"e{e}:{K.vnames[t]}->{K.vnames[h]}" for e, (t, h) in enumerate(K.edges)))
    out("FACES (attaching chains, (edge,sign)):")
    for f, cyc in enumerate(K.faces):
        out(f"   F{f}: " + " ".join(f"{'+' if s > 0 else '-'}e{e}" for (e, s) in cyc))
    out(f"gamma_F = " + " ".join(f"{'+' if s > 0 else '-'}e{e}" for (e, s) in K.gF))
    out(f"gamma_C = " + " ".join(f"{'+' if s > 0 else '-'}e{e}" for (e, s) in K.gC))
    out()

    b0, b1, b2, r1, r2 = K.betti()
    chi = K.nV - K.nE + K.nF
    gauge = K.nV - b0
    inv = K.nE - gauge
    curv, flat = r2, b1
    d1d2 = K.d1d2_max()

    out(f"  d1 . d2 max|entry| = {d1d2}        {check(tag, d1d2, 0, 'd1.d2=0')}")
    out(f"  regular (attaching maps injective on the boundary): {K.regular()}")
    out(f"  V={K.nV} E={K.nE} F={K.nF}   rank d1={r1} rank d2={r2}")
    out(f"     V     got {K.nV:3d}  S4 {R['V']:3d}   {check(tag, K.nV, R['V'], 'V')}")
    out(f"     E     got {K.nE:3d}  S4 {R['E']:3d}   {check(tag, K.nE, R['E'], 'E')}")
    out(f"     F     got {K.nF:3d}  S4 {R['F']:3d}   {check(tag, K.nF, R['F'], 'F')}")
    out(f"     chi   got {chi:3d}  S4 {R['chi']:3d}   {check(tag, chi, R['chi'], 'chi')}")
    out(f"     b0    got {b0:3d}  S4 {R['b0']:3d}   {check(tag, b0, R['b0'], 'b0')}")
    out(f"     b1    got {b1:3d}  S4 {R['b1']:3d}   {check(tag, b1, R['b1'], 'b1')}")
    out(f"     b2    got {b2:3d}  S4 {R['b2']:3d}   {check(tag, b2, R['b2'], 'b2')}")
    out(f"     gauge got {gauge:3d}  S4 {R['gauge']:3d}   {check(tag, gauge, R['gauge'], 'gauge')}")
    out(f"     inv   got {inv:3d}  S4 {R['inv']:3d}   {check(tag, inv, R['inv'], 'inv')}")
    out(f"     curv  got {curv:3d}  S4 {R['curv']:3d}   {check(tag, curv, R['curv'], 'curv')}")
    out(f"     flat  got {flat:3d}  S4 {R['flat']:3d}   {check(tag, flat, R['flat'], 'flat')}")
    out(f"     chi = b0-b1+b2 : {chi} = {b0-b1+b2}   {check(tag, chi, b0-b1+b2, 'euler')}")
    out(f"     curv+flat = inv: {curv+flat} = {inv}   {check(tag, curv+flat, inv, 'split')}")

    # loops verified against the boundary maps, not asserted
    cF, cC = K.chain(K.gF), K.chain(K.gC)
    dF, dC = K.d1_times(cF), K.d1_times(cC)
    cols = K.d2_columns()
    F_cycle = all(x == 0 for x in dF)
    C_cycle = all(x == 0 for x in dC)
    F_bounds = in_column_span(cols, cF)
    C_bounds = in_column_span(cols, cC)
    indep = exact_rank([[cF[e], cC[e]] for e in range(K.nE)]) == 2
    out(f"  gamma_F: d1(gF)=0 ? {F_cycle}   bounds ? {F_bounds}  "
        f"{check(tag, F_bounds, R['gF_bounds'], 'gF bounds')}")
    out(f"  gamma_C: d1(gC)=0 ? {C_cycle}   bounds ? {C_bounds}  "
        f"{check(tag, C_bounds, R['gC_bounds'], 'gC bounds')}")
    out(f"  independent (rank[gF|gC] = 2) ? {indep}  "
        f"{check(tag, indep, R['independent'], 'independent')}")
    check(tag, F_cycle, True, "gF is a cycle")
    check(tag, C_cycle, True, "gC is a cycle")

    ms = K.class_multiset()
    cl = K.classes()
    out("  VERTEX CLASSES (a_v = in gamma_F, b_v = in gamma_C):")
    out("     " + "  ".join(f"{K.vnames[v]}:{CLASS_NAME[cl[v]]}" for v in range(K.nV)))
    got_ms = {k: ms[k] for k in [(0, 0), (0, 1), (1, 0), (1, 1)]}
    out(f"  CLASS MULTISET  got {{00:{got_ms[(0,0)]}, 01:{got_ms[(0,1)]}, "
        f"10:{got_ms[(1,0)]}, 11:{got_ms[(1,1)]}}}")
    out(f"                   S4 {{00:{R['multiset'][(0,0)]}, 01:{R['multiset'][(0,1)]}, "
        f"10:{R['multiset'][(1,0)]}, 11:{R['multiset'][(1,1)]}}}   "
        f"{check(tag, got_ms, R['multiset'], 'class multiset')}")
    occ = [c for c in CLASSES if ms[c] > 0]
    out(f"  OCCUPIED CLASSES: {[CLASS_NAME[c] for c in occ]}   (|occ| = {len(occ)})")
    piU = K.pushforward_uniform()
    out(f"  SENSE U pushforward pi = (p00,p10,p01,p11) = "
        f"({piU[0]}, {piU[1]}, {piU[2]}, {piU[3]})   sum = {sum(piU)}")
    out(f"  support difference lattice L_S basis (EXACT) = {L_S(occ)}   rank {len(L_S(occ))}")
    # four-class-only structure: does P factor?
    p00, p10, p01, p11 = piU
    out(f"  factorisation test  p00*p11 - p10*p01 = {p00*p11 - p10*p01}   "
        f"(P factors as (a+bx)(c+dy) iff this is 0)")
    out()

out("=" * 100)
out("VALIDATION SUMMARY")
out("=" * 100)
if FAILS:
    out(f"**{len(FAILS)} MISMATCHES AGAINST S4'S PUBLISHED ROWS — READ THESE LOUDLY**")
    for f in FAILS:
        out("   " + f)
else:
    out("ALL 4 CARRIERS REPRODUCE S4's PUBLISHED ROWS EXACTLY:")
    out("  V, E, F, chi, b0, b1, b2, gauge params, invariants, curvature, flat,")
    out("  gF bounds / gC does not / independent, and the CLASS MULTISET — 0 mismatches.")
    out("  S4 was never audited; these rows are audited now and they hold.")
out()

out("-" * 100)
out("THE FOUR-CLASS FACT THAT IS THE POINT OF THIS SCRIPT")
out("-" * 100)
out("B1 (K1)  occupies {01,10,11}  — a vertex in BOTH loops, none in NEITHER.")
out("B1q      occupies {00,01,10}  — a vertex in NEITHER, none in BOTH.")
out("B0b, B4  occupy  {00,01,10,11} — BOTH a pinch and a spectator.  These are the only two")
out("         carriers of S4's ten that do, exactly as W-09 recorded, and they are built here.")
out()
out("A carrier enters Z_k = sum_ab p_ab (u^a v^b)^k ONLY through pi.  So the ONLY thing a")
out("carrier can decide about any Z-derived result is WHICH pi ARE REALIZABLE, i.e. which")
out("face of the 3-simplex the pushforward is confined to.  On B0b and B4 that face is the")
out("WHOLE simplex; on every other carrier the corpus owns it is a proper face.")
out()
out("EXACT: the factorisation locus p00*p11 = p10*p01 is EMPTY on any three-class carrier —")
out("p00 = 0 forces p10*p01 = 0, i.e. a second empty class.  Symmetrically p11 = 0 forces")
out("p10*p01 = 0.  So 'P factors into a curvature factor times a flat factor' is a")
out("FOUR-CLASS-ONLY phenomenon.  It is exercised in script 3.")

with open("w10a_1_build.OUT.txt", "w") as fh:
    fh.write("\n".join(LOG) + "\n")
sys.exit(1 if FAILS else 0)
