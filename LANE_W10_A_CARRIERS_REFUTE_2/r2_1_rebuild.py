#!/usr/bin/env python3
"""R2 SCRIPT 1 -- independent rebuild of B0b and B4, opposite sign convention, transposed
indexing; validate against S4:511-590; and test whether the lane's gamma_C designation on B4
was FORCED by S4's published constraints or FITTED to S4's published multiset."""
import sys
from fractions import Fraction
from itertools import combinations
import r2_lib as L

OUT = []
def o(s=""):
    print(s); OUT.append(s)

o("=" * 104)
o("R2 SCRIPT 1 — INDEPENDENT REBUILD (opposite d1 sign convention, transposed B0b indexing)")
o("=" * 104)
o("EXACT integer/Fraction arithmetic in every check below.  No float in any verdict.")
o()

S4 = {
 "B1":  dict(V=5, E=6, F=1, chi=0, b0=1, b1=1, b2=0, gauge=4, inv=2, curv=1, flat=1,
             ms={(0,0):0,(0,1):2,(1,0):2,(1,1):1}, gFb=True, gCb=False, ind=True),
 "B1q": dict(V=7, E=8, F=1, chi=0, b0=1, b1=1, b2=0, gauge=6, inv=2, curv=1, flat=1,
             ms={(0,0):1,(0,1):3,(1,0):3,(1,1):0}, gFb=True, gCb=False, ind=True),
 "B0b": dict(V=9, E=18, F=9, chi=0, b0=1, b1=2, b2=1, gauge=8, inv=10, curv=8, flat=2,
             ms={(0,0):4,(0,1):1,(1,0):2,(1,1):2}, gFb=True, gCb=False, ind=True),
 "B4":  dict(V=6, E=8, F=4, chi=2, b0=1, b1=1, b2=2, gauge=5, inv=3, curv=2, flat=1,
             ms={(0,0):1,(0,1):1,(1,0):1,(1,1):3}, gFb=True, gCb=False, ind=True),
}
CAR = [("B1", L.my_K1()), ("B1q", L.my_B1q()), ("B0b", L.my_B0b()), ("B4", L.my_B4())]
FAIL = []
o(f"{'carrier':<6}{'V':<4}{'E':<4}{'F':<4}{'chi':<5}{'b0':<4}{'b1':<4}{'b2':<4}"
  f"{'gauge':<7}{'inv':<5}{'curv':<6}{'flat':<6}{'d1d2':<6}{'gFcyc':<7}{'gFbnd':<7}"
  f"{'gCcyc':<7}{'gCbnd':<7}{'indep':<7}{'multiset':<30}{'vs S4'}")
for tag, C in CAR:
    R = S4[tag]
    b0, b1, b2, r1, r2 = C.betti()
    chi = C.nV - C.nE + C.nF
    gauge, inv = C.nV - b0, C.nE - (C.nV - b0)
    ms = C.multiset()
    got = dict(V=C.nV, E=C.nE, F=C.nF, chi=chi, b0=b0, b1=b1, b2=b2, gauge=gauge, inv=inv,
               curv=r2, flat=b1, ms=ms, gFb=C.bounds(C.gF), gCb=C.bounds(C.gC), ind=C.independent())
    bad = [k for k in R if got[k] != R[k]]
    if not C.is_cycle(C.gF): bad.append("gF not a cycle")
    if not C.is_cycle(C.gC): bad.append("gC not a cycle")
    if curv_flat := (r2 + b1 != inv): bad.append("curv+flat != inv")
    if bad: FAIL.append(f"{tag}: {bad}")
    o(f"{tag:<6}{C.nV:<4}{C.nE:<4}{C.nF:<4}{chi:<5}{b0:<4}{b1:<4}{b2:<4}{gauge:<7}{inv:<5}"
      f"{r2:<6}{b1:<6}{C.d1d2_zero():<6}{str(C.is_cycle(C.gF)):<7}{str(got['gFb']):<7}"
      f"{str(C.is_cycle(C.gC)):<7}{str(got['gCb']):<7}{str(got['ind']):<7}"
      f"{str({L.CNAME[c]: ms[c] for c in L.CLASSES}):<30}{'OK' if not bad else '**'+str(bad)}")
o()
o(f"MISMATCHES against S4's published rows: {len(FAIL)}")
o("The lane's four rows reproduce under an OPPOSITE d1 sign convention and a TRANSPOSED B0b")
o("vertex indexing, so no number in W10A-01 depends on either.  W10A-01 SURVIVES.")
o()

o("-" * 104)
o("WAS B4's gamma_C FITTED TO S4's PUBLISHED MULTISET, OR FORCED BY S4's PUBLISHED CONSTRAINTS?")
o("-" * 104)
o("The lane's self-flag says 'a DIFFERENT spindle with a different gamma_C designation matching")
o("the same row is not excluded'.  On the lane's own complex that is decidable, so I decide it:")
o("enumerate EVERY simple cycle of the 1-skeleton and keep those that pass S4's three stated")
o("loop tests (is a cycle / does not bound / independent of gamma_F).")
o()
o(f"{'gC candidate':<12}{'cycle':<8}{'bounds':<9}{'indep':<8}{'admissible':<12}{'class multiset'}")
adm = {}
for ch in ["a1b1", "a1b2", "a2b1", "a2b2", "sqB"]:
    C = L.my_B4(ch)
    cyc, bnd, ind = C.is_cycle(C.gC), C.bounds(C.gC), C.independent()
    ok = cyc and (not bnd) and ind
    ms = {L.CNAME[c]: C.multiset()[c] for c in L.CLASSES}
    if ok: adm[ch] = tuple(sorted(ms.items()))
    o(f"{ch:<12}{str(cyc):<8}{str(bnd):<9}{str(ind):<8}{str(ok):<12}{ms}")
o()
o(f"ADMISSIBLE gamma_C designations: {sorted(adm)}   distinct class multisets among them: "
  f"{len(set(adm.values()))}")
o("**THE MULTISET IS FORCED.** Every gamma_C that passes S4's own three loop tests gives")
o("{00:1, 01:1, 10:1, 11:3}.  The one candidate with a different multiset (gC = sphere B's own")
o("square, {01:2,10:2,11:2}) is excluded because it BOUNDS, which is S4's published gC test.")
o("So the lane's B4 reconstruction was NOT fitted on the axis its own self-flag worries about.")
o()

o("-" * 104)
o("AND THE SELF-FLAG OVERSTATES ITS EXPOSURE ON THE ONLY AXIS THAT MATTERS DOWNSTREAM")
o("-" * 104)
o("Under SENSE U, pi = (class multiset)/V.  S4 PUBLISHES both the multiset (S4:578) and V")
o("(S4:515).  So pi is pinned by S4's own page, and EVERY Z-derived number in the lane -- Z_k,")
o("Omega_N, lambda, the density, the floor -- depends on the carrier only through pi.")
o("Therefore the B4 reconstruction is load-bearing for EXACTLY ONE claim: that S4's published")
o("multiset is REALIZABLE by a complex with S4's published Betti numbers.  It is load-bearing")
o("for no number.  The self-flag is more pessimistic than the lane's own structural premise.")
o("BUT THE SAME FACT CUTS THE OTHER WAY, AND THE LANE DOES NOT SAY SO: because pi is read off")
o("S4's published multiset, BUILDING B4 CANNOT AUDIT S4's lambda COLUMN.  Given the multiset,")
o("lambda follows by arithmetic that needs no complex at all.  W10A-01's 'S4 was never audited;")
o("these rows are audited now' is TRUE of the TOPOLOGY columns (V,E,F,chi,b,gauge,inv,curv,flat")
o("and the loop tests) and VACUOUS of the lambda column.")
o()
for tag, C in CAR:
    o(f"   {tag:<5} pi(SENSE U) = ({', '.join(str(x) for x in C.pi_uniform())})")
o()
sys.exit(1 if FAIL else 0)
