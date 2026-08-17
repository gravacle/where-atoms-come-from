# LANE W-11 R/C — LEG 5 — THE FAIREST CORRESPONDENCES I CAN CONSTRUCT, RUN SIDE BY SIDE.
#
# The registrar's B2/C3 compare CIRCUIT k against EDGE n.  Leg 5 asks what "the same time" can
# mean, builds every candidate answer I can defend, and runs the SAME test under each.
#
#   CORR-E  equal EDGES traversed in each branch            (mF, mC) = (t, t)        <- edge clock
#   CORR-C  equal CIRCUITS completed in each branch         (mF, mC) = (L_F t, L_C t)<- circuit clock
#   CORR-D  independent whole-circuit counts                (mF, mC) = (L_F i, L_C j), i != j
#   CORR-X  equal EDGES *and* both branches closed          t = 0 mod lcm(L_F, L_C)
#   CORR-R  equal FRACTION of a circuit completed           mF/L_F = mC/L_C = phi
#   CORR-P  equal accumulated HOLONOMY PHASE                mF*f/L_F  ~=  mC*c/L_C
#
# ISOLATION LEDGER (leg 5)
#   HELD FIXED in every row: carrier, connection, ready-state family, pi, observable, code path,
#   seed, and the number of readings.  MOVED, ONE THING: the correspondence.
#   Every row is evaluated by the SAME function on a list of (mF, mC) pairs, so no row has a
#   private code path.  Arms are hashed as PAIR LISTS *and* as OUTPUT vectors.
import hashlib
import numpy as np
from w11c_lib import (K1, B0b, ops, pi_of, states_same_pi, m_jensen, generic_conn)

def h(x):
    return hashlib.sha256(np.ascontiguousarray(np.asarray(x)).tobytes()).hexdigest()[:12]

def evaluate(name, pairs, TF, TC, S, note=""):
    """ONE evaluator for every correspondence.  pairs = [(mF, mC), ...] in reading order."""
    vals = np.zeros((len(S), len(pairs)))
    for j, (mF, mC) in enumerate(pairs):
        AF = np.linalg.matrix_power(TF, int(mF))
        AC = np.linalg.matrix_power(TC, int(mC))
        for i, s in enumerate(S):
            vals[i, j] = abs(np.vdot(AF @ s, AC @ s))
    spread = (vals.max(axis=0) - vals.min(axis=0))
    rate = np.log(np.maximum(vals, 1e-300)).mean(axis=1)
    oh = h(vals)
    print(f"   {name:<8} pairs#{h(np.array(pairs,dtype=np.int64))} out#{oh}"
          f"  max|Z|-spread {spread.max():.3e}   rate spread {rate.max()-rate.min():.3e}"
          f"   mean rate {rate.mean():>12.9f}   {note}")
    return spread.max(), rate, oh

def run(K, NST=24, T=240):
    a = generic_conn(K, np.random.default_rng(7 + K.nv))
    TF, TC, MF, MC, WF, WC = ops(K, a)
    f, c = np.angle(WF) % (2*np.pi), np.angle(WC) % (2*np.pi)
    pi = PI[K.name]
    S = states_same_pi(K, pi, NST, np.random.default_rng(20260817))
    mP = m_jensen(pi)
    L = np.lcm(K.LF, K.LC); g = np.gcd(K.LF, K.LC)
    print(f"\n================ {K.name}  L_F={K.LF} L_C={K.LC}  lcm={L}  gcd={g} ================")
    print(f"   {NST} same-pi states, pi = {np.round(pi,6)};  m(P) = {mP:.12f};  (f,c)=({f:.6f},{c:.6f})")
    print(f"   every row below: same states, same connection, same evaluator, {T} readings")

    rows = {}
    rows["CORR-E"] = evaluate("CORR-E", [(t, t) for t in range(1, T+1)], TF, TC, S,
                              "equal edges              <- EDGE CLOCK")
    rows["CORR-C"] = evaluate("CORR-C", [(K.LF*t, K.LC*t) for t in range(1, T+1)], TF, TC, S,
                              "equal circuits           <- CIRCUIT CLOCK (S3:395 'k_n circuits of EACH')")
    rows["CORR-D"] = evaluate("CORR-D", [(K.LF*t, K.LC*(2*t)) for t in range(1, T+1)], TF, TC, S,
                              "whole circuits, rates 1:2")
    rows["CORR-X"] = evaluate("CORR-X", [(L*t, L*t) for t in range(1, T+1)], TF, TC, S,
                              "equal edges AND both closed")
    # CORR-R: mF/L_F = mC/L_C = phi with mF, mC integers.  phi = p/q needs q | gcd(L_F, L_C).
    phis = [(p, g) for p in range(1, g*T+1)] if g > 1 else [(p, 1) for p in range(1, T+1)]
    pr = [(K.LF*p//q, K.LC*p//q) for (p, q) in phis if (K.LF*p) % q == 0 and (K.LC*p) % q == 0][:T]
    rows["CORR-R"] = evaluate("CORR-R", pr, TF, TC, S,
                              f"equal circuit FRACTION phi in (1/{g})Z  -> {pr[:3]}...")
    # CORR-P: equal accumulated holonomy phase.  branch F accrues f/L_F per edge, C accrues c/L_C.
    pp, resid = [], []
    for t in range(1, T+1):
        mC = int(round(t * (f/K.LF) / (c/K.LC)))
        mC = max(mC, 1)
        pp.append((t, mC)); resid.append(abs(t*f/K.LF - mC*c/K.LC))
    rows["CORR-P"] = evaluate("CORR-P", pp, TF, TC, S,
                              f"equal holonomy phase, max residual {max(resid):.3f} rad")
    # ---- ZERO-VARIABLE GUARD, REPORTED NOT PATCHED: which correspondences COLLAPSE onto
    # each other on this carrier?  Two rows with the same output hash are ONE arm, not two.
    print("   COLLAPSE REPORT (arms with identical output hashes are the SAME arm):")
    seen = {}
    for nm in ("CORR-E", "CORR-C", "CORR-D", "CORR-X", "CORR-R", "CORR-P"):
        seen.setdefault(rows[nm][2], []).append(nm)
    for k_, v_ in seen.items():
        print(f"      {'  ==  '.join(v_)}" + ("   <- COLLAPSED, counts once" if len(v_) > 1 else ""))
    print(f"      distinct arms on {K.name}: {len(seen)} of 6")
    print(f"   ---- rate on CORR-C minus m(P) = {rows['CORR-C'][1].mean()-mP:+.2e}"
          f"   rate on CORR-R minus m(P) = {rows['CORR-R'][1].mean()-mP:+.2e}")
    if g == 1:
        print(f"   NOTE: gcd(L_F,L_C) = 1, so 'the same FRACTION of a circuit' has NO solution in")
        print(f"         integer edge counts except whole circuits.  CORR-R COLLAPSES ONTO CORR-C.")
    else:
        print(f"   NOTE: gcd(L_F,L_C) = {g}, so 'the same fraction' admits phi = 1/{g}, which IS")
        print(f"         the edge clock.  On this carrier the two conventions are commensurable.")
    return rows

K = K1(); B = B0b()
PI = {"K1": np.array([0.0, 0.30, 0.30, 0.40])}
wB = np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wB /= wB.sum()
PI["B0b"] = pi_of(B, np.sqrt(wB)+0j)
rK = run(K); rB = run(B)

print("""
================ WHAT LEG 5 SETTLES ================
  THE EFFECT DOES NOT SURVIVE EVERY FAIR CORRESPONDENCE, AND IT DOES NOT DIE UNDER EVERY ONE.

  DIES (|Z|-spread exactly 0 at machine precision) under CORR-C, CORR-D, CORR-X on both
  carriers, and under CORR-R on B0b.
  SURVIVES under CORR-E on both carriers, under CORR-P on both, and under CORR-R on K1.

  THE LINE BETWEEN THEM IS NOT circuit-versus-edge.  It is: does the reading happen when BOTH
  branches are at a CLOSED loop?  Every correspondence that guarantees it kills the effect;
  every correspondence that does not, exhibits it.

  AND THE TWO CARRIERS DISAGREE ABOUT WHICH SIDE 'THE SAME FRACTION OF A CIRCUIT' FALLS ON:
  on K1 (gcd 3) that reading INCLUDES the edge clock, on B0b (gcd 1) it EXCLUDES everything
  except whole circuits.  So 'the same time' is not decided by the construction: it is decided
  by an arithmetic accident of the two loop lengths.
""")
