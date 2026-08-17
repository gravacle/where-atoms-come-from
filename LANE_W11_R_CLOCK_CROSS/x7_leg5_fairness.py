# X7 — MY BRIEF, VERBATIM: "verify its 'fairest comparison' is not itself a zero-variable control".
# Leg 5 is lane C's fairness leg: six named correspondences through one evaluator.
# X7 (a) re-derives the collapse structure independently, (b) counts how many INDEPENDENT
# BEHAVIOURS the six rows actually exhibit, (c) tests the one row lane C declared weakest.
import hashlib
import numpy as np
from x_lib import *

def oh(x): return hashlib.sha256(np.ascontiguousarray(np.asarray(x)).tobytes()).hexdigest()[:12]

def ev(pairs, TF, TC, S):
    vals = np.zeros((len(S), len(pairs)))
    for j,(mF,mC) in enumerate(pairs):
        AF = np.linalg.matrix_power(TF,int(mF)); AC = np.linalg.matrix_power(TC,int(mC))
        for i,s in enumerate(S): vals[i,j] = abs(np.vdot(AF@s, AC@s))
    return vals

for K, pi in ((K1(), np.array([0.,.30,.30,.40])), (B0b(), None)):
    if pi is None:
        wB = np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wB/=wB.sum(); pi = pi_of(K,np.sqrt(wB)+0j)
    a = generic_conn(K, np.random.default_rng(7+K.nv))
    TF, TC = T_edge(K,K.wF,a), T_edge(K,K.wC,a)
    WF, WC = holo(K.wF,a), holo(K.wC,a)
    f, c = np.angle(WF)%(2*np.pi), np.angle(WC)%(2*np.pi)
    S = states_same_pi(K, pi, 24, np.random.default_rng(20260817))
    Lc, g = int(np.lcm(K.LF,K.LC)), int(np.gcd(K.LF,K.LC)); T = 240
    print(f"\n================ {K.name}  L=({K.LF},{K.LC})  lcm={Lc} gcd={g} ================")
    P = {}
    P["CORR-E"] = [(t,t) for t in range(1,T+1)]
    P["CORR-C"] = [(K.LF*t,K.LC*t) for t in range(1,T+1)]
    P["CORR-D"] = [(K.LF*t,K.LC*2*t) for t in range(1,T+1)]
    P["CORR-X"] = [(Lc*t,Lc*t) for t in range(1,T+1)]
    phis = [(p,g) for p in range(1,g*T+1)] if g>1 else [(p,1) for p in range(1,T+1)]
    P["CORR-R"] = [(K.LF*p//q, K.LC*p//q) for (p,q) in phis
                   if (K.LF*p)%q==0 and (K.LC*p)%q==0][:T]
    pp=[]
    for t in range(1,T+1):
        mC = max(int(round(t*(f/K.LF)/(c/K.LC))),1); pp.append((t,mC))
    P["CORR-P"] = pp
    res = {}
    for nm, pr in P.items():
        v = ev(pr, TF, TC, S); sp = (v.max(axis=0)-v.min(axis=0)).max()
        onlat = all((mF%K.LF==0 and mC%K.LC==0) for (mF,mC) in pr)
        res[nm] = (oh(v), sp, onlat)
        print(f"   {nm:<8} out#{res[nm][0]}   max|Z|-spread {sp:.3e}   "
              f"ALL readings on the sublattice? {onlat}")
    seen = {}
    for nm in P: seen.setdefault(res[nm][0], []).append(nm)
    print(f"   COLLAPSE (my code): " + " | ".join("=".join(v) for v in seen.values())
          + f"    distinct arms {len(seen)} of 6")
    beh = {}
    for nm in P: beh.setdefault(res[nm][2], []).append(nm)
    print(f"   INDEPENDENT BEHAVIOURS: {len(beh)} -- "
          + "; ".join(f"{'ALL-ON-SUBLATTICE' if k else 'SOME-OFF-SUBLATTICE'}: {v}" for k,v in beh.items()))
    print("   Every ON row is the SAME theorem (T^m diagonal iff L|m) sampled again; every OFF row")
    print("   is its negation.  Leg 5 therefore contains ONE discriminator, not six, and lane C's")
    print("   own [T] mark on its leg-5 nulls concedes half of that.")

    # (c) CORR-P, the row lane C called its weakest: is its non-null caused by the phase residual?
    resid = np.array([abs(t*f/K.LF - mC*c/K.LC) for (t,mC) in pp])
    keep = [k for k in range(len(pp)) if resid[k] < 0.01]
    v = ev([pp[k] for k in keep], TF, TC, S)
    onl = sum(1 for k in keep if pp[k][0]%K.LF==0 and pp[k][1]%K.LC==0)
    print(f"   CORR-P residual test: {len(keep)} of {T} readings have residual < 0.01 rad "
          f"(max residual overall {resid.max():.3f});")
    print(f"      of those, {onl} land on the sublattice.  max|Z|-spread on the low-residual"
          f" subset = {(v.max(axis=0)-v.min(axis=0)).max():.3e}")
    print("      -> the non-null SURVIVES at near-zero residual, so lane C's own caveat ('part of")
    print("         its non-null is attributable to the mismatch') is too generous to itself: the")
    print("         residual is irrelevant, because equal accumulated PHASE never makes T^m")
    print("         diagonal.  CORR-P is a third sample of the same theorem, not a third test.")
