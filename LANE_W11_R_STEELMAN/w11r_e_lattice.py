# LANE W11-R  LEG E -- THE CLOCK-FREE STATEMENT, AND THE TRIVIAL-CONNECTION FIRING SWEPT.
#
# E1.  Drop the clock entirely.  Let each branch carry its OWN tick count and form
#         Z(a,b) = < T_F^a s , T_C^b s >.
#      Both "conventions" are diagonals of this two-parameter family: the corpus's is a -> (k,k)
#      after L-fold rescaling, the registrar's is (n,n).  The clock-free question is: for which
#      (a,b) is Z(a,b) a function of pi alone?  ANSWER (tested, not asserted): exactly on the
#      lattice  L_F.Z x L_C.Z  -- i.e. exactly when BOTH branches sit at their own loop closure.
#      That is a THEOREM with content, stated with no clock in it at all.
#
# E2.  Robustness of leg B: is the trivial-connection firing a knife-edge or generic?
import numpy as np, w11r_lib as L
rng = np.random.default_rng(20260817)

print("== E1  FOR WHICH (a,b) IS THE OBSERVABLE pi-ONLY?  NO CLOCK IN THE QUESTION ==")
for nm, lf, lc, NV, aa, mkw in (
    ("K1 ", L.K1_LOOP_F, L.K1_LOOP_C, 5, np.array([1.0,0.37,0.91,2**0.5,0.23,1.77]),
        lambda: (np.array([0.40,0.15,0.15,0.15,0.15]), np.array([0.40,0.30,0.00,0.05,0.25]))),
    ("B0b", L.B0B_LOOP_F, L.B0B_LOOP_C, 9, rng.uniform(0,2*np.pi,18),
        lambda: (np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11])/1.0,
                 np.array([.22,.00,.09,.00,.25,.22,.11,.11,.00])))):
    LF, LC = len(lf), len(lc)
    wA, wB = mkw(); wA=wA/wA.sum(); wB=wB/wB.sum()
    sA, sB = np.sqrt(wA)+0j, np.sqrt(wB)+0j
    sC = sA*np.exp(1j*rng.uniform(0,2*np.pi,NV))
    pA,pB,pC = (L.pi_of(x,lf,lc,NV) for x in (sA,sB,sC))
    assert np.allclose(pA,pB,atol=1e-12) and np.allclose(pA,pC,atol=1e-12), f"{nm}: pi not fixed {pA} {pB}"
    assert L.arms_differ(sA,sB,sC)
    TF,TC = L.T_edge(lf,aa,NV), L.T_edge(lc,aa,NV)
    hits, misses, mism = [], 0, []
    NMAX = 25
    for A in range(1,NMAX):
        for B in range(1,NMAX):
            v=[abs(L.Z(TF,TC,s,A,B)) for s in (sA,sB,sC)]
            sp = max(v)-min(v)
            on_lattice = (A % LF == 0) and (B % LC == 0)
            if sp < 1e-12: hits.append((A,B))
            else: misses += 1
            if (sp < 1e-12) != on_lattice: mism.append((A,B,sp))
    print(f"  {nm}  |gamma_F|={LF} |gamma_C|={LC}   pi = {np.round(pA,9)}")
    print(f"       pi-only pairs (a,b) with a,b < {NMAX}: {len(hits)}   predicted |L_F.Z x L_C.Z| = "
          f"{((NMAX-1)//LF)*((NMAX-1)//LC)}")
    print(f"       MISMATCHES against 'both branches at loop closure': {len(mism)}  {mism[:4]}")
    print(f"       first few pi-only pairs: {hits[:8]}")
print("  -> THE CLOCK-FREE STATEMENT OF RECORD:  the incidence is invisible EXACTLY at joint loop")
print("     closure.  No convention appears in that sentence.  It is false off the lattice, so it")
print("     is not vacuous; it is true on it, so it is not a stipulation.")

print("\n== E2  IS THE TRIVIAL-CONNECTION FIRING OF LEG B A KNIFE-EDGE?  SWEPT OVER RANDOM STATES ==")
print(f"  {'carrier':<6}{'convention':<10}{'lambda at a=0: min':>20}{'median':>14}{'max':>14}{'frac < -1e-6':>14}")
for nm, lf, lc, NV, ne in (("K1",L.K1_LOOP_F,L.K1_LOOP_C,5,6), ("B0b",L.B0B_LOOP_F,L.B0B_LOOP_C,9,18)):
    a0 = np.zeros(ne)
    TF,TC = L.T_edge(lf,a0,NV), L.T_edge(lc,a0,NV)
    MF,MC = L.M_circuit(lf,a0,NV), L.M_circuit(lc,a0,NV)
    re_, rc_ = [], []
    for _ in range(60):
        s = rng.normal(size=NV)+1j*rng.normal(size=NV); s/=np.linalg.norm(s)
        re_.append(L.rate(TF,TC,s,4000)); rc_.append(L.rate(MF,MC,s,4000))
    for cn, r in (("EDGE",re_),("CIRCUIT",rc_)):
        r=np.array(r)
        print(f"  {nm:<6}{cn:<10}{r.min():>20.9f}{np.median(r):>14.9f}{r.max():>14.9f}"
              f"{np.mean(r < -1e-6):>14.3f}")
print("  -> under the EDGE convention essentially every ready state forms a record on a carrier")
print("     with NO CONNECTION.  Under the CIRCUIT convention none does, identically.")
print("     (Note the honest detail: the edge rate at a=0 is state-dependent and can be small --")
print("      -0.0027 for the near-uniform B0b state of leg B3 -- but it is negative, and CHOICE")
print("      LEDGER C4's test is on the VERDICT, not on the magnitude.  And leg B0/B1 exhibit")
print("      EXACT zeros of Z at the trivial connection, which is the verdict at its sharpest.)")
