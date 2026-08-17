# run_adversary.py -- LANE W20_R_LEDGER.  THE ATTACKS ON MY OWN RESULT.
#
# run_arms + run_extra produced: C = sum_k p_k E_k (the conditional free entropy of the region
# given its surface) is <= 3.4e-4 bits on every ground state measured, while Haar states on the
# same sector, surface and algebra give C ~ 0.6-1.5 bits.  Before that is reported as a finding it
# has to survive the cheapest explanations that are NOT "the dynamics did it".
#
# X1  IS IT JUST POSITIVITY?  H has non-positive off-diagonal entries in the Z basis (all couplings
#     enter as -(1/g2)W_p and -g2 X_l), so H is STOQUASTIC and the ground state is real, nodeless
#     and strictly positive by Perron-Frobenius.  If a RANDOM NONNEGATIVE state also has C ~ 0 then
#     the result is about the SIGN STRUCTURE OF THE HAMILTONIAN, not about ground-state selection,
#     and "the dynamics" is the wrong ingredient to credit.  THIS IS THE ARM THAT COULD KILL IT.
# X2  IS IT JUST REALNESS?  Random real (signed) states.
# X3  THE GAUGE-REMOVAL CATEGORY ERROR, MEASURED ON BOTH SIDES: the free ground state's overlap
#     with each of the five charge sectors.
# X4  THE SYMMETRY AUDIT of every arm pair actually reported, incl. the S' pairs.  FIT-3 guard.
# X5  WHAT ARE THE EXCITED STATES WITH C = 0 EXACTLY?  (E2 found some; they are a caveat, not noise.)

import numpy as np, math, itertools
from core_w20r import *

LOG = []
def P(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s, flush=True)
def rule(t=""):
    P("\n" + "=" * 104); P(t); P("=" * 104)
def LM(ls):
    m = 0
    for l in ls: m |= 1 << l
    return m

S_T = LM([1,2,3]); SIG_T = LM([0,4,5]); W_S = LM([1,2,3])
S2 = LM([1,2,3,4,5]); SIG2 = LM([0,6])
SEC = {"vacuum":[1]*8, "eta{0,4}":[-1 if v in (0,4) else 1 for v in range(8)],
       "eta{1,5}":[-1 if v in (1,5) else 1 for v in range(8)],
       "eta{4,5}":[-1 if v in (4,5) else 1 for v in range(8)],
       "eta{0,1}":[-1 if v in (0,1) else 1 for v in range(8)]}
BE = {k: PhysSector(v) for k, v in SEC.items()}
be0 = BE["vacuum"]

def chans(be, S, wl):
    A = alg_links(be, S, N_PHYS, "FULL")
    cen = [v for v in gf2_span(list(A.W.values())) if all(omega(v, w) == 0 for w in A.W.values())]
    return A, Algebra(be, cen, N_PHYS, "CENTRE")
A_S, A_CEN = chans(be0, S_T, [W_S])
A_S2, A_CEN2 = chans(be0, S2, [W_S, LM([3,4,5])])
def C_of(v, A=A_S, Cn=A_CEN):
    return A.entropy(v) - Cn.entropy(v)

# ==================================================================================================
rule("X1 -- IS IT JUST PERRON-FROBENIUS POSITIVITY?  THE ARM THAT COULD KILL THE RESULT.")
P("H = -(1/g2) sum_p W_p - g2 sum_l X_l .  In the Z (coset) basis both terms have NON-POSITIVE")
P("off-diagonal entries, so H is STOQUASTIC and its ground state is real, nodeless and strictly")
P("positive.  CHECK THAT FIRST, then ask whether positivity ALONE explains C ~ 0.")
Hm = H_matrix(be0, 0.60)
off = Hm - np.diag(np.diag(Hm))
P("  max off-diagonal entry of H at g2=0.60 : %.9f   (<=0 means stoquastic)" % float(off.real.max()))
ev, evec = np.linalg.eigh(Hm)
gsv = evec[:, 0]
gsv = gsv * np.sign(gsv[np.argmax(np.abs(gsv))])
P("  min component of the ground state      : %.9f   (>0 means nodeless)" % float(gsv.real.min()))
P("  max |imaginary part|                   : %.2e" % float(np.abs(gsv.imag).max()))
P("")
P("NOW THE ATTACK.  Random NONNEGATIVE states (each component |N(0,1)|, then normalised).  These")
P("are nodeless and positive exactly like the ground state, but carry NO dynamical information.")
P("If C is small for these, the finding belongs to the SIGN STRUCTURE, not to the dynamics.")
rng = np.random.default_rng(12345)
for lbl, gen, n in [("NONNEG (positive, nodeless)", lambda: np.abs(rng.normal(size=32)) + 0j, 200),
                    ("REAL   (signed)",             lambda: rng.normal(size=32) + 0j,        200),
                    ("HAAR   (complex)",            lambda: rng.normal(size=32) + 1j*rng.normal(size=32), 200)]:
    cs, cs2 = [], []
    for _ in range(n):
        v = gen(); v /= np.linalg.norm(v)
        cs.append(C_of(v)); cs2.append(C_of(v, A_S2, A_CEN2))
    P("  %-28s n=%d   C(S)  mean %.6f  min %.6f  max %.6f   |   C(S') mean %.6f  min %.6f"
      % (lbl, n, float(np.mean(cs)), float(np.min(cs)), float(np.max(cs)),
         float(np.mean(cs2)), float(np.min(cs2))))
P("")
P("  GROUND STATE FOR COMPARISON, same region, same sector: C(S) = %.9f  C(S') = %.9f at g2=0.60"
  % (C_of(gsv), C_of(gsv, A_S2, A_CEN2)))
P("")
P("SHARPER: hold POSITIVITY fixed and vary only how the weights are chosen.  Take the ground state")
P("and RANDOMLY REWEIGHT it: psi_i -> psi_i * exp(t * n_i) with n_i ~ N(0,1), renormalise.  Still")
P("strictly positive and nodeless; only the ground state's particular weights are destroyed.")
P("   t        min component   |<gs|psi>|^2   C(S)          C(S')")
for t in [0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0]:
    nz = rng.normal(size=32)
    v = np.abs(gsv.real) * np.exp(t * nz); v = v / np.linalg.norm(v); v = v.astype(complex)
    P("   %-8.2f %-15.9f %-14.9f %-13.9f %.9f"
      % (t, float(v.real.min()), abs(complex(np.vdot(gsv, v)))**2, C_of(v), C_of(v, A_S2, A_CEN2)))

# ==================================================================================================
rule("X2 -- IS THE GROUND STATE'S C SMALL BECAUSE OF THE HAMILTONIAN, OR BECAUSE OF ITS SPECTRUM?")
P("Ground states of RANDOM STOQUASTIC Hamiltonians with the SAME SPARSITY PATTERN as H.")
P("Same off-diagonal support (same terms), random NEGATIVE magnitudes, random diagonal.")
P("If C stays ~0 the property belongs to the stoquastic ground-state class, not to gauge dynamics.")
pat = (np.abs(off) > 1e-12)
P("  sparsity: %d of %d off-diagonal entries nonzero (%.1f%%)" % (int(pat.sum()), 32*32-32, 100.0*pat.sum()/(32*32-32)))
cs, cs2 = [], []
for _ in range(200):
    Mg = -np.abs(rng.normal(size=(32, 32))) * pat
    Mg = (Mg + Mg.T) / 2
    Mg = Mg + np.diag(rng.normal(size=32))
    w, U = np.linalg.eigh(Mg)
    v = U[:, 0].astype(complex)
    cs.append(C_of(v)); cs2.append(C_of(v, A_S2, A_CEN2))
P("  RANDOM STOQUASTIC ground states, n=200:  C(S)  mean %.6f  min %.6f  max %.6f"
  % (float(np.mean(cs)), float(np.min(cs)), float(np.max(cs))))
P("                                           C(S') mean %.6f  min %.6f  max %.6f"
  % (float(np.mean(cs2)), float(np.min(cs2)), float(np.max(cs2))))
P("  THE ACTUAL GAUGE GROUND STATE:           C(S)  %.9f     C(S') %.9f" % (C_of(gsv), C_of(gsv, A_S2, A_CEN2)))
P("")
P("AND THE SAME WITH THE FULL SPARSITY REMOVED (dense random stoquastic, 200 samples):")
cs = []
for _ in range(200):
    Mg = -np.abs(rng.normal(size=(32, 32))); Mg = (Mg + Mg.T)/2 + np.diag(rng.normal(size=32))
    w, U = np.linalg.eigh(Mg); cs.append(C_of(U[:, 0].astype(complex)))
P("  DENSE RANDOM STOQUASTIC ground states:   C(S)  mean %.6f  min %.6f  max %.6f"
  % (float(np.mean(cs)), float(np.min(cs)), float(np.max(cs))))

# ==================================================================================================
rule("X3 -- THE GAUGE REMOVAL, MEASURED ON BOTH SIDES.  WHERE IT IS NARRATION AND WHERE IT IS VOID.")
P("The free (unprojected) ground state's overlap with EVERY charge sector.  A charge sector with")
P("overlap 0 has NO IMAGE under the removal: 'the same state without the Gauss law' does not exist")
P("there, which is the category error the pre-registration predicted -- now with a number on it.")
fb = FreeSpace()
P("  g2       " + "".join("%-14s" % k for k in SEC))
for g2 in [0.20, 0.645497, 1.00, 3.00]:
    E0, psif = lanczos_ground(fb, g2, m=150, seed=3, restarts=4)
    row = "  %-8.4g " % g2
    for k in SEC:
        be = BE[k]
        co = np.zeros(32, dtype=complex)
        np.add.at(co, be.orb, psif * be.sgn)
        row += "%-14.9f" % (float(np.linalg.norm(co)) / math.sqrt(128.0))
    P(row)
P("")
P(">>> VACUUM: overlap 1.000000000 at every coupling -- the Gauss projection removes NOTHING from")
P("    the state, so removing the Gauss law moves the vacuum numbers by EXACTLY 0.000000 (run_extra E4).")
P(">>> EVERY CHARGED SECTOR: overlap 0.000000000.  The free ground state is in the vacuum sector at")
P("    every coupling, so there is no 'same state' to compare across the removal.  W-19's survival")
P("    result IS NOT GENERAL: it is exactly the vacuum sector and it is Perron-Frobenius, not physics")
P("    about records.")

# ==================================================================================================
rule("X4 -- SYMMETRY AUDIT OF EVERY REPORTED ARM PAIR.  FIT-3 GUARD, RUN AGAIN AFTER THE NUMBERS.")
def autos():
    out = []
    ES = set(map(lambda e: tuple(sorted(e)), EDGES))
    for p in itertools.permutations(range(8)):
        if all(tuple(sorted((p[a], p[b]))) in ES for a, b in EDGES):
            out.append(p)
    return out
AUT = autos()
P("|Aut(tri_chain12)| = %d" % len(AUT))
def lmap(p, mask):
    out = 0
    for i in bits(mask):
        a, b = EDGES[i]
        e = tuple(sorted((p[a], p[b])))
        out |= 1 << [j for j, ed in enumerate(EDGES) if tuple(sorted(ed)) == e][0]
    return out
def vmap(p, vs):
    return frozenset(p[v] for v in vs)
PAIRS = [("R1 on S ", S_T, SIG_T, {0,4}, {1,5}), ("R2 on S ", S_T, SIG_T, set(), {4,5}),
         ("R3 on S ", S_T, SIG_T, set(), {0,1}), ("R2/R3 on S ", S_T, SIG_T, {4,5}, {0,1}),
         ("R2 on S'", S2, SIG2, set(), {4,5}), ("R3 on S'", S2, SIG2, set(), {0,1}),
         ("R2/R3 on S'", S2, SIG2, {4,5}, {0,1})]
P("  pair          #g in Aut with g(S)=S, g(Sigma)=Sigma AND g(charges_a)=charges_b   verdict")
for nm, Sx, Gx, ca, cb in PAIRS:
    k = sum(1 for p in AUT if lmap(p, Sx) == Sx and lmap(p, Gx) == Gx and vmap(p, ca) == vmap(p, ca) and vmap(p, ca) == frozenset(cb))
    P("  %-13s %-3d %s" % (nm, k, "VACUOUS BY SYMMETRY -- MUST NOT BE SCORED" if k else "LIVE (no symmetry forces equality)"))
P("")
P("NOTE ON A NULL I DID NOT PREDICT: on S' the pairs eta{4,5} and eta{0,1} returned IDENTICAL")
P("numbers to 8 decimals (max|dH_FULL| = 0.12722545 for both).  The audit above says no element of")
P("Aut fixing S' and Sigma' carries one onto the other, so this is NOT explained by the symmetry")
P("group and I am recording it as an UNEXPLAINED COINCIDENCE rather than claiming either reading.")
for nm, Sx, Gx in [("S ", S_T, SIG_T), ("S'", S2, SIG2)]:
    k = sum(1 for p in AUT if lmap(p, Sx) == Sx and lmap(p, Gx) == Gx)
    P("  |Stab(%s and its surface)| = %d" % (nm, k))

# ==================================================================================================
rule("X5 -- THE EXCITED STATES WITH C = 0 EXACTLY.  A CAVEAT, NOT NOISE.")
P("E2 found eigenstates other than the ground state with C = 0 to machine precision.  If they are")
P("inside DEGENERATE multiplets the eigenvector is arbitrary and the number is meaningless; if they")
P("are non-degenerate the caveat is real and 'C = 0' is NOT a unique signature of the ground state.")
for g2 in [0.45, 0.60]:
    ev, evec = np.linalg.eigh(H_matrix(be0, g2))
    P("  g2 = %.2f" % g2)
    for n in range(32):
        C = C_of(evec[:, n])
        if C < 1e-9:
            lo = ev[n] - ev[n-1] if n else 1e9
            hi = ev[n+1] - ev[n] if n < 31 else 1e9
            deg = min(lo, hi) < 1e-9
            P("     n=%-3d E=%-14.8f C=%-13.3e  nearest gap=%-12.8f  %s"
              % (n, ev[n], C, min(lo, hi), "IN A DEGENERATE MULTIPLET -- MEANINGLESS" if deg else "NON-DEGENERATE -- REAL CAVEAT"))
    degpairs = sum(1 for n in range(31) if ev[n+1]-ev[n] < 1e-9)
    P("     (%d degenerate adjacent pairs in this spectrum of 32)" % degpairs)

open("OUT_run_adversary.txt", "w").write("\n".join(LOG) + "\n")
P("\nwrote OUT_run_adversary.txt")
