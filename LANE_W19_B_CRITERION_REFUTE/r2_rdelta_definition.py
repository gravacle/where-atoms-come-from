"""
r2 -- IS THE LANE'S R_delta THE PUBLISHED R_delta?

THE THREE STATEMENTS, SIDE BY SIDE.

 (Z) PUBLISHED (Blume-Kohout & Zurek 2006; Zurek, Nature Phys. 5, 181 (2009)).  Build the PARTIAL
     INFORMATION PLOT: for each fragment SIZE m, the AVERAGE of I(S:F) over all fragments of that
     size.  Let m* be the smallest size at which the average reaches (1-delta) H(S).  Then
              R_delta  =  N / m*        (equivalently 1/f_delta with f_delta = m*/N).
     R_delta is a function of (state, system, environment, delta).  NO PARTITION IS CHOSEN.

 (B) AS THE BRIEF OF THIS PROGRAM STATES IT: "partition E into DISJOINT fragments ... R_delta
     counts how many disjoint fragments independently reach it."  Computed here EXACTLY, as the
     maximum number of pairwise-disjoint fragments that clear the bar (a set-packing DP, not a
     floor(N/m*) proxy -- the proxy over-reports and is not used).

 (L) AS THE TARGET LANE COMPUTES IT (b1_discrimination.py: `single = [mi_ext(...) for e in ENV]`
     then `Rdelta = sum(1 for v in single if v >= (1-DELTA)*HS)`; b3 the same over ENVF; b4 PART 3
     the same over plaquette fragments):  fix ONE partition by hand, count its members that clear
     the bar.  R_delta is then a function of (state, cut, delta, AND THE CHOSEN PARTITION).

(L) is not (Z) and is not (B).  This file measures how much that costs, and whether the target's
findings survive being re-scored under the other two.

ISOLATION LEDGER
  held fixed : carrier, state, algebra rule, system region, delta = 0.1, and the entropy code
               (lib_b, imported unmodified from the target lane).
  moved      : PART B -- THE ENVIRONMENT PARTITION, and nothing else.
               PART C -- THE STATE, against a Haar control on the same carrier and cut.
               PART D -- THE SCORING RULE, and nothing else.
"""
import itertools
import numpy as np
from rlib import *

np.set_printoptions(precision=6, suppress=True)
DELTA = 0.1

hr("r2  R_delta: THE LANE'S COUNT versus THE PUBLISHED REDUNDANCY")

# =====================================================================================
hr("PART B  --  ONE STATE, ONE CARRIER, ONE ALGEBRA, ONE delta.  ONLY THE PARTITION MOVES.")
print("""
  Carrier theta_10.  State: the Wilson-loop GHZ of b3 PART C (the target's own BRD_W).
  System region S = the plaquette {0,1}; system algebra = the target's own declared pointer
  alg{W_01} = alg{Z_0 Z_1} (finding B-10's repair).  Fragment algebras = C3_FULL, the target's own.
  delta = 0.1.  THE ONLY THING THAT DIFFERS BETWEEN THE TWO ARMS BELOW IS HOW E IS PARTITIONED.""")
L = 10; car = theta(L)
a = sym_basis_state(car, 0)
b = sym_basis_state(car, sum(1 << i for i in range(1, L, 2)))
psiW = project_physical((a + b) / np.linalg.norm(a + b), car)
E = pauli_table(psiW, L, "direct")
POINTER = [SP(0, 0b11)]
HW = algebra_entropy(E, POINTER)[0]
ENV = list(range(2, L))
print(f"  H(W_01) = {HW:.6f} bits.   environment E = links {ENV} (N = {len(ENV)} links)")

for nm, part in [("ARM 1: partition into 8 single links", [[e] for e in ENV]),
                 ("ARM 2: partition into 4 plaquettes  ", [[2, 3], [4, 5], [6, 7], [8, 9]])]:
    vals = [mutual_information(E, POINTER, gens_FULL(car, f))[0] for f in part]
    print(f"  {nm}  I(W_01:F_k) = {['%.6f' % v for v in vals]}")
    print(f"  {'':<36}  lane R_delta = {rdelta_count(vals, HW, DELTA)} of {len(part)}")
print("""
  SAME STATE. SAME CARRIER. SAME ALGEBRA RULE. SAME delta.  R_delta = 0 or R_delta = 4 according
  to a partition that appears NOWHERE in the target's PUBLISHED_CONVENTIONS, its isolation ledger,
  or its four-item next_step.  b3 PART B partitions into single links (and reports 0); b4 PART 3
  partitions into plaquettes (and reports 5).  THE PARTITION IS A FIFTH UNDECLARED CONVENTION AND
  IT IS WORTH THE WHOLE VERDICT.""")

print("\n  THE OTHER TWO DEFINITIONS ON THE SAME OBJECT (neither chooses a partition by hand):")
def pip_alg(Efun, gsys, car, ENV, maxsize):
    mean, best, qual = {}, {}, []
    for m in range(1, maxsize + 1):
        vals = []
        for F in itertools.combinations(ENV, m):
            v = mutual_information(Efun, gsys, gens_FULL(car, list(F)))[0]
            vals.append(v)
            qual.append((sum(1 << (ENV.index(x)) for x in F), v))
        mean[m] = float(np.mean(vals)); best[m] = float(np.max(vals))
    return mean, best, qual
mean, best, qual = pip_alg(E, POINTER, car, ENV, 4)
print(f"  {'size m':>8}{'mean I(W_01:F)':>18}{'max I(W_01:F)':>16}   (all C(8,m) fragments)")
for m in sorted(mean):
    print(f"  {m:>8}{mean[m]:>18.6f}{best[m]:>16.6f}")
Rz, mz = rdelta_zurek(mean, len(ENV), HW, DELTA)
bar = (1 - DELTA) * HW - 1e-9
Rb = rdelta_packing([q for q, v in qual if v >= bar], len(ENV))
print(f"  PUBLISHED (Z):  m* = {mz}, R_delta = N/m* = {Rz:.2f}")
print(f"  BRIEF'S   (B):  exact max disjoint packing = {Rb}   (fragment sizes searched: 1..4)")
print("""  Both recover redundancy on a state the lane's ARM 1 scores at zero.  The zero is an
  artefact of the hand-chosen partition -- not of the state, the algebra, delta, or the criterion.""")

# =====================================================================================
hr("PART C  --  A STATE THE LANE'S R_delta SCORES AT ZERO AND WHICH IS A PERFECT BROADCAST")
print("""
  CONSTRUCTION (electric-basis parity broadcast, on the target's own b1 carrier theta_8).
  In the X (electric) basis, link 0 carries a bit s; the three pairs {1,2} {3,4} {5,6} each carry
  a uniformly random label whose PARITY equals s; link 7 is a |+> spectator.  Total X-parity =
  s + 3s + 0 = 0 mod 2, so every basis element is already in the physical sector.
     - every SINGLE link is uniformly random and statistically independent of s;
     - every MATCHED PAIR determines s with certainty.
  THIS IS THE SHAPE A HOLONOMY RECORD MUST HAVE: a Wilson loop is a PRODUCT over >= girth links,
  never a single link.  A criterion that only ever reads single links is blind to it by
  construction -- which is the situation of b3 PART B, whose verdict is finding B-7.""")
L = 8; car = theta(L); ENV = list(range(1, L))
def parity_broadcast(car, sysl, pairs, spectators=()):
    Lc = car["L"]; phi = np.zeros(1 << Lc, dtype=complex)
    for s in (0, 1):
        for ch in itertools.product((0, 1), repeat=len(pairs)):
            lab = s << sysl
            for k, (u, v) in enumerate(pairs):
                lab |= (ch[k] << u) | ((ch[k] ^ s) << v)
            if bin(lab).count("1") % 2: continue          # keep only the physical (even-X) sector
            phi[lab] += 1.0
    phi /= np.linalg.norm(phi)
    Hd = np.array([[1, 1], [1, -1]]) / np.sqrt(2); U = np.array([[1.0]])
    for _ in range(Lc): U = np.kron(U, Hd)
    return U @ phi
psiP = parity_broadcast(car, 0, [(1, 2), (3, 4), (5, 6)])
print(f"  ||psi - P_phys psi|| = {np.linalg.norm(psiP - project_physical(psiP, car)):.3e}"
      f"   (0 => already physical; no projection applied)")

def pip_ext(psi, Lc, S, ENV, maxsize):
    mean, best, qual = {}, {}, []
    for m in range(1, maxsize + 1):
        vals = []
        for F in itertools.combinations(ENV, m):
            v = mi_ext(psi, Lc, S, list(F)); vals.append(v)
            qual.append((sum(1 << ENV.index(x) for x in F), v))
        mean[m] = float(np.mean(vals)); best[m] = float(np.max(vals))
    return mean, best, qual

for lbl, psi in [("PAR  parity broadcast", psiP),
                 ("SCR  Haar physical seed 101 (control, same carrier & cut)",
                  haar_physical(car, 101))]:
    HS = vn_entropy(reduce_links(psi, L, [0]))
    mn, bs, ql = pip_ext(psi, L, [0], ENV, len(ENV))
    bar = (1 - DELTA) * HS - 1e-9
    RL = rdelta_count([mi_ext(psi, L, [0], [e]) for e in ENV], HS, DELTA)
    Rz, mz = rdelta_zurek(mn, len(ENV), HS, DELTA)
    Rb = rdelta_packing([q for q, v in ql if v >= bar], len(ENV))
    print(f"\n  {lbl}      H(S) = {HS:.6f}")
    print(f"     mean I(S:F) by size 1..7 : {' '.join('%.4f' % mn[m] for m in sorted(mn))}")
    print(f"     max  I(S:F) by size 1..7 : {' '.join('%.4f' % bs[m] for m in sorted(bs))}")
    print(f"     (L) lane's count over the SINGLE-LINK partition ....... {RL} of 7")
    print(f"     (L) lane's count over the MATCHED-PAIR partition ...... "
          f"{rdelta_count([mi_ext(psi, L, [0], list(p)) for p in [(1,2),(3,4),(5,6)]], HS, DELTA)}"
          f" of 3")
    print(f"     (Z) published N/m* .................................... {Rz:.2f}  (m* = {mz})")
    print(f"     (B) exact maximum disjoint packing .................... {Rb}")
print("""
  READING, AND IT CUTS BOTH WAYS.
  On the PAR state three disjoint fragments each determine the system bit with certainty.  The
  exact disjoint packing (B) -- the statement THIS PROGRAM'S BRIEF makes -- returns 3 against 1
  for the Haar control on the same carrier and cut, so it separates them and the test had a way
  to fail.  The lane's (L) returns 0, a FALSE NEGATIVE of the implementation.
  BUT THE PUBLISHED (Z) IS ALSO BLIND HERE: 1.75 for PAR and 1.75 for the Haar control, identical
  to two decimals, because averaging over ALL fragments of a size washes out a record that lives
  in 3 fragments out of C(7,2) = 21.  So this is NOT 'the lane paraphrased and the literature is
  fine'.  A HOLONOMY-SHAPED RECORD -- one carried by a PRODUCT over >= girth links and by no
  single link -- IS INVISIBLE TO BOTH THE LANE'S COUNT AND THE PUBLISHED PARTIAL-INFORMATION PLOT,
  and visible only to the disjoint-packing reading.  That is a defect of quantum Darwinism as
  usually written when it is carried into a gauge theory, and it is not in the lane's list.""")

# =====================================================================================
hr("PART D  --  RE-SCORING THE TARGET'S OWN b1 TABLE.  ONLY THE SCORING RULE MOVES.")
print("""  Carrier theta_8, S = {link 0}, E = links 1..7, extended-Hilbert-space estimator,
  delta = 0.1 -- b1_discrimination.py exactly.  Every arm below is the target's own state.""")
def bell_localised(car):
    Lc = car["L"]; psi = np.zeros(1 << Lc, dtype=complex)
    for c in range(1 << (Lc - 2)):
        bits = sum(((c >> i) & 1) << (i + 2) for i in range(Lc - 2))
        psi[bits] += 1.0; psi[bits | 0b11] += 1.0
    return project_physical(psi / np.linalg.norm(psi), car)
ARMS = [("BRD_M magnetic GHZ", sym_basis_state(car, 0)),
        ("BRD_E electric GHZ", elec_ghz(car)),
        ("BRD_Ew elec GHZ a=0.99", elec_ghz(car, 0.99)),
        ("LOC   Bell(0,1) x |+>^6", bell_localised(car)),
        ("VAC   |+>^8", project_physical(np.ones(1 << L, dtype=complex), car)),
        ("SCR   Haar seed 101", haar_physical(car, 101)),
        ("SCR   Haar seed 202", haar_physical(car, 202)),
        ("PAR   parity broadcast [new]", psiP)]
print(f"  {'ARM':<30}{'H(S)':>10}{'(L) count':>11}{'(Z) N/m*':>10}{'m*':>5}{'(B) packing':>13}")
print("-" * 104)
for nm, psi in ARMS:
    HS = vn_entropy(reduce_links(psi, L, [0]))
    sing = [mi_ext(psi, L, [0], [e]) for e in ENV]
    mn, bs, ql = pip_ext(psi, L, [0], ENV, len(ENV))
    bar = (1 - DELTA) * HS - 1e-9
    Rz, mz = rdelta_zurek(mn, len(ENV), HS, DELTA)
    Rb = rdelta_packing([q for q, v in ql if v >= bar], len(ENV))
    print(f"  {nm:<30}{HS:>10.6f}{rdelta_count(sing,HS,DELTA):>11}{Rz:>10.2f}{str(mz):>5}{Rb:>13}")
print("""
  WHAT SURVIVES THE RE-SCORING.  B-1 (the criterion discriminates): the broadcast/scrambled
  separation is 7 vs 0 under (L), 7.00 vs 1.75 under (Z), 7 vs 1 under (B).  SURVIVES under all
  three.  B-2 (H(S) = 0 is trivially satisfied): VAC scores maximally under all three -- so it is
  a defect of the PUBLISHED criterion, correctly identified.  B-3 (blind to how much is recorded):
  BRD_Ew scores maximally under all three at H(S) = 0.14 -- also correctly identified.
  WHAT DOES NOT SURVIVE: the identification of the lane's hand-partitioned count WITH R_delta.
  The LOC row moves 1 -> 1.75 -> 1 and the PAR row moves 0 -> 1.75 -> 3.""")
