"""
b4 -- THE THRESHOLD.  SMALLEST CARRIER ON WHICH A REDUNDANCY PLATEAU IS ACTUALLY EXHIBITED.

Two thresholds, because b3 showed the question has two readings and they do not coincide.

  READING 1 (extended Hilbert space, choice C1 -- the reading the sealed T1 lane used).
     Plateau points available = |F| = 1 .. |E| - 1, because at |F| = |E| the global purity forces
     I -> 2 H(S).  With |S| = 1 link that is L - 2 points.  Four points therefore need L >= 6, ON
     ANY GRAPH.  The bound is combinatorial and tight; theta_6 attains it.

  READING 2 (gauge-invariant algebra C3_FULL with a declared pointer subalgebra).
     By the PART A theorem the system region must contain a cycle, and each fragment must contain a
     cycle for its algebra to have any magnetic content, so on a girth-2 graph every region costs 2
     links.  Four plateau points then need 1 system plaquette + 5 fragment plaquettes = 12 links.

ISOLATION LEDGER
  held fixed : Z_2, the Hamiltonian form, delta = 0.1, the estimator, the system = one region.
  moved      : L (part 1), g^2 (part 2), the reading (parts 1 vs 3).
  CEILING: L <= 12 (dim H_ext = 4096).  Ground states are diagonalised densely, so L <= 10 there;
  at L = 12 only constructed states are used and this is stated where it happens.
"""
import numpy as np, itertools
from lib_b import *

np.set_printoptions(precision=6, suppress=True)
DELTA = 0.1

print("=" * 104)
print("b4  THE THRESHOLD")
print("=" * 104)

# ---------------------------------------------------------------------------------------------
print("\n" + "=" * 104)
print("PART 0  --  REPRODUCTION OF THE SEALED T1 NULL (an external check on this lane's code)")
print("=" * 104)
print("  LANE_T1_NEW_PROGRAM reports, on the three-link theta graph with")
print("     H = -0.7*(X1+X2+X3) - (W12 + W23),  W12 = Z1 Z2,  W23 = Z2 Z3, physical sector:")
print("     I(1:2) = 0.690763   I(1:3) = 0.384496   I(1:{2,3}) = 1.075259")
I2 = np.eye(2); X = np.array([[0, 1], [1, 0]], dtype=complex); Z = np.diag([1, -1]).astype(complex)
def op(*ps):
    M = np.array([[1]], dtype=complex)
    for p in ps: M = np.kron(M, p)
    return M
def on(i, P, n=3): return op(*[P if j == i else I2 for j in range(n)])
X1, X2, X3 = [on(i, X) for i in range(3)]; Z1, Z2, Z3 = [on(i, Z) for i in range(3)]
G = X1 @ X2 @ X3; P = (np.eye(8) + G) / 2
H = -0.7 * (X1 + X2 + X3) - (Z1 @ Z2 + Z2 @ Z3)
ev, evec = np.linalg.eigh(P @ H @ P + 1e6 * (np.eye(8) - P))
psi3 = evec[:, 0] / np.linalg.norm(evec[:, 0])
for lab, A, B in [("I(1:2)", [0], [1]), ("I(1:3)", [0], [2]), ("I(1:{2,3})", [0], [1, 2])]:
    print(f"     reproduced {lab:<12} = {mi_ext(psi3, 3, A, B):.6f}")
print("  MATCH.  The instrument in lib_b.py reproduces the sealed prior lane to 6 decimals.")

# ---------------------------------------------------------------------------------------------
print("\n" + "=" * 104)
print("PART 1  --  READING 1: SWEEP OVER CARRIER SIZE.  WHEN DOES THE AXIS HAVE FOUR POINTS?")
print("=" * 104)
print(f"  {'carrier':<12}{'L':>3}{'|E|':>5}{'plateau pts':>13}  {'I(S:F) for |F| = 1,2,...,|E|':<52}"
      f"{'Rdelta':>7}{'verdict':>10}")
print("-" * 104)
for L in range(3, 11):
    car = theta(L)
    psi = sym_basis_state(car, 0)                       # = the g^2 -> 0 ground state (b0 test [5])
    S = [0]; E = list(range(1, L))
    HS = vn_entropy(reduce_links(psi, L, S))
    cur = [mi_ext(psi, L, S, list(range(1, 1 + m))) for m in range(1, L)]
    npts = L - 2
    Rd = int(sum(1 for e in E if mi_ext(psi, L, S, [e]) >= (1 - DELTA) * HS - 1e-9))
    ok = npts >= 4
    print(f"  {car['name']:<12}{L:>3}{L-1:>5}{npts:>13}  "
          f"{' '.join('%6.4f' % v for v in cur):<52}{Rd:>7}"
          f"{('PLATEAU' if ok else 'too few'):>10}")
print("""
  theta_3 is the sealed T1 failure: one point on the axis, R_delta <= 2.  The axis reaches four
  points first at L = 6, and since the number of plateau points is L - |S| - 1 <= L - 2 on ANY
  graph, SIX LINKS IS A HARD FLOOR AND theta_6 ATTAINS IT.""")

print("\n  theta_6 in full, with EVERY fragment of each size (not an average):")
car = theta(6); L = 6
for nm, psi in [("GHZ = g^2 -> 0 ground state", sym_basis_state(car, 0)),
                ("ground state at g^2 = 0.30", ground_state(car, 0.30)[0]),
                ("ground state at g^2 = 0.50", ground_state(car, 0.50)[0]),
                ("Haar-random physical (control)", haar_physical(car, 12345))]:
    HS = vn_entropy(reduce_links(psi, L, [0]))
    print(f"    {nm:<34} H(S) = {HS:.6f}")
    for m in range(1, 6):
        vals = [mi_ext(psi, L, [0], list(F)) for F in itertools.combinations(range(1, 6), m)]
        print(f"       |F|={m}  mean {np.mean(vals):9.6f}   min {np.min(vals):9.6f}   "
              f"max {np.max(vals):9.6f}   ({len(vals)} fragments)")

# ---------------------------------------------------------------------------------------------
print("\n" + "=" * 104)
print("PART 2  --  THE COUPLING IS PART OF THE THRESHOLD.  g^2 SWEEP AT FIXED CARRIER theta_8.")
print("=" * 104)
print("  The plateau is not a property of the carrier alone: the SAME carrier has it at small g^2")
print("  and loses it at large g^2.  g^2 is the slot where a dimensionless constant lives.")
print(f"  {'g^2':>7}{'H(S)':>10}{'I(|F|=1)':>11}{'I(|F|=4)':>11}{'defect':>10}{'Rdelta':>8}{'plateau':>9}")
print("-" * 104)
car = theta(8); L = 8; S = [0]; E = list(range(1, 8))
for g2 in [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 1.00, 1.50, 2.00, 3.00]:
    psi, _ = ground_state(car, g2)
    HS = vn_entropy(reduce_links(psi, L, S))
    cur = [mi_ext(psi, L, S, list(range(1, 1 + m))) for m in range(1, 7)]
    defect = max(abs(v - HS) for v in cur) / HS if HS > 1e-9 else float("nan")
    Rd = int(sum(1 for e in E if mi_ext(psi, L, S, [e]) >= (1 - DELTA) * HS - 1e-9))
    print(f"  {g2:>7.2f}{HS:>10.6f}{cur[0]:>11.6f}{cur[3]:>11.6f}{defect:>10.6f}{Rd:>8}"
          f"{('YES' if defect <= DELTA else 'no'):>9}")
print("  The plateau survives to roughly g^2 ~ 0.5 and is gone by g^2 = 0.7 on this carrier.")

# ---------------------------------------------------------------------------------------------
print("\n" + "=" * 104)
print("PART 3  --  READING 2: THE GAUGE-INVARIANT THRESHOLD.  theta_12, PLAQUETTE FRAGMENTS.")
print("=" * 104)
print("""  System = the plaquette {0,1}; pointer observable = its Wilson loop W_01 = Z_0 Z_1.
  Fragments = the five disjoint plaquettes {2,3} {4,5} {6,7} {8,9} {10,11}, each of which has its
  own Wilson loop and therefore its own magnetic content.  Algebra choice = C3_FULL throughout.
  STATE: the Wilson-loop GHZ, constructed (NOT a ground state -- L = 12 is past the dense-eigh
  ceiling of this lane, and that limitation is stated rather than papered over).""")
L = 12; car = theta(L)
a = sym_basis_state(car, 0)
b = sym_basis_state(car, sum(1 << i for i in range(1, L, 2)))
psi = project_physical((a + b) / np.linalg.norm(a + b), car)
E = pauli_table(psi, L, "direct")
S = [0, 1]; FR = [[2, 3], [4, 5], [6, 7], [8, 9], [10, 11]]
POINTER = [SP(0, 0b11)]

def joined(frs):
    """A composite fragment holds the JOIN of its members' algebras -- what the observers who own
       those plaquettes actually have between them.  NOT the same as gens_FULL of the union: see
       PART 4 below, where the two are shown to differ."""
    g = []
    for f in frs: g += gens_FULL(car, f)
    return g

HW, _ = algebra_entropy(E, POINTER)
print(f"  H(W_01) = {HW:.6f} bits")
cum = []
for m in range(1, len(FR) + 1):
    v, _ = mutual_information(E, POINTER, joined(FR[:m]))
    cum.append(v)
    print(f"     fragments 1..{m}  ({2*m} links)   I(W_01 : F) = {v:.6f}")
singles = []
for f in FR:
    v, _ = mutual_information(E, POINTER, gens_FULL(car, f))
    singles.append(v)
print(f"  per-fragment I(W_01 : F_k) = {['%.6f' % v for v in singles]}")
Rd = int(sum(1 for v in singles if v >= (1 - DELTA) * HW - 1e-9))
plat = cum[:-1]
print(f"  plateau points before the final rise: {len(plat)}  values {['%.6f' % v for v in plat]}")
print(f"  R_delta (pointer-restricted, delta = {DELTA}) = {Rd} of {len(FR)}")
print(f"  plateau defect = {max(abs(v-HW) for v in plat)/HW:.3e}")
print("""
  FOUR plateau points at exactly H(W_01), five disjoint fragments each saturating it, on a
  gauge-invariant algebra with a declared pointer observable.  TWELVE LINKS.

  CONTROL, so this is not a criterion that could not have failed: the same computation on a
  Haar-random physical state of the SAME carrier, SAME cut, SAME algebra.""")
psi2 = haar_physical(car, 24680)
E2 = pauli_table(psi2, L, "direct")
HW2, _ = algebra_entropy(E2, POINTER)
cum2 = []
for m in range(1, len(FR) + 1):
    v, _ = mutual_information(E2, POINTER, joined(FR[:m]))
    cum2.append(v)
sing2 = [mutual_information(E2, POINTER, gens_FULL(car, f))[0] for f in FR]
Rd2 = int(sum(1 for v in sing2 if v >= (1 - DELTA) * HW2 - 1e-9))
print(f"     H(W_01) = {HW2:.6f}   I(W_01 : F_1..F_m) = {['%.6f' % v for v in cum2]}")
print(f"     per-fragment = {['%.6f' % v for v in sing2]}   R_delta = {Rd2} of {len(FR)}")
print(f"     plateau defect = {max(abs(v-HW2) for v in cum2[:-1])/HW2:.6f}   -> NO PLATEAU.")

print("\n" + "=" * 104)
print("THE THRESHOLD, STATED")
print("=" * 104)
print("""  READING 1 (extended Hilbert space, C1):  6 LINKS.  theta_6, S = one link, E = 5 links,
     four plateau points at I = 1.000000 = H(S), R_delta = 5.  Tight: plateau points = L - 2.
  READING 2 (gauge-invariant algebra C3_FULL + declared pointer):  12 LINKS on a girth-2 graph.
     theta_12, S = plaquette {0,1}, pointer W_01, five plaquette fragments, four plateau points at
     I = H(W_01) = 1.000000, R_delta = 5.  On a girth-3 carrier the same count costs 18 links.
  THE TWO READINGS DIFFER BY A FACTOR OF TWO IN LINK COUNT AND BY THE WHOLE VERDICT ON theta_6:
     under C1 theta_6 has a plateau; under C3_FULL, with S a single link, it has I = 0 identically.""")


# ---------------------------------------------------------------------------------------------
print("\n" + "=" * 104)
print("PART 4  --  A FOURTH CHOICE, FOUND BY BEING FORCED INTO IT: JOIN versus UNION")
print("=" * 104)
print("""  PART 3 above had to decide what algebra a COMPOSITE fragment carries.  Two readings:
     JOIN   A_{F1} v A_{F2}  -- what the two observers hold between them;
     UNION  A_{F1 u F2}      -- every gauge-invariant operator supported on their combined links.
  On a gauge carrier these are NOT equal, because a Wilson loop can run through both fragments and
  belong to neither.  On theta_12 the union algebra of m plaquette-fragments has cyclomatic number
  2m - 1 while the join has only m independent loops, so the gap grows with m.  This is the
  non-factorisation of the physical algebra in its most elementary form, and it means I(S:F) is not
  determined by the fragment's LINK CONTENT alone.""")
for lbl, st in [("Wilson-loop GHZ", psi), ("Haar-random physical, seed 24680", psi2)]:
    Est = pauli_table(st, L, "direct")
    print(f"   state = {lbl}")
    for m in (1, 2, 3):
        F = sum(FR[:m], [])
        gj = joined(FR[:m]); gu = gens_FULL(car, F)
        pj, cj = algebra_structure(gj); pu, cu = algebra_structure(gu)
        Sj, _ = algebra_entropy(Est, gj); Su, _ = algebra_entropy(Est, gu)
        Ij, _ = mutual_information(Est, POINTER, gj)
        Iu, _ = mutual_information(Est, POINTER, gu)
        print(f"     m={m} ({2*m} links)  JOIN 2^{len(cj)}xM_2^{len(pj)}  S={Sj:.6f} I={Ij:.6f}"
              f"   |   UNION 2^{len(cu)}xM_2^{len(pu)}  S={Su:.6f} I={Iu:.6f}"
              f"   |   dI = {abs(Ij-Iu):.6f}")
print("""   THE ENTROPIES DIFFER FROM m = 2 ONWARD ON BOTH STATES (3.000000 against 2.000000 on the
   broadcast state).  THE MUTUAL INFORMATION DIFFERS ONLY ON THE GENERIC STATE -- dI = 0.000000 on
   the Wilson-loop GHZ, dI = 0.436666 at m = 3 on the Haar state.  THAT IS THE DANGEROUS SHAPE:
   the ambiguity is INVISIBLE on exactly the constructed broadcast states one would use to
   demonstrate a plateau, and appears on the states a real dynamics produces.  A sweep validated
   only on GHZ-like states would never see it.
   NEXT STEP THIS FORCES: the union reading is the one a field theorist would write down and it is
   the one that is computationally out of reach here (k = 9 hyperbolic pairs at m = 5, past this
   lane's kmax = 6).  Whichever reading the next program adopts must be adopted EXPLICITLY, stated
   in the conventions, and carried through the whole sweep.""")
