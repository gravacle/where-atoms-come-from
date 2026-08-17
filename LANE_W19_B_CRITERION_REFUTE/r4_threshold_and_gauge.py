"""
r4 -- THE THRESHOLD, AND HOW MUCH GAUGE THEORY IS IN IT.

Four questions, in order of how much they cost the target lane's headline.

  PART 1  IS THE GAUSS LAW DOING ANY WORK ON THE STATES THIS LANE USES?
          A theorem (Perron-Frobenius) plus a measurement: for H = -(1/g^2) sum_p W_p - g^2 sum_l X_l
          on ANY of these carriers, the UNCONSTRAINED ground state is already in the physical
          sector, so lib_b.ground_state's projection is inert.  Combined with the fact that the
          READING-1 estimator (mi_ext) never consults the graph, every reading-1 number in b1,
          b3 PART B and b4 PARTS 1-2 is a computation with ZERO GAUGE INPUT.
          CONTROL THAT CAN FAIL: flip the sign of the electric term on an odd-L carrier.

  PART 2  WHAT CARRIES THE READING-1 PLATEAU?  Decompose C1's one bit of I(S:F) into the
          gauge-invariant part alg{X_0} and the gauge-VARIANT part alg{Z_0}.  This adjudicates the
          question the lane's last self-flag says it deliberately did not adjudicate.

  PART 3  IS THE THRESHOLD TIGHT?  B-12 is proved only for PURE global states and the lane records
          the caveat as untested.  I test it.  Reading 1's floor falls from 6 to 5 only for a
          gauge-invariant state that is NOT in the physical sector, so 6 stands.  Reading 2's 12
          falls to 10 once a spurious purity discard is removed.

  PART 4  AND 10 IS STILL WRONG.  Reading 2's threshold is SIX links.  b4 PART 3's derivation
          requires every FRAGMENT to contain a cycle; that does not follow from B-8, which is about
          where a magnetic OPERATOR can be localised, not where information about one can be
          recorded.  A constructed physical stabiliser state puts a full bit of holonomy entropy
          into four PURELY ELECTRIC single-link fragments.  PART 3's own 10 is left in place, with
          a warning, because it repeats half of the same error.

ISOLATION LEDGER
  held fixed : Z_2, the Hamiltonian form, delta = 0.1, the estimator, the system region, and the
               entropy code (lib_b, unmodified).
  moved      : PART 1 -- whether the Gauss projection is applied (and, in the control, the SIGN of
               the electric term).  PART 2 -- the system ALGEBRA.  PART 3 -- global PURITY only:
               same carrier, same cut, same pointer, same delta.
  CEILING    : dense eigh to L = 10 (2^10); constructed and mixed states to L = 12; algebra
               entropies to k = 6 hyperbolic pairs.
"""
import itertools
import numpy as np
from rlib import *

np.set_printoptions(precision=6, suppress=True)
DELTA = 0.1

hr("r4  THE THRESHOLD, AND HOW MUCH GAUGE THEORY IS IN IT")

# =====================================================================================
hr("PART 1  --  IS THE GAUSS LAW DOING ANY WORK ON THE GROUND STATES THIS LANE USES?")
print("""
  THEOREM.  Let H = -(1/g^2) sum_p W_p - g^2 sum_l X_l with g^2 > 0 on any of these carriers.
  In the Z (link-configuration) basis, H has a diagonal part and OFF-DIAGONAL ENTRIES ALL EQUAL
  TO -g^2 < 0 (the plaquette term is diagonal for Z_2).  So cI - H has non-negative entries and is
  irreducible, and by Perron-Frobenius the ground state of H is UNIQUE with STRICTLY POSITIVE
  entries.  Every Gauss operator G_v is a permutation of basis states and commutes with H, so the
  ground state is a G_v eigenvector; a strictly positive vector cannot have eigenvalue -1 under a
  permutation.  Hence the unconstrained ground state is ALREADY PHYSICAL, for every v, every
  carrier, every g^2 > 0.  THE PROJECTION IN lib_b.ground_state IS INERT.
  MEASURED: |<psi_free | psi_projected>| and the Gauss expectations of the UNPROJECTED state.""")
print(f"  {'carrier':<10}{'g^2':>6}{'|<free|proj>|':>15}{'min <G_v> (free)':>19}"
      f"{'||psi_free - P psi_free||':>27}")
print("-" * 104)
for car in [theta(5), theta(6), theta(7), theta(8), fan(3), fan(4), wheel(4)]:
    Lc = car["L"]
    for g2 in (0.30, 1.00, 3.00):
        H = hamiltonian(car, g2)
        w, v = np.linalg.eigh(H)
        free = v[:, 0]
        proj = ground_state(car, g2)[0]
        Ef = pauli_table(free, Lc, "direct")
        gmin = min(sp_expect(Ef, SP(g, 0)).real for g in indep_gauss(car))
        print(f"  {car['name']:<10}{g2:>6.2f}{abs(np.vdot(free, proj)):>15.12f}{gmin:>19.12f}"
              f"{np.linalg.norm(free - project_physical(free, car)):>27.3e}")
print("""
  CONTROL THAT CAN FAIL, AND DOES.  Flip the SIGN of the electric term: H' = -(1/g^2) sum_p W_p
  + g^2 sum_l X_l.  Now the off-diagonals are POSITIVE, Perron-Frobenius applies to the TOP of the
  spectrum instead, and on an ODD number of links the ground state lands in the G = -1 sector.""")
print(f"  {'carrier':<10}{'L parity':>10}{'g^2':>6}{'<G> of free ground state':>28}{'physical?':>12}")
for car in [theta(5), theta(7), theta(6), theta(8)]:
    Lc = car["L"]
    for g2 in (1.00,):
        Hp = hamiltonian(car, g2)
        idx = np.arange(1 << Lc)
        for l in range(Lc):                      # undo the -g^2 and add +g^2
            Hp[idx ^ (1 << l), idx] += 2 * g2
        w, v = np.linalg.eigh(Hp)
        free = v[:, 0]
        Ef = pauli_table(free, Lc, "direct")
        gexp = min(sp_expect(Ef, SP(g, 0)).real for g in indep_gauss(car))
        print(f"  {car['name']:<10}{'odd' if Lc % 2 else 'even':>10}{g2:>6.2f}{gexp:>28.12f}"
              f"{str(gexp > 0.5):>12}")
print("""
  CONSEQUENCE.  Under READING 1 the estimator is lib_b.mi_ext, which takes (psi, L, A, B) and
  never touches car['gauss'] or car['plaq'].  The state it is fed is, by the theorem above, the
  state the SAME Hamiltonian would have without any Gauss law.  So b1 in its entirety, b3 PART B's
  C1 column, and b4 PARTS 1 and 2 are computations in which the gauge structure enters NOWHERE --
  neither the state nor the estimator.  They are redundancy measurements on the ground state of a
  transverse-field Ising model on a star (theta) or on fan/wheel graphs.  That does not make the
  numbers wrong.  It makes 'SIX LINKS is the smallest GAUGE CARRIER' an overstatement: six is the
  smallest carrier of ANY kind, and the demonstration has ZERO VARIABLES MOVED with respect to
  gauge invariance.  Shown directly:""")
print(f"  {'L':>3}{'theta_L, GHZ = g^2->0 GS':>28}{'L bare qubits, GHZ, no Gauss law':>36}{'identical?':>12}")
for L in range(3, 9):
    car = theta(L)
    p1 = sym_basis_state(car, 0)
    p2 = np.zeros(1 << L, dtype=complex); p2[0] = p2[(1 << L) - 1] = 1 / np.sqrt(2)
    c1 = [mi_ext(p1, L, [0], list(range(1, 1 + m))) for m in range(1, L)]
    c2 = [mi_ext(p2, L, [0], list(range(1, 1 + m))) for m in range(1, L)]
    print(f"  {L:>3}{' '.join('%5.3f' % v for v in c1):>28}{' '.join('%5.3f' % v for v in c2):>36}"
          f"{str(np.allclose(c1, c2)):>12}")

# =====================================================================================
hr("PART 2  --  WHAT CARRIES THE READING-1 PLATEAU?  THE LANE DECLINED TO ADJUDICATE; THIS DOES.")
print("""
  The lane's last self-flag: 'either the extended-Hilbert-space plateau is a real record of a
  variable no gauge-invariant measurement can read ... or the null means the magnetic record is
  simply not localised in one link ... no number here can settle it.'
  A number can.  A_S under C1 is the full 2x2 matrix algebra on link 0.  Split it into subalgebras
  and ask which one carries the bit.  Same state, same carrier, same cut, same code; only the
  SYSTEM ALGEBRA and the matching FRAGMENT ALGEBRA move, and both are printed.""")
L = 8; car = theta(L); psi = sym_basis_state(car, 0); E = pauli_table(psi, L)
print(f"  {'A_S':<16}{'A_F (link 1)':<16}{'gauge-inv?':<12}{'S_A(S)':>9}{'S_A(F)':>9}{'I(S:F)':>9}")
print("-" * 104)
tests = [("alg{X_0}", [SP(1, 0)], "alg{X_1}", [SP(2, 0)]),
         ("alg{Z_0}", [SP(0, 1)], "alg{Z_1}", [SP(0, 2)]),
         ("alg{Y_0}", [SP(1, 1)], "alg{Y_1}", [SP(2, 2)]),
         ("B(H_0) = C1", gens_EXT(car, [0]), "B(H_1) = C1", gens_EXT(car, [1]))]
for na, ga, nb, gb in tests:
    gi = is_gauge_invariant(car, ga) and is_gauge_invariant(car, gb)
    Sa = algebra_entropy(E, ga)[0]; Sb = algebra_entropy(E, gb)[0]
    I = mutual_information(E, ga, gb)[0]
    print(f"  {na:<16}{nb:<16}{str(gi):<12}{Sa:>9.6f}{Sb:>9.6f}{I:>9.6f}")
print("""
  ADJUDICATED.  The whole of C1's one bit is the Z_0-Z_1 correlation.  alg{Z_0} is NOT gauge
  invariant -- Z_0 anticommutes with G = prod_l X_l -- and the gauge-invariant part alg{X_0}
  carries EXACTLY ZERO.  So the reading-1 plateau on theta_L is the redundancy of an operator that
  no measurement in this theory can perform.  B-7's two readings are NOT symmetric: one of them
  is measuring a quantity that is not an observable.  That is a disposition the lane could have
  reached with three lines of its own instrument and did not.
  THE HONEST RESIDUE, WHICH IS THE LANE'S REAL POINT AND SURVIVES: the gauge-invariant null is
  about the REGION, not about the record.  b3's own PART C shows it -- enlarge S to a plaquette
  and the gauge-invariant algebra sees a flat I = 1.000000.""")

# =====================================================================================
hr("PART 3  --  IS THE THRESHOLD TIGHT?  B-12's PURITY CAVEAT, TESTED.")
print("""
  B-12: 'plateau points = L - |S| - 1 <= L - 2 for a pure global state on ANY graph, so four points
  force L >= 6'.  Recorded caveat: 'a mixed global state (no final jump) could reach four points at
  L = 5 -- untested here'.  Tested here, on both readings.

  3a. READING 1 at L = 5.  rho = 1/2(|0^5><0^5| + |1^5><1^5|), the classical GHZ.  ONLY PURITY
      MOVES: same carrier, same cut, same estimator, same delta.""")
L = 5; car = theta(L)
psis = [np.zeros(1 << L, dtype=complex), np.zeros(1 << L, dtype=complex)]
psis[0][0] = 1.0; psis[1][(1 << L) - 1] = 1.0
rho = rho_from_mixture(psis, [0.5, 0.5])
P = np.eye(1 << L); idx = np.arange(1 << L)
for g in indep_gauss(car):
    Pg = np.zeros((1 << L, 1 << L)); Pg[idx, idx] += 0.5; Pg[idx ^ g, idx] += 0.5
    P = P @ Pg
g0 = indep_gauss(car)[0]
Gperm = np.zeros((1 << L, 1 << L)); Gperm[idx ^ g0, idx] = 1.0    # G = prod_l X_l as a permutation
print(f"      gauge invariant?  ||G rho G - rho|| = "
      f"{np.linalg.norm(Gperm @ rho @ Gperm - rho):.3e}")
print(f"      IN THE PHYSICAL SECTOR?  Tr(rho P_+) = {np.trace(P @ rho).real:.6f}   <-- NOT 1")
HS = vn_entropy(reduce_rho(rho, L, [0]))
print(f"      H(S) = {HS:.6f}")
for m in range(1, 5):
    vals = [mi_ext_rho(rho, L, [0], list(F)) for F in itertools.combinations(range(1, 5), m)]
    print(f"      |F| = {m}   I(S:F) mean {np.mean(vals):.6f}  min {np.min(vals):.6f}  "
          f"max {np.max(vals):.6f}   ({len(vals)} fragments)")
print("""      => FOUR flat points at I = H(S) = 1 on FIVE links.  B-12's caveat is real.
      BUT I MUST SAY WHAT THIS STATE IS: it is gauge-invariant as an operator (G rho G = rho) and
      it is an EQUAL MIXTURE OF THE TWO GAUSS SECTORS, Tr(rho P_+) = 0.5, so it is not a state of
      the physical Hilbert space.  For READING 1 the correction is therefore CONDITIONAL: the
      bound is 5 for gauge-invariant density matrices on H_ext and I have NOT exhibited a
      physical-sector state achieving 5.  A physical-sector thermal state is tested next.""")

print("""
  3b. READING 1 at L = 5, PHYSICAL-SECTOR GIBBS STATE -- the dynamical version of the same test.
      rho_beta = P exp(-beta H) P / Z at g^2 = 0.30.  Same carrier, same cut; beta moves.""")
Hm = hamiltonian(car, 0.30)
w, v = np.linalg.eigh(Hm)
print(f"      {'beta':>7}{'H(S)':>10}{'I |F|=1':>10}{'I |F|=2':>10}{'I |F|=3':>10}{'I |F|=4':>10}"
      f"{'4 pts within 10%?':>20}")
for beta in (0.5, 1.0, 2.0, 4.0, 8.0, 20.0):
    ex = v @ np.diag(np.exp(-beta * (w - w[0]))) @ v.conj().T
    r = P @ ex @ P
    r = r / np.trace(r).real
    HSb = vn_entropy(reduce_rho(r, L, [0]))
    cur = [np.mean([mi_ext_rho(r, L, [0], list(F))
                    for F in itertools.combinations(range(1, 5), m)]) for m in range(1, 5)]
    ok = HSb > 1e-9 and all(abs(c - HSb) <= DELTA * HSb for c in cur)
    print(f"      {beta:>7.1f}{HSb:>10.6f}" + "".join(f"{c:>10.6f}" for c in cur) + f"{str(ok):>20}")
print("""      => no beta gives four points inside 10% on theta_5: at high beta the state is the
      pure GHZ and the last point jumps to 2H(S); at low beta the plateau height itself decays.
      SO: READING 1's SIX-LINK FLOOR STANDS FOR PHYSICAL-SECTOR STATES OF THIS THEORY, and falls
      to five only if one admits gauge-invariant mixtures across Gauss sectors.  That is a
      correction to the word 'tight', not to the number.""")

print("""
  3c. READING 2's THRESHOLD IS TEN LINKS, NOT TWELVE -- AND THE TARGET'S OWN SEALED OUTPUT
      ALREADY SHOWS IT.  b4 PART 3 announces TWELVE links after computing five values of
      I(W_01 : F_1..F_m) that are ALL 1.000000, then reporting `cum[:-1]` under the label
      'plateau points before the final rise: 4'.  THERE IS NO FINAL RISE.  The purity jump
      I -> 2H(S) is a theorem about a FULL FACTOR and its commutant; with the pointer restricted
      to the one-bit algebra alg{W_01} and the fragments restricted to C3_FULL, neither is a full
      factor and no jump occurs.  b4 imported reading 1's discard rule into a setting where it
      does not hold, and threw away a good point.
      The target's OWN b3 PART C2 line, on theta_10, reads
          BRD_W ... H(W_01) = 1.00000   I(W_01 : F) = 1.00000 1.00000 1.00000 1.00000  Rdelta = 4
      -- four flat points at TEN links.  b4 PART 3 and b3 PART C2 are inconsistent, and PART C2
      is the correct one.  Reproduced below as ARM 1, with a mixed-state ARM 2 to show the answer
      is not a purity artefact in either direction.  ONLY PURITY MOVES between the two arms.""")
L = 10; car = theta(L)
pa = sym_basis_state(car, 0)
pb = sym_basis_state(car, sum(1 << i for i in range(1, L, 2)))
pure = project_physical((pa + pb) / np.linalg.norm(pa + pb), car)
POINTER = [SP(0, 0b11)]
FR = [[2, 3], [4, 5], [6, 7], [8, 9]]
for lbl, Eobj in [("ARM 1  PURE  (|psi_a>+|psi_b>)/sqrt2 -- the target's BRD_W",
                   pauli_table(pure, L, "direct")),
                  ("ARM 2  MIXED 1/2(|psi_a><psi_a| + |psi_b><psi_b|), both branches physical",
                   MixE([pa, pb], [0.5, 0.5], L, "direct"))]:
    HW = algebra_entropy(Eobj, POINTER)[0]
    cum = [mutual_information(Eobj, POINTER,
                              sum([gens_FULL(car, f) for f in FR[:m]], []))[0]
           for m in range(1, 5)]
    sing = [mutual_information(Eobj, POINTER, gens_FULL(car, f))[0] for f in FR]
    R = rdelta_count(sing, HW, DELTA)
    flat = [c for c in cum if abs(c - HW) <= DELTA * HW]
    print(f"  {lbl}")
    print(f"     H(W_01) = {HW:.6f}   I(W_01 : F_1..F_m) for m = 1,2,3,4 = "
          f"{['%.6f' % c for c in cum]}")
    print(f"     per-fragment = {['%.6f' % s for s in sing]}   R_delta = {R} of 4   "
          f"plateau points within 10% = {len(flat)}")
print("""
  => FOUR plateau points at I = H(W_01) = 1.000000, R_delta = 4 of 4, on TEN links, in the
  gauge-invariant algebra with the declared pointer -- for the PURE state as well as the mixed one.
  READING 2's THRESHOLD IS AT MOST 10, NOT 12, AND THE CORRECTION HAS NOTHING TO DO WITH PURITY.
  *** WARNING TO THE READER: 10 IS STILL WRONG.  PART 4 BELOW TAKES IT TO 6 AND SHOWS THAT THIS
  *** SECTION REPEATS HALF OF THE TARGET'S OWN LOGICAL ERROR.  THE WRONG STEP IS LEFT IN PLACE.

  AND 10 WOULD BE TIGHT, IF the lane's own reason for the fragment cost held: four plateau points
  need four fragments; by theorem B-8 every fragment must contain a cycle to have any magnetic
  content, so each costs >= girth links, and so must the system region -- 2 + 4*2 = 10 on girth 2.
  THAT REASON IS FALSE AND PART 4 BREAKS IT.  A fragment does not need magnetic content of its own
  to be correlated with a holonomy.  This paragraph is the target's error reproduced, and it is
  left standing so the correction in PART 4 has something to correct.
  THE NEGATIVE SIDE OF TIGHTNESS -- theta_8, three plaquette fragments, everything else identical:""")
car8 = theta(8)
pa8 = sym_basis_state(car8, 0)
pb8 = sym_basis_state(car8, sum(1 << i for i in range(1, 8, 2)))
pure8 = project_physical((pa8 + pb8) / np.linalg.norm(pa8 + pb8), car8)
E8 = pauli_table(pure8, 8, "direct")
FR8 = [[2, 3], [4, 5], [6, 7]]
HW8 = algebra_entropy(E8, POINTER)[0]
cum8 = [mutual_information(E8, POINTER, sum([gens_FULL(car8, f) for f in FR8[:m]], []))[0]
        for m in range(1, 4)]
print(f"     theta_8, S = {{0,1}}, pointer W_01, 3 plaquette fragments:  H(W_01) = {HW8:.6f}   "
      f"I = {['%.6f' % c for c in cum8]}")
print(f"     -> only {len(cum8)} points on the fragment-count axis.  A four-point plateau is not"
      f" available at L = 8.")
print("""
  CONTROL, SAME CARRIER / CUT / ALGEBRA / delta, only the state moved -- a Haar physical state and
  a classical mixture of two Haar physical states, on theta_10:""")
for lbl, Eobj in [("Haar physical seed 24680 (pure)", pauli_table(haar_physical(car, 24680), L, "direct")),
                  ("mixture of Haar seeds 111 and 222 (mixed)",
                   MixE([haar_physical(car, 111), haar_physical(car, 222)], [0.5, 0.5], L, "direct"))]:
    HW = algebra_entropy(Eobj, POINTER)[0]
    cum = [mutual_information(Eobj, POINTER,
                              sum([gens_FULL(car, f) for f in FR[:m]], []))[0] for m in range(1, 5)]
    sing = [mutual_information(Eobj, POINTER, gens_FULL(car, f))[0] for f in FR]
    print(f"     {lbl:<44} H = {HW:.6f}  I = {['%.6f' % c for c in cum]}  "
          f"R_delta = {rdelta_count(sing, HW, DELTA)} of 4")
print("     -- neither purity nor mixedness manufactures a plateau: both Haar arms are flat-free.")
print("     What the BRD_W arms exhibit at L = 10 is a record, not an artefact.")

# =====================================================================================
hr("PART 4  --  READING 2's THRESHOLD IS SIX LINKS.  THE DERIVATION OF 12 HAS A LOGICAL ERROR,")
hr("            AND SO DOES MY OWN 10 OF PART 3c.  I REPORT BOTH.", "-")
print("""
  b4 PART 3's derivation: 'by the PART A theorem the system region must contain a cycle, and EACH
  FRAGMENT MUST CONTAIN A CYCLE FOR ITS ALGEBRA TO HAVE ANY MAGNETIC CONTENT, so on a girth-2
  graph every region costs 2 links'.  The second clause does not follow from theorem B-8.
  B-8 says where a magnetic OPERATOR can be LOCALISED.  A record is not a localisation: a fragment
  needs no magnetic content of its own to be CORRELATED with a holonomy.  Z_0 Z_1 X_l is a
  perfectly good gauge-invariant operator, so a state can pin the ELECTRIC operator X_l of a single
  link to the HOLONOMY W_01 of a distant plaquette.  A purely electric single-link fragment then
  records a magnetic bit.  My PART 3c inherited half of the same error and said 10.

  CONSTRUCTION.  Take the stabiliser group generated by the independent Gauss operators together
  with  S_l = W_S X_l  for each fragment link l, where W_S is the system region's Wilson loop.
  Every generator is gauge invariant, they mutually commute, and W_S itself is NOT in the group, so
  the stabilised space splits into W_S = +1 and W_S = -1 halves.  Their equal superposition is a
  PHYSICAL state with H(W_S) = 1 bit in which every fragment's electric operator reads it.
  ONE VARIABLE MOVES between the two rows of each block: the state (broadcast versus Haar), with
  carrier, cut, pointer, fragment algebras (C3_FULL) and delta held fixed.""")

def apply_sp(psi, a, Lc):
    D = 1 << Lc; idx = np.arange(D); popc = popc_table(D)
    out = ((-1.0) ** popc[idx & a.z]) * psi
    out = out[idx ^ a.x]
    return a.c * (1j) ** (pc(a.x & a.z) % 4) * out

def holonomy_broadcast(car, S, frags):
    """Physical state with H(W_S) = 1 bit recorded in every single-link fragment.  None if the
       stabiliser count leaves no room (see the wheel_4 row)."""
    Lc = car["L"]; Szz = sum(1 << l for l in S)
    stab = [SP(g, 0) for g in indep_gauss(car)] + [SP(sum(1 << x for x in f), Szz) for f in frags]
    W = SP(0, Szz); found = {}
    for seed in range(1 << Lc):
        p = np.zeros(1 << Lc, dtype=complex); p[seed] = 1.0
        for g in stab: p = 0.5 * (p + apply_sp(p, g, Lc))
        n = np.linalg.norm(p)
        if n < 1e-9: continue
        p /= n
        w = round(float(np.vdot(p, apply_sp(p, W, Lc)).real), 6)
        found.setdefault(w, p)
        if 1.0 in found and -1.0 in found: break
    if 1.0 not in found or -1.0 not in found:
        return None, len(stab)
    v = found[1.0] + found[-1.0]
    return v / np.linalg.norm(v), len(stab)

CASES = [(theta(6), [0, 1], [[2], [3], [4], [5]], "girth 2"),
         (fan(3), [0, 1, 2], [[3], [4], [5], [6]], "girth 3, interior vertices"),
         (wheel(4), [0, 1, 4], [[2], [3], [5], [6]], "girth 3, interior vertices")]
for car, S, FRs, note in CASES:
    Lc = car["L"]; Szz = sum(1 << l for l in S)
    print(f"\n  {car['name']}  L = {Lc}  ({note}).  S = {S}, cyclomatic(S) = {cyclomatic(car,S)}, "
          f"pointer W_S = {sp_str(SP(0,Szz), Lc)}")
    print(f"     fragments = {FRs}, each a SINGLE link with a purely electric algebra alg{{X_l}}")
    psi, nstab = holonomy_broadcast(car, S, FRs)
    if psi is None:
        print(f"     NO SUCH STATE: {len(indep_gauss(car))} independent Gauss operators + "
              f"{len(FRs)} fragment stabilisers = {nstab} = L, so the stabilised space is")
        print(f"     one-dimensional and W_S is deterministic.  The construction needs")
        print(f"     L >= (independent Gauss operators) + (fragments) + 1.  RECORDED AS A NULL.")
        continue
    print(f"     ||psi - P_phys psi|| = {np.linalg.norm(psi - project_physical(psi, car)):.2e}"
          f"   (physical, no projection applied)")
    for lbl, st in [("BROADCAST (constructed)", psi),
                    ("HAAR physical control  ", haar_physical(car, 31415))]:
        Es = pauli_table(st, Lc, "direct"); P = [SP(0, Szz)]
        HW = algebra_entropy(Es, P)[0]
        sing = [mutual_information(Es, P, gens_FULL(car, f))[0] for f in FRs]
        cum = [mutual_information(Es, P, sum([gens_FULL(car, f) for f in FRs[:m]], []))[0]
               for m in range(1, len(FRs) + 1)]
        uni = [mutual_information(Es, P, gens_FULL(car, sum(FRs[:m], [])))[0]
               for m in range(1, len(FRs) + 1)]
        print(f"     {lbl}  H(W_S) = {HW:.6f}")
        print(f"        per-fragment I(W_S:F_k) = {['%.6f' % v for v in sing]}"
              f"   R_delta = {rdelta_count(sing, HW, DELTA)} of {len(FRs)}")
        print(f"        cumulative, JOIN reading  = {['%.6f' % v for v in cum]}")
        print(f"        cumulative, UNION reading = {['%.6f' % v for v in uni]}")
print("""
  RESULT.  FOUR plateau points at I = H(W_S) = 1.000000 bit, R_delta = 4 of 4, plateau defect 0,
  in the FULL GAUGE-INVARIANT ALGEBRA C3_FULL with a declared Wilson-loop pointer, on a PHYSICAL
  state, at SIX LINKS on a girth-2 carrier and SEVEN on a girth-3 carrier.  The JOIN and UNION
  readings agree here, so B-11's ambiguity does not rescue the larger numbers.
  The Haar control on the identical carrier, cut, pointer and algebra returns R_delta = 0 and no
  flat curve, so the test had a way to fail.

  THE THRESHOLD TABLE, CORRECTED:
     the target says   reading 1 = 6 links,   reading 2 = 12 links (18 on girth 3)
     PART 3c said      reading 1 = 6 links,   reading 2 = 10 links   <-- MY OWN INTERMEDIATE ERROR
     this part says    reading 1 = 6 links,   reading 2 =  6 links (7 on girth 3)
  AND 6 IS TIGHT under BOTH readings for the same reason: four plateau points need four disjoint
  fragments, each at least one link, plus a system region that must contain a cycle -- at least
  girth links -- so L >= girth + 4.  theta_6 attains it.
  WHAT THIS COSTS THE TARGET'S HEADLINE: 'the two readings differ by a factor of two in link count'
  is REFUTED.  They agree, at six.  The thing that genuinely differs between the readings is the
  VERDICT ON A GIVEN STATE (B-7), not the carrier size -- and B-7 survives this file untouched.""")
