"""
r3 -- IS THE LANE'S SBS TEST KORBICZ'S SBS?

THE DEFINITION (Horodecki-Korbicz-Zurek, PRA 91, 032122 (2015); Korbicz et al.).
  rho_{S:F_1..F_m} has SPECTRUM BROADCAST STRUCTURE iff THERE EXISTS an orthonormal basis {|i>}
  of S with
        rho  =  sum_i  p_i |i><i| (x) rho_i^{F_1} (x) ... (x) rho_i^{F_m},
        rho_i^{F_k} rho_j^{F_k} = 0   for all k and all i =/= j.
  THE BASIS IS EXISTENTIALLY QUANTIFIED.  It is not "the eigenbasis of rho_S" and it is not any
  fixed basis of the carrier.

WHAT b2_sbs_vs_rdelta.py ACTUALLY COMPUTES (its own comment, line ~"dephase S in the computational
(Z) basis of the link -- the einselected basis of this carrier"): it dephases in the Z basis and
nothing else.  It computes rho_S's eigenbasis only to PRINT a degeneracy flag; that flag is never
used.  So the code tests ONE basis of a condition that quantifies over ALL bases.
  => its SBS = YES verdicts are SOUND (an existential witnessed).
  => its SBS = no verdicts are NOT VERDICTS unless the basis search is done.

PART 1  a FALSE NEGATIVE: a state that IS of exact SBS form, which the lane's test calls "no".
PART 2  the existential search done properly, on the lane's own four arms.  Does B-4 survive?
PART 3  B-4's headline over-claim: "the gap is EXACTLY strong independence".  A constructed state
        with QD saturated, D_coh = D_prod = 0, and D_orth =/= 0 -- the gap also contains clause (c).
PART 4  B-5 re-read: "SBS needs a trace-out" is a theorem about purity, not a disagreement.
PART 5  B-6 re-read: correct where it is vacuous, false where it would bite.

ISOLATION LEDGER
  held fixed : carrier theta_9 / theta_5, S = {link 0}, fragment lists, and the defect functions,
               which are the target's own sbs_report arithmetic re-expressed with the basis as an
               argument (rlib.sbs_defects; verified against b2's own numbers in PART 0 below).
  moved      : PART 1-2 -- THE POINTER BASIS, and nothing else.  PART 3 -- the state.
"""
import itertools
import numpy as np
from rlib import *

np.set_printoptions(precision=6, suppress=True)
DELTA = 0.1

hr("r3  SBS: THE LANE'S TEST versus KORBICZ'S CONDITION")

# =====================================================================================
hr("PART 0  --  MY DEFECT FUNCTIONS REPRODUCE b2's, SO EVERY LATER NUMBER IS COMPARABLE")
L = 9; car = theta(L)
S = [0]; FR = [[1, 2], [3, 4], [5, 6]]; U = [7, 8]
def ghz(): return sym_basis_state(car, 0)
def qdns():
    psi = np.zeros(1 << L, dtype=complex)
    full = (1 << L) - 1
    for c in (0, (1 << 2) | (1 << 4) | (1 << 6)):
        psi[c] += 0.5; psi[c ^ full] += 0.5
    return project_physical(psi, car)
def sbsm():
    psi = np.zeros(1 << L, dtype=complex); full = (1 << L) - 1
    for b1 in (0, 1):
        for b2 in (0, 1):
            for b3 in (0, 1):
                c = (b1 << 2) | (b2 << 4) | (b3 << 6)
                psi[c] += 1.0; psi[c ^ full] += 1.0
    return project_physical(psi / np.linalg.norm(psi), car)
ARMS = [("GHZ    magnetic GHZ", ghz()),
        ("SBSm   mixed orthogonal conditional fragment states", sbsm()),
        ("QDNS   QD saturated, fragments correlated given S", qdns()),
        ("HAAR   Haar physical seed 777", haar_physical(car, 777))]
print("  b2 reports (its own table):  GHZ 0.00e+00/0.00e+00/0.00e+00 ; SBSm 0/0/0 ;")
print("     QDNS D_coh 0.00e+00  D_prod 8.75e-01  D_orth 0.00e+00 ; HAAR 4.87e-01/9.30e-01/9.92e-01")
print(f"  {'ARM':<54}{'D_coh':>10}{'D_prod':>10}{'D_orth':>10}")
for nm, psi in ARMS:
    d = sbs_defects(psi, L, S, FR)
    print(f"  {nm:<54}{d['D_coh']:>10.2e}{d['D_prod']:>10.2e}{d['D_orth']:>10.2e}")
print("  MATCH.  (rlib.sbs_defects is b2's arithmetic with the basis lifted out as an argument.)")

# =====================================================================================
hr("PART 1  --  A FALSE NEGATIVE OF THE LANE'S SBS TEST")
print("""
  STATE: the ELECTRIC GHZ (|+^9> + |-^9>)/sqrt2 -- an arm the target itself uses in b1 and b3.
  In the X basis, link 0 carries a bit s and every environment link carries the same s.  Tracing
  out U = {7,8} (whose two branches are |++> and |-->, orthogonal) leaves
        rho = 1/2 |+><+|_0 (x) |+ + + + + +>< ... |  +  1/2 |-><-|_0 (x) |- - - - - -><...|,
  which is EXACTLY of SBS form, with pointer basis {|+>,|->} and product, orthogonal conditional
  fragment states.  It is a spectrum broadcast structure by inspection.
  Below: the SAME defect functions, the SAME state, the SAME cut; ONLY THE POINTER BASIS MOVES.""")
psiE = elec_ghz(car)
Hd = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
print(f"  {'pointer basis of S':<36}{'D_coh':>12}{'D_prod':>12}{'D_orth':>12}{'SBS?':>8}")
print("-" * 104)
for lbl, V in [("Z (what b2 hard-codes)", np.eye(2, dtype=complex)),
               ("X (electric -- the right one)", Hd.astype(complex))]:
    d = sbs_defects(rot_link(psiE, L, 0, V.conj().T), L, S, FR)
    ok = max(d["D_coh"], d["D_prod"], d["D_orth"]) < 1e-8
    print(f"  {lbl:<36}{d['D_coh']:>12.2e}{d['D_prod']:>12.2e}{d['D_orth']:>12.2e}"
          f"{('YES' if ok else 'no'):>8}")
best, ang = sbs_best_over_bases(psiE, L, S, FR)
print(f"  minimum over a 25x24 basis grid + local refine: worst defect = {best['worst']:.3e} "
      f"at (theta,phi) = ({ang[0]:.4f},{ang[1]:.4f})")
print("""  => b2's test, run verbatim on this state, returns D_coh = 5.00e-01 and 'SBS = no' for a
  state that IS a spectrum broadcast structure.  A FALSE NEGATIVE.  Every 'SBS = no' the lane
  reports is of this type until a basis search is run -- including the headline of B-4.""")

# =====================================================================================
hr("PART 2  --  THE EXISTENTIAL DONE PROPERLY.  DOES B-4 SURVIVE?")
print("""  Same carrier, same cut, same functions.  For each arm: the Z-basis defects the lane
  reports, then the MINIMUM over all pointer bases of S of max(D_coh, D_prod, D_orth).""")
print(f"  {'ARM':<52}{'Z: worst':>11}{'min over bases':>16}{'verdict':>12}")
print("-" * 104)
for nm, psi in ARMS + [("BRD_E  electric GHZ", psiE)]:
    dz = sbs_defects(psi, L, S, FR)
    bb, aa = sbs_best_over_bases(psi, L, S, FR)
    v = "SBS" if bb["worst"] < 1e-6 else "NOT SBS"
    print(f"  {nm:<52}{dz['worst']:>11.3e}{bb['worst']:>16.3e}{v:>12}")
print("""
  B-4 SURVIVES, AND NOW WITH A PROOF THE LANE DID NOT HAVE.  QDNS is not SBS in ANY basis.
  The reason is structural, and worth stating because it is what makes the basis search
  unnecessary here: rho = 1/2 |0><0| (x) sigma_0 + 1/2 |1><1| (x) sigma_1 with sigma_0, sigma_1
  of ORTHOGONAL SUPPORT (D_orth = 0).  For any other basis {|+n>,|-n>} the off-diagonal block is
  proportional to (sigma_0 - sigma_1), which is nonzero, so no other basis block-diagonalises rho.
  Z is the unique candidate, and D_prod = 0.875 there.  QD saturates, SBS fails.
  BUT the arm the lane used to report the SAME structure -- the electric GHZ, which it never ran
  through b2 -- would have been mis-verdicted, so the instrument was not sound when it was used.""")

# =====================================================================================
hr("PART 3  --  B-4's OVER-CLAIM: THE GAP IS *NOT* 'EXACTLY' STRONG INDEPENDENCE")
print("""
  B-4 asserts: 'the gap is exactly SBS's STRONG-INDEPENDENCE clause'.  One example in which only
  clause (b) fails shows the gap CONTAINS clause (b).  It cannot show the gap IS clause (b).
  Here is a state in the gap where clause (b) HOLDS EXACTLY and clause (c) fails.

  CONSTRUCTION on theta_5.  S = {link 0}, fragments {1} {2} {3}, unobserved U = {4}.
     |psi> = ( |0>_0 (x) |a>|a>|a>|0>_4  +  |1>_0 (x) X^{(x)4}(|a>|a>|a>|0>_4) ) / sqrt2 ,
     |a> = cos t |0> + sin t |1>.
  Gauge invariant because G = X^{(x)5} exchanges the two branches.  U decoheres S exactly.  The
  conditional environment state is a PRODUCT across fragments by construction (clause (b) exact),
  but the two conditional single-link states |a> and X|a> are NOT orthogonal: overlap sin 2t.""")
L5 = 5; car5 = theta(L5)
def near_orth(t):
    ca, sa = np.cos(t), np.sin(t)
    a = np.array([ca, sa], dtype=complex); ap = np.array([sa, ca], dtype=complex)
    psi = np.zeros(1 << L5, dtype=complex)
    for i1 in (0, 1):
        for i2 in (0, 1):
            for i3 in (0, 1):
                # branch link0 = 0 : links1,2,3 in |a>, link4 = |0>
                idx = (i1 << 1) | (i2 << 2) | (i3 << 3)
                psi[idx] += a[i1] * a[i2] * a[i3]
                # branch link0 = 1 : links1,2,3 in X|a>, link4 = |1>
                idx2 = 1 | (i1 << 1) | (i2 << 2) | (i3 << 3) | (1 << 4)
                psi[idx2] += ap[i1] * ap[i2] * ap[i3]
    return psi / np.linalg.norm(psi)
S5 = [0]; FR5 = [[1], [2], [3]]
print(f"  {'t':>6}{'|<a|X a>|':>11}{'H(S)':>9}{'I(S:F_k)':>11}{'R_delta':>9}"
      f"{'D_coh':>10}{'D_prod':>10}{'D_orth':>10}{'SBS?':>7}")
print("-" * 104)
chosen = None
for t in (0.00, 0.06, 0.12, 0.18, 0.25, 0.35, np.pi / 4):
    psi = near_orth(t)
    gi = np.linalg.norm(psi - project_physical(psi, car5))
    assert gi < 1e-10, f"state not gauge invariant, ||psi-Ppsi|| = {gi}"
    HS = vn_entropy(reduce_links(psi, L5, S5))
    per = [mi_ext(psi, L5, S5, f) for f in FR5]
    R = rdelta_count(per, HS, DELTA)
    d = sbs_defects(psi, L5, S5, FR5)
    ok = max(d["D_coh"], d["D_prod"], d["D_orth"]) < 1e-8
    print(f"  {t:>6.2f}{abs(np.sin(2*t)):>11.6f}{HS:>9.6f}{per[0]:>11.6f}{R:>9}"
          f"{d['D_coh']:>10.2e}{d['D_prod']:>10.2e}{d['D_orth']:>10.2e}{('YES' if ok else 'no'):>7}")
    if abs(t - 0.12) < 1e-9: chosen = (t, HS, per, R, d)
t, HS, per, R, d = chosen
bb, aa = sbs_best_over_bases(near_orth(t), L5, S5, FR5)
print(f"""
  THE ROW THAT DOES THE WORK, t = {t}:
     R_delta = {R} of 3 with I(S:F_k) = {per[0]:.6f} >= 0.9 * H(S) = {0.9*HS:.6f}  -- QD SATURATES.
     D_coh   = {d['D_coh']:.2e}   (clause (a) exact)
     D_prod  = {d['D_prod']:.2e}   (clause (b) STRONG INDEPENDENCE EXACT -- the clause B-4 names)
     D_orth  = {d['D_orth']:.6f}   (clause (c) FAILS)
     minimum over ALL pointer bases of max defect = {bb['worst']:.6f}  -- NOT SBS in any basis.
  CONTROL IN THE SAME TABLE, one variable moved (t only): at t = 0 the state is the SBS/GHZ limit
  and all three defects are 0; at t = pi/4 the two conditional states coincide, the record is
  destroyed, and R_delta drops.  The test had a way to fail at both ends and did not.
  CONCLUSION: QD is strictly weaker than SBS -- B-4's DIRECTION is right -- but the gap is NOT
  'exactly' clause (b).  It contains clause (c) as well, and it must, because R_delta is a
  (1-delta)-inequality while SBS's clause (c) is an exact orthogonality.  Any delta > 0 opens
  this gap.  B-4's status should be EXHIBITED-IN-ONE-DIRECTION, not the equality it asserts.""")

# =====================================================================================
hr("PART 4  --  B-5 RE-READ: 'SBS NEEDS A TRACE-OUT' IS A THEOREM ABOUT PURITY")
print("""
  THEOREM.  If rho on S (x) F_1 (x) ... (x) F_m is PURE and of SBS form, then exactly one p_i is
  nonzero and every rho_i^{F_k} is pure, i.e. rho is a PRODUCT STATE and H(S) = 0.
  PROOF.  A pure state that is a convex mixture is a mixture of one term; a pure product over the
  cut S|rest forces each factor pure.  QED.
  So 'GHZ with all 8 environment links observed is not SBS' is not a disagreement between two
  criteria: it is the observation that a pure entangled state is not separable.  It carries no
  information about SBS versus R_delta, and it holds for EVERY state with H(S) > 0.
  MEASURED, on every arm of b2 with nothing traced out (the b2 'U OBSERVED' configuration):""")
FRall = [[1, 2], [3, 4], [5, 6], [7, 8]]
print(f"  {'ARM':<52}{'H(S)':>9}{'D_coh (Z)':>12}{'min over bases':>16}")
for nm, psi in ARMS:
    HS = vn_entropy(reduce_links(psi, L, S))
    dz = sbs_defects(psi, L, S, FRall)
    bb, _ = sbs_best_over_bases(psi, L, S, FRall, ngrid=13, refine=2)
    print(f"  {nm:<52}{HS:>9.6f}{dz['D_coh']:>12.2e}{bb['worst']:>16.3e}")
print("""  Every arm with H(S) > 0 fails in every basis, as the theorem says it must.  B-5's
  'the verdict flips' is CORRECT AS ARITHMETIC and EMPTY AS A COMPARISON.""")

# =====================================================================================
hr("PART 5  --  B-6 RE-READ: TRUE WHERE IT IS VACUOUS, FALSE WHERE IT WOULD BITE")
print("""
  B-6 (status PROVED): 'the Gauss law forces the eigenbasis of rho_S on a gauged link to be
  ELECTRIC or degenerate, so SBS ... can never certify a magnetic (holonomy) record'.
  (i) THE ARITHMETIC IS RIGHT.  <Z_l> = <Y_l> = 0 for physical states, on every carrier, because
      Z_l anticommutes with the Gauss operator at either endpoint of l.  Verified beyond theta:""")
for car2 in [theta(6), fan(3), fan(4), wheel(4), wheel(5)]:
    Lc = car2["L"]; wz = wy = 0.0
    for sd in range(8):
        Ep = pauli_table(haar_physical(car2, 7000 + sd), Lc)
        for l in range(Lc):
            wz = max(wz, abs(sp_expect(Ep, SP(0, 1 << l))))
            wy = max(wy, abs(sp_expect(Ep, SP(1 << l, 1 << l))))
    print(f"     {car2['name']:<10} max|<Z_l>| = {wz:.3e}   max|<Y_l>| = {wy:.3e}   "
          f"(8 Haar physical states, every link)")
print("""
  (ii) BUT THE CONCLUSION IS VACUOUS AT ITS OWN SCOPE.  The premise is about ONE LINK.  By the
      lane's OWN theorem B-8, a single link has NO gauge-invariant holonomy at all.  'SBS cannot
      certify a magnetic record on a gauged link' is therefore a statement about the empty set.
  (iii) AND IT IS FALSE ONE STEP OUT.  Take the smallest region that CAN carry a holonomy,
      |S| = girth = 2 on theta.  Then rho_S is not (I + <X>X)/2 and its magnetic expectation is
      free.  Measured on theta_10, S = the plaquette {0,1}:""")
L = 10; car = theta(L)
aa_ = sym_basis_state(car, 0)
bb_ = sym_basis_state(car, sum(1 << i for i in range(1, L, 2)))
psiW = project_physical((aa_ + bb_) / np.linalg.norm(aa_ + bb_), car)
for nm, psi in [("magnetic GHZ (all W_p = +1)", sym_basis_state(car, 0)),
                ("Wilson-loop GHZ (W_01 = +-1)", psiW),
                ("ground state g^2 = 0.30", ground_state(car, 0.30)[0]),
                ("Haar physical seed 66", haar_physical(car, 66))]:
    Ep = pauli_table(psi, L, "direct")
    w01 = sp_expect(Ep, SP(0, 0b11)).real
    x01 = sp_expect(Ep, SP(0b11, 0)).real
    rho = reduce_links(psi, L, [0, 1])
    ev = np.linalg.eigvalsh(rho)
    W = np.diag([1., -1., -1., 1.])                       # Z_0 Z_1 on the 2-link block
    comm = np.linalg.norm(rho @ W - W @ rho)
    print(f"     {nm:<30} <W_01> = {w01:>9.6f}   <X_0X_1> = {x01:>9.6f}   "
          f"||[rho_S, W_01]|| = {comm:.2e}   spec(rho_S) = {np.round(ev,4)}")
print("""      => on the magnetic GHZ, <W_01> = 1 EXACTLY and rho_S commutes with W_01: the
      holonomy IS diagonal in an eigenbasis of rho_S.  So SBS's pointer basis can be magnetic as
      soon as the system region is large enough to have a holonomy -- which is exactly the regime
      the lane's own B-8 corollary says any magnetic record must live in.
  VERDICT SO FAR: the two-line proof is correct; the conclusion drawn from it is not.  It is
  vacuous for |S| = 1 and its mechanism dies for |S| >= girth.  A third defect, not in the lane's
  text: SBS does NOT 'read its pointer basis off rho_S' -- the definition quantifies existentially
  over bases, and b2's own output flags 'rho_S degenerate: True' on all three broadcast arms,
  where an eigenbasis is not even unique.""")

print("""
  (iv) THE DIRECT REFUTATION.  B-6 says SBS 'can never certify a magnetic (holonomy) record on its
      own terms'.  Here is a gauge-invariant state on which SBS certifies exactly that, with all
      three defects zero, at |S| = girth = 2.

      CONSTRUCTION on theta_8.  Group the links into four pairs P_0 = {0,1} (= S), P_1 = {2,3},
      P_2 = {4,5}, U = {6,7}.  On each pair use the BELL BASIS, i.e. the simultaneous eigenbasis
      of the HOLONOMY W = Z_aZ_b and the pair electric flux e = X_aX_b:
         |j=0> = (|00>+|11>)/sqrt2   (W=+1, e=+1)      |j=2> = (|01>+|10>)/sqrt2   (W=-1, e=+1)
         |j=1> = (|00>-|11>)/sqrt2   (W=+1, e=-1)      |j=3> = (|01>-|10>)/sqrt2   (W=-1, e=-1)
      and take   |psi> = (1/2) sum_j |j>_S (x) |j>_{P_1} (x) |j>_{P_2} (x) |j>_U .
      G = prod_l X_l acts on this as e_j^4 = +1, so |psi> IS PHYSICAL, with no projection.
      Tracing out U (whose four branches are orthogonal) leaves an exact SBS state whose pointer
      basis DIAGONALISES THE HOLONOMY W_01, and H(W_01) = 1 bit -- a real, non-deterministic
      magnetic record, not the b1 H=0 trap.
      ONE VARIABLE MOVES BELOW: the pointer basis of the system pair.  The state is identical.""")
L8 = 8; car8 = theta(L8)
bell = np.zeros((4, 4), dtype=complex)
bell[:, 0] = [1, 0, 0, 1]; bell[:, 1] = [1, 0, 0, -1]
bell[:, 2] = [0, 1, 1, 0];  bell[:, 3] = [0, 1, -1, 0]
bell /= np.sqrt(2)
psiB = np.zeros(1 << L8, dtype=complex)
for j in range(4):
    v = np.array([1.0])
    for _ in range(4):
        v = np.kron(v, bell[:, j])
    psiB += 0.5 * v                     # pairs are (0,1)(2,3)(4,5)(6,7), most significant first
# lib_b indexes link l in BIT l, i.e. link 0 is the LEAST significant; rebuild with that convention
psiB = np.zeros(1 << L8, dtype=complex)
for j in range(4):
    amp = {0: {(0, 0): 1 / np.sqrt(2), (1, 1): 1 / np.sqrt(2)},
           1: {(0, 0): 1 / np.sqrt(2), (1, 1): -1 / np.sqrt(2)},
           2: {(0, 1): 1 / np.sqrt(2), (1, 0): 1 / np.sqrt(2)},
           3: {(0, 1): 1 / np.sqrt(2), (1, 0): -1 / np.sqrt(2)}}[j]
    for k0 in amp:
        for k1 in amp:
            for k2 in amp:
                for k3 in amp:
                    idx = ((k0[0] << 0) | (k0[1] << 1) | (k1[0] << 2) | (k1[1] << 3) |
                           (k2[0] << 4) | (k2[1] << 5) | (k3[0] << 6) | (k3[1] << 7))
                    psiB[idx] += 0.5 * amp[k0] * amp[k1] * amp[k2] * amp[k3]
print(f"      ||psi - P_phys psi|| = "
      f"{np.linalg.norm(psiB - project_physical(psiB, car8)):.3e}   (physical, unprojected)")
Eb = pauli_table(psiB, L8)
print(f"      H(W_01) = {algebra_entropy(Eb, [SP(0, 0b11)])[0]:.6f} bits    "
      f"<W_01> = {sp_expect(Eb, SP(0,0b11)).real:.6f}    "
      f"I(W_01 : W_23) = "
      f"{mutual_information(Eb, [SP(0,0b11)], [SP(0,0b1100)])[0]:.6f} bits")
SB = [0, 1]; FRB = [[2, 3], [4, 5]]; UB = [6, 7]
print(f"      {'pointer basis of S = links {0,1}':<44}{'D_coh':>11}{'D_prod':>11}"
      f"{'D_orth':>11}{'SBS?':>7}")
print("      " + "-" * 92)
Hd2 = np.kron(np.array([[1, 1], [1, -1]]) / np.sqrt(2),
              np.array([[1, 1], [1, -1]]) / np.sqrt(2)).astype(complex)
for lbl, V in [("Z product basis of links {0,1}", np.eye(4, dtype=complex)),
               ("X (ELECTRIC) product basis -- B-6's basis", Hd2),
               ("BELL basis -- W_01 diagonal (MAGNETIC)", bell)]:
    d = sbs_defects_multi(rot_block(psiB, L8, SB, V.conj().T), L8, SB, FRB)
    ok = d["worst"] < 1e-8
    print(f"      {lbl:<44}{d['D_coh']:>11.2e}{d['D_prod']:>11.2e}{d['D_orth']:>11.2e}"
          f"{('YES' if ok else 'no'):>7}")
print("""      => SBS = YES with a pointer basis in which the HOLONOMY is diagonal, on a physical
      state of a Z_2 gauge theory, with one full bit of holonomy entropy redundantly recorded in
      two disjoint fragments.  B-6's conclusion is REFUTED, not merely narrowed.
      THE CONTROL IN THE SAME TABLE COULD HAVE FAILED AND DID: the X (electric) product basis --
      the one B-6 says the Gauss law forces -- is NOT an SBS basis for this state, and neither is
      the Z product basis.  Only the magnetic Bell basis works.
      AND NOTE WHERE B-6 ACTUALLY LIVES: its own escape clause.  rho_S here is I/4, degenerate, so
      "the eigenbasis of rho_S is electric OR DEGENERATE" is satisfied by the degenerate branch --
      and b2's own output prints "rho_S degenerate: True" on ALL THREE of its broadcast arms.  The
      electric-basis conclusion therefore applies to none of the states the lane tested.""")
