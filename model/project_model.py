"""THE UNIVERSAL RECORD MODEL (URM) — named by the principal, 2026-08-20.
One object that represents the entire program, and the overall representation of the proof:
PROOF_V002 is this model's narration. "Universal" is operational — one model, any record surface
(T-33: six mechanisms at machine precision; D-25: world-tier surfaces provenance-pinned) — and is a
title, not a claim.

The principal, 2026-08-20: "We need a model to work on that represents the full project."

Layers, each resting on registered rows:
  DEFINITION   clauses (i')-(v') with tolerances (O-51; census GR1-GR4)  -> ProjectModel.clauses
  LAWS         lifetime and steady value from the record's own Liouvillian mode
               (C-69, C-70 -- FORMAL, six mechanisms)                    -> .lifetime, .steady_value
  FORMATION    written vs unwritten configuration; the five source
               standards (C-71 PROVED; C-72 PARTIAL)                     -> .configuration
  CORNER       t_m -> infinity, W -> 0 reduces to the exact commutant
               records of record_model (C-75, T-28)                      -> .corner
  GEOMETRY     the C-77 increments as a first-class layer (T-46):
               earned distance d_W (C-78), interface rank + the
               boundary law, both tiers (C-79), the two-region
               taxonomy and w_min = d (C-80), certifiability and
               the locality zero (C-81), the certification window
               law (C-83), encoding-level formation (C-71/C-72)          -> .distance,
               .interface_rank, .world_interface, .mutual_interface,
               .winding_interface, .coupling_cost, .certifiability,
               .cert_window, .formation_orientation, .formation_occupancy
               (machinery ported from the sealed lanes in model/geometry.py)
  ARROW        record threshold, entanglement-without-record, inside
               invariance, history, and fragment redundancy (F-17..F-21) -> .arrow_*
  COUNT LAW    the surviving-record census and its two widths (C-86)      -> .census, .count_widths
  CLASSES      venue-earned coupling classes and critical kernel
               (C-87, C-90)                                              -> .coupling_venue, .reachable_class, ...
  WRITING      conserving/critical writer kernels and the surface gap
               (C-91)                                                    -> .writing_*
  ROLES        the three role sentences, honestly scoped                 -> .roles

A RecordSurface is what a physicist supplies: their own device or specimen's constants, in SI.
New observations enter through the URM gate: a new surface carries provenance; a new law is a
layer method with a validator gate and claim row; a new external number is a pinned comparison
with stated tolerance, semantics, and a power control.  validate_urm.py is the integrated
conjunction; validate_project.py preserves the public-entry and backward-compatibility gates."""
import numpy as np
import grounded as G

_SZ = np.array([[1, 0], [0, -1]], dtype=complex)
_SP = np.array([[0, 1], [0, 0]], dtype=complex)
_SM = _SP.conj().T


class RecordSurface:
    """A real record surface: named constants in SI units. dE and E_b in J, T in K, f0 in Hz."""
    def __init__(self, name, mechanism, dE, E_b, T, f0, thermal=True):
        self.name, self.mechanism = name, mechanism
        self.dE, self.E_b, self.T, self.f0, self.thermal = dE, E_b, T, f0, thermal

    def open_system(self):
        """(H, {L_k}) for the two-state record. CONVENTION (corrected per the solidity review):
           E_b is the ACTIVATION ENERGY FROM THE METASTABLE (upper) WELL -- the standard Arrhenius
           meaning. Escape from upper: f0 exp(-E_b/kT); from lower: f0 exp(-(E_b+dE)/kT); detailed
           balance holds automatically. The earlier midpoint convention silently reduced the escape
           barrier by dE/2, producing a 5-6 order contradiction with the measured azobenzene cis
           lifetime. Returns None where the surface is NOT thermally activated -- the model DECLINES
           rather than returning a number."""
        if not self.thermal:
            return None
        kT = G.KB * self.T
        gu = self.f0 * np.exp(-self.E_b / kT)              # escape from the metastable well
        gd = self.f0 * np.exp(-(self.E_b + self.dE) / kT)  # reverse, from the lower well
        if not (np.isfinite(gd) and np.isfinite(gu)) or (gd < 1e-300 and gu < 1e-300):
            return None
        H = -(self.dE / 2) * _SZ
        return H, [np.sqrt(gd) * _SM, np.sqrt(gu) * _SP], _SZ


class ProjectModel:
    # ---------------------------------------------------------------- DEFINITION (O-51)
    def clauses(self, s, t_m, C_v, W_write):
        """The five amended clauses on a RecordSurface, with the tolerances declared."""
        sys = s.open_system()
        if sys is None:
            return dict(applies=False, why="not a thermally activated two-state record; the model declines")
        H, Ls, R = sys
        return dict(applies=True,
                    i=dict(spectral_family=True),
                    ii=G.clause_ii(H, Ls, R, t_m),
                    iii=G.clause_iii(s.dE, s.T, C_v),
                    iv=G.clause_iv(W_write, s.dE, s.T),
                    v=G.clause_v(s.E_b, s.T))

    # ---------------------------------------------------------------- LAWS (C-69, C-70; FORMAL)
    def lifetime(self, s):
        """tau from the RECORD'S OWN Liouvillian mode. Closed form exp(E_b/kT)/f0 is the check,
           never the source."""
        sys = s.open_system()
        if sys is None:
            return None
        H, Ls, R = sys
        return G.lifetime(H, Ls, R)

    def steady_value(self, s):
        """<R>_ss from the Liouvillian's stationary state. Closed form tanh(dE/2kT) is the check."""
        sys = s.open_system()
        if sys is None:
            return None
        H, Ls, R = sys
        if not Ls:
            return None    # no dissipator: a steady state is not defined; decline, never an arbitrary vector
        Lv = G.liouvillian(H, Ls)
        w, V = np.linalg.eig(Lv)
        j = int(np.argmin(np.abs(w)))
        rho = V[:, j].reshape(2, 2, order='F'); rho = rho / np.trace(rho)
        return float(np.real(np.trace(rho @ R)))

    # ---------------------------------------------------------------- FORMATION (C-71/C-72 machinery)
    def configuration(self, per_record, pattern):
        """The accumulated configuration quantity of a medium: per_record = the one-record magnitude
           (charge, moment-energy, ...), pattern = array of record values (0/±1). Returns the sum,
           the absolute sum, and |sum|/sum|.| -- the accumulation-vs-screening discriminator (C-46).
           The claim this machinery carries: the UNWRITTEN medium's accumulated quantity is null and
           WRITING creates an extensive, sign-definite value. Scoped pending the solidity review."""
        p = np.asarray(pattern, dtype=float)
        ssum, sabs = per_record * p.sum(), per_record * np.abs(p).sum()
        return dict(sum=ssum, abs_sum=sabs,
                    ratio=(abs(ssum) / sabs if sabs > 0 else None))

    # ---------------------------------------------------------------- CORNER (C-75; T-28)
    def corner(self, H, Ls=()):
        """DEF-A: at t_m -> infinity with no dissipation the slow modes are exactly the commutant,
           whose dimension is sum over eigenspaces of m_E^2. Both are returned so the reduction is
           CHECKED, never assumed."""
        rates, obs = G.slow_modes(H, list(Ls), np.inf)
        w = np.linalg.eigvalsh(np.asarray(H, dtype=complex)); mult = {}
        for x in w:
            key = round(float(x), 9); mult[key] = mult.get(key, 0) + 1
        return dict(slow_dim=len(rates), commutant_dim=sum(m * m for m in mult.values()))

    # ---------------------------------------------------------------- GEOMETRY (T-46; C-77 increments)
    # Every method below delegates to model/geometry.py, where the sealed-lane machinery
    # lives (exact F_2 bitmasks, writer searches, cut-ranks -- ported, not reimplemented).
    # Each geometry function's docstring names its claim row, sealed source, and owners.
    def distance(self, carrier, s, s2):
        """C-78: earned distance d_W between configurations s, s2 -- minimal admissible-
           writer weight, searched exhaustively (LANE_T42_A_DISTANCE, sealed).  carrier:
           the toric torus side L (the sealed venue).  On the torus d_W = L*a + L*b, the
           minimal boundary-crossing cost.  s, s2: label pairs in {0,1}^2."""
        import geometry as GE
        return GE.dW(int(carrier), s, s2)

    def distance_matrix(self, carrier):
        """C-78: the full 4x4 d_W class matrix on the toric carrier (side L = carrier),
           plus class minima and the computed code distance."""
        import geometry as GE
        return GE.dW_class_matrix(int(carrier))

    def interface_rank(self, L, s):
        """C-79 corner: cut-rank IR2 and content capacity capP of the s x s block at
           side L; thick scope: IR2 = 8s-10, capP = 2s^2+2s-5, and
           32*capP == (IR2+10)^2 + 8*(IR2+10) - 160 (LANE_T42_C_BOUNDCAP /
           LANE_T42_D_DERIVE, sealed)."""
        import geometry as GE
        return GE.interface_rank(L, s)

    def world_interface(self, n):
        """C-79 world: counted access aggregates of the n^3 barrier-record block --
           IFACE = 6n^2, CAP_d1 = n^3-(n-2)^3 = 6n^2-12n+8 (the C-82 per-epoch law),
           CAP_total = n^3, with IFACE^3 == 216*CAP_total^2 (LANE_T42_C_BOUNDCAP world
           tier / LANE_T42_D_DERIVE, sealed)."""
        import geometry as GE
        return GE.world_counts(n)

    def mutual_interface(self, L, s, g):
        """C-80: mutual interface I(A:B) of two s x s blocks at separation g -- zero
           identically when separated, seam law I_IR = 2(s-2) at single contact
           (LANE_O54_A_CORNER, sealed)."""
        import geometry as GE
        return GE.two_region_blocks(L, s, g)

    def winding_interface(self, L, w, r0):
        """C-80: the non-contractible pair's exact constant w at every distance
           (LANE_O54_A_CORNER control (ii), sealed)."""
        import geometry as GE
        return GE.winding_pair(L, w, r0)

    def coupling_cost(self, Lx, Ly, u, v):
        """C-80: the coupling-writer law w_min = d_gen between two hole records,
           exhaustive coset scan (LANE_O54_C_ATTEMPT, sealed; defect construction
           Bravyi-Kitaev)."""
        import geometry as GE
        return GE.coupling_cost(Lx, Ly, u, v)

    def certifiability(self, L, s):
        """C-81: the certify/write ledger of the s x s block -- CERT computed two sealed
           routes and equal to the cut-rank law 8s-10; WRITE0 = 0, the locality theorem
           (LANE_T43_A_CORNER, sealed)."""
        import geometry as GE
        return GE.certifiability(L, s)

    def cert_window(self, n, W):
        """C-83: CERT_W(n) = min(n^3, W*(6n^2-12n+8)); W=1 is C-82's per-epoch law,
           W=inf is DEF-A's immortal-record volume corner (LANE_T45_CLOCK, sealed)."""
        import geometry as GE
        return GE.cert_window(n, W)

    def formation_orientation(self):
        """C-71/C-72, orientation encoding (T-32, sealed seed 7): the accumulation table
           for DC-saturated / random / DC-free / AC-erased patterns, through the same
           configuration instrument as the FORMATION layer.  Returns
           dict[(name, N)] -> configuration dict (per_record = the CoCrPt grain moment)."""
        import geometry as GE
        return {key: self.configuration(GE.M_GRAIN, pat)
                for key, pat in GE.orientation_patterns().items()}

    def formation_occupancy(self):
        """C-71/C-72, occupancy encoding (T-34, sealed seed 11): written pages are
           one-signed for EVERY pattern (electron-only injection, mechanism-constitutive);
           the unwritten page is null within the declared +-5e/cell tolerance.  Returns
           dict(written_ratios, all_programmed, unwritten, null_within_tolerance)."""
        import geometry as GE
        pats = GE.occupancy_patterns()
        q = GE.N_E * GE.E_CHARGE
        # electron injection: each programmed cell stores -q; the sign was open in design
        # space and is fixed once per device by the write (the lane's premise, kept)
        ratios = [self.configuration(q, -p)['ratio'] for p in pats['written']
                  if p.sum() > 0]
        allp = self.configuration(q, -pats['all_programmed'])
        unw = self.configuration(GE.E_CHARGE, pats['unwritten_e'])
        tol = pats['tol_e'] * GE.E_CHARGE
        return dict(written_ratios=ratios, all_programmed=allp, unwritten=unw,
                    tolerance=tol, null_within_tolerance=bool(abs(unw['sum']) <= tol))

    # ---------------------------------------------------------------- ARROW (T-54/T-55; F-17..F-21)
    # Every method below delegates to model/arrow.py -- machinery ported from the sealed
    # arrow lanes or delegated to record_model, per that module's header.
    def arrow_threshold(self, lam=0.8, weights=(1, 2), coupling=None):
        """F-17: THE ARROW CARRIES THE RECORD'S OWN THRESHOLD.  chi(O:B) swept over every
           observable of each listed weight in the mean-force state under the weight-d
           coupling -- sealed: 0.00000000 on all 24 weight-1 observables, 0.11448276 at
           weight 2 = d, closed form from Z_B(+-1) exact (LANE_F1_ARROW part 4, sealed).
           Owner: ORIGINAL.  Scope: SINGLE-CARRIER, toric-2x2 (sealed T-9 audit)."""
        import arrow as AW
        return AW.arrow_threshold(lam=lam, weights=weights, coupling=coupling)

    def arrow_ledger(self, lam=0.8):
        """F-18: THE RECORD'S ARROW IS NOT AMBIENT DECOHERENCE.  The four-row coupling
           ledger (weight, I(S:B), chi(record:B)) -- sealed: the weight-1 coupling entangles
           I(S:B) = 0.04549256 yet transfers ZERO record bits (LANE_F1_ARROW/f1b part (c),
           sealed).  Owner: ORIGINAL.  Scope: SINGLE-CARRIER."""
        import arrow as AW
        return AW.arrow_ledger(lam=lam)

    def arrow_invariance(self, n_unitaries=12, seed=5, lam=0.8):
        """F-19: IRREVERSIBILITY FROM INSIDE.  System-only unitaries cannot reduce I(S:B)
           (sealed 3.686e-14 over 12 unitaries; covariance check passes; chi about the FIXED
           label moves 1.145e-01) -- the copy is relocatable, never erasable, from inside
           (LANE_F1_ARROW/f1b parts (a,b), sealed).  Owner: BORROWED (textbook invariance);
           the relative-arrow reading for records ours (T-III.6)."""
        import arrow as AW
        return AW.arrow_invariance(n_unitaries=n_unitaries, seed=seed, lam=lam)

    def arrow_history(self, times, coupling=None, lam=0.8, env=None, keep_states=()):
        """The arrow as a HISTORY (F-20-adjacent): chi(record:B)(t) from a product state
           with <record> exactly conserved -- the record is READ, not written; negative
           times are the sealed reversal control (LANE_PF2_DYNAMICAL, sealed; one
           eigendecomposition serves every time).  RecordModel.formation is the single-time
           twin (validate_formation.py).  Owner: ASSEMBLED.  F-20's mechanism TWO-CARRIER
           (toric-2x2; bouquet); this venue toric."""
        import arrow as AW
        return AW.arrow_history(times, coupling=coupling, lam=lam, env=env,
                                keep_states=keep_states)

    def arrow_redundancy(self, coupling=None, lam=0.8, t=4.0, env=None):
        """F-21: REDUNDANCY CARRIES THE RECORD'S THRESHOLD.  Whole-bath and per-fragment
           chi through RecordModel.redundancy -- sealed: fragments 0.789366/0.048377/
           0.678602 under the weight-d coupling, EXACTLY ZERO under weight-1
           (LANE_PF2_DYNAMICAL parts 3-4, sealed).  Owner: ASSEMBLED (quantum-Darwinism-
           style apparatus, Zurek/Blume-Kohout; the threshold-in-fragments finding ours).
           Scope: fragment bits SINGLE-CARRIER toric-only -- the T-9 battery replicated
           only the whole-bath weight-1 null on [[8,1,2]]/[[4,2,2]]."""
        import arrow as AW
        return AW.arrow_redundancy(coupling=coupling, lam=lam, t=t, env=env)

    def arrow_observation(self, env, coupling, record=None, model=None, lam=0.8, t=4.0,
                          tier="world", provenance=None):
        """OBSERVATION ENTRY for the arrow family (T-54/T-55): score a NEW bath/fragment
           observation through the family's own instruments -- I(S:B), whole-bath chi,
           per-fragment chi, and the verdicts (holds_record_bits, entangled_without_record
           = the F-18 class, redundant_fragments).  A custom RecordModel must bring its own
           explicit record; the toric default is used only with the default model.  D-25 AT
           THE GATE: world-tier baths require provenance; corner baths must self-declare
           'DEF-A'.  Every outcome registers -- entangled-without-record is a RESULT, not
           a failure."""
        import arrow as AW
        return AW.score_bath_observation(env, coupling, record=record, model=model, lam=lam,
                                         t=t, tier=tier, provenance=provenance)

    # ---------------------------------------------------------------- COUNT LAW (T-54; C-86)
    # Delegates to model/countlaw.py -- machinery ported from LANE_T47_A_WIDTH and
    # LANE_T47_B_STAIRCASE (sealed), registered via LANE_T47_D_REGISTER. The corner law
    # k = min_E v2(m_E) remains homed at model/count_law.py (C-14), referenced not duplicated.
    def census(self, surfaces, t_m):
        """C-86: THE SURVIVING-RECORD COUNT LAW k(t_m) -- the URM's wholly-owned
           falsifiable census (LANE_T47_A_WIDTH sec. 6, LANE_T47_B_STAIRCASE secs.
           C/D/G, sealed; registered via LANE_T47_D_REGISTER).  surfaces: a LIST of
           RecordSurface objects, world-tier ones built through URM.surface() so D-25
           provenance rides on every row; t_m: the retention spec in s.  Counts a
           record only while BOTH its values are durable -- clause (ii') on the
           record's OWN Liouvillian mode, |lambda_record| <= 1/t_m -- so each record
           dies at its SHALLOWER value's escape, at the parameter-free drop time
           t*_i = f0^-1 exp((B_i-dE_i)/kT)/(1+e^{-dE_i/kT}).  Returns k by the
           instrument AND k_formula by the derived width delta_pop(t_m) =
           kT ln(expm1(B/kT - ln f0 t_m)) so the agreement is CHECKED, never assumed;
           plus the dated drop schedule, delta_coh = hbar/t_m, the departure term
           sum_over_dead tanh(dE_i/2kT) (remanence persists while records die), and
           the declined list (non-thermal surfaces are declared, never silently
           counted).  THE SIGNATURE IS THE C-76 GATE: (surfaces, t_m) and nothing
           else -- no width, tolerance, or clustering parameter exists on this path;
           checks_countlaw.py gates the unreachability.  Owners per C-86 (ownership
           PARTIAL): Neel/Street-Woolley/Sharrock own the activation window as a
           remanence-decay device; Charap-Lu-He and Weller-Moser the dE = 0 corner;
           Preisach-Neel rate-level asymmetry; ours the derived width, the two-face
           unification, the margin-free integer census, the departure term."""
        import countlaw as CLW
        return CLW.census(surfaces, t_m)

    def count_widths(self, s, t_m):
        """C-86: both faces of clause (ii') on one RecordSurface at one retention
           spec -- delta_pop from the record's own constants (B = E_b + dE, the exact
           convention map; None when t_m is beyond the symmetric bound exp(B/kT)/(2 f0),
           which IS the no-crossing condition), delta_coh = hbar/t_m, and the record's
           own instrument drop time beside its closed form (the check, never the
           source).  LANE_T47_A_WIDTH, sealed."""
        import countlaw as CLW
        B = s.E_b + s.dE
        return dict(delta_pop=CLW.delta_pop(B, s.T, s.f0, t_m),
                    delta_coh=CLW.delta_coh(t_m),
                    t_star=CLW.drop_time(s),
                    t_star_formula=CLW.drop_time_formula(s))

    # ---------------------------------------------------------------- CLASSES (T-54; C-87/C-90)
    # Every method below delegates to model/classes.py, where the sealed-lane machinery
    # lives (exact ints/Fractions on the measurement path; certified tails; computed
    # booleans -- ported, not reimplemented).  Each docstring names its claim row, sealed
    # source, and owners.
    def coupling_venue(self, name, adj, provenance=None, tier="world", sector=None):
        """C-87 observation entry: a NEW venue graph (a record surface's access geometry)
           enters as an adjacency structure through the D-25 provenance gate -- world tier
           REFUSED without a pinned source, corner tier must self-declare 'DEF-A'
           (LANE_T44_B_WORLD S0 / LANE_T44_A_CORNER S0, sealed).  adj: list over nodes of
           [(neighbor, multiplicity), ...]; the gate also refuses empty adjacency,
           non-integer/out-of-range neighbor indices, and multiplicities that are not
           strictly positive integers.  sector optionally declares the venue limit
           ('Z3'/'Z2'/'Z1') for the evidence instruments."""
        import classes as CC
        return CC.venue(name, adj, provenance=provenance, tier=tier, sector=sector)

    def world_coupling_venue(self, n):
        """C-90: the world venue (earned D=3) -- the n^3 census-grain torus with the PINNED
           provenance (GR1 grains, face adjacency, T42_C/T43_B lineage; one walk step = one
           grain-boundary crossing = one unit of writer weight, the Gamma price C-80/O-54).
           Returns (VenueGraph, cells, idx).  (LANE_T44_B_WORLD S0, sealed.)"""
        import classes as CC
        return CC.world_venue(n)

    def corner_coupling_venue(self, Lx, Ly):
        """C-87: the corner venue (earned D=2) -- the dual lattice of the model's OWN ported
           carrier: plaquettes of geometry.Torus adjacent iff their supports share a carrier
           edge, computed from the plaquette masks alone, multiplicity kept; DEF-A
           self-declared.  (LANE_T44_A_CORNER/t44a_lib.py plaquette_adjacency,
           verbatim-in-substance; sealed row sums exactly 4 on (4,6) and (3,7).)"""
        import classes as CC
        return CC.corner_venue(Lx, Ly)

    def chain_coupling_venue(self, L):
        """C-87: the D=1 control venue (cycle C_L, DEF-A).  (t44a_lib.py cycle_adjacency,
           verbatim; sealed: row sums exactly 2, mu_c = 1/2 the venue's own.)"""
        import classes as CC
        return CC.chain_venue(L)

    def critical_price(self, venue):
        """C-87/C-90: the venue's OWN critical coupling price mu_c, LOCATED by computation,
           never a literal (D-8): Perron row-sum candidate 1/deg, (I - mu_c A) annihilates
           the constant vector, exact-rational resolvent SINGULAR at mu_c (the D-15 zero)
           and SOLVABLE at mu_c*(19/20), mu_c*(21/20) beside it (the positive controls; on
           the sealed deg-6 venue these are exactly the sealed 19/120 and 7/40).  Sealed
           values: 1/6 world D=3, 1/4 corner D=2, 1/2 chain D=1.  DECLINES on a venue that
           is not degree-regular.  Owners: Perron-Frobenius/Gershgorin standard.
           (LANE_T44_B_WORLD S0/S2; LANE_T44_A_CORNER S2, sealed.)"""
        import classes as CC
        return CC.mu_c_of(venue)

    def reachable_class(self, venue, mu, evidence=False):
        """C-87: THE REACHABLE-CLASS VERDICT -- the computed boolean triple (exponential,
           critical, divergent) for a declared coupling price mu, by exact rational
           comparison against the venue's own computed mu_c; exactly one True whenever mu_c
           is located.  evidence=True additionally runs the class's own instrument on the
           declared sector: subcritical booleans (ratios <= 1 - 1/20, Cauchy, power
           exclusion), the critical 1/d kernel signature (INV window, deepening-stabilized),
           or the divergence witness.  A price must be nonnegative; negative mu is refused
           before any verdict or evidence is emitted.  (LANE_T44_B_WORLD taxonomy +
           S2/S3/S4, sealed.)"""
        import classes as CC
        return CC.class_verdict(venue, mu, evidence=evidence)

    def coupling(self, mu, target, K):
        """C-87 class (1): the Gamma-priced coupling G_mu(d) = sum over admissible strings
           of mu^weight on the D=3 sector -- exact partial sum (Fraction) with the exact
           geometric tail (6mu)^{K+1}/(1-6mu); leading term N_min mu^d with w_min = d the
           confinement cost (C-80/O-54 standing).  Owners: walk generating functions,
           Spitzer/Lawler; comparison tier only.  Negative mu is refused rather than
           producing a signed 'tail bound'.  (t44b_lib.py series_3d, verbatim.)"""
        import classes as CC
        return CC.series_3d(mu, target, K)

    def critical_kernel(self, targets, M):
        """C-90: the regularized critical kernel a_M(x) = sum (N_2m(0)-N_2m(x))/36^m at the
           COMPUTED mu_c = 1/6, exact rationals, with certified tails beside it
           (classes.diff_tail_bound / abs_tail_bound; honest for M >~ 1000).  The D=3
           critical member is the 1/d POWER LAW: the sealed M=1400 exponent bracket contains
           1; separately, the M=2800 coefficient bracket [0.476369, 0.487321] contains owner
           3/(2 pi) (comparison).  See the fixed-distance caveat below.  Owners: Polya
           1921 transience, Watson 1939 G(0), Spitzer P26.1 coefficient -- comparisons only.
           (LANE_T44_B_WORLD S4 + register row C-90, sealed.)"""
        import classes as CC
        return CC.crit_kernel_3d(targets, M)

    def class_discriminator(self, K2=6000, K1=80000):
        """C-87: the cross-dimension discriminator -- ONE instrument (doubling-increment
           ratios of the regularized critical kernel), three venues, three pairwise-DISJOINT
           declared windows: D=1 LINEAR [9/5,11/5], D=2 LOG [4/5,5/4], D=3 INV [2/5,3/5].
           The critical class is the earned dimension's own; mu = 1/6 is critical on D=3 and
           subcritical on D=2.  (LANE_T44_B_WORLD S5, sealed.)"""
        import classes as CC
        return CC.discriminator(K2=K2, K1=K1)

    # ---------------------------------------------------------------- WRITING (T-54/T-55; C-91)
    # Every method below delegates to model/writing.py, where the sealed-lane machinery
    # lives (exact Fraction kernels, venues from carrier supports, the coset instrument
    # reused from geometry.py -- ported, not reimplemented).
    def writing_kernel_verdict(self, venue, c=0):
        """C-91 kernel tier: the invariant lazy-family kernel c*I + ((1-c)/deg)*A on a
           kernel venue ('C8' | 'T3' | 'T4' | 'Z27'), with conservation (double
           stochasticity), criticality (exact det(I - K)), link amplitudes, and the
           CTRL-LEAK control beside it.  The sealed identity: conserving <=> critical,
           per-crossing amplitude 1/deg identically for every nonzero crossing share
           (LANE_T48_A_DERIVATION, sealed).  The identity endpoint c=1 is conserving and
           critical but has no crossing, so per_crossing is None.  Owners: unital channels /
           Birkhoff (re-verified on explicit operators in-lane), Perron, Gershgorin."""
        import writing as WW
        from fractions import Fraction as Fr
        try:
            c = Fr(c)
        except (TypeError, ValueError, ZeroDivisionError) as exc:
            raise ValueError("writing_kernel_verdict requires a rational 0 <= c <= 1") from exc
        if c < 0 or c > 1:
            raise ValueError("writing_kernel_verdict requires 0 <= c <= 1")
        adj = {"C8": WW.ring_venue(8)["adj"], "T3": WW.plaquette_venue(3)["adj"],
               "T4": WW.plaquette_venue(4)["adj"], "Z27": WW.grain_venue(3)["adj"]}[venue]
        K = WW.kernel_uniform(adj, c)
        deg = WW.venue_degree(adj)
        conserving = WW.is_doubly_stochastic(K)
        return dict(deg=deg, conserving=conserving,
                    critical=bool(conserving and WW.crit_det(K) == 0),
                    link_amplitudes=WW.link_amplitudes(adj, K),
                    per_crossing=(Fr(1, deg) if c < 1 else None),
                    leak_det=WW.crit_det(WW.leak_kernel(adj, Fr(9, 10))))

    def writing_uniformity(self, Lx, Ly):
        """C-91 corner tier: uniformity EARNED from the writer algebra on the (Lx, Ly)
           plaquette venue rebuilt from carrier supports -- the elementary writer's
           invariant tuple on every link, mu_c located in-lane (Perron row sums + exact
           resolvent singular at 1/deg, nonsingular beside), and the unique conserving
           member t* = 1/deg = mu_c (LANE_T48_B_CORNER, sealed at (4,6), (3,7), (5,5)).
           Owners: carrier Kitaev quant-ph/9707021; coset instrument o54c lineage
           (geometry.py); Perron/Gershgorin/Feller as in-lane."""
        import writing as WW
        from fractions import Fraction as Fr
        cv = WW.corner_venue(Lx, Ly)
        inv = WW.writer_invariants(cv)
        return dict(invariant_tuple=inv[0], identical=all(v == inv[0] for v in inv),
                    n_links=len(inv),
                    mu_c=WW.mu_c_locate(cv["rows"], 0, beside=(Fr(1, 8), Fr(23, 100))),
                    conserving_member=WW.conserving_member(cv["rows"]))

    def writing_transport(self, n, a):
        """C-91 world tier, E1 TRANSPORT (iv' literal): the energy-conserving writer on the
           n^3 census venue at per-attempt amplitude a -- conserving AND critical
           (mu = 1/deg = mu_c) at EVERY dE and barrier; verdicts computed, never narrated
           after enforcing the probability domain 0 <= a <= 1/deg (invalid rows refuse);
           (LANE_T48_C_WORLD V1, sealed).  Owners: Stinespring (why row sums 1 is the
           dilation's structural property; the sums themselves computed)."""
        import writing as WW
        _c, _i, nbr = WW.torus3(n)
        return WW.transport_verdict(nbr, WW.ensemble_transport(nbr, a))

    def writing_trail_retreat(self, n, u, b):
        """C-91 world tier, E2 TRAIL WITH RETREAT: the raw two-rate writer with
           backtracking kept -- conserving at every dE, uniform (1/deg) exactly at dE = 0,
           split (b, 1)/(5b + 1) beside.  The constructor computes deg from the venue and
           refuses u < 0, b < 0, or (deg-1)*u*b + u > 1 before a verdict is formed
           (LANE_T48_C_WORLD V2, sealed).  Owners: Goldstein 1951 / Kac 1974 named for
           the persistence remark, comparison-only."""
        import writing as WW
        _c, _i, nbr = WW.torus3(n)
        return WW.retreat_verdict(nbr, WW.ensemble_trail_retreat(nbr, u, b))

    def writing_trail_decay(self, n, u, b, counting="H1"):
        """C-91 world tier, E3 TRAIL WITH DECAY (the model's own erase channel): NEVER
           critical -- mu = b/(deg*b + 1), f0 and E_b dropping out exactly; the COMPUTED
           mass ratio mu_c/mu is returned; the closed form ln(1 + e^{dE/kT}/l) is checked
           against it by the validator, never sourced from it (LANE_T48_C_WORLD V3/V4,
           sealed; the DONE_WHEN control).  counting='NB' uses the venue's own
           directed-edge criticality reference (Hashimoto 1989, earned by row sums).
           Both counting rules refuse negative or overfull probability rows; any other
           counting label is refused rather than silently treated as NB."""
        import writing as WW
        _c, _i, nbr = WW.torus3(n)
        W = WW.ensemble_trail_decay(nbr, u, b, counting)
        return (WW.decay_verdict(nbr, W) if counting == "H1"
                else WW.decay_verdict_nb(nbr, W))

    def writing_gap(self, s, n=4, den=10 ** 9):
        """C-91 + D-25, THE OBSERVATION ENTRY: a real record surface's written-trail mass
           gap.  REFUSES without provenance (build s through URM.surface); declines
           non-thermal surfaces and dials whose exponential underflows or whose declared
           rational denominator cannot resolve a positive lower bound.  Otherwise the
           surface's b = exp(-dE/kT) is bracketed by exact rationals, the E3 kernel is
           built and measured exactly at both brackets (u independence re-computed at
           entry), the closed form is checked against every computed ratio, and the float
           gap is certified INSIDE the computed bracket."""
        import writing as WW
        return WW.surface_gap(s, n=n, den=den)

    # ---------------------------------------------------------------- ROLES (honestly scoped)
    def roles(self):
        return {
            "EM": "supplies the carrier: in the corner, both boundary maps of the chain complex "
                  "(A-EM, FORMAL); in the world, every terrestrial census record's mechanism is "
                  "electromagnetic at bottom (census GR1). Restatement for the world = TD-2/T-15.",
            "GAMMA": "in the corner: record space, writer, protection via homology (A-GR, FORMAL; "
                     "C-74 exact on the torus). In the world: protection is by BARRIERS, of which "
                     "homology is the E_b -> infinity idealisation (census GR2). The world-level "
                     "emergence claim has no falsifiable statement yet = TD-1/T-37.",
            "ALPHA": "in the corpus: a generic coupling, nothing dimensionful (T-22, FORMAL). Under "
                     "the amended laws every terrestrial barrier is an EM energy scale, so alpha's "
                     "VALUE may enter record durability through E_b = TD-3/T-38, undecided.",
        }


# ---------------------------------------------------------------------- D-25: provenance registry
# The model is grounded in real record data, never the toy category (the principal, 2026-08-20).
# World-tier surfaces constructed through the URM's public gate MUST carry provenance; corner
# carriers must self-declare DEF-A. Sealed lanes predate this gate and use the raw class.
PROVENANCE = {
    "CoCrPt grain": "magnetic anisotropy; K_u ~ 2e5 J/m^3 class; stability rule K_uV/kT onset 35+/-2, "
                    "~60 for 10-yr: Weller & Moser, IEEE Trans. Magn. 35, 4423 (1999)",
    "CoCrPt HDD grain": "magnetic anisotropy; K_u ~ 2e5 J/m^3 class; stability rule K_uV/kT onset "
                        "35+/-2, ~60 for 10-yr: Weller & Moser, IEEE Trans. Magn. 35, 4423 "
                        "(1999); display-name alias of CoCrPt grain",
    "NAND floating gate": "trapped charge; detrapping E_a ~ 1.0 eV across generations; JEDEC 1-yr/30C "
                          "EOL spec (LANE_T41_EXTERNAL/CITATIONS.md)",
    "Azobenzene": "photoisomerisation; t1/2 = 1.4 d in benzene at 35 C, dH = 21.1 kcal/mol = 0.915 eV "
                  "(ACS Cent. Sci., PMC9951306)",
    "Azobenzene cis/trans": "photoisomerisation; t1/2 = 1.4 d in benzene at 35 C, dH = 21.1 "
                            "kcal/mol = 0.915 eV (ACS Cent. Sci., PMC9951306); display-name "
                            "alias of Azobenzene",
    "DNA base tautomer": "chemical bond; census GR1 values, source class: quantum-chemistry literature",
    "Fe(II) spin crossover": "spin-state transition; census GR1 values",
    "Alanine enantiomer": "parity violation; PVED ~ 1e-19..1e-17 eV literature class; instrument declines",
}


class URM(ProjectModel):
    """The public gate of the Universal Record Model.

       D-25 is enforced at both observation constructors: surface() REFUSES a world-tier
       surface without provenance, and corner() REFUSES an exact carrier that does not
       self-declare DEF-A.  New laws enter as ProjectModel layer methods plus validator
       gates and claim rows; new external numbers enter only as pinned comparisons with
       stated tolerance, semantics, and a power control."""

    @staticmethod
    def surface(name, mechanism, dE, E_b, T, f0, thermal=True, provenance=None, tier="world"):
        if tier == "corner":
            if provenance != "DEF-A":
                raise ValueError(
                    "URM REFUSES: a corner carrier must self-declare provenance='DEF-A' — the exact "
                    "idealisation may never silently pose as the world (D-25).")
        else:
            p = provenance or PROVENANCE.get(name)
            if not p or not str(p).strip():
                raise ValueError(
                    "URM REFUSES: a world-tier RecordSurface requires PROVENANCE — the real record it "
                    "models and its constants' pinned sources (D-25, the principal 2026-08-20: the "
                    "model is grounded in real record data, never the toy category). Register the "
                    "source in model.project_model.PROVENANCE or pass provenance=... explicitly.")
        s = RecordSurface(name, mechanism, dE, E_b, T, f0, thermal=thermal)
        s.provenance = provenance or PROVENANCE.get(name) or ("DEF-A" if tier == "corner" else None)
        s.tier = tier
        return s

    def corner(self, H, Ls=(), provenance=None):
        """D-25-gated DEF-A entry to the exact commutant corner."""
        if provenance != "DEF-A":
            raise ValueError(
                "URM REFUSES: an exact corner carrier must self-declare provenance='DEF-A' "
                "(D-25).")
        return super().corner(H, Ls=Ls)
