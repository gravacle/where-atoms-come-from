"""THE UNIVERSAL RECORD MODEL (URM) — named by the principal, 2026-08-20.
One object that represents the entire program, and the overall representation of the proof:
PROOF_V002 is this model's narration. "Universal" is operational — one model, any record surface
(T-33: six mechanisms at machine precision; D-25: world-tier surfaces provenance-pinned) — and is a
title, not a claim.

The principal, 2026-08-20: "We need a model to work on that represents the full project."

Layers, each resting on registered rows:
  DEFINITION   clauses (i')-(v') with tolerances (O-51; census GR1-GR4)  -> ProjectModel.clauses
  LAWS         lifetime and steady value from the record's own Liouvillian mode
               (C-69, C-70 -- the two PROVED rows, six mechanisms)       -> .lifetime, .steady_value
  FORMATION    written vs unwritten configuration; the five source
               standards (C-71, C-72 -- promotion pending solidity)      -> .configuration
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
  ROLES        the three role sentences, honestly scoped                 -> .roles

A RecordSurface is what a physicist supplies: their own device or specimen's constants, in SI.
validate_project.py reproduces every sealed headline number through this one entry point."""
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

    # ---------------------------------------------------------------- LAWS (C-69, C-70; PROVED)
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
    "NAND floating gate": "trapped charge; detrapping E_a ~ 1.0 eV across generations; JEDEC 1-yr/30C "
                          "EOL spec (LANE_T41_EXTERNAL/CITATIONS.md)",
    "Azobenzene": "photoisomerisation; t1/2 = 1.4 d in benzene at 35 C, dH = 21.1 kcal/mol = 0.915 eV "
                  "(ACS Cent. Sci., PMC9951306)",
    "DNA base tautomer": "chemical bond; census GR1 values, source class: quantum-chemistry literature",
    "Fe(II) spin crossover": "spin-state transition; census GR1 values",
    "Alanine enantiomer": "parity violation; PVED ~ 1e-19..1e-17 eV literature class; instrument declines",
}


class URM(ProjectModel):
    """The public gate of the Universal Record Model. D-25 enforced here: surface() REFUSES a
       world-tier surface without provenance, and corner carriers must self-declare DEF-A."""

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
