"""THE FULL-PROJECT MODEL. One object that represents the entire program.

The principal, 2026-08-20: "We need a model to work on that represents the full project."

Layers, each resting on registered rows:
  DEFINITION   clauses (i')-(v') with tolerances (O-51; census GR1-GR4)  -> ProjectModel.clauses
  LAWS         lifetime and steady value from the record's own Liouvillian mode
               (C-69, C-70 -- the two PROVED rows, six mechanisms)       -> .lifetime, .steady_value
  FORMATION    written vs unwritten configuration; the five source
               standards (C-71, C-72 -- promotion pending solidity)      -> .configuration
  CORNER       t_m -> infinity, W -> 0 reduces to the exact commutant
               records of record_model (C-75, T-28)                      -> .corner
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
        """(H, {L_k}) for the two-state record, rates by detailed balance over the barrier.
           Returns None where the surface is NOT thermally activated -- the model DECLINES rather
           than returning a number (T-33's zircon and CMB controls)."""
        if not self.thermal:
            return None
        kT = G.KB * self.T
        gd = self.f0 * np.exp(-(self.E_b + self.dE / 2) / kT)
        gu = self.f0 * np.exp(-(self.E_b - self.dE / 2) / kT)
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
