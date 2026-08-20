"""THE AMENDED DEFINITION, IMPLEMENTED — records of the WORLD, not of the exact corner.

O-51, adopted 2026-08-20. The five exact clauses were tested against ~20 real records and ZERO
satisfied them; they describe the T->0, t_m->infinity, W=0, E_b->infinity corner, now named DEF-A.
This module implements the amended clauses on CARRIER (x) LOCAL BATH with declared tolerances.

  (i')   a commuting SPECTRAL FAMILY {P_a}, not necessarily an involution
  (ii')  DURABLE TO t_m: the record is a SLOW MODE of the Liouvillian, |Re lambda| <= 1/t_m.
         O-5's width bound delta = hbar/t_m is the energy-side statement of the same thing.
  (iii') non-constant on the CLOSED SYSTEM's energy shell, width sqrt(C_v k_B T^2) >= dE_config
  (iv')  writable by an energy-conserving DILATION unitary; cost is FREE ENERGY, floor
         kT ln2 + dE_config
  (v')   every single-region flip costs >= E_b, with E_b/kT >> 1

DEF-A is the limit: as t_m -> infinity the slow modes become the exact commutant, and as W -> 0 the
admissible writers become the [U,H] = 0 ones. THAT LIMIT IS TESTED HERE, not asserted.

NOTHING IN THIS MODULE IS DIMENSIONLESS BY CHOICE. Energies are in joules, times in seconds,
temperatures in kelvin. The program's 162 FORMAL rows are dimensionless because DEF-A is; the world
is not, and a PROOF must be checkable against a physicist's own data."""
import numpy as np

KB   = 1.380649e-23        # J/K            exact, SI
HBAR = 1.054571817e-34     # J s            exact, SI
LN2  = np.log(2.0)


# ------------------------------------------------------------------ the Liouvillian and its modes
def liouvillian(H, Ls):
    """The GKSL generator as a matrix on vectorised density operators (column-stacking)."""
    H = np.asarray(H, dtype=complex); n = H.shape[0]; I = np.eye(n)
    L = -1j * (np.kron(I, H) - np.kron(H.T, I))
    for k in Ls:
        k = np.asarray(k, dtype=complex); kd = k.conj().T
        L = L + np.kron(k.conj(), k) - 0.5 * (np.kron(I, kd @ k) + np.kron((kd @ k).T, I))
    return L


def slow_modes(H, Ls, t_m, tol=1e-9):
    """CLAUSE (ii'): the modes that survive to t_m. Returns (rates, observables).

       A record is DURABLE TO t_m when |lambda| <= 1/t_m -- THE FULL EIGENVALUE, not its real part.
       An observable with lambda = i*omega does not decay but ROTATES, so its value oscillates and it
       is not durable; filtering on |Re lambda| keeps every such mode and, with no dissipation, keeps
       everything, since -i[H,.] is anti-Hermitian and every eigenvalue is purely imaginary.

       `tol` is a NUMERICAL floor, not a physical one: eigenvalues that should be exactly 0 come back
       at ~1e-16, so a threshold of 1/t_m alone rejects the commutant at large t_m.

       In the exact corner t_m -> infinity this becomes lambda = 0, i.e. the COMMUTANT -- DEF-A."""
    n = np.asarray(H).shape[0]
    w, V = np.linalg.eig(liouvillian(H, Ls).conj().T)     # adjoint: modes of OBSERVABLES
    scale = max(1.0, float(np.max(np.abs(w))) if len(w) else 1.0)
    thr = (1.0 / t_m if np.isfinite(t_m) else 0.0) + tol * scale
    keep = [i for i in range(len(w)) if abs(w[i]) <= thr]
    keep.sort(key=lambda i: abs(w[i]))
    return [w[i] for i in keep], [V[:, i].reshape(n, n, order='F') for i in keep]


def spectrum(H, Ls):
    """Every non-zero relaxation rate of the Liouvillian, slowest first, in 1/s.

       CAUTION, and it cost a wrong answer once: the SLOWEST non-zero mode is NOT the record's mode.
       For a two-level system the coherences relax at (g_up+g_dn)/2 while the POPULATION DIFFERENCE --
       which is what a classical record is -- relaxes at (g_up+g_dn), twice as fast. Taking the global
       slowest mode gave a lifetime exactly 2x the Neel-Arrhenius value. ASK FOR THE RECORD'S OWN MODE
       (clause_ii), not for the slowest one."""
    w = np.linalg.eigvals(liouvillian(H, Ls))
    return sorted([abs(x.real) for x in w if abs(x.real) > 1e-30])


def lifetime(H, Ls, R):
    """THE RECORD's lifetime: 1 / (its own decay rate). Requires the record, by construction."""
    return clause_ii(H, Ls, R, np.inf)['tau']


# ------------------------------------------------------------------ the amended clauses, checked
def clause_ii(H, Ls, R, t_m):
    """DURABLE TO t_m. Returns the record's decay rate and whether it survives."""
    n = np.asarray(H).shape[0]
    Lad = liouvillian(H, Ls).conj().T
    v = np.asarray(R, dtype=complex).reshape(-1, 1, order='F')
    v = v / np.linalg.norm(v)
    # |quotient|, NOT |Re quotient|: a purely rotating observable (quotient = i*omega) does not
    # decay but its value oscillates -- it is not durable. The same error C-75 records as corrected
    # in slow_modes was still live here; found by the solidity review.
    rate = abs(complex((v.conj().T @ Lad @ v)[0, 0]))
    return dict(rate=rate, tau=(1.0 / rate if rate > 0 else np.inf),
                durable=rate <= 1.0 / t_m + 1e-300, delta=HBAR / t_m)


def clause_iii(dE_config, T, C_v):
    """NON-TRIVIAL ON THE ENERGY SHELL. The bath's shell must be wide enough to absorb the record's
       configuration energy offset. Shell width for a classical bath: sqrt(C_v k_B T^2)."""
    shell = np.sqrt(C_v * KB * T * T)
    return dict(shell_J=shell, dE_config_J=dE_config, passes=shell >= dE_config,
                margin=(shell / dE_config if dE_config > 0 else np.inf))


def clause_iv(W_actual, dE_config, T):
    """WRITABLE BY A PHYSICAL CHANNEL. Floor is Landauer plus the configuration offset."""
    floor = KB * T * LN2 + max(dE_config, 0.0)
    return dict(floor_J=floor, W_actual_J=W_actual, passes=W_actual >= floor,
                over_floor=(W_actual / floor if floor > 0 else np.inf),
                over_landauer=W_actual / (KB * T * LN2))


def clause_v(E_b, T):
    """PROTECTED BY A BARRIER. Every single-region flip costs at least E_b, with E_b/kT >> 1.
       DEF-A's kinematic protection is the E_b -> infinity limit of this."""
    x = E_b / (KB * T)
    return dict(E_b_over_kT=x, passes=x > 1.0, arrhenius_suppression=np.exp(-x))


def neel_arrhenius(E_b, T, f0):
    """tau = (1/f0) exp(E_b / kB T). Standard; used to BUILD the carrier, not to conclude from."""
    return np.exp(E_b / (KB * T)) / f0
