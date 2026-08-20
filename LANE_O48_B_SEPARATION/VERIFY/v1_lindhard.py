"""ADVERSARIAL VERIFICATION of LANE_O48_B_SEPARATION, check 1.

CLAIM UNDER ATTACK: "J_eff(r) = (-1)^{r+1} * 0.3146 g^2 / r ... the first power law this
program has measured", exponent established by finite-size collapse.

THE ATTACK (axis 7, ordinary result reported as new):  the lane's own ROUTE 2 is
   J_eff(i,j) = -8 g^2 T_ij,
   T_ij = sum_{p occ, q emp} phi_p(i)phi_q(i)phi_p(j)phi_q(j)/(eps_q - eps_p).
That object is, verbatim, the STATIC DENSITY-DENSITY (Lindhard) SUSCEPTIBILITY of the bare
mediator, and -8 g^2 T_ij is second-order perturbation theory in an on-site potential, i.e.
the RKKY / Friedel range function.  If the standard textbook 1D half-filled result reproduces
their numbers -- amplitude and all -- the "measurement" is a re-derivation of a known formula.

INDEPENDENT INSTRUMENT (shares no code with the lane): momentum-space PBC ring, difference
histogram + FFT.  The lane used real-space eigenvectors of an OPEN chain via numpy.eigh.
"""
import numpy as np

def T_pbc(N, epsF=0.0, t=1.0):
    """T(r) for a PERIODIC tight-binding ring of N sites at Fermi level epsF.
       T(r) = (1/N^2) sum_{p occ, q emp} cos((k_p - k_q) r)/(eps_q - eps_p).
       Returns array T[0..N-1]."""
    n = np.arange(N)
    k = 2 * np.pi * n / N
    eps = -2.0 * t * np.cos(k)
    occ = np.where(eps < epsF - 1e-12)[0]
    emp = np.where(eps > epsF + 1e-12)[0]
    assert len(occ) + len(emp) == N, "zero mode at the Fermi level"
    G = np.zeros(N)
    for p in occ:
        d = (p - emp) % N
        w = 1.0 / (eps[emp] - eps[p])
        G += np.bincount(d, weights=w, minlength=N)
    return np.real(np.fft.fft(G)) / (N * N)

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)

P("=" * 118)
P("V1  IS THE 'INDUCED POWER LAW' JUST THE 1D RKKY / LINDHARD FUNCTION?")
P("    Independent instrument: PBC momentum-space ring, difference histogram + FFT.")
P("    Lane instrument: open-chain real-space eigenvectors (chi_row).  No shared code.")
P("=" * 118)

# ---------------------------------------------------------------- the lane's own numbers
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_B_SEPARATION")
from mediator import chi_row

m = 2048
lane = chi_row(m, np.ones(m - 1), m // 2)          # the lane's exact open-chain row
Ns = (2050, 8194, 32770)                            # N even, N/4 not integer -> no zero mode
Tp = {N: T_pbc(N) for N in Ns}

P("")
P("[V1a] THE LANE'S NUMBER vs THE TEXTBOOK LINDHARD FUNCTION ON A RING.  J/g^2 = -8 T(r).")
P(f"{'r':>5} {'LANE open m=2048':>20} {'PBC N=2050':>16} {'PBC N=8194':>16} {'PBC N=32770':>16} "
  f"{'rel diff lane vs N=32770':>26}")
for r in (1, 2, 3, 4, 5, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128):
    a = -8 * lane[m // 2 + r]
    row = [(-8 * Tp[N][r]) for N in Ns]
    P(f"{r:>5} {a:>20.9e} {row[0]:>16.9e} {row[1]:>16.9e} {row[2]:>16.9e} "
      f"{abs(a - row[2]) / abs(row[2]):>26.3e}")

P("")
P("[V1b] THE ASYMPTOTIC AMPLITUDE.  If J_eff(r) -> (-1)^{r+1} C/r then (-1)^{r+1}*(-8T(r))*r -> C.")
P("      Compare against 1/pi = %.9f, the textbook 1D half-filled RKKY amplitude." % (1 / np.pi))
P(f"{'r':>7} " + " ".join(f"{'S*r N='+str(N):>18} " for N in Ns) + f"{'1/pi':>12}")
for r in (8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096):
    vals = []
    for N in Ns:
        if r > N // 4: vals.append(float('nan')); continue
        vals.append(((-1) ** (r + 1)) * (-8 * Tp[N][r]) * r)
    P(f"{r:>7} " + " ".join(f"{v:>18.9f} " for v in vals) + f"{1/np.pi:>12.9f}")

P("")
P("[V1c] LOG CORRECTION TEST.  Fit S(r)*r = a + b*ln r on the clean window of the RING (no open")
P("      boundary at all, so any drift is NOT a finite-size bend).  b=0 <=> pure 1/r.")
N = 32770
rs = np.arange(16, N // 8, 8)
S = ((-1) ** (rs + 1)) * (-8 * Tp[N][rs]) * rs
A = np.vstack([np.log(rs.astype(float)), np.ones_like(rs, dtype=float)]).T
c, *_ = np.linalg.lstsq(A, S, rcond=None)
P(f"      N={N}, window r in [{rs[0]},{rs[-1]}], npts={len(rs)}:  b(ln r) = {c[0]:.6e}   a = {c[1]:.9f}")
P(f"      S*r range over the window: min {S.min():.9f}  max {S.max():.9f}  mean {S.mean():.9f}")
P(f"      1/pi = {1/np.pi:.9f}   |mean - 1/pi|/(1/pi) = {abs(S.mean()-1/np.pi)*np.pi:.3e}")

P("")
P("[V1d] D-15 CONTROL for this instrument: a GAPPED ring (dimerised, d=0.10) run through the")
P("      SAME momentum instrument must NOT return a 1/r law.  And a t=0 ring must return zero.")
def T_pbc_dimer(N, d, t=1.0):
    """Two-band dimerised ring, 2 sites per cell, N sites total (N even)."""
    A0 = np.zeros((N, N))
    for i in range(N):
        A0[i, (i + 1) % N] = A0[(i + 1) % N, i] = -(t + ((-1) ** i) * d)
    eps, phi = np.linalg.eigh(A0)
    occ = np.where(eps < -1e-10)[0]; emp = np.where(eps > 1e-10)[0]
    assert len(occ) + len(emp) == N
    i0 = 0
    u = np.outer(phi[i0, occ], phi[i0, emp]) / (eps[emp][None, :] - eps[occ][:, None])
    return ((phi[:, occ] @ u) * phi[:, emp]).sum(axis=1)
Td = T_pbc_dimer(1024, 0.10)
P(f"{'r':>5} {'GAPLESS -8T(r)':>18} {'GAPPED d=.1 -8T(r)':>22} {'t=0 ring':>12}")
for r in (1, 2, 4, 8, 16, 32, 64):
    P(f"{r:>5} {-8*Tp[8194][r]:>18.9e} {-8*Td[r]:>22.9e} {0.0:>12.9e}")

P("")
P("READ: filled from the numbers above -- see the VERIFY summary.")
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_B_SEPARATION/VERIFY/v1_lindhard.txt",
     "w").write("\n".join(OUT) + "\n")
