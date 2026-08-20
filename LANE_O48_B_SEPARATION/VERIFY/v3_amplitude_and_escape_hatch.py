"""ADVERSARIAL VERIFICATION of LANE_O48_B_SEPARATION, check 3.

THREE CLAIMS UNDER ATTACK.

 (A) "J_eff(r) = (-1)^{r+1} * 0.3146 g^2 / r, C = 0.3146"  and the log-log amplitudes
     "0.370346, 0.368511, 0.369040, 0.371385, 0.375254".  Both cannot be right, and on a RING
     (no open boundary at all) the amplitude is checkable against a closed form.

 (B) The finite-size collapse in [B3] is the lane's ONLY fit-free evidence for exponent 1.  Its
     loop is  r = int(round(frac*m)); r += (r % 2)  -- r is forced EVEN at every point, so the
     collapse never samples an ODD separation, i.e. never samples an opposite-sublattice pair.

 (C) THE ESCAPE HATCH: "if the records occupy only every SECOND mediator site ... the same
     induced 1/r interaction ACCUMULATES instead of cancelling", |sum J|/sum|J| = 1.000000000.
     That ratio is taken over the COUPLING CONSTANTS J(r).  The thing that would have to
     accumulate is the record-dependent ENERGY, whose two-body part is sum_{i<j} J_ij z_i z_j,
     and z_i z_j takes BOTH SIGNS.  Clause (iv) is exactly the statement that each record is
     balanced, so both values occur.  This script computes the ratio for the actual energy.
"""
import numpy as np

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)

def T_pbc(N, epsF=0.0, t=1.0):
    n = np.arange(N); k = 2 * np.pi * n / N; eps = -2.0 * t * np.cos(k)
    occ = np.where(eps < epsF - 1e-12)[0]; emp = np.where(eps > epsF + 1e-12)[0]
    assert len(occ) + len(emp) == N
    G = np.zeros(N)
    for p in occ:
        G += np.bincount((p - emp) % N, weights=1.0 / (eps[emp] - eps[p]), minlength=N)
    return np.real(np.fft.fft(G)) / (N * N)

P("=" * 122)
P("V3  AMPLITUDE, THE EVEN-ONLY COLLAPSE, AND WHETHER THE ESCAPE HATCH ACCUMULATES")
P("=" * 122)

N = 32770
T = T_pbc(N)
J = -8.0 * T                      # J_eff(r)/g^2 on the ring, EXACT at O(g^2)

P("")
P("[V3a] THE AMPLITUDE ON A RING (no open boundary, so no finite-size bend to blame).")
P("      LOCAL log-log slope of the envelope, and the envelope times r, deep in the clean window.")
P("-" * 122)
P(f"{'r':>7} {'S(r)=(-1)^(r+1)J':>20} {'S(r)*r':>12} {'local slope':>13} {'1/pi':>12} "
  f"{'S*r / (1/pi)':>14}")
def S(r): return ((-1) ** (r + 1)) * J[r]
for r in (8, 16, 32, 64, 128, 256, 512, 1024):
    sl = -(np.log(S(2 * r)) - np.log(S(r))) / np.log(2.0)
    P(f"{r:>7} {S(r):>20.9e} {S(r)*r:>12.9f} {sl:>13.6f} {1/np.pi:>12.9f} {S(r)*r*np.pi:>14.9f}")
P("")
P("      the lane reports C = 0.3146 (from f(r/m -> 0)) and C = 0.368-0.375 (from its log-log fits).")
P("      1/pi = %.9f." % (1 / np.pi))

P("")
P("[V3b] DOES THE COLLAPSE HOLD FOR ODD r?  The lane's [B3] loop forces r EVEN")
P("      (`r = int(round(frac*m)); r += (r % 2)`), so only same-sublattice pairs were sampled.")
P("      Here: S(r)*r at fixed r/N for ODD r and for EVEN r, on rings of four sizes.")
P("-" * 122)
rings = (2050, 8194, 32770)
Ts = {n: -8.0 * T_pbc(n) for n in rings}
def Sn(n, r): return ((-1) ** (r + 1)) * Ts[n][r]
P(f"{'r/N':>8} {'parity':>7} " + " ".join(f"{'N='+str(n):>14}" for n in rings) + f" {'spread %':>10}")
for frac in (0.01, 0.02, 0.05, 0.10):
    for par in (0, 1):
        vals = []
        for n in rings:
            r = int(round(frac * n))
            if r % 2 != par: r += 1
            vals.append(Sn(n, r) * r)
        P(f"{frac:>8.3f} {'even' if par==0 else 'odd':>7} " + " ".join(f"{v:>14.6f}" for v in vals)
          + f" {100*(max(vals)-min(vals))/np.mean(vals):>10.4f}")

P("")
P("[V3c] IS THE COEFFICIENT RECORD-BLIND?  (attack axis 5)")
P("-" * 122)
P("      T_ij is a functional of the BARE mediator's spectrum alone -- it is evaluated at g=0")
P("      background and contains no reference to the record configuration z.  The record content")
P("      of the energy is entirely the factor z_i z_j.  Stated plainly so [V3d] can be read.")

P("")
P("[V3d] THE ESCAPE HATCH, TESTED ON THE ENERGY RATHER THAN ON THE COUPLING CONSTANTS.")
P("      Records on every SECOND ring site (the lane's k=2 placement), K records, all separations")
P("      even, so every J_ij < 0 and the lane's ratio over J(r) is 1.000000000 by construction.")
P("      Here: ratio = |sum_{i<j} J_ij z_i z_j| / sum_{i<j} |J_ij|, for real record configurations.")
P("      D-15 CONTROLS in the same table: the ALIGNED configuration (the lane's implicit best case,")
P("      must give 1.000000) and the lane's own coupling-only ratio (must give 1.000000).")
P("-" * 122)
NR = 32770
TR = -8.0 * T_pbc(NR)
rng = np.random.default_rng(2026)
P(f"{'K records':>10} {'sum|J_ij|':>14} {'ALIGNED z':>13} {'ratio ALIGNED':>14} "
  f"{'median |E| over 64 random z':>29} {'ratio RANDOM (median)':>22} {'ratio, coupling-only':>21}")
for K in (16, 64, 256, 1024):
    pos = 2 * np.arange(K)
    d = np.abs(pos[:, None] - pos[None, :])
    Jm = np.where(d > 0, TR[np.clip(d, 0, NR - 1)], 0.0)
    iu = np.triu_indices(K, 1)
    Jp = Jm[iu]
    tot = np.abs(Jp).sum()
    z1 = np.ones(K)
    Ealign = abs((Jm @ z1 * z1).sum() - np.trace(Jm)) / 2.0
    Ealign = abs(0.5 * (z1 @ Jm @ z1))
    es = []
    for _ in range(64):
        z = rng.choice([-1.0, 1.0], size=K)
        es.append(abs(0.5 * (z @ Jm @ z)))
    med = float(np.median(es))
    # the lane's statistic: signs of the coupling constants alone, over separations
    seps = 2 * np.arange(1, K)
    Jr = TR[seps]
    lane_ratio = abs(Jr.sum()) / np.abs(Jr).sum()
    P(f"{K:>10} {tot:>14.6f} {Ealign:>13.6f} {Ealign/tot:>14.6f} {med:>29.6f} "
      f"{med/tot:>22.6f} {lane_ratio:>21.9f}")
P("")
P("      DOUBLING TEST (criterion (a) shape, applied to the same two quantities):")
P(f"{'K':>8} {'sum|J_ij|':>14} {'S(2K)/S(K) unsigned':>21} {'|E| median random':>19} {'S(2K)/S(K) random':>19}")
prev_tot = prev_med = None
for K in (64, 128, 256, 512, 1024, 2048):
    pos = 2 * np.arange(K)
    d = np.abs(pos[:, None] - pos[None, :])
    Jm = np.where(d > 0, TR[np.clip(d, 0, NR - 1)], 0.0)
    tot = np.abs(Jm[np.triu_indices(K, 1)]).sum()
    es = [abs(0.5 * (rng.choice([-1.0, 1.0], size=K) @ Jm @ rng.choice([-1.0, 1.0], size=K)))
          for _ in range(32)]
    es = []
    for _ in range(32):
        z = rng.choice([-1.0, 1.0], size=K); es.append(abs(0.5 * (z @ Jm @ z)))
    med = float(np.median(es))
    a = f"{tot/prev_tot:>21.6f}" if prev_tot else f"{'-':>21}"
    b = f"{med/prev_med:>19.6f}" if prev_med else f"{'-':>19}"
    P(f"{K:>8} {tot:>14.6f} {a} {med:>19.6f} {b}")
    prev_tot, prev_med = tot, med

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_B_SEPARATION/VERIFY/"
     "v3_amplitude_and_escape_hatch.txt", "w").write("\n".join(OUT) + "\n")
