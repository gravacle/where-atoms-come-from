# W19-A REFUTER / step 1.  IS THE INSTRUMENT WHAT IT CLAIMS?
#
# Lane A builds the physical sector by GAUGE FIXING to a spanning tree: it never touches the full
# N^L space when it diagonalises.  The lane's own validation (out_00 block 2) only checks that the
# LIFTED vector is Gauss-invariant -- which is guaranteed by construction of the lift and therefore
# COULD NOT HAVE FAILED.  It is not a check that the HAMILTONIAN was projected correctly.
#
# This script builds the physical sector a SECOND, INDEPENDENT WAY:
#   * enumerate the full N^L Z-basis,
#   * partition it into gauge ORBITS by explicit application of the Gauss permutations G_v,
#   * take the orbit-uniform vectors as an orthonormal basis of the +1 eigenspace of every G_v
#     (G_v are permutations, so their common +1 eigenspace is exactly the functions constant on orbits),
#   * build H in the FULL space from the brief's formula, compress it to that basis, diagonalise.
# Then compare the FULL SPECTRUM with lane A's gauge-fixed spectrum.  A wrong projection, a wrong
# plaquette set, a mis-signed cycle, or a missing electric term would all show up as a spectrum gap.
#
# Also: the eigenvector residual  ||H_full Psi - E0 Psi||  on the 21-link carrier itself (2^21 = 2097152),
# computed by matvec without ever forming a matrix.  THAT is the check the lane did not run.
import numpy as np, sys, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_A_SWEEP")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_A_SWEEP_REFUTE")
from zn_gauge import ZNGauge, S_of, mutual_information, level_cuts, nested_fragments, girth_through
from carriers import theta, ladder, cube, grid, petersen, heawood, theta_sub
from refute_carriers import tri_chain12, dbl_chain9

np.set_printoptions(precision=9, suppress=True)
print("=" * 110)
print("W19-A REFUTER / 1  INDEPENDENT RECONSTRUCTION OF THE PHYSICAL SECTOR")
print("=" * 110)

def digits(N, L):
    idx = np.arange(N ** L, dtype=np.int64)
    return np.stack([(idx // (N ** (L - 1 - e))) % N for e in range(L)], axis=1)   # link 0 = MSD

def orbit_basis(N, V, edges):
    """Independent physical sector: orbits of the Z-basis under z -> z + d.lambda."""
    L = len(edges); D = N ** L
    z = digits(N, L)
    pw = np.array([N ** (L - 1 - e) for e in range(L)], dtype=np.int64)
    # canonical orbit label: minimum index over the whole gauge group N^V (brute force, small L only)
    lab = np.arange(D, dtype=np.int64)
    for lam in itertools.product(range(N), repeat=V):
        lam = np.array(lam, dtype=np.int64)
        shift = np.array([(lam[b] - lam[a]) % N for (a, b) in edges], dtype=np.int64)
        znew = (z + shift) % N
        j = znew @ pw
        lab = np.minimum(lab, j)          # NOT a full union-find, but the group acts, so min over
                                          # the whole group from every point IS the orbit canonical rep
    reps = np.unique(lab)
    return lab, reps

def H_full_dense(N, V, edges, plaq, mag, elec):
    L = len(edges); D = N ** L
    z = digits(N, L)
    H = np.zeros((D, D))
    diag = np.zeros(D)
    for p in plaq:
        ph = (z @ np.array(p, dtype=np.int64)) % N
        diag += -(mag / 2.0) * 2.0 * np.cos(2 * np.pi * ph / N)
    np.fill_diagonal(H, diag)
    pw = np.array([N ** (L - 1 - e) for e in range(L)], dtype=np.int64)
    for e in range(L):
        for s in (+1, -1):
            zz = z.copy(); zz[:, e] = (zz[:, e] + s) % N
            j = zz @ pw
            H[j, np.arange(D)] += -(elec / 2.0)
    return H

print("\n[1a] FULL SPECTRUM COMPARISON.  arm 1 = lane A's gauge-fixed orbit-basis H (zn_gauge.ZNGauge).")
print("     arm 2 = this script's full N^L H compressed onto the orbit-uniform basis.  SAME g^2, SAME")
print("     plaquette set (taken FROM lane A's object, so a plaquette-set error would be shared, not hidden).")
print(f"     {'carrier':<16}{'N':>3}{'L':>4}{'C':>4}{'dimP(lane)':>12}{'dimP(orbits)':>14}"
      f"{'max|dE|':>12}{'E0 lane':>15}{'E0 indep':>15}")
CASES = [("theta", 2, theta()), ("theta_subdiv2", 2, theta_sub(2)), ("ladder_2sq", 2, ladder(2)),
         ("ladder_2sq", 3, ladder(2)), ("dbl_chain9", 2, dbl_chain9()),
         ("cube_Q3", 2, cube()), ("tri_chain12", 2, tri_chain12()), ("grid_3x3_open", 2, grid(3, 3, False))]
for nm, N, (V, E) in CASES:
    L = len(E)
    if N ** L > 5000:
        print(f"     {nm:<16}{N:>3}{L:>4}  -- {N}^{L} too big for the DENSE independent check, skipped here"); continue
    g = ZNGauge(nm, V, E, N)
    Hlane = g.hamiltonian(2.0, 2.0)
    wl = np.linalg.eigvalsh(Hlane)
    lab, reps = orbit_basis(N, V, E)
    D = N ** L
    B = np.zeros((D, len(reps)))
    for k, r in enumerate(reps):
        m = (lab == r); B[m, k] = 1.0 / np.sqrt(m.sum())
    Hf = H_full_dense(N, V, E, g.plaq, 2.0, 2.0)
    Hp = B.T @ Hf @ B
    wi = np.linalg.eigvalsh(Hp)
    # the two spectra must agree as multisets
    n = min(len(wl), len(wi))
    dmax = float(np.abs(np.sort(wl)[:n] - np.sort(wi)[:n]).max()) if len(wl) == len(wi) else float("nan")
    print(f"     {nm:<16}{N:>3}{L:>4}{g.C:>4}{g.dimP:>12}{len(reps):>14}{dmax:>12.3e}{wl[0]:>15.9f}{wi[0]:>15.9f}")
    # also: is the orbit basis really Gauss-invariant, and is B^T B = I ?
    assert abs((B.T @ B - np.eye(len(reps))).max()) < 1e-12

print("\n[1b] IS THE COMPRESSION LOSSLESS?  Does [H_full, P_phys] = 0, i.e. is the physical sector")
print("     actually H-invariant?  If it were not, lane A's whole spectrum would be a fiction.")
for nm, N, (V, E) in [("theta", 2, theta()), ("ladder_2sq", 2, ladder(2)), ("dbl_chain9", 2, dbl_chain9())]:
    g = ZNGauge(nm, V, E, N); L = len(E); D = N ** L
    lab, reps = orbit_basis(N, V, E)
    B = np.zeros((D, len(reps)))
    for k, r in enumerate(reps):
        m = (lab == r); B[m, k] = 1.0 / np.sqrt(m.sum())
    P = B @ B.T
    Hf = H_full_dense(N, V, E, g.plaq, 2.0, 2.0)
    print(f"     {nm:<16} ||[H_full,P]||_inf = {np.abs(Hf@P - P@Hf).max():.3e}   "
          f"rank P = {int(round(np.trace(P)))} (= dim_phys {g.dimP})")

print("\n[1c] THE CHECK LANE A DID NOT RUN: is the lifted 21-link ground state an EIGENVECTOR of the")
print("     FULL 2^21 Hamiltonian with lane A's reported E0?  matvec only, no matrix formed.")
def H_apply(g, Psi, mag, elec):
    N, L = g.N, g.L
    T = Psi.reshape([N] * L)
    out = np.zeros_like(T)
    # magnetic: diagonal phase
    idx = np.arange(N ** L, dtype=np.int64)
    dig = [((idx // (N ** (L - 1 - e))) % N).astype(np.int64) for e in range(L)]
    ph = np.zeros(N ** L, dtype=np.int64)
    diag = np.zeros(N ** L)
    for p in g.plaq:
        ph[:] = 0
        for e in range(L):
            if p[e]: ph = (ph + int(p[e]) * dig[e]) % N
        diag += -(mag / 2.0) * 2.0 * np.cos(2 * np.pi * ph / N)
    out += (diag.reshape([N] * L) * T)
    for e in range(L):
        for s in (+1, -1):
            out += -(elec / 2.0) * np.roll(T, s, axis=e)
    return out.ravel()

for nm, (V, E), gsq in [("heawood", heawood(), 0.5), ("heawood", heawood(), 1.0), ("petersen", petersen(), 1.0)]:
    g = ZNGauge(nm, V, E, 2)
    psi, E0, gap = g.ground(2.0 / gsq, 2.0 * gsq)
    Psi = g.full_vector(psi)
    R = H_apply(g, Psi, 2.0 / gsq, 2.0 * gsq) - E0 * Psi
    print(f"     {nm:<10} g^2={gsq:<4} L={g.L}  2^L={2**g.L}  E0={E0:.9f}  gap={gap:.9f}  "
          f"||H_full Psi - E0 Psi||_2 = {np.linalg.norm(R):.3e}  <Psi|H_full|Psi> = "
          f"{float(Psi @ H_apply(g, Psi, 2.0/gsq, 2.0*gsq)):.9f}")

print("\n[1d] ORBIT COUNT / MULTIPLICITY.  Each gauge orbit must have exactly N^(V-1) members for")
print("     Psi = psi(Gamma z)/N^((V-1)/2) to be normalised.  Counting them on the small carriers:")
for nm, N, (V, E) in [("theta", 2, theta()), ("ladder_2sq", 2, ladder(2)), ("ladder_2sq", 3, ladder(2)),
                      ("dbl_chain9", 2, dbl_chain9())]:
    lab, reps = orbit_basis(N, V, E)
    sizes = np.bincount(np.searchsorted(reps, lab))
    print(f"     {nm:<14} Z_{N}: orbits={len(reps)}  orbit sizes: min={sizes.min()} max={sizes.max()} "
          f"expected N^(V-1)={N**(V-1)}  -> uniform: {sizes.min()==sizes.max()==N**(V-1)}")
print("\nDONE 1.")
