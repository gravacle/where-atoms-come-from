"""
b0 -- SELF-TEST OF THE INSTRUMENT, RUN BEFORE ANY ARM.
Nothing below is a result about gauge theory.  It is the check that the algebra-entropy machinery
reproduces objects computable a second, independent way.  If this fails, every later number is void.
"""
import numpy as np, itertools
from lib_b import *

np.set_printoptions(precision=6, suppress=True)
fails = 0
def chk(name, a, b, tol=1e-8):
    global fails
    ok = abs(a - b) < tol
    if not ok: fails += 1
    print(f"  [{'ok ' if ok else 'FAIL'}] {name:<58} {a:.12f}   vs   {b:.12f}")

print("=" * 100)
print("b0  SELF-TEST")
print("=" * 100)

# ---- 1. algebra_entropy on a FULL matrix algebra must equal the reduced-density-matrix entropy
print("\n[1] full matrix algebra on a link set  ==  von Neumann entropy of the reduced state")
for car, seed in [(theta(6), 1), (theta(7), 2), (fan(3), 3)]:
    L = car["L"]
    psi = haar_physical(car, seed)
    E = pauli_table(psi, L)
    for R in [[0], [1], [0, 1], [0, 1, 2], [1, 2, 3]]:
        if max(R) >= L: continue
        S_alg, dims = algebra_entropy(E, gens_EXT(car, R))
        S_rdm = vn_entropy(reduce_links(psi, L, R))
        chk(f"{car['name']} R={R} (c,k)={dims}", S_alg, S_rdm)

# ---- 2. abelian algebra entropy must equal the Shannon entropy of the measured distribution
print("\n[2] abelian algebra alg{X_l} == Shannon entropy of the joint X-measurement distribution")
car = theta(6); L = 6
psi = haar_physical(car, 11)
E = pauli_table(psi, L)
for R in [[0], [0, 1], [0, 1, 2]]:
    S_alg, dims = algebra_entropy(E, gens_ELEC(car, R))
    # independent route: rotate to the X basis, take diagonal of the reduced density matrix
    Hd = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    U = np.array([[1.0]])
    for _ in range(L): U = np.kron(U, Hd)
    phi = U @ psi
    p = np.diag(reduce_links(phi, L, R)).real
    chk(f"R={R} (c,k)={dims}", S_alg, shannon(p))

# ---- 3. mutual information via algebras with C1 == ordinary extended-Hilbert-space MI
print("\n[3] I_{C1}(S:F) == ordinary mutual information of the reduced states")
for car, seed in [(theta(7), 5), (fan(3), 6)]:
    L = car["L"]; psi = haar_physical(car, seed); E = pauli_table(psi, L)
    for S_, F_ in [([0], [1]), ([0], [1, 2]), ([0], [1, 2, 3])]:
        if max(F_) >= L: continue
        I_alg, _ = mutual_information(E, gens_EXT(car, S_), gens_EXT(car, F_))
        chk(f"{car['name']} S={S_} F={F_}", I_alg, mi_ext(psi, L, S_, F_))

# ---- 4. GHZ sanity: I(S:F) = 1 bit exactly on every proper fragment
print("\n[4] theta_7 GHZ (= g2->0 ground state):  I_ext(link0 : F) for |F| = 1..6")
car = theta(7); L = 7
ghz = sym_basis_state(car, 0)
print("      GHZ overlap with (|0..0>+|1..1>)/sqrt2 :",
      float(abs(ghz[0]) ** 2 + abs(ghz[(1 << L) - 1]) ** 2))
for m in range(1, L):
    F = list(range(1, 1 + m))
    print(f"      |F|={m}  I_ext = {mi_ext(ghz, L, [0], F):.12f}")

# ---- 5. the g2 -> 0 ground state really is the GHZ state (checks the Hamiltonian code)
print("\n[5] ground state of H at small g2 vs GHZ, theta_7")
psi0, e0 = ground_state(car, 0.05)
print(f"      |<GHZ|psi_0>|^2 = {abs(np.vdot(ghz, psi0))**2:.12f}   E0 = {e0:.6f}")
print(f"      physical?  ||psi - P psi|| = {np.linalg.norm(psi0 - project_physical(psi0, car)):.3e}")

# ---- 6. physical-sector dimensions
print("\n[6] physical sector dimensions")
for car in [theta(3), theta(6), theta(7), theta(9), fan(3), fan(4)]:
    print(f"      {car['name']:<10} L={car['L']:<3} indep Gauss ops={len(indep_gauss(car)):<3} "
          f"dim H_ext={1<<car['L']:<6} dim H_phys={phys_dim(car)}")

print("\n" + "=" * 100)
print(f"b0 SELF-TEST  {'PASSED' if fails == 0 else 'FAILED with ' + str(fails) + ' mismatches'}")
print("=" * 100)
