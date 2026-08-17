"""
r0 -- SELF-TEST OF THE REFUTATION INSTRUMENT, RUN BEFORE ANY CLAIM IS BELIEVED.
Nothing here is a result.  It is the check that the four things this lane ADDS to the target's
own lib_b -- mixed-state expectations, mixed-state partial traces, exhaustive subspace
enumeration, exact disjoint-fragment packing, and basis rotations -- reproduce objects computable
a second, independent way.  If this fails, every number in r1-r4 is void.

The target's own instrument is imported UNMODIFIED and is not re-tested here; it was validated in
its own b0 and, separately, all five of its scripts were re-run in this lane and reproduce their
sealed outputs byte-identically (recorded in REPRODUCTION.txt).
"""
import itertools
import numpy as np
from rlib import *

fails = 0
def chk(name, a, b, tol=1e-9):
    global fails
    ok = abs(a - b) < tol
    if not ok: fails += 1
    print(f"  [{'ok ' if ok else 'FAIL'}] {name:<62} {a:.12f}  vs  {b:.12f}")

def chkint(name, a, b):
    global fails
    ok = (a == b)
    if not ok: fails += 1
    print(f"  [{'ok ' if ok else 'FAIL'}] {name:<62} {a}  vs  {b}")

hr("r0  SELF-TEST OF THE REFUTATION INSTRUMENT")

# ---- 1  mixed-state partial trace vs the target's pure-state partial trace
print("\n[1] reduce_rho(|psi><psi|) == lib_b.reduce_links(psi)   (my trace vs theirs)")
for car, seed in [(theta(6), 1), (theta(7), 2), (fan(3), 3)]:
    L = car["L"]; psi = haar_physical(car, seed); rho = np.outer(psi, psi.conj())
    for R in [[0], [0, 1], [1, 2, 3]]:
        chk(f"{car['name']} R={R}", vn_entropy(reduce_rho(rho, L, R)),
            vn_entropy(reduce_links(psi, L, R)))

# ---- 2  MixE: algebra entropies of a MIXED state against a direct density-matrix route
print("\n[2] algebra_entropy(MixE, full matrix algebra on R) == vn_entropy of the mixed rho_R")
for car, s1, s2, w in [(theta(6), 11, 12, 0.5), (theta(7), 21, 22, 0.3), (fan(3), 31, 32, 0.75)]:
    L = car["L"]
    p1 = haar_physical(car, s1); p2 = haar_physical(car, s2)
    rho = rho_from_mixture([p1, p2], [w, 1 - w])
    Em = MixE([p1, p2], [w, 1 - w], L)
    for R in [[0], [0, 1], [0, 1, 2]]:
        chk(f"{car['name']} w={w} R={R}", algebra_entropy(Em, gens_EXT(car, R))[0],
            vn_entropy(reduce_rho(rho, L, R)))

print("\n[2b] MixE on an ABELIAN algebra == Shannon entropy of the X-measurement distribution")
car = theta(6); L = 6
p1 = haar_physical(car, 41); p2 = haar_physical(car, 42)
rho = rho_from_mixture([p1, p2], [0.5, 0.5])
Em = MixE([p1, p2], [0.5, 0.5], L)
Hd = np.array([[1, 1], [1, -1]]) / np.sqrt(2); Ufull = np.array([[1.0]])
for _ in range(L): Ufull = np.kron(Ufull, Hd)
rhoX = Ufull @ rho @ Ufull.conj().T
for R in [[0], [0, 1], [0, 1, 2]]:
    p = np.diag(reduce_rho(rhoX, L, R)).real
    chk(f"abelian alg{{X_l}} R={R}", algebra_entropy(Em, gens_ELEC(car, R))[0], shannon(p))

print("\n[2c] MixE with a single component == the target's pure-state PauliExpect")
psi = haar_physical(theta(7), 55)
E1 = pauli_table(psi, 7); E2 = MixE([psi], [1.0], 7)
worst = max(abs(E1(x, z) - E2(x, z)) for x in range(128) for z in range(128))
chk("max |MixE - PauliExpect| over all 4^7 strings", worst, 0.0, 1e-12)

# ---- 3  exhaustive subspace enumeration against the Galois numbers
print("\n[3] all_subspaces(n) counts == Galois numbers G_n (complete enumeration, not a sample)")
for n in range(1, 7):
    chkint(f"n={n}", len(all_subspaces(n)), GALOIS[n])
print("     and every enumerated set is genuinely closed under XOR:")
bad = 0
for n in (3, 4):
    for s in all_subspaces(n):
        for u in s:
            for v in s:
                if (u ^ v) not in s: bad += 1
chkint("closure violations over n=3,4", bad, 0)

# ---- 4  exact disjoint packing, against hand-computable cases
print("\n[4] rdelta_packing: maximum number of pairwise-disjoint qualifying fragments")
chkint("N=6, all three matched pairs {0,1}{2,3}{4,5} qualify",
       rdelta_packing([0b000011, 0b001100, 0b110000], 6), 3)
chkint("N=6, only ONE fragment qualifies", rdelta_packing([0b000011], 6), 1)
chkint("N=6, overlapping pairs {0,1}{1,2}{2,3} -> at most 2 disjoint",
       rdelta_packing([0b000011, 0b000110, 0b001100], 6), 2)
chkint("N=7, all seven singletons qualify", rdelta_packing([1 << i for i in range(7)], 7), 7)
chkint("N=7, nothing qualifies", rdelta_packing([], 7), 0)
chkint("N=8, four disjoint plaquettes qualify",
       rdelta_packing([0b11, 0b1100, 0b110000, 0b11000000], 8), 4)

# ---- 5  basis rotations are unitary and act where they say they act
print("\n[5] rot_link / rot_block: unitarity, and agreement with an explicit kron")
psi = haar_physical(theta(6), 77); L = 6
V = bloch_basis(0.7, 1.3)
r = rot_link(psi, L, 2, V)
chk("rot_link preserves norm", float(np.linalg.norm(r)), 1.0)
I2 = np.eye(2, dtype=complex); M = np.array([[1.0]])
for l in range(L - 1, -1, -1):                       # link l is bit l => leftmost kron factor is L-1
    M = np.kron(M, V if l == 2 else I2)
chk("rot_link == explicit kron on link 2", float(np.max(np.abs(r - M @ psi))), 0.0, 1e-12)
psi8 = haar_physical(theta(8), 88)
B = np.linalg.qr(np.random.default_rng(0).normal(size=(4, 4)))[0].astype(complex)
r2 = rot_block(psi8, 8, [0, 1], B)
chk("rot_block preserves norm", float(np.linalg.norm(r2)), 1.0)
chk("rot_block is an involution with its inverse",
    float(np.max(np.abs(rot_block(r2, 8, [0, 1], B.conj().T) - psi8))), 0.0, 1e-12)

# ---- 6  sbs_defects reduces to sbs_defects_multi for a one-qubit system
print("\n[6] sbs_defects (2-dim system) == sbs_defects_multi (general) on the same input")
L = 9; car = theta(L)
for nm, p in [("GHZ", sym_basis_state(car, 0)), ("Haar 777", haar_physical(car, 777))]:
    a = sbs_defects(p, L, [0], [[1, 2], [3, 4], [5, 6]])
    b = sbs_defects_multi(p, L, [0], [[1, 2], [3, 4], [5, 6]])
    for k in ("D_coh", "D_prod", "D_orth"):
        chk(f"{nm} {k}", a[k], b[k], 1e-12)

# ---- 7  the two R_delta scorers agree where they must
print("\n[7] rdelta_count and rdelta_packing agree when the chosen partition IS the best packing")
car = theta(8); L = 8; ENV = list(range(1, L))
psi = sym_basis_state(car, 0)
HS = vn_entropy(reduce_links(psi, L, [0]))
sing = [mi_ext(psi, L, [0], [e]) for e in ENV]
qual = []
for m in range(1, len(ENV) + 1):
    for F in itertools.combinations(ENV, m):
        if mi_ext(psi, L, [0], list(F)) >= 0.9 * HS - 1e-9:
            qual.append(sum(1 << ENV.index(x) for x in F))
chkint("GHZ on theta_8: count over single links vs exact packing",
       rdelta_count(sing, HS, 0.1), rdelta_packing(qual, len(ENV)))

hr(f"r0 SELF-TEST  {'PASSED' if fails == 0 else 'FAILED with ' + str(fails) + ' mismatches'}")
