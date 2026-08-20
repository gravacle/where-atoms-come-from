"""ADVERSARIAL CHECK of T-31: is the clustering width the amended definition's width, or an ad hoc
tolerance? The lane claims 'clustering within clause (ii')'s stated durability width recovers the
count' (C-76). Clause (ii')'s width is delta = hbar/t_m = 3.341e-43 J (its own table prints it).
The width actually used is wdt = 1e-2 in the carrier's DIMENSIONLESS units, while the table labels
the splittings '(J)' and prints splitting/delta ~ 1e36..1e41 -- i.e. by the lane's own column the
splittings are NEVER within delta.

Test 1: rerun the clustered count with the lane's own delta (taking the printed '(J)' label at
        face value). If k = 0, the demonstrated fix is NOT the durability width.
Test 2: give the carrier a PHYSICAL energy scale from this program's own census (T-29 grain:
        splitting 0.007..7 kT at 300 K; census range 0.2 eV..MJ). Convert delta into carrier units
        and cluster. If k = 0 for every physical scale, no real record's count is recovered by the
        amended definition's width.
Test 3: how small would the carrier's unit energy have to be for delta to cover a 1e-6 splitting,
        and what would E_b/kT then be? (clause (v') requires E_b/kT >> 1)."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'model'))
import grounded as G

I2 = np.eye(2); X = np.array([[0,1],[1,0]], dtype=complex); Z = np.array([[1,0],[0,-1]], dtype=complex)
def word(n, s):
    M = np.array([[1]], dtype=complex)
    for c in s: M = np.kron(M, {'I':I2,'X':X,'Z':Z}[c])
    return M
def v2(m):
    k = 0
    while m % 2 == 0 and m > 0: m //= 2; k += 1
    return k
def count_law(H, width):
    w = np.sort(np.linalg.eigvalsh(H)); shells = []
    for x in w:
        if shells and abs(x - shells[-1][0]) <= width: shells[-1][1] += 1
        else: shells.append([float(x), 1])
    return min(v2(m) for _, m in shells), [m for _, m in shells]

n = 4
Hsym = -(word(n,'X'*n) + word(n,'Z'*n))
def generic(eps):
    return eps*sum((1.0+0.6180339887*j)*word(n, ''.join('Z' if i==j else 'I' for i in range(n)))
                   for j in range(n))
t_m = 10.0*3.156e7
delta = G.HBAR/t_m
print(f"delta = hbar/t_m = {delta:.3e} J   (the lane's own printed durability width)")
print()
print("TEST 1 — cluster with delta itself, taking the lane's '(J)' column label at face value:")
print(f"{'eps':>10}{'width used':>14}{'multiplicities':>40}{'k':>4}")
for eps in (1e-6, 1e-3, 1e-1):
    k, m = count_law(Hsym + generic(eps), delta)
    print(f"{eps:>10.0e}{delta:>14.3e}{str(m):>40}{k:>4}")
k_lane, m_lane = count_law(Hsym + generic(1e-6), 1e-2)
print(f"the lane's own result at eps=1e-6 with its ACTUAL width 1e-2: k={k_lane}, {m_lane}")
print()
print("TEST 2 — physical energy scales for the carrier unit, census range; delta in carrier units:")
EV = 1.602176634e-19
print(f"{'carrier unit':>16}{'delta (carrier units)':>24}{'covers eps=1e-6?':>18}{'k at eps=1e-6':>15}")
for label, unitJ in (("0.2 eV (census min)", 0.2*EV), ("1 kT @300K", G.KB*300),
                     ("0.007 kT (T-29 min dE)", 0.007*G.KB*300), ("1 J", 1.0)):
    d_units = delta/unitJ
    k, m = count_law(Hsym + generic(1e-6), d_units)
    print(f"{label:>22}{d_units:>18.2e}{str(d_units>=1e-6):>18}{k:>15}")
print()
print("TEST 3 — for delta to cover a splitting of eps (carrier units), the carrier unit must be")
print("         <= delta/eps in J; the barrier (2 units on this carrier) is then:")
for eps in (1e-6, 1e-3):
    unit_max = delta/eps
    Eb = 2*unit_max
    print(f"  eps={eps:.0e}: unit <= {unit_max:.2e} J -> E_b <= {Eb:.2e} J = {Eb/(G.KB*300):.2e} kT at 300 K")
print("  clause (v') requires E_b/kT >> 1: impossible by ~40 orders of magnitude.")
print()
print("CROSS-CHECK against the lane's own table: it prints splitting/delta = 2.99e+36..2.99e+41 and,")
print("two lines later, 'within delta? True' for the same splittings. Both cannot be true.")
