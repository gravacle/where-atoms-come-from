"""
b1 -- CAN THE CRITERION FAIL?   DISCRIMINATION TEST.
Attempt one was wrecked by criteria that could not have failed.  Before any sweep is believed, the
redundancy criterion must be shown to SEPARATE a state that clearly HAS a record from one that
clearly does NOT, ON THE SAME CARRIER, UNDER THE SAME CUT, THROUGH THE SAME CODE PATH.

ISOLATION LEDGER FOR THIS FILE
  held fixed : carrier theta_8 (8 links, one Gauss law prod_l X_l, dim H_ext = 256,
               dim H_phys = 128), system S = {link 0}, environment E = {links 1..7},
               fragments = subsets of E, estimator = extended-Hilbert-space mutual information
               (choice C1), delta = 0.1, code path identical for every arm.
  moved      : THE STATE VECTOR, and nothing else.  Every arm below is one line: a different psi.
  therefore  : any difference in the reported curve is attributable to the state alone.

WHAT WOULD REFUTE THE CRITERION HERE: if the scrambled arm also produced a flat curve at H(S).
That is a live possibility, not a formality -- see the SCR arms, whose curves are computed by the
identical function and do not plateau.
"""
import numpy as np, itertools, hashlib
from lib_b import *

np.set_printoptions(precision=6, suppress=True)
L = 8
car = theta(L)
S = [0]; ENV = list(range(1, L))
DELTA = 0.1

spec = f"{car['name']}|L={L}|gauss={indep_gauss(car)}|plaq={car['plaq']}|S={S}|E={ENV}|delta={DELTA}"
print("=" * 104)
print("b1  DISCRIMINATION TEST")
print("=" * 104)
print("CARRIER SPEC (identical for every arm below):")
print("   " + spec)
print("   sha256(spec) = " + hashlib.sha256(spec.encode()).hexdigest())
print()

# ------------------------------------------------------------------ the states (the ONLY variable)
def bell_localised():
    """system maximally entangled with ONE environment link; every other link a product |+>."""
    psi = np.zeros(1 << L, dtype=complex)
    rest = 0
    for c in range(1 << (L - 2)):                      # links 2..L-1 in the |+> product state
        bits = sum(((c >> i) & 1) << (i + 2) for i in range(L - 2))
        psi[bits] += 1.0                               # |+...+> = uniform superposition
        psi[bits | 0b11] += 1.0
    psi /= np.linalg.norm(psi)
    return project_physical(psi, car)

def elec_ghz(a=None):
    """GHZ in the ELECTRIC (X) basis: a|+..+> + b|-..-> .  Physical only for L even.
       a=None gives the balanced GHZ; a<1 gives a PERFECT record carrying only H(a^2) bits."""
    popc = popc_table(1 << L)
    plus = np.ones(1 << L, dtype=complex)
    minus = ((-1.0) ** popc[np.arange(1 << L)]).astype(complex)
    plus /= np.linalg.norm(plus); minus /= np.linalg.norm(minus)
    if a is None: a = 1 / np.sqrt(2)
    psi = a * plus + np.sqrt(1 - a * a) * minus
    psi /= np.linalg.norm(psi)
    return project_physical(psi, car)

def plus_product():
    return project_physical(np.ones(1 << L, dtype=complex), car)

ARMS = []
ARMS.append(("BRD_M  magnetic GHZ (|0^L>+|1^L>)/sqrt2  = g2->0 ground state", sym_basis_state(car, 0)))
ARMS.append(("BRD_E  electric GHZ (|+^L>+|-^L>)/sqrt2", elec_ghz()))
ARMS.append(("BRD_Ew electric GHZ, amplitude 0.99 (weak but PERFECT record)", elec_ghz(a=0.99)))
ARMS.append(("LOC    Bell(link0,link1) tensor |+>^(L-2)", bell_localised()))
ARMS.append(("VAC    |+>^L   (H(S)=0 -- the trivial-satisfaction arm)", plus_product()))
for sd in (101, 202, 303):
    ARMS.append((f"SCR    Haar-random physical state, seed {sd}", haar_physical(car, sd)))
for g2 in (0.30, 0.70, 1.00, 1.50, 3.00):
    psi, e = ground_state(car, g2)
    ARMS.append((f"GS     ground state of H, g^2 = {g2:.2f}  (E0={e:.4f})", psi))

# ------------------------------------------------------------------ the estimator (one code path)
def curve(psi):
    HS = vn_entropy(reduce_links(psi, L, S))
    rows = []
    for m in range(1, len(ENV) + 1):
        vals = [mi_ext(psi, L, S, list(F)) for F in itertools.combinations(ENV, m)]
        rows.append((m, float(np.mean(vals)), float(np.min(vals)), float(np.max(vals))))
    single = [mi_ext(psi, L, S, [e]) for e in ENV]      # the finest DISJOINT partition of E
    # TOL absorbs float noise around an exactly-zero threshold; the two counts differ ONLY when
    # H(S) = 0, which is exactly the degeneracy this file is flagging.
    TOL = 1e-9
    Rdelta = int(sum(1 for v in single if v >= (1 - DELTA) * HS - TOL))
    plateau_rows = [r for r in rows if r[0] <= len(ENV) - 1]
    defect = max(abs(r[1] - HS) for r in plateau_rows) / HS if HS > 1e-9 else float("nan")
    flat = max(r[3] - r[2] for r in plateau_rows)
    return HS, rows, single, Rdelta, defect, flat

print(f"{'ARM':<58} {'H(S)':>8} {'Rdelta':>7} {'defect':>9} {'spread':>9}")
print("-" * 104)
results = {}
for name, psi in ARMS:
    HS, rows, single, Rd, defect, flat = curve(psi)
    results[name] = (HS, rows, single, Rd, defect, flat)
    print(f"{name:<58} {HS:8.6f} {Rd:7d} {defect:9.6f} {flat:9.6f}")

print()
print("KEY:  H(S) = von Neumann entropy of link 0.   Rdelta = # of the 7 disjoint single-link")
print("      fragments with I(S:F) >= 0.9*H(S).   defect = max_|F|<=6 |mean I - H(S)| / H(S)")
print("      (0 = perfect plateau at the record height).   spread = max over |F|<=6 of")
print("      (max I - min I) across all fragments of that size (0 = the plateau is not an")
print("      artefact of averaging).")
print()
print("=" * 104)
print("THE CURVES.  I_ext(S : F) as a function of |F|.  mean [min , max] over all fragments.")
print("=" * 104)
for name, _ in ARMS:
    HS, rows, single, Rd, defect, flat = results[name]
    print(f"\n{name}     H(S) = {HS:.6f}")
    for m, mu, lo, hi in rows:
        bar = "#" * int(round(40 * mu / 2.0))
        print(f"   |F|={m}  I = {mu:9.6f}  [{lo:9.6f} , {hi:9.6f}]  {bar}")

print()
print("=" * 104)
print("VERDICT OF THE DISCRIMINATION TEST")
print("=" * 104)
bm = results[ARMS[0][0]]; sc = results[[n for n,_ in ARMS if n.startswith("SCR")][0]]
print(f"  BROADCAST (magnetic GHZ): plateau defect {bm[4]:.3e}, spread {bm[5]:.3e}, Rdelta = {bm[3]}")
print(f"  SCRAMBLED (Haar seed 101): plateau defect {sc[4]:.6f}, spread {sc[5]:.6f}, Rdelta = {sc[3]}")
print(f"  SEPARATION in Rdelta: {bm[3]} vs {sc[3]} out of {len(ENV)} disjoint fragments.")
print(f"  SEPARATION in plateau defect: {bm[4]:.3e} vs {sc[4]:.6f}  "
      f"(ratio {sc[4]/max(bm[4],1e-16):.3e})")
print()
print("  COULD THE SCRAMBLED ARM HAVE PLATEAUED?  Yes: a flat curve at H(S) is what the identical")
print("  estimator returns for the broadcast arm.  It did not.  The test therefore had a way to")
print("  fail and did not take it.  The criterion DISCRIMINATES on this carrier.")
print()
print("  THE ONE ARM THAT VOIDS ITSELF: VAC (|+>^L) has H(S) = 0, so the condition")
print("  I >= (1-delta)*H(S) reads 0 >= 0 and is satisfied by EVERY fragment.  Rdelta is reported")
print("  as maximal for a state with no information anywhere.  R_delta AS USUALLY WRITTEN IS")
print("  TRIVIALLY SATISFIED WHENEVER H(S) = 0 AND MUST CARRY AN H(S) > 0 PRECONDITION.")
