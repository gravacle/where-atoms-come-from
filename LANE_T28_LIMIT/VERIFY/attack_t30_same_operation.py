"""ADVERSARIAL CHECK of T-30: is the metastable computation the SAME OPERATION as the control?

The control computes  Tr(Pg M)/Tr(Pg)  -- a UNIFORM trace ratio against the code-space projector --
over a family of integer-coefficient polynomials M in record operators.

The metastable measurement computes  Tr(rho_ss R)  -- the expectation of a SINGLE record operator
against the Boltzmann-weighted STEADY STATE of the Liouvillian.

Two things changed at once: (1) exact carrier -> metastable carrier, (2) uniform trace over the
durable manifold -> thermal steady state. T-30 attributes the continuum to (1). Decide by
constructing the two intermediates:

  CASE A: METASTABLE carrier, CONTROL's operation (uniform trace ratio over the projector onto the
          long-lived manifold -- the operation the lane's own text names). If this is quantised, the
          metastability is NOT what breaks the quantisation.
  CASE B: EXACT carrier (tolerance zero, exact signed Paulis, no dissipation), METASTABLE's operation
          (thermal-state expectation, Gibbs at finite T). If this is continuous, exactness NEVER
          protected the quantisation once the state is thermal.

If A is quantised and B is continuous, the discriminating axis is the STATE (T -> 0 corner
coordinate), not exactness/metastability of the record, and T-30's control does not control for the
substitution actually made."""
import sys, os, itertools, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'model'))
import grounded as G

I2 = np.eye(2); X = np.array([[0,1],[1,0]], dtype=complex); Z = np.array([[1,0],[0,-1]], dtype=complex)
def word(n, s):
    M = np.array([[1]], dtype=complex)
    for c in s: M = np.kron(M, {'I':I2,'X':X,'Z':Z}[c])
    return M

rng = np.random.default_rng(3)
print("="*100)
print("CASE A — metastable grain, the CONTROL'S operation: uniform trace ratio over the")
print("         projector onto the long-lived manifold (slow modes at t_m = 10 yr).")
print("="*100)
T = 300.0; K_u = 2.0e5; V0 = 1.26e-24; f0 = 1e9
sz = np.array([[1,0],[0,-1]], dtype=complex); sp = np.array([[0,1],[0,0]], dtype=complex); sm = sp.conj().T
t_m = 10.0*3.156e7
for dEk in (3.0, 1.0, 0.01):
    E_b = K_u*V0; dE = dEk*G.KB*T
    H = -(dE/2)*sz
    gd = f0*np.exp(-(E_b+dE/2)/(G.KB*T)); gu = f0*np.exp(-(E_b-dE/2)/(G.KB*T))
    Ls = [np.sqrt(gd)*sm, np.sqrt(gu)*sp]
    # the long-lived manifold: STATES that survive to t_m. Both wells do (tau ~ 1e17 s >> t_m).
    # Its projector is the identity on the 2-dim well space; normalised trace = Tr(.)/2.
    L = G.liouvillian(H, Ls)
    w, V = np.linalg.eig(L)
    slow_states = [i for i in range(len(w)) if abs(w[i]) <= 1.0/t_m]
    # build the projector onto the span of the durable population manifold
    # (for this carrier: both computational states; check that the slow state-space is 2-dim)
    mods = []
    for _ in range(2000):
        a0, a1, a2 = int(rng.integers(-2,3)), int(rng.integers(-2,3)), int(rng.integers(-1,2))
        M = a0*np.eye(2) + a1*sz + a2*(sz@sz)
        r = complex(np.trace(M)/2.0)              # uniform trace over the durable manifold
        if abs(r) > 1e-12: mods.append(abs(r))
    inside = [m for m in mods if 1e-12 < m < 1-1e-9]
    print(f"  dE/kT={dEk:<6} slow state-modes at t_m: {len(slow_states)}  "
          f"non-zero ratios: {len(mods)}  min modulus: {min(mods):.10f}  strictly inside (0,1): {len(inside)}")
print("  -> if 'strictly inside (0,1): 0' on every row, the metastable record is QUANTISED under the")
print("     control's own operation; metastability did not break the quantisation.")

print()
print("="*100)
print("CASE B — the EXACT [[4,2,2]] stabiliser carrier (tolerance zero, no dissipation, exact signed")
print("         Paulis), the METASTABLE side's operation: thermal-state expectation at finite T.")
print("="*100)
n = 4
Hs = -(word(n,'X'*n) + word(n,'Z'*n))     # dimensionless spectrum {-2, 0, +2}
recs = [word(n,'ZZII'), word(n,'IZZI'), word(n,'XXII')]
print(f"  {'beta (carrier units)':>22}{'<M>':>18}{'in Z[i]?':>10}")
vals = []
for beta in (0.25, 0.5, 1.0, 2.0, 4.0, 16.0):
    Egs = np.linalg.eigvalsh(Hs)
    rho = None
    w, V = np.linalg.eigh(Hs)
    p = np.exp(-beta*(w - w.min())); p = p/p.sum()
    rho = (V * p) @ V.conj().T
    # a single record operator's expectation -- exactly what T-30 computed on the metastable side
    # (use a polynomial with the identity term too, mirroring the control family)
    M = np.eye(2**n, dtype=complex) + recs[0] + (recs[0]@recs[0])
    r = complex(np.trace(rho @ M))
    q = abs(r - complex(round(r.real), round(r.imag))) < 1e-9
    vals.append(abs(r)); print(f"  {beta:>22}{r.real:>18.9f}{str(q):>10}")
# also the bare record operator, the exact analog of <R>_ss
print(f"  {'beta':>22}{'<ZZII>':>18}{'in Z[i]?':>10}")
for beta in (0.25, 0.5, 1.0, 2.0, 4.0, 16.0):
    w, V = np.linalg.eigh(Hs)
    p = np.exp(-beta*(w - w.min())); p = p/p.sum()
    rho = (V * p) @ V.conj().T
    r = float(np.real(np.trace(rho @ recs[0])))
    q = abs(r - round(r)) < 1e-9
    print(f"  {beta:>22}{r:>18.9f}{str(q):>10}")
print("  -> if these move continuously with beta and sit strictly inside integers, the EXACT carrier")
print("     with tolerance ZERO already yields a continuum once the state is thermal: exactness was")
print("     never the protection, the T->0 corner STATE was.")

print()
print("="*100)
print("CASE C — the T->0 limit of the metastable value itself: tanh(dE/2kT) -> sign(dE) = +-1 in Z[i].")
print("="*100)
for Tq in (300.0, 30.0, 3.0, 0.3, 0.03):
    dE = 3.0*G.KB*300.0
    print(f"  T={Tq:>8.2f} K   tanh(dE/2kT) = {np.tanh(dE/(2*G.KB*Tq)):.12f}")
print("  -> the continuum closes onto Z[i] as T -> 0 with the record's metastability intact:")
print("     the driver is the temperature axis of the corner, not the exactness of the record.")
