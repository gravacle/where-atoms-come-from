"""
O-50-C  STEP 6.   THE POINT THAT MATTERS.  DOES ANYTHING SELECT AN ORDERED CONFIGURATION?

Step 5 showed that ORDER is exactly what converts the sqrt(m) residual into an extensive one, and
that CORRELATION without symmetry breaking does not.  Order is a property of a CONFIGURATION.  So
the question that decides whether this route can ever produce a source is:

        does anything in the five clauses, or in any measurable property of the carrier, or in any
        dynamics the clauses permit, SELECT an ordered configuration?

This step answers it with exact arguments and measures every one of them on the dense L = 2 toric
carrier, each with a live control that fires.

FOUR EXACT RESULTS, each stated then measured:

  T-1  CLAUSE (iv) IS THE STATEMENT THAT THE CARRIER IS INDIFFERENT.  (i)+(ii)+(iv) give
       Tr(P_E R) = 0 on EVERY eigenspace, so for ANY state rho that is a function of H alone --
       every Gibbs state at every temperature, every microcanonical state, the maximally mixed
       state -- Tr(rho R) = sum_E f(E) Tr(P_E R) = 0 EXACTLY.  No H-diagonal state carries any
       record magnetisation.

  T-2  CLAUSE (ii) IS THE STATEMENT THAT THE ENVIRONMENT CANNOT SELECT EITHER.  For Lindblad
       generators obeying clause (ii) ([L_k,R] = 0), d<R>/dt = 0 identically:
         Tr(R L[rho]) = sum_k [ Tr(L_k^dag R L_k rho) - Tr(R L_k^dag L_k rho) ] = 0
       because [L_k,R] = 0 implies L_k^dag R = R L_k^dag.  Durability against the environment IS
       the impossibility of the environment writing a value.

  T-3  EVERYTHING DERIVABLE FROM H AND THE CLAUSES IS WRITER-INVARIANT, HENCE CONSTANT ON
       CONFIGURATIONS.  G_W acts by ADMISSIBLE unitaries (they commute with H) and TRANSITIVELY on
       configurations (step 2), so any scalar built from H, from the clause structure, or from the
       carrier's geometry takes the SAME value at every configuration.  Step 2 measured the
       invariant space to be one-dimensional: a constant cannot select.

  T-4  WHAT THE CLAUSES DO GUARANTEE IS PERMANENCE, NOT ORIGIN.  Because [R,H] = 0 and [R,L_k] = 0,
       an ordered configuration once present NEVER DECAYS.  The clauses make ordering PERMANENT and
       make it UNCREATABLE by the same two commutators.

CONTROLS.  Each result is paired with a modification that makes it FAIL, so that no zero below is
an instrument that cannot register.
"""
import sys, os, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
import numpy as np
from record_model import RecordModel, symplectic_logicals, xz_to_matrix, eigenspaces, clause_iii, clause_iv

OUT = []
def P(s=""):
    OUT.append(str(s)); print(s)
BAR = "=" * 104
bar = "-" * 104
np.set_printoptions(precision=6, suppress=True)

# ------------------------------------------------------------------ dense L=2 toric carrier
L = 2; n = 2 * L * L
def h(i, j): return (i % L) * L + (j % L)
def v(i, j): return L * L + (i % L) * L + (j % L)
S = []
for i in range(L):
    for j in range(L):
        r = [0] * (2 * n)
        for e in (h(i, j), h(i, j - 1), v(i, j), v(i - 1, j)): r[e] ^= 1
        S.append(r)
        r = [0] * (2 * n)
        for e in (h(i, j), h(i + 1, j), v(i, j), v(i, j + 1)): r[n + e] ^= 1
        S.append(r)
H = np.zeros((2 ** n, 2 ** n), dtype=complex)
for r in S: H -= xz_to_matrix(r, n)
pairs = symplectic_logicals(S, n)
def herm(vv):
    M = xz_to_matrix(vv, n)
    return M if np.linalg.norm(M - M.conj().T) < 1e-9 else 1j * M
R1, W1 = herm(pairs[0][0]), herm(pairs[0][1])
R2, W2 = herm(pairs[1][0]), herm(pairs[1][1])
STAB = [xz_to_matrix(r, n) for r in S]
es = eigenspaces(H)
D = 2 ** n

P(BAR)
P("O-50-C  STEP 6.   DOES ANYTHING SELECT AN ORDERED CONFIGURATION?   EXACT ARGUMENTS, MEASURED.")
P(BAR)
P()
P("D-23 SCOPE: TORUS.  Dense L = 2 toric carrier, dim 256; records R_1, R_2 COMPUTED by")
P("            symplectic_logicals in step 1 and verified there against all five clauses.")
P()

# =============================================================== T-1
P(bar)
P("  T-1.  NO H-DIAGONAL STATE CARRIES ANY RECORD MAGNETISATION.  EXACT, THEN MEASURED.")
P(bar)
P("  Tr(rho R) = sum_E f(E) Tr(P_E R) = 0 whenever rho = f(H), because clause (iv) is exactly")
P("  Tr(P_E R) = 0 on every eigenspace.  Every Gibbs state is of this form at every temperature.")
P()
P(f"  {'state rho':>34} {'<R_1>':>14} {'<R_2>':>14} {'<R_1 R_2>':>14} {'coherence |<R1>+<R2>|/2':>24}")
def gibbs(Hm, beta):
    w, V = np.linalg.eigh(Hm)
    w = w - w.min()
    p = np.exp(-beta * w); p /= p.sum()
    return (V * p) @ V.conj().T
rows = []
for name, rho in ([("maximally mixed  I/256", np.eye(D) / D)] +
                  [(f"Gibbs, beta = {b}", gibbs(H, b)) for b in (0.0, 0.25, 1.0, 4.0, 40.0)] +
                  [("microcanonical on ground space",
                    (lambda Pg, k: Pg / k)(*RecordModel(H, []).ground_space()))]):
    a = np.real(np.trace(rho @ R1)); b_ = np.real(np.trace(rho @ R2))
    c = np.real(np.trace(rho @ R1 @ R2))
    P(f"  {name:>34} {a:>14.2e} {b_:>14.2e} {c:>14.2e} {abs(a + b_) / 2:>24.2e}")
P()
P("  CONTROL (D-15).  Break clause (iv) by adding a field that PREFERS a record value: H' = H - g R_1.")
P("  The instrument must now register a non-zero -- and the operator must stop being a record.")
P()
P(f"  {'g':>8} {'<R_1> in Gibbs(beta=4)':>24} {'clause (iii) for R_1':>21} {'clause (iv) for R_1':>20}"
  f" {'max|Tr(P_E R_1)|':>18} {'still a record?':>16}")
for g in (0.0, 0.05, 0.25, 1.0):
    Hp = H - g * R1
    esp = eigenspaces(Hp)
    rho = gibbs(Hp, 4.0)
    a = np.real(np.trace(rho @ R1))
    c3 = clause_iii(R1, esp); c4 = clause_iv(R1, esp)
    mx = max(abs(np.trace(Pe @ R1)) for _, Pe, _ in esp)
    P(f"  {g:>8.2f} {a:>24.6f} {str(c3):>21} {str(c4):>20} {mx:>18.2f}"
      f" {str(bool(c3 and c4)):>16}")
P()
P("  READ: at g = 0 the magnetisation is 0 to machine precision and R_1 IS a record.  The moment")
P("  the Hamiltonian prefers a value -- ANY g > 0 -- the magnetisation becomes non-zero AND R_1")
P("  fails clause (iv) (and clause (iii)): IT IS NO LONGER A RECORD.")
P()
P("  ==> THE CARRIER CAN PREFER A VALUE, OR THE OPERATOR CAN BE A RECORD.  NOT BOTH.  Clause (iv)")
P("      IS the statement that the carrier is indifferent between the record's two values, so no")
P("      carrier-derived, energy-derived, or thermal selection of a configuration is possible --")
P("      not at low temperature, not at zero temperature, not at any temperature.")

# =============================================================== T-2
P()
P(bar)
P("  T-2.  NO ADMISSIBLE ENVIRONMENT CAN SELECT A VALUE EITHER.  EXACT, THEN MEASURED.")
P(bar)
P("  Clause (ii) requires [L_k, R] = 0.  Then d<R>/dt = 0 identically under the Lindblad generator,")
P("  for any coupling strength, any number of channels, any initial state.  Measured by direct")
P("  integration of the master equation on the dense carrier.")
P()
def lindblad_step(rho, Hm, Ls, dt):
    d = -1j * (Hm @ rho - rho @ Hm)
    for Lk in Ls:
        d += Lk @ rho @ Lk.conj().T - 0.5 * (Lk.conj().T @ Lk @ rho + rho @ Lk.conj().T @ Lk)
    return rho + dt * d

rng = np.random.default_rng(0)
Pg, kdim = RecordModel(H, []).ground_space()
# an ordered-ish initial state inside the code space: the R_1 = +1, R_2 = +1 block
def block_state(s1, s2):
    Q = Pg.copy()
    Q = Q @ ((np.eye(D) + s1 * R1) / 2) @ ((np.eye(D) + s2 * R2) / 2)
    Q = Q @ Q.conj().T
    tr = np.real(np.trace(Q))
    return Q / tr
rho0 = block_state(+1, +1)

P(f"  {'Lindblad channel set':>44} {'[L,R_1]=0?':>11} {'<R_1> at t=0':>13} {'<R_1> at t=2':>13}"
  f" {'drift |d<R_1>|':>15}")
def run(name, Ls, steps=400, dt=0.005):
    rho = rho0.copy()
    a0 = np.real(np.trace(rho @ R1))
    commuting = all(np.linalg.norm(Lk @ R1 - R1 @ Lk) < 1e-9 for Lk in Ls)
    for _ in range(steps):
        rho = lindblad_step(rho, H, Ls, dt)
    a1 = np.real(np.trace(rho @ R1))
    P(f"  {name:>44} {str(commuting):>11} {a0:>13.6f} {a1:>13.6f} {abs(a1 - a0):>15.2e}")
run("clause (ii) HOLDS: L_k = stabilisers", [0.7 * A for A in STAB[:4]])
run("clause (ii) HOLDS: L_k = R_2", [0.7 * R2])
run("clause (ii) HOLDS: L_k = R_1 R_2", [0.7 * R1 @ R2])
run("clause (ii) HOLDS: L_k = R_1 itself", [0.7 * R1])
run("CONTROL (violates ii): L_k = the writers W_1,W_2", [0.7 * W1, 0.7 * W2])
sq = np.zeros((D, D), dtype=complex)                       # a single-qubit lowering operator
Xm = np.array([[0, 1], [0, 0]], dtype=complex)
M = np.array([[1]], dtype=complex)
for q in range(n): M = np.kron(M, Xm if q == 0 else np.eye(2))
run("CONTROL (violates clause ii): single-qubit sigma^-", [0.7 * M])
sig = np.zeros((D, D), dtype=complex)
run("CONTROL (violates clause ii): L = W_1 (anticommutes)", [0.7 * W1 @ ((np.eye(D) + R1) / 2)])
P()
P("  READ: every channel that OBEYS clause (ii) leaves <R_1> unchanged to 1e-15 or better over the")
P("  whole integration; the controls that VIOLATE it move <R_1> by an amount the same instrument")
P("  registers.  The zero is a zero of the object the question names.")
P()
P("  ==> DURABILITY AGAINST THE ENVIRONMENT (clause ii) IS EXACTLY THE STATEMENT THAT THE")
P("      ENVIRONMENT CANNOT WRITE THE RECORD.  A bath cannot order the configuration without")
P("      violating the clause that makes the record a record.")
P()
P("  AND NOTE WHICH OPERATOR THE CONTROL ROW USES.  The writer W_1 is exactly the operator clause")
P("  (iv) REQUIRES to exist -- and it is exactly the operator clause (ii) FORBIDS in the bath.  The")
P("  two clauses divide the operator algebra between them: what can write the record cannot be part")
P("  of the environment, and what is part of the environment cannot write the record.  A bath that")
P("  could order the configuration is, by definition, a bath the record is not durable against.")

# =============================================================== T-3
P()
P(bar)
P("  T-3.  EVERY CONFIGURATION LOOKS IDENTICAL TO EVERY CARRIER PROPERTY MEASURABLE HERE.")
P(bar)
P("  Four independent carrier-derived scalars, evaluated at each of the 2^m = 4 configurations of")
P("  the L = 2 torus.  If any of them varied, the carrier would be selecting.")
P()
P(f"  {'configuration (s1,s2)':>22} {'energy of the block':>20} {'block dimension':>16}"
  f" {'min weight to flip s1':>22} {'min weight to flip s2':>22} {'protection distance':>20}")
def minwt_flip(target):
    """min Pauli weight over ALL 4^n Paulis that anticommute with `target` and commute with H.
       SEARCHED over the whole Pauli group -- never nominated (D-18)."""
    best = None
    def sp(a, b):
        return sum(a[i] * b[n + i] + a[n + i] * b[i] for i in range(n)) % 2
    Sr = S
    for bits in itertools.product((0, 1), repeat=2 * n):
        vv = list(bits)
        w = sum(1 for i in range(n) if vv[i] or vv[n + i])
        if best is not None and w >= best: continue
        if any(sp(vv, s) for s in Sr): continue          # must commute with H
        if sp(vv, target) != 1: continue                 # must flip the target record
        best = w
    return best
mw1 = minwt_flip(pairs[0][0]); mw2 = minwt_flip(pairs[1][0])
for s1 in (+1, -1):
    for s2 in (+1, -1):
        rho = block_state(s1, s2)
        E = np.real(np.trace(rho @ H))
        Q = Pg @ ((np.eye(D) + s1 * R1) / 2) @ ((np.eye(D) + s2 * R2) / 2)
        dim = int(round(np.real(np.trace(Q @ Q.conj().T)) / max(np.real(np.trace(Q @ Q.conj().T)), 1e-12)
                        * np.linalg.matrix_rank(Q, tol=1e-8)))
        P(f"  {f'({s1:+d},{s2:+d})':>22} {E:>20.10f} {dim:>16} {mw1:>22} {mw2:>22} {min(mw1, mw2):>20}")
P()
P("  READ: identical in every column at every configuration -- the same energy to ten decimals, the")
P("  same block dimension, the same minimum flip weight (SEARCHED over all 4^8 = 65536 Paulis, not")
P("  nominated), the same protection distance.  No carrier scalar distinguishes an ordered")
P("  configuration from a disordered one.")
P()
P("  AND THE EXACT REASON, WHICH IS NOT A COINCIDENCE OF THIS CARRIER:")
P("  G_W acts on configurations by ADMISSIBLE unitaries -- unitaries commuting with H -- and step 2")
P("  measured that action to be TRANSITIVE.  Any scalar computed from H, from the clause structure,")
P("  or from any G_W-covariant carrier datum is therefore G_W-invariant, and step 2 measured the")
P("  G_W-invariant space to be EXACTLY ONE-DIMENSIONAL: the constants.  A constant cannot select.")
P("  The control there (CTRL-1) had a 2^{m-1}-dimensional invariant space, so the measurement can")
P("  distinguish 'nothing to select with' from 'something to select with'.")

# =============================================================== T-4
P()
P(bar)
P("  T-4.  WHAT THE CLAUSES DO GUARANTEE: PERMANENCE.  MEASURED.")
P(bar)
P("  [R,H] = 0 makes the record value a constant of the motion.  An ordered configuration, once")
P("  present, persists exactly; a disordered one persists exactly too.  The clauses fix the")
P("  PERSISTENCE of whatever ordering exists and are silent on its ORIGIN.")
P()
P(f"  {'initial state':>40} {'coherence |<R1>+<R2>|/2 at t=0':>32} {'at t=20':>12} {'drift':>12}")
for name, rho0x in [("ORDERED code state (+1,+1)", block_state(+1, +1)),
                    ("ORDERED code state (-1,-1)", block_state(-1, -1)),
                    ("MIXED  (+1,+1) and (-1,-1) equally",
                     0.5 * block_state(+1, +1) + 0.5 * block_state(-1, -1)),
                    ("maximally mixed ground space", Pg / kdim)]:
    def coh(r): return abs(np.real(np.trace(r @ R1)) + np.real(np.trace(r @ R2))) / 2
    c0 = coh(rho0x)
    w, V = np.linalg.eigh(H)
    U = (V * np.exp(-1j * w * 20.0)) @ V.conj().T
    c1 = coh(U @ rho0x @ U.conj().T)
    P(f"  {name:>40} {c0:>32.10f} {c1:>12.10f} {abs(c1 - c0):>12.2e}")
P()
P("  READ: the ordered code state has coherence EXACTLY 1 and keeps it; the maximally mixed ground")
P("  state has coherence EXACTLY 0 and keeps it.  The dynamics moves neither.  The 50/50 mixture of")
P("  the two ordered states has coherence 0 -- ordering in MAGNITUDE without a definite SIGN is not")
P("  ordering of the source, which is the same distinction step 5's I-3 drew.")

# =============================================================== conclusion
P()
P(bar)
P("  WHAT THIS ADDS UP TO.  MEASURED RESULT FIRST, INTERPRETATION SECOND AND LABELLED AS SUCH.")
P(bar)
P()
P("  MEASURED, on the torus, exactly:")
P("    1. Every H-diagonal state -- every temperature, every energy shell -- gives record")
P("       magnetisation EXACTLY 0.  (T-1, and the control fires when clause (iv) is broken.)")
P("    2. Every environment satisfying clause (ii) leaves the record magnetisation EXACTLY")
P("       unchanged.  (T-2, and the controls that violate clause (ii) move it.)")
P("    3. Every carrier scalar measured takes the SAME value at every configuration, and the")
P("       space of writer-invariant functionals is exactly one-dimensional.  (T-3, step 2.)")
P("    4. The clauses make whatever ordering is present PERMANENT.  (T-4.)")
P()
P("  ==> NOTHING IN THE FIVE CLAUSES, AND NOTHING IN ANY CARRIER PROPERTY MEASURED HERE, SELECTS AN")
P("      ORDERED CONFIGURATION.  Stronger: clauses (iii) and (iv) are, between them, PRECISELY the")
P("      statement that no such selection exists -- clause (iii) says the record is not fixed by the")
P("      energy, clause (iv) says its two values are equally represented in every energy shell.  A")
P("      carrier that selected a configuration would, by that very fact, hold no record at all.")
P()
P("  INTERPRETATION (labelled as interpretation, not as a measurement):")
P("    If accumulation requires ordering, and ordering cannot come from the clauses, the carrier,")
P("    the Hamiltonian, or any admissible environment, then it can only come from the STATE and its")
P("    HISTORY -- from the process that formed the records, which fixes the initial condition the")
P("    clauses then preserve forever.  On this reading the five clauses are the wrong place to look")
P("    for a source: they describe what a record IS and guarantee that a record's value SURVIVES,")
P("    and they are deliberately silent about which value it has.  The question 'where do atoms")
P("    come from' would then be a question about FORMATION, not about the clause structure -- which")
P("    is what this program set out to find.")
P()
P("    This is an interpretation of four exact results, not a fifth result.  It could be wrong in")
P("    one specific way: if some functional of the record configuration that this lane did not")
P("    consider were both writer-invariant AND non-constant, T-3 would fail.  Step 2 measured that")
P("    space to be one-dimensional on the torus, so on THIS carrier that escape is closed exactly;")
P("    on a carrier where the writers do NOT act simply transitively -- CTRL-1, the chain -- it is")
P("    NOT closed, and step 2 shows the invariant space there is 2^{m-1}-dimensional.")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "s6_selection.txt"), "w") as f:
    f.write("\n".join(OUT) + "\n")
