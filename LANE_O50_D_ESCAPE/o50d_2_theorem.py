"""O-50 D  PART 2 -- THE THEOREM, STATED EXACTLY, WITH ITS ONE FALSE CLAUSE REMOVED.

The candidate as handed to this lane:
   IF records are independently writable, G_W acts simply transitively on configurations,
   THEN every G_W-invariant functional is CONSTANT,
   AND every non-constant functional is G_W-ODD, hence has mean exactly zero, hence CANCELS.

The middle line is TRUE.  The last line is FALSE as written, and a one-line counterexample
kills it.  What survives is stronger in one way and weaker in another, and the difference is
exactly where the escapes have to live.
"""
import sys, itertools, math
from fractions import Fraction
import numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_D_ESCAPE")
from o50d_common import *

say("=" * 104)
say("O-50 D  PART 2   THE THEOREM, EXACTLY")
say("=" * 104)

# ---------------------------------------------------------------- 1. the counterexample
say("")
say("1. THE LITERAL STATEMENT IS FALSE.  ONE LINE KILLS IT.")
k = 3
cfg = list(itertools.product((1, -1), repeat=k))
f = {s: Fraction(1 + s[0], 2) for s in cfg}                       # the indicator of s_1 = +1
mean = sum(f.values()) / len(cfg)
resp = max(abs(f[s] - f[(-s[0],) + s[1:]]) for s in cfg)
say(f"   f(s) = (1 + s_1)/2 on k={k} records.  values {sorted(set(f.values()))}")
say(f"   non-constant: True     responsive to writing record 1: max|df| = {resp}")
say(f"   mean over configurations = {mean}  -- NOT ZERO.  It does not cancel.")
say("   So 'non-constant => G_W-odd => mean zero' is REFUTED.  A non-constant functional is a")
say("   sum of a CONSTANT and an odd part; only the odd part is forced to vanish.")

# ---------------------------------------------------------------- 2. what is actually true
say("")
say("2. WHAT IS ACTUALLY TRUE -- AND IT NEEDS NO HYPOTHESIS AT ALL.")
say("   LEMMA (orbit averaging).  Let ANY group G act on ANY finite configuration set X, and let")
say("   f : X -> R.  Define (Pf)(x) = (1/|Gx|) sum_{y in Gx} f(y), the orbit average.  Then")
say("       (1) Pf is G-INVARIANT, hence UNCHANGED by every admissible write;")
say("       (2) f - Pf has mean EXACTLY ZERO on every orbit;")
say("       (3) f is responsive to writing  <=>  f != Pf  <=>  f - Pf != 0.")
say("   Therefore: THE PART OF ANY FUNCTIONAL THAT RESPONDS TO WRITING AVERAGES TO EXACTLY ZERO")
say("   OVER THE WRITE ORBIT.  Proof: (1) and (2) are the definition of an average over a group")
say("   orbit; (3) is immediate.  NO simple transitivity, NO independent writability, NO")
say("   abelianness, NO structure on X is used.  EXACT.")
say("")
say("   SIMPLE TRANSITIVITY ADDS EXACTLY ONE THING: the orbit is all of X, so Pf is a GLOBAL")
say("   constant and the invariant space is 1-dimensional (Part 1 section 6 measured dim = 1).")
say("   Without it the invariant space is larger -- one dimension per orbit -- and that is the")
say("   only door the hypothesis opens.  Section 5 walks through it.")

# ---------------------------------------------------------------- 3. Fourier form on Z_2^k
say("")
say("3. THE Z_2^k FORM, WHICH IS WHAT THE TORUS GIVES, AND ITS EXACT CANCELLATION LAW")
say("   f(s) = sum_{T subset [k]} fhat(T) chi_T(s), chi_T(s) = prod_{i in T} s_i.")
say("   A write by generator g multiplies chi_T by (-1)^{|T and supp(g)|}.  chi_empty is the only")
say("   invariant character, so Pf = fhat(empty) = the mean, and every responsive piece is a sum")
say("   of NON-TRIVIAL characters, each of which takes +1 and -1 EQUALLY OFTEN.")
say("")
say(f"   {'k':>3}{'#configs':>10}{'mean of chi_T (T nonempty), exact':>36}{'#sign-definite characters':>28}")
for k in range(1, 13):
    cfg = list(itertools.product((1, -1), repeat=k))
    worst = 0; sd = 0
    for T in itertools.chain.from_iterable(itertools.combinations(range(k), r) for r in range(1, k + 1)):
        vals = [math.prod(s[i] for i in T) for s in cfg]
        worst = max(worst, abs(sum(vals)))
        if min(vals) * max(vals) > 0: sd += 1
    say(f"   {k:>3}{len(cfg):>10}{('max |sum| = ' + str(worst)):>36}{sd:>28}")
say("   EVERY non-trivial character sums to EXACTLY 0 and NONE is sign-definite.  Exact integers,")
say("   not a numerical trend.  A responsive functional therefore FAILS standard (d) identically.")

# ---------------------------------------------------------------- 4. the accumulation law
say("")
say("4. HOW BADLY IT CANCELS AS RECORDS ARE ADDED -- and the control that shows this is not")
say("   about records at all.")
def coherence(weights, rng, samples=200000):
    m = len(weights); w = np.asarray(weights, dtype=float)
    S = rng.integers(0, 2, size=(samples, m)) * 2 - 1
    v = S @ w
    return float(np.abs(v).mean() / np.abs(w).sum())
rng = np.random.default_rng(7)
say(f"   {'m terms':>9}{'coherence |sum|/sum| |':>24}{'sqrt(2/pi)/sqrt(m)':>22}{'ratio':>9}")
for m in (4, 16, 64, 256, 1024, 4096):
    c = coherence(np.ones(m), rng)
    pred = math.sqrt(2 / math.pi) / math.sqrt(m)
    say(f"   {m:>9}{c:>24.6f}{pred:>22.6f}{c / pred:>9.4f}")
say("   Reproduces C-62's closed form sqrt(2/pi) m^(-1/2) out of sample.  D-20: this is a")
say("   PREDICTION checked, not a fit.")
say("")
say("   CONTROL (C-61 discipline).  The same law holds for a sum over the +-1 labels of ANY")
say("   register, record or not:")
Hc, J, h, diag = control_carrier(10, seed=3)
say(f"   zero-record control carrier (C-61): H' = sum J Z Z + sum h Z, n=10, "
    f"{len(set(np.round(diag,9)))} distinct energies of {len(diag)} -> NON-DEGENERATE -> "
    f"ZERO records by P-1")
lab = np.array([[1 - 2 * ((b >> (9 - q)) & 1) for q in range(10)] for b in range(1024)])
Ez = np.array([sum(J[q] * lab[b, q] * lab[b, q + 1] for q in range(9)) for b in range(1024)])
say(f"   its bond energy sum J_i s_i s_i+1 over ALL 1024 label states: mean = {Ez.mean():+.3e}, "
    f"coherence = {np.abs(Ez).mean() / np.abs(J).sum():.6f}")
say(f"   the same quantity on a genuine 10-record register: identical by construction -- the")
say("   cancellation law is a fact about +-1 LABELS, not about records.  IT IS RECORD-BLIND.")
say("   That cuts both ways: it is why no responsive functional survives, and it is why")
say("   cancellation alone can never be evidence that records are involved.")

# ---------------------------------------------------------------- 5. the hypothesis failing
say("")
say("5. ESCAPE (5a): WHAT HAPPENS WHEN THE HYPOTHESIS FAILS -- NON-INDEPENDENT WRITABILITY.")
say("   C-65's carrier: H = sum J_i Z_i Z_i+1, open chain.  Records SEARCHED, writers SEARCHED")
say("   over the whole Pauli group -- nothing nominated (D-18).")
I2n = np.eye(2, dtype=complex)
def pauli_str(ss):
    M = np.array([[1]], dtype=complex)
    for c in ss:
        M = np.kron(M, {'I': I2n, 'X': Xm, 'Y': Ym, 'Z': Zm}[c])
    return M
for n in (3, 4, 5):
    J = [1.0 + 0.3 * i for i in range(n - 1)]
    Zs = [pauli_str(''.join('Z' if q == i else 'I' for q in range(n))) for i in range(n)]
    H = sum(J[i] * Zs[i] @ Zs[i + 1] for i in range(n - 1))
    es = eigenspaces(H)
    mults = [m for _, _, m in es]
    # SEARCH every Pauli for admissibility and for its action on the sign vector of the Z_i
    adm = []
    for ss in itertools.product('IXYZ', repeat=n):
        P = pauli_str(ss)
        if np.linalg.norm(P @ H - H @ P) > 1e-9: continue
        flips = tuple(0 if np.linalg.norm(P @ Zs[i] - Zs[i] @ P) < 1e-9 else 1 for i in range(n))
        adm.append((''.join(ss), flips))
    writers = sorted({f for _, f in adm if any(f)})
    cfg = list(itertools.product((1, -1), repeat=n))
    orb = {}
    for s in cfg:
        reach = frozenset(tuple(s[i] * (-1) ** f[i] for i in range(n)) for f in [(0,) * n] + writers)
        # close the orbit under the group
        changed = True
        while changed:
            changed = False
            new = set(reach)
            for t in reach:
                for f in writers:
                    u = tuple(t[i] * (-1) ** f[i] for i in range(n))
                    if u not in new: new.add(u); changed = True
            reach = frozenset(new)
        orb[s] = reach
    norb = len(set(orb.values())); osize = sorted({len(o) for o in orb.values()})
    say(f"   n={n}: eigenvalue multiplicities {mults}; admissible Paulis {len(adm)}; distinct "
        f"config-flips realised {sorted(set(writers))}")
    say(f"        orbits: {norb} of size {osize}  ->  dim(invariant functionals) = {norb} "
        f"of {2 ** n}   SIMPLY TRANSITIVE: {norb == 1}")
    if norb > 1:
        s0 = cfg[0]; g = orb[s0]
        inv = {s: s[0] * s[1] for s in cfg}
        const_on_orbits = all(len({inv[t] for t in o}) == 1 for o in set(orb.values()))
        say(f"        NON-CONSTANT INVARIANT EXHIBITED: chi_{{1,2}}(s) = s_1 s_2 ; constant on every "
            f"orbit: {const_on_orbits}; distinct values {sorted(set(inv.values()))}")
        say(f"        and BY CONSTRUCTION it is INVARIANT, i.e. NO admissible write changes it.")
say("   THE ESCAPE IS REAL AND IT IS EMPTY.  Dropping independent writability does create")
say("   non-constant invariants -- but 'invariant' MEANS 'does not respond to writing'.  The")
say("   door opens onto the same wall: C-60's dE=0 and C-65's scalar-on-every-eigenspace are two")
say("   instances of exactly this.")
say("=" * 104)
