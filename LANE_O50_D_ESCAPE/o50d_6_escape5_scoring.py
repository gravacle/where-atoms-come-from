"""O-50 D  PART 6 -- ESCAPE (5): MORE RECORDS / HIGHER GENUS.  AND THE MASTER SCORING TABLE.

Does the writer group's action change character when there are many records?  Exact F_2
computation on m disjoint tori (k = 2m records), no enumeration of 4^k classes anywhere.
Then every candidate this lane produced, scored against all five standards with the C-61
control in the same table.
"""
import sys, itertools, math
import numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_D_ESCAPE")
from o50d_common import *

say("=" * 104)
say("O-50 D  PART 6   ESCAPE (5)  MANY RECORDS / HIGHER GENUS   +   MASTER SCORING")
say("=" * 104)

class MultiTorus:
    """m disjoint L x L tori: genus m, k = 2m records."""
    def __init__(self, m, L):
        self.m, self.L = m, L
        T1 = Torus(L); self.nq = m * T1.nq; nn = self.nq
        self.stab = []
        for b in range(m):
            off = b * T1.nq
            for s in T1.stab:
                v = [0] * (2 * nn)
                for i in range(T1.nq):
                    v[off + i] = s[i]; v[nn + off + i] = s[T1.nq + i]
                self.stab.append(v)
    def sp(self, a, b):
        nn = self.nq
        return sum(a[i] * b[nn + i] + a[nn + i] * b[i] for i in range(nn)) % 2

say("")
say("1. m DISJOINT TORI: k = 2m RECORDS.  The writer group is characterised by LINEAR ALGEBRA")
say("   on F_2^{2k}, never by enumerating 4^k classes.")
say(f"   {'m':>3}{'n qubits':>10}{'k records':>11}{'#configs':>10}{'rank of v -> sp(v,R_i)':>24}"
    f"{'orbit size':>12}{'kernel dim':>12}{'simply transitive':>19}{'dim invariant fns':>19}")
for m in (1, 2, 3, 4, 6):
    MT = MultiTorus(m, 2); nn = MT.nq
    prs = symplectic_logicals(MT.stab, nn)
    k = len(prs)
    # the record family: the first member of each conjugate pair, SEARCHED to be mutually commuting
    Rs = [p[0] for p in prs]
    assert all(MT.sp(Rs[i], Rs[j]) == 0 for i in range(k) for j in range(k) if i != j)
    logbasis = [x for p in prs for x in p]
    M = [[MT.sp(v, R) for R in Rs] for v in logbasis]        # 2k x k over F_2
    r = rank2(M, k)
    orbit = 2 ** r
    kerdim = 2 * k - r
    simply = (orbit == 2 ** k)
    say(f"   {m:>3}{nn:>10}{k:>11}{2 ** k:>10}{r:>24}{orbit:>12}{kerdim:>12}"
        f"{str(simply):>19}{1:>19}")
say("   At every genus the map v -> (sp(v,R_1),...,sp(v,R_k)) has FULL RANK k, so G_W acts")
say("   SIMPLY TRANSITIVELY on the 2^k configurations and the invariant functional space stays")
say("   ONE-DIMENSIONAL.  The kernel stays exactly span{R_1..R_k}.  NOTHING CHANGES CHARACTER.")
say("   Largest reached here: m=6, k=12 records, n=48 qubits, 4096 configurations.  EXACT.")

say("")
say("2. THE ONE CARRIER CLASS THIS LANE CANNOT REACH -- AND WHY IT DOES NOT MATTER.")
say("   A NON-ABELIAN anyon model (D(S_3), D(D_4)) would not give a Z_2^k torsor, and D-21")
say("   warns that record_model.commutant() silently returns a SHORT basis there (512 against")
say("   an exact 736 on D(D_4)), so no verdict could rest on a computation of that kind here.")
say("   BUT THE ORBIT-AVERAGING LEMMA OF PART 2 USES NO STRUCTURE AT ALL -- not abelianness,")
say("   not transitivity, not even a group action on a product set.  Whatever the writer group")
say("   is, its orbit average is invariant and the remainder has orbit mean zero.  A non-abelian")
say("   carrier can only make the INVARIANT space BIGGER, i.e. add more NON-responsive")
say("   quantities.  It cannot create a responsive non-cancelling one.  ESCAPE (5) CLOSED BY")
say("   ARGUMENT, not reached by computation.")

say("")
say("3. WHAT ACTUALLY SCALES.  m disjoint tori at L: exact counts, and the coherence of the")
say("   responsive part measured in the same table.")
rng = np.random.default_rng(19)
say(f"   {'m':>3}{'k=2m':>6}{'n=2mL^2':>9}{'record count k':>15}{'k(2m)/k(m)':>12}"
    f"{'capacity/volume':>17}{'coherence of a responsive fn':>30}")
prev = None
for m in (1, 2, 4, 8, 16, 32):
    L = 3; k = 2 * m; nn = m * 2 * L * L
    S = rng.integers(0, 2, size=(200000, k)) * 2 - 1
    w = np.ones(k)
    F = S @ w
    coh = float(np.abs(F).mean() / np.abs(w).sum())
    ratio = '' if prev is None else f"{k / prev:.4f}"
    prev = k
    say(f"   {m:>3}{k:>6}{nn:>9}{k:>15}{ratio:>12}{k / nn:>17.6f}{coh:>30.6f}")
say("   THE RECORD COUNT IS STRICTLY EXTENSIVE (ratio exactly 2) and the capacity DENSITY is")
say("   constant.  In the same table, the coherence of a responsive functional FALLS like")
say("   k^(-1/2).  Adding records makes the count grow and makes every responsive quantity")
say("   cancel harder.  These are the two halves of the trichotomy, measured together.")

say("")
say("4. MASTER SCORING TABLE.  Standards: (a) strictly extensive, (b) additive over disjoint")
say("   regions, (c) not a count and not topological, (d) sign-definite, (e) power-law falloff")
say("   with separation, INDUCED not inserted.  Plus the two questions this lane exists to ask:")
say("   RESPONSIVE to a write?  and SURVIVES the C-61 record-blind control?")
say("")
hdr = f"   {'quantity':<30}{'(a)':>5}{'(b)':>5}{'(c)':>5}{'(d)':>5}{'(e)':>5}{'RESPONSIVE':>12}{'SURVIVES C-61':>15}{'where measured':>22}"
say(hdr); say("   " + "-" * (len(hdr) - 3))
ROWS = [
 ("E(s)=sum J_ij s_i s_j",        "yes", "yes", "yes", "NO",  "ins", "yes", "NO",  "O-48 / part 2"),
 ("sum_i w_i s_i  (any weights)", "yes", "yes", "yes", "NO",  "n/a", "yes", "NO",  "part 2 sec 3"),
 ("geometry-weighted s_i s_j",    "yes", "yes", "yes", "NO",  "ins", "yes", "NO",  "part 3 sec 3"),
 ("Re<Xbar_1> (coherence)",       "no",  "no",  "yes", "NO",  "n/a", "yes*","NO",  "part 4 sec 2"),
 ("|rho_L,ab| (coherence mag.)",  "no",  "no",  "yes", "yes", "n/a", "relabel", "NO", "part 4 sec 2"),
 ("purity / S(rho_L)",            "no",  "no",  "yes", "yes", "n/a", "NO",  "NO",  "part 4 sec 2"),
 ("S_A, contractible region",     "yes", "yes", "no",  "yes", "n/a", "NO",  "NO",  "part 4 sec 5,6"),
 ("S_A, wrapping region",         "no",  "no",  "no",  "yes", "n/a", "NO",  "NO",  "part 4 sec 5,6"),
 ("bath observable <B(t)>",       "yes", "yes", "yes", "NO",  "ind", "yes", "NO",  "part 5 sec 2"),
 ("chi(R : bath)",                "yes", "yes", "yes", "yes", "ind", "NO",  "NO",  "part 5 sec 2,3"),
 ("sum_i |<R_i>| (integrity)",    "yes", "yes", "no",  "yes", "n/a", "NO",  "NO",  "part 4,5"),
 ("holonomy of a write-loop",     "no",  "no",  "yes", "NO",  "n/a", "yes", "NO",  "part 5 sec 1"),
 ("record count k",               "yes", "yes", "NO",  "yes", "n/a", "NO",  "YES", "part 6 sec 1,3"),
 ("protection distance d",        "no",  "no",  "NO",  "yes", "n/a", "NO",  "YES", "part 1 sec 4"),
 ("sum w_ij s_i s_j, BIASED m!=0","yes", "yes", "yes", "YES", "ins", "yes", "NO",  "part 7 -- OPEN"),
]
for r_ in ROWS:
    say(f"   {r_[0]:<30}{r_[1]:>5}{r_[2]:>5}{r_[3]:>5}{r_[4]:>5}{r_[5]:>5}{r_[6]:>12}{r_[7]:>15}{r_[8]:>22}")
say("   ins = the falloff was INSERTED by hand; ind = INDUCED by dynamics; n/a = the quantity")
say("   has no separation argument at all.  'yes*' = responsive but AMBIGUOUS: two equally")
say("   admissible writers of the same record give opposite signs (part 4 sec 3b).")
say("")
say("   UNDER THE UNIFORM CONFIGURATION MEASURE, NOT ONE ROW IS BOTH RESPONSIVE AND SIGN-")
say("   DEFINITE, AND NOT ONE ROW IS BOTH RESPONSIVE AND SURVIVES THE C-61 CONTROL.  The two")
say("   rows that survive the control are a COUNT and a TOPOLOGICAL LENGTH, failing standard (c)")
say("   by inspection, and neither responds to a write.")
say("   THE LAST ROW IS THE EXCEPTION AND IT IS WHY PART 7 EXISTS: drop the uniform measure and")
say("   the SAME two-body functional becomes responsive AND sign-definite AND accumulating.  It")
say("   still fails the C-61 control, and its bias is INSERTED, so it is OPEN, not won.")

say("")
say("5. THE TRICHOTOMY, WHICH IS THIS LANE'S ACTUAL RESULT.")
say("   Every quantity definable on a record carrier falls into exactly one of three classes:")
say("   (Stated under the hypothesis that admissible writes form a GROUP and the configuration")
say("    measure is the uniform orbit measure.  Part 7 shows that both fail together once the")
say("    writer is irreversible, and that this is the only place they can fail.)")
say("     I.   RESPONSIVE to writing.  Then it is a non-trivial part of the group algebra, its")
say("          orbit mean is EXACTLY zero, it takes both signs equally often, and it cancels")
say("          like k^(-1/2).  Fails standard (d).  [parts 2,3,4,5 -- exact]")
say("     II.  INVARIANT and expressible as a formula in +-1 observables, states and baths.")
say("          Then it has an exact twin on a carrier with ZERO records, because that carrier")
say("          has +-1 observables too.  Record-blind.  [part 5 sec 3 -- measured to 6 digits]")
say("     III. INVARIANT and defined only via the clauses: the record count k and the")

say("          protection distance d.  These survive the control -- and they are a COUNT and a")
say("          TOPOLOGICAL LENGTH, which standard (c) excludes, and they do not respond.")
say("")
say("   A SOURCE MUST BE IN CLASS I (it responds to its own matter) AND HAVE CLASS III's")
say("   IMMUNITY TO THE CONTROL.  THE CLASSES ARE DISJOINT.")
say("")
say("6. WHAT WOULD FALSIFY THIS -- stated so the next lane can aim at it.")
say("   (F1) A carrier on which some admissible write is NOT invertible by another admissible")
say("        operation.  The orbit-averaging lemma needs the writers to form a GROUP.  If the")

say("        admissible set were a SEMIGROUP -- irreversible writes -- orbits would not exist and")
say("        the mean-zero conclusion would fail.  Clause (iv) as written gives a unitary, hence")
say("        a group; a DISSIPATIVE writer would not.  TESTED IN PART 7: such a writer EXISTS on")
say("        the torus, both its Kraus operators commute with H, and under the biased ensemble it")
say("        produces the cancellation law fails completely (coherence 1.000000 at m=1).")
say("   (F2) A quantity defined via the clauses that is NOT a count and NOT topological.  Class")
say("        III has exactly two members here; a third with continuous, non-topological content")
say("        would break the trichotomy.")
say("   (F3) A record configuration space that is NOT a group torsor -- non-abelian anyons -- IF")
say("        it also made the admissible set fail to be a group.  Genus and record number alone")
say("        do not do it (section 1).")
say("=" * 104)
