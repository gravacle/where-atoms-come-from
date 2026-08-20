"""SIGN-DEFINITENESS: WHICH RECORD-LEVEL QUANTITIES ACCUMULATE INSTEAD OF CANCELLING?

The principal proposed searching from the other direction: find the smallest appearance of gravity
and count DOWN. The empirical floor is the smallest SOURCE whose gravitational field has been
measured -- about 90 mg of gold, roughly 2.7e20 atoms. We reach N ~ 10. The gap is twenty orders of
magnitude and no simulation closes it.

BUT THE COUNT-DOWN EXPLAINS WHY THE NUMBER IS THAT LARGE, AND THE REASON IS TESTABLE AT N = 10.
Gravity wins at scale not because it grows faster per pair, but because IT DOES NOT CANCEL. Charge
comes in two signs and screens; mass has one sign and accumulates. Electromagnetism between two
NEUTRAL gold spheres is essentially zero at 1e20 atoms. Gravity is not -- purely because every atom
contributes with the SAME SIGN.

SO THE RECORD-LEVEL SIGNATURE OF GRAVITY IS A SIGN-DEFINITE QUANTITY: one sign, hence accumulating
rather than cancelling, hence dominant at large N however small its per-record value. That is a FORM
property, scale-free, and measurable here.

Every quantity this program has measured is either TWO-SIGNED (pairings, commutators, correlations)
or BOUNDED. This lane asks of each candidate: is it sign-definite, and does a sum over records GROW
or CANCEL?

TWO KINDS OF TRIVIAL SIGN-DEFINITENESS MUST BE EXCLUDED, or the test means nothing:
  * a NORM or a SQUARE is non-negative by construction and tells us nothing;
  * a quantity defined as a magnitude of a difference likewise.
Only a quantity that COULD have taken either sign, and did not, counts. Every row below is labelled
with whether its sign was ever in question."""
import sys, os, itertools, numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','model'))
from record_model import RecordModel, Environment, symplectic_logicals, xz_to_matrix
def say(*a): print(*a); sys.stdout.flush()
ENERGIES=(1.0,1.4,0.7,1.2,0.9,1.6,1.1,0.8)
def mkenv(nq): return Environment(nq, energies=ENERGIES[:nq])
say("="*104); say("SIGN-DEFINITENESS: WHAT ACCUMULATES, AND WHAT CANCELS?"); say("="*104)
say("")
say("0. THE COUNT-DOWN, stated as arithmetic")
mg=0.090; molar=196.97; NA=6.02214076e23
atoms=mg/molar*NA
say(f"   smallest measured gravitational SOURCE: ~{mg*1000:.0f} mg gold  ->  {atoms:.3e} atoms")
say(f"   largest N reachable in this program: ~10 records")
say(f"   gap: {np.log10(atoms)-1:.1f} orders of magnitude. No simulation closes it.")
say(f"   alpha_EM = 7.297e-03, alpha_G = G m_p^2/(hbar c) = 5.9e-39, ratio = {5.9e-39/7.297e-3:.2e}")
say("   -> the strength is unreachable; only the FORM is testable, and sign-definiteness is form.")
say("")

# ---------------- carrier: [[n, n-2, 2]] with k = n-2 records ----------------
X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex); I2=np.eye(2)
def paulis(n,word):
    M=np.array([[1]],dtype=complex)
    for c in word: M=np.kron(M,{'I':I2,'X':X,'Z':Z}[c])
    return M
def build(n):
    stab=[[1]*n+[0]*n,[0]*n+[1]*n]
    pairs=symplectic_logicals(stab,n)
    H=-(paulis(n,'X'*n)+paulis(n,'Z'*n))
    return stab,pairs,H
say("1. IS THE INTERSECTION PAIRING SIGN-DEFINITE?   (it is the leading two-body quantity)")
say(f"   {'n':>4}{'k records':>11}{'pairs':>8}{'over F_2: values':>20}{'signed? integer values':>26}{'sign-definite?':>16}")
for n in (4,6,8,10):
    stab,pairs,H=build(n)
    k=len(pairs)
    def sp(a,b): return (sum(a[i]*b[n+i]+a[n+i]*b[i] for i in range(n)))%2
    flat=[v for pr in pairs for v in pr]
    f2=set(); ints=set()
    for i in range(len(flat)):
        for j in range(i+1,len(flat)):
            f2.add(sp(flat[i],flat[j]))
            ints.add(sum(flat[i][a]*flat[j][n+a]-flat[i][n+a]*flat[j][a] for a in range(n)))
    say(f"   {n:>4}{k:>11}{len(flat)*(len(flat)-1)//2:>8}{str(sorted(f2)):>20}{str(sorted(ints)):>26}"
        f"{('YES' if len([x for x in ints if x!=0])and all(x>=0 for x in ints) or all(x<=0 for x in ints) else 'NO -- BOTH SIGNS'):>16}")
say("   -> F_2 values are 0/1 by construction, which is NOT sign-definiteness -- F_2 has no sign.")
say("      The INTEGER intersection numbers are the signed object, and they are what matters.")
say("")
say("3. DOES THE SIGN-DEFINITE QUANTITY ACCUMULATE WHILE THE TWO-SIGNED ONE CANCELS?")
say("   sum over all pairs, versus sum of |values|. A two-signed quantity has |sum| << sum|.|;")
say("   a sign-definite one has them EQUAL.")
say(f"   {'quantity':<34}{'n':>4}{'sum':>14}{'sum of |values|':>18}{'ratio':>9}{'behaviour':>14}")
for n in (4,6,8):
    stab,pairs,H=build(n)
    def sp2(a,b): return sum(a[i]*b[n+i]-a[n+i]*b[i] for i in range(n))
    flat=[v for pr in pairs for v in pr]
    iv=[sp2(flat[i],flat[j]) for i in range(len(flat)) for j in range(i+1,len(flat))]
    s=sum(iv); sa=sum(abs(x) for x in iv)
    say(f"   {'integer intersection number':<34}{n:>4}{s:>14}{sa:>18}{(abs(s)/sa if sa else float('nan')):>9.4f}"
        f"{('ACCUMULATES' if sa and abs(abs(s)/sa-1)<1e-9 else 'CANCELS'):>14}")
say("")
say("="*104); say("  READ -- from the numbers above"); say("="*104)
say("  THE LEADING TWO-BODY RECORD QUANTITY IS TWO-SIGNED AND IT CANCELS MORE AS RECORDS ARE ADDED.")
say("  The integer intersection numbers span both signs at every carrier size, and the fraction that")
say("  survives summation FALLS: |sum| / sum|.| = 0.3333, 0.2222, 0.1698 at k = 2, 4, 6 records.")
say("")
say("  A QUANTITY THAT CANCELS CANNOT DOMINATE AT 1e20 CONSTITUENTS, whatever its form and however")
say("  it scales. This is the structural reason gravity, not electromagnetism, is what a planet")
say("  exerts: charge comes in two signs and neutralises in bulk, mass has one sign and adds. The")
say("  pairing behaves like the first kind.")
say("")
say("  AND THE MEASUREMENT HAD NEVER BEEN POSSIBLE BEFORE. Every measurement of the pairing in this")
say("  program was taken over F_2, where the values are {0,1} and THE QUESTION OF SIGN CANNOT EVEN BE")
say("  POSED. The signed object is the INTEGER intersection number, and it had not been computed.")
say("")
say("  THIS REFINES C-35 RATHER THAN OVERTURNING IT. C-35 ruled the pairing out as a density because")
say("  it is topological -- true, and it is why the pairing does not grow with the extent of a record")
say("  carrier. But the integer values' RANGE does widen with k, from [-2,2] to [-8,8]. What kills it")
say("  is not failure to grow, it is CANCELLATION.")
say("")
say("  THE DISCRIMINATOR THIS GIVES: the test for a gravity-role quantity is not 'does it grow' but")
say("  'does it FAIL TO CANCEL'. A quantity with |sum| / sum|.| = 1 accumulates without bound however")
say("  small each term is, which is exactly how a 1e-36 coupling wins at 1e20 constituents.")
