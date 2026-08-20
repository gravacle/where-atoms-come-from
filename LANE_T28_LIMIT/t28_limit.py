"""T-28: DOES THE AMENDED DEFINITION REDUCE TO DEF-A IN THE LIMIT?

O-51 adopted the amended clauses on carrier (x) local bath with tolerances, and named DEF-A -- the
exact five clauses -- as the T->0, t_m->infinity, W=0, E_b->infinity CORNER. That relationship was
ASSERTED. If it is false, the 162 FORMAL rows are not theorems about a corner of the amended
definition; they are theorems about something else.

THE LIMIT TO TEST. With no dissipation the Liouvillian is -i[H, .], whose zero modes are exactly the
operators commuting with H -- the commutant. So as t_m -> infinity, slow_modes must return the
commutant, and the amended clause (ii') must return DEF-A's clause (ii).

CONTROLS IN THE SAME TABLE (D-15): a carrier where the commutant is known exactly, and a DISSIPATIVE
run at finite t_m where the slow set must be STRICTLY SMALLER -- if it is not, the instrument cannot
tell the limit from the general case."""
import sys, os, itertools, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'model'))
import grounded as GR
from record_model import RecordModel
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
def word(n,s):
    M=np.array([[1]],dtype=complex)
    for c in s: M=np.kron(M,{'I':I2,'X':X,'Z':Z}[c])
    return M
say("="*100); say("T-28   DOES THE AMENDED DEFINITION REDUCE TO DEF-A?"); say("="*100)
say("")
say("1. WITH NO DISSIPATION, DO THE SLOW MODES BECOME THE COMMUTANT AS t_m -> infinity?")
say(f"   {'carrier':<22}{'dim':>6}{'slow modes at t_m=inf':>24}{'dim commutant (exact)':>24}{'agree':>8}")
CARRIERS=[("[[4,2,2]]", -(word(4,'XXXX')+word(4,'ZZZZ'))),
          ("[[6,4,2]]", -(word(6,'X'*6)+word(6,'Z'*6))),
          ("3-qubit Ising", -(word(3,'ZZI')+word(3,'IZZ')))]
ok1=True
for name,H in CARRIERS:
    d=H.shape[0]
    rates,obs = GR.slow_modes(H,[],np.inf)
    # exact commutant dimension: sum of squares of eigenvalue multiplicities of H
    w=np.linalg.eigvalsh(H); mult={}
    for x in w:
        key=round(float(x),9); mult[key]=mult.get(key,0)+1
    exact=sum(m*m for m in mult.values())
    ok1 &= (len(rates)==exact)
    say(f"   {name:<22}{d:>6}{len(rates):>24}{exact:>24}{str(len(rates)==exact):>8}")
say(f"   -> {'the slow modes ARE the commutant in the limit' if ok1 else 'THEY ARE NOT — DEF-A is not the limit of the amended definition'}")
say("")
say("2. CONTROL — WITH DISSIPATION THE SLOW SET MUST GROW MONOTONICALLY AS t_m SHRINKS, and must be")
say("   STRICTLY SMALLER than the commutant when t_m is long enough to resolve the decay rates.")
say(f"   {'carrier':<22}{'t_m (s)':>12}{'slow modes':>13}{'commutant':>12}{'strictly smaller?':>19}")
sz=np.array([[1,0],[0,-1]],dtype=complex); sp=np.array([[0,1],[0,0]],dtype=complex); sm=sp.conj().T
H2=-(0.5)*sz
ok2=True; prev=None; shrank=False
for tm in (1e9, 1e3, 1.0, 1e-3):
    Ls=[np.sqrt(1e-2)*sm, np.sqrt(1e-3)*sp]
    r,_=GR.slow_modes(H2,Ls,tm)
    r0,_=GR.slow_modes(H2,[],np.inf)
    if prev is not None: ok2 &= (len(r)>=prev)
    if len(r)<len(r0): shrank=True
    prev=len(r)
    say(f"   {'2-level dissipative':<22}{tm:>12.0e}{len(r):>13}{len(r0):>12}{str(len(r)<len(r0)):>19}")
ok2 = ok2 and shrank
say(f"   -> {'monotone in t_m, and strictly smaller than the commutant at long t_m' if ok2 else 'CONTROL FAILED'}")
say("")
say("3. DOES CLAUSE (ii') BECOME CLAUSE (ii)? A record durable to every t_m must commute with H.")
say(f"   {'carrier':<22}{'record':<14}{'rate (no dissip.)':>19}{'||[H,R]||':>13}{'both zero?':>12}")
ok3=True
for name,H in CARRIERS[:2]:
    n=int(np.log2(H.shape[0]))
    R=word(n,'ZZ'+'I'*(n-2))
    c=GR.clause_ii(H,[],R,np.inf)
    cm=np.linalg.norm(H@R-R@H)
    both = c['rate']<1e-12 and cm<1e-9
    ok3 &= both
    say(f"   {name:<22}{'ZZI..':<14}{c['rate']:>19.2e}{cm:>13.2e}{str(both):>12}")
say(f"   -> {'clause (ii-prime) reduces to clause (ii)' if ok3 else 'IT DOES NOT'}")
say("")
say("="*100); say("  READ — from the numbers above"); say("="*100)
if ok1 and ok2 and ok3:
    say("  DEF-A IS THE LIMIT OF THE AMENDED DEFINITION, VERIFIED RATHER THAN ASSERTED. With no")
    say("  dissipation the slow modes at t_m -> infinity are exactly the commutant, on every carrier")
    say("  tested, matching the exact dimension sum(m_E^2) computed independently. A dissipative")
    say("  control at finite t_m returns a strictly smaller set, so the instrument distinguishes the")
    say("  limit from the general case. And clause (ii') reduces to clause (ii).")
    say("")
    say("  CONSEQUENCE: the 162 FORMAL rows are theorems about a genuine corner of the amended")
    say("  definition, not about an unrelated object.")
else:
    say("  THE LIMIT DOES NOT HOLD AS TESTED. The relationship between DEF-A and the amended")
    say("  definition is not established, and the FORMAL rows' scope is unsettled.")
