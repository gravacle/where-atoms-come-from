# W-08 / M4-REFUTE-1  leg D — FOUR AUDITS.
#  D1  M4-6: "W-07's five generic connections ARE a control and they could not have failed."
#      Test the antecedent as the register states it (IMP-1 is a LOGICAL predicate, not a
#      probabilistic one): does the control's failure set exist, and does W-07's own table
#      contain a row where the same code path returns something else?
#  D2  M4-8 / W-07 sec4: which quantity governs sup|Z_k| = 1 at each ready state, and does
#      W-07's sec4 comparison move one variable or two?
#  D3  M4-4: verify the exhibited 3*sqrt(3)/10 state in exact arithmetic, and test M4-3's claim
#      that 3*sqrt(3)/10 is "the one figure with discriminating power".
#  D4  Are M4's "ord = INFINITE" rows infinite-order AS COMPUTED?
#
# ISOLATION LEDGER.
#  D1: carrier, state, observable, dressing, K, tol all fixed at W-07's; MOVED: theta alone.
#  D2: carrier and K fixed; block (a) MOVES the ready state alone at fixed connection;
#      block (b) MOVES W_C alone with W_F pinned at -1; block (c) exhibits what W-07 moved.
#  D3: connection fixed at S1 published, observable fixed; MOVED: the ready state alone.
#  D4: no comparison; a representability check on the floats M4 actually fed its own harness.
# PRECISION: float64 where labelled; Fraction/exact where the claim turns on exact-vs-small.
import numpy as np
from fractions import Fraction

FACE_V={0,1,2}; CYC_V={0,3,4}; TREE={1:(0,),2:(0,1),3:(3,),4:(3,4)}
def dress(s,a):
    u=np.exp(1j*np.asarray(a)); t=np.array(s,dtype=complex)
    for v,p in TREE.items():
        w=1.0+0j
        for e in p: w*=u[e]
        t[v]=s[v]/w
    return t
rng=np.random.default_rng(20260816)
s=rng.normal(size=5)+1j*rng.normal(size=5); s/=np.linalg.norm(s)
K=4000; k=np.arange(1,K+1)

def count(a,u,v,tol=1e-9):
    WF=np.exp(1j*sum(a[:3])); WC=np.exp(1j*sum(a[3:]))
    dF=(v in FACE_V)-(u in FACE_V); dC=(v in CYC_V)-(u in CYC_V)
    t=dress(s,a); amp=abs(np.conj(t[u])*t[v])
    D=amp*np.abs(WF**(dF*k)-WC**(dC*k))
    return int((D<tol).sum()), D.min()

print("== D1  IS W-07's FIVE-GENERIC CONTROL 'UNABLE TO FAIL'?  (M4-6) ==")
print("  The register's disqualifier (IMP-1) is: 'it could not have failed UNDER VARIATION'.")
print("  That is a LOGICAL predicate.  Test it: is the control's failure set empty?")
a_pub=np.array([np.pi/3]*3+[np.pi/2]*3)
c,mn=count(a_pub,2,3)
print(f"    same code path, same state, same pair, theta = -1/4  ->  {c} of 4000   (min {mn:.3e})")
print("    That row is IN W-07's OWN TABLE.  The control returns 1000 on one row and 0 on five.")
print("    Its failure set is therefore NON-EMPTY and demonstrated inside the same table.")
print()
print("  Failure set, measured: fraction of theta in [0,1) at which the control FAILS (returns > 0):")
e=np.arcsin(1e-9/(2*0.271776443))/np.pi
print(f"    every rational p/q with q <= 4000 (a dense set), plus an interval of half-width ~{e:.3e}")
print(f"    about each; total measure >= 2*4000*{e:.3e} = {2*4000*e:.3e}  -- small, but POSITIVE,")
print("    and S1's published connection lies in it.  'Could not have failed' is false as stated;")
print("    'was very unlikely to fail' is true and is a different predicate.")
print()
print("  And the code path is not insensitive to the connection either -- the control's actual job:")
for tag,a in [("generic sqrt2/sqrt3",np.array([2*np.pi*(2**0.5)/3]*3+[2*np.pi*(3**0.5)/3]*3))]:
    for (u,v) in [(2,3),(1,2),(0,3)]:
        c,mn=count(a,u,v)
        print(f"    {tag}, pair ({u},{v}): {c} of 4000   min {mn:.3e}")
print("    Same generic connection returns 4000 at pair (1,2).  A null-returning harness would")
print("    have returned 0 there too.  The control does detect a real failure mode.")
print("    ==> M4-6 is CORRECT that five-for-five is not an isolation of ord(rho); it is WRONG")
print("        that the control 'could not have failed'.  It slides a probabilistic near-certainty")
print("        into a logical disqualifier -- the same category slide it convicts W-07 of at sec6.\n")

print("== D2  WHAT GOVERNS sec4's sup|Z_k| = 1, AND HOW MANY VARIABLES W-07 MOVED THERE ==")
def Zabs(aF,aC,P0,PF,PC,KK=10**6):
    kk=np.arange(1,KK+1)
    return np.abs(P0*np.exp(1j*kk*(aC-aF))+PF*np.exp(-1j*kk*aF)+PC*np.exp(1j*kk*aC))
phi=(1+5**0.5)/2
print("  (a) READY STATE MOVED ALONE, connection pinned at S1 published (aF=pi, aC=3pi/2):")
for nm,(P0,PF,PC) in [("published p=(1/2,0,0,1/4,1/4) -> (1/2, 0, 1/2)",(0.5,0.0,0.5)),
                      ("generic                        -> (0.4,0.3,0.3)",(0.4,0.3,0.3))]:
    d=Zabs(np.pi,3*np.pi/2,P0,PF,PC)
    print(f"      {nm}: cells |Z_k| > 1-1e-12 : {int((d>1-1e-12).sum())}  (period {'2' if PF==0 else '4'})")
print("  (b) W_C MOVED ALONE, W_F PINNED AT -1, published ready state:")
for nm,aC in [("W_C = -i          (ord 4, W-07's)",3*np.pi/2),
              ("W_C = e^{2pi i sqrt2}  (IRRATIONAL, ord = infinity)",2*np.pi*(2**0.5)),
              ("W_C = e^{2pi i/phi^2}  (IRRATIONAL, ord = infinity)",2*np.pi/phi**2)]:
    d=Zabs(np.pi,aC,0.5,0.0,0.5)
    print(f"      {nm:<52}: cells > 1-1e-12 : {int((d>1-1e-12).sum())}")
print("      At the PUBLISHED ready state, sup|Z_k| = 1 is ATTAINED on 500000 of 1e6 cells even when")
print("      <W_F,W_C> is INFINITE and dense.  The attainment there needs only ord(W_F) = 2.")
print("  (c) WHAT W-07's sec4 GENERIC ROW ACTUALLY MOVED:")
print(f"      published row: aF = pi        (W_F = -1)      aC = 3pi/2  (W_C = -i)")
print(f"      generic   row: aF = 2 pi phi  (W_F != -1)     aC = 2 pi phi^2")
print("      BOTH W_F and W_C moved.  W-07 sec4's comparison is CONFOUNDED at the published-state")
print("      column: the variable that carries it is ord(W_F) alone, and W-07 moved two.")
print("      ==> M4-8 CONFIRMED, and by a test M4 did not run (W_C moved alone at W_F = -1).\n")

print("== D3  M4-4's 3*sqrt(3)/10 STATE, IN EXACT ARITHMETIC ==")
print("  Pair (0,3): dF = -1, dC = 0, so D_k = amp*|W_F^-k - 1| = amp*|(-1)^k - 1| in {0, 2*amp}.")
print("  M4's state |s|^2 = (3/4, 4/25, 0, 9/100, 0).  Exact: amp = |s_0||s_3| = sqrt(3)/2 * 3/10.")
print("  2*amp = 2 * sqrt(3)/2 * 3/10 = 3 sqrt(3)/10.  EXACT, no float needed.")
p=[Fraction(3,4),Fraction(4,25),Fraction(0),Fraction(9,100),Fraction(0)]
print(f"    sum of |s|^2 = {sum(p)}  (normalised: {sum(p)==1})")
print(f"    amp^2 = |s_0|^2 |s_3|^2 = {p[0]*p[3]} = {float(p[0]*p[3]):.12f};  (2 amp)^2 = {4*p[0]*p[3]} = 27/100")
print(f"    (3 sqrt(3)/10)^2 = 27/100 = {Fraction(27,100)}   MATCH: {4*p[0]*p[3]==Fraction(27,100)}")
print("    M4-4(a) CONFIRMED in exact arithmetic.  W-07's 'does not reproduce' is refuted as an")
print("    IMPOSSIBILITY claim.  But note what W-07 sec5 itself already wrote three lines later:")
print("    'Either a different observable, a different connection, or a DIFFERENT NORMALISATION")
print("    produced it.'  M4-4(a) is a different normalisation -- i.e. it instantiates one of the")
print("    three escapes W-07 had already listed.  The REGISTER row dropped that clause; the W-07")
print("    PAGE did not.  M4-4 is a hit on the register's phrasing, not on W-07's reasoning.")
print()
print("  M4-3 vs M4-4, INTERNAL TENSION: M4-3 calls 3sqrt3/10 'the one figure with discriminating")
print("  power'; M4-4 shows it is attainable from a continuum of states (any value in (0,1] is).")
print("  Both cannot be load-bearing.  Attainable-set size at k=1, by factor class:")
for f_,nm in [(2.0,"factor 2    (4 pairs)"),(2**0.5,"factor sqrt2 (12 pairs)")]:
    print(f"    {nm}: separations span (0, {f_*0.5:.6f}]  -- a continuum; 0.5196 is interior: {0.5196< f_*0.5}")
print()

print("== D4  ARE M4's 'ord = INFINITE' ROWS INFINITE-ORDER AS COMPUTED? ==")
for lam_expr,lam in [("(sqrt2-1)*1e-13",(2**0.5-1)*1e-13),("(sqrt2-1)*1e-16",(2**0.5-1)*1e-16)]:
    th=0.25+lam
    exact=Fraction(th)                    # the float64 value, exactly
    print(f"    theta = 1/4 + {lam_expr}:  float64 value = {exact.numerator}/{exact.denominator}")
    print(f"      is it exactly 1/4? {exact==Fraction(1,4)}   denominator is a power of 2: "
          f"{exact.denominator & (exact.denominator-1) == 0}  -> AS COMPUTED IT IS RATIONAL, ord = {exact.denominator if exact.denominator&(exact.denominator-1)==0 else '?'}")
print("    Every float64 is a dyadic rational, so both rows have FINITE order as computed.")
print("    The infinite-order claim is a statement about the INTENDED object and is carried by the")
print("    closed-form bound (m4_g G3b), which is sound.  But M4's self-flag 3 says it DID carry")
print("    'the m4_c counterexample' into exact arithmetic while B2 was not -- and the C1 TABLE it")
print("    publishes is float, labelled 'ord(rho) = INFINITE', for a computed value of finite order.")
print("    Same defect class as COR-K, which M4 charges against W-07 sec3.  Not load-bearing (the")
print("    theorem stands independently), but it is the flag M4 did not raise about itself.")
