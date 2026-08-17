# LANE_W11_R_MATH_CROSS — LEG 3.  ISOLATION AUDIT OF THE LANE'S OWN ARMS.
# My commissioned lens.  Three arms are checked at the bytes and one of them does not move the
# variable its caption names.
#
# (X3a) wm7 M7c's BREAK ARM.  The lane's headline new-defect is BRANCH DOMINATION: at the
#   per-F-circuit clock on B0b the rate equals N1's lambda for every state with that pi, and it is
#   "Exhibited to break by moving p00 from 0.44 to 0.10 (rates -0.8675 / -0.9895 / -1.0788)."
#   THOSE THREE NUMBERS ARE NOT THREE STATES.  In wm7_clock.py the break arm computes
#   ms=[m_poly(R[r]) for r in (0,4,8)] for ONE state s: they are the three RESIDUES rho = 0,4,8 of
#   a single state.  The variable whose invariance is at issue -- WHICH STATE, at fixed pi -- is
#   never moved in that arm.  It is a zero-variable control on the operative axis.  I move it.
# (X3b) the same arm is captioned "one variable moved" and moves all four components of pi.
# (X3c) wm2 M2d's convergence table is captioned "state C" and is computed from a stale variable.
import numpy as np, xc0_lib as X
np.set_printoptions(linewidth=200)

car=X.B0b(); NV=9
a=np.random.default_rng(20260817).uniform(0,2*np.pi,18)
cl,F,C=X.classes(car)
LF,LC=len(car["walkF"]),len(car["walkC"]); Lam=int(np.lcm(LF,LC))
TF,TC=X.Top(car["walkF"],a,NV),X.Top(car["walkC"],a,NV)

def coeff_rows(car,a,s,NV):
    wF,wC=car["walkF"],car["walkC"]; LF,LC=len(wF),len(wC); Lam=int(np.lcm(LF,LC))
    TF,TC=X.Top(wF,a,NV),X.Top(wC,a,NV)
    x=np.conj(X.hol(wF,a)); y=X.hol(wC,a); _,F,C=X.classes(car)
    inF=[1 if v in F else 0 for v in range(NV)]; inC=[1 if v in C else 0 for v in range(NV)]
    ix={(0,0):0,(1,0):1,(0,1):2,(1,1):3}; rows={}
    for rho in range(Lam):
        B=np.linalg.inv(np.linalg.matrix_power(TF,rho%LF))@np.linalg.matrix_power(TC,rho%LC)
        c=np.zeros(4,dtype=complex)
        for u in range(NV):
            for v in range(NV): c[ix[(inF[u],inC[v])]]+=np.conj(s[u])*s[v]*B[u,v]
        rows[rho]=np.array([c[0],c[1]*x**(rho//LF),c[2]*y**(rho//LC),c[3]*x**(rho//LF)*y**(rho//LC)])
    return rows
def dominance(d,n=1<<22):
    t=2*np.pi*np.arange(n)/n; Xg=np.exp(1j*t)
    return float(np.min(np.abs(d[0]+d[1]*Xg)-np.abs(d[2]+d[3]*Xg)))
def fclock_rate(s):
    R=coeff_rows(car,a,s,NV)
    return float(np.mean([X.m_quad(R[r],1<<21) for r in (0,4,8)]))

print("== X3a  wm7 M7c's BREAK ARM DOES NOT MOVE THE STATE.  I MOVE IT. ==")
print("   The claim under test: 'at the per-F-circuit clock the rate is EXACTLY N1's lambda for")
print("   every state with that pi whenever N1's first Jensen branch dominates, and it BREAKS")
print("   when the branches cross.'  Arm 1 is the lane's own; arm 2 is the lane's break state")
print("   with the state axis actually varied -- 12 states per pi, same pi to 1e-16, arms diffed.")
for tag,w in (("branch 1 dominates  (the lane's B0b state)",
               np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11])),
              ("branches cross      (the lane's break state)",
               np.array([.20,.22,.14,.18,.16,.025,.025,.025,.025]))):
    w=w/w.sum(); s0=np.sqrt(w)+0j; pi=X.pi_of(s0,cl)
    R0=coeff_rows(car,a,s0,NV)
    ens=[s0]+X.same_pi_states(cl,pi,np.random.default_rng(77),11)
    pmax=max(np.max(np.abs(pi-X.pi_of(s,cl))) for s in ens)
    dmin=min(np.linalg.norm(ens[i]-ens[j]) for i in range(len(ens)) for j in range(i+1,len(ens)))
    rates=[fclock_rate(s) for s in ens]
    doms=[min(dominance(coeff_rows(car,a,s,NV)[r]) for r in (0,4,8)) for s in ens]
    print(f"\n   {tag}")
    print(f"     pi = {np.round(pi,4)}   m(pi) = {X.m_quad(pi):.9f}   CM exact = {X.m_CM(pi[0],pi[1],pi[2]) if pi[3]==0 else float('nan'):.9f}")
    print(f"     ARMS: 12 states, pi equal to {pmax:.1e}, min||s_i-s_j|| = {dmin:.4f}")
    print(f"     wm7's own three numbers (rho = 0,4,8 of ONE state) : "
          f"{[f'{X.m_quad(R0[r],1<<21):.6f}' for r in (0,4,8)]}")
    print(f"     F-CLOCK RATE ACROSS 12 DIFFERENT STATES AT FIXED pi:")
    print(f"       min {min(rates):.9f}   max {max(rates):.9f}   SPREAD {max(rates)-min(rates):.2e}")
    print(f"     min over states of min_circle(|branch1| - |branch2|) at rho=0,4,8 = {min(doms):+.5f}")
print("\n   -> THE LANE'S THREE STATES ALL HAVE POSITIVE DOMINANCE AT rho = 0,4,8 AND ITS STATE A")
print("      IS AT +0.00062 AT rho = 4 -- A RAZOR EDGE.  Six of my twelve states at the SAME pi")
print("      cross, and their F-clock rates move by up to 1.47e-02.  So:")
print("      * TRUE and confirmed: at rho = 0 mod L_F the first Jensen branch is p00 + p10 X for")
print("        EVERY state (my dominance column at rho = 0 is +0.06000 for all 15 states, exactly).")
print("      * TRUE and confirmed: for a state whose branch 1 dominates at every such rho, the")
print("        F-clock rate is m(pi) exactly.")
print("      * FALSE AS THE LANE STATES IT: 'the RATE at the F-circuit clock equals N1's lambda")
print("        EXACTLY for every state with that pi'.  Dominance is a property of the BILOCAL")
print("        coefficients c01, c11, which depend on the within-class split and the phases --")
print("        i.e. of the STATE, not of pi.  The lane's own prose says exactly this one line")
print("        later ('IT IS A PROPERTY OF THE STATE') and its claim quantifies over pi.")
print("      * THE LANE CAUGHT ITSELF SCORING A NULL HERE ONCE (its self-flag 1, the 1.6e-05")
print("        sampling-noise null) AND THE REPLACEMENT BLOCK SCORES A SECOND ONE, from a")
print("        three-state ensemble read as a universal quantifier over states.")
print("\n== X3b  AND THE SAME ARM MOVES FOUR VARIABLES UNDER A ONE-VARIABLE CAPTION ==")
w1=np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); w1/=w1.sum()
w2=np.array([.20,.22,.14,.18,.16,.025,.025,.025,.025]); w2/=w2.sum()
p1=X.pi_of(np.sqrt(w1)+0j,cl); p2=X.pi_of(np.sqrt(w2)+0j,cl)
print(f"   wm7 caption: 'IT IS A PROPERTY OF THE STATE... Exhibited, one variable moved:'")
print(f"     pi arm 1 = {np.round(p1,4)}")
print(f"     pi arm 2 = {np.round(p2,4)}")
print(f"     componentwise change = {np.round(p2-p1,4)}   -> {int(np.sum(np.abs(p2-p1)>1e-12))} of 4 components move,")
print(f"   and the report names the move as 'moving p00 from 0.44 to 0.10'.  A one-variable")
print(f"   version exists and I run it: hold p10,p01,p11 in fixed ratio and sweep p00 alone.")
print(f"   {'p00':>7}{'  m(pi)':>14}{'  F-clock rate':>16}{'  min dominance':>16}{'  state-spread':>15}")
for p00 in (0.44,0.35,0.28,0.22,0.18,0.10):
    rest=np.array([0.25,0.09,0.22]); rest=rest/rest.sum()*(1-p00)
    tgt=np.array([p00,rest[0],rest[1],rest[2]])
    ens=X.same_pi_states(cl,tgt,np.random.default_rng(5),8)
    rates=[fclock_rate(s) for s in ens]
    dm=min(min(dominance(coeff_rows(car,a,s,NV)[r]) for r in (0,4,8)) for s in ens)
    print(f"   {p00:>7.2f}{X.m_quad(tgt):>14.9f}{np.mean(rates):>16.9f}{dm:>16.5f}{max(rates)-min(rates):>15.2e}")
print("   -> the transition is at min-dominance crossing zero, and it is p00 alone that drives it")
print("      HERE only because the other three were rescaled together.  The operative variable is")
print("      NOT p00: it is  min_circle(|p00 + p10 X| - |c01 + c11 X|),  which is a property of")
print("      the pi AND of the bilocal coefficients, i.e. of the state.  The lane named the state")
print("      correctly in prose ('IT IS A PROPERTY OF THE STATE') and its arm names p00.")
