# LANE_W11_R_MATH — LEG 7.  THE STRONGEST SURVIVING PRO-READING-A MOVE, RUN.
#
# The move: the intermediate ticks may not be OBSERVATION times.  S3's directed system writes one
# record slot per CELL, and CHOICE LEDGER A2 sets "time = number of CIRCUITS of the loop".  So a
# COR-F-faithful system might still only ever evaluate Z at circuit-completion times -- in which
# case COR-F's transport changes nothing and Reading B loses its exhibit.
#
# RESULT.  The move SUCCEEDS COMPLETELY ON K1 and FAILS COMPLETELY ON B0b, and the reason is
# structural, not numerical: on K1 both loops have length 3, so "one circuit each" and "one edge
# each" are the same clock every third tick.  On B0b, |gamma_F| = 4 and |gamma_C| = 3, so the
# corpus's own operator pair (M_F^k, M_C^k) advances branch F by 4k edges and branch C by 3k --
# THE CIRCUIT CLOCK RUNS AT DIFFERENT SPEEDS IN THE TWO BRANCHES.  Invisibility then holds at
# neither the F-circuit times nor the C-circuit times, only on the lcm sublattice, which is not
# "one write per circuit" for either branch.
#
# CONSEQUENCE FOR THE REGISTRAR'S REPORT: three of its four legs (A, B, D) run on K1, and K1
# CANNOT DECIDE THE QUESTION.  Only leg C can, and leg C was run at one unexamined connection,
# with no rate, and with its headline C2 vacuous (wm6).
import numpy as np, wm0_lib as L
rng=np.random.default_rng(20260817)
def spr(oF,oC,n,sts):
    v=[abs(np.vdot(np.linalg.matrix_power(oF,n)@s,np.linalg.matrix_power(oC,n)@s)) for s in sts]
    return max(v)-min(v)
def timeavg(oF,oC,s,N,step=1):
    xF=s.copy(); xC=s.copy(); tot=0.0; cnt=0
    for i in range(1,N*step+1):
        xF=oF@xF; xC=oC@xC
        if i%step==0:
            z=abs(np.vdot(xF,xC)); tot+= np.log(z) if z>1e-300 else -700.0; cnt+=1
    return tot/cnt

print("== M7a  ON K1 THE TWO CONVENTIONS ARE THE SAME OBJECT AT CIRCUIT-COMPLETION TIMES ==")
car=L.K1(); NV=5; a=np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])
TF,TC=L.Top(car["walkF"],a,NV),L.Top(car["walkC"],a,NV)
MF,MC=L.Mop(car["walkF"],a,NV),L.Mop(car["walkC"],a,NV)
sA=np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
sB=np.sqrt(np.array([0.40,0.30,0.00,0.05,0.25]))+0j
sC=sA*np.exp(1j*np.array([0.0,1.3,-0.7,2.2,0.4]))
sts=(sA,sB,sC)
print(f"   |gamma_F| = |gamma_C| = 3, so T^3 = M in BOTH branches at the SAME tick.")
print(f"   max_k |Z_edge(3k) - Z_circuit(k)|, k <= 40, three states = "
      f"{max(abs(np.vdot(np.linalg.matrix_power(TF,3*k)@s,np.linalg.matrix_power(TC,3*k)@s)-np.vdot(np.linalg.matrix_power(MF,k)@s,np.linalg.matrix_power(MC,k)@s)) for s in sts for k in range(1,41)):.2e}")
print(f"   spread at per-circuit times n = 3,6,9,12 : "
      f"{[f'{spr(TF,TC,n,sts):.1e}' for n in (3,6,9,12)]}")
print(f"   edge rate SAMPLED at circuit times (N=2e5): "
      f"{[f'{timeavg(TF,TC,s,200000,step=3):.9f}' for s in sts]}")
print(f"   circuit rate                       (N=2e5): "
      f"{[f'{timeavg(MF,MC,s,200000):.9f}' for s in sts]}")
print("   -> IDENTICAL.  On K1, adopting COR-F's transport and keeping the corpus's per-circuit")
print("      clock changes NOTHING -- not one figure, not the rate, not the invisibility.")
print("      K1 CANNOT DECIDE THE CONVENTION QUESTION.  It is a carrier on which the two")
print("      conventions coincide wherever the corpus ever looks.\n")

print("== M7b  ON B0b THE CORPUS'S OWN CLOCK RUNS THE TWO BRANCHES AT DIFFERENT EDGE SPEEDS ==")
carB=L.B0b(); NB=9; aB=np.random.default_rng(20260817).uniform(0,2*np.pi,18)
TFb,TCb=L.Top(carB["walkF"],aB,NB),L.Top(carB["walkC"],aB,NB)
MFb,MCb=L.Mop(carB["walkF"],aB,NB),L.Mop(carB["walkC"],aB,NB)
w=np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); w/=w.sum()
wB=w.copy(); wB[0],wB[1]=w[0]+w[1],0.0; wB[3],wB[4]=0.0,w[3]+w[4]; wB[5],wB[8]=w[5]+w[8],0.0
sAb=np.sqrt(w)+0j; sBb=np.sqrt(wB)+0j
sCb=sAb*np.exp(1j*np.random.default_rng(7).uniform(0,2*np.pi,9))
stb=(sAb,sBb,sCb)
clb,_,_=L.classes(carB)
assert max(np.max(np.abs(L.pi_of(sAb,clb)-L.pi_of(s,clb))) for s in (sBb,sCb))<1e-15
print(f"   M_F^k, M_C^k advance branch F by 4k edges and branch C by 3k edges.  At k = 1 that is")
print(f"   4 edges against 3.  The EDGE convention is the one that holds the edge count equal.")
print(f"   {'clock':<44}{'spread at 1st':>15}{'2nd':>12}{'3rd':>12}{'4th':>12}")
for nm,step in (("edge tick, every n","1"),):
    pass
rows=[("per F-circuit  (n = 4,8,12,16)",[4,8,12,16]),
      ("per C-circuit  (n = 3,6,9,12)",[3,6,9,12]),
      ("per lcm        (n = 12,24,36,48)",[12,24,36,48]),
      ("corpus's circuit convention  M^k",None)]
for nm,ns in rows:
    if ns is None:
        v=[spr(MFb,MCb,k,stb) for k in (1,2,3,4)]
    else:
        v=[spr(TFb,TCb,n,stb) for n in ns]
    print(f"   {nm:<44}{v[0]:>15.2e}{v[1]:>12.2e}{v[2]:>12.2e}{v[3]:>12.2e}")
print("   -> the ONLY edge clock that restores invisibility is the lcm one, and 12 ticks is")
print("      3 F-circuits AND 4 C-circuits -- not 'one write per circuit' for either branch.")
print("      The pro-Reading-A move dies on B0b.\n")
print(f"   And the RATES on B0b, which the registrar never computed (N = 2e5):")
print(f"   {'clock':<44}{'state A':>15}{'state B':>15}{'state C':>15}{'spread':>10}")
for nm,step,op in (("edge tick, per tick",1,(TFb,TCb)),
                   ("edge, sampled per F-circuit (step 4)",4,(TFb,TCb)),
                   ("edge, sampled per C-circuit (step 3)",3,(TFb,TCb)),
                   ("edge, sampled per lcm       (step 12)",12,(TFb,TCb)),
                   ("corpus's circuit convention",1,(MFb,MCb))):
    v=[timeavg(op[0],op[1],s,200000//max(1,step//1) if step>1 else 200000,step=step) for s in stb]
    print(f"   {nm:<44}{v[0]:>15.9f}{v[1]:>15.9f}{v[2]:>15.9f}{max(v)-min(v):>10.1e}")
print(f"   m(pi) = N1's registered value on this pi = {L.m_poly(L.pi_of(sAb,clb)):.12f}")

print()
print("== M7c  THE 1.6e-05 AT THE PER-F-CIRCUIT CLOCK WAS SAMPLING NOISE, AND WHAT IS UNDER IT ==")
print("   IS THE DEFECT CLASS THE REGISTER HAS MISNAMED SIX TIMES (W-10 N-7, N-3).")
print("   At 5e4 samples the noise is itself O(1e-5), so the time average cannot decide it.")
print("   Settle it with LEG 2's closed form restricted to the residue subset.")
def coeff_rows(car,a,s):
    NV=car["NV"]; wF,wC=car["walkF"],car["walkC"]; LF,LC=len(wF),len(wC); Lam=int(np.lcm(LF,LC))
    TF_,TC_=L.Top(wF,a,NV),L.Top(wC,a,NV)
    x=np.conj(L.hol(wF,a)); y=L.hol(wC,a); _,F,C=L.classes(car)
    inF=[1 if v in F else 0 for v in range(NV)]; inC=[1 if v in C else 0 for v in range(NV)]
    ix={(0,0):0,(1,0):1,(0,1):2,(1,1):3}; rows={}
    for rho in range(Lam):
        B=np.linalg.inv(np.linalg.matrix_power(TF_,rho%LF))@np.linalg.matrix_power(TC_,rho%LC)
        c=np.zeros(4,dtype=complex)
        for u in range(NV):
            for v in range(NV): c[ix[(inF[u],inC[v])]]+=np.conj(s[u])*s[v]*B[u,v]
        eF,eC=rho//LF,rho//LC
        rows[rho]=np.array([c[0],c[1]*x**eF,c[2]*y**eC,c[3]*x**eF*y**eC])
    return rows
def dominance(d,n=1<<16):
    t=2*np.pi*np.arange(n)/n; X=np.exp(1j*t)
    return float(np.min(np.abs(d[0]+d[1]*X)-np.abs(d[2]+d[3]*X)))
CLOCKS={"every tick (rho = 0..11)":list(range(12)),
        "per F-circuit (rho = 0,4,8)":[0,4,8],
        "per C-circuit (rho = 0,3,6,9)":[0,3,6,9],
        "per lcm (rho = 0)":[0]}
print(f"   {'clock':<34}{'state A':>16}{'state B':>16}{'state C':>16}{'spread':>11}")
for nm,rs in CLOCKS.items():
    vals=[np.mean([L.m_poly(coeff_rows(carB,aB,s)[r],1<<21) for r in rs]) for s in stb]
    print(f"   {nm:<34}{vals[0]:>16.9f}{vals[1]:>16.9f}{vals[2]:>16.9f}{max(vals)-min(vals):>11.2e}")
print("   The per-F-circuit spread is EXACTLY 0, not 1.6e-05.  The time average was noise and I")
print("   was about to score a null.  THE REASON IT IS ZERO IS NOT INVISIBILITY:")
print()
print("   STRUCTURE.  At n = 0 mod L_F the tick B = T_F^0 T_C^{r_C} = T_C^{r_C} is the IDENTITY off")
print("   gamma_C, and classes 00 and 10 both lie OUTSIDE gamma_C.  So c00 = p00 and |c10| = p10")
print("   EXACTLY, for every state and every such n: THE FIRST JENSEN BRANCH IS N1's OWN, always.")
print("   All the state dependence is pushed into the second branch (c01, c11).  Hence:")
print("     whenever N1's first branch DOMINATES the second everywhere on the circle -- W-10 N-3's")
print("     exact mechanism, the one that makes B0b's lambda = log(4/9) closed-form -- the rate at")
print("     the F-circuit clock is EXACTLY N1's lambda for every state with that pi, even though")
print("     |Z_n| itself is visibly state-dependent (2.26e-01) at every one of those very ticks.")
print("   IT IS A PROPERTY OF THE STATE, NOT OF THE CONSTRUCTION.  Exhibited, one variable moved:")
for tag,w in (("registrar's state, p00 = 0.44 (branch 1 dominates)",
               np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11])),
              ("same carrier, p00 = 0.10 (branches cross)         ",
               np.array([.20,.22,.14,.18,.16,.025,.025,.025,.025]))):
    w=w/w.sum(); s=np.sqrt(w)+0j; R=coeff_rows(carB,aB,s); pi=L.pi_of(s,clb)
    ms=[L.m_poly(R[r],1<<21) for r in (0,4,8)]
    print(f"     {tag}")
    print(f"       pi = {np.round(pi,4)}   m(pi) = {L.m_poly(pi):.9f}")
    print(f"       min(|branch1| - |branch2|) at rho = 0,4,8 : "
          f"{dominance(R[0]):+.5f} {dominance(R[4]):+.5f} {dominance(R[8]):+.5f}")
    print(f"       m at rho = 0,4,8 : {ms[0]:.9f} {ms[1]:.9f} {ms[2]:.9f}   mean {np.mean(ms):.9f}")
print("   -> at p00 = 0.44 all three equal m(pi) and the F-circuit clock looks invisible.")
print("      At p00 = 0.10 the branches cross and the same clock gives -0.8675 / -0.9895 / -1.0788.")
print("      A NULL THAT READS TWO WAYS, AND I SCORE IT AS NEITHER: it is not a restoration of")
print("      invisibility and it is not nothing.  It is BRANCH DOMINATION -- the second of the two")
print("      mechanisms W-10 N-7 says no proposed name covers -- reappearing one convention down.")
print("      And it is invisible on K1, where L_F = L_C makes the F-circuit clock the lcm clock.")
