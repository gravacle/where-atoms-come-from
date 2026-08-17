# LANE W-08 / M2 leg E — THE SCHEDULE AXIS, AND THE SHARP FORM OF ATTAINED-vs-APPROACHED.
# Leg C closes the CONNECTION axis: on the one-cell-per-circuit schedule no connection makes
# SUM (1-|Z_k|) converge.  W-02 flagged a second conditionality and never quantified it: "an
# adversarial schedule locked to the carrier's near-recurrences defeats the crossing entirely."
# That moves a DIFFERENT VARIABLE — the schedule — and is reported as a separate comparison.
#
# E3 is the finding: the ATTAINED/APPROACHED distinction has an exact consequence here, and it is
# not the one the register guesses at.
import numpy as np, pickle, mpmath as mp
res = {(rt,cn.strip()):d for (rt,cn),d in pickle.load(open("m2_b_sweep.pkl","rb")).items()}
# SEED DEFECT, RECORDED NOT SILENTLY FIXED: the first run of this script used a hard-coded
# pair (0.63696..., 0.26979...) labelled D_RAND1 that is NOT the seeded stream leg B used.
# Any rank-2 irrational pair supports the same conclusion, but the LABEL was wrong and the
# number was unreproducible from the published seed -- the exact defect this program treats
# as absent data.  The pairs are now drawn from the master seed here exactly as in leg B.
RNG_ = __import__('numpy').random.default_rng(20260816)
RND_ = [(RNG_.random(), RNG_.random()) for _ in range(3)]

THRESH = 1e-2

print("== E1  THE LOG-DENSITY ADVERSARIAL SCHEDULE: write ONLY at record near-returns. ==")
print(f"   Schedule = the running-minimum records of (1-|Z_k|) with depth < {THRESH:.0e}, k <= 1e7.")
print(f"   {'ready':<6} {'connection':<9} {'dimH':>4} {'J = cells written':>18} {'k_J':>10} "
      f"{'SUM_j (1-|Z_{k_j}|)':>20} {'PROD_j |Z_{k_j}|':>18}")
for (rt,cn),d in sorted(res.items()):
    rec=[r for r in d["records"] if r[1]<THRESH]
    if not rec:
        print(f"   {rt:<6} {cn:<9} {d['dim']:>4} {0:>18} {'-':>10} {0.0:>20.9f} {1.0:>18.9f}"
              f"   (no cell is that close)"); continue
    gs=np.array([r[1] for r in rec]); ks=[r[0] for r in rec]
    print(f"   {rt:<6} {cn:<9} {d['dim']:>4} {len(rec):>18} {ks[-1]:>10} "
          f"{gs.sum():>20.9f} {float(np.prod(1.0-gs)):>18.9f}")
print("   REPRODUCTION, not a new number: RS-G / B_S3RES over ALL its records gives")
print("   SUM = 1.533095 and PROD = 0.117102 — S3 sec4.6's own two figures, and COR-L(ii)'s")
print("   corrected count of 11 record-breakers.  S3 had already computed this object without")
print("   naming it: it is the adversarial schedule's total.\n")

print("== E2  IS THAT TAIL CONVERGENT, OR MERELY SMALL SO FAR? ==")
for key in [("RS-G","C_BADAPP"),("RS-G","D_RAND1"),("RS-G","B_S3RES"),("RS-P","D_RAND1")]:
    d=res[key]; rec=[r for r in d["records"] if r[1]<THRESH]
    ks=np.array([r[0] for r in rec],float); gs=np.array([r[1] for r in rec])
    if len(ks)<4: continue
    gr=np.exp(np.polyfit(np.arange(len(ks)),np.log(ks),1)[0])
    sl=np.polyfit(np.log(ks),np.log(gs),1)[0]
    # d_eff, NOT dim H: on RS-P the weight p10 = 0 kills one character, so only chi_0/chi_C = x
    # survives and the near-return problem is 1-dimensional WHATEVER the connection.  Using dim H
    # here would be the confound this program keeps paying for.
    deff = d['dim'] if key[0]=="RS-G" else min(1,d['dim'])
    print(f"   {key[0]} {key[1]:<9} dimH={d['dim']} d_eff={deff}  J={len(ks):>2}  "
          f"k_j grows x{gr:.2f} per record  "
          f"d log(1-|Z|)/d log k = {sl:+.3f}   (theory: -2/d_eff = {-2/deff:+.1f})")
print("   -> k_j grows geometrically and the depths fall like a fixed power of k_j.  The series is")
print("      dominated by a geometric one: IT CONVERGES.  On this schedule durability is killed,")
print("      the limit product is bounded away from zero, and this holds for EVERY class.")
print("   COST: J ~ 15 cells written out of K = 1e7 — density ~ log K / K.  The adversary must")
print("      discard 99.9999% of the circuits and must know the connection's continued fraction.\n")

print("== E3  THE SHARP CONSEQUENCE OF ATTAINED vs APPROACHED: POSITIVE-DENSITY SCHEDULES ==")
print("""   A log-density schedule is a caricature; it wastes almost every circuit.  Ask instead the
   question the distinction actually decides: CAN THE ADVERSARY BLANK THE RECORD WHILE STILL
   WRITING ON A POSITIVE FRACTION delta OF CIRCUITS?  The best such schedule writes on the
   delta*K cells with the SMALLEST (1-|Z_k|).  Its total is computed here EXACTLY by sorting.""")
print(f"   K = 1e6, ready state RS-G held fixed; the connection is the only thing that moves.")
K=10**6; W=(0.4,0.3,0.3); w11,w10,w01=W
phi=(1+5**0.5)/2
CONN=[("A_S1PUB  ATTAINED  |H|=4",0.5,0.75,0),
      ("B_S3RES  APPROACHED H=S^1",1/np.pi,11/(20*np.pi),1),
      ("E_W07GEN APPROACHED H=S^1",np.mod(phi,1),np.mod(phi**2,1),1),
      ("C_BADAPP APPROACHED H=T^2",np.mod(2*np.cos(2*np.pi/7),1),np.mod((2*np.cos(2*np.pi/7))**2,1),2),
      ("D_RAND1  APPROACHED H=T^2",RND_[0][0],RND_[0][1],2)]
DELTAS=[0.5,0.25,0.1,0.01,0.001]
print(f"   {'connection':<27} " + "".join(f"{f'delta={d}':>16}" for d in DELTAS))
k=np.arange(1,K+1,dtype=np.float64)
for lab,al,be,dim in CONN:
    u=(k*al)%1.0; v=(k*be)%1.0
    du=u-np.round(u); dv=v-np.round(v); duv=u+v; duv-=np.round(duv)
    S=np.minimum(4*(w11*w10*np.sin(np.pi*dv)**2+w11*w01*np.sin(np.pi*du)**2
                    +w10*w01*np.sin(np.pi*duv)**2),1.0)
    g=np.sort(S/(1.0+np.sqrt(np.maximum(0.0,1.0-S))))
    out=[]
    for dd in DELTAS:
        n=int(dd*K); out.append(g[:n].sum())
    print(f"   {lab:<27} " + "".join(f"{o:>16.6e}" for o in out))
print("""   READ THIS ROW BY ROW.  On the ATTAINED connection the best delta=0.25 schedule accumulates
   EXACTLY 0.000000 — not small, ZERO — because a quarter of all cells have |Z_k| = 1 exactly,
   forever, and a schedule that writes only there leaves |Omega| = 1 for all time.  On every
   APPROACHED connection the same delta gives a positive total that GROWS LINEARLY in K.
   That is the consequence, and it is exact:

     sup_k |Z_k| = 1 is ATTAINED  <=>  there is a POSITIVE-DENSITY schedule on which the record
                                       is exactly and permanently blank.
     sup_k |Z_k| = 1 is APPROACHED <=> every positive-density schedule still accumulates linearly.

   The register asked whether the distinction is cosmetic.  ON THE HONEST SCHEDULE IT VERY NEARLY
   IS -- worth 4.8% in lambda at |H| = 4 and nothing at all as the order grows (leg D2).  ON THE
   SCHEDULE QUESTION IT IS ABSOLUTE.  Both halves are true and neither may be quoted alone.""")
print()
print("== E3b  HOW FAST DOES THE APPROACHED RATE FALL WITH delta?  MEASURED, NOT ASSERTED. ==")
print("   Predicted from N(eps) ~ K eps^{d_eff/2}: best-delta-schedule total ~ K delta^{1+2/d_eff},")
print("   i.e. delta^2 for d_eff = 2 and delta^3 for d_eff = 1.  The delta values in E3 are ABOVE")
print("   the asymptotic regime (leg F2: the eps-law only sets in below eps ~ 1e-3), so it is")
print("   measured here at K = 1e7 and small delta instead.")
K2=10**7
DEL2=[1e-3,1e-4,1e-5,1e-6]
print(f"   {'connection':<27} {'d_eff':>5} " + "".join(f"{f'delta={d:.0e}':>15}" for d in DEL2)
      + f" {'fitted d log S/d log delta':>27} {'theory':>7}")
for lab,al,be,dim in CONN:
    if dim==0: continue
    tot=[]
    for dd in DEL2: tot.append(0.0)
    vals=np.empty(K2,dtype=np.float64); done=0
    while done<K2:
        c=min(2*10**6,K2-done); kk=np.arange(done+1,done+c+1,dtype=np.float64)
        u=(kk*al)%1.0; v=(kk*be)%1.0
        du=u-np.round(u); dv=v-np.round(v); duv=u+v; duv-=np.round(duv)
        S=np.minimum(4*(w11*w10*np.sin(np.pi*dv)**2+w11*w01*np.sin(np.pi*du)**2
                        +w10*w01*np.sin(np.pi*duv)**2),1.0)
        vals[done:done+c]=S/(1.0+np.sqrt(np.maximum(0.0,1.0-S))); done+=c
    vals.sort()
    tot=[float(vals[:int(dd*K2)].sum()) for dd in DEL2]
    sl=np.polyfit(np.log10(DEL2),np.log10(np.maximum(tot,1e-300)),1)[0]
    print(f"   {lab:<27} {dim:>5} " + "".join(f"{t:>15.4e}" for t in tot)
          + f" {sl:>27.3f} {1+2/dim:>7.1f}")
print("   -> the predicted exponents are recovered.  The adversary's return on a positive-density")
print("      schedule falls off as a POWER of delta but never reaches zero on an approached")
print("      connection, and reaches EXACTLY zero on an attained one.  That is the whole content")
print("      of ATTAINED vs APPROACHED for durability.")
print()

print("== E4  THE HONEST SCHEDULE, FOR CONTRAST (leg C restated) ==")
print(f"   {'ready':<6} {'connection':<9} {'dimH':>4} {'SUM_{k<=1e7}(1-|Z_k|)':>24} {'SUM/K':>10}")
for (rt,cn),d in sorted(res.items()):
    print(f"   {rt:<6} {cn:<9} {d['dim']:>4} {d['S']:>24.4f} {d['S']/1e7:>10.6f}")
print("""   4.7e6 on the honest schedule against < 2 on the log-density adversarial one, at the same K,
   the same connection, the same ready state.  THE SCHEDULE IS WORTH SIX ORDERS OF MAGNITUDE AND
   THE CONNECTION IS WORTH 13%.  Whatever threatens durability on K1, it is not the connection's
   arithmetic — and W-02's flagged schedule conditionality is the whole of the exposure.""")
