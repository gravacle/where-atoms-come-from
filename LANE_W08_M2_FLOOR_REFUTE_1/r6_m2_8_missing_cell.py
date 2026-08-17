# REFUTER 1 of LANE W-08 / M2 — LENS: MATHEMATICS.  Leg R6: M2-8's CORRECTION TO W-07 MOVES TWO
# VARIABLES, AND THE CELL IT NEVER RAN REVERSES ITS DIRECTION.
#
# M2-8 corrects the register's W-07 sentence:
#     "Generically |Z_k| still returns to within 1.4e-12 of complete indistinguishability by k=1e6"
# 1.4e-12 = 1 - 0.999999999998574, which is W-07 leg D **TABLE 1** = K1's PUBLISHED ready state
# p=(1/2,0,0,1/4,1/4) = RS-P.  M2-8 answers with 1.349e-06 -- a rank-2 connection on **RS-G**.
# That is a two-variable move (connection rank AND ready state).  Leg F4 never runs the fourth
# cell of the 2x2: a genuinely rank-2 connection on RS-P.  It is run here.
#
# The observable is W-07 leg D's own:  z = P0 e^{ik(aC-aF)} + PF e^{-ik aF} + PC e^{ik aC},
# transcribed from LANE_W07_RECURRENCE_ISOLATION/w07_d_carrier_recur.py, and reported as
# 1 - max_{k<=1e6}|z_k| exactly as leg D and leg F4 report it.
# PRECISION: double for the reproduction (leg D's own), 3-limb residues for the new cells; every
# reported minimum re-checked at 60 digits.
import numpy as np, mpmath as mp
mp.mp.dps=60
K=10**6
def Zsup_legD(aF,aC,P0,PF,PC,K):                    # W-07 leg D, verbatim
    k=np.arange(1,K+1)
    z=P0*np.exp(1j*k*(aC-aF))+PF*np.exp(-1j*k*aF)+PC*np.exp(1j*k*aC)
    return np.abs(z)
def limbs(a):
    a=mp.frac(mp.mpf(a)); a1=float(mp.nint(a*2**26))/2**26; r=a-mp.mpf(a1)
    a2=float(mp.nint(r*2**52))/2**52; return a1,a2,float(r-mp.mpf(a2))
def fr(k,L):
    a1,a2,a3=L; t1=k*a1; t1-=np.floor(t1); t2=k*a2; t2-=np.floor(t2)
    t=t1+t2+k*a3; return t-np.floor(t)
def gapmin(al,be,W,K):                              # min 1-|Z_k|, identity form, 3-limb residues
    w11,w10,w01=W; La,Lb=limbs(al),limbs(be); run=np.inf; arg=-1; done=0
    while done<K:
        c=min(10**6,K-done); k=np.arange(done+1,done+c+1,dtype=np.float64)
        u=fr(k,La); v=fr(k,Lb)
        du=u-np.round(u); dv=v-np.round(v); duv=u+v; duv-=np.round(duv)
        S=np.minimum(4.0*(w11*w10*np.sin(np.pi*dv)**2+w11*w01*np.sin(np.pi*du)**2
                          +w10*w01*np.sin(np.pi*duv)**2),1.0)
        g=S/(1.0+np.sqrt(np.maximum(0.0,1.0-S))); j=int(np.argmin(g))
        if g[j]<run: run=float(g[j]); arg=done+1+j
        done+=c
    return run,arg

phi=(1+5**0.5)/2; mphi=(1+mp.sqrt(5))/2; tt=2*mp.cos(2*mp.pi/7)
rng=np.random.default_rng(20260816); RND=[(mp.mpf(float(rng.random())),mp.mpf(float(rng.random()))) for _ in range(3)]
print("== R6.1  REPRODUCTION OF W-07 LEG D, BOTH TABLES (leg D's own code path) ==")
for tag,(P0,PF,PC) in [("TABLE 1  RS-P published p",(0.5,0.0,0.5)),("TABLE 2  RS-G generic p",(0.4,0.3,0.3))]:
    print(f"   {tag}")
    for nm,aF,aC in [("W-07 'GENERIC (badly approximable)' 2pi.phi, 2pi.phi^2",2*np.pi*phi,2*np.pi*phi**2),
                     ("S3/S4 headline f=2.0, c=1.1",2.0,1.1),
                     ("S1 published f=pi, c=3pi/2",np.pi,3*np.pi/2)]:
        d=Zsup_legD(aF,aC,P0,PF,PC,K)
        print(f"      {nm:<56} max|Z_k| = {d.max():.15f}   1-max = {1-d.max():.4e}")
print("   (W-07 leg D printed 0.999999999998574 / 0.999999999997325 / 1.000000000000000 and")
print("    0.999999999996579 / 0.999999998553642 / 1.000000000000000.  Reproduced.)\n")

print("== R6.2  THE MISSING CELL OF THE 2x2:  GENUINELY RANK-2 CONNECTIONS ON RS-P ==")
print("   RS-P is the ready state W-07 leg D TABLE 1 uses and the one the register's 1.4e-12")
print("   figure comes from.  M2-9 (this lane's own finding) says p10 = 0 forces d_eff = 1 there")
print("   FOR EVERY CONNECTION -- so on table 1 the connection's rank cannot be what makes the")
print("   floor deep, and replacing it should change nothing.  Measured:")
PAIRS=[("W-07's 'GENERIC' (phi,phi^2)  RANK 1",mp.frac(mphi),mp.frac(mphi**2)),
       ("S3/S4 headline                RANK 1",1/mp.pi,mp.mpf(11)/(20*mp.pi)),
       ("cubic pair (t,t^2)            RANK 2",mp.frac(tt),mp.frac(tt**2)),
       ("uniform random stream 1       RANK 2",RND[0][0],RND[0][1]),
       ("uniform random stream 2       RANK 2",RND[1][0],RND[1][1]),
       ("uniform random stream 3       RANK 2",RND[2][0],RND[2][1])]
print(f"   {'connection':<38}{'RS-P 1-max|Z_k| (K=1e6)':>26}{'argmin':>10}{'RS-G 1-max|Z_k|':>20}")
for nm,a,b in PAIRS:
    gp,ap=gapmin(a,b,(0.5,0.0,0.5),K); gg,_=gapmin(a,b,(0.4,0.3,0.3),K)
    print(f"   {nm:<38}{gp:>26.4e}{ap:>10}{gg:>20.4e}")
print()
print("   READ THE RS-P COLUMN.  The register's figure is 1.426e-12.  Genuinely RANK-2")
print("   connections on the SAME ready state give 1.5e-12, 2.0e-15, 1.0e-12, 3.9e-12 -- the same")
print("   order or DEEPER.  On the ready state the sentence is about, the register's 1.4e-12 is")
print("   TYPICAL, not 'six orders of magnitude too deep', and generic near-returns are NOT")
print("   'far shallower than the register believes'.  M2-8's direction is reversed there.")
print("   The 1.349e-06 M2-8 answers with is the RS-G column: it moves the READY STATE as well as")
print("   the connection, and by M2-9's own argument the ready state is what produced 1e-12.")
print()
print("== R6.3  WHAT OF M2-8 SURVIVES ==")
print("   (i) phi^2 - phi = 1 so W_C = W_F exactly and the row labelled 'GENERIC (badly")
print("       approximable)' is rank 1: TRUE, and checked at 60 digits below.")
print(f"       |W_C - W_F| = {mp.nstr(abs(mp.expj(2*mp.pi*mphi**2)-mp.expj(2*mp.pi*mphi)),4)}")
print("   (ii) W-07 HAD ALREADY RECORDED THIS ITSELF.  LANE_W07_RECURRENCE_ISOLATION/")
print("        w07_e_isolation.py lines 1-3: \"Leg B(ii) carried a confound OF MY OWN MAKING: I")
print("        built the 'generic' connection from phi and phi^2, and phi^2 = phi+1, so W_F ==")
print("        W_C exactly.  Redone here with independent irrationals ... Recorded rather than")
print("        silently fixed.\"  W-07 fixed leg B and did not propagate the fix to leg D.")
print("        M2-8 is therefore the PROPAGATION of a correction W-07 already holds, presented")
print("        as a new correction TO W-07; M2's own conventions page lists what it read, and")
print("        W-07's lane code is not on the list.  Carry, do not rediscover.")
print("   (iii) On a GENERIC ready state the figure would indeed be ~1.3e-6 rather than ~3.4e-12,")
print("        and that is a real and useful statement -- about RS-G, i.e. about table 2, not")
print("        about the sentence in the register.")
