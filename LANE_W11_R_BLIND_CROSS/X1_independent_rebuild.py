"""X1 -- INDEPENDENT REBUILD.  Does the blind lane's arithmetic hold under a third
implementation with a different internal representation (T = P.D, not entrywise)?
Everything here is re-derived from S1 sec.1, COR-F (S3 audit :160-215) and W-01's register row.
double precision (float64) is the default; the exact-arithmetic checks are in X6."""
import numpy as np, sys
sys.path.insert(0, '/Users/bgm/MB Work/where-atoms-come-from/LANE_W11_R_BLIND_CROSS')
from xcore import *
np.set_printoptions(precision=12, linewidth=200)

print("X1.0  P.D FACTORISATION AGREES WITH AN ENTRYWISE BUILD (my own internal control)")
rng = np.random.default_rng(90210)
worst = 0.0
for _ in range(2000):
    k = Cx('rand', 5, [0,1,2], [0,3,4], rng.uniform(-np.pi,np.pi,3), rng.uniform(-np.pi,np.pi,3))
    for w in 'FC':
        worst = max(worst, np.linalg.norm(k.T(w) - k.T_entrywise(w)))
print("   max || P.D - entrywise T ||, 2000 connections = %.2e" % worst)

print("\nX1.1  COR-F's SEALED EXHIBIT (S3_THE_CROSSING_AUDIT_V001.md:180-190), a4,a5,a6=0.7,1.3,-0.4")
k = Cx('K1/CORF', 5, [0,1,2], [0,3,4], [0.0,0.0,0.0], [0.7,1.3,-0.4])
T = k.T('C')
WC = np.exp(1j*(0.7+1.3-0.4))
rho = np.diag([0.40,0.15,0.15,0.15,0.15]).astype(complex)
print("   || T*T - I ||                 = %.2e        [COR-F prints 0.00e+00]" % np.linalg.norm(T.conj().T@T-np.eye(5)))
print("   T diagonal?                     %s                    [COR-F: False]" % (np.linalg.norm(T-np.diag(np.diag(T)))<1e-14))
print("   T^3                           = %s" % np.array2string(np.diag(np.linalg.matrix_power(T,3)), precision=6))
print("   W_C                           = %.6f%+.6fj   [COR-F: -0.029200+0.999574j]" % (WC.real, WC.imag))
print("   || T^3 - diag(WC,1,1,WC,WC) || = %.2e" % np.linalg.norm(np.linalg.matrix_power(T,3)-np.diag([WC,1,1,WC,WC])))
print("   diag(T rho T*)                = %s   [COR-F: 0.15 0.15 0.15 0.40 0.15]" % np.array2string(np.real(np.diag(T@rho@T.conj().T)), precision=2))

print("\nX1.2  T UNITARY, T^L = M, T GAUGE-COVARIANT, M GAUGE-INVARIANT.  2000 connections, seed 90210")
eu = eF = eC = gT = gM = 0.0
rng = np.random.default_rng(90210)
for _ in range(2000):
    pf, pc = rng.uniform(-np.pi,np.pi,3), rng.uniform(-np.pi,np.pi,3)
    k = Cx('K1', 5, [0,1,2], [0,3,4], pf, pc)
    th = rng.uniform(0,2*np.pi,5); G = np.diag(np.exp(1j*th))
    sh = lambda loop, ph: np.array([ph[j]+th[loop[(j+1)%len(loop)]]-th[loop[j]] for j in range(len(loop))])
    kg = Cx('K1g', 5, [0,1,2], [0,3,4], sh([0,1,2],pf), sh([0,3,4],pc))
    for w, Lp in (('F',3), ('C',3)):
        A = k.T(w)
        eu = max(eu, np.linalg.norm(A.conj().T@A-np.eye(5)))
        d = np.linalg.norm(np.linalg.matrix_power(A,Lp) - k.M(w))
        if w=='F': eF = max(eF,d)
        else: eC = max(eC,d)
        gT = max(gT, np.linalg.norm(kg.T(w) - G@k.T(w)@np.linalg.inv(G)))
        gM = max(gM, np.linalg.norm(kg.M(w) - k.M(w)))
print("   max ||T*T - I||           = %.2e     [blind 3.85e-16 | registrar exhibit only]" % eu)
print("   max ||T_F^3 - M_F||       = %.2e     [blind 1.59e-15 | registrar 4.64e-15]" % eF)
print("   max ||T_C^3 - M_C||       = %.2e     [blind 1.74e-15 | registrar 3.25e-15]" % eC)
print("   max ||T[a^g] - g T[a] g*||= %.2e     [blind 1.53e-15 | registrar 4.78e-15]" % gT)
print("   max ||M[a^g] - M[a]||     = %.2e     [blind 3.87e-15 | registrar did not run]" % gM)

print("\nX1.3  THE DECISIVE TEST, MY OWN CONNECTION AND MY OWN THREE ARMS")
print("   connection: f = 1.0, c = sqrt(2)  -- S4:603 / W-10 N-4, the ONLY generic connection")
print("   the corpus publishes.  (The registrar's lane used f = 2.28, c = 2+sqrt(2), which is")
print("   published nowhere; genericity of that pair is checked in X6.)")
k = K1(1.0, np.sqrt(2.0))
sA, sB, sC = equal_pi_triple(k, [0.40,0.15,0.15,0.15,0.15],
                             moves=[(1,2,0.09),(3,4,-0.13)],       # inside class 10 and class 01
                             phases=[0.0,-1.1,0.35,2.9,-0.62])
St = [sA,sB,sC]
print("   ARMS DIFF (guard against the ZERO-VARIABLE control):")
for nm,s in zip('ABC',St):
    print("      s%s |s|^2 = %s   arg = %s" % (nm, np.array2string(np.abs(s)**2, precision=6),
                                               np.array2string(np.angle(s), precision=4)))
print("      ||sA-sB|| = %.6f  ||sA-sC|| = %.6f  ||sB-sC|| = %.6f" %
      (np.linalg.norm(sA-sB), np.linalg.norm(sA-sC), np.linalg.norm(sB-sC)))
print("      pi identical to %.2e ; pi = %s" % (max(np.abs(k.pi(sA)-k.pi(x)).max() for x in (sB,sC)),
                                               np.array2string(k.pi(sA), precision=12)))
N = 2000
circ = np.array([traj(k, s, N, k.T('F'), k.T('C'), lambda t:(3*t,3*t)) for s in St])
edge = np.array([traj(k, s, N, k.T('F'), k.T('C'), lambda t:(t,t)) for s in St])
print("   CIRCUIT advance (3k,3k), spread over k<=%d : %.3e" % (N, np.ptp(circ,axis=0).max()))
print("   EDGE    advance (n,n),   spread over n<=%d : %.3e" % (N, np.ptp(edge,axis=0).max()))
sp = np.ptp(edge,axis=0)
print("      on 3Z: max %.3e     off 3Z: min %.3e" % (sp[2::3].max(), np.delete(sp, np.s_[2::3]).min()))

print("\nX1.4  BOTH B0b RECONSTRUCTIONS.  The registrar's and the blind lane's B0b are DIFFERENT")
print("      objects (shared edge traversed the SAME way vs OPPOSITE ways).  Every B0b claim")
print("      is run on both; if a claim depends on which, it is not a claim about B0b.")
for B in (B0b_registrar(1.0, np.sqrt(2.0)), B0b_blind(1.0, np.sqrt(2.0))):
    cl = B.classes()
    mult = {''.join(map(str,o)): len(cl[o]) for o in ORDER}
    print("   %-10s gamma_F=%s (L=%d)  gamma_C=%s (L=%d)  classes %s  [S4:575 {00:4,01:1,10:2,11:2}]"
          % (B.name, B.loopF, B.LF, B.loopC, B.LC, mult))
    print("      ||T_F^4 - M_F|| = %.2e   ||T_C^3 - M_C|| = %.2e" %
          (np.linalg.norm(np.linalg.matrix_power(B.T('F'),4)-B.M('F')),
           np.linalg.norm(np.linalg.matrix_power(B.T('C'),3)-B.M('C'))))
    sU = np.sqrt(np.ones(9)/9).astype(complex)
    print("      SENSE U pi = %s   m(P) = %.12f   log(4/9) = %.12f" %
          (np.array2string(B.pi(sU), precision=6), mahler4(*B.pi(sU)), np.log(4/9)))
    wu = np.ones(9)/9.0
    a,b,c3 = equal_pi_triple(B, wu, moves=auto_moves(B, wu, 0.45),
                             phases=[0.,0.5,-1.2,2.0,0.3,-0.8,1.7,0.9,-2.4])
    print("      arm B within-class moves (derived from THIS carrier's classes): %s" % auto_moves(B, wu, 0.45))
    print("      |sB|^2 = %s" % np.array2string(np.abs(b)**2, precision=6))
    Sb=[a,b,c3]
    cc = np.array([traj(B,s,600,B.T('F'),B.T('C'), lambda t:(B.LF*t,B.LC*t)) for s in Sb])
    ee = np.array([traj(B,s,600,B.T('F'),B.T('C'), lambda t:(t,t)) for s in Sb])
    spe = np.ptp(ee,axis=0)
    n12 = (np.arange(1,601) % 12 == 0)
    print("      CIRCUIT spread <=600 : %.2e    EDGE spread <=600 : %.2e" % (np.ptp(cc,axis=0).max(), spe.max()))
    print("      EDGE spread on 12Z: max %.2e   off 12Z: min %.2e" % (spe[n12].max(), spe[~n12].min()))
    same=[n for n in range(1,20001) if n % B.LF==0 and n % B.LC==0 and n//B.LF==n//B.LC]
    print("      n<=20000 with n/L_F == n/L_C : %s" % same)
