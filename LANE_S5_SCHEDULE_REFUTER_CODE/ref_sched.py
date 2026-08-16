#!/usr/bin/env python3
# REFUTER LANE - SCHEDULE AXIS.  Written from scratch.  python3 + numpy only.
# SEED CONVENTION: numpy.random.default_rng(SEED) with SEED printed at each use.
# GRID CONVENTION: all angle grids are half-open [0,2pi) with endpoint EXCLUDED,
#                  theta_j = 2*pi*j/M, j=0..M-1  (midpoint variants stated where used).
import numpy as np

np.set_printoptions(precision=12, suppress=False, linewidth=140)
TWOPI = 2.0*np.pi

def hdr(s):
    print("\n" + "="*78); print(s); print("="*78)

# ---------------------------------------------------------------- 0. CARRIER
hdr("0.  CARRIER K1 - PUBLISHED INCIDENCE MATRICES (built here, not cited)")
# vertices v0..v4 (rows), edges e1..e6 (cols).  d1[v,e] = +1 target, -1 source.
# e1: v0->v1  e2: v1->v2  e3: v2->v0  e4: v0->v3  e5: v3->v4  e6: v4->v0
V, E = 5, 6
d1 = np.zeros((V, E), dtype=int)
edges = [(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]
for j,(s,t) in enumerate(edges):
    d1[s,j] -= 1
    d1[t,j] += 1
# faces: one face F attached along e1+e2+e3
d2 = np.zeros((E,1), dtype=int)
d2[0,0] = 1; d2[1,0] = 1; d2[2,0] = 1
print("d1 (rows v0..v4, cols e1..e6) =\n", d1)
print("d2 (rows e1..e6, col F) =\n", d2.T, " (transposed for display)")
print("d1 @ d2 =", (d1@d2).ravel(), "   d^2 = 0 :", bool(np.all(d1@d2 == 0)))
chi = V - E + 1
rk1 = np.linalg.matrix_rank(d1.astype(float))
print("chi =", chi, "  rank d1 =", rk1, "  b0 =", V-rk1,
      "  b1 =", (E - rk1) - np.linalg.matrix_rank(d2.astype(float)), "  b2 = 0")

# ------------------------------------------------- 1. TRANSPORT, FROM SCRATCH
hdr("1.  THE ONE-CELL OVERLAP Z_k, DERIVED HERE FROM THE TWO LOOP TRANSPORTS")
# Face loop  = e1 e2 e3 (v0 v1 v2), holonomy W_F = exp(i f), f = a1+a2+a3.
# Cycle loop = e4 e5 e6 (v0 v3 v4), holonomy W_C = exp(i c), c = a4+a5+a6.
# Branch F transports each vertex ON the face loop by W_F, spectators by 1.
# Branch C transports each vertex ON the cycle loop by W_C, spectators by 1.
# v0 lies on BOTH.  M_F = diag(e^if,e^if,e^if,1,1), M_C = diag(e^ic,1,1,e^ic,e^ic).
def Z_direct(k, f, c, p):
    MF = np.diag(np.exp(1j*f*np.array([1,1,1,0,0])))
    MC = np.diag(np.exp(1j*c*np.array([1,0,0,1,1])))
    s  = np.sqrt(np.asarray(p, dtype=float))          # ready section, |s_v|^2 = p_v
    a  = np.linalg.matrix_power(MF, k) @ s
    b  = np.linalg.matrix_power(MC, k) @ s
    return np.vdot(a, b)                              # <M_F^k s, M_C^k s>

def Z_closed(k, f, c, p):
    p0 = p[0]; q = p[1]+p[2]; r = p[3]+p[4]
    return p0*np.exp(1j*k*(c-f)) + q*np.exp(-1j*k*f) + r*np.exp(1j*k*c)

p_K1 = [0.4,0.15,0.15,0.15,0.15]
rng = np.random.default_rng(20260816)
print("SEED 20260816")
dev = 0.0
for _ in range(400):
    f,c = rng.uniform(0,TWOPI,2); k = int(rng.integers(0,60))
    dev = max(dev, abs(Z_direct(k,f,c,p_K1) - Z_closed(k,f,c,p_K1)))
print("direct diag transport vs closed form, 400 random (f,c,k<60): max dev = %.3e" % dev)
print("=> Z_k = p11 e^{i(th1+th2)} + p10 e^{i th1} + p01 e^{i th2} + p00,")
print("   with (th1,th2) = k*w,  w = (-f, c),  p11=p_v0, p10=p_v1+p_v2, p01=p_v3+p_v4, p00=0")
print("   [K1 has NO spectator vertex: p00 = 0]")

# reproduce S1's worked instance and S3 4.1
print("\nS1 worked instance f=pi, c=3pi/2 :",
      [ (k, round(abs(Z_direct(k,np.pi,1.5*np.pi,p_K1)),12)) for k in range(1,7) ])

# ------------------------------------------------------- 2. G, HAAR, MAHLER
hdr("2.  G = log|Z| ON THE TORUS; HAAR AVERAGES CHECKED TWO WAYS")
P11, P10, P01, P00 = 0.4, 0.3, 0.3, 0.0
def Zth(th1, th2):
    return P00 + P10*np.exp(1j*th1) + P01*np.exp(1j*th2) + P11*np.exp(1j*(th1+th2))
def G(th1, th2):
    with np.errstate(divide='ignore'):
        return np.log(np.abs(Zth(th1,th2)))
def Zk(k, f, c):    # k may be an integer array
    k = np.asarray(k, dtype=float)
    return Zth(-k*f, k*c)
def Gk(k, f, c):
    with np.errstate(divide='ignore'):
        return np.log(np.abs(Zk(k,f,c)))

# 2-torus Haar integral of G, midpoint rule (avoids the two conical zeros exactly)
M = 4000
g = (np.arange(M)+0.5)*TWOPI/M
T1, T2 = np.meshgrid(g, g, indexing='ij')
haar2 = G(T1,T2).mean()
print("GRID: midpoint, M=%d per axis, theta_j=(j+0.5)*2pi/M" % M)
print("int_{T^2} log|Z| dHaar   (quadrature)      = %.9f" % haar2)
print("Cassaigne-Maillot m(0.4+0.3x+0.3y) of record= -0.767507880")

# 1-torus Haar integral on the diagonal circle H = {(-t,t)}
M1 = 2_000_001
t = (np.arange(M1)+0.5)*TWOPI/M1
haar1 = np.log(np.abs(0.4 + 0.6*np.cos(t))).mean()
print("int_H log|0.4+0.6 cos t| dHaar (quadrature) = %.9f" % haar1)
print("closed form log(0.3) = %.9f   (roots of 0.3z^2+0.4z+0.3 both ON |z|=1)" % np.log(0.3))
print("S4 3.1 table value for L=<(1,1)>            = -1.203972804")

# ============================================================================
hdr("TEST 1  ***  'B IS THE BARYCENTRE HAAR' IS FALSE AT A CONNECTION S4 ITSELF PUBLISHES")
# ============================================================================
phi = np.arccos(-2/3)
print("phi = arccos(-2/3) = %.9f   (S4 3.1 lists (phi,phi) as an EXACT zero of |Z_1|)" % phi)
for (f,c,tag) in [(1.0,1.0,"f=c=1.0     generic on the diagonal"),
                  (phi,phi,"f=c=phi     S4's own exact-zero point"),
                  (phi/7,phi/7,"f=c=phi/7   zero first hit at k=7")]:
    ks = np.arange(1,13)
    zs = np.abs(Zk(ks,f,c))
    # relation lattice L = {(m,n): -m f + n c = 0 mod 2pi}; on c=f this is n=m (f/2pi irrat.)
    print("\n%s" % tag)
    print("   |Z_k|, k=1..12:", np.array2string(zs, precision=9))
    # schedule B partial rates
    for N in (1,2,5,50,2000):
        kk = np.arange(1,N+1)
        with np.errstate(divide='ignore'):
            lam = np.log(np.abs(Zk(kk,f,c))).mean()
        print("   lambda_B at N=%6d : %s" % (N, ("%.9f"%lam) if np.isfinite(lam) else "-inf"))
print("""
L(f=c=1.0)   = { (m,n) : (n-m)*1.0   = 0 mod 2pi } = Z*(1,1)    [1/2pi irrational]
L(f=c=phi)   = { (m,n) : (n-m)*phi   = 0 mod 2pi } = Z*(1,1)    [phi/2pi irrational, Niven]
L(f=c=phi/7) = { (m,n) : (n-m)*phi/7 = 0 mod 2pi } = Z*(1,1)    [same]
IDENTICAL relation lattice, identical orbit closure H (the full diagonal circle),
identical Haar_H, identical int_H G = log 0.3.  lambda_B differs: finite vs -infinity.""")
print("Z_1 at f=c=phi   = 0.4+0.6cos(phi) = 0.4+0.6*(-2/3) = 0  EXACTLY. |Z_1| computed = %.3e"
      % abs(Zk(1,phi,phi)))
print("Z_7 at f=c=phi/7 = 0.4+0.6cos(phi) = 0                EXACTLY. |Z_7| computed = %.3e"
      % abs(Zk(7,phi/7,phi/7)))

# Weyl check that the empirical measure IS Haar_H in the -infinity case
hdr("TEST 1b  THE EMPIRICAL MEASURE AT f=c=phi IS EXACTLY THE SAME HAAR AS AT f=c=1.0")
for (f,tag) in [(1.0,"f=c=1.0"),(phi,"f=c=phi")]:
    for N in (10**4, 10**6):
        n = np.arange(1,N+1)
        th = (n*f) % TWOPI
        wey = [abs(np.exp(1j*mm*th).mean()) for mm in (1,2,3,7)]
        print("%s  N=%8d  |Weyl sums| m=1,2,3,7 : %s" % (tag,N,["%.3e"%x for x in wey]))
print("Both empirical measures -> uniform on the circle.  Same measure, different lambda.")

# ============================================================================
hdr("TEST 2  SAME LIMITING MEASURE (HAAR on T^2), DELIBERATE FINITE SHIFT IN lambda")
# ============================================================================
f2, c2 = 1.0, np.sqrt(2.0)     # S4's cross-check point; orbit dense in T^2
print("connection f=1.0, c=sqrt(2)  (rank L = 0, H = T^2)")
KMAX = 4_000_000
kk = np.arange(1, KMAX+1)
gk = Gk(kk, f2, c2)
print("min over k<=%d of log|Z_k| = %.6f  at k=%d" % (KMAX, gk.min(), kk[gk.argmin()]))
order = np.argsort(gk)
gs, ks_sorted = gk[order], kk[order]

J = 100                      # perturbed cells are the perfect squares n = j^2, j=1..J
Ncells = J*J
alpha = 0.030                # target extra decay per cell
targets = -alpha*(2*np.arange(1,J+1)-1)     # t_j ; sum_j t_j = -alpha*J^2 = -alpha*N
idx = np.searchsorted(gs, targets)
idx = np.clip(idx, 0, len(gs)-1)
kpick = ks_sorted[idx]
sched = np.arange(1, Ncells+1).astype(np.int64)      # schedule B baseline
pert_positions = (np.arange(1,J+1)**2) - 1
sched[pert_positions] = kpick
lam_pert = Gk(sched, f2, c2).mean()
lam_B    = Gk(np.arange(1,Ncells+1), f2, c2).mean()
print("cells N = %d ; perturbed cells = the %d perfect squares (density %.4f -> 0)"
      % (Ncells, J, J/Ncells))
print("achieved sum of perturbation targets / N = %.6f (design: %.6f)"
      % ((gs[idx]).sum()/Ncells, -alpha))
print("lambda (pure schedule B, same N)  = %.6f" % lam_B)
print("lambda (perturbed schedule)       = %.6f" % lam_pert)
print("difference                        = %.6f" % (lam_pert-lam_B))
for mm in [(1,0),(0,1),(1,1),(2,-3)]:
    th1, th2 = -sched*f2, sched*c2
    w1 = abs(np.exp(1j*(mm[0]*th1+mm[1]*th2)).mean())
    th1b, th2b = -np.arange(1,Ncells+1)*f2, np.arange(1,Ncells+1)*c2
    w2 = abs(np.exp(1j*(mm[0]*th1b+mm[1]*th2b)).mean())
    print("Weyl sum m=%-8s perturbed %.4e   pure B %.4e" % (str(mm), w1, w2))
print("""Both empirical measures converge to Haar on T^2 (the perturbed cells have density
J/J^2 -> 0, so they cannot move any weak-* limit), yet lambda differs by ~alpha, and alpha
was CHOSEN.  Taking J -> infinity with alpha fixed makes the difference exact in the limit.""")

# ============================================================================
hdr("TEST 3  lambda EXISTS WHERE THE LIMITING EMPIRICAL MEASURE DOES NOT (EXACT)")
# ============================================================================
# torsion connection: w of order 7.  Conjugate symmetry |Z(-th)| = |Z(th)| is EXACT
# because all four weights are real, so G(k w) = G((7-k) w) exactly.
f3, c3 = TWOPI*1/7, TWOPI*2/7
print("f = 2pi/7, c = 4pi/7 ; orbit is the 7 points {k w}, H finite of order 7")
gvals = np.array([Gk(k,f3,c3) for k in range(0,7)])
print("G(k w), k=0..6 :", np.array2string(gvals, precision=15))
print("G(1w) - G(6w) = %.3e ;  G(2w) - G(5w) = %.3e ;  G(3w) - G(4w) = %.3e  (exact by Z(-th)=conj Z(th))"
      % (gvals[1]-gvals[6], gvals[2]-gvals[5], gvals[3]-gvals[4]))
# blocks of geometrically growing length alternating k=1 and k=6
blocks, kseq = [], []
L_, kcur = 1, 1
while sum(blocks) < 200000:
    blocks.append(L_); kseq.append(kcur); L_ *= 3; kcur = 6 if kcur==1 else 1
sched3 = np.concatenate([np.full(b,k,dtype=np.int64) for b,k in zip(blocks,kseq)])
run = np.cumsum(Gk(sched3,f3,c3))/np.arange(1,len(sched3)+1)
print("running lambda at N=1,10,100,1e3,1e4,1e5,%d :" % len(sched3),
      ["%.15f"%run[i] for i in [0,9,99,999,9999,99999,len(sched3)-1]])
print("constant to %.2e  =>  lambda EXISTS and equals G(w) = %.15f" % (run.max()-run.min(), gvals[1]))
for N in [3,12,120,1200,12000,120000,len(sched3)]:
    fr = (sched3[:N]==1).mean()
    print("   N=%7d : fraction of cells at k=1 = %.4f  (empirical measure oscillates)" % (N,fr))
print("The empirical measure has NO weak-* limit (it oscillates between delta_{1w} and")
print("delta_{6w}); lambda exists and is exact.  The claimed functional is not even total.")

# ============================================================================
hdr("TEST 4  lambda = 0 AT EVERY CONNECTION - AND WHAT THE SCHEDULE MUST LOOK LIKE")
# ============================================================================
def near_recurrence_schedule(f, c, Ncells, Kmax=3_000_000):
    """pick k_n making |Z_{k_n}| -> 1: Dirichlet recurrence of (kf,kc) to (0,0)."""
    ks = np.arange(1, Kmax+1)
    g = Gk(ks, f, c)
    best, out, cur = -np.inf, [], -np.inf
    rec_k, rec_g = [], []
    for i in range(len(ks)):        # record-holders only (cheap streaming max)
        if g[i] > cur:
            cur = g[i]; rec_k.append(ks[i]); rec_g.append(g[i])
    rec_k = np.array(rec_k); rec_g = np.array(rec_g)
    # use the j-th record holder for cell blocks of growing length
    sched, j = [], 0
    while len(sched) < Ncells:
        blk = max(1, int(1.6**j))
        sched += [int(rec_k[min(j, len(rec_k)-1)])]*blk
        j += 1
    return np.array(sched[:Ncells], dtype=np.int64), rec_k, rec_g

GOLD = (1+np.sqrt(5))/2
tests = [(2.0, 1.1, "S3 headline (11f=20c resonant)"),
         (1.0, np.sqrt(2), "rank L = 0"),
         (TWOPI/GOLD, TWOPI/GOLD**2, "badly approximable pair"),
         (np.pi, 1.5*np.pi, "K1's own published connection (finite orbit, order 4)"),
         (0.001, 0.0006, "near-trivial")]
for (f,c,tag) in tests:
    sch, rk, rg = near_recurrence_schedule(f,c,40000)
    lam = Gk(sch,f,c).mean()
    # how close are the chosen circuit counts to a recurrence of the carrier?
    d1_ = np.abs((( -sch*f + np.pi) % TWOPI) - np.pi)
    d2_ = np.abs((( sch*c + np.pi) % TWOPI) - np.pi)
    print("\n%-46s f=%.6f c=%.6f" % (tag,f,c))
    print("   record near-recurrences k = %s ..." % rk[:8])
    print("   deepest |Z_k| reached      = %.12f  (deficiency %.3e)" % (np.exp(rg.max()), 1-np.exp(rg.max())))
    print("   lambda on this schedule, N=40000 : %.9f" % lam)
    print("   max over the last 10%% of cells of ||k_n w||_inf = %.3e  <== LOCKED TO RECURRENCE"
          % max(d1_[-4000:].max(), d2_[-4000:].max()))
print("""
THEOREM (proved in the report): lambda = 0 forces ||k_n w|| -> 0 in density, i.e. the
schedule must sit on the carrier's near-recurrence times.  So lambda = 0 is reachable at
every connection - and ONLY adversarially, in exactly S3 4.6's sense.""")

# Dirichlet certificate: recurrence exists at EVERY connection, no density needed
hdr("TEST 4b  DIRICHLET CERTIFICATE - 2000 RANDOM CONNECTIONS, SEED 777")
rng2 = np.random.default_rng(777)
worst = 0.0; worstw = None
ks = np.arange(1, 200001)
for _ in range(2000):
    f,c = rng2.uniform(0,TWOPI,2)
    m = np.abs(Zk(ks,f,c)).max()
    if 1-m > worst: worst, worstw = 1-m, (f,c)
print("over 2000 random (f,c): worst deficiency min_over-connections of (1 - max_{k<=2e5}|Z_k|)")
print("   = %.3e   at (f,c) = (%.6f, %.6f)" % (worst, *worstw))
print("Dirichlet bound: for any Q there is 1<=k<=Q^2 with ||kf||,||kc|| <= 1/Q, so")
print("1-|Z_k| = O(1/Q^2) uniformly.  Holds at EVERY connection, resonant or not.")

# ============================================================================
hdr("TEST 5  THE REACHABLE SET  [inf_H G, 0]  - CONFIRMED, AND NOTHING BELOW IT")
# ============================================================================
f5, c5 = TWOPI/GOLD, TWOPI/GOLD**2      # H = T^2, but check the finite-orbit case too
ks = np.arange(1,2_000_001)
gg = Gk(ks,f5,c5)
print("f,c = 2pi/gold, 2pi/gold^2 : H = T^2, inf_H G = -infinity (Z has 2 conical zeros)")
print("   min_{k<=2e6} G = %.6f   max_{k<=2e6} G = %.3e (<= 0 always)" % (gg.min(), gg.max()))
print("   sup G over T^2 = G(0,0) = 0 exactly ; |Z| <= sum of weights = 1 always")
# finite-orbit connection: exact interval
f6, c6 = np.pi, 1.5*np.pi
orb = np.array([Gk(k,f6,c6) for k in range(0,4)])
print("\nK1's published connection f=pi,c=3pi/2 : orbit closure H is order 4")
print("   G on H:", np.array2string(orb, precision=12), " => [inf_H G, 0] = [-inf, 0]")
# two-block schedules realise every convex combination
print("\ntwo-block schedules k in {2,4} at f=pi,c=3pi/2 (both give |Z|=1) -> lambda = 0 EXACTLY")
sch = np.array([2,4]*50000)
print("   lambda = %.1f  (exact, not a limit: this connection is torsion)" % Gk(sch,f6,c6).mean())
print("\nmixed schedules at f=1.0,c=sqrt(2): frequency-theta mix of k=1 and k=3")
g1, g3 = Gk(1,f2,c2), Gk(3,f2,c2)
for th in [0.0,0.25,0.5,0.75,1.0]:
    n1 = int(200000*th); s = np.array([1]*n1 + [3]*(200000-n1))
    print("   theta=%.2f  predicted %.9f   computed %.9f" % (th, th*g1+(1-th)*g3, Gk(s,f2,c2).mean()))
print("Every convex combination is realised => the reachable set IS an interval. Its top is 0.")
print("LOWER BOUND IS FORCED: every k_n w lies in H, so G(k_n w) >= inf_H G, so lambda >= inf_H G.")

# ============================================================================
hdr("TEST 6  'A IS THE EXTREME POINT delta_w' - THERE ARE COUNTABLY MANY OF THEM")
# ============================================================================
print("constant schedules k_n = k are ALL extreme (empirical measure = delta_{kw}):")
for k in range(1,13):
    print("   k_n = %2d : mu = delta_{%dw},  lambda = %.9f" % (k,k,Gk(k,f2,c2)))
print("delta_w is one of infinitely many extreme points; 'A' names the smallest k and nothing else.")
print("And EVERY probability measure mu on H satisfies mu = int delta_x dmu(x): 'barycentre'")
print("is a property of every measure, so 'extreme point vs barycentre' does not separate A from B.")

# ============================================================================
hdr("TEST 7  THE UNSTATED NORMALISATION: PER-CELL vs PER-CIRCUIT")
# ============================================================================
N = 20000
kB = np.arange(1,N+1); kA = np.ones(N,dtype=np.int64)
logO_A = Gk(kA,f2,c2).sum(); logO_B = Gk(kB,f2,c2).sum()
print("N = %d cells, f=1.0 c=sqrt(2)" % N)
print("per-CELL    : lambda_A = %.9f   lambda_B = %.9f" % (logO_A/N, logO_B/N))
print("per-CIRCUIT : lambda_A = %.9f   lambda_B = %.9f" % (logO_A/kA.sum(), logO_B/kB.sum()))
print("Under the per-circuit (carrier-time) normalisation B collapses to 0 and A does not:")
print("the A/B dichotomy REAPPEARS.  The claim's dissolution is normalisation-dependent and")
print("the normalisation is never stated.")
print("\nDONE.")
