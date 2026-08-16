import numpy as np, math, sys
sys.path.insert(0,'/private/tmp/claude-501/-Users-bgm-MB-Work/549fc5c6-445d-490b-bf95-ac9313727b33/scratchpad')
from rx import *
rng = np.random.default_rng(RNG_SEED)

print("="*78)
print("A4  'FOR ANY CHARGE'  --  a CLASS-COMPATIBLE collapse that moves lambda")
print("="*78)
K, L = K1(), K1_collapse_e2()
print("  phi : K1 -> K1/e2  identifies v1 with v2, both class (1,0).  CLASS-COMPATIBLE.")
print("  L2 d1 (rows v0,w,v3,v4 ; cols e1,e3,e4,e5,e6) =\n", L.d1())
print("  L2 classes:", L.classes(), " betti(b0,b1,b2,chi):", L.betti())

def lam_charge(carrier, p, q, f, c, N):
    """vectorised: Z_k = sum_v p_v u^{k a_v q_v} v^{k b_v q_v},  u=e^{-if}, v=e^{ic}"""
    cls = carrier.classes()
    ks = np.arange(1, N+1)
    Z = np.zeros(N, dtype=complex)
    for v in carrier.V:
        a, b = cls[v]
        Z += p[v]*np.exp(1j*ks*q[v]*(-a*f + b*c))
    return float(np.mean(np.log(np.abs(Z)))), Z

f, c = 1.0, math.sqrt(2)     # generic: 1, sqrt2, 2pi rationally independent
p  = {'v0':0.25,'v1':0.20,'v2':0.15,'v3':0.20,'v4':0.20}
q  = {'v0':1,'v1':1,'v2':2,'v3':1,'v4':1}
pL = {'v0':0.25,'w':0.35,'v3':0.20,'v4':0.20}
print(f"\n  charge q = {q}   (q(v1)=1 but q(v2)=2)")
print(f"  class weights source {class_weights(K,p)}")
print(f"  class weights target {class_weights(L,pL)}  -- IDENTICAL: the class pushforward R is preserved")
lamK, ZK = lam_charge(K, p, q, f, c, 2000000)
print(f"  lambda_B  source = {lamK:.9f}")
for qw in (1,2):
    qL = {'v0':1,'w':qw,'v3':1,'v4':1}
    lamL, ZL = lam_charge(L, pL, qL, f, c, 2000000)
    print(f"  lambda_B  target with q(w)={qw} : {lamL:.9f}   GAP = {abs(lamK-lamL):.3e}   |Z_1| {abs(ZL[0]):.9f} vs source {abs(ZK[0]):.9f}")
print("""  No charge on the target reproduces the source: the source's single class (1,0) emits TWO
  characters u^1 and u^2, and one merged vertex can carry only one.  The invariant that is
  actually preserved is the (class x charge) partition -- strictly finer than loop-incidence.""")

print("\n  --- and the NET sentence dies on ONE carrier, with no map at all ---")
pS = {'v0':0.4,'v1':0.15,'v2':0.15,'v3':0.15,'v4':0.15}
q1 = {v:1 for v in K.V}; q2 = {'v0':1,'v1':2,'v2':2,'v3':2,'v4':2}
l1,_ = lam_charge(K,pS,q1,f,c,4000000); l2,_ = lam_charge(K,pS,q2,f,c,4000000)
wS = {(0,0):0.0,(1,0):0.3,(0,1):0.3,(1,1):0.4}
print(f"  K1: SAME carrier, SAME loop-incidence partition, SAME state pushforward (0,.3,.3,.4)")
print(f"    q=(1,1,1,1,1)  lambda_B direct = {l1:.9f}    exact m(.4+.3x+.3y) = {lambdaB_exact(wS):.12f}")
print(f"    q=(1,2,2,2,2)  lambda_B direct = {l2:.9f}    exact = log(0.3) = {math.log(0.3):.12f}")
print("      derivation: Z_k = .4 u^k v^k + .3 u^2k + .3 v^2k ; on the generic torus")
print("      |Z| = |0.4 e^{i(A+B)} + 0.3 e^{2iA} + 0.3 e^{2iB}| = |0.4 + 0.6 cos(A-B)| ,")
print("      whose Mahler measure is m(0.3 + 0.4T + 0.3T^2) = log 0.3 (both roots on |T|=1).")
print(f"    W-03 records -1.200555 for this row; the exact value is log(0.3) = {math.log(0.3):.12f}.")
print("    Deviation 3.4e-03 -- consistent with a finite-N schedule-B estimate.  FLAG R-F1.")

print()
print("="*78)
print("A5  'FOR ANY STRUCTURE GROUP VIA THE CLASS DENSITY PUSHFORWARD R'  --  SU(2)")
print("="*78)
def su2(r):
    a = r.normal(size=4); a /= np.linalg.norm(a)
    return np.array([[a[0]+1j*a[1], a[2]+1j*a[3]], [-a[2]+1j*a[3], a[0]-1j*a[1]]])
def unit(z): return z/np.linalg.norm(z)
I2 = np.eye(2, dtype=complex)
A_,B_,C_,D_,E_ = [su2(rng) for _ in range(5)]    # U_e1,U_e3,U_e4,U_e5,U_e6 ; U_e2 = I (pullback)
HF  = {'v0':B_@A_, 'v1':A_@B_, 'v2':A_@B_, 'v3':I2, 'v4':I2}
HC  = {'v0':E_@D_@C_, 'v1':I2, 'v2':I2, 'v3':C_@E_@D_, 'v4':D_@C_@E_}
HFL = {'v0':B_@A_, 'w':A_@B_, 'v3':I2, 'v4':I2}
HCL = {'v0':E_@D_@C_, 'w':I2, 'v3':C_@E_@D_, 'v4':D_@C_@E_}

N = 200000
def powers(H, N):
    out = np.empty((N,2,2), dtype=complex); M = np.eye(2, dtype=complex)
    for k in range(N):
        M = M @ H; out[k] = M
    return out
PF = {v: powers(HF[v], N) for v in HF}
PC = {v: powers(HC[v], N) for v in HC}
PFL = {'w': PF['v1']}                       # identical matrix, reuse

sK = {'v0':math.sqrt(0.25)*unit(rng.normal(size=2)+1j*rng.normal(size=2)),
      'v1':math.sqrt(0.20)*unit(rng.normal(size=2)+1j*rng.normal(size=2)),
      'v2':math.sqrt(0.15)*unit(rng.normal(size=2)+1j*rng.normal(size=2)),
      'v3':math.sqrt(0.20)*unit(rng.normal(size=2)+1j*rng.normal(size=2)),
      'v4':math.sqrt(0.20)*unit(rng.normal(size=2)+1j*rng.normal(size=2))}
def contrib(Pf, Pc, s):
    x = np.einsum('kij,j->ki', Pf, s); y = np.einsum('kij,j->ki', Pc, s)
    return np.einsum('ki,ki->k', x.conj(), y)
Zsrc = sum(contrib(PF[v], PC[v], sK[v]) for v in sK)
Zfix = contrib(PF['v0'],PC['v0'],sK['v0']) + contrib(PF['v3'],PC['v3'],sK['v3']) + contrib(PF['v4'],PC['v4'],sK['v4'])
lamK = float(np.mean(np.log(np.abs(Zsrc))))
print(f"  source lambda_B (SU(2), N={N}) = {lamK:.9f}    class weights (0, 0.35, 0.40, 0.25)")
print("  |s_v1>, |s_v2> overlap |<s1,s2>|/(|s1||s2|) = %.6f  (NOT collinear)" %
      (abs(np.vdot(sK['v1'],sK['v2']))/(np.linalg.norm(sK['v1'])*np.linalg.norm(sK['v2']))))
print("  target: the merged vertex w carries ONE vector of norm^2 = 0.35 -- the class DENSITY")
print("  pushforward R fixes the norm and NOTHING ELSE.  Scan the whole Bloch sphere of directions:")
IDN = np.broadcast_to(I2,(N,2,2))
lams, zdev, best = [], [], None
for th in np.linspace(0, math.pi, 41):
    for ph in np.linspace(0, 2*math.pi, 81):
        t = math.sqrt(0.35)*np.array([math.cos(th/2), math.sin(th/2)*np.exp(1j*ph)])
        Ztg = Zfix + contrib(PF['v1'], IDN, t)
        lam = float(np.mean(np.log(np.abs(Ztg))))
        lams.append(lam); zdev.append(float(np.max(np.abs(Ztg[:4000]-Zsrc[:4000]))))
        if best is None or abs(lam-lamK) < best[0]: best = (abs(lam-lamK), th, ph)
lams = np.array(lams)
print(f"    lambda_target ranges over [{lams.min():.9f}, {lams.max():.9f}]  -- SPREAD {lams.max()-lams.min():.6f}")
print(f"    lambda_source = {lamK:.9f}.  R does NOT determine lambda: it leaves an interval of width {lams.max()-lams.min():.3f}.")
print(f"    min over directions of max_{{k<=4000}} |Z_k^L - Z_k^K| = {min(zdev):.6e}  -- NEVER 4.5e-16, for ANY direction.")
print(f"    (a two-parameter family of directions meeting a one-dimensional target does contain")
print(f"     accidental lambda matches: best gap in the scan {best[0]:.2e}.  That is coincidence, not transport.)")
# two source states, identical class density, different lambda
sK2 = dict(sK)
sK2['v1'] = math.sqrt(0.20)*unit(rng.normal(size=2)+1j*rng.normal(size=2))
sK2['v2'] = math.sqrt(0.15)*unit(rng.normal(size=2)+1j*rng.normal(size=2))
Z2 = sum(contrib(PF[v], PC[v], sK2[v]) for v in sK2)
lamK2 = float(np.mean(np.log(np.abs(Z2))))
print(f"\n    AND ON ONE CARRIER, NO MAP: two K1 states with IDENTICAL class density (0,.35,.40,.25)")
print(f"      lambda_B = {lamK:.9f}   vs   {lamK2:.9f}     |difference| = {abs(lamK-lamK2):.6e}")
print( "      => at SU(2) the class density pushforward R is not even a FUNCTION of lambda's argument.")
print("""  STRUCTURAL, not numerical.  The source contributes
      <H^k s_v1, s_v1> + <H^k s_v2, s_v2> = tr( (H^k)^dag rho ),  rho = s1 s1* + s2 s2*  (rank 2)
  and the target can contribute only tr( (H^k)^dag t t* )  (rank 1, same trace).  Agreement for
  all k forces rho rank one, i.e. s_v1 || s_v2.  So a CLASS-COMPATIBLE map fails at SU(2) for
  every non-collinear pair: the class DENSITY pushforward R is not enough data once the holonomy
  is not a scalar.  (W-03 already records the underlying break; the map claim ignores it.)""")
