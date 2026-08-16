import numpy as np, math, sys
sys.path.insert(0,'/private/tmp/claude-501/-Users-bgm-MB-Work/549fc5c6-445d-490b-bf95-ac9313727b33/scratchpad')
from rx import *
rng = np.random.default_rng(RNG_SEED + 1)

print("="*78)
print("A5b  SU(2): a CLASS-COMPATIBLE merge of two class-(1,1) vertices -- NO pushforward works")
print("="*78)
# THETA carrier: gamma_F = f1 f2 f3 (filled), gamma_C = g1 g2 f3^{-1}... use 4 vertices a,b,c,d
# f1: a->c, f2: c->b, f3: b->a   (face)          g1: a->d, g2: d->b, g3: b->a
Vt = ['a','b','c','d']
Et = [('a','c'),('c','b'),('b','a'),('a','d'),('d','b'),('b','a')]
Th = Carrier('THETA', Vt, Et, [[(0,1),(1,1),(2,1)]], [(0,1),(1,1),(2,1)], [(3,1),(4,1),(5,1)])
print("  THETA carrier  d1 =\n", Th.d1())
print("  betti(b0,b1,b2,chi) =", Th.betti(), "  classes:", Th.classes())
print("  a and b are BOTH class (1,1).  Collapse f3 (b->a): CLASS-COMPATIBLE merge.")

def su2(r):
    x = r.normal(size=4); x /= np.linalg.norm(x)
    return np.array([[x[0]+1j*x[1], x[2]+1j*x[3]], [-x[2]+1j*x[3], x[0]-1j*x[1]]])
def unit(z): return z/np.linalg.norm(z)
I2 = np.eye(2, dtype=complex)
U1,U2,U4,U5,U6 = [su2(rng) for _ in range(5)]     # U_f3 = I  (pulled back through the collapse)
HF = {'a':U2@U1, 'b':U1@U2, 'c':U1@U2, 'd':I2}     # wrong-order guard fixed below
HF = {'a':I2@U2@U1, 'b':U1@I2@U2, 'c':U2@U1@I2, 'd':I2}
HC = {'a':U6@U5@U4, 'b':U4@U6@U5, 'c':I2, 'd':U5@U4@U6}
# after collapsing f3, a and b become one vertex 'A' with H_F(A)=U2 U1, H_C(A)=U6 U5 U4
HFm, HCm = U2@U1, U6@U5@U4
N = 100000
def powers(H,N):
    out = np.empty((N,2,2), dtype=complex); M = np.eye(2, dtype=complex)
    for k in range(N): M = M@H; out[k] = M
    return out
PF = {v: powers(HF[v],N) for v in HF}; PC = {v: powers(HC[v],N) for v in HC}
PFm, PCm = powers(HFm,N), powers(HCm,N)
def contrib(Pf,Pc,s):
    x = np.einsum('kij,j->ki',Pf,s); y = np.einsum('kij,j->ki',Pc,s)
    return np.einsum('ki,ki->k',x.conj(),y)
s = {'a':math.sqrt(0.20)*unit(rng.normal(size=2)+1j*rng.normal(size=2)),
     'b':math.sqrt(0.25)*unit(rng.normal(size=2)+1j*rng.normal(size=2)),
     'c':math.sqrt(0.30)*unit(rng.normal(size=2)+1j*rng.normal(size=2)),
     'd':math.sqrt(0.25)*unit(rng.normal(size=2)+1j*rng.normal(size=2))}
Zsrc = sum(contrib(PF[v],PC[v],s[v]) for v in s)
lamS = float(np.mean(np.log(np.abs(Zsrc))))
Zfix = contrib(PF['c'],PC['c'],s['c']) + contrib(PF['d'],PC['d'],s['d'])
print(f"  source lambda_B (N={N}) = {lamS:.9f}   class weights: (00:0, 10:0.30, 01:0.25, 11:0.45)")
print("  the merged vertex A must carry ONE vector with |t|^2 = 0.45 (that is ALL R specifies).")
lams, zd, best = [], [], None
for th in np.linspace(0, math.pi, 61):
    for ph in np.linspace(0, 2*math.pi, 121):
        t = math.sqrt(0.45)*np.array([math.cos(th/2), math.sin(th/2)*np.exp(1j*ph)])
        Zt = Zfix + contrib(PFm,PCm,t)
        lam = float(np.mean(np.log(np.abs(Zt))))
        lams.append(lam); zd.append(float(np.max(np.abs(Zt[:5000]-Zsrc[:5000]))))
        if best is None or abs(lam-lamS) < best[0]: best = (abs(lam-lamS), th, ph)
lams = np.array(lams)
print(f"    lambda_target over 7381 directions: [{lams.min():.9f}, {lams.max():.9f}]")
print(f"    MIN |lambda_target - lambda_source|            = {best[0]:.6e}")
print(f"    MIN over directions of max_k<=5000 |Z_k^L-Z_k^K| = {min(zd):.6e}   (claim asserts <= 4.5e-16)")
print("""  REASON.  The merged vertex is seen only through  Tr[ (H_F^k)^dag H_C^k  rho ].  When BOTH
  holonomies act and generate a non-abelian group, {(H_F^k)^dag H_C^k} spans M_2(C), so rho is
  fully visible: the source's rank-2 rho = s_a s_a* + s_b s_b* (3 real parameters at fixed trace)
  cannot be matched by any rank-1 t t* (2 parameters).  Deficiency 1: NO pushforward exists.""")

print()
print("="*78)
print("A6  'CODIMENSION-1 PULLBACK IMAGE => THE FAMILY SEES ONLY A MEASURE-ZERO SLICE'")
print("="*78)
# K_A = K1 + extra unfilled triangle at v1  (v5,v6; e7:v1->v5, e8:v5->v6, e9:v6->v1)
KA = Carrier('K_A', ['v0','v1','v2','v3','v4','v5','v6'],
             [('v0','v1'),('v1','v2'),('v2','v0'),('v0','v3'),('v3','v4'),('v4','v0'),
              ('v1','v5'),('v5','v6'),('v6','v1')],
             [[(0,1),(1,1),(2,1)]], [(0,1),(1,1),(2,1)], [(3,1),(4,1),(5,1)])
KB = Carrier('K_B', ['v0','v1','v2','v3','v4','v5'],
             [('v0','v1'),('v1','v2'),('v2','v0'),('v0','v3'),('v3','v4'),('v4','v0'),('v1','v5')],
             [[(0,1),(1,1),(2,1)]], [(0,1),(1,1),(2,1)], [(3,1),(4,1),(5,1)])
print("  K_A = K1 + a third (undesignated) unfilled triangle hung at v1 ;  K_B = K1 + a whisker.")
print("  phi: v6->v5, e7->e7', e8->constant, e9->reverse(e7').  v5,v6 are class (0,0) in K_A and")
print("  v5 is class (0,0) in K_B  => CLASS-COMPATIBLE.")
print(f"  betti K_A {KA.betti()}   gauge-invariant dim = {KA.gauge_invariant_dim()}   (f, c, h)")
print(f"  betti K_B {KB.betti()}   gauge-invariant dim = {KB.gauge_invariant_dim()}   (f, c)")
hit = []
for _ in range(2000):
    a = rng.uniform(0, 2*math.pi, size=7)          # a random connection on K_B
    pb = np.array([a[0],a[1],a[2],a[3],a[4],a[5], a[6], 0.0, -a[6]])   # pullback to K_A
    f_ = (pb[0]+pb[1]+pb[2]) % (2*math.pi); c_ = (pb[3]+pb[4]+pb[5]) % (2*math.pi)
    h_ = (pb[6]+pb[7]+pb[8]) % (2*math.pi)
    hit.append((f_, c_, min(h_, 2*math.pi-h_)))
hit = np.array(hit)
print(f"  2000 random K_B connections pulled back:  max |h_A| = {hit[:,2].max():.3e}  -> image = {{h=0}},")
print(f"     a CODIMENSION-1 subtorus of K_A's 3-torus: measure zero, exactly as the claim says.")
print(f"  BUT the projection onto the pair lambda can see is onto: (f,c) covers T^2 --")
H2,_,_ = np.histogram2d(hit[:,0], hit[:,1], bins=12, range=[[0,2*math.pi]]*2)
print(f"     12x12 grid of (f,c): {int((H2>0).sum())}/144 cells occupied, min count {int(H2.min())}.")
pA = {'v0':0.20,'v1':0.10,'v2':0.10,'v3':0.15,'v4':0.15,'v5':0.15,'v6':0.15}
pB = {'v0':0.20,'v1':0.10,'v2':0.10,'v3':0.15,'v4':0.15,'v5':0.30}
wA, wB = class_weights(KA,pA), class_weights(KB,pB)
print(f"  class weights  source {wA}\n                 target {wB}")
print(f"  lambda_B  source {lambdaB_exact(wA):.12f}   target {lambdaB_exact(wB):.12f}   GAP = {abs(lambdaB_exact(wA)-lambdaB_exact(wB)):.3e}")
print("""  So a codimension-1 pullback image coexists with EXACT preservation of lambda on the whole
  of the relevant torus.  lambda and G read the connection through exactly TWO coordinates,
  (W_F, W_C).  Codimension in the other 23 is codimension in coordinates the functional is blind
  to -- S4 already established that ('topology is inert', d2 enters nowhere).  'Measure-zero
  slice of the connection space' is therefore not a wall; it is a statement about gauge-invariant
  directions that carry no formation datum.""")

print()
print("="*78)
print("A7  THE TERMINAL WALL: b1 = 0 AND gamma_C -> ZERO CHAIN, YET FORMATION OCCURS")
print("="*78)
D = Carrier('D2', ['v0','v1','v2','v5'],
            [('v0','v1'),('v1','v2'),('v2','v0'),('v0','v5')],
            [[(0,1),(1,1),(2,1)]], [(0,1),(1,1),(2,1)], [])     # gamma_C = the ZERO chain
print("  psi : K1 -> D2 sends v3,v4 -> v5, e4 -> e7'(v0->v5), e5 -> constant, e6 -> reverse(e7').")
print("  psi_*(gamma_C) = e7' + 0 - e7' = THE ZERO CHAIN.   D2 d1 =\n", D.d1())
print(f"  betti(b0,b1,b2,chi) = {D.betti()}   ->  b1 = 0, NO free cycle, gamma_C cannot be re-designated.")
pD = {'v0':0.40,'v1':0.15,'v2':0.15,'v5':0.30}
wD = class_weights(D, pD)
print(f"  vertex classes on D2: {D.classes()}     class weights {wD}")
formed, r, S = formation(wD)
print(f"  S = {S}   |S| = {len(S)}   rank G = {r}   FORMATION = {formed}   (G = <u> = <conj(W_F)>)")
print(f"  lambda_B = {lambdaB_exact(wD):.12f}      exact: Jensen gives log max(p00,p10) = {math.log(max(wD[(0,0)],wD[(1,0)])):.12f}")
print("""  A terminal carrier with b1 = 0 on which the pushed-forward formation datum is ALIVE:
  G is non-trivial and lambda is finite and negative.  What dies at b1 = 0 is the DESIGNATION
  of a fresh free cycle (choice C4), not formation.  The claim reads a designation wall as a
  mathematical one, and generalises from one terminal object (S4's B5, on which every class
  collapses to a single one) to every terminal object.""")
