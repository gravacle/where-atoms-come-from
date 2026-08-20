"""ADVERSARIAL VERIFY for T42_B_WORLDDIST. Independent re-implementations.
Probe A: part (i) axioms re-verified with brute-force (no Dijkstra reuse): BFS-free direct
         all-paths check via independent Bellman-Ford-style relaxation + numpy triangle scan.
Probe B: part (ii) key numbers re-derived: tri-viol counts at coarse p's, spectrum at 1/3,
         shuffle/chain controls, Procrustes.
Probe C: ATTACK the dimension pin: the 1% eigencount threshold is an IMPORTED instrument
         parameter. Recompute dim2 p-range at thresholds 0.5%, 1%, 2%, 5% and finer p grid
         (step 0.005). Also p_max at step 0.001 near 1/3 (near-field slack could push it above).
Probe D: is the p=1/3 'exact' point special, or does e.g. p=0.335 also earn metric+2D?
"""
import numpy as np, itertools

# ---------- regenerate U exactly as lane (parameters cross-checked vs LANE_T34_NAND) ----------
E_CH, EPS0 = 1.602176634e-19, 8.8541878128e-12
q, h, pitch, eps_r = 100*E_CH, 10e-9, 40e-9, 3.9
kq2 = q*q/(4*np.pi*eps_r*EPS0)
L = 8
pos = np.array([(x*pitch, y*pitch) for y in range(L) for x in range(L)])
n = len(pos)
def U_of_r(r): return kq2*(1.0/r - 1.0/np.sqrt(r*r + 4*h*h))
def build_U(P):
    m=len(P); Um=np.zeros((m,m))
    for i in range(m):
        for j in range(i+1,m):
            r=float(np.hypot(*(P[i]-P[j]))); Um[i,j]=Um[j,i]=U_of_r(r)
    return Um
U = build_U(pos)

def delta_of(p):
    d=U.copy(); np.fill_diagonal(d,1.0); d=d**(-float(p)); np.fill_diagonal(d,0.0)
    return d/np.max(d)
def tri_viol(delta, rtol=1e-9):
    v=0
    for k in range(n):
        v+=int(np.count_nonzero(delta > delta[:,[k]]+delta[[k],:]+rtol*np.max(delta)))
    return v
def spectrum(delta):
    J=np.eye(n)-np.ones((n,n))/n
    B=-0.5*J@(delta**2)@J
    return np.linalg.eigvalsh(B)[::-1]

print("PROBE B: coarse-table re-derivation (independent code):")
for p in (0.15,0.20,0.25,0.30,1/3,0.35,0.40,0.50,1.00):
    d=delta_of(p); ev=spectrum(d); l1=ev[0]
    dim=int(np.count_nonzero(ev>0.01*l1)); neg=abs(min(ev.min(),0.0))/l1; tv=tri_viol(d)
    print(f"  p={p:.4f} tv={tv} dim@1%={dim} neg={neg:.4f} l3/l1={ev[2]/l1:.4f}")
dlog = U.copy(); np.fill_diagonal(dlog,np.max(U)); dlog=-np.log(dlog/np.max(U)); np.fill_diagonal(dlog,0.0)
dlog/=np.max(dlog); evL=spectrum(dlog)
print(f"  -logU tv={tri_viol(dlog)} dim@1%={int(np.count_nonzero(evL>0.01*evL[0]))}")

ev3=spectrum(delta_of(1/3))
print(f"  spectrum at 1/3 top6: {[round(float(v/ev3[0]),4) for v in ev3[:6]]} min={ev3.min()/ev3[0]:+.5f}")

print("\nPROBE C: threshold-dependence of the DIMENSION pin (the imported 1% choice):")
pgrid=np.round(np.arange(0.02,1.0001,0.005),4)
res={}
for p in pgrid:
    d=delta_of(p); ev=spectrum(d)
    res[p]=(tri_viol(d), ev)
metric_ps=[p for p in pgrid if res[p][0]==0]
print(f"  p_max (step 0.005) = {max(metric_ps)}")
for thr in (0.005,0.01,0.02,0.05):
    dims={p:int(np.count_nonzero(res[p][1]>thr*res[p][1][0])) for p in metric_ps}
    dmin=min(dims.values()); span=[p for p in metric_ps if dims[p]==dmin]
    print(f"  threshold {thr*100:.1f}%: dim_min={dmin} on p in [{min(span)}, {max(span)}]"
          f"  (count {len(span)})")

print("\nPROBE C2: p_max at step 0.001 near 1/3 (near-field slack attack):")
for p in np.round(np.arange(0.330,0.3451,0.001),4):
    d=delta_of(p)
    print(f"  p={p:.3f} tv={tri_viol(d)}")

print("\nPROBE D: does exactly p=1/3 earn, and nearby p? dim@1% just below/above:")
for p in (0.325,0.330,1/3,0.334,0.336):
    d=delta_of(p); ev=spectrum(d)
    print(f"  p={p:.4f} tv={tri_viol(d)} dim@1%={int(np.count_nonzero(ev>0.01*ev[0]))} l3/l1={ev[2]/ev[0]:.5f}")

# ---------- controls ----------
print("\nPROBE B2: controls, independent rerun:")
rng=np.random.default_rng(42)
iu=np.triu_indices(n,1); vals=U[iu].copy(); rng.shuffle(vals)
Ush=np.zeros_like(U); Ush[iu]=vals; Ush+=Ush.T
dsh=Ush.copy(); np.fill_diagonal(dsh,1.0); dsh=dsh**(-1/3); np.fill_diagonal(dsh,0.0); dsh/=np.max(dsh)
evS=spectrum(dsh)
print(f"  shuffle: tv={tri_viol(dsh)} dim={int(np.count_nonzero(evS>0.01*evS[0]))} neg={abs(evS.min())/evS[0]:.3f}")
chain=np.array([(x*pitch,0.0) for x in range(n)]); Uch=build_U(chain)
dch=Uch.copy(); np.fill_diagonal(dch,1.0); dch=dch**(-1/3); np.fill_diagonal(dch,0.0); dch/=np.max(dch)
evC=spectrum(dch)
print(f"  chain: tv={tri_viol(dch)} dim={int(np.count_nonzero(evC>0.01*evC[0]))} l2/l1={evC[1]/evC[0]:.5f}")

# Procrustes
d3=delta_of(1/3)*np.max((U.copy()**0))  # normalized; scale-free claim anyway
J=np.eye(n)-np.ones((n,n))/n
B=-0.5*J@(d3**2)@J
w,V=np.linalg.eigh(B); order=np.argsort(w)[::-1]
X=V[:,order[:2]]*np.sqrt(np.maximum(w[order[:2]],0.0))
P0=pos-pos.mean(0); X0=X-X.mean(0)
sc=np.linalg.norm(P0)/np.linalg.norm(X0)
Uu,_,Vt=np.linalg.svd(P0.T@(X0*sc)); R=Uu@Vt
err=np.linalg.norm((X0*sc)@R.T-P0)/np.linalg.norm(P0)
print(f"  Procrustes rel RMS = {err:.6f}")

print(f"\n  U(2a)/U(4a) = {U_of_r(2*pitch)/U_of_r(4*pitch):.3f}")
