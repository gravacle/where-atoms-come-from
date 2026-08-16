import numpy as np, math, sys
sys.path.insert(0, '/private/tmp/claude-501/-Users-bgm-MB-Work/549fc5c6-445d-490b-bf95-ac9313727b33/scratchpad')
from rx import *

print("="*78); print("R1  VALIDATION -- rebuild the corpus from scratch"); print("="*78)

K = K1()
print("\n-- PUBLISHED INCIDENCE, K1 (rows v0..v4, cols e1..e6)")
print("d1 =\n", K.d1())
print("d2 (rows e1..e6, col F) =\n", K.d2().T)
print("d1 @ d2 =", (K.d1() @ K.d2()).ravel(), " -> d^2 = 0:", np.all(K.d1() @ K.d2() == 0))
b0,b1,b2,chi = K.betti()
print("b0,b1,b2,chi =", b0,b1,b2,chi, "  gauge-invariant dim =", K.gauge_invariant_dim())
print("classes:", K.classes())

# ---- S3 test point
f, c = 2.0, 1.1
p = {'v0':0.4,'v1':0.15,'v2':0.15,'v3':0.15,'v4':0.15}
w = class_weights(K, p)
print("\n-- S3 test point f=2.0 c=1.1 p=(0.4,0.15,0.15,0.15,0.15)")
print("class weights (00,10,01,11) =", tuple(w[k] for k in [(0,0),(1,0),(0,1),(1,1)]))

dev = max(abs(Z_k_direct(K,f,c,p,k) - Z_k_class(w,f,c,k)) for k in range(1,201))
print("max |Z_k(direct matrix) - Z_k(class formula)|, k<=200 :", f"{dev:.3e}")

Zs = [abs(Z_k_direct(K,f,c,p,k)) for k in range(1,4001)]
print("min |Z_k| k<=400  =", f"{min(Zs[:400]):.6f}", " at k =", int(np.argmin(Zs[:400]))+1, "   (S3: 0.024654 at 42)")
print("sup |Z_k| k<=4000 =", f"{max(Zs):.6f}", " at k =", int(np.argmax(Zs))+1, "   (S3: 0.999941 at 377)")
print("|Z_1| =", f"{Zs[0]:.12f}", "  (S3: 0.411271)")

# ---- lambda_B  : exact torus average vs direct schedule B
print("\n-- lambda_B at the S3 test point")
print("EXACT (2-var Mahler, generic-torus value) =", f"{lambdaB_exact(w):.12f}",
      "   [S3 corrected generic value -0.767507880]")
for N in (4000, 200000):
    print(f"   direct schedule B, f=2.0 c=1.1, N={N:>7} : {lambdaB_direct_class(w,f,c,N):.9f}")
print("   NOTE: (f,c)=(2.0,1.1) is EXACTLY RESONANT (-11f+20c=0); direct B converges to the")
print("   SUBTORUS value, not the generic one.  check:  -11*2.0 + 20*1.1 =", -11*2.0+20*1.1)
f2, c2 = 1.0, math.sqrt(2)
for N in (200000, 2000000):
    print(f"   direct schedule B, f=1.0 c=sqrt2, N={N:>7} : {lambdaB_direct_class(w,f2,c2,N):.9f}")

# ---- reproduce the S4 SENSE U / SENSE C table rows I depend on
print("\n-- S4 family rows re-derived independently (SENSE U = uniform on vertices)")
for cr in (K1(), K1_subdivided()):
    pu = {v: 1.0/len(cr.V) for v in cr.V}
    wu = class_weights(cr, pu)
    pc = None
    print(f"  {cr.name:24s} V={len(cr.V):2d} classes={ {k:round(v,6) for k,v in wu.items() if v>0} }")
    print(f"      lambda_B SENSE U (exact) = {lambdaB_exact(wu):.12f}")
wC = {(0,0):0.0,(1,0):0.3,(0,1):0.3,(1,1):0.4}
print(f"      lambda_B SENSE C (0.4,0.3,0.3)  = {lambdaB_exact(wC):.12f}   [S4: -0.767507880358]")
print(f"      m(0.4+0.4x+0.2y)                = {mahler2(0.4,0.4,0.2,0.0):.12f}   [S4 B1 SENSE U: -0.756573585640]")
print(f"      m(5/11+5/11x+1/11y)             = {mahler2(5/11,5/11,1/11,0.0):.12f}   [S4 B1s SENSE U: -0.724759919461]")
print(f"      B0b 4-term m(4/9+2/9x+1/9y+2/9xy) = {mahler2(4/9,2/9,1/9,2/9):.12f}  [erratum: log(4/9) = {math.log(4/9):.12f}]")
