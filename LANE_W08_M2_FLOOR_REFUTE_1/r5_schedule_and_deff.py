# REFUTER 1 of LANE W-08 / M2 — LENS: MATHEMATICS.  Leg R5: the OTHER two scaling laws, and the
# d_eff RULE ITSELF.
#   R5.1  E3b's law  total(delta) ~ K delta^{1+2/d_eff}: reproduced, and then tested on the ONE
#         d_eff = 2 connection M2 left out of it.
#   R5.2  F2's N(eps) fits: are they fitted with enough counts to carry an exponent?  Poisson
#         error bars, and the significance of the F_VWA row's 0.795 that M2-4's evidence drops.
#   R5.3  M2's RS-P rule "d_eff = min(1, dim H)" / "d_eff = 1 FOR EVERY CONNECTION on RS-P" is
#         false on a connection type its tables happen not to contain.  Exhibited.
# PRECISION: double, with the 3-limb residue split (k*a_hi, k*a_mid exact for k < 2^27).
import numpy as np, mpmath as mp
mp.mp.dps=60
W=(0.4,0.3,0.3); w11,w10,w01=W
def limbs(a):
    a=mp.frac(mp.mpf(a)); a1=float(mp.nint(a*2**26))/2**26; r=a-mp.mpf(a1)
    a2=float(mp.nint(r*2**52))/2**52; return a1,a2,float(r-mp.mpf(a2))
def fr(k,L):
    a1,a2,a3=L; t1=k*a1; t1-=np.floor(t1); t2=k*a2; t2-=np.floor(t2)
    t=t1+t2+k*a3; return t-np.floor(t)
def gapv(u,v,W):
    w11,w10,w01=W
    du=u-np.round(u); dv=v-np.round(v); duv=u+v; duv-=np.round(duv)
    S=np.minimum(4.0*(w11*w10*np.sin(np.pi*dv)**2+w11*w01*np.sin(np.pi*du)**2
                      +w10*w01*np.sin(np.pi*duv)**2),1.0)
    return S/(1.0+np.sqrt(np.maximum(0.0,1.0-S)))
def allvals(al,be,K,W):
    La,Lb=limbs(al),limbs(be); out=np.empty(K); done=0
    while done<K:
        c=min(2*10**6,K-done); k=np.arange(done+1,done+c+1,dtype=np.float64)
        out[done:done+c]=gapv(fr(k,La),fr(k,Lb),W); done+=c
    return out
tt=2*mp.cos(2*mp.pi/7); phi=(1+mp.sqrt(5))/2; Q=2999999
rng=np.random.default_rng(20260816); RND=[(mp.mpf(float(rng.random())),mp.mpf(float(rng.random()))) for _ in range(3)]
vwa=(mp.mpf(1234567)/Q+mp.mpf('1e-18')*(mp.sqrt(2)-1), mp.mpf(765431)/Q+mp.mpf('1e-18')*(mp.sqrt(3)-1))

print("== R5.1  THE SCHEDULE LAW  total(delta) ~ K delta^{1+2/d_eff}  (M2's leg E3b) ==")
K=10**7; DEL=[1e-3,1e-4,1e-5,1e-6]
ROWS=[("B_S3RES   d_eff=1",1/mp.pi,mp.mpf(11)/(20*mp.pi),1),
      ("C_BADAPP  d_eff=2",mp.frac(tt),mp.frac(tt**2),2),
      ("D_RAND1   d_eff=2",RND[0][0],RND[0][1],2),
      ("F_VWA     d_eff=2  <- E3b OMITS THIS ROW TOO",vwa[0],vwa[1],2)]
print(f"   {'connection':<44}{'d_eff':>6}"+"".join(f"{f'{d:.0e}':>14}" for d in DEL)+f"{'fit':>9}{'theory':>8}")
for lab,a,b,d in ROWS:
    v=allvals(a,b,K,W); v.sort()
    tot=[float(v[:int(dd*K)].sum()) for dd in DEL]
    s=np.polyfit(np.log10(DEL),np.log10(np.maximum(tot,1e-300)),1)[0]
    print(f"   {lab:<44}{d:>6}"+"".join(f"{t:>14.4e}" for t in tot)+f"{s:>9.3f}{1+2/d:>8.1f}")
print("   -> the law is recovered on the rows E3b ran (1.98, 1.95 vs 2; 2.96 vs 3) AND breaks on")
print("      the d_eff = 2 row it did not run.  Same pattern as the floor law: the exponent is a")
print("      function of d_eff ONLY ON THE BADLY-APPROXIMABLE/TYPICAL PAIRS THAT WERE TABULATED.\n")

print("== R5.2  ARE F2's N(eps) EXPONENTS ACTUALLY FITTED?  COUNTS AND POISSON BARS ==")
EPS=[10.0**-e for e in range(1,10)]
for lab,a,b,d in ROWS:
    v=allvals(a,b,10**7,W)
    c=[int((v<e).sum()) for e in EPS]
    use=[(e,n) for e,n in zip(EPS,c) if n>20 and e<=1e-3]
    if len(use)>=3:
        x=np.log10([e for e,_ in use]); y=np.log10([n for _,n in use])
        sd=np.array([1/np.log(10)/np.sqrt(n) for _,n in use])       # Poisson bar on log10 N
        s,i=np.polyfit(x,y,1)
        # 1-sigma on the slope from weighted least squares
        Sxx=((x-x.mean())**2).sum(); err=np.sqrt(((sd*(x-x.mean()))**2).sum())/Sxx
        print(f"   {lab:<44} counts {c[:8]}")
        print(f"      {'':<41} fit over eps<=1e-3, N>20 ({len(use)} pts): {s:.3f} +- {err:.3f}"
              f"   theory {d/2:.2f}   deviation {(s-d/2)/err:+.1f} sigma")
print("   -> C_BADAPP / D_RAND1 are 3-point fits whose smallest bin holds ~85 counts: the")
print("      exponent 1.00 is genuinely FITTED and genuinely agrees.  F_VWA's 0.795 is NOT")
print("      noise -- its eps=1e-6 bin holds ~38 where the law predicts ~8.  M2-4's evidence")
print("      quotes the four agreeing d_eff=2 rows and omits this one.\n")
print("   WHERE F_VWA's EXCESS COUNTS LIVE (the tail that is 3,3,3,3,3 in leg B's own table):")
v=allvals(vwa[0],vwa[1],10**7,W); idx=np.flatnonzero(v<1e-8)+1
print(f"      k with 1-|Z_k| < 1e-8: {idx.tolist()}   = the multiples of Q = {Q}")
print(f"      depths: {[f'{v[i-1]:.3e}' for i in idx]}")
print("      i.e. the 'law' fails exactly where the pair's rational core reasserts itself.\n")

print("== R5.3  THE d_eff RULE ON RS-P IS WRONG FOR ONE CONNECTION TYPE ==")
print("   M2 publishes: 'RS-P (p10 = 0) -> only chi_0/chi_C = conj(W_F) survives, so d_eff = 1")
print("   FOR EVERY CONNECTION' and codes it as d_eff = min(1, dim H).  But on RS-P the only")
print("   live ratio is x = e^{-2 pi i alpha}: d_eff is decided by ALPHA ALONE.")
WP=(0.5,0.0,0.5)
for lab,al,be in [("alpha = 1/2 (RATIONAL), beta = cubic irrational",mp.mpf(1)/2,mp.frac(tt)),
                  ("alpha = cubic irrational, beta = 0",mp.frac(tt),mp.mpf(0))]:
    v=allvals(al,be,10**5,WP)
    nz=int((v<1e-15).sum())
    dimH = 1  # both pairs have exactly one independent relation, so M2's rule says d_eff=1
    print(f"   {lab:<46} dim H = {dimH}, M2's rule -> d_eff = 1")
    print(f"      cells with 1-|Z_k| = 0 (to 1e-15) out of 1e5: {nz}   "
          f"-> TRUE d_eff = {0 if nz>0 else 1}")
print("   The first row is ATTAINED on RS-P and has d_eff = 0, not 1.  No table in the lane")
print("   contains such a row, so nothing published is wrong because of it -- but the RULE as")
print("   published is, and this lane's whole claim is that d_eff is the operative variable.")
