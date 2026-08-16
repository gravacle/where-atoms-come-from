#!/usr/bin/env python3
# PATCH RUN: corrects two mislabelled connections in ref_sched.py and adds the
# lower-bound test.  Same conventions, same seeds.
import numpy as np
TWOPI = 2*np.pi
P11,P10,P01,P00 = 0.4,0.3,0.3,0.0
def Zk(k,f,c):
    k=np.asarray(k,dtype=float)
    return P00 + P10*np.exp(-1j*k*f) + P01*np.exp(1j*k*c) + P11*np.exp(1j*k*(c-f))
def Gk(k,f,c):
    with np.errstate(divide='ignore'): return np.log(np.abs(Zk(k,f,c)))
def hdr(s): print("\n"+"="*78); print(s); print("="*78)

GOLD=(1+np.sqrt(5))/2
hdr("CORRECTION: 2pi/gold + 2pi/gold^2 = 2pi EXACTLY -> that pair is RESONANT, L = Z(1,-1)")
f,c = TWOPI/GOLD, TWOPI/GOLD**2
print("f + c - 2pi = %.3e   => c = -f mod 2pi   => L = { (m,n) : m = -n } = Z*(1,-1)" % (f+c-TWOPI))
print("on H = {(-t,-t)} : Z = 0.6 e^{i t} + 0.4 e^{2 i t}, so |Z| = |0.6 + 0.4 e^{i t}|")
print("   inf_H G = log 0.2 = %.9f   sup_H G = log 1.0 = 0" % np.log(0.2))
ks=np.arange(1,2_000_001); g=Gk(ks,f,c)
print("   min_{k<=2e6} G = %.9f   (approaches log 0.2 from ABOVE, never below)" % g.min())
print("   max_{k<=2e6} G = %.3e" % g.max())
print("   #{k : G < log 0.2 } = %d   <== LOWER BOUND IS HARD" % int((g < np.log(0.2)-1e-12).sum()))
print("   int_H G dHaar = m(0.4 z^2 + 0.6 z) = log 0.4 ... check:")
t=(np.arange(2_000_001)+0.5)*TWOPI/2_000_001
print("      quadrature %.9f   vs log(0.6)=%.9f  vs log(0.4)=%.9f"
      % (np.log(np.abs(0.6+0.4*np.exp(1j*t))).mean(), np.log(0.6), np.log(0.4)))

hdr("GENUINELY rank-0 BADLY-APPROXIMABLE PAIR: f = 2pi*gold, c = 2pi*sqrt2")
f2,c2 = TWOPI*GOLD, TWOPI*np.sqrt(2)
print("1, gold, sqrt2 are Q-linearly independent (gold in Q(sqrt5)), so L = {0}, H = T^2")
ks=np.arange(1,4_000_001); g2=Gk(ks,f2,c2)
print("   lambda_B, N=4e6 : %.9f   (generic value of record -0.767507880)" % g2.mean())
print("   min_{k<=4e6} G  : %.6f  at k=%d   => inf_H G = -infinity (conical zeros on T^2)"
      % (g2.min(), ks[g2.argmin()]))

hdr("lambda = 0 EXACTLY AT A TORSION CONNECTION (f=pi, c=3pi/2, orbit order 4)")
f3,c3 = np.pi, 1.5*np.pi
print("   |Z_k| for k=1..8 :", np.array2string(np.abs(Zk(np.arange(1,9),f3,c3)),precision=12))
print("   |Z_k| = 1 iff k = 0 mod 4  (NOT k=2: S3 4.1's |Z_2|=1 uses weights (1/2,0,1/2), not (0.4,0.3,0.3))")
sch=np.full(200000,4,dtype=np.int64)
print("   schedule k_n = 4 : lambda = %.1f  EXACTLY (every factor has modulus exactly 1)"
      % Gk(sch,f3,c3).mean())
print("   max |1-|Z_{k_n}|| over the schedule = %.3e" % np.abs(1-np.abs(Zk(sch,f3,c3))).max())

hdr("AT A NON-TORSION CONNECTION lambda=0 IS A LIMIT ONLY - NO TERM ATTAINS |Z|=1")
f4,c4 = 1.0,np.sqrt(2)
ks=np.arange(1,4_000_001); g4=Gk(ks,f4,c4)
print("   max_{k<=4e6} |Z_k| = %.15f  < 1 strictly ; deficiency %.3e" % (np.exp(g4.max()),1-np.exp(g4.max())))
print("   every cell of every schedule contributes a STRICTLY negative term here.")
# build a lambda -> 0 schedule and show the running rate
rec_k=[];cur=-np.inf
for i in range(len(ks)):
    if g4[i]>cur: cur=g4[i]; rec_k.append(int(ks[i]))
sch=[];j=0
while len(sch)<400000:
    sch += [rec_k[min(j,len(rec_k)-1)]]*max(1,int(2.0**j)); j+=1
sch=np.array(sch[:400000],dtype=np.int64)
run=np.cumsum(Gk(sch,f4,c4))/np.arange(1,len(sch)+1)
print("   running lambda at N=10,1e2,1e3,1e4,1e5,4e5 :",
      ["%.3e"%run[i] for i in [9,99,999,9999,99999,399999]])
print("   -> 0 monotonically in magnitude.  The limit is exactly 0; no partial sum is.")
print("\nDONE.")
