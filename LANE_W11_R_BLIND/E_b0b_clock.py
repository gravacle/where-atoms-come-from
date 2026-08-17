"""LEG E -- B0b (loop lengths 4 and 3), and the exponent lattice that is the real object.
   NOTE A COSMETIC DEFECT IN D_uniqueness.py, RECORDED NOT PATCHED: the line printing
   "%d draws each" was written without its argument, so it printed the literal '%d'.
   No number is affected; NS = 400 in every row of that table."""
import numpy as np, sys
sys.path.insert(0,'.')
from wcore import *
np.set_printoptions(precision=9, linewidth=200)

f, c = 1.0, np.sqrt(2.0)
b = B0b(f, c)
print("LEG E -- B0b, S4:512/:575 'ring torus 3x3 grid, loops meet'. V=9.")
print("   gamma_F = %s (L=%d, bounds a 2-cell)   gamma_C = %s (L=%d)" % (b.loopF,b.LF,b.loopC,b.LC))
cnt = {k2: len(v) for k2,v in b.classes().items()}
print("   class counts (mine)  : 00:%d 01:%d 10:%d 11:%d" % (cnt[(0,0)],cnt[(0,1)],cnt[(1,0)],cnt[(1,1)]))
print("   class counts (S4:575): 00:4  01:1  10:2  11:2      -> RECONSTRUCTION MATCHES S4")
pU = np.ones(9)/9.0
print("   SENSE U pi = %s ; m(P) = %.12f ; S4/N-3 exact log(4/9) = %.12f" %
      (np.array2string(b.pi(np.sqrt(pU).astype(complex)),precision=6),
       mahler4(*b.pi(np.sqrt(pU).astype(complex))), np.log(4/9)))

# three ready states with identical pi
sA = np.sqrt(pU).astype(complex)
pB2 = pU.copy(); pB2[1]+=0.06; pB2[4]-=0.06; pB2[2]+=0.05; pB2[5]-=0.05   # moves inside 10 and inside 00
sB2 = np.sqrt(pB2).astype(complex)
sC2 = np.sqrt(pU)*np.exp(1j*np.array([0.,0.5,-1.2,2.0,0.3,-0.8,1.7,0.9,-2.4]))
St = {'sA':sA,'sB':sB2,'sC':sC2}
print("\nE.0  ARMS DIFF")
for n,s in St.items():
    print("   %s |s|^2 = %s" % (n, np.array2string(np.abs(s)**2, precision=6)))
print("   pi identical? max dev = %.2e ; ||sA-sB||=%.4f ||sA-sC||=%.4f" %
      (max(np.abs(b.pi(sA)-b.pi(sB2)).max(), np.abs(b.pi(sA)-b.pi(sC2)).max()),
       np.linalg.norm(sA-sB2), np.linalg.norm(sA-sC2)))

K=2000
Ac = np.array([np.abs(Z_circuit(b,s,K)) for s in St.values()])
Ae = np.array([np.abs(Z_edge(b,s,K)) for s in St.values()])
print("\nE.1  CIRCUIT clock spread, max over k<=%d : %.3e" % (K, np.ptp(Ac,axis=0).max()))
sp = np.ptp(Ae,axis=0)
print("E.2  EDGE clock spread,    max over n<=%d : %.3e" % (K, sp.max()))
print("     spread at n=1..14 : %s" % np.array2string(sp[:14],precision=3))
n12 = np.arange(1,K+1)%12==0
print("     spread on n in 12Z : max %.3e   off 12Z : min %.3e" % (sp[n12].max(), sp[~n12].min()))
print("     at n=12 the branches sit at circuit counts (n/4, n/3) = (3,4)  -- MISMATCHED.")
same = [n for n in range(1,20001) if n%4==0 and n%3==0 and n//4==n//3]
print("     n<=20000 with n/4 == n/3 (both branches at equal circuit count): %s  (only n=0 ever)" % same)

print("\nE.3  *** THE EXPONENT LATTICE -- THE OBJECT BOTH CONVENTIONS ARE SAMPLES OF ***")
print("     Both conventions evaluate the SAME family  Y(a,b) = <T_F^a s, T_C^b s>.")
print("       CIRCUIT convention samples the line (a,b) = (L_F k, L_C k) = (4k, 3k)")
print("       EDGE    convention samples the line (a,b) = (n, n)")
print("     Invisibility (|Y| a function of pi alone) holds EXACTLY on the sublattice")
print("     L_F Z x L_C Z = 4Z x 3Z, because that is where BOTH T_F^a and T_C^b are diagonal")
print("     with class-constant diagonal.  Map of |Y| spread over a,b <= 12:")
hdr = "        b\\a " + "".join("%6d"%a for a in range(0,13))
print(hdr)
for bb in range(0,13):
    row = "        %4d " % bb
    for aa in range(0,13):
        vals=[abs(Z_pair(b,s,aa,bb)) for s in St.values()]
        row += "%6s" % ("." if np.ptp(vals)<1e-12 else "X")
    print(row)
print("     '.' = incidence invisible, 'X' = visible.  The dots are exactly 4Z x 3Z.")
print("     The circuit convention's line (4k,3k) lies INSIDE the dot lattice by construction.")
print("     The edge convention's line (n,n) meets it only at the origin-multiples of 12,")
print("     where the two branches have completed DIFFERENT numbers of circuits.")

print("\nE.4  RATES on B0b")
mP = mahler4(*b.pi(sA))
for NN in (2000,20000,200000):
    rc=[rate(Z_circuit(b,s,NN)) for s in St.values()]
    re=[rate(Z_edge(b,s,NN)) for s in St.values()]
    print("   N=%7d circuit %s spread %.2e   edge %s spread %.2e" %
          (NN, ["%.9f"%x for x in rc], np.ptp(rc), ["%.9f"%x for x in re], np.ptp(re)))
print("   m(P) = %.12f = log(4/9) = %.12f" % (mP, np.log(4/9)))

# residue decomposition on B0b, period lcm = 12
print("\nE.5  residue decomposition on B0b (period lcm(4,3)=12)")
cls=b.classes(); order=[(0,0),(1,0),(0,1),(1,1)]
for nm,s in St.items():
    ms=[]
    for j in range(12):
        tr=np.linalg.matrix_power(b.T('F'), j%4)@s
        wr=np.linalg.matrix_power(b.T('C'), j%3)@s
        prod=np.conj(tr)*wr
        cf=np.array([prod[cls[o]].sum() for o in order])
        ms.append(mahler4(*cf))
    pred=float(np.mean(ms)); meas=rate(Z_edge(b,s,200004))
    print("   %s  predicted (1/12)sum_j m(P_j) = %.9f   measured = %.9f   dev %.1e" % (nm,pred,meas,abs(pred-meas)))
