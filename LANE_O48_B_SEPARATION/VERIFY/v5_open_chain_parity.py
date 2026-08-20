"""VERIFY 5: the lane's collapse [B3] forces r EVEN.  Does the ODD branch collapse to the SAME
function?  Run on the LANE'S OWN open-chain instrument (chi_row), not on my ring."""
import sys, numpy as np
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_B_SEPARATION")
from mediator import chi_row
OUT=[]
def P(*a):
    s=" ".join(str(x) for x in a); print(s); OUT.append(s)
MS=(512,1024,2048,4096)
row={m:chi_row(m,np.ones(m-1),m//2) for m in MS}
def S(m,r): return ((-1)**(r+1))*(-8*row[m][m//2+r])
P("="*112)
P("V5  THE LANE'S OWN [B3] COLLAPSE, RE-RUN SEPARATELY FOR EVEN r AND ODD r (open chain, chi_row)")
P("    lane's loop:  r = int(round(frac*m)); r += (r % 2)   -> r is ALWAYS EVEN")
P("="*112)
P(f"{'r/m':>7} {'parity':>7} " + " ".join(f"{'m='+str(m):>12}" for m in MS)
  + f" {'spread %':>9} {'even vs odd gap %':>18}")
for frac in (0.02,0.05,0.10,0.15,0.20,0.25,0.30):
    got={}
    for par in (0,1):
        vals=[]
        for m in MS:
            r=int(round(frac*m))
            if r%2!=par: r+=1
            vals.append(S(m,r)*r)
        got[par]=vals
        P(f"{frac:>7.3f} {'even' if par==0 else 'odd':>7} " + " ".join(f"{v:>12.6f}" for v in vals)
          + f" {100*(max(vals)-min(vals))/np.mean(vals):>9.4f}"
          + (f" {'':>18}" if par==0 else
             f" {100*abs(np.mean(got[1])-np.mean(got[0]))/np.mean(got[0]):>18.3f}"))
P("")
P("D-15 CONTROL: at small r/m both branches must agree (they share the 1/pi asymptote).")
for frac in (0.002,0.005):
    v=[]
    for par in (0,1):
        vv=[]
        for m in MS:
            r=int(round(frac*m))
            if r<2: r=2
            if r%2!=par: r+=1
            vv.append(S(m,r)*r)
        v.append(np.mean(vv))
    P(f"   r/m={frac}:  even {v[0]:.6f}   odd {v[1]:.6f}   gap {100*abs(v[1]-v[0])/v[0]:.3f} %")
open("v5_open_chain_parity.txt","w").write("\n".join(OUT)+"\n")
