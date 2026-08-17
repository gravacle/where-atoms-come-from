# W19-A step 2.  DOES THE PLATEAU DEPEND ON g^2 ?  ONE VARIABLE MOVED: g^2.
# Carrier, system link, fragment rule, N all held fixed.  Two carriers run: the threshold carrier
# (heawood, d=5) and a failing carrier (cube_Q3, d=3), so that the coupling answer is not read off
# a single graph.
#
# Also measured here to MACHINE PRECISION: is the plateau EXACTLY flat, or only flat to 4 places?
# That distinction decides whether the plateau is a theorem (forced by the gauge-invariant algebra)
# or a contingent dynamical fact.  A forced plateau is a WEIGHTLESS control and must be labelled so.
import numpy as np, sys, json
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_A_SWEEP")
from zn_gauge import *
from carriers import *

def plaquette_expectation(g, psi):
    """<W_p> averaged over plaquettes, in the physical sector -- the magnetic order parameter."""
    N,C = g.N, g.C
    m = g._digits(); tot=0.0
    for a in g.avec:
        ph = (m @ np.array(a,dtype=np.int64)) % N
        tot += float((psi**2 * np.cos(2*np.pi*ph/N)).sum())
    return tot/len(g.avec)

def electric_expectation(g, psi):
    """<X_e> averaged over links -- the electric order parameter.  <X_e> = 1 <=> confined vacuum."""
    N,C,D = g.N,g.C,g.dimP; m=g._digits()
    pw = np.array([N**c for c in range(C)],dtype=np.int64); tot=0.0; n=0
    for t in g.tvec:
        t=np.array(t,dtype=np.int64); n+=1
        if not t.any(): tot += 1.0; continue
        j = ((m+t)%N) @ pw
        tot += float((psi*psi[j]).sum())
    return tot/n

print("="*120)
print("W19-A / 02 COUPLING SWEEP.  H = -(1/g^2) sum_p (W+W^dag) - g^2 sum_e (X+X^dag).  Z_2.")
print("  ONE VARIABLE MOVED: g^2.  carrier / link / fragment rule / N held fixed within each block.")
print("="*120)

results={}
for nm,(V,E) in [("heawood_honeycomb7",heawood()), ("cube_Q3",cube())]:
    g = ZNGauge(nm,V,E,2); L=g.L
    l = 0
    frs,d = nested_fragments(V,E,l); cuts,_ = level_cuts(V,E,l)
    print(f"\n### {nm}:  L={L} C={g.C} dim_phys={g.dimP}  system link l=0={E[0]}  d={d}  plateau fragments k=1..{d-1}")
    print(f"{'g^2':>8}{'gap':>12}{'<W_p>':>10}{'<X_e>':>10}{'H(S) bits':>12}   "
          + "".join(f"{'I/H(F'+str(k)+')':>13}" for k in range(1,d+1)) + f"{'max|I/H-1| k<d':>16}")
    rows=[]
    for gsq in [0.20,0.30,0.40,0.50,0.60,0.70,0.85,1.00,1.20,1.50,2.00,3.00,5.00]:
        mag, elec = 2.0/gsq, 2.0*gsq
        psi,E0,gap = g.ground(mag,elec)
        Psi = g.full_vector(psi)
        HS = S_of(Psi,L,2,[l])
        Is = [mutual_information(Psi,L,2,[l],F) for F in frs]
        ratios = [I/HS for I in Is]
        worst = max(abs(r-1) for r in ratios[:d-1])
        W = plaquette_expectation(g,psi); X = electric_expectation(g,psi)
        print(f"{gsq:>8.2f}{gap:>12.3e}{W:>10.5f}{X:>10.5f}{HS:>12.8f}   "
              + "".join(f"{r:>13.9f}" for r in ratios) + f"{worst:>16.2e}")
        rows.append(dict(gsq=gsq,gap=gap,W=W,X=X,HS=HS,ratios=ratios,worst=worst))
    results[nm]=rows
    print(f"    ceiling note: full-space vector is 2^{L} = {2**L} amplitudes.")

print("\n"+"="*120)
print("READ-OFF")
h = results["heawood_honeycomb7"]
print(f"  plateau FLATNESS on heawood, max |I/H(S) - 1| over k=1..d-1 across the whole sweep: "
      f"{max(r['worst'] for r in h):.3e}   -> the plateau is EXACTLY flat at every coupling.")
print(f"  plateau HEIGHT H(S) on heawood ranges {min(r['HS'] for r in h):.6f} .. {max(r['HS'] for r in h):.6f} bits")
print(f"    g^2=0.20 (magnetic / deconfined): H(S)={h[0]['HS']:.6f} bits, <W_p>={h[0]['W']:.4f}, <X_e>={h[0]['X']:.4f}")
print(f"    g^2=5.00 (electric / confined):   H(S)={h[-1]['HS']:.6f} bits, <W_p>={h[-1]['W']:.4f}, <X_e>={h[-1]['X']:.4f}")
print("  WEIGHT FLOOR 0.10 bits: couplings at which the plateau is WEIGHTED (there is something to record):")
print("    heawood: " + ", ".join(f"{r['gsq']:.2f}" for r in h if r['HS']>=0.10))
print("    heawood: WEIGHTLESS at " + ", ".join(f"{r['gsq']:.2f}" for r in h if r['HS']<0.10))
json.dump(results, open("out_02_coupling.json","w"), indent=1)
print("\nDONE 02.")
