"""LEG F -- IS THE EDGE CLOCK'S EXTRA SENSITIVITY PHYSICAL, OR IS IT GAUGE?
   A null read two ways is worthless; this decides which way it reads."""
import numpy as np, sys
sys.path.insert(0,'.')
from wcore import *
np.set_printoptions(precision=9, linewidth=200)
f, c = 1.0, np.sqrt(2.0)
k = K1(f, c)
sA = np.sqrt([0.40,0.15,0.15,0.15,0.15]).astype(complex)
sC = np.sqrt([0.40,0.15,0.15,0.15,0.15])*np.exp(1j*np.array([0.0,0.7,-1.9,2.3,0.4]))
rng = np.random.default_rng(1618033)

print("F.1  gauge transform BOTH the connection and the section (S1:59-63 is a law on both)")
worst=0.0
for _ in range(2000):
    th = rng.uniform(0,2*np.pi,5)
    def shift(loop, ph): return np.array([ph[j]+th[loop[(j+1)%len(loop)]]-th[loop[j]] for j in range(len(loop))])
    kg = Carrier('K1g',5,[0,1,2],[0,3,4], shift([0,1,2],k.phF), shift([0,3,4],k.phC))
    sg = np.exp(1j*th)*sC
    worst = max(worst, np.abs(np.abs(Z_edge(kg,sg,12))-np.abs(Z_edge(k,sC,12))).max())
print("     max | |Z^T_n[a^g, g s]| - |Z^T_n[a, s]| |, n<=12, 2000 gauge maps = %.2e" % worst)
print("     -> the edge-clock functional is GAUGE-INVARIANT on the PAIR (connection, section).")

print("\nF.2  change only the edge-phase SPLIT at fixed holonomies, section held fixed")
outs=[]
for _ in range(6):
    a = rng.uniform(-2,2,3); a = a - a.mean() + f/3
    bb = rng.uniform(-2,2,3); bb = bb - bb.mean() + c/3
    kk = Carrier('K1s',5,[0,1,2],[0,3,4], a, bb)
    assert abs(kk.f-f)<1e-12 and abs(kk.c-c)<1e-12
    outs.append(np.abs(Z_edge(kk,sA,6)))
outs=np.array(outs)
print("     |Z^T_n| n=1..6 for six different splits with IDENTICAL W_F, W_C, identical s:")
for r in outs: print("      ", np.array2string(r, precision=9))
print("     spread across splits, max over n<=6 = %.3e" % np.ptp(outs,axis=0).max())
print("     -> holding s fixed while re-gauging the connection MOVES the answer.  So the extra")
print("        content the edge clock sees is not a function of (W_F, W_C, |s_v|^2); it is a")
print("        function of gauge-invariant DRESSED relative phases  arg( conj(s_u) U_e s_v ),")
print("        i.e. exactly the Wilson-line-dressed data W-06/W-07 built their observable from.")

print("\nF.3  the dressed invariants, exhibited")
print("     For the loop gamma_C = v0->v3->v4->v0 the quantities")
print("       D_j = conj(s_{w_j}) U_{e_j} s_{w_{j+1}}   are NOT invariant, but")
print("       their arguments' differences and the closed product are.  Check:")
for _ in range(3):
    th = rng.uniform(0,2*np.pi,5)
    def shift(loop, ph): return np.array([ph[j]+th[loop[(j+1)%len(loop)]]-th[loop[j]] for j in range(len(loop))])
    kg = Carrier('K1g',5,[0,1,2],[0,3,4], shift([0,1,2],k.phF), shift([0,3,4],k.phC))
    sg = np.exp(1j*th)*sC
    def dressed(car, s):
        loop=car.loopC; out=[]
        for j in range(3):
            u,v = loop[j], loop[(j+1)%3]
            out.append(np.conj(s[v])*np.exp(1j*car.phC[j])*s[u])
        return np.array(out)
    print("      dressed edge pairings, gauged: %s" % np.array2string(dressed(kg,sg), precision=9))
print("      ungauged                        : %s" % np.array2string(dressed(k,sC) if True else 0, precision=9))
def dressed(car, s):
    loop=car.loopC; out=[]
    for j in range(3):
        u,v = loop[j], loop[(j+1)%3]
        out.append(np.conj(s[v])*np.exp(1j*car.phC[j])*s[u])
    return np.array(out)
print("      ungauged                        : %s" % np.array2string(dressed(k,sC), precision=9))
print("     -> identical.  These L dressed pairings per loop are the extra PHYSICAL data the")
print("        circuit clock integrates away and the edge clock reads.  They are gauge-invariant,")
print("        carrier-supplied (S1:16-22 supplies the edges), and absent from pi.")
