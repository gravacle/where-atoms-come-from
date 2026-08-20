"""V5 -- IS THE NOISE FLOOR ADEQUATE, AND WHAT WOULD IT ACTUALLY HAVE DETECTED?

The lane's floor is measured three ways, and ALL THREE compare configurations that an exact
symmetry maps onto each other (permutation replicas, block-exchange replicas, and a
separate-site control that factorises exactly).  Such a comparison measures ROUND-OFF, not
sensitivity.  The honest question is: how large would a real separation dependence have to be
for this pipeline to see it?  Measured here by modulating the partner's coupling by (1-eps),
which is what any metric fall-off would look like at leading order.

Also: the reported numbers are averages over ONE time window, linspace(1,13,25).  If the window
moves, the 12-digit numbers move.  Then '0.385118611788 at every separation' is a within-window
identity, which is what the symmetry argument already guarantees.
"""
from mpmath import mp, mpf, cos, sin, sqrt, log, tanh
mp.dps = 40
BETA = mpf(2); Z0 = -tanh(BETA)
def bl(K, lam, t):
    r = sqrt(1 + (lam*K)**2); a = lam*K/r; c = mpf(1)/r; th = 2*t*r
    C = cos(th); S = sin(th)
    return (a*c*Z0*(1-C), -a*Z0*S, Z0*C + c*c*Z0*(1-C))
def H2(x):
    if x <= 0 or x >= 1: return mpf(0)
    return -(x*log(x)+(1-x)*log(1-x))/log(2)
def Sv(v): return H2((1+sqrt(v[0]**2+v[1]**2+v[2]**2))/2)
def chi2(lam, w, times):
    """readout A; partner B present with coupling weight w (w = 1 is the lane's case).
       sectors: A = +-1 couples lam*sA, B = +-1 couples lam*w*sB -> effective K is real."""
    acc = mpf(0)
    for t in times:
        vp = [mpf(0)]*3; vm = [mpf(0)]*3
        for sB in (1, -1):
            b = bl(1 + w*sB, lam, t)
            for i in range(3): vp[i] += b[i]/2
            b = bl(-1 + w*sB, lam, t)
            for i in range(3): vm[i] += b[i]/2
        av = [(vp[i]+vm[i])/2 for i in range(3)]
        d = Sv(av) - (Sv(vp)+Sv(vm))/2
        acc += d if d > 0 else mpf(0)
    return acc/len(times)
def chi1(lam, times):
    acc = mpf(0)
    for t in times:
        vp = bl(1, lam, t); vm = bl(-1, lam, t)
        av = [(vp[i]+vm[i])/2 for i in range(3)]
        d = Sv(av) - (Sv(vp)+Sv(vm))/2
        acc += d if d > 0 else mpf(0)
    return acc/len(times)

T_lane = [mpf(1) + (mpf(12)*i)/24 for i in range(25)]
lam = mpf('0.8')
I0 = chi1(lam, T_lane) - chi2(lam, 1, T_lane)
print("="*100)
print("V5(a)  WHAT SIZE OF SEPARATION DEPENDENCE WOULD THE PIPELINE HAVE REGISTERED?")
print("="*100)
print(f"  lane interaction I(w=1) = {mp.nstr(I0,14)}   (lane prints 0.385118611788)")
print(f"  {'eps':>10}{'I(1-eps) - I(1)':>26}{'vs float64 floor 1e-15':>26}")
for k in (3, 6, 9, 12, 15):
    eps = mpf(10)**(-k)
    d = (chi1(lam, T_lane) - chi2(lam, 1-eps, T_lane)) - I0
    print(f"  1e-{k:<8}{mp.nstr(d,10):>26}{('DETECTABLE' if abs(d)>mpf('1e-15') else 'below floor'):>26}")
print("  => the float layer excludes a separation-dependent modulation of the coupling only down")
print("     to ~1e-15 relative.  A gravity-strength one (1e-36) is twenty-one orders below that.")
print("     Only the SYMMETRY argument excludes it absolutely -- and that argument is a statement")
print("     about the carrier's symmetry group, not about records (see V3a).")
print()
print("="*100)
print("V5(b)  THE 12-DIGIT NUMBERS ARE WINDOW ARTEFACTS")
print("="*100)
print(f"  {'time window':>34}{'chi alone':>20}{'interaction':>20}")
for lo, hi, k in ((1,13,25), (1,13,51), (1,26,25), (5,17,25), (0.5,6.5,25)):
    T = [mpf(lo) + (mpf(hi)-mpf(lo))*i/(k-1) for i in range(k)]
    a = chi1(lam, T); c = chi2(lam, 1, T)
    print(f"  {('[%s,%s] x%d'%(lo,hi,k)):>34}{mp.nstr(a,12):>20}{mp.nstr(a-c,12):>20}")
print("  The quantity is window-dependent at the 1e-2 level.  Reporting agreement to 1e-16")
print("  BETWEEN configurations inside one window is therefore a symmetry identity, not evidence")
print("  that a physical constant has been measured to 16 digits.")
