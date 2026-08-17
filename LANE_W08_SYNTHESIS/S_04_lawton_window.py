# S_04 — IS M1_07's (2584,1597) HEADLINE +4.716e-09 CONVERGED, OR A QUADRATURE WINDOW?
# M1_07:50 caps the np.roots branch at degree 1200, so that row alone goes through a 2^23-point
# trapezoid.  Moving ONE variable — NQ — with everything else identical to M1_07's own code.
import numpy as np
P10, P01, P11 = 0.3, 0.3, 0.4
mP = -0.767507880357775871645874051819   # 50-dp value, independently confirmed by two refuters

def m_quad(mm, nn, NQ):
    terms = [(nn, P10), (-mm, P01), (nn-mm, P11)]
    th = np.arange(NQ) * (2*np.pi/NQ)
    val = np.zeros(NQ, dtype=complex)
    for e, c in terms:
        val += c*np.exp(1j*e*th)
    return float(np.mean(np.log(np.maximum(np.abs(val), 1e-300))))

def m_roots(mm, nn):
    terms = [(nn,P10),(-mm,P01),(nn-mm,P11)]
    shift = -min(e for e,_ in terms); deg = max(e+shift for e,_ in terms)
    coef = np.zeros(deg+1)
    for e,c in terms: coef[e+shift] += c
    nz = np.nonzero(coef)[0]; coef = coef[nz[0]:nz[-1]+1]
    r = np.roots(coef[::-1])
    return float(np.log(abs(coef[-1])) + np.sum(np.log(np.maximum(np.abs(r),1.0))))

print("== S4a  THE ROW M1_07 PUBLISHES AS THE ACCUMULATION HEADLINE:  (2584,1597), degree 4181 ==")
print("   M1_07 prints deviation +4.716e-09 from a single NQ = 2^23 evaluation.")
print(f"   {'NQ':>8} {'lambda':>18} {'deviation from m(P)':>22}")
for e in (22,23,24,25,26,27):
    v = m_quad(2584,1597,1<<e)
    print(f"   2^{e:<6} {v:>18.12f} {v-mP:>+22.4e}")
print("   ROOTS route (degree 4181 -- outside M1_07's own cap, run here anyway):")
v = m_roots(2584,1597); print(f"   {'np.roots':>8} {v:>18.12f} {v-mP:>+22.4e}")

print()
print("== S4b  THE ROW BELOW IT, (610,377), degree 987 -- INSIDE the cap, computed by ROOTS ==")
vr = m_roots(610,377); print(f"   np.roots  {vr:>18.12f}  deviation {vr-mP:+.4e}   (M1_07 prints +5.251e-07)")
for e in (22,23,24,25,26):
    v = m_quad(610,377,1<<e); print(f"   quad 2^{e:<3} {v:>18.12f}  deviation {v-mP:+.4e}")
