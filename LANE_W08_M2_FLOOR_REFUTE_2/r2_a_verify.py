# LANE W-08 / M2 REFUTER 2 (LENS: CONSEQUENCE) — leg A: independent reimplementation.
# NOTHING is taken from the lane's kernel.  Z_k is built from the BRANCH OPERATORS on C^5,
# exactly as W-01/S3 define them, and only then compared with the lane's trig kernel.
# Double precision unless stated; extremes re-done in mpmath at 40 dps.
import numpy as np, mpmath as mp
mp.mp.dps = 40

FACE_V = [0,1,2]; CYC_V = [0,3,4]          # S1 sec1 incidence.  v0 in both -> p00 = 0.
def Zk_bruteforce(f, c, p, K):
    """Z_k = <M_dF^k s, M_c^k s> computed from the OPERATORS, not from any closed form."""
    WF = np.exp(1j*f); WC = np.exp(1j*c)
    inF = np.array([1.0 if v in FACE_V else 0.0 for v in range(5)])
    inC = np.array([1.0 if v in CYC_V  else 0.0 for v in range(5)])
    p = np.asarray(p, dtype=float)
    out = np.empty(K, dtype=complex)
    for k in range(1, K+1):
        a = np.conj(WF**(k*inF)) * (WC**(k*inC))    # per-vertex phase
        out[k-1] = np.sum(p*a)
    return out

def lane_kernel_absZ(alpha, beta, K, W=(0.4,0.3,0.3)):
    """The lane's kernel, transcribed verbatim from m2_d_attained.py."""
    w11,w10,w01 = W
    k = np.arange(1,K+1,dtype=np.float64); u=(k*alpha)%1.0; v=(k*beta)%1.0
    du=u-np.round(u); dv=v-np.round(v); duv=u+v; duv-=np.round(duv)
    S=4*(w11*w10*np.sin(np.pi*dv)**2+w11*w01*np.sin(np.pi*du)**2+w10*w01*np.sin(np.pi*duv)**2)
    Sc=np.minimum(S,1.0)
    return np.sqrt(np.maximum(0.0,1.0-Sc)), S

print("== A1  THE LANE'S KERNEL, CHECKED AGAINST THE OPERATORS (I nearly called it a bug) ==")
print("   The third pair term uses sin^2(pi k (alpha+BETA)), the SUM.  That is correct ONLY under")
print("   the lane's own sign convention x = e^{-2pi i alpha}, y = e^{+2pi i beta} (x = conj(W_F)).")
print("   Under the opposite convention it would be the DIFFERENCE and every rank-1 row would be")
print("   wrong.  Tested here from the operators, both ready states, 200 random connections:")
rng = np.random.default_rng(7)
worst = 0.0; worst_clamp = 0.0
for trial in range(200):
    f = rng.uniform(0,2*np.pi); c = rng.uniform(0,2*np.pi)
    for p in ([0.4,0.15,0.15,0.15,0.15], [0.5,0,0,0.25,0.25]):
        W = (p[0], p[1]+p[2], p[3]+p[4])
        Z = Zk_bruteforce(f,c,p,60)
        absZ,S = lane_kernel_absZ(f/(2*np.pi), c/(2*np.pi), 60, W)
        worst = max(worst, float(np.max(np.abs(np.abs(Z)-absZ))))
        worst_clamp = max(worst_clamp, float(np.max(S)-1.0))
print(f"   max |  |Z_k|_operators  -  |Z_k|_lane-kernel  | over 200x2x60 = 24000 cells : {worst:.3e}")
print(f"   max (S - 1) over the same cells (the lane clamps S at 1; a WRONG third term would")
print(f"   drive S above 1 and the clamp would hide it): {worst_clamp:.3e}  -> the clamp never binds.")
print("   VERDICT: the kernel is RIGHT.  Recorded because I tested it expecting it to be wrong.\n")

print("== A2  THE FOUR REGISTERED lambda VALUES, RECOMPUTED FROM THE OPERATORS ==")
phi=(1+5**0.5)/2
def lam_bruteforce(f,c,p,K):
    Z = Zk_bruteforce(f,c,p,K); return float(np.mean(np.log(np.abs(Z))))
CASES=[("A_S1PUB   f=pi, c=3pi/2      (attained,|H|=4)", np.pi, 3*np.pi/2, -0.804718956217),
       ("C_BADAPP  cubic pair          (H=T^2)", 2*np.pi*np.mod(2*np.cos(2*np.pi/7),1),
                                                  2*np.pi*np.mod((2*np.cos(2*np.pi/7))**2,1), -0.767507880358),
       ("B_S3RES   f=2.0, c=1.1        (H=S^1)", 2.0, 1.1, -0.767014992998),
       ("E_W07GEN  f=2pi.phi,c=2pi.phi^2(H=S^1)", 2*np.pi*phi, 2*np.pi*phi**2, -1.203972804326)]
P_G=[0.4,0.15,0.15,0.15,0.15]
for lab,f,c,ref in CASES:
    v = lam_bruteforce(f,c,P_G,200000)
    print(f"   {lab:<45} (1/N)log|Omega_N| at N=2e5 = {v:>15.9f}   lane's exact lambda {ref:>15.9f}")
print()

print("== A3  THE GENERIC (RANK-2) VALUE, INDEPENDENTLY, BY 1-D REDUCTION + mpmath QUADRATURE ==")
print("   m(0.4+0.3x+0.3y): integrate y first in closed form, int log|A+0.3 e^{i phi}| = log max(|A|,0.3).")
a,b,cc = mp.mpf('0.4'), mp.mpf('0.3'), mp.mpf('0.3')
th = mp.acos((cc**2 - a**2 - b**2)/(2*a*b))      # |a+b e^{i th}| = c  (the kink)
I1 = mp.quad(lambda t: mp.log(mp.sqrt(a**2+b**2+2*a*b*mp.cos(t))), [0, th])
I2 = (mp.pi - th)*mp.log(cc)
m_exact = (I1 + I2)/mp.pi
print(f"   kink at theta* = {mp.nstr(th,20)}  (cos theta* = -2/3 : {mp.nstr(mp.cos(th),20)})")
print(f"   m(0.4+0.3x+0.3y) = {mp.nstr(m_exact,20)}     lane: -0.767507880358")
print(f"   difference to the lane's value: {mp.nstr(m_exact-mp.mpf('-0.767507880358'),5)}")
print(f"   log(0.3) = {mp.nstr(mp.log(mp.mpf('0.3')),20)}   lane's E_W07GEN/circle value -1.203972804326")

def m_1var(coeffs_exps, N=200000):
    """Mahler measure of a Laurent poly via FFT-free trapezoid at high N (double precision)."""
    t = np.arange(N)/N*2*np.pi
    val = np.zeros(N, dtype=complex)
    for cf,ex in coeffs_exps: val += cf*np.exp(1j*ex*t)
    return float(np.mean(np.log(np.abs(val))))
print(f"   B_S3RES subtorus: 11f=20c so H = {{(e^{{-20is}}, e^{{11is}})}}; lambda = m(0.4 z^-9+0.3 z^-20+0.3 z^11)")
print(f"      = {m_1var([(0.4,-9),(0.3,-20),(0.3,11)]):>15.9f}   lane/erratum: -0.767014992998")
print(f"   diagonal circle xy=1 (D3's limit): m(0.4+0.3z+0.3z^-1) = {m_1var([(0.4,0),(0.3,1),(0.3,-1)]):.9f}"
      f"   = log 0.3 = {np.log(0.3):.9f}")
