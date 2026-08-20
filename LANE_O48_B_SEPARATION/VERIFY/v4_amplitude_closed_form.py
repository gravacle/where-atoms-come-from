"""VERIFY 4: the closed form the lane did not identify.

DERIVATION (textbook 1D Lindhard):  T(r) = (1/2pi) int dQ Pi(Q) cos(Qr) over the BZ, and at half
filling Pi(Q) has the perfect-nesting log divergence Pi(Q) ~ (1/(2 pi v_F)) ln(Lambda/|Q-2k_F|)
at 2k_F = pi, with v_F = 2 t sin(k_F) = 2t.  Since int ln(Lambda/|u|) cos(ur) du = pi/|r| and
cos((pi+u)r) = (-1)^r cos(ur), the single BZ singularity gives
      T(r) -> (-1)^r / (8 pi t r),        J_eff(r)/g^2 = -8 T(r) -> (-1)^{r+1} / (pi t r).
So C = 1/pi = 0.3183098862 EXACTLY at t = 1, with NO free parameter and NO fit.
This checks that prediction against rings of growing size, then against the lane's own numbers.
"""
import numpy as np
OUT=[]
def P(*a):
    s=" ".join(str(x) for x in a); print(s); OUT.append(s)
def T_pbc(N, epsF=0.0, t=1.0):
    n=np.arange(N); k=2*np.pi*n/N; eps=-2.0*t*np.cos(k)
    occ=np.where(eps<epsF-1e-12)[0]; emp=np.where(eps>epsF+1e-12)[0]
    assert len(occ)+len(emp)==N
    G=np.zeros(N)
    for p in occ:
        G+=np.bincount((p-emp)%N, weights=1.0/(eps[emp]-eps[p]), minlength=N)
    return np.real(np.fft.fft(G))/(N*N)
P("="*110)
P("V4  THE CLOSED FORM:  J_eff(r)/g^2 -> (-1)^{r+1} / (pi r)   [1D RKKY at half filling]")
P("="*110)
P("")
P("[V4a] PEAK OF S(r)*r on rings of growing size, against 1/pi = %.10f" % (1/np.pi))
P(f"{'N':>8} {'argmax r':>9} {'max S(r)*r':>14} {'/(1/pi)':>10} {'S*r at r=N/2048':>17}")
for N in (2050, 8194, 32770, 65538):
    T=T_pbc(N); J=-8.0*T
    rs=np.arange(4, N//8)
    S=((-1)**(rs+1))*J[rs]*rs
    i=int(np.argmax(S))
    rr=max(4,N//2048)
    P(f"{N:>8} {rs[i]:>9} {S[i]:>14.9f} {S[i]*np.pi:>10.6f} "
      f"{(((-1)**(rr+1))*J[rr]*rr):>17.9f}")
P("")
P("[V4b] THE LANE'S OWN REPORTED AMPLITUDES, against the closed form:")
P("      lane 'C = 0.3146' (collapse f(r/m->0))           ->  %.4f %% of 1/pi" % (100*0.3146*np.pi))
P("      lane log-log fits 0.370346 ... 0.375254          ->  %.2f - %.2f %% of 1/pi"
  % (100*0.370346*np.pi, 100*0.375254*np.pi))
P("      closed form 1/pi                                 ->  100.0000 %")
open("v4_amplitude_closed_form.txt","w").write("\n".join(OUT)+"\n")
