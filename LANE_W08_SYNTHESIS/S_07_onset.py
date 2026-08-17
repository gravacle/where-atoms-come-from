# S_07 — THE ONE THING THAT IS NOT UNIFORM: the ONSET SCALE.  lambda < 0 for every non-trivial
# connection, but lambda -> 0 continuously at the trivial connection, and the linear bound of
# M2-1 only bites after K0 ~ 1/|sin(tau/2)| circuits.  ONE VARIABLE: the connection's size.
import numpy as np, math
w=(0.4,0.3,0.3)
def lam_and_onset(f,c,K=10**6):
    k=np.arange(1,K+1,dtype=np.float64)
    z=np.minimum(np.abs(w[0]*np.exp(1j*k*(c-f))+w[1]*np.exp(-1j*k*f)+w[2]*np.exp(1j*k*c)),1.0)
    s=np.log(np.maximum(z,1e-300))
    csum=np.cumsum(1.0-z)
    # first K at which SUM(1-|Z|) exceeds 1 -- "one nat of record written"
    idx=int(np.searchsorted(csum,1.0))+1
    return s.sum()/K, idx
print("   f = c = t on the rank-1 locus u v = 1 (W_C = W_F): tau = the character separation.")
print(f"   {'t':>10} {'lambda (K=1e6)':>16} {'1/|sin(t/2)|':>14} {'K to write 1 nat':>18}")
for t in (1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5):
    L,idx = lam_and_onset(t,t)
    print(f"   {t:>10.0e} {L:>16.9f} {1/abs(math.sin(t/2)):>14.3e} {idx:>18d}")
print()
print("   READ: lambda < 0 at EVERY non-trivial t (durability never fails off G={1}), but the")
print("   number of circuits needed to write one nat of record grows like ~1/t^2.  DURABILITY IS")
print("   A LIMIT STATEMENT WITH NO UNIFORM RATE.  This is the honest residue of the founding")
print("   obstruction and it is a statement about the CONNECTION's size, not about recurrence.")
