# W-14 leg D — the two proofs the inventory actually rests on, supplied late.
# Written after the principal asked what "by theorem rather than by exhaustion" means, which
# exposed that leg B was a CONTROL THAT COULD NOT HAVE FAILED, presented as the finding.
import numpy as np
rng=np.random.default_rng(20260822)
print("== D1  PROOF (two lines) — the state enters ONLY through pi, given a diagonal convention ==")
print("  Z_k = <s, Q_k s> = SUM_{u,v} conj(s_u) Q_k[u,v] s_v.")
print("  (a) If Q has ANY off-diagonal entry, the sum contains a cross-vertex product")
print("      conj(s_u) s_v with u != v, which pi does not determine: hold |s| fixed and rotate")
print("      the phase of s_v alone and the term moves. So pi-dependence forces Q DIAGONAL.")
print("  (b) For diagonal Q,  Z_k = SUM_v d_v |s_v|^2 , which is a function of the CLASS SUMS")
print("      iff d_v is constant on each class. QED.")
print("  W-01 DEFINES M_gamma diagonal, so Q_k is diagonal with entry conj(W_F)^{ka} W_C^{kb} at")
print("  class (a,b) -- constant on classes BY CONSTRUCTION. The state side is closed by (a)+(b),")
print("  not by any run.\n")
print("  The proof is CONSTRUCTIVE, so it predicts exactly where a run WOULD see something.")
print("  Leg B could not have failed; this can. One off-diagonal entry, everything else fixed:")
def sweep(offdiag):
    V=5; CLS=[(1,1),(1,0),(1,0),(0,1),(0,1)]
    Q=np.diag(np.exp(1j*np.array([0.3,1.1,1.1,2.0,2.0])))     # class-constant diagonal
    if offdiag: Q[1,2]+=0.4                                    # ONE entry, within a class
    w=np.array([.4,.15,.15,.15,.15]); vals=[]
    for _ in range(200):
        ph=rng.uniform(0,2*np.pi,V); s=np.sqrt(w)*np.exp(1j*ph)
        vals.append(abs(np.vdot(s,Q@s)))
    return max(vals)-min(vals)
print(f"    Q class-constant diagonal        -> spread over 200 phase draws = {sweep(False):.3e}")
print(f"    Q + ONE off-diagonal entry (1,2) -> spread over 200 phase draws = {sweep(True):.3e}")
print("  -> the test has a live failure mode and finds it exactly where (a) says it is.\n")
print("== D2  THE CONNECTION SIDE IS A THEOREM ALREADY (W-12), RESTATED FOR THE INVENTORY ==")
print("  a |-> (<gF,a>, <gC,a>) is a continuous homomorphism of compact CONNECTED groups, so its")
print("  image is a connected closed subgroup of T^2 -- {1}, a circle, or T^2 -- of dimension the")
print("  R-rank of the 2 x E matrix. Two DISTINCT simple cycles are R-dependent only if they share")
print("  support, so the rank is 2 and the image is T^2, with Haar pushforward. Hence any")
print("  absolutely continuous measure on connections stays absolutely continuous on T^2, where")
print("  the resonant set is null. NO SAMPLING IS INVOLVED ANYWHERE IN THAT ARGUMENT.")
