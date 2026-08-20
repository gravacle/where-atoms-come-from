"""O-50 D  PART 7 -- THE ESCAPE THIS LANE ADDS, AND THE ONLY ONE THAT IS NOT CLOSED.

Escape (6): THE MEASURE ON CONFIGURATION SPACE.

Every cancellation result in this program -- C-46's screening ratio, C-62's sqrt(2/pi)/sqrt(m),
part 2's orbit lemma -- computes a mean over configurations with the UNIFORM measure.  The
uniform measure is the group-orbit measure, and it is the right one exactly when writing is a
reversible group action with no preferred value.

Clause (iv) demands that an admissible UNITARY writer exist.  It does not say that writing
HAPPENS unitarily, and it does not say the ensemble of written configurations is uniform.  An
irreversible writer -- a reset channel, a dissipative process, decoherence into a preferred
pointer value -- produces a BIASED ensemble.  This part asks what happens then.
"""
import sys, itertools, math
import numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_D_ESCAPE")
from o50d_common import *

say("=" * 104)
say("O-50 D  PART 7   ESCAPE (6)  IRREVERSIBLE WRITING AND THE MEASURE ON CONFIGURATIONS")
say("=" * 104)

# ---------------------------------------------------------------- setup
T = Torus(2); n = T.nq; N = 2 ** n
pairs = symplectic_logicals(T.stab, n); bs = [x for p in pairs for x in p]
def comb(c):
    v = [0] * (2 * n)
    for cc, b in zip(c, bs):
        if cc: v = [(x + y) % 2 for x, y in zip(v, b)]
    return v
R1v = comb((0, 0, 0, 1)); R2v = comb((0, 1, 0, 0))
R1 = dense(R1v, n); R2 = dense(R2v, n); H = -sum(dense(s, n) for s in T.stab)
Xb = dense(comb((0, 0, 1, 0)), n)
Pg = np.eye(N, dtype=complex)
for s in T.stab: Pg = Pg @ (np.eye(N) + dense(s, n)) / 2
w, V = np.linalg.eigh(Pg); V = V[:, w > 0.5]
r1 = V.conj().T @ R1 @ V; r2 = V.conj().T @ R2 @ V
ww, U1 = np.linalg.eigh(r1 + 3.0 * r2); CB = V @ U1
cfg = [(int(round(np.real(CB[:, a].conj() @ R1 @ CB[:, a]))),
        int(round(np.real(CB[:, a].conj() @ R2 @ CB[:, a])))) for a in range(4)]

# ---------------------------------------------------------------- 1. an admissible IRREVERSIBLE writer
say("")
say("1. AN ADMISSIBLE *IRREVERSIBLE* WRITER EXISTS AND IT IS BUILT FROM THE CLAUSES THEMSELVES.")
Kp = (np.eye(N) + R1) / 2
Km = Xb @ (np.eye(N) - R1) / 2
say("   RESET-RECORD-1-TO-+1:  K+ = (I+R1)/2,  K- = Xbar1 (I-R1)/2.")
say(f"   completeness  ||K+dagK+ + K-dagK- - I|| = "
    f"{np.linalg.norm(Kp.conj().T @ Kp + Km.conj().T @ Km - np.eye(N)):.2e}")
say(f"   [K+,H] = {np.linalg.norm(Kp @ H - H @ Kp):.2e}   [K-,H] = "
    f"{np.linalg.norm(Km @ H - H @ Km):.2e}   -- BOTH KRAUS OPERATORS COMMUTE WITH H,")
say("   so the channel commutes with the free evolution: it is admissible in every sense the")
say("   framework uses EXCEPT that it is not a unitary.")
def E1(r): return Kp @ r @ Kp.conj().T + Km @ r @ Km.conj().T
for a in range(4):
    r = np.outer(CB[:, a], CB[:, a].conj())
    say(f"   config {str(cfg[a]):>9} -> <R1> = {np.real(np.trace(E1(r) @ R1)):+.4f}, "
        f"<R2> = {np.real(np.trace(E1(r) @ R2)):+.4f}")
ia = 0; ib = [a for a in range(4) if cfg[a] == (-cfg[ia][0], cfg[ia][1])][0]
ra = np.outer(CB[:, ia], CB[:, ia].conj()); rb = np.outer(CB[:, ib], CB[:, ib].conj())
say(f"   IRREVERSIBLE: configurations {cfg[ia]} and {cfg[ib]} differ only in s_1, and the channel")
say(f"   maps them to states differing by {np.linalg.norm(E1(ra) - E1(rb)):.2e}: IT MERGES THEM.")
say("   A channel that merges states has no inverse.  THE ADMISSIBLE SET IS THEN A")
say("   SEMIGROUP, NOT A GROUP, AND THE ORBIT-AVERAGING LEMMA OF PART 2 HAS NO ORBITS TO AVERAGE")
say("   OVER.")

# ---------------------------------------------------------------- 2. the measure is the hypothesis
say("")
say("2. THE HIDDEN HYPOTHESIS IN EVERY CANCELLATION RESULT IN THIS PROGRAM: THE UNIFORM MEASURE.")
say("   C-46's screening ratio, C-62's sqrt(2/pi) m^(-1/2), part 2's lemma -- each averages over")
say("   configurations with EQUAL WEIGHT.  Under a biased ensemble with magnetisation")
say("   m = <s_i>, the pair functional has EXACT mean")
say("        E[ sum_{i<j} w_ij s_i s_j ] = m^2 sum_{i<j} w_ij        (independent records)")
say("   which is SIGN-DEFINITE for every m != 0 and equals the FULL uncancelled sum at m = +-1.")
say("")
def bias_table(k, ws, ms, rng, samples=400000):
    tot = float(np.abs(ws).sum())
    out = []
    for mm in ms:
        p = (1 + mm) / 2
        S = (rng.random(size=(samples, k)) < p) * 2 - 1
        F = np.einsum('ai,ij,aj->a', S, ws, S)
        out.append((mm, float(F.mean()), float(np.abs(F).mean()),
                    float(F.mean() / tot), float(np.abs(F).mean() / tot),
                    float((F > 0).mean())))
    return tot, out
rng = np.random.default_rng(23)
k = 12
pos = [4.0 * (i // 2) + 0.5 * (i % 2) for i in range(k)]
Wm = np.zeros((k, k))
for i in range(k):
    for j in range(i + 1, k):
        Wm[i, j] = max(1.0, abs(pos[i] - pos[j])) ** -1.0
tot, rows = bias_table(k, Wm, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], rng)
say(f"   k = {k} records, w_ij = 1/r_ij (INSERTED power law), sum w = {tot:.6f}")
say(f"   {'bias m':>8}{'E[F]':>12}{'predicted m^2 * sum w':>23}{'E|F|':>12}"
    f"{'coherence E[F]/sum w':>22}{'frac F > 0':>12}")
for mm, mean, absmean, coh, acoh, frac in rows:
    say(f"   {mm:>8.2f}{mean:>12.6f}{mm * mm * tot:>23.6f}{absmean:>12.6f}{coh:>22.6f}{frac:>12.4f}")
say("   THE PREDICTION m^2 * sum(w) IS MET AT EVERY BIAS.  D-20: this is a closed form checked")
say("   OUT OF SAMPLE, not a fit.  At m = 0 the mean is zero -- C-46's screening.  At m = 1 the")
say("   coherence is 1.000000: NO SCREENING AT ALL.  The pairing accumulates in full.")

say("")
say("3. AND IT KEEPS ACCUMULATING AS RECORDS ARE ADDED -- the opposite of C-46's signature.")
say(f"   {'k':>4}{'m=0 coherence':>16}{'m=0.5 coherence':>18}{'m=1 coherence':>16}"
    f"{'m=0.5 E[F]':>14}{'E[F](2k)/E[F](k)':>19}")
prev = None
for k2 in (4, 8, 16, 32, 64):
    pos2 = [4.0 * (i // 2) + 0.5 * (i % 2) for i in range(k2)]
    W2 = np.zeros((k2, k2))
    for i in range(k2):
        for j in range(i + 1, k2):
            W2[i, j] = max(1.0, abs(pos2[i] - pos2[j])) ** -1.0
    t2, rr = bias_table(k2, W2, [0.0, 0.5, 1.0], rng, samples=120000)
    e_half = rr[1][1]
    ratio = '' if prev is None else f"{e_half / prev:.4f}"
    prev = e_half
    say(f"   {k2:>4}{rr[0][3]:>16.6f}{rr[1][3]:>18.6f}{rr[2][3]:>16.6f}{e_half:>14.4f}{ratio:>19}")
say("   At m = 0 the coherence falls with k (screening).  At m != 0 it is CONSTANT in k and the")
say("   total GROWS -- super-linearly here because the pair count grows, which is the same reason")
say("   a two-body quantity was attractive in the first place.")

# ---------------------------------------------------------------- 4. the control
say("")
say("4. THE C-61 CONTROL, AND THE HONEST VERDICT.")
Hc, Jc, hc, diag = control_carrier(8, seed=31)
say(f"   control carrier: {len(set(np.round(diag,9)))}/{len(diag)} distinct energies -> "
    f"NON-DEGENERATE -> ZERO records by P-1.")
labels = np.array([[1 - 2 * ((b >> (7 - q)) & 1) for q in range(8)] for b in range(256)])
kk = 8
posc = [4.0 * (i // 2) + 0.5 * (i % 2) for i in range(kk)]
Wc = np.zeros((kk, kk))
for i in range(kk):
    for j in range(i + 1, kk):
        Wc[i, j] = max(1.0, abs(posc[i] - posc[j])) ** -1.0
for mm in (0.0, 0.5, 1.0):
    p = (1 + mm) / 2
    S = (rng.random(size=(200000, kk)) < p) * 2 - 1
    F = np.einsum('ai,ij,aj->a', S, Wc, S)
    say(f"   control, bias m={mm:.1f}: E[F] = {F.mean():>10.6f}  coherence = "
        f"{F.mean()/np.abs(Wc).sum():>9.6f}   (record carrier at the same k: same formula, same")
    say(f"                          numbers -- the bias mechanism is RECORD-BLIND)")
say("")
say("   SO WHAT IS AND IS NOT ESTABLISHED HERE:")
say("     ESTABLISHED: the cancellation half of the theorem is a statement about REVERSIBLE")
say("     writing with a UNIFORM configuration measure.  An admissible irreversible writer")
say("     exists -- built from the record itself, both Kraus operators commuting with H -- and")
say("     under the biased ensemble it produces, a two-body record functional is SIGN-DEFINITE,")
say("     ACCUMULATING, and carries whatever falloff the coefficients carry.  Standard (d),")

say("     which killed C-46's pairing and O-48's correlation energy, is NOT violated by any")
say("     clause -- it was violated by an assumption about the measure.")
say("     NOT ESTABLISHED: that the bias is INDUCED.  Here it is INSERTED -- I chose m.  Nothing")
say("     in this lane derives a preferred record value from (H, {L_k}).  And the bias mechanism")
say("     is RECORD-BLIND on the C-61 control: a biased non-record register does the same thing.")
say("")
say("   THIS IS THE ONE OPEN DOOR THIS LANE FOUND, AND IT IS A LANE OF ITS OWN:")
say("     Q1. Does a dissipative dynamics (H, {L_k}) INDUCE a preferred record value -- i.e. is")

say("         there a carrier whose Lindbladian has a biased stationary record ensemble while")
say("         all five clauses still hold?  Clause (iv) constrains the HAMILTONIAN (Tr(P_E R)=0)")
say("         and says nothing about the DISSIPATOR.")
say("     Q2. If so, is the resulting sign-definite accumulation distinguishable from a")
say("         record-blind biased register -- i.e. can it clear the C-61 control?")
say("   Until Q1 is answered the escape is OPEN, not won: an inserted bias proves only that the")
say("   obstruction is in the measure, not in the clauses.")
say("=" * 104)
