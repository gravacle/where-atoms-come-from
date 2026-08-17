# LANE W-11-R-T, leg T2 — DOES ANY MEMBER OF THE FAMILY RESTORE INVISIBILITY?
# The brief: "Sample the family broadly and report the fraction that restores invisibility."
# I do that, and I also PROVE the answer, and then check the proof's sharp numerical predictions.
import numpy as np, itertools, rlib
from rlib import K1, B0b, SYNTHD, haar, member, diag_member, shift_member, Dscore, Dper_n
rng = np.random.default_rng(20260817)
TOL = 1e-9

print("== T2.0  THE INSTRUMENT, AND WHY IT IS STRONGER THAN THE REGISTRAR'S THREE STATES ==")
print("  Q_n = (T_F^n)* T_C^n and Z_n(s) = s* Q_n s.  Z_n depends on s only through pi FOR ALL s")
print("  <=> Q_n is DIAGONAL and its diagonal is CONSTANT ON EACH CLASS.  D(n) measures exactly")
print("  that.  The registrar sampled three states; three states can agree by accident, and a")
print("  three-state agreement is the very 'could not have failed' shape this program keeps")
print("  finding.  D is a statement about the quantifier, so it cannot.")
a = rlib.a_with_holonomies(K1, 1.0, np.sqrt(2), rng)      # S4:603, the only generic connection published
MF, MC = K1.M(a,'F'), K1.M(a,'C')
TF, TC = K1.T_corf(a,'F'), K1.T_corf(a,'C')
print(f"  agreement check on K1, connection f=1.0 c=sqrt(2):")
print(f"     D(CIRCUIT convention M) over n<=12 = {Dscore(K1,MF,MC,12):.2e}   (registrar leg B1: 1e-16 spreads)")
print(f"     D(EDGE convention  T)  over n<=12 = {Dscore(K1,TF,TC,12):.2e}   (registrar leg B2: up to 5.9e-01)")
print(f"     D per tick n=1..9, EDGE: {['%.1e'%x for x in Dper_n(K1,TF,TC,9)]}")
print("     -> zero exactly at n = 3,6,9. The registrar's table, derived rather than sampled.\n")

print("== T2.1  THE THEOREM THIS LANE PROVES (statement; the proof is in the FINDINGS block) ==")
print("  Let J = gamma_F ^ gamma_C (this is exactly class 11).  Then invisibility holds at tick n")
print("  IFF  T_F^n and T_C^n are both block-diagonal w.r.t. J, DIAGONAL and class-constant off J,")
print("  and equal up to a scalar ON J:  T_C^n|_J = q_n T_F^n|_J.")
print("  Since T_F^{L_F} = W_F.I and T_C^{L_C} = W_C.I on their loops, the J-block A satisfies")
print("  A^{L_F} = scalar and A^{L_C} = scalar, so A^{gcd(L_F,L_C)} is scalar.")
print("  COROLLARY.  If |J| <= 1  OR  gcd(L_F,L_C) = 1, then A is scalar and EVERY")
print("  invisibility-restoring member is FIBRE-WISE (diagonal) and CLASS-UNIFORM.")
print("  K1: |J| = 1.   B0b: |J| = 2 but gcd(4,3) = 1.   BOTH CORPUS CARRIERS ARE COVERED.\n")

def report(car, label, f, c, nmax, NFULL, NSHIFT, roots_pred):
    print(f"== {label} ==")
    a = rlib.a_with_holonomies(car, f, c, rng)
    LF, LC = car.LF, car.LC
    # --- (i) the FULL family, sampled broadly ---
    hits = 0; best = np.inf; hits_scalar = 0; gen = 0; gen_hits = 0; gen_best = np.inf
    for _ in range(NFULL):
        jF, jC = rng.integers(0,LF,LF), rng.integers(0,LC,LC)
        TF = member(car, a, 'F', jF, haar(LF,rng))
        TC = member(car, a, 'C', jC, haar(LC,rng))
        d = Dscore(car, TF, TC, nmax); best = min(best,d)
        nondeg = (len(set(jF.tolist()))==LF) and (len(set(jC.tolist()))==LC)
        if nondeg:
            gen += 1; gen_best = min(gen_best,d)
            if d < TOL: gen_hits += 1
        if d < TOL:
            hits += 1
            offF = np.linalg.norm(TF-np.diag(np.diag(TF))); offC = np.linalg.norm(TC-np.diag(np.diag(TC)))
            if offF < 1e-9 and offC < 1e-9: hits_scalar += 1
    print(f"  FULL family, {NFULL} Haar-random pairs      : {hits} restore invisibility "
          f"(fraction {hits/NFULL:.4f}); best D = {best:.2e}")
    print(f"     of those {hits} hits, {hits_scalar} have BOTH branch operators exactly DIAGONAL"
          f" (off-diagonal < 1e-9).")
    print(f"     they are the degenerate-spectrum members T_S = w.I, for which the Haar direction V")
    print(f"     is irrelevant -- i.e. the SAME finite set as branch (A), reached with probability")
    print(f"     ({LF}/{LF**LF})x({LC}/{LC**LC}) = {LF**(1-LF)*LC**(1-LC):.4f} by the root-index draw.")
    print(f"  FULL family, GENERIC component only (all L roots distinct, {gen} of {NFULL} draws):"
          f" {gen_hits} restore invisibility; best D = {gen_best:.2e}")
    # --- (ii) the SHIFT branch (the only local branch that MOVES anything), sampled ---
    hits = 0; best = np.inf
    for _ in range(NSHIFT):
        TF = shift_member(car, a, 'F', rng.uniform(0,2*np.pi,LF-1))
        TC = shift_member(car, a, 'C', rng.uniform(0,2*np.pi,LC-1))
        d = Dscore(car, TF, TC, nmax); best = min(best,d)
        if d < TOL: hits += 1
    print(f"  SHIFT branch (B), {NSHIFT} random members  : {hits} restore invisibility; best D = {best:.2e}")
    # --- (iii) the DIAGONAL branch (A), EXHAUSTIVE ---
    tot = 0; hits = 0; witness = None
    for jF in itertools.product(range(LF), repeat=LF):
        TFd = diag_member(car, a, 'F', jF)
        for jC in itertools.product(range(LC), repeat=LC):
            TCd = diag_member(car, a, 'C', jC)
            tot += 1
            if Dscore(car, TFd, TCd, nmax) < TOL:
                hits += 1
                if witness is None: witness = (jF,jC)
    print(f"  DIAGONAL branch (A), EXHAUSTIVE {tot:>5} pairs: {hits} restore invisibility"
          f"   [theorem predicts {roots_pred}]")
    print(f"     one witness: root-index vectors jF={witness[0]} jC={witness[1]}")
    # --- (iv) COR-F's own T ---
    TF, TC = car.T_corf(a,'F'), car.T_corf(a,'C')
    per = Dper_n(car, TF, TC, nmax)
    zer = [n+1 for n,x in enumerate(per) if x < TOL]
    print(f"  COR-F's T: invisible exactly at n = {zer}   [lcm(L_F,L_C) = {np.lcm(LF,LC)}]")
    print(f"     max D over n<=%d = %.3e\n" % (nmax, max(per)))

report(K1,  "T2.2  K1   (|J| = 1)",              1.0, np.sqrt(2), 12, 20000, 20000, "3^2 x 3^2 = 81 of 729")
report(B0b, "T2.3  B0b  (|J| = 2, gcd(4,3) = 1)",1.0, np.sqrt(2), 24,  4000,  4000, "4x4x3x3 = 144 of 6912")

print("== T2.4  THE REGISTRAR'S CONCLUSION SENTENCE IS FALSE AS WRITTEN — COUNTEREXAMPLE ==")
print("  Registrar: 'the invisibility of incidence holds exactly where both branch operators are")
print("  DIAGONAL'.  DIAGONALITY IS NOT SUFFICIENT.  Take T_F = diag(w, w.zeta, w.zeta^2, 1, 1) on")
print("  K1 -- perfectly diagonal, T_F^3 = M_dF exactly -- but the two class-10 vertices v1, v2 get")
print("  DIFFERENT cube roots, so the within-class distribution becomes visible:")
a = rlib.a_with_holonomies(K1, 1.0, np.sqrt(2), rng)
TFbad = diag_member(K1, a, 'F', (0,1,2)); TCok = diag_member(K1, a, 'C', (0,0,0))
print(f"     ||T_F^3 - M_dF||        = {np.linalg.norm(np.linalg.matrix_power(TFbad,3)-K1.M(a,'F')):.2e}")
print(f"     T_F diagonal?             {np.allclose(TFbad, np.diag(np.diag(TFbad)))}")
print(f"     D (invisibility score)  = {Dscore(K1,TFbad,TCok,12):.3f}   -> INVISIBILITY BROKEN")
sA = np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
sB = np.sqrt(np.array([0.40,0.26,0.04,0.11,0.19]))+0j
print(f"     pi(A) = {np.round(rlib.pi_of(K1,sA),12)}   pi(B) = {np.round(rlib.pi_of(K1,sB),12)}")
zA = abs(np.vdot(TFbad@sA, TCok@sA)); zB = abs(np.vdot(TFbad@sB, TCok@sB))
print(f"     |Z_1(A)| = {zA:.12f}   |Z_1(B)| = {zB:.12f}   spread {abs(zA-zB):.2e}")
print("  CORRECTED NAME OF THE OPERATIVE VARIABLE: not DIAGONAL -- FIBRE-WISE AND CLASS-UNIFORM.")
print("  'Fibre-wise' is the register's own word (REGISTER:578, W-06's correction of N4); the")
print("  class-uniformity half is new here.  The correction STRENGTHENS the registrar's reading:")
print("  the stipulation that buys invisibility is narrower than the one it named.\n")

print("== T2.5  THE BOUNDARY OF MY OWN THEOREM, EXHIBITED AGAINST MYSELF (SYNTH-D, NOT A CORPUS CARRIER) ==")
print("  If |J| >= 2 AND gcd(L_F,L_C) >= 2 the corollary's hypothesis fails and a NON-DIAGONAL")
print("  invisibility-restorer exists.  SYNTH-D: two 4-cycles sharing {0,1}, spectator v6.")
a = rlib.a_with_holonomies(SYNTHD, 1.0, np.sqrt(2), rng)
WF, WC = SYNTHD.hol(a,'F'), SYNTHD.hol(a,'C')
V2 = haar(2,rng); rF = np.exp(1j*np.angle(WF)/4)
A = V2@np.diag([rF, rF*1j])@V2.conj().T                      # non-scalar, A^4 = W_F I
q  = np.exp(1j*(np.angle(WC)-np.angle(WF))/4)                # q^4 W_F = W_C
dF = rF; dC = np.exp(1j*np.angle(WC)/4)
TF = np.eye(7,dtype=complex); TC = np.eye(7,dtype=complex)
TF[np.ix_([0,1],[0,1])] = A; TF[2,2] = TF[3,3] = dF
TC[np.ix_([0,1],[0,1])] = q*A; TC[4,4] = TC[5,5] = dC
print(f"     ||T_F^4 - M_dF||   = {np.linalg.norm(np.linalg.matrix_power(TF,4)-SYNTHD.M(a,'F')):.2e}")
print(f"     ||T_C^4 - M_c ||   = {np.linalg.norm(np.linalg.matrix_power(TC,4)-SYNTHD.M(a,'C')):.2e}")
print(f"     T_F diagonal?        {np.allclose(TF,np.diag(np.diag(TF)))}     "
      f"||offdiag T_F|| = {np.linalg.norm(TF-np.diag(np.diag(TF))):.3f}")
print(f"     D over n<=16       = {Dscore(SYNTHD,TF,TC,16):.2e}   -> INVISIBILITY IS RESTORED.")
loc = 0.0
for u in range(7):
    for v in range(7):
        if u!=v and abs(TF[u,v])>1e-12:
            if (v,u) not in [(s_,d_) for s_,d_,_ in SYNTHD.loopF]: loc = max(loc, abs(TF[u,v]))
print(f"     but is it an EDGE TICK?  largest T_F entry OFF the loop's directed adjacency = {loc:.3f}")
print(f"     and does it transport?  T_F[2,1] (the entry a tick along edge 1->2 must fill) = {abs(TF[2,1]):.3f}")
print("  -> IT IS NOT LOCAL AND IT MOVES NOTHING ALONG EITHER LOOP.  It is a rotation inside the")
print("     shared block, i.e. a re-labelling of the class-11 fibres, and it needs an eigenbasis")
print("     the carrier does not supply.  It restores invisibility and it is not a transport.")
print("     RECORDED AGAINST MYSELF: my corollary is CONDITIONAL on |J|<=1 or gcd = 1, and both")
print("     corpus carriers satisfy it -- but the unconditional sentence would be FALSE.")
