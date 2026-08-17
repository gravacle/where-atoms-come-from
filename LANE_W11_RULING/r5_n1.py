# RULING LEG 5 -- N1, THE RESULT PROPOSED FOR PUBLICATION.
#   N1: lambda = m(p00 + p10 x + p01 y + p11 xy), the logarithmic Mahler measure.
# Three separate questions, kept separate on purpose:
#   (Q1) is the CARRIER-INDEPENDENCE of lambda a restatement?   (legs R2/R3/R4: yes)
#   (Q2) is the MAHLER IDENTIFICATION a restatement?            (this leg: NO)
#   (Q3) what hypothesis does N1 need stated?                   (this leg: equidistribution)
import numpy as np, rlib
rng = np.random.default_rng(rlib.SEED)

print("== R5.1  m(P) TWO WAYS: Jensen reduction vs a Jensen-FREE 2-D quadrature ==")
def m2d(p, n=3000):
    a,b,c,d = p
    t = 2*np.pi*np.arange(n)/n
    X = np.exp(1j*t)[:,None]; Y = np.exp(1j*t)[None,:]
    return float(np.log(np.abs(a + b*X + c*Y + d*X*Y)).mean())
for p in ([0.0,0.30,0.30,0.40], [0.44,0.25,0.09,0.22], [4/9,2/9,2/9,1/9]):
    print(f"  pi={np.round(p,4)}  Jensen(2^20) = {rlib.mahler4(np.array(p)):.12f}   2-D quad(3000^2) = {m2d(p):.12f}")
print(f"  register N1 value for pi=(0,.3,.3,.4): -0.767507880357   [W-10 / ERRATUM-vs-W-02]")
print(f"  W-10 N-3: B0b SENSE U is log(4/9) = {np.log(4/9):.12f} exactly\n")

print("== R5.2  (Q3) N1 NEEDS AN EQUIDISTRIBUTION HYPOTHESIS, AND THE CORPUS ALREADY KNOWS IT ==")
K = rlib.K1(); PI = np.array([0.0,0.30,0.30,0.40])
w = np.array([0.40,0.15,0.15,0.15,0.15]); s0 = np.sqrt(w)+0j
for (f,c,tag) in ((1.0, 2**0.5, "GENERIC   (S4:603, the corpus's only published generic pair)"),
                  (2.0, 1.1,    "RESONANT  (S3's headline; -11f+20c = 0, W-10 N-4)"),
                  (np.pi, 3*np.pi/2, "ORDER 4   (S1 sec6's own published connection, W-07)")):
    a = rlib.a_generic(K, rng, f, c)
    MF = rlib.Mcirc(K,K.VF,rlib.holon(K.walkF,a)); MC = rlib.Mcirc(K,K.VC,rlib.holon(K.walkC,a))
    r = rlib.rate(MF, MC, s0, 200000)
    print(f"  {tag}\n     circuit rate at N=2e5 = {r:.9f}   m(P) = {rlib.mahler4(PI):.9f}   |diff| = {abs(r-rlib.mahler4(PI)):.2e}")
print("  -> the identification lambda = m(P) is TRUE at the generic connection and FALSE at the")
print("     resonant and finite-order ones, where the orbit averages over a proper subtorus.")
print("     That is not a defect of N1: it is a HYPOTHESIS N1 must carry.  The corpus's own")
print("     ERRATUM AGAINST W-02 (REGISTER:162-175) is exactly this correction, one round early.\n")

print("== R5.3  (Q2) THE MAHLER STEP IS *NOT* A RESTATEMENT.  It survives changing the operator ==")
print("   Under the fibre-wise root D (a FINER clock on the corpus's own kind of operator) the")
print("   polynomial identity holds at EVERY tick with FRACTIONAL winding, and the rate is m(P).")
for C in (rlib.K1(), rlib.B0b()):
    a = rlib.a_generic(C, rng, 1.0, 2**0.5)
    WF, WC = rlib.holon(C.walkF,a), rlib.holon(C.walkC,a)
    DF, DC = rlib.Droot(C,C.VF,WF,C.LF), rlib.Droot(C,C.VC,WC,C.LC)
    base = rng.dirichlet(np.ones(C.nv)); s = np.sqrt(base)+0j; pi = C.pi_of(s)
    worst = 0.0
    for n in range(1,25):
        uu = np.conj(np.exp(1j*np.angle(WF)/C.LF))**n; vv = np.exp(1j*np.angle(WC)/C.LC)**n
        pred = pi[0] + pi[1]*uu + pi[2]*vv + pi[3]*uu*vv
        worst = max(worst, abs(rlib.Z(DF,DC,s,n)-pred))
    rD = rlib.rate(DF,DC,s,200000)
    print(f"  {C.name:4s} max|Z_n^D - (p00+p10 u^(n/L_F)+p01 v^(n/L_C)+p11 ...)| over n<=24 = {worst:.2e}")
    print(f"       rate under D at N=2e5 = {rD:.9f}   m(pi) = {rlib.mahler4(pi):.9f}   |diff| = {abs(rD-rlib.mahler4(pi)):.1e}")
print("  -> N1's polynomial and its Mahler measure are a statement about the pair (pi, characters),")
print("     NOT about M_gamma.  Changing the operator inside the fibre-wise class relabels the")
print("     characters and leaves the identification standing.  N1 SURVIVES READING B INTACT.\n")

print("== R5.4  UNDER AN ADMISSIBLE EDGE TICK THE RATE IS STATE-DEPENDENT.  12 states, not 3 ==")
for C in (rlib.K1(), rlib.B0b()):
    a = rlib.a_generic(C, rng, 1.0, 2**0.5)
    WF,WC = rlib.holon(C.walkF,a), rlib.holon(C.walkC,a)
    MF,MC = rlib.Mcirc(C,C.VF,WF), rlib.Mcirc(C,C.VC,WC)
    TF,TC = rlib.Tedge(C,C.walkF,a), rlib.Tedge(C,C.walkC,a)
    base = rng.dirichlet(np.ones(C.nv)); S = rlib.same_pi_states(C, rng, base, 12)
    pi = C.pi_of(S[0])
    rc = [rlib.rate(MF,MC,s,20000) for s in S]
    re = [rlib.rate(TF,TC,s,20000) for s in S]
    dmin = min(np.linalg.norm(S[i]-S[j]) for i in range(12) for j in range(i+1,12))
    print(f"  {C.name:4s} pi = {np.round(pi,6)}   m(pi) = {rlib.mahler4(pi):.9f}")
    print(f"       ARMS DIFFED: min||s_i - s_j|| = {dmin:.4f}   max|pi_i-pi_0| = {np.max(np.abs(np.array([C.pi_of(s) for s in S])-pi)):.1e}")
    print(f"       CIRCUIT rate  N=2e4, 12 states: spread = {max(rc)-min(rc):.2e}   value = {rc[0]:.9f}")
    print(f"       EDGE    rate  N=2e4, 12 states: spread = {max(re)-min(re):.2e}   range = [{min(re):.6f}, {max(re):.6f}]")
    print(f"       L_F x EDGE rate (per gamma_F circuit): [{C.LF*min(re):.6f}, {C.LF*max(re):.6f}]  vs m(pi) = {rlib.mahler4(pi):.6f}")
print("  -> clause 1 of the registrar's leg D (the edge rate is STATE-DEPENDENT while the circuit")
print("     rate is not) reproduces on 12 states, both carriers, my connection.  It is the strong")
print("     clause and it is not a units question.  Clause 2 ('not m(P)/L') is well posed but")
print("     weaker, and on B0b it is ILL-POSED: the two loops have different lengths, so there is")
print("     no single L to rescale by -- see R1.4's empty list of common circuit times.")
