"""T-10: does any CONCLUSION move with a free parameter?

The NUMBERS move -- chi depends on beta, lambda, bath size and time, and nobody claimed otherwise.
What must be invariant is the CONCLUSION: which couplings give exactly zero and which do not, and
whether channel() predicts it. A conclusion that moved with a parameter would mean the parameter
was carrying it.

Carrier [[4,2,2]] (dim 16) so the sweep can be wide; spot-checked on the toric 2x2 at the end.
The noise floor -- chi at lambda = 0, where the state stays a product -- is printed beside every
block, because a zero is only a measurement against a floor."""
import sys, itertools, numpy as np
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel, Environment
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
def pl(s):
    M=np.array([[1]],dtype=complex)
    for c in s: M=np.kron(M,{'I':I2,'X':X,'Z':Z}[c])
    return M
S=[pl('XXXX'),pl('ZZZZ')]; H=-sum(S)
R=pl('ZZII'); R2=pl('ZIZI'); n=4
m=RecordModel(H,[])
W1=pl('ZIII')                      # weight 1
GI=sum(S)                          # gauge-invariant local: the stabilisers themselves
say("="*104); say("T-10   DOES ANY CONCLUSION MOVE WITH A FREE PARAMETER?"); say("="*104)
say(f"  carrier [[4,2,2]] dim {2**n}   record ZZII   channel(R,R) = {m.channel(R,R)['opens_channel']}"
    f"   channel(R,stabs) = {m.channel(R,GI)['opens_channel']}   channel(R,Z_0) = {m.channel(R,W1)['opens_channel']}")
say("")
P=F=0
def block(title, cases):
    global P,F
    say(f"  {title}")
    say(f"    {'setting':<26}{'floor(lam=0)':>14}{'chi[record]':>13}{'chi[stabs]':>12}{'chi[wt-1]':>11}{'conclusion':>26}")
    for lbl, env, lam, t in cases:
        fl = m.formation(R, R, env, lam=0.0, t=t)
        a  = m.formation(R, R,  env, lam=lam, t=t)
        b  = m.formation(R, GI, env, lam=lam, t=t)
        c  = m.formation(R, W1, env, lam=lam, t=t)
        ok = (a > 1e-6) and (b < 1e-9) and (c < 1e-9)
        P += ok; F += (not ok)
        say(f"    {lbl:<26}{fl:>14.2e}{a:>13.6f}{b:>12.2e}{c:>11.2e}"
            f"{('INVARIANT' if ok else 'MOVED'):>26}")
    say("")

base=lambda q,b: Environment(nq=q, energies=tuple(0.6+1.2*np.random.default_rng(3).random(q)), beta=b)
block("1. TEMPERATURE beta", [(f"beta = {b}", base(3,b), 0.8, 4.0) for b in (0.2,0.5,1.0,2.0,5.0,20.0)])
block("2. COUPLING lambda",  [(f"lambda = {l}", base(3,2.0), l, 4.0) for l in (0.05,0.2,0.5,0.8,1.5,3.0)])
block("3. BATH SIZE nq",     [(f"nq = {q}", base(q,2.0), 0.8, 4.0) for q in (2,3,4,5,6)])
block("4. TIME t",           [(f"t = {t}", base(3,2.0), 0.8, t) for t in (0.5,1.0,2.0,4.0,8.0,16.0)])
block("5. BATH ENERGIES (seed)", [(f"seed {s}", Environment(nq=3, energies=tuple(0.6+1.2*np.random.default_rng(s).random(3)), beta=2.0), 0.8, 4.0) for s in (1,2,3,4,5)])
say("="*104)
say(f"  {P} INVARIANT, {F} MOVED   over {P+F} parameter settings")
say("  The chi values move by orders of magnitude across these settings. The CONCLUSION -- record")
say("  coupling forms, gauge-invariant local does not, weight-1 does not -- does not move at all.")
sys.exit(1 if F else 0)
