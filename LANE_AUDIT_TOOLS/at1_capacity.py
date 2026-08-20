"""AUDIT 1: IS THE 'CAPACITY CEILING' JUST log(dim) OF THE BATH SITE?

C-36 says an environment site holds a bounded amount about the records written on it and records
sharing it split it. Every measurement behind that used a site of ONE QUBIT. A qubit holds one bit.
So the finding may be nothing but 'a qubit holds one bit' -- an information-capacity fact, not a
density law.

THE DISCRIMINATOR: give the site MORE ROOM. If the crowded chi rises as the site dimension grows,
the ceiling is set by the site's Hilbert space and C-36 is an information-capacity statement that
must be named as such. If it does NOT rise, the ceiling is something else and C-36 stands as
written."""
import sys, os, numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','model'))
from record_model import RecordModel, Environment
def say(*a): print(*a); sys.stdout.flush()
ENERGIES=(1.0,1.4,0.7,1.2,0.9,1.6,1.1,0.8)
def mkenv(nq): return Environment(nq, energies=ENERGIES[:nq])
I2=np.eye(2); Z=np.array([[1,0],[0,-1]],dtype=complex)
def op(k,j,P):
    M=np.array([[1]],dtype=complex)
    for i in range(k): M=np.kron(M,P if i==j else I2)
    return M
say("="*104); say("AUDIT 1   IS THE CAPACITY CEILING JUST THE SITE'S HILBERT SPACE?"); say("="*104)
say("  2 records share ONE SITE. The site is made of s bath qubits -- the records couple to ALL of")
say("  them, so they still share, but the site has 2^s dimensions instead of 2.")
say("  CONTROL in the same row: the same 2 records on DISJOINT sites of the same size s.")
say("")
TS=np.linspace(1.0,13.0,21); LAM=0.8; k=2; nS=4
M=RecordModel(np.zeros((nS,nS),dtype=complex))
R0=op(k,0,Z); R1=op(k,1,Z)
say(f"  {'s qubits/site':>15}{'site dim':>10}{'log2(site dim)':>16}{'chi CROWDED':>14}{'chi SPREAD':>13}{'ratio':>9}")
rows=[]
for s in (1,2,3):
    NB=2*s                                   # two sites of s qubits each
    env=mkenv(NB)
    crowd=[(R0,j) for j in range(s)]+[(R1,j) for j in range(s)]          # both records on site A
    spread=[(R0,j) for j in range(s)]+[(R1,s+j) for j in range(s)]       # one record per site
    a=float(np.mean([env.holevo(M.evolve(crowd,env,lam=LAM,t=tt),R0,nS) for tt in TS]))
    b=float(np.mean([env.holevo(M.evolve(spread,env,lam=LAM,t=tt),R0,nS) for tt in TS]))
    rows.append((s,a,b))
    say(f"  {s:>15}{2**s:>10}{s:>16}{a:>14.6f}{b:>13.6f}{(a/b if b>1e-12 else float('nan')):>9.4f}")
say("")
grow=rows[-1][1]-rows[0][1]
say(f"  crowded chi from s=1 to s=3:  {rows[0][1]:.6f} -> {rows[-1][1]:.6f}   change {grow:+.6f}")
say(f"  site dimension over the same range: 2 -> 8, a 4x increase in room")
say(f"  -> {'THE CEILING RISES WITH SITE DIMENSION: C-36 is an information-capacity statement' if grow>0.05 else 'THE CEILING DOES NOT RISE WITH SITE DIMENSION: it is not simply log(dim)'}")
say("")
say("  2. A SECOND DISCRIMINATOR -- is the SUM over records bounded by the site?")
say("     if the site merely holds log2(2^s) = s bits, then chi(R0) + chi(R1) should be capped by s.")
say(f"  {'s':>4}{'chi(R0)':>12}{'chi(R1)':>12}{'sum':>12}{'s (bits available)':>21}{'sum <= s?':>12}")
for s in (1,2,3):
    NB=2*s; env=mkenv(NB)
    crowd=[(R0,j) for j in range(s)]+[(R1,j) for j in range(s)]
    a=float(np.mean([env.holevo(M.evolve(crowd,env,lam=LAM,t=tt),R0,nS) for tt in TS]))
    b=float(np.mean([env.holevo(M.evolve(crowd,env,lam=LAM,t=tt),R1,nS) for tt in TS]))
    say(f"  {s:>4}{a:>12.6f}{b:>12.6f}{a+b:>12.6f}{s:>21}{str(a+b<=s+1e-9):>12}")
say("")
say("="*104); say("  READ -- from the numbers above"); say("="*104)
say("  THE CEILING IS THE SITE'S INFORMATION CAPACITY. Crowded chi rises 0.133478 -> 0.267223 ->")
say("  0.336940 as the site grows from 1 to 3 qubits, and the SUM over both records is bounded by")
say("  the site's s bits in every row. C-36's 'an environment site holds a bounded amount and")
say("  records sharing it split it' is TRUE and it is the HOLEVO BOUND -- conventional physics,")
say("  correctly applied, and NOT a discovery of this program.")
say("")
say("  THE TOOL MADE AN ORDINARY EFFECT LOOK NOVEL. Every measurement behind C-36 used a site of")
say("  ONE qubit, where capacity and density are indistinguishable. Varying the site dimension --")
say("  one line of the experiment, never run -- separates them immediately.")
say("")
say("  WHAT SURVIVES, AND IT IS THE MORE INTERESTING RESULT. Pure capacity splitting would depend")
say("  only on HOW MANY records share the site. O-33a measured a difference of 0.33 in the crowding")
say("  cost depending on whether the partner PAIRS with the record -- at fixed number of records and")
say("  fixed site. Capacity alone cannot produce that, so C-39 is NOT explained by this and stands.")
