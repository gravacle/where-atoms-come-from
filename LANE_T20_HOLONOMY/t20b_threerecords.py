"""T-20b: the diagnostic needs THREE records, and v1 had two defects.

DEFECT 1: v1's formation arm gave 4.409e-01 and 7.867e-01 on BOTH carriers -- identical to four
digits across dim 256 and dim 16, different codes. A quantity that does not notice the carrier is
measuring the bath, not the records. Dropped.
DEFECT 2: v1's writer arm mixed a NOMINATED Xbar with a CONSTRUCTED one and measured their
commutator. 22.70 is the disagreement between two differently-obtained operators. Both writers are
now computed the same way and each is VERIFIED against its own record before use.

AND THE CRITERION NEEDS THREE. 'One record cannot show gravity, two give relation, several make
path comparison possible.' With two records there is one route and nothing to compare. A->{B,C}->D
needs three transports.  Carrier: [[4,2,2]] has k=2. The toric 2x2 has k=2. A 2x3 torus has k=2.
For THREE independent records we need k>=3: the toric 3x3 is k=2 as well, so we use [[8,3,2]] --
three logical qubits on 8 physical -- built here and verified."""
import sys, itertools, numpy as np
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
def pl(s):
    M=np.array([[1]],dtype=complex)
    for c in s: M=np.kron(M,{'I':I2,'X':X,'Z':Z}[c])
    return M
say("="*100); say("T-20b   RECORD-LEVEL HOLONOMY ON THREE RECORDS"); say("="*100)

# [[8,3,2]]: stabilisers XXXXXXXX and ZZZZZZZZ plus ... use the standard k=3 code:
# S = {XXXXXXXX, ZZZZZZZZ, ZZZZIIII, ZZIIZZII, ZIZIZIZI}  -> n=8, 5 stabilisers, k=3, d=2
S=[pl('XXXXXXXX'),pl('ZZZZZZZZ'),pl('ZZZZIIII'),pl('ZZIIZZII'),pl('ZIZIZIZI')]
H=-sum(S); n=8
# logical Z's and their conjugate X's, all COMPUTED the same way and verified below
ZL=[pl('ZZIIIIII'),pl('ZIZIIIII'),pl('ZIIZIIII')]
XL=[pl('XXIIXXII'),pl('XXIIIIXX'),pl('XIXIXIXI')]
m=RecordModel(H,[]); Pg,k=m.ground_space()
say(f"  carrier [[8,3,2]]  dim {2**n}  ground-space dim {k}  (expect 8 = 2^3)")
say("")
say("  VERIFICATION -- every operator checked before use, both writers obtained the same way")
ok=True
for i in range(3):
    a=np.linalg.norm(ZL[i]@H-H@ZL[i]); b=np.linalg.norm(XL[i]@H-H@XL[i])
    c=np.linalg.norm(ZL[i]@XL[i]+XL[i]@ZL[i])
    cross=[np.linalg.norm(ZL[i]@XL[j]-XL[j]@ZL[i]) for j in range(3) if j!=i]
    good = a<1e-9 and b<1e-9 and c<1e-9 and max(cross)<1e-9
    ok &= good
    say(f"    record {i}: [Z_i,H]={a:.1e}  [X_i,H]={b:.1e}  {{Z_i,X_i}}={c:.1e}  "
        f"max|[Z_i,X_j!=i]|={max(cross):.1e}   {'PASS' if good else 'FAIL'}")
if not ok:
    say("\n  VERIFICATION FAILED -- not reporting a holonomy on operators that are not what they claim")
    sys.exit(1)
say("")
say("  THE DIAGNOSTIC:  H_ABCD = (U_A U_B)(U_C U_D)^-1  for the two chains A->B->D and A->C->D,")
say("  with transport = the record's own admissible writer.")
say("")
say(f"  {'chain 1':<16}{'chain 2':<16}{'||H - I||':>13}{'verdict':>10}")
Ident=np.eye(2**n)
names=['A','B','C']
tot=curved=0
for (i,j),(p,q) in itertools.product(itertools.permutations(range(3),2),repeat=2):
    if (i,j)==(p,q): continue
    Hol=(XL[i]@XL[j]) @ np.linalg.inv(XL[p]@XL[q])
    dev=float(np.linalg.norm(Hol-Ident)); tot+=1; curved += dev>1e-9
    if tot<=6:
        say(f"  {names[i]+'->'+names[j]:<16}{names[p]+'->'+names[q]:<16}{dev:>13.3e}{('CURVED' if dev>1e-9 else 'FLAT'):>10}")
say(f"\n  {tot-curved} of {tot} path pairs CLOSE exactly;  {curved} do not")
say("")
say("  AND THE COMMUTATOR TABLE -- does the ORDER of writing matter?")
say(f"    {'':<6}" + "".join(f"{n2:>12}" for n2 in names))
for i in range(3):
    say(f"    {names[i]:<6}" + "".join(f"{np.linalg.norm(XL[i]@XL[j]-XL[j]@XL[i]):>12.1e}" for j in range(3)))
