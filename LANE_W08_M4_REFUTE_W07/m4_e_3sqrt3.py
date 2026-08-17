# W-08 / M4 leg E — IS  3*sqrt(3)/10 = 0.5196152422706632  REALLY UNREACHABLE ON S1's CONNECTION?
#
# W-07 sec5 / register: "a factor sqrt(3) requires an element of order 3, and <W_F,W_C> = Z_4
# contains none.  Unresolvable while the code does not exist."
# Two independent refutations, and the second one is the one that bites.
#
# ISOLATION LEDGER.  Held fixed throughout: carrier K1, S1's PUBLISHED connection
# a = (pi/3,pi/3,pi/3, pi/2,pi/2,pi/2) — the exact bytes of S1 sec6 — and the dressed sesquilinear
# form conj(t_u) t_v.  Moved: E1 moves the READY STATE only.  E2 moves the TRANSPORT SCHEDULE only
# (whole-circuit M_gamma -> edge-by-edge T, COR-F's object), with state and connection fixed.
# Double precision for the searches; every exhibited hit is re-checked in exact closed form.
import numpy as np
from fractions import Fraction

TARGET = 3*np.sqrt(3)/10
print(f"  TARGET  3*sqrt(3)/10 = {TARGET:.16f}\n")

FACE_V = {0,1,2}; CYC_V = {0,3,4}; TREE = {1:(0,),2:(0,1),3:(3,),4:(3,4)}
a_pub = np.array([np.pi/3]*3+[np.pi/2]*3)
WF = np.exp(1j*sum(a_pub[:3])); WC = np.exp(1j*sum(a_pub[3:]))

print("== E1  THE SEPARATION IS  amp * |group element - 1|.  amp IS CONTINUOUS AND FREE. ==")
print("  |t_v| = |s_v| exactly (dressing has unit modulus), so amp = |s_u||s_v| ranges over")
print("  (0, 1/2] for a normalised state.  The quantised factor takes values {0, sqrt2, 2}.")
print("  Attainable separations at k=1 are therefore the FULL INTERVAL (0, 1].  3*sqrt(3)/10 is in it.")
for factor, name in [(2.0, "|W_F^-1 - 1| = 2   (pairs (0,3),(0,4),(3,0),(4,0))"),
                     (np.sqrt(2), "|rho - 1| = sqrt2  (the 12 pairs of leg D)")]:
    amp_needed = TARGET/factor
    a2 = 0.9; b2 = amp_needed**2/a2                    # |s_u|^2 = a2, |s_v|^2 = b2
    print(f"    factor {name}")
    print(f"      needs amp = {amp_needed:.16f}   e.g. |s_u|^2 = {a2}, |s_v|^2 = {b2:.16f}  (sum {a2+b2:.4f} <= 1)")

print("\n  EXHIBITED, EXACTLY.  Pair (u,v) = (0,3), factor 2, |s_0|^2 = 3/4, |s_3|^2 = 9/100:")
s = np.zeros(5, dtype=complex)
s[0] = np.sqrt(3.0/4.0); s[3] = np.sqrt(9.0/100.0)
s[1] = np.sqrt(1 - 0.75 - 0.09)                      # remainder parked on v1 (in F only)
s /= np.linalg.norm(s)
print(f"      |s|^2 by vertex = {np.round(np.abs(s)**2,12)}   (sums to {float((np.abs(s)**2).sum()):.12f})")
def dress(s,a):
    u=np.exp(1j*np.asarray(a)); t=np.array(s,dtype=complex)
    for v,p in TREE.items():
        w=1.0+0j
        for e in p: w*=u[e]
        t[v]=s[v]/w
    return t
def sep(s,a,u,v,kk=1):
    t=dress(s,a); dF=(v in FACE_V)-(u in FACE_V); dC=(v in CYC_V)-(u in CYC_V)
    return abs(np.conj(t[u])*t[v])*abs(WF**(dF*kk)-WC**(dC*kk))
got = sep(s,a_pub,0,3)
print(f"      |A_03[M_dF s] - A_03[M_c s]| = {got:.16f}")
print(f"      target                        = {TARGET:.16f}      |difference| = {abs(got-TARGET):.3e}")
print("      closed form: 2 * sqrt(3)/2 * 3/10 = 3*sqrt(3)/10.  EXACT, on S1's published connection,")
print("      with the observable W-07 itself reconstructed.  The impossibility claim is false as stated.\n")

print("== E2  AND THE FACTOR sqrt(3) ITSELF IS AVAILABLE — VIA COR-F, THE CORRECTION W-07 CARRIES ==")
print("  COR-F (S3 audit :794, read and listed as load-bearing in W-07 sec0): loop transport is NOT")
print("  the whole-circuit scalar; the EDGE-transport operator T has T^3 = W-01's operator.")
print("  On S1's PUBLISHED connection the edge holonomies are:")
print(f"    e1,e2,e3 : exp(i pi/3)  order {12//np.gcd(2,12)}", end="")
u1 = np.exp(1j*np.pi/3); u4 = np.exp(1j*np.pi/2)
def ordr(z):
    w=z
    for n in range(1,200):
        if abs(w-1)<1e-12: return n
        w*=z
print(f"  -> ord = {ordr(u1)}")
print(f"    e4,e5,e6 : exp(i pi/2)  -> ord = {ordr(u4)}")
grp = sorted({round(np.angle(u1**i*u4**j)/(np.pi/6)) % 12 for i in range(12) for j in range(12)})
print(f"    <U_e : e in E> = <exp(i pi/6)> = Z_12   (angles at multiples of pi/6: {grp})")
print("    Z_12 CONTAINS ELEMENTS OF ORDER 3.  W-07's premise 'the group is Z_4' is the")
print("    whole-circuit stipulation COR-F corrects, applied as if it were the whole story.")
print("    |zeta^j - 1| for zeta = exp(2 pi i/12):")
for j in range(1,7):
    z = np.exp(2j*np.pi*j/12)
    nm = {1:"", 2:"= 1", 3:"= sqrt2", 4:"= SQRT(3)  <-- order-3 element, on S1's own edges", 6:"= 2"}.get(j,"")
    print(f"      j={j:>2}  ord={ordr(z):>2}  |zeta^j - 1| = {abs(z-1):.12f}   {nm}")
print()
print("  TWO EDGES OF THE FACE, transported one at a time (T^2, not T^3 = M_gamma), accumulate")
print("  exp(2 i pi/3): an element of order 3, of modulus-difference EXACTLY sqrt(3).")
print(f"  With amp = 3/10 the separation is EXACTLY 3*sqrt(3)/10 = {3*np.sqrt(3)/10:.16f}")
amp_check = 0.3
print(f"    check: 3/10 * |exp(2 i pi/3) - 1| = {amp_check*abs(np.exp(2j*np.pi/3)-1):.16f}")
print("  amp = |s_u||s_v| = 3/10 is realised by, e.g., |s_u|^2 = |s_v|^2 = 3/10 — a normalised")
print("  state with two equally weighted vertices at weight 3/10.  Nothing exotic is required.\n")

print("== E3  WHAT SURVIVES OF W-07 sec5 ==")
print("  SURVIVES: W-06 left no code, and 3*sqrt(3)/10 does not arise from the SPECIFIC pair")
print("            (whole-circuit transport, W-07's reconstructed A_uv, W-07's own random state).")
print("  FALSE  : 'a factor sqrt(3) requires an element of order 3 and Z_4 has none' — as an")
print("            impossibility claim it fails twice: amp is a free continuous parameter (E1),")
print("            and the edge-transport group on the same published connection is Z_12,")
print("            which has order-3 elements and yields sqrt(3) exactly (E2).")
print("  The register states the impossibility flatly.  It is not one.")
