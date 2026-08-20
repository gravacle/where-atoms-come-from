"""V1 -- IS ARM A A MEASUREMENT, OR AN ALGEBRAIC IDENTITY?

ARM A is the finding's headline evidence: "record count x16 at fixed G/Z leaves phi_fix
unchanged to every digit".  EVERY step of ARM A is of the form  G -> G x Z_2.

CLAIM UNDER TEST (mine, adversarial): for ANY finite G and ANY abelian A, the minimal-torus
carrier of G x A is  H(GxA) = H(G) (x) I_{|A|^2}  with the transport action
rho_{GxA} = rho_G (x) 1^{(+)|A|^2}.  Hence every eigenspace dimension and every isotypic
multiplicity is multiplied by EXACTLY |A|^2, and phi_fix = sum m^2 / sum d^2 is invariant
IDENTICALLY -- no experiment could have returned anything else.  If that is right, ARM A
carries zero discriminating information about record count.

TESTED, NOT ARGUED:
  (a) eigenvalue multiplicities of H(GxA) vs |A|^2 * (those of H(G))          [matrices]
  (b) transport character chi'_{E}((h,alpha)) vs |A|^2 * chi_E(h)             [matrices]
  (c) isotypic multiplicities m'_rho vs |A|^2 * m_rho                          [characters]
  (d) phi(GxA) - phi(G) == 0 exactly
POSITIVE CONTROL IN THE SAME TABLE (D-15): a step that is NOT a product with an abelian
group -- D_4 -> D_8, D_8 -> D_16, D_4 -> ES_2^(1+4) -- where the same computation MUST
return a non-zero difference, otherwise the comparison discriminates nothing.
"""
import sys, numpy as np
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_A_TRANSPORT"
sys.path.insert(0, LANE)
import glib
from carriers import census, isotypic, phi, minimal_torus, eigblocks
def say(*a): print(*a); sys.stdout.flush()

Z2 = glib.cyclic(2)
def named(n):
    for G in glib.ladder(64):
        if G.name == n: return G
    raise KeyError(n)

say("="*128)
say("V1   ARM A UNDER TEST: IS 'phi_fix UNCHANGED WHEN THE RECORD COUNT x4' A RESULT OR AN IDENTITY?")
say("="*128)
say("")
say("  (A) MATRIX TEST of  H(G x Z_2) = H(G) (x) I_4  and  chi'((h,a)) = 4 chi(h).")
say(f"  {'base G':<12}{'dim(G)':>8}{'dim(GxZ_2)':>12}{'mult(G)':>22}{'mult(GxZ_2)':>26}"
    f"{'mult ratio':>12}{'max|chi_prime-4chi|':>21}")
for base in ("D_4", "Q_8", "D_8", "D_4xZ_2"):
    G = named(base); GA = glib.direct(G, Z2, base + "xZ_2_v")
    ceG, ceA = census(G), census(GA)
    mG = [ceG['dims'][v] for v in (-2,-1,0)]
    mA = [ceA['dims'][v] for v in (-2,-1,0)]
    ratio = sorted({(mA[i]/mG[i]) for i in range(3) if mG[i] > 0})
    # transport character comparison, exact closed form from census
    err = 0.0
    for v in (-2,-1,0):
        for i, g in enumerate(GA.el):
            hidx = G.idx[g[0]]
            err = max(err, abs(ceA['chis'][v][i] - 4.0*ceG['chis'][v][hidx]))
    say(f"  {base:<12}{G.n**2:>8}{GA.n**2:>12}{str(mG):>22}{str(mA):>26}{str(ratio):>12}{err:>21.3e}")

say("")
say("  numeric confirmation on actual matrices (D_4 only, dim 64 -> 256):")
G = named("D_4"); GA = glib.direct(G, Z2, "D_4xZ_2_v")
H1, p1, D1 = minimal_torus(G); H2, p2, D2 = minimal_torus(GA)
w1 = np.round(np.linalg.eigvalsh(H1), 9); w2 = np.round(np.linalg.eigvalsh(H2), 9)
from collections import Counter
c1, c2 = Counter(w1.tolist()), Counter(w2.tolist())
say(f"    spectrum(H(D_4))       = {sorted(c1.items())}")
say(f"    spectrum(H(D_4xZ_2))   = {sorted(c2.items())}")
say(f"    every multiplicity exactly x4 : {all(abs(c2[k]-4*c1[k])==0 for k in c1) and set(c1)==set(c2)}")

say("")
say("  (B) ISOTYPIC MULTIPLICITIES AND phi.   'x|A|^2 exact' = every m_rho multiplied by 4.")
say(f"  {'base G':<12}{'records G':>11}{'records GxZ_2':>15}{'phi(G)':>12}{'phi(GxZ_2)':>14}"
    f"{'phi difference':>17}{'m_rho all x4':>14}")
for base in ("D_4", "Q_8", "D_8", "D_16", "D_4xZ_2", "M_4(2)", "Pauli16"):
    G = named(base); GA = glib.direct(G, Z2, base+"xZ_2_v")
    ceG, ceA = census(G), census(GA)
    isoG, isoA = isotypic(G, ceG), isotypic(GA, ceA)
    pG = phi(isoG, ceG['dims'])[2]; pA = phi(isoA, ceA['dims'])[2]
    ok = True
    for v in (-2,-1,0):
        a = sorted([m for _, m in isoG[v]]); b = sorted([m for _, m in isoA[v]])
        if len(b) != 2*len(a) and len(b) != len(a): ok = ok  # irrep count may double
        if sorted(4*x for x in a) != sorted(b): ok = False
    say(f"  {base:<12}{ceG['dims'][-2]:>11}{ceA['dims'][-2]:>15}{pG:>12.10f}{pA:>14.10f}"
        f"{pA-pG:>17.2e}{str(ok):>14}")

say("")
say("  (C) POSITIVE CONTROL IN THE SAME TABLE -- steps that are NOT 'x abelian'.  The same")
say("      quantity MUST move here, or the ARM-A zero discriminates nothing.")
say(f"  {'step':<26}{'records':>18}{'phi before':>13}{'phi after':>12}{'phi difference':>17}{'|G/Z| before->after':>21}")
pairs = [("D_4","D_8"), ("D_8","D_16"), ("D_4","ES_2^(1+4)"), ("D_8xZ_2","ES_2^(1+4)"),
         ("D_4","M_4(2)"), ("M_4(2)","M_5(2)")]
for a, b in pairs:
    Ga, Gb = named(a), named(b)
    ca, cb = census(Ga), census(Gb)
    pa = phi(isotypic(Ga, ca), ca['dims'])[2]; pb = phi(isotypic(Gb, cb), cb['dims'])[2]
    say(f"  {a+' -> '+b:<26}{str(ca['dims'][-2])+' -> '+str(cb['dims'][-2]):>18}{pa:>13.4f}{pb:>12.4f}"
        f"{pb-pa:>17.2e}{str(Ga.n//ca['Z'])+' -> '+str(Gb.n//cb['Z']):>21}")

say("")
say("="*128); say("  READ -- filled from the numbers above"); say("="*128)
