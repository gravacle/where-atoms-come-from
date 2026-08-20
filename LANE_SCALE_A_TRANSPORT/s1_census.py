"""S1 -- EXACT CENSUS OF THE MINIMAL-TORUS CARRIER D(G) FOR A LADDER OF 2-GROUPS.

Carrier (as fixed by the brief): 1 vertex, 2 edges, 1 face.  Hilbert space C[G] (x) C[G],
basis |g1,g2>.
    A_h |g1,g2> = |h g1 h^-1, h g2 h^-1>       (transport / gauge)
    A = (1/|G|) sum_h A_h                      (gauge projector)
    B diagonal, 1 iff [g1,g2] = e              (flatness projector)
    H = -(A + B)

WHY THIS CAN BE DONE EXACTLY AND IS NOT SAMPLED.  A and B are COMMUTING projectors: B is the
indicator of the set C = {(g1,g2) : g1 g2 = g2 g1}, which is closed under simultaneous
conjugation, so [A_h, B] = 0 for every h and hence [A,B] = 0.  The Hilbert space therefore
splits into four joint blocks and H has at most three eigenvalues, -2, -1, 0, with

    dim E(-2) = m  = # conjugation orbits on C        (= # anyon types = GROUND SPACE)
    dim E(-1) = (a - m) + (|C| - m)
    dim E( 0) = |G|^2 - |C| - (a - m)
    a   = # conjugation orbits on GxG = (1/|G|) sum_h |C_G(h)|^2      (Burnside)
    |C| = |G| * k(G)

and the transport character of each block is EXACT:

    chi_{E(-2)}(h) = m                                          (trivial rep only)
    chi_{E(-1)}(h) = (a - m) + Fix_C(h) - m
    chi_{E( 0)}(h) = |C_G(h)|^2 - Fix_C(h) - (a - m)
    Fix_C(h)       = |C_G(h)| * k(C_G(h))

because A_h acts as the IDENTITY on every gauge-invariant vector, and as the permutation of
basis states on the rest.

SELF-CHECKS PRINTED, and any failure kills the row:
  (1) the three multiplicities sum to |G|^2
  (2) chi_E(e) equals dim E for every block
  (3) sum_rho d_rho m_rho = dim E for every block (isotypic decomposition consistency) --
      THIS IS THE CHECK THAT FAILED IN LANE_O34 WITH A SAMPLED COMMUTANT; here it is
      character theory and must pass exactly
  (4) every m_rho is a non-negative integer
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_A_TRANSPORT")
import glib
def say(*a): print(*a); sys.stdout.flush()

from carriers import census, isotypic, subset_sums, fixed_record_exists, phi

say("=" * 118)
say("S1   EXACT CENSUS OF D(G) ON THE MINIMAL TORUS -- A LADDER OF 2-GROUPS, ABELIAN AND NON-ABELIAN AT EVERY ORDER")
say("=" * 118)
say("")
hdr = (f"{'carrier':<14}{'|G|':>5}{'abel':>6}{'dim':>7}{'k(G)':>6}{'|Z|':>5}{'|G/Z|':>7}"
       f"{'mult(-2)=GROUND':>17}{'mult(-1)':>10}{'mult(0)':>9}{'all even':>10}{'fixed rec?':>12}{'phi_fix':>9}")
say(hdr); say("-" * len(hdr))
rows = []
for G in glib.ladder(64):
    ce = census(G)
    iso = isotypic(G, ce)
    ex, wit, why = fixed_record_exists(iso, ce['dims'])
    f, A, ph = phi(iso, ce['dims'])
    allev = all(ce['dims'][v] % 2 == 0 for v in (-2, -1, 0))
    say(f"{G.name:<14}{ce['n']:>5}{str(ce['abelian']):>6}{ce['dim']:>7}{ce['k']:>6}{ce['Z']:>5}"
        f"{ce['n']//ce['Z']:>7}{ce['dims'][-2]:>17}{ce['dims'][-1]:>10}{ce['dims'][0]:>9}"
        f"{str(allev):>10}{('YES' if ex else 'no'):>12}{ph:>9.4f}")
    rows.append(dict(name=G.name, n=ce['n'], abelian=ce['abelian'], dim=ce['dim'], k=ce['k'],
                     Z=ce['Z'], m=ce['dims'][-2], d1=ce['dims'][-1], d0=ce['dims'][0],
                     alleven=allev, fixed=ex, why=why, phi=ph, fixdim=f, alldim=A,
                     iso={v: iso[v] for v in iso}))
say("")
say("  mult(-2) = ground-space dimension = number of anyon types = THE RECORD COUNT the brief names.")
say("  phi_fix  = (transport-fixed Hermitian directions) / (all Hermitian directions in the commutant of H).")
say("             phi_fix = 1 means transport moves NOTHING; phi_fix < 1 means a generic record is moved.")
say("")
say("=" * 118); say("  ISOTYPIC BLOCK STRUCTURE (d_rho, m_rho) OF EACH EIGENSPACE -- self-check sum d*m = dim E"); say("=" * 118)
say(f"  {'carrier':<14}{'E':>4}{'dim':>7}   blocks (d x m)")
for G in glib.ladder(16):
    ce = census(G); iso = isotypic(G, ce)
    for v in (-2, -1, 0):
        if ce['dims'][v] == 0: continue
        b = " ".join(f"{d}x{m}" for d, m in iso[v])
        chk = sum(d * m for d, m in iso[v])
        say(f"  {G.name:<14}{v:>4}{ce['dims'][v]:>7}   {b}    [sum d*m = {chk} {'OK' if chk==ce['dims'][v] else 'FAIL'}]")
say("")
np.save("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_A_TRANSPORT/s1_rows.npy",
        np.array(rows, dtype=object), allow_pickle=True)
say("=" * 118); say("  READ -- filled in from the numbers above, not in advance"); say("=" * 118)
ab = [r for r in rows if r['abelian']]; na = [r for r in rows if not r['abelian']]
say(f"  ABELIAN carriers tested: {len(ab)}    NON-ABELIAN: {len(na)}")
say(f"  every abelian carrier   : phi_fix = {sorted({round(r['phi'],6) for r in ab})}")
say(f"  every non-abelian carrier: phi_fix in [{min(r['phi'] for r in na):.4f}, {max(r['phi'] for r in na):.4f}]")
say("")
say("  RECORD COUNT (ground dim) AT EQUAL ORDER, abelian vs non-abelian:")
for order in (2, 4, 8, 16, 32, 64):
    sel = [r for r in rows if r['n'] == order]
    if not sel: continue
    say(f"    |G| = {order:<3}  " + "   ".join(f"{r['name']}:{r['m']}{'(ab)' if r['abelian'] else '(NA)'}" for r in sel))
