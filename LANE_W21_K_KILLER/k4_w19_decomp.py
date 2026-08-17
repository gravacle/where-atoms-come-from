"""
k4_w19_decomp.py -- LANE W21-K, BLOCK 4.  WHAT W-19'S ONE FULL BIT ACTUALLY IS.

The brief says W-19 measured "the disagreement between the extended and the gauge-invariant
channel" as ONE FULL BIT in the magnetic phase and EXACTLY ZERO in the electric phase, and
offers it as the finite instance of the Casini-Huerta-Rosabal centre ambiguity, which the
type III -> type II literature is then supposed to resolve.

THIS BLOCK SHOWS THAT IS A MISIDENTIFICATION, AND SAYS WHAT THE BIT IS INSTEAD.

  4A  THEOREM + EXHIBIT: for Z_N pure gauge theory on any graph, any physical state, any set of
      links R, the EXTENDED-HILBERT-SPACE entropy of R equals the entropy of the FULL
      gauge-invariant algebra of R, exactly.  The extended construction adds NOTHING.  So no
      "extended vs gauge-invariant" gap can be a gauge-variance effect.
  4B  W-19's numbers reproduced from scratch on theta_6 / theta_8, and the bit decomposed.
  4C  THE ISOLATING DIFF: one variable, "does the union of the two regions carry a cycle".
  4D  The centre ambiguity of CHR, measured on the same page, so the two are visibly different
      quantities of different sizes.
"""
import numpy as np, itertools
from klib import (PauliAlgebra, Z2Gauge, pack, unpack, rref_basis, span_elements, in_span,
                  gauge_invariant_subspace, sympl_perp, intersect, pauli_matrix)

np.set_printoptions(precision=9, suppress=True)
def line(c="="): print(c * 100)
rng = np.random.default_rng(419)


def vn(rho):
    w = np.clip(np.linalg.eigvalsh((rho + rho.conj().T) / 2).real, 0, None)
    return float(-sum(x * np.log2(x) for x in w if x > 1e-14))


def ptrace_links(psi, L, keep):
    pos = sorted(L - 1 - l for l in keep)
    rest = [i for i in range(L) if i not in pos]
    t = psi.reshape([2] * L)
    t = np.transpose(t, pos + rest).reshape(1 << len(pos), 1 << len(rest))
    return t @ t.conj().T


def ginv_algebra(G, links, name=""):
    """The FULL gauge-invariant algebra of a set of links: all X-strings on those links, plus
    Z-strings whose support is a cycle inside them."""
    n = G.L
    m = sum(1 << l for l in links)
    gens = [pack((1 << l, 0), n) for l in links]
    for c in span_elements(G.cycles):
        if c and (c & ~m) == 0:
            gens.append(pack((0, c), n))
    return PauliAlgebra(gens, n, name)


def elec_algebra(G, links, name=""):
    n = G.L
    return PauliAlgebra([pack((1 << l, 0), n) for l in links], n, name)


# ============================================================================ 4A
line()
print("BLOCK 4A -- THEOREM: THE EXTENDED HILBERT SPACE ADDS NOTHING TO Z_N PURE GAUGE THEORY.")
line()
print("""STATEMENT.  Let R be any set of links, |psi> any state in the Gauss-physical sector.  Let
S_ext(R) = von Neumann entropy of Tr_{R^c}|psi><psi| in the UNCONSTRAINED tensor product, and let
A_gi(R) be the gauge-invariant algebra of R.  Then S_ext(R) = S(rho|_{A_gi(R)}) EXACTLY.

PROOF.  (i) A Pauli supported on R is gauge-invariant iff its Z-part is a cycle of R; any other
Pauli anticommutes with some Gauss operator G_v, so <psi|O|psi> = <psi|G_v O G_v|psi> = -<psi|O|psi>
= 0.  Hence rho_R, expanded in the Pauli basis of R, is supported ENTIRELY on A_gi(R), so the
Hilbert-Schmidt projection of the state onto A_gi(R) is rho_R itself: omega is the same matrix in
both computations.
(ii) The block data also agree.  A_gi(R) = {(x,z) : x in F_2^{|R|}, z in Cyc(R)}, of F_2-dimension
|R| + c_R where c_R = dim Cyc(R).  Its symplectic radical is {(x,0) : x . z = 0 for all z in
Cyc(R)} = the CUT SPACE of R, of dimension |R| - c_R.  So r = |R| - c_R and k = c_R, giving
r + k = |R|, which is exactly the k of the full |R|-link matrix algebra (r=0).  The multiplicity
m = 2^{L-r-k} is therefore the same for both.  Same omega, same m, same entropy.  QED.

CONSEQUENCE, AND IT IS THE ONE THAT MATTERS FOR THIS ROUND: no gap between an "extended" channel
and a "gauge-invariant" channel can be attributed to gauge-variance.  If a gap is measured, the
gauge-invariant side was a PROPER SUBALGEBRA of the region's own algebra, and the gap is the
content of the generators that were dropped.

The exhibit below is a CODE CHECK on a proved statement.  It could not have failed and is void
as evidence; it is here to certify the code, and because this program has three times published
a forced identity as a result.""")

carriers = [
    ("theta_6", 2, [(0, 1)] * 6),
    ("ladder3", 6, [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]),
    ("k4", 4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]),
]
worst = 0.0
ncheck = 0
print(f"\n  {'carrier':>9} {'g^2':>6} {'region R':>22} {'|R|':>4} {'c_R':>4} "
      f"{'S_ext':>14} {'S(A_gi(R))':>14} {'diff':>10}")
for nm, nv, E in carriers:
    G = Z2Gauge(nv, E, nm)
    for g2 in (0.10, 1.00, 3.00):
        psi, pdim, e0 = G.ground_state(g2)
        for trial in range(4):
            k = rng.integers(1, min(4, G.L) + 1)
            R = sorted(rng.choice(G.L, size=int(k), replace=False).tolist())
            m = sum(1 << l for l in R)
            cR = len([1 for c in span_elements(G.cycles) if c and (c & ~m) == 0])
            cR = int(np.log2(cR + 1)) if cR else 0
            A = ginv_algebra(G, R)
            se = vn(ptrace_links(psi, G.L, R))
            sa = A.entropy(psi)
            worst = max(worst, abs(se - sa)); ncheck += 1
            print(f"  {nm:>9} {g2:6.2f} {str(R):>22} {len(R):>4} {cR:>4} "
                  f"{se:14.9f} {sa:14.9f} {abs(se-sa):10.2e}")
print(f"\n  >>> {ncheck} checks, worst |S_ext - S(A_gi)| = {worst:.3e}.  [FORCED -- a code check]")


# ============================================================================ 4B
line()
print("BLOCK 4B -- W-19'S ONE FULL BIT, REPRODUCED FROM SCRATCH AND DECOMPOSED.")
line()
for nm, nlinks in (("theta_6", 6), ("theta_8", 8)):
    G = Z2Gauge(2, [(0, 1)] * nlinks, nm)
    n = G.L
    # the magnetic GHZ = the exact g^2 -> 0 ground state
    ghz = np.zeros(1 << n, dtype=complex)
    ghz[0] = ghz[(1 << n) - 1] = 1 / np.sqrt(2)
    psi_small, pdim, e0 = G.ground_state(0.02)
    ov = abs(np.vdot(ghz, psi_small))
    print(f"\n  carrier {nm}: V=2 L={n} cycle_dim={G.cycle_dim} physical dim={pdim}")
    print(f"    magnetic GHZ vs the g^2=0.02 ground state: |<GHZ|psi_0>| = {ov:.12f}   "
          f"Gauss residual {G.gauss_residual(ghz):.1e}")

    S_, F_ = [0], [1]
    A_ext_S = PauliAlgebra([pack((1, 0), n), pack((0, 1), n)], n)
    A_gi_S = ginv_algebra(G, S_); A_gi_F = ginv_algebra(G, F_)
    A_gi_SF = ginv_algebra(G, S_ + F_)
    A_join = PauliAlgebra(A_gi_S.basis + A_gi_F.basis, n)
    A_el_SF = elec_algebra(G, S_ + F_)

    sS = vn(ptrace_links(ghz, n, S_)); sF = vn(ptrace_links(ghz, n, F_))
    sSF = vn(ptrace_links(ghz, n, S_ + F_))
    I_ext = sS + sF - sSF
    I_gi = A_gi_S.entropy(ghz) + A_gi_F.entropy(ghz) - A_gi_SF.entropy(ghz)
    I_elec = A_gi_S.entropy(ghz) + A_gi_F.entropy(ghz) - A_join.entropy(ghz)
    print(f"    S_ext(S)={sS:.9f}  S_ext(F)={sF:.9f}  S_ext(S u F)={sSF:.9f}")
    print(f"    I_EXTENDED  (W-19's 'EXT'  channel)                          = {I_ext:.9f} bits")
    print(f"    I with the FULL gauge-invariant algebra of S u F             = {I_gi:.9f} bits")
    print(f"    I with alg{{X_l}} on each side (W-19's gauge-invariant channel) = {I_elec:.9f} bits")
    print(f"    W-19 reported 1.000000 and 0.000000 for these two.  Reproduced independently.")
    print(f"    ALGEBRA DIMENSIONS: A_gi(S) v A_gi(F) has dim_F2 = {A_join.dimV}; "
          f"A_gi(S u F) has dim_F2 = {A_gi_SF.dimV}")
    missing = [b for b in A_gi_SF.basis if not in_span(b, A_join.basis)]
    for b in missing:
        x, z = unpack(b, n)
        gv = all(bin(g[0] & z).count("1") % 2 == 0 for g in G.gauss)
        print(f"    THE GENERATOR THE JOIN IS MISSING: X-part {bin(x)} Z-part {bin(z)}  "
              f"-> the Wilson loop of the 2-link cycle {{0,1}}.  GAUGE-INVARIANT: {gv}")
    print(f"    S(alg{{X_0,X_1}}) = {A_join.entropy(ghz):.9f}   "
          f"S(A_gi({{0,1}})) = {A_gi_SF.entropy(ghz):.9f}   "
          f"DIFFERENCE = {A_join.entropy(ghz)-A_gi_SF.entropy(ghz):.9f} bits")
    print(f"    >>> THE ENTIRE BIT IS THAT ONE WILSON LOOP.  It is gauge-INVARIANT, so the bit is")
    print(f"        NOT gauge-variance; and A_gi(S u F) is UNIQUELY determined by the link set, so")
    print(f"        the bit is NOT a centre choice either.  It is the failure of the region-to-")
    print(f"        algebra map to be a join: A(S u F) strictly contains A(S) v A(F).")

    # larger fragments, same structure
    print(f"    same structure at larger |F|:")
    for fs in (2, 3):
        F2 = list(range(1, 1 + fs))
        gS, gF = ginv_algebra(G, S_), ginv_algebra(G, F2)
        gSF = ginv_algebra(G, S_ + F2)
        jn = PauliAlgebra(gS.basis + gF.basis, n)
        ie = (vn(ptrace_links(ghz, n, S_)) + vn(ptrace_links(ghz, n, F2))
              - vn(ptrace_links(ghz, n, S_ + F2)))
        print(f"      |F|={fs}: I_EXT={ie:.9f}  I_full-gi={gS.entropy(ghz)+gF.entropy(ghz)-gSF.entropy(ghz):.9f}"
              f"  I_join(no loops)={gS.entropy(ghz)+gF.entropy(ghz)-jn.entropy(ghz):.9f}"
              f"  loops missing from the join: {gSF.dimV-jn.dimV}")


# ============================================================================ 4C
line()
print("BLOCK 4C -- THE ISOLATING DIFF.  ONE VARIABLE: DOES S u F CARRY A CYCLE.")
line()
print("Two link pairs on the SAME carrier, same state, same coupling, same algebras rule.")
G = Z2Gauge(6, [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)], "ladder3")
n = G.L
G2 = Z2Gauge(2, [(0, 1)] * 6, "theta_6")
print(f"\n  {'carrier':>9} {'S':>5} {'F':>5} {'cycle in SuF?':>14} {'g^2':>6} "
      f"{'I_EXT':>12} {'I_full-gi':>12} {'I_join':>12} {'EXT - join':>12}")
for GG, pairs in ((G, [([0], [1]), ([0], [2])]), (G2, [([0], [1]), ([0], [1, 2])])):
    nn = GG.L
    for g2 in (0.10, 1.00, 3.00):
        psi, pdim, e0 = GG.ground_state(g2)
        for Ss, Fs in pairs:
            m = sum(1 << l for l in Ss + Fs)
            hascyc = any(c and (c & ~m) == 0 for c in span_elements(GG.cycles))
            gS, gF = ginv_algebra(GG, Ss), ginv_algebra(GG, Fs)
            gSF = ginv_algebra(GG, Ss + Fs)
            jn = PauliAlgebra(gS.basis + gF.basis, nn)
            ie = (vn(ptrace_links(psi, nn, Ss)) + vn(ptrace_links(psi, nn, Fs))
                  - vn(ptrace_links(psi, nn, Ss + Fs)))
            ig = gS.entropy(psi) + gF.entropy(psi) - gSF.entropy(psi)
            ij = gS.entropy(psi) + gF.entropy(psi) - jn.entropy(psi)
            print(f"  {GG.name:>9} {str(Ss):>5} {str(Fs):>5} {str(hascyc):>14} {g2:6.2f} "
                  f"{ie:12.9f} {ig:12.9f} {ij:12.9f} {ie-ij:12.9f}")
print("""
  THE COLUMN 'EXT - join' IS ZERO EXACTLY WHEN S u F CARRIES NO CYCLE AND NON-ZERO WHEN IT DOES.
  That is the whole of W-19's finding restated with its mechanism visible.  The coupling
  dependence W-19 read as a PHASE INDEX is the coupling dependence of <W_loop>: at small g^2 the
  ground state is a Wilson-loop eigenstate and the loop carries a full bit; at large g^2 it is an
  electric eigenstate and the loop carries nothing.  Nothing in it is about von Neumann type, and
  nothing in it is about a choice of centre.""")


# ============================================================================ 4D
line()
print("BLOCK 4D -- THE TWO QUANTITIES SIDE BY SIDE, SO THEY ARE NOT CONFUSED AGAIN.")
line()
print("""  QUANTITY 1 -- W-19's bit.  A(S u F) minus A(S) v A(F).  Region-to-algebra assignment is
    FIXED and unique; nothing is chosen.  The gap is the region's own Wilson loops.  Size here:
    exactly 1.000000000 bits on theta_6 in the magnetic phase, 0 in the electric phase.
    THIS IS NOT AN AMBIGUITY.  It is a determinate structural fact about gauge theory:
    the map region -> algebra is not a lattice homomorphism.

  QUANTITY 2 -- the Casini-Huerta-Rosabal centre ambiguity.  Region FIXED; the algebra assigned
    to it is CHOSEN from every subalgebra between "what is strictly inside" and "the commutant of
    what is strictly inside the complement".  Measured exhaustively in k2b on grid3x3, region
    (0,1,2,4): 67 admissible algebras, 41 distinct entropies, spread 3.998925600 bits at
    g^2 = 0.10 and 2.851355785 bits at g^2 = 0.50.
    THIS IS THE AMBIGUITY.  It is larger than W-19's bit, it is not phase-indexed in the same
    way, and W-19 did not measure it.

  Both are finite, exact, and fully computable.  NEITHER involves a trace that fails to exist.""")

line()
print("END BLOCK 4")
