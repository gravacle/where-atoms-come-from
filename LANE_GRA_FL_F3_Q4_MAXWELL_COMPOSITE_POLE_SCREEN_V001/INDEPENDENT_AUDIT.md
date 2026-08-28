# Independent hostile audit -- Gaussian-Maxwell composite pole screen

**Lane:** `GRA-FL-F3-Q4-MCPS-V001`

**Audit date:** 2026-08-27

## Verdict

`ACCEPT_AFTER_THRESHOLD_ONSET_TD_CUSTODY_AND_RESIDUE_REPAIR__EXACT_ICE_PARITY_MAP_AND_CONDITIONAL_GAUSSIAN_NO_GO_SURVIVE__NO_HELICITY2_OR_GRAVITY_CLOSURE`

The exact six-state ice projection, incidence-oriented flux map, complement
parity, and finite-group decompositions all survive independent replay.  The
conditional Maxwell calculation also survives: the odd one-link channel has
the transverse spin-one pole, while the even centered-pair channel has no
one-photon matrix element and begins with a two-photon continuum.  No isolated
helicity-two pole is present at the Gaussian fixed point.

The hostile audit repaired four scoped precision defects without changing
that conclusion.  It removed an unjustified universal threshold-suppression
claim, made the full-tetrahedral `T_d \simeq S_4` polar convention explicit,
separated microscopic parity from the conditional symmetric-state selection
rule, and stated the positive-frequency electric-field residue and the
electric/magnetic dual naming convention exactly.

## 1. Frozen dependencies and primary-source custody

The verifier independently recomputes and matches the four load-bearing
dependency hashes printed in `THEOREM.md`:

- FJ theorem: `05f4a619a6f80aa40c48570ab4035ab874426502a31468a08f435e66610bd769`;
- FJ hostile audit: `44690ad431c85af7a4947a431a4c57ad4cd8b19a346e3d26720b341f77256f90`;
- CROSS-CW theorem: `5e68e4a8c62ad89cff309781a5cb54e071092e82594e60c04e7992414bc18dbe`;
- CROSS-CW primary-source ledger: `4ee84b4f9b78003cdc5ce80a86cba6cbab618feb1fcd78d25903cb5e97c42a62`.

I also checked the cited primary paper itself, Shannon et al., arXiv
`1105.4196v3`, against the APS DOI record.  Its Eq. (1) is the source-free
Maxwell action, Eq. (2) is the fully-packed-loop quantum-ice Hamiltonian,
Fig. 5(b) is the direct `mu=0` ED/GFMC flux-scaling test on 80-, 320-, 640-,
and 1280-site clusters, and Fig. 6 is the `mu=0` small-wave-vector structure
factor.  The paper expressly places `mu=0` inside the reported U(1) liquid and
describes the excitations by a `3+1`-dimensional Maxwell action.

This is numerical phase evidence for the published pure-ice model.  It is not
an all-orders thermodynamic theorem for the complete F3 parent.  CROSS-CW
internally earns only the sixth-order fixed-graph Hamiltonian and imports the
published infinite-system phase identification.  FL correctly keeps that
boundary in the named `MAXWELL-IR` premise.

The source calls the divergence-free microscopic arrow field `B`, whereas FL
uses the canonical lattice-Hamiltonian name electric flux.  In the sourceless
Gaussian phase these are electric-magnetic dual conventions.  Both have the
same transverse rank-two polarization projector, linear cone, complement
parity, and pole-versus-continuum classification.  This does not identify the
emergent U(1) with visible electromagnetism.

## 2. Exact ice, orientation, and pair-sector replay

Direct enumeration gives exactly six words in

\[
 \Omega_2=\{\varepsilon\in\{-1,+1\}^4:\sum_a\varepsilon_a=0\}.
\]

Their one-link value matrix has rank three.  The six pair columns have total
rank three, but after subtracting their state average they have rank two.  On
every word,

\[
 j_{12}=j_{34},\qquad j_{13}=j_{24},\qquad j_{14}=j_{23},
 \qquad \sum_{a<b}j_{ab}=-2.
\]

The exact edge projectors have ranks

\[
 \operatorname{rank}(P_{A_1},P_E,P_{T_2})=(1,2,3).
\]

The enumerated pair matrix obeys `pairs P_T2=0`; its centered part obeys
`centered pairs P_A1=0` and `centered pairs P_E=centered pairs`.  Thus the
ice-projected pair response is constant `A1` plus nonconstant `E`, not the
unprojected `A1+E+T2` module.

The endpoint repair is also exact.  With `eta_v=+1` on `V_+`, `eta_v=-1` on
`V_-`, and `epsilon_(v,e)=eta_v Z_e`, both `epsilon` and the outward bond
vector reverse across a shared edge.  Their product is endpoint independent,
and

\[
 \mathcal E_v={3\over4}\sum_a e_{v,a}\varepsilon_{v,a},
 \qquad \varepsilon_{v,a}=e_{v,a}\cdot\mathcal E_v
\]

holds on all six ice words.  The resulting local sum-zero condition is the
incidence-correct Gauss law.  No raw same-sign `Z_e` is silently treated as an
outward flux.

Complement maps each ice word to its negative.  It makes every one-link
observable odd and every pair observable even, and it commutes with the
symmetric pure-kinetic ring Hamiltonian.  These operator statements are
internal.  The zero one-photon pair matrix element additionally requires the
complement-even vacuum and odd photon supplied by the complement-symmetric
`MAXWELL-IR` state; the custody table now says so.

## 3. Full tetrahedral typing is not continuum spin

All 24 permutations preserve the exact `A1/E/T2` projectors.  Using the
geometric tetrahedral frame, every permutation has a unique polar orthogonal
matrix `R` satisfying `R e_a=e_(pi(a))`.  Even permutations have determinant
`+1`; odd permutations are the improper classes of the full point group
`T_d`.  The proper-rotation subgroup alone is `A_4`, so writing only
`ell downarrow S4` without the polar/full-`T_d` convention was ambiguous.

In the cycle-class order `1^4, 2 1^2, 2^2, 3 1, 4`, the sum-zero vertex module
has character

\[
 \chi_{T_2}=(3,1,-1,0,-1).
\]

For polar symmetric-traceless rank two,

\[
 \chi_{\mathrm{ST}^2}(R)
 ={(\operatorname{tr}R)^2+\operatorname{tr}(R^2)\over2}-1
 =(5,1,1,-1,-1)=\chi_E+\chi_{T_2}.
\]

The verifier checks these identities for every group element, not only by
dimension count.  They establish the finite `T_d` decomposition.  They do not
turn the one-link `T2` pole into continuum helicity two.

## 4. Maxwell one-link response

For each transverse momentum mode, the imported quadratic action is a
harmonic oscillator of mass `chi` and frequency `omega_k=c|k|`.  Canonical
normalization gives

\[
 |\langle0|E_i(\mathbf k)|\gamma(\mathbf k,\lambda)\rangle|^2
 ={c|\mathbf k|\over2\chi}
 \epsilon_i^{(\lambda)}\epsilon_i^{(\lambda)*}.
\]

Consequently the propagating retarded tensor is proportional to

\[
 P^T_{ij}(\mathbf k)
 {c^2|\mathbf k|^2/\chi\over
  (\omega+i0)^2-c^2|\mathbf k|^2}.
\]

Its positive-frequency complex residue is
`c|k| P_T/(2 chi)`: it vanishes linearly, while the normalized polarization
matrix has rank two.  Replacing the numerator by `omega^2` changes the
expression by the analytic equal-time term `P_T/chi`, not by a propagating
pole.  Pullback through the complete four-link tetrahedral frame preserves
rank two for every nonzero momentum.  This is a transverse vector/helicity
`+/-1` pole.

## 5. Pair parity, two-particle cut, and threshold repair

In a complement-even vacuum, an even pair operator has zero matrix element to
the odd one-photon sector.  The finite Fock replay verifies that an odd field
reaches one particle, while its normal-ordered square reaches two particles
but not one.  Operator renormalization cannot violate this exact selection
rule without breaking complement or expanding about a nonzero odd background.

At the Gaussian point the remaining gapless pair response is a Wick
convolution.  For two photons with energies `c|p|` and `c|k-p|`,

\[
 c|\mathbf p|+c|\mathbf k-\mathbf p|\ge c|\mathbf k|.
\]

The lower edge is attained by collinear co-propagating photons.  Taking
`p=-t khat` with `t>=0` gives energy `c|k|+2ct`, so continuously varying an
internal momentum fills all energies above the edge.  With nonzero allowed
vertex overlap this is a two-particle branch cut, not a delta-function pole.

The provisional theorem additionally claimed that every field-strength
bilinear is power-suppressed by invariant distance from threshold.  That is
too strong component by component.  For example, collinear photons moving
along `z` with `x` polarization have a nonzero threshold matrix element for
`:E_xE_x:`.  Phase-space and tensor projections determine the actual onset;
some channels vanish by powers and others need not share that exponent.  The
repaired theorem claims no universal threshold exponent.  This correction
does not weaken the branch-cut/no-delta conclusion.

Local composite counterterms and equal-time terms are analytic at the
propagating singularity.  Gapped microscopic states may add higher cuts or
poles, but cannot manufacture the massless Gaussian one-particle tensor pole
being screened here.

## 6. TT comparison and exact no-go scope

The spatial TT projector in the theorem is symmetric, idempotent,
transverse, traceless, and rank two.  Under a rotation by `theta` about
nonzero momentum, the Maxwell polarization plane has character
`2 cos(theta)`, while the TT plane has `2 cos(2 theta)`.  Equal polarization
counts therefore do not identify the representations.

At the Gaussian fixed point the odd link block contains the spin-one pole and
the even pair block contains the two-photon cut.  Their cross block vanishes
by complement parity.  A local `S4`-equivariant rearrangement cannot change
that particle content, parity, denominator, or little-group action.  Hence no
linear combination of the admitted FJ link and pair probes obtains an
isolated helicity-two pole by the direct Gaussian composite route.

The theorem does not exclude a same-parent non-Gaussian bound state generated
by inherited higher-order F3 interactions, nor a distinct rank-two constrained
phase.  Either route still owes a thermodynamic pole with nonzero residue,
the TT projector and tensor Ward/constraint identities, a common cone,
universal stress coupling, and stability.  No fitted `j-j` attraction may be
inserted, and this packet proves neither visible electromagnetism nor gravity.

## 7. Reproduction and final result

Run:

```text
python3 LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/verify_maxwell_composite_pole_screen.py
```

The strengthened replay reports `108 passed, 0 failed`.  It includes the six
ice words, all exact projectors, endpoint gluing, all 24 full-tetrahedral
polar actions, the `E+T2` rank-two character, vector and TT projectors,
Maxwell residue/contact identity, Fock parity witness, two-photon threshold,
and every dependency and claim ceiling.

Final disposition:

```text
ODD ONE-LINK: CONDITIONAL MASSLESS TRANSVERSE SPIN-1 POLE
EVEN CENTERED PAIR: NO ONE-PHOTON POLE; TWO-PHOTON CUT + CONTACTS
HELICITY 2: NO ISOLATED POLE AT THE GAUSSIAN FIXED POINT
SCOPE: DIRECT COMPOSITE ROUTE ONLY; ALL-ORDERS F3, VISIBLE EM, AND GRAVITY OPEN
```

