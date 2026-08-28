# Independent hostile audit -- ice-projected q4/F3 hybrid tensor response

**Lane:** `GRA-FK-F3-Q4-IHTR-V001`

**Date:** 2026-08-27

## Verdict

`ACCEPT_AFTER_GLOBAL_SUPPORT_AND_FLIPPABILITY_DOMAIN_REPAIR__LOCAL_ICE_REPRESENTATION_FISHER_AND_SIXTH_ORDER_RING_RESULTS_SURVIVE__PHYSICAL_T2_METRIC_SOLDER_MASSLESS_TENSOR_AND_GRAVITY_OPEN`

The hostile audit independently reproduced the local representation ranks,
the symmetric-query Fisher derivatives, the sixth-order ring coefficient,
the Walsh commutators, and the compressed spectral response.  No error was
found in those calculations.  One real promotion-boundary defect was found
and repaired: the original result chain blurred the authenticated interior
q4 links supplied by `FH` with the separately supplied regular global domain
required by `FE/CW`.  `FH` itself proves that the raw finite q4 slab has an
empty global degree-two sector.  The repaired theorem now requires a
compatible boundary/periodic completion or controlled infinite-support
definition for every global use of the ice Hamiltonian, without calling any
added boundary links q4-authenticated.  It also states the two premises
needed by every nonzero response claim: `h\ne0` and an admitted hexagon with a
compatible ice state in which that hexagon is flippable.

The original `SELF_AUDIT.md` is preserved unchanged as builder custody.  This
independent audit governs the repaired theorem.

## 1. Dependency and physical-domain audit

The audit checked the exact bytes of all five theorem dependencies.  It also
pins the final `FH` re-audit and the final `FJ` independent audit, because
both contain load-bearing type/domain qualifications:

- `FH` distinguishes authenticated support memory `K_e` from active link
  occupation `n_e`, requires the raw ungated flip to be absent during
  quarantine, and proves `Omega_2(E_N)=emptyset` on the raw finite slab.
- `FJ` realizes the six Walsh pairs on the unprojected four-link factor, but
  does not identify them with a prepared PMMDC family or qualified records.
- `FE/CW` supply the degree-two diamond-ice Hamiltonian only on a compatible
  coordination-four global domain.  Its finite raw-slab use requires added
  boundary physics; its alternative is a controlled infinite-support
  definition.

Thus `FK` lawfully combines authenticated **interior** q4/F3 links with an
ice calculation only under the now-explicit compatible-global-domain
antecedent.  The local algebra is not a theorem that BQ4 lineage selects or
authenticates a periodic boundary completion.

The inherited coefficient is

\[
 J_6={63h^6\over8U_d^5}.
\]

For `U_d>0` it is nonnegative and is strictly positive exactly when
`h\ne0`.  The amended theorem no longer silently includes `h=0` in claims of
a nonzero commutator or response.

## 2. Exact local representation replay

The six local ice states are the sign vectors with two plus and two minus
entries.  Direct integer enumeration gives

\[
 \sum_as_a=0,\qquad
 j_{12}=j_{34},\quad j_{13}=j_{24},\quad j_{14}=j_{23},\qquad
 \sum_{a<b}j_{ab}=-2.
\]

The one-link evaluation matrix has rank three.  The pair evaluation matrix
has rank three, while its centered matrix has rank two.  The exact characters
on the class order `1,(12),(12)(34),(123),(1234)` are

```text
six-state functions  (6,2,2,0,0) = A1 + E + T2
one-link span         (3,1,-1,0,-1) = T2
pair span             (3,1,3,0,1) = A1 + E
centered pair span    (2,0,2,-1,0) = E
```

The constant, three one-link directions, and two centered matching
directions have joint rank six as functions on the six states.  Hence the
five nonconstant statistics give the full normalized tangent space; this is
a finite mathematical exponential family, not a physical preparation
theorem.

The PMMDC maps were replayed independently in a three-dimensional
orthonormal basis of `V=1^perp`.  The three tensors

\[
 P\operatorname{diag}(x)P\big|_V
\]

have rank three and are trace-free.  The two complementary-edge-symmetric,
zero-sum pair tensors

\[
 \sum_{a<b}y_{ab}(v_av_b^{\mathsf T}+v_bv_a^{\mathsf T})
\]

have rank two and are trace-free.  Their joint span is rank five, and adding
`I_V` gives rank six.  This verifies only the multiplicity-one abstract
isomorphism

\[
 A_1\oplus E\oplus T_2\cong\operatorname{Sym}^2(V).
\]

It does not determine a physical scalar owner, the relative sector
normalizations, or a metric-response map.

## 3. Exact Fisher-query replay

Under the uniform ice measure, exact rational moments give

\[
 \mathbb E[s_as_b]
 =\begin{cases}1,&a=b,\\-1/3,&a\ne b,\end{cases}
 \qquad
 \operatorname{Cov}(s)={4\over3}P.
\]

All cubic moments vanish exactly under global sign reversal.  Therefore the
first derivative of `Cov(s)` in every one-link score is zero; this confirms
the claimed first-order `T2` no-go without relying on floating-point rank
tests.

For the two independent `E` sources, chosen as duplicated matching weights
`(1,-1,0)` and `(1,1,-2)`, exact enumeration of the six states gives

\[
 D_yF_{\rm ice}(0,0)={8\over3}M(y).
\]

Linearity then proves the formula for every `y` in `E_pair`.  The resulting
metric tangent has rank two.  A uniform pair source is exactly `-2` on every
ice state, so it changes only the partition-function normalization and
cannot supply an `A1` probability tangent.  These facts confirm the packet's
most important ceiling: representation availability is not physical metric
availability, and the direct symmetric query is still missing a nonzero
first-order `T2` response.

## 4. Ring coefficient, commutators, and propagation replay

The audit independently enumerated all `6!=720` virtual flip orders.  After
sorting the five intermediate gap coefficients, the multiplicities are

```text
(2,2,2,2,2)  96
(2,2,2,2,4) 144
(2,2,2,4,4) 216
(2,2,4,4,4) 192
(2,2,4,4,6)  72
```

and their exact rational denominator sum is `63/8`.  The local diamond
three-label hexagons expose all six unordered pairs of incident append
labels, so their allowed local moves span the rank-three one-link and
rank-two centered-pair sectors as claimed.

For all 64 Walsh subsets of one alternating hexagon, direct matrix replay
gives

\[
 B_CW_A=(-1)^{|A\cap C|}W_AB_C.
\]

It follows algebraically that

\[
 [H,W_A]=2J_6\sum_{C:\,|A\cap C|\ \mathrm{odd}}W_AB_C
\]

and, for two commuting diagonal observables odd on the same ring,

\[
 [[H,A],D]=-4J_6ADB_C.
\]

The original text inferred a nonzero full-sum matrix element without first
stating that the selected ring is flippable in a compatible global ice
state.  The repaired theorem now conditions on such a state.  Then the
matrix element from `|n>` to `|n triangle C>` is nonzero, and no distinct
hexagon can cancel it because a different edge set gives a different
symmetric difference.

## 5. Compressed response and exact ceiling

For the two states related by one admitted flippable ring, compression gives

\[
 H_Q=E_CI-J_6\sigma_x,
 \qquad O_Q=o_O\sigma_z
\]

for every odd one-link or crossing-pair observable.  Direct diagonalization
reproduces

\[
 \chi^{R,Q_C}_{OD}(z)
 =o_Oo_D\left({1\over z-2J_6}-{1\over z+2J_6}\right).
\]

At `z=i kappa`, `kappa>0`, this equals
`-4o_Oo_D J_6/(kappa^2+4J_6^2)` and is nonzero for `h\ne0`.  This is exactly
the response of the compressed two-state Hamiltonian.  Other plaquettes can
leave that subspace, so it is not the resolvent of the full many-ring model.
Its poles are finite-gap poles at `+/-2J_6`, not a massless helicity-two pole.

Accordingly the strongest audited result is an exact ice-projected local
module, an exact abstract hybrid tensor representation candidate, an exact
pair-`E` but zero first-order one-link-`T2` Fisher response at the symmetric
point, and an exact finite ring-mediated operator response.  Physical
`T2` metric solder, scalar same-parent ownership, sector calibration,
preparation/record qualification, thermodynamic tensor spectrum, helicity
two, refinement/gluing, universal stress coupling, RGRL-B, gravity, and `G`
remain open.

## 6. Reproduction and frozen hashes

Run:

```text
python3 LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/verify_ice_hybrid_tensor_response.py
```

Expected result: `SUMMARY 74/74 PASS`.

Canonical content hashes before adding this audit to the manifest:

```text
THEOREM.md                          cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98
verify_ice_hybrid_tensor_response.py 7897f92ea2b21b219ed4ae99cf2cbeef0aece599a3a82d2514b7d6a267a551d8
SELF_AUDIT.md                       2e1b5c41343574334d98e448b3446456c1d6b57f46c2d8f21c1294144e69af29
```

The manifest pins the final theorem, verifier, preserved self-audit,
independent audit, and regenerated verification transcript.
