# Builder self-audit -- ice-projected hybrid tensor response

**Lane:** `GRA-FK-F3-Q4-IHTR-V001`

**Verdict:**
`PASS__EXACT_ICE_REPRESENTATION_AND_HYBRID_REPRESENTATION_CANDIDATE__SYMMETRIC_FISHER_T2_NO_GO__EXACT_INHERITED_RING_RESPONSE__PHYSICAL_METRIC_MASSLESS_TENSOR_AND_GRAVITY_OPEN`

## 1. Representation result checked

The local `d_*=2` constraint leaves the six two-plus/two-minus sign states.
It imposes the exact operator identities

\[
 \sum_as_a=0,
 \qquad j_{12}=j_{34},\quad j_{13}=j_{24},\quad j_{14}=j_{23},
 \qquad \sum_{a<b}j_{ab}=-2I.
\]

The one-link span is therefore the three-dimensional tetrahedral `T2`
module.  The pair span is only `A1+E`: opposite-pair differences vanish, and
its `A1` observable is a fixed constant.  The normalized pair tangent is
exactly `E`, of dimension two.  Direct character replay gives

```text
six ice states: (6,2,2,0,0) = A1 + E + T2
one links:      (3,1,-1,0,-1) = T2
pairs:          (3,1,3,0,1) = A1 + E
centered pairs: (2,0,2,-1,0) = E
```

This corrects, rather than hides, the projection issue: FJ's six pair
operators are independent on the unconstrained four-qubit factor, but not
after compression to ice.

## 2. Tensor-rank result checked

The five nonconstant ice statistics—three one-link contrasts and two
centered matching/pair contrasts—have rank five, positive uniform Fisher
covariance, and together with the constant span all six diagonal functions.
Thus their finite exponential family is locally saturated on the six-state
simplex.  This is only a mathematical preparation family.

The frozen PMMDC one-body intertwiner carries one-link `T2` isomorphically to
the `T2` tensor summand. Its edge intertwiner restricted to
opposite-symmetric, zero-sum coefficients carries pair `E` isomorphically to
tensor `E`. Adding one external scalar multiple of `I_V` gives ranks
`3+2+1=6` and an exact **representation-isomorphism candidate** onto
`Sym^2(V)`. This does not identify those intertwiners with an actual physical
metric response. Three nonzero sector normalizations, their relative signs,
and their common physical calibration remain open.

The scalar is not manufactured from a uniform pair coupling.  Inside ice,
that coupling is `lambda sum j_ab=-2lambda I` and disappears into the
normalization. Same-parent physical ownership and calibration of a genuine
accumulation/scale scalar remain open.

The most direct symmetric ice Fisher query was checked separately. With
`s in V` as the localization statistic, its uniform covariance is exactly
`(4/3) I_V`. Differentiation in a one-link `T2` source gives zero at first
order by global sign symmetry. Differentiation in a six-edge pair-`E` source
is exactly `(8/3) M(y)` and has rank two. Thus the actual symmetric query
realizes the pair `E` metric tangent but has a first-order `T2` no-go. The
saturated state-family rank and the nonzero one-link operator dynamics do not
erase that physical-metric ceiling.

## 3. Ring dynamics checked

The inherited sixth-order coefficient was independently re-summed over all
720 virtual orders and reproduced `63/8`.  For every one of the 64 Walsh
subsets on one hexagon, direct finite matrix replay verifies

\[
 B_CW_A=(-1)^{|A\cap C|}W_AB_C.
\]

Thus one-link observables on the ring and pair observables containing exactly
one ring edge are nonconserved without adding an interaction.  Local allowed
ring moves span all one-link `T2` and pair `E` directions.  Two odd
observables on a common ring have the exact nested commutator

\[
 [[H,A],D]=-4J_6ADB_C
\]

for that ring contribution, and the corresponding full sum is nonzero as an
operator.

Compression to one alternating configuration and its ring-flipped partner
gives `H_Q=-J_6 sigma_x`.  Odd observables become signed `sigma_z`, so the
cross response is exactly `o_A o_D R_(2J_6)(z)` and is nonzero on the positive
imaginary-energy axis.  The audit keeps the crucial qualifier: this is the
compressed linked-ring resolvent.  Other plaquettes can leave the block, so
it is not silently promoted to the full many-ring Green function.

## 4. Scientific ceilings checked

- The finite poles are at `+/-2J_6`; they are not massless.
- The imported diamond-ice `U(1)` phase is spin one, not a helicity-two mode.
- The actual symmetric Fisher query has no first-order `T2` metric tangent.
- The independent scalar has no ring dynamics or common-port bind yet, and
  the three representation-sector normalizations are not physically fixed.
- Active link/pair observables are not called records solely because their
  support word `K_e` is authenticated.
- No arbitrary PMMDC six-pair preparation is claimed inside ice.
- Support selection, compatible boundary completion, refinement, all-orders
  control, tensor constraints, universal stress coupling, RGRL-B, Einstein
  response, gravity, and `G` remain open.

## 5. Reproduction

Run:

```text
python3 LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/verify_ice_hybrid_tensor_response.py
```

Expected result: `SUMMARY 65/65 PASS`.
