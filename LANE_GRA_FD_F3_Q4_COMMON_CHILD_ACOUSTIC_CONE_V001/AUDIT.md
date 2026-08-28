# Independent hostile audit -- F3/q4 common-child acoustic cone

**Lane:** `GRA-FD-F3-Q4-CCMAC-V001`  
**Audit date:** 2026-08-27  
**Audited theorem SHA-256:**
`776b079858525b9c4e035925cd5f2df599a3f8062b84e612a3323c95d2a9d764`  
**Audited verifier SHA-256:**
`48d5f591fc0212553cfec87111189fd346265b7874bd0e3510e2cf6db30caff7`

**Disposition:**
`FAIL_PENDING_CUSTODY_AND_SCOPE_CORRECTIONS__FINITE_INCIDENCE_SCHUR_A3_SYMBOL_VOLUME_AND_SUPPLIED_ACOUSTIC_ACTION_MATHEMATICS_PASS__CURRENT_F3_DETUNING_COLLECTIVE_PHASE_PAIR_FIELD_LIFT_COMMON_CONE_AND_PHYSICAL_RGRL_B_REMAIN_OPEN`

## 1. Replay and exact finite algebra

The supplied verifier replays `33/33 PASS`.  Independent calculation confirms
the central finite identities.

1. Every column of `B_N` has four ones.  Two distinct columns have inner
   product one exactly when the two parents share one child, and otherwise
   have inner product zero.  A distinct pair cannot share two children.
   Therefore

   \[
   B_N^\dagger B_N=4I+A_N
   \]

   for every finite `N`, including the nonregular boundary of `S_N`.

2. With child amplitude `eta`, the lower block of `(H_N-zI)(psi,eta)=0`
   gives

   \[
   \eta={t\over\Delta-z}B_N\psi .
   \]

   Substitution into the upper block gives exactly

   \[
   \left[-zI-{t^2\over\Delta-z}B_N^\dagger B_N\right]\psi=0.
   \]

   Thus the sign, operator order, and denominator in (FD06) are correct.

3. Put `K=B_N^dagger B_N`.  Since every parent has degree four and every
   child has at most four parents, `||B_N||<=4` and `0<=K<=16I`.  Under the
   natural spectral identification of the low block with the parent space,
   its exact operator is

   \[
   f_{\Delta,t}(K)
   ={\Delta I-\sqrt{\Delta^2I+4t^2K}\over2}
   =-{t^2\over\Delta}K+{t^4\over\Delta^3}K^2
     +O_{\|\cdot\|}\!\left({t^6\over\Delta^5}\right).
   \]

   Consequently the order asserted in (FD07) is correct, and can actually be
   made uniform in `N`, but the theorem must define `H_(N,low)`, state that the
   remainder is in operator norm, and freeze a smallness domain.  As written,
   the perturbative formula is mathematically under-specified rather than
   algebraically wrong.

4. The three-mode spectrum (FD08) is correct for the separately truncated
   incidence row `B=(1,1)`.  It is not literally a reduction of the full
   `B_N` block: each full q4 parent has three additional children, and the
   two-parent principal block of `4I+A_N` is `[[4,1],[1,4]]`, not
   `[[1,1],[1,1]]`.  The text must call (FD08) a one-common-child comparator
   or per-child contribution.  It may not say that the isolated bright energy
   agrees with all of (FD07)'s displayed diagonal shift.

## 2. `A3`/FCC symbol, covolume, and refinement

These parts pass.

- The tetrahedral frame has unit diagonal Gram entries and off-diagonal
  entries `-1/3`.
- The six roots obey
  `sum_(a<b) alpha_ab alpha_ab^T=(16/3)I_V`.
- On the infinite interior diamond/FCC pair, the exact symbols satisfy
  `|b(k)|^2=4+A(k)` and
  `A(k)=12-(16/3)a_*^2|k|^2+O(a_*^4|k|^4)`.
- The scalar parent band is even at the symmetric `k=0` point and has no
  first-order Weyl/Dirac principal part there.
- The primitive `A3` cell covolume is
  `16 a_*^3/(3 sqrt(3))`; the tetrahedral-domain volume, its ratio
  `N^3/6`, and the exact stars-and-bars boundary corrections are correct.
- Under `a_* -> 0`, `N -> infinity`, and `Na_* -> L`, every fixed compact
  subset of the limiting tetrahedron's interior is a diverging number of
  lattice steps from the boundary.  The corrected condition
  `min_a m_a -> infinity` is the right lattice-unit statement.  The result is
  a refining flat mathematical affine atlas on the interior, not a physical
  spacetime or a global topology theorem.

The Bloch formulas are exact for the translation-invariant infinite interior
operator.  They are not the exact spectrum of a finite simplex `A_N`; the
current phrase “translation-invariant interior symbol” preserves that
distinction and should be retained.

## 3. Hydrodynamic action and raw-versus-collective cone

For the **supplied** complete quadratic action (FD15), Fourier transformation
does give

\[
 \omega^2={2\kappa\over\chi}\sum_{a<b}
 [1-\cos(a_*k\cdot\alpha_{ab})],
 \qquad
 c_*^2={16\kappa a_*^2\over3\chi}.
\]

For positive `chi,kappa`, its long-wavelength principal polynomial is
`-omega^2+c_*^2|k|^2`.  The distinction in (FD19) is exact and important: the
raw q4 append cone remains four-ray/polyhedral, while a separately supplied
collective wave operator can have a smooth quadratic characteristic cone.
There is no contradiction.

Two claim corrections are nevertheless required.

1. Stability alone does not produce (FD15).  The antecedent assumes a
   gapless, shift-symmetric real collective coordinate, a positive kinetic
   coefficient, and precisely the nearest-root stiffness action.  Neither
   the Schur kernel nor present F3 derives that collective coordinate,
   `chi`, `kappa`, or the absence of a mass term.  Lines 33--38 and the summary
   at lines 443--447 must not imply that the current scalar transfer itself
   has already supplied the physical stiffness.  What is proved is that its
   common-child algebra supplies the same `A3` support and that **if** the
   displayed massless collective action is generated, its principal cone is
   isotropic.

2. The geometric refinement does not by itself give a nondegenerate
   dynamical continuum.  If `a_* -> 0` while `kappa/chi` is fixed, then
   `c_* -> 0`.  A finite physical-speed refinement must additionally scale
   `kappa/chi=3c_phys^2/(16a_*^2)` (or equivalently rescale controller time).
   Under that scaling the dispersion correction is
   `O(c_phys^2 a_*^2 |k|^4)`, not a fixed-coefficient
   `O(a_*^4|k|^4)`.  This is an open scale/dynamics binding, not a defect in
   the fixed-lattice infrared expansion.

## 4. Six-sector screen

The abstract representation statement is correct: the six unordered edges of
a tetrahedron carry `A1 + E + T2`, with dimensions `1+2+3`, and the six edge
dyads span `Sym^2(V)`.  In the explicitly factorized quadratic class, positive
sector kinetic/stiffness coefficients give

\[
 \omega_X^2={m_X^2\over\chi_X}
 +{\kappa_X\over\chi_X}\lambda_{A3}(k).
\]

Thus the retained massless sectors have one common cone exactly when their
positive ratios `kappa_X/chi_X` agree.

The present wording does not yet lawfully bind that screen to q4/F3 physics.

1. PMMDC proves an invertible Jacobian for one supplied finite joint
   statistical pair family.  PMMDC explicitly does **not** identify its four
   binary ports with the four reusable `Q4-MERGE` operations and does not
   provide compactly supported pair-memory fields, gluing, propagation, or an
   F3 constitutive law.  Reusing its result here therefore requires an
   explicit `Q4-PAIR-FIELD-LIFT`/type join.  Until that join is derived, the
   six-component action is a prospective representation-theoretic screen,
   not “the q4 pair-memory fields” of the current parent.
2. `P_A,P_E,P_T` are not defined in the theorem, although the verifier builds
   one particular opposite-edge realization.  Define them as the exact
   orthogonal projectors for the edge representation (or cite the frozen
   inherited definitions).  Use the conventional irrep label `E`, not the
   unexplained `E_2` in (FD29).
3. State `chi_X>0` and `kappa_X>0` for every retained propagating sector.
   Otherwise equality of ratios does not establish a healthy hyperbolic cone.
4. The “exactly when” in (FD32) must either assume that all three sectors are
   retained or require equality only among the retained sectors.  If a later
   constraint removes a sector, its ratio is irrelevant to the common cone of
   the physical retained modes.

The theorem correctly leaves the constraint/Ward reduction, helicity-two
pole, universal stress vertex, and general orientation-coupled principal
symbol open.

## 5. Source custody and physical binding

The strongest custody guard in the packet is sound: `Q4-CARRIER-LIFT` is
explicitly new and open.  The bounded q4 witness earns count-front states and
blind retained word custody, not coexisting physical carrier sites or a
global coherent walk.  FOCST earns only a prospective local common-child
SWAP.  The packet does not silently erase that gap.

One additional load-bearing term is not yet owned.  Present F3 declares the
content-symmetric transfer amplitude `t`, but its carrier onsite term is
uniform.  The relative parent/child detuning `Delta` in (FD05) is not supplied
by that transfer and is essential to the low-energy Schur reduction.  It must
be added explicitly to the prospective interface, with its source-off
preparation and controller/work/boundary ownership.  Therefore lines 106 and
443--444 cannot currently attribute the whole low-energy kernel to “the same
scalar transfer declared in F3” under `Q4-CARRIER-LIFT` alone.

The remaining physical ceilings are correctly stated and must remain:

- controller time is not proper time without a common clock;
- one scalar acoustic cone is not a common matter/probe cone;
- coordinate covolume is not physical event volume without cell/scale binding;
- no current parent derives the collective phase, common-cone sector ratios,
  constraints, universal stress response, Einstein dynamics, `G`, or gravity.

The assertion that all ports “remain owned” in `Q4-CARRIER-LIFT` is an
antecedent obligation.  The symbolic F3 port slot does not itself discharge
that obligation.

## 6. Verifier ceiling

The `33/33` replay correctly checks sampled finite incidence identities,
tetrahedral/root moments, dyad rank, covolume and tetrahedron constants, the
three-mode spectrum, bulk symbol identity, quadratic convergence, and a
particular `1+2+3` projector decomposition.

It does not check the general Schur sign/order, define the low-block operator
or bound its remainder, distinguish the isolated three-mode truncation from
the full q4 incidence block, establish source custody for `Delta`, derive
(FD15), test dynamical refinement scaling, type-join PMMDC pair variables to
q4/F3 fields, prove physical ports, or establish gravity.  Its passing count
must not be used as a full promotion gate.

## 7. Required corrections before promotion

1. Add the relative child detuning and its complete physical ownership to the
   prospective carrier-lift antecedent; narrow every “declared F3 gives the
   kernel” sentence accordingly.
2. Define the exact low spectral block and operator-norm remainder/smallness
   domain in (FD07).
3. Relabel (FD08) as a separately truncated one-child comparator and restrict
   its agreement claim to the per-common-child contribution.
4. State explicitly that the massless collective action, not merely a stable
   phase or the one-carrier Schur kernel, is the antecedent of CCMAC-2; keep
   `chi,kappa` un-derived.  Add the finite-speed scaling condition if the
   acoustic dynamics is composed with the `a_* -> 0` refinement.
5. Add the pair-field type join, define the three projectors, require positive
   sector coefficients, and scope (FD32) to the actually retained sectors.

After those narrow corrections, the promotable result is:

`EXACT_FINITE_Q4_COMMON_CHILD_INCIDENCE__EXACT_PROSPECTIVE_DETUNED_F3_FORM_SCHUR_SIBLING_KERNEL__EXACT_A3_BULK_SECOND_MOMENT_COVOLUME_AND_REFINING_AFFINE_ATLAS__EXACT_SUPPLIED_MASSLESS_ACOUSTIC_ACTION_TO_SMOOTH_IR_PRINCIPAL_CONE__RAW_APPEND_CONE_REMAINS_POLYHEDRAL__CARRIER_LIFT_DETUNING_OWNERSHIP_COLLECTIVE_PHASE_PAIR_FIELD_JOIN_COMMON_PROBES_CONSTRAINTS_STRESS_RGRL_B_AND_GRAVITY_OPEN`

