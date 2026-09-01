# Independent post-repair hostile audit — GL6AY

**Frozen repaired theorem SHA-256:**
`5b86ab5eb2998eb719dffd09e05add131863fd2a3290d87fb749dc8aebc1891c`

**Frozen repaired author-manifest SHA-256:**
`e81ec1cfd4bdcdc43b4709b8f90f9eceac3dfba82be80701dd4a2a7e08de089b`

**Frozen repaired author-seal-file SHA-256:**
`740f051b3347d7387e481a9991f536bd61a7e47ad51d80b680475dff394e5cbb`

**Disposition:**
`PASS__FINITE_COUPLING_PRETHERMAL_BRIDGE__LOCAL_COLLAR_TYPED__GLOBAL_PROJECTOR_CLAIMS_REMOVED`

No author byte was modified by this audit.  The original fail audit remains
frozen at the hashes in `PRIOR_AUDIT_CUSTODY.sha256`.

## 1. Custody and source replay

All twelve repaired author files are pinned in `AUDITED_TARGETS.sha256`.
The author manifest and seal resolve.  Both author verifiers pass in normal
and optimized modes.  The independent replay imports no author module.

The exact external source mapping was replayed against:

- ADHH, arXiv `1509.05386v3`, Theorems 3.1--3.3;
- Else--Fendley--Kemp--Nayak, arXiv `1704.08703v2`, Appendix A; and
- Bravyi--DiVincenzo--Loss, arXiv `1105.0675v1`, Section 4.

The primary constants, remainder factor, local-observable horizon, strong-
support extension, and whole-band proof boundary remain correctly stated.

## 2. F3 strong-support mapping — pass

Each coarse cell owns four link qubits and has dimension sixteen.  The
commuting constraints `q_v^2` have spectrum `{0,1,4}`, so `N_def` is integer
spectrum.  A link flip is strongly supported on the complete four-cell child
star

```text
S_(x,a)={x+d_a-d_b:b=0,1,2,3}.
```

That star also contains the parent cell; all other constraints commute with
the flip.  The set is connected and a fixed cell lies in exactly sixteen
translated labeled flip supports.  The potential-norm constants `16` and
`32` therefore survive.

The exact pinching retains charged-sector resonances.  Only finite-volume
`P_LD_0P_L=0` is asserted.

## 3. Finite local collar — pass

For each strongly supported term the repair defines

```text
N_S=sum_(v:supp(q_v^2) subset S)q_v^2,
P_S^0=chi(N_S=0),
Phi(S)=P_S^0D_hat(S)P_S^0.                               (RA.1)
```

This is sufficient.  Strong support gives individual commutation with every
constraint whose full star is not contained in `S`.  Termwise `N_def`
pinching gives `[D_hat(S),N_def]=0`.  Splitting the formal sum into the finite
contained part and individually commuting outside part gives

```text
[D_hat(S),N_S]=0.                                       (RA.2)
```

Every summand of `N_S` is nonnegative and commuting, so its zero projection
is exactly the product of the individual zero projections in (RA.1).
Starting from a globally locked configuration, `D_hat(S)` preserves the
contained zero sector by (RA.2) and preserves every noncontained constraint
individually by strong support.  Therefore `Phi(S)` agrees with `D_hat(S)`
on that input and maps it to globally locked output without any global
infinite-volume projector.

`P_S^0` is supported in `S`, commutes with every noncontained constraint, and
does not enlarge the strong support.  Projection is contractive, so the
exponential potential norm survives.

The finite-volume compression is correctly limited to

```text
P_LD_hat_LP_L=P_L[sum_S Phi_L(S)]P_L.                   (RA.3)
```

No `P_L` is promoted to an infinite-volume algebra element.

## 4. Locked-endpoint port conservation — pass

Every matrix element of a finite `Phi(S)` between globally locked
configurations has a finite, contractible symmetric difference.  GL6AX then
forces zero change in all four port totals.  The author now states exactly
this locked-endpoint/local-locked-algebra property and explicitly denies a
full-Hilbert port symmetry of each collar term.

Wrapping supports on a torus remain outside that implication and enter the
explicit tail `T_L`.  Exponential strong-support decay controls the tail.

## 5. Second twist moment — pass

The local collar and twist generator are both diagonal in link occupation,
so `[A_S,P_S^0]=0`.  Consequently

```text
[A_S,[A_S,Phi(S)]]
 =P_S^0[A_S,[A_S,D_hat(S)]]P_S^0.                       (RA.4)
```

This is an exact identity, not only a norm estimate.  With
`||A_S||<=C_geo|S|^2`, projection contractivity gives the displayed quartic
bound.  The exponential strong-support norm dominates that polynomial, so
`D_2` is uniform in volume at every fixed theorem-admissible finite ratio.
The quasi-local GL6AX anisotropic dichotomy therefore applies to the locally
typed locked interaction, with its original centered-sector and GNS
ceilings.

## 6. Finite-order contact and finite-horizon dynamics — pass

The finite-volume identity (RA.3) makes the GL6AO order-six coefficient
comparison well typed.  Lower locked-sector orders are scalar, so the first
nonscalar coefficient is invariant under the allowed near-identity
within-`P_L` gauge.  No convergence claim is added.

ADHH Theorem 3.3 is reproduced with

```text
0<r_1<ln(3/2)/4,
||tau_t^H(O)-tau_t^(U_dN_def+D_hat)(O)||<=K_3(O)/U_d,
t<=exp(r_1n_*).
```

The repaired text restricts the dressing conclusion to declared local
potentials and fixed local observables.  It explicitly states that neither
the source nor AY.20 controls

```text
||Y_L^*P_LY_L-P_L||
```

uniformly in volume and that no global dressed spectral subspace exists in
the infinite quasi-local algebra.

## 7. Leakage, winding, and whole-band boundary — pass

All projection symbols in these sections are now finite-volume `P_L,Q_L`.
Local `P_L->Q_L` leakage requires no winding.  A port-changing return to
`P_L` needs a noncontractible symmetric difference of at least `2L_min`
links.

For the simple winding row of length `r=2L_j`, every proper nonempty subset
has positive even boundary.  No intermediate return to `P_L` occurs, folded
terms cannot enter at the first endpoint-changing order, and every direct
word has the same sign.  The nonzero coefficient and factorial bounds are
correct.

The unperturbed `P_L`--`Q_L` gap is `2U_d`, while
`||W_L||=h|E_L|`.  The repaired Feshbach resolvent is explicitly

```text
Q_L(E-Q_LH_LQ_L)^(-1)Q_L.
```

The resulting volume collapse is used only as the boundary of the standard
global whole-band proof route, not as an impossibility theorem.

## 8. Ceiling attack — pass

The repaired packet does not claim:

- `V_hat=0` or an exact invariant locked space for the full Hamiltonian;
- uniform all-orders Schrieffer--Wolff/Kato/Feshbach convergence;
- a volume-uniform global dressed-projector norm;
- centered-sector selection or selected-GNS gaplessness;
- an isotropic physical mode or common cone;
- record authentication of the effective mode;
- gravity, an Einstein law, or numerical `G`.

The promoted result is exactly a small-but-finite-coupling, volume-uniform
prethermal normal form plus a well-typed local locked interaction and a
conditional finite-observation-horizon corollary.

**Hostile verdict: PASS.**
