# Distinct hostile audit — GL6AW anisotropic Følner twist closure

**Target:** `LANE_CROSS_RFT_GRA_GL6AW_ANISOTROPIC_FOLNER_TWIST_CLOSURE_V001/`  
**Frozen theorem SHA-256:** `ad315f2434b1042f183eeeb0244f9d323214717bebf01643c46eccd110373f9d`  
**Frozen author-manifest SHA-256:** `73b0055d4fca77c044282fe80c9de78e3976c9bc77b37c31037d75f8e1ecb234`  
**Frozen author-seal-file SHA-256:** `382b10d28301c536f007fe79241167a1a830b68ce47135dbd5dfca37850b93e0`  
**Disposition:** `PASS__EXACT_FINITE_CENTERED_SECTOR_DICHOTOMY__ANISOTROPIC_FOLNER_CLOSURE__GNS_AND_PHYSICAL_PROMOTION_EXPLICITLY_OPEN`

## 1. Independence and custody

All eleven author files are pinned in `AUDITED_TARGETS.sha256`.  The
independent replay imports no author module.  Both author verifiers and both
independent audit verifiers pass in normal and optimized Python modes.  All
eighteen pinned GL6AR, GL6AS, and GL6AU author/audit dependencies resolve, as
do their manifests and seals.

## 2. Rectangular geometry and nonempty centered sector

The audit independently reconstructs every elementary six-link loop on a
large family of rectangular quotients.  For periods at least four, the
`4V` loop terms are distinct; each loop has six distinct links, uses three
ports twice at opposite alternating parity, and has degree two on each of
its three parent and three child vertices.  Toggling an alternating loop
therefore preserves the lock.

For even `L1`, the union of the alternating port-`0`/port-`3` matching and
the constant port-`1` matching has degree two at every parent and child and

```text
N_0=V/2.
```

Thus the centered sector used by the theorem is constructively nonempty.
The witness does not select a controlling component, and the author does not
claim that it does.

## 3. Translation character and component dichotomy

With

```text
q=2pi/L1,
U_0=exp(iq sum_x x1 n_(x,0)),
```

the wrap-plane correction is an integer multiple of `L1`, so translation
along `x1` gives

```text
Y U_0 Y^-1=exp(-2pi i N_0/L1)U_0.
```

In the centered sector this phase is `exp(-pi i L0L2)=-1` exactly when the
transverse area is odd.  The sign convention for the direction of `Y` can
invert the displayed exponent, but cannot change the decisive character
`-1`.

The diagonal twist preserves each flip-component span.  If a controlling
component is translation stable, its unique positive PF ground vector has
translation eigenvalue `+1`; the twisted vector has eigenvalue `-1` and is
therefore orthogonal in the same connected component.  If the component is
not stable, translation maps it to a distinct orthogonal component with the
same controlling energy.  This proves the stated finite centered-sector
dichotomy without assuming flip-graph connectivity or component stability.

## 4. Exact term count and variational energy

For the four unordered port triples, the port-zero `x1` twist changes are

```text
{0,1,2}: +/-1,
{0,1,3}: +/-1,
{0,2,3}: 0,
{1,2,3}: no port zero.
```

Wrap changes `+/-(L1-1)` have the same unit phase with the opposite sign.
Because the `4V` loops are distinct, exactly `2V` Hamiltonian terms are
affected.

For an affected partial flip, direct two-configuration matrix algebra gives

```text
<Upsi,-J T_c Upsi>-<psi,-J T_c psi>
  =J[1-cos(q)] <psi,T_c psi>.
```

PF positivity and the partial-permutation norm imply
`0<=<T_c><=1`; equal amplitudes and graph regularity are not used.  Summing
the affected terms and applying min--max inside the stable component yields

```text
Delta_C <=2JV[1-cos(2pi/L1)]
        <=4pi^2J L0L2/L1.
```

The coefficient, factor of two, and direction-pairing survive independent
numeric and symbolic attacks.

## 5. Anisotropic Følner conclusion

For odd `m>=5`, the periods `(m,2m^3,m)` satisfy every parity and minimum-
period premise.  All three periods and the injectivity radius diverge, while

```text
boundary/volume = 2(1/m+1/(2m^3)+1/m) -> 0.
```

The associated boxes are therefore a genuine three-dimensional Følner/van
Hove sequence despite their diverging aspect ratio.  The stable-component
branch obeys

```text
Delta_C(m)<=2pi^2J/m -> 0.
```

The other branch is exact translation-related ground-component degeneracy.
The correct result is consequently an obstruction to a *unique uniformly
isolated* centered-sector ground, not a proof that the positive gap above a
possibly degenerate ground subspace must close.

## 6. Counterexample and scope attacks

The following stronger readings fail and are not made by the author:

1. With even transverse area the twist character is `+1`, so orthogonality
   is absent.
2. With isotropic side lengths the elementary upper bound scales as a length
   and does not close.  GL6AW proves no isotropic finite-size law.
3. A translation-related zero-energy or topological sector need not become
   an internal excitation of a selected pure infinite-volume phase.  The
   compatible-state/full-zero-energy-projection GNS bridge remains open.
4. The centered port sector need not contain the all-sector physical ground.
5. The theorem controls the exact pure order-six Hamiltonian only; it proves
   no stability under unknown higher-order finite-`h/U_d` corrections.
6. The large-twist translation character is not calibrated physical
   momentum.  No low-character pole, dispersion, cone, stress/Ricci law,
   gravity, or `G` follows.

The phrase "LSM-like" is descriptive only.  The proof uses the native lock,
port conservation, finite PF positivity, translation algebra, exact local
counting, and min--max; no conventional gauge, photon, or graviton theorem is
imported.

**Hostile verdict: PASS.**
