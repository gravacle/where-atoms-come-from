# Distinct hostile audit — GL6AS native hexagon collective response

**Target:** `LANE_CROSS_RFT_GRA_GL6AS_NATIVE_HEXAGON_COLLECTIVE_RESPONSE_V001/`  
**Frozen theorem SHA-256:** `bfe36071a24ccc7d6d7a16afeeea1b5554a95562ae91ac59c709db478000db9f`  
**Frozen author-manifest SHA-256:** `0b7c12e51ff2892cc44e0e4e39d68b6939c5fdb47fce8c018db1ad09191e0f0f`  
**Frozen author-seal-file SHA-256:** `9e7f058f68e989a00e16a17611246392fb4cdc9957a24a0f770ac530d029a957`  
**Disposition:** `PASS__PURE_HEXAGON_PORT_T2_CONSERVATION_AND_CYCLE_CONTINUITY_EXACT__SMA_REQUIRES_STRUCTURE_FACTOR__HARMONIC_TWO_CHARACTER_MODES_CONDITIONAL__PAIR_E_IS_COMPLEMENT_EVEN_QUADRATIC_WITH_CONDITIONAL_TWO_T2_CHANNEL__RETAINED_SOURCE_A1_PLUS_T2_ZERO_CHARACTER_E_CROSS__SYM2_TRACELESS_T2_EQUALS_E_PLUS_T2_ALGEBRA_ONLY__NO_POLE_PHYSICAL_CONE_STRESS_GRAVITY_OR_G`

## 1. Independence and custody

The GL6AS author bytes were frozen before this audit.  All eleven author
files, including the manifest and seal, are pinned in
`AUDITED_TARGETS.sha256`.  The audit imports no author Python module.

The independent replay passes `11059/11059` checks in normal and optimized
Python modes.  The author replay separately passes `1704/1704`, and the
author packet verifier passes `125/125`.  All seventeen dependency hashes
resolve: six objects each from the sealed GL6AO and GL6AP author/audit
chains, and the five actually existing GL6AQ author/audit objects.  The
missing GL6AQ author seal is not invented.

## 2. Exact port conservation and its strict model scope

The replay reconstructs all 256 elementary hexagons on the declared `Q_4`.
Every hexagon contains three port labels twice, once on each alternating
half.  Toggling an alternating cycle therefore removes and inserts exactly
one occupation for each used port.  Consequently

```text
[H_hex,N_a]=0,
sum_a N_a=2|X|,
```

and the three centered totals transform as the centered four-port `T2`.
This is exact for the declared pure projected hexagon Hamiltonian.  GL6AS
correctly does not assign it to the unprojected one-link finite-`h`
Hamiltonian or claim that the all-orders expansion truncates.

The complement automorphism reverses the centered density and preserves
both the locked configurations and every hexagon toggle.  Thus the density
channel is complement odd while the pair read and loop source are
complement even.  Any spectral parity claim still requires an invariant
state or a definite-parity vector; GL6AS includes that premise.

## 3. Cycle complex, continuity, and normalization

For every exact nonconstant rational four-tuple `z`, the independent replay
finds

```text
B(z) C(z)=0,
rank B(z)=2,
rank C(z)=2.
```

Therefore `im C(z)=ker B(z)` at every nontrivial character.  Direct
reconstruction from the six parent-cell link changes reproduces the column

```text
(z_b-z_c, z_c-z_a, z_a-z_b, 0)^T
```

up to one common nonzero factor, which is a unit phase on the actual
character torus as declared in the theorem.  This also establishes that the
operator continuity equation is the boundary map of the actual local move,
not an inserted gauge condition.

For centered real `theta`, the exact leading symbol obeys

```text
C_1 C_1^T = 4[(theta.theta) P_T2 - theta theta^T].
```

Hence its restriction to the two-dimensional plane orthogonal to both
`1` and `theta` is `4(theta.theta)I`.  The audit replays the factor four over
all nonzero centered integer vectors in a bounded exhaustive family and
checks the normalization directly against the small-character exponential.
At the trivial character `C(1)=0`; fixed `Q_4` itself has no infrared limit.

## 4. Oscillator strength proves a diagnostic, not a mode

The two-state double-commutator replay gives the exact per-transition factor

```text
f=(J <tau>/2)|delta rho|^2.
```

Summing translations and orientations therefore yields GL6AS's
`(J t_hex/2) u^* C C^* u`, and the leading transverse value is
`2J t_hex I_2(theta)`.  The single-mode quotient is only the
energy-weighted average

```text
Delta_T2 <= f/S_T2^+.
```

The positive-frequency structure factor is indispensable.  The audit gives
an explicit algebraic counterexample in which both `f` and `S^+` vanish
quadratically while their ratio remains a fixed nonzero gap.  Thus exact
conservation and `f=O(theta^2)` do not prove a gapless mode or pole.

The phrase “only charge-supported candidate” is correctly confined to the
nontrivial link/pair sectors compared.  GL6AS expressly excludes the always
conserved scalar energy density and possible nonlocal homology labels from
that quantifier.

## 5. Conditional harmonic result and character ceiling

Given the additional projected canonical variables and positive coherent
Hessians, the linear equations have squared frequencies equal to the two
nonzero eigenvalues of

```text
C G_0 C^* K_0.
```

Under the further isotropic conditions `K_0=kappa P_T2` and `G_0=gI`, the
replay independently recovers

```text
omega_1^2=omega_2^2=4g kappa I_2(theta)+O(|theta|^3).
```

The author correctly withholds this result when the coherent state,
canonical pair, or either positive Hessian is absent.  In particular, the
pure order-six operator supplies neither `kappa` nor a normalization that
would set `g=J`.

The variable `theta` is only an additive translation character.  GL6AS
does not call it physical momentum and requires an embedding, length/time
calibration, scaling limit, and cross-sector agreement before any physical
velocity or cone could be asserted.

## 6. Pair-E composite, parity, and continuum scope

On each strict degree-two node,

```text
Z_a=-2e_a,
O(c)=4 sum_(a<b) c_ab e_a e_b,
c in ker R=E.
```

The replay reconstructs `rank R=4`, the two-dimensional `E` kernel, and

```text
(1/6) sum_(k=2) O(c)^2 = (8/3)c.c.
```

This is an exact locked mean-square operator overlap, not a pole residue or
universal dynamical response.  Complement sends `e` to `-e` and leaves
`O(c)` invariant.  Accordingly the one-density odd `T2` matrix element
vanishes in a complement-parity state.  If complement symmetry is absent,
the zero-character rule still follows in the retained `S4`-invariant setting
from `Hom_S4(T2,E)=0`; a small-character `E` amplitude is symmetry allowed.

Under the explicitly stronger Gaussian premise, direct expansion of the
quadratic operator gives the two-mode form factor

```text
F_c^(lambda mu)=4 sum_(a<b)c_ab
  [r_(a lambda)r_(b mu)+r_(b lambda)r_(a mu)].
```

The `1/2` in the spectral measure is an identical-particle counting
convention.  The linear `v|k|` and quadratic `Dk^2/2` expressions are merely
kinematic continuum minima in character coordinates.  They become support
thresholds only when the form factor remains nonzero near the minimizers;
they are never promoted to isolated poles.

The isolated flippable-hexagon doublet also replays exactly:

```text
|<+|O(delta M)|->|=8,
E_- - E_+=2J.
```

The audit confirms that this block is not invariant in the full lattice and
therefore supplies no bulk pole or dispersion.

## 7. Retained source and exact zero-character selection

The retained sixth-order coefficient is the six-support product
`J product_e kappa_e`.  Linearizing uniform port perturbations gives the
four-by-four map with eigenvalue `6` on `A1` and `-2` on centered `T2`.
Thus the four unoriented loop-source components transform as `A1+T2` and
have no `E` projection at zero character in a translation- and
`S4`-invariant state:

```text
K_(E<-loop)(omega,1)=0.
```

The scalar zero-character source is proportional to `H_hex` itself and has
no finite-frequency off-diagonal energy-basis matrix elements.  At small
character, `T2 x T2` contains one `E`, while a scalar first reaches `E`
quadratically.  GL6AS correctly says these coefficients are allowed, not
forced to be nonzero.  Complement parity forbids an even loop source from a
one-density odd channel; its Gaussian two-mode interpretation remains
conditional on the separately assumed coherent construction.

## 8. Algebraic tensor channel is not gravity

The full `S4` character replay gives

```text
Sym^2(T2)=A1+E+T2,
Sym^2_0(T2)=E+T2.
```

This establishes only a five-component finite-group composite channel.  At
one strict locked node the authenticated pair read supplies its `E` part;
the local pair-`T2` part vanishes by the three opposite-pair equalities.  A
nonzero companion `T2` therefore requires bilocal, derivative, or
coarse-grained data that has not been authenticated here.

The representation coincidence supplies no `SO(3)` enhancement, common
pole, residue or speed, helicity reduction, stress conservation, Ward
identity, Ricci/Einstein response, gravity, or `G`.  A conditional linear
two-particle threshold is neither a graviton nor a physical cone.

## 9. State and universal-response attack

The normalized product trace has positive pair-`E` mean-square weight and,
by tracial cyclicity, zero expectation for every retarded commutator.  The
independent replay checks cyclicity on a complete local Pauli-word basis.
This is a lawful counterexample to a universal nonzero-response inference.
It is not a selected locked ground-state witness and does not refute the
existence of another state with useful response.  GL6AS preserves these
quantifiers and makes no trace-to-ground-state promotion.

## 10. Verdict

GL6AS is sound at the pinned snapshot.  It identifies the exact native
charge/continuity variable of the declared pure hexagon model, proves its
rank-two cycle geometry, and determines the exact representation routes by
which pair and retained reads can access composite channels.  It also proves
why none of those algebraic facts alone supplies a mode, pole, physical cone,
stress law, or gravity.

**Hostile verdict: PASS.**
