# Distinct hostile audit — GL6AQ authenticated E-loop selection obstruction

**Target:** `LANE_CROSS_RFT_GRA_GL6AQ_AUTHENTICATED_E_LOOP_SELECTION_OBSTRUCTION_V001/`  
**Frozen theorem SHA-256:** `1d1b01380ec8fd7ce83c69d45b68d9bde36bbe1dacdd32e3a5909ee6723a5ace`  
**Frozen author-manifest SHA-256:** `1adb601b5c0a957e9182ecfb89fd24fc2e63a318a7006759826911592d979afd`  
**Disposition:** `PASS__AUTHENTICATED_PAIR_QUERY_HAS_EXACT_NONZERO_LOCAL_LOCKED_E_OVERLAP__TRANSVERSE_K_HAS_ZERO_DIRECT_LOCKED_AND_ONE_CELL_S4_LINEAR_E_PROJECTION__SIX_RETAINED_SUPPORTS_GATE_MINUS63_OVER8_E_CHANGING_LOOP__TRACE_REFUTES_ONLY_UNIVERSAL_NONZERO_STATIONARY_CONTRAST__EXISTENTIAL_SELECTED_LOCKED_BULK_WITNESS_OPEN__NO_STATE_POLE_CONE_GRAVITY_OR_G`

## 1. Independence and custody

The GL6AQ author bytes were frozen before this review.  The ten exact target
files, including the author manifest, are pinned in `AUDITED_TARGETS.sha256`.
The author packet does not self-seal.

The independent replay imports no author module and performs 307976 checks
in normal and optimized Python modes.  The author replay separately passes
`1010/1010`; its packet/custody verifier passes `95/95`.  The pinned GL6AM
and GL6AN hostile-audit packets pass `54/54` and `58/58`.

All eleven direct dependency hashes resolve.  They are confined to the
frozen GL6AM and GL6AN author/audit chains, including both hostile-audit
seals.  No mutable reconnaissance result and no later effective-Hamiltonian
lane is a premise.

## 2. Authenticated pair query/read and exact locked E overlap

The replay reconstructs the unsigned four-port/six-pair incidence `R` over
the rationals.  It has rank four, so `ker R` is exactly two-dimensional.  On
the six local degree-two configurations the pair vector `M(z)` obeys

```text
R M(z)=(-1,-1,-1,-1).
```

The exact uniform locked covariance `C` satisfies

```text
P_E=(3/8)C,
P_E^2=P_E,
rank(P_E)=2,
C|_E=(8/3) I_E.
```

Therefore, for every real `c in E`,

```text
(1/6) sum_{k(z)=2} (c.M(z))^2 = (8/3)||c||^2.
```

This proves `P_x^(2) O_x(c) P_x^(2)` is nonzero whenever `c` is nonzero.
The six pair Pauli strings have the identity normalized-trace Gram matrix,
so the product trace separately gives `tr_0(O_x(c)^2)=||c||^2`.

The source/read language is correctly bounded.  `O_x(c)` is an already
authenticated finite real combination of GL6AM pair queries and pulses; it
is not identified with the retained-`K` Hamiltonian source and does not
select a state.

## 3. Direct and one-cell linear K obstruction

One transverse `X_e` flip changes degree two to degree one or three at each
endpoint.  Hence the exact finite-regulator identity

```text
P_Q X_e P_Q = 0
```

and its finite-word sum are sound.  Pair-`Z` reads and single-link-`X`
sources are distinct Pauli words, so their direct normalized-trace overlap
also vanishes.

The replay reconstructs all 24 `S4` actions.  The four-port permutation
representation has zero character overlap with the two-dimensional pair
`E` representation.  Reynolds averaging every elementary `6 x 4` map
confirms that each equivariant one-cell port-to-pair kernel satisfies

```text
P_E K(t)=0.
```

GL6AQ does not extend this one-cell statement to a generic spatial defect
word.  Such a word need not be `S4`-closed and receives no homogeneous
sector decomposition from GL6AM.

## 4. Six-support nonlinear loop gate

The audit exhausts all choices of two local loop ports, both alternating
orientations, and both compatible external occupations.  Every loop exchange
preserves degree two and has pair displacement

```text
R delta M=0,
delta M != 0,
||delta M||^2=16.
```

Thus the loop changes the local `E` coordinate, and the authenticated choice
`c=delta M` distinguishes its endpoints by exactly 16.

For the sealed alternating hexagon, every proper nonempty flip subset has
positive lock energy.  Independent enumeration of all `6!=720` orders gives

```text
-(63/8) h^6/U_d^5.
```

Every contributing order uses each of the six differing links once.  Giving
them retained coefficients therefore multiplies the result by

```text
product_{e in C} kappa_e.
```

The replay verifies this for independent rational coefficient vectors and
all 64 binary support words.  Removing any support kills this particular
leading off-diagonal entry.

This does not turn the pair query into a physical transverse source.  The
effect is a nonlinear six-support gate, while the direct and one-cell linear
projections remain zero.  GL6AQ also correctly withholds claims about
nonuniform lower-order diagonal shifts and the complete order-six bulk
generator.

## 5. Stationarity, trace obstruction, and quantifier scope

For a chosen stationary homogeneous GL6AM state, contraction with `c`
produces a positive correlation measure.  Its total mass is variance, but
positivity does not imply a nonzero retarded commutator or positive
dissipative measure.

The normalized product trace makes the distinction exact.  It is the unique
trace on the UHF spin algebra and is therefore invariant under the
homogeneous and every finite-defect automorphism.  Exact Pauli-basis replay
confirms tracial cyclicity locally.  Consequently it has positive pair-`E`
correlation mass for nonzero `c`, while

```text
tr_0([tau_t(A),B])=0,
tr_0(gamma_t^kappa(A)-gamma_t^kappa'(A))=0.
```

This is a valid counterexample to a claim that **every** lawful stationary
bulk state has nonzero response or nonzero defect contrast.  It is not a
locked-state witness and does not prove nonexistence of a useful selected
locked state.  GL6AQ states the quantifiers correctly: the universal nonzero
claim is refuted, while existence of some independently selected stationary
locked bulk state with nonzero contrast remains open.

The original homogeneous state is generally not stationary under a defect
dynamics, so GL6AM's positive homogeneous spectral measure and scalar
`A1/E/T2` split cannot be transferred to a generic defect comparison.

## 6. Promotion attacks

The audit attacked and rejected the following hidden promotions:

1. **Query/source confusion:** pair pulse/read access is not identified with
   the physical retained-`K` transverse perturbation.
2. **Finite/bulk confusion:** `P_Q` is a finite regulator projector, not an
   infinite locked-sector state or projector.
3. **Linear/nonlinear confusion:** the order-six product gate is not called a
   nonzero linear `K -> E` projection.
4. **Universal/existential confusion:** the trace refutes only a universal
   nonzero assertion; it does not refute existence of a selected locked
   witness.
5. **Stationary/nonequilibrium confusion:** homogeneous spectral positivity
   is not assigned to a generic defect contrast.
6. **Lower-order promotion:** nonuniform diagonal terms are left open.
7. **Gravity promotion:** no state selection, pole, physical momentum, cone,
   stress, Ricci/Einstein response, gravity, or `G` is inferred.

## 7. Verdict

GL6AQ is sound at the pinned snapshot.  It proves a real authenticated local
`E` read, an exact direct/linear obstruction for the physical transverse
`K` source, and the first sealed nonlinear bridge between retained supports
and an `E`-changing collective operation.  It also proves that these facts
alone cannot guarantee a nonzero stationary bulk response.  The remaining
gate is an independently justified stationary locked-state selection and a
controlled thermodynamic effective dynamics.

**Hostile verdict: PASS.**
