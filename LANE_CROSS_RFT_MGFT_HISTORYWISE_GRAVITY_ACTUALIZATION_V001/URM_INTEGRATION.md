# URM integration contract

## Purpose

The Universal Record Model may expose this lane only as a sealed **formal
discriminant**. It is not a gravity solver, an outcome selector, a GARH-D
admission engine, or a source of empirical proof.

The public surface is deliberately zero-input:

```python
URM.historywise_gravity_discriminant()
URM.historywise_gravity_discriminant_certificate()
```

No caller-supplied Boolean, manifest, lane root, physical packet, or claim label
can change the result. Scientific weight from caller input is therefore exactly
zero.

## Positive formal content

The certificate may report only these sealed results:

- HGA1: exact finite-group orbit nonselection;
- HGA1 feedback: exact equivariant fixed-point orbit closure;
- HGA1a: a unique deterministic equivariant flow maps symmetry-fixed complete
  input to a symmetry-fixed history;
- HGA1b: a unique equivariant fixed probability law on a finite transitive
  history orbit is invariant and uniform, but does not supply a realized
  sample;
- HGA2: conditional finite affine mean-field nonselection;
- HGA3: exact transitive finite-G-set stabilizer criterion
  (K\subseteq L), up to conjugacy; and
- the two-cell negative and orienting-input witnesses as exact finite algebraic
  examples, not physical GR.

The discriminator may state that endogenous equivariant data cannot select one
member of a fixed-point-free history orbit. It may also state that a nominated
input orbit has the mathematical transformation capacity for a covariant
deterministic map precisely when HGA3 holds. It may not infer that nature
supplies such an input.

## Fixed scientific ceiling

Every certificate must retain all of the following statuses:

```text
claim_class = FORMAL_FINITE_GROUP_DISCRIMINANT_ONLY
physical_GARH_D = NOT_ADMITTED_BY_THIS_DISCRIMINANT
GARH_Q = NOT_DERIVED_NOT_FORCED_BY_THIS_DISCRIMINANT
GARH_D_Q_DECISION = NOT_MADE
objective_actualization = OPEN_IN_THIS_LANE
physical_gravity = NO_PROOF_OUTPUT
record_causes_gravity = NO_PROOF_OUTPUT
Born_law = NO_PROOF_OUTPUT
general_relativity = NO_PROOF_OUTPUT
empirical_validation = NONE
caller_input_scientific_weight = ZERO
```

All claim-authorization flags are fixed false. In particular:

- failure of a GARH-D proposal cannot promote GARH-Q;
- the positive boundary table cannot promote its algebraic vector to a physical
  orienting field;
- the rational Lorentzian-signature proxy cannot promote itself to a physical
  metric or general relativity; and
- a successful custody or executable check cannot promote a formal theorem to
  empirical validation.

## Custody and refusal behavior

The URM implementation must pin this lane's closed manifest and independently
accepted audit. It must verify the manifest, every listed artifact, the exact
64/64 transcript and verdict, and a fresh execution of the sealed verifier.
The general finite-group proof remains a human-checkable analytic proof; the
64 checks reproduce finite witnesses and claim ceilings and do not enumerate
all finite groups.

The public gateway must refuse if:

- the lane, manifest, accepted audit, verifier, or verification transcript is
  absent, symlinked, malformed, incomplete, or hash mismatched;
- fresh verifier execution differs from the sealed transcript or verdict;
- any positional or keyword argument is supplied; or
- a caller attempts to substitute another lane, root, gate list, or physical
  packet.

The dedicated validator must test custody failure, zero-argument signatures,
direct URM delegation, immutable fresh certificate output, false authorization
flags, nonpromotion across both the positive and negative boundaries, and
propagation of any failure into the integrated URM conjunction.

## Integration boundary

This gateway is additive. It must not broaden or reinterpret the existing
world-observation, formation-input, gamma-flow, or proof-frontier contracts. Its
only job is to make the sealed theorem and its exact claim ceiling inspectable
through the URM.
