# Typed clarification of the adopted U-DCL predicate

**Clarification ID:** `UDCL-TYPED-CLARIFICATION-V002`

**Date:** 2026-08-25

**Status:** faithful theorem-level clarification of the adopted V001 postulate;
no empirical confirmation and no enlargement of the physical domain

## 1. Why this clarification is required

The adopted decision stated D1--D4 in physical prose.  A hostile reread found
four places where an exact theorem could otherwise exploit an unintended weak
reading:

1. the physical existence of a DCL realization and the prospective evidence
   used to certify it were not separately typed;
2. the K/W alternatives in D2--D4 were not bound by one displayed branch
   quantifier;
3. ``directed'' did not expressly say that the external macro-incidence is
   acyclic after permitted bounded internal modules are contracted; and
4. a family indexed by pre-frontier histories could be mistaken for the one
   joint state that must contain the history register and its distribution.

This file removes those weak readings.  It does not add a new representation
category, relax a falsifier, shrink the record domain, infer a joint state from
marginals, or turn adoption into evidence.

## 2. Ontic predicate and evidential certificate

`DCL_phys(r)` is an ontic predicate.  It says that the actual episode has an
outcome-independent physical realization satisfying the typed clauses below.
The realization and its physical referents are fixed by the episode and its
law; they are not selected as a function of the observed record outcome,
contrast, or successful Coverage-U score.

`Cert_DCL(r;P)` is an evidential predicate.  It says that a protocol packet `P`
was prospectively frozen before scored-data access and soundly authenticates
one such physical realization.  Under the packet's soundness and custody
premises,

\[
 \operatorname{Cert}_{DCL}(r;P)\Longrightarrow DCL_{\rm phys}(r).
\]

Failure to obtain a certificate does **not** imply
\(\neg DCL_{\rm phys}(r)\).  This distinction is required because the sealed
`C`, `S`, and `J` predicates are ontic, while preregistration is evidential.

In this lane, the original symbol is retained as an alias:

\[
 DCL(r):=DCL_{\rm phys}(r).
\]

Thus the adopted global statement is read without changing its quantifier:

\[
 \operatorname{U\!DCL}\Longleftrightarrow
 \forall r\in\mathfrak R^{\rm actual,bf}_{\rm FM},
 \ DCL_{\rm phys}(r).
\]

## 3. One common typed realization

For a fixed record, a physical DCL witness has the form

\[
 \mathfrak d_r=
 (b,G,{\cal I},{\cal P},H_G,Z_G,\Omega^b_{x,z},
  \Phi^b_{<G},\Phi^b_{>G},{\cal J}^b,U^b),
 \qquad b\in\{K,W\}.                                  \tag{C1}
\]

The same `b`, frontier `G`, physical incidence \({\cal I}\), port census
\({\cal P}\), history register \(H_G\), separator object, maps, instruments,
and future-input object must satisfy all four clauses.  A classical separator
cannot be sewn to otherwise untyped quantum maps, or conversely.

Then

\[
 DCL_{\rm phys}(r)\Longleftrightarrow
 \exists\mathfrak d_r\;[D1(\mathfrak d_r)\land
 D2_b(\mathfrak d_r)\land D3_b(\mathfrak d_r)\land
 D4_b(\mathfrak d_r)].                                \tag{C2}
\]

The existential is over physically admissible realizations, not arbitrary
endpoint fits.

## 4. Exact readings of D1--D4

### D1 -- complete acyclic external frontier

After every bounded internal-feedback, continuum, or higher-order object
admitted by D3 is contracted to its actual exact ordinary K/W macro-map, the
complete **external** operational incidence of the registered ancestor mission
has an outcome-independent finite physical bound and is a directed acyclic
graph.  `G` is one
branch-independent topological frontier in that graph.  Every actual relevant
source/write occurrence and its complete registered ancestry is pre-`G`; every
registered query is post-`G`; every actual source-to-query route crosses `G`;
and the complete physical census contains no return or bypass.

An unresolved physical loop is not made acyclic by drawing a box around it.
Internal feedback is allowed only when the actual theory already supplies a
unique normalized/positive exact K/W macro-map for that bounded module.  This
is the finite-operational acyclic-frontier specialization of the sealed `C`
predicate.  It is an explicit wording restriction on the ambiguous V001 word
``directed'', not a new complement architecture.

### D2 -- one joint state containing history

For every positive-probability source arm `x` and fixed exogenous context `z`,
D2 supplies one normalized apparatus-inclusive separator object
\(\Omega^b_{x,z}\):

- in branch `K`, one normalized joint probability measure on a
  standard-Borel crossing space that includes the classical history register
  \(H_G\); or
- in branch `W`, one normalized positive state on one unital joint C*-algebra
  containing the declared classical history-register subalgebra and the
  component embeddings.

The object includes the actual history probabilities and all correlations with
the remaining separator degrees.  Conditional objects for positive-probability
histories are conditionals of this supplied joint object; they are not an
unrelated family and are not used to infer joint existence.  The object also
contains every crossing carrier, apparatus, controller, environment, reference,
boundary/charge/dressing datum, memory, and external port.  It is never
conditioned on a post-`G` registered outcome.

### D3 -- well-posed typed two-sided sewing

The maps, instruments, and separator in (C1) share branch `b` and matching
domains/codomains.  Their finite sequential/parallel contraction gives unique,
well-defined, normalized stochastic macro-maps in `K`, or UCP macro-maps with
complete CP instruments in `W`, on both sides of `G`, for every permitted
arm/context/history.  Identity, normalization, all adaptive and failure
outcomes, and every crossing memory are retained.  A continuum, feedback, or
higher-order object qualifies only through its actual exact contracted macro-map;
an endpoint-fitted table does not.

### D4 -- locality, provenance, and joint future inputs

All typed modules respect the same outcome-independently fixed physical
incidence and port census.
No remote free setting or old source label enters through an undeclared path.
After conditioning on the complete separator/history object and fixed context,
the **joint** law/state/process for all later fresh inputs is common to the
source arms; equality of separate marginals is insufficient.  Post-`G`
adaptation can depend on the source only through the separator/history wires.
Classical incomparable laboratories factor conditionally on the separator and
their local settings.  Quantum incomparable local algebras commute or are
explicitly tensor separated, and their localized maps obey causal
factorization.

D4 deliberately contains T5-grade provenance.  The formal result must not
advertise the common-future condition as physics obtained without this strong
custody premise.

## 5. Clarification versus strengthening

The common-branch quantifier, joint-history state, composed-map typing, and
ontic/evidential split clarify content already intended by “one existing
sufficient category,” “complete joint separator,” “exact two-sided K/W sewing,”
and the earlier audit's description of T5-grade provenance.

External macro-acyclicity and joint rather than marginal fresh-input
independence are stronger **wording** than V001's shortest sentences.  They are
not a newly chosen physical sector: they make explicit the exact finite-
operational `C` and common-future requirements that V001 claimed to instantiate.
Without those readings, the displayed implication to sealed `C,S,J` is not an
exact theorem.  This lane therefore uses the clarified predicate in every
theorem and labels the clarification openly.

The historical V001 decision and adoption record remain unchanged.  The
adoption record should receive a hash-pinned clarification note stating that
its symbol `DCL` is implemented as `DCL_phys` in this file; if the program
principal rejects either explicit wording restriction, the conditional theorem
must revert to pending rather than silently weakening this definition.
