# Same-Sector Alpha Inheritance and Non-Embedding Theorem

**Claim class:** exact conditional effective-action theorem, joined to the
already proved RFT alpha-nonselection theorem

**Official short name:** SAI

**Not claimed:** a derivation of electromagnetism, a numerical prediction of
\(\alpha_{\rm obs}\), a proof that records choose a vacuum, or the destruction of
every possible foreign-sector object

## 1. The two quantifiers must not be confused

Freeze the complete comparison context

\[
 \chi=({\cal O},\text{eigenmode/current},\text{canonical normalization},
 \text{charge unit},{\cal S},\mu,\text{vacuum/background},
 \text{spectrum/mixing},\text{thresholds},\text{medium}).
\]

Let

\[
 \mathfrak M_\alpha({\cal A}_{\rm RFT};\chi)
 :=\{M\in\mathfrak M({\cal A}_{\rm RFT}):
 \alpha_\chi(M)\text{ is physically defined}\}.
\]

The existing nonselection theorem concerns a class of possible models:

\[
 \left|\left\{\alpha_\chi(M):
 M\in\mathfrak M_\alpha({\cal A}_{\rm RFT};\chi)
 \right\}\right|\ge2                                   \tag{1}
\]

It proves that the current portable RFT predicates do not themselves project
the admitted model class onto one electromagnetic coupling.  In particular, a
causally decoupled electromagnetic spectator can be changed without changing
any registered record law, and finite continuous strict formation margins
normally admit an open coupling neighborhood.

The present theorem asks a different, within-model question.  Once one parent
electromagnetic sector has been fixed, can two records made only from that same
sector possess independently different fundamental alphas?  Under the premises
below, the answer is no.

There is no contradiction between these statements:

\[
 \underbrace{\text{alpha not uniquely fixed across possible RFT models}}
 _{\text{nonselection}}
 \quad\land\quad
 \underbrace{\text{alpha common to all same-sector records in one model}}
 _{\text{inheritance}}.                                 \tag{2}
\]

## 2. Domain and load-bearing premises

All occurrences of “the same alpha” below are typed by \(\chi\).  Values quoted
at different contexts are first transported to a common context using a supplied
physical map.

### SAI1 — one connected selected phase

The domain \(D\) is one connected causal/gauge component of one realized vacuum
or phase with one unbroken compact \(U(1)\) connection \(A\).  A phase or domain
wall is not silently included in \(D\).

### SAI2 — one mediator

At the declared order there is no second \(U(1)\), kinetic-mixing matrix, hidden
photon, or independent connection used by a record subsystem.

### SAI3 — one common effective action

At reference scale \(\mu_0\) in scheme \({\cal S}\), the parent has

\[
 \Gamma_{\mu_0}=
 -{Z_F(\mu_0)\over4}\int_D F_{\mu\nu}F^{\mu\nu}
 +\sum_i n_iq_0(\mu_0)\int_D J_i^\mu A_\mu
 +\Gamma_{\rm rest},                                   \tag{3}
\]

where \(Z_F>0\), \(q_0\ne0\), and \(n_i\) belongs to one fixed charge lattice.
The coefficients are global action parameters on \(D\), not independently
assignable state variables.  No unlisted \(f(\phi,x)F^2\), varying-charge
spurion, or interface action is admitted.

### SAI4 — gauge custody

Gauge-preserving renormalization and the relevant Ward/BRST identities hold.
The integers or representation labels \(n_i\) specify charges under the same
generator; they are not independent gauge couplings.

Throughout any claimed RG/matching interval, the massless gauge eigenmode remains
uniquely identified, the fixed charge lattice and gauge custody persist, and the
retained pole/minimal-coupling sector retains the form (3).  Higher, nonminimal,
and process-dependent operators remain in \(\Gamma_{\rm rest}\).  Any mixing,
symmetry breaking, or phase transition is handled by an explicitly supplied
unique matching map; otherwise SAI-1 is asserted only at \(\mu_0\).

### SAI5 — typed base coupling

At fixed scheme and scale, after canonical pole/two-derivative normalization,
define

\[
 \alpha_D^{\cal S}(\mu)
 ={q_0^{\cal S}(\mu)^2\over4\pi Z_F^{\cal S}(\mu)}.       \tag{4}
\]

The ratio is invariant under gauge-field coordinate rescaling.  An alleged exact
duality is identified only when a supplied bijection maps the complete
observable, current, and charge-lattice data; an unspecified dual description is
not automatically quotiented.

### SAI6 — common RG ancestry

The complete coupling vector has one boundary value at \(\mu_0\), a unique RG
initial-value solution on the claimed interval, and one declared matching map at
each threshold.  A subsystem is not assigned a separately fitted boundary
constant.

### SAI7 — same-sector record

A record \(R\) in scope is a state, history, region, device, or subsystem
algebra obtained by a coefficient-preserving inclusion/restriction of (3), or by
an EFT equipped with a supplied exact matching map back to the same parent base
coupling.  Preparing, writing, holding, and reading \(R\) can change states,
sources, controls, and boundary data; those operations do not replace the parent
action coefficients.  An unmatched material or open-system effective response
is not \(\alpha_R\) in this theorem.

### SAI8 — aligned comparison

Two inferred couplings are declared different only after aligning the target
observable, physical eigenmode/current, scale, scheme, charge representation,
spectrum and mixing, vacuum/background, thresholds, medium response, and
estimator uncertainty.

## 3. SAI-1 — fixed-scale inheritance

Define the typed record class

\[
 {\cal R}^{\rm same}_\alpha(M;\chi)
 :=\{R:R\text{ satisfies SAI7 and its matched base-alpha target at }\chi
 \text{ is defined}\}.
\]

For every pair \(R_1,R_2\in{\cal R}^{\rm same}_\alpha(M;\chi)\) satisfying
SAI1–SAI8,

\[
 \boxed{
 \alpha_{R_1}^{\cal S}(\mu)
 =\alpha_D^{\cal S}(\mu)
 =\alpha_{R_2}^{\cal S}(\mu)}.                           \tag{5}
\]

Here \(\alpha_R\) means the ideal parent-matched base action coefficient governing
the declared record interaction.  It is not an intrinsic scalar stored inside a
record, a raw medium coupling, or a noisy finite estimator.  Data test statistical
compatibility with (5); a point estimate need not equal the ideal coefficient.

### Proof

Canonically normalize the one parent gauge field once:

\[
 A_\mu^c=\sqrt{Z_F}\,A_\mu,
 \qquad e={q_0\over\sqrt{Z_F}}.                          \tag{6}
\]

Equation (3) becomes

\[
 \Gamma_{\mu}=
 -{1\over4}\int_D(F^c)^2
 +\sum_i n_i e(\mu)\int_DJ_i\mathbin{\cdot}A^c
 +\Gamma_{\rm rest}.                                   \tag{7}
\]

There is one base minimal-vertex coefficient \(e(\mu)\) at the declared order.
A coefficient-preserving restriction cannot replace it, and an exactly matched
subsystem carries it through the supplied map.  For charged species \(i\) with
\(n_i\ne0\), the minimal one-photon vertex is \(e_i=n_ie\), so

\[
 {e_i^2\over4\pi n_i^2}=\alpha_D,                       \tag{8}
\]

For the common massless current pole/minimal vertex at this order, charge factors
are therefore fixed by \(n_i\).  Generic amplitudes can also contain loops,
thresholds, form factors, other couplings, and nonminimal operators.  Ward
identities preserve the common gauge generator and charge relations; those extra
terms are not independent zero-momentum base couplings.  Therefore every member
of \({\cal R}^{\rm same}_\alpha(M;\chi)\) inherits (4).  QED.

## 4. SAI-2 — trajectory inheritance

All same-sector records share one RG trajectory, more precisely one
scheme-equivalence class of trajectories.  Records sampling two scales may
correctly quote two numerical values, but both values lie on the one transported
curve fixed by the common parent.

### Proof

By SAI6 the complete RG initial-value problem has one solution from the one
boundary vector.  Equation (4) is a projection of that solution.  Unique matching
maps preserve uniqueness inductively across the declared thresholds.  A
subsystem restriction does not introduce another boundary condition.  QED.

## 5. SAI-3 — different-alpha non-embedding

Suppose a proposed device requires an invariant base coupling
\(\alpha'\ne\alpha_D\) after the complete SAI8 alignment.  The device cannot be
embedded as an ordinary SAI7 subsystem while retaining both its \(\alpha'\)
dynamics and SAI1–SAI6.

### Proof by contradiction

If the device were an SAI7 subsystem, SAI-1 would give
\(\alpha'=\alpha_D\).  This contradicts the aligned premise
\(\alpha'\ne\alpha_D\).  Therefore at least one premise has changed.  QED.

Representative ways a SAI premise can fail include:

1. another gauge connection or kinetically mixed sector;
2. a dynamical or prescribed spacetime-dependent gauge-kinetic field, with its
   gradients, background, and interface data;
3. another vacuum, phase, or disconnected domain;
4. a charge assignment outside the frozen lattice/normalization or a different
   gauge generator;
5. an explicit gauge-breaking or portal/interface interaction; or
6. a merely different scale, scheme, medium, charge factor, or estimator that
   ceases to differ after SAI8 alignment; or
7. an exact dual presentation related by a supplied full observable/current/
   charge-lattice bijection, which is physically equivalent rather than an
   inequivalent alpha.

The exhaustive exact conclusion is only that at least one of SAI1–SAI8 fails;
the list is not asserted to exhaust every enlarged theory.

Thus the host sector “forces alpha” in one exact sense: anything independently
established as a coefficient-preserving restriction or correctly matched EFT of
its fields and action is governed by its coupling.  Actual recordhood,
`REC(r)`, and `DCL_phys(r)` do not by themselves establish this same-visible-
\(U(1)\) ancestry.  If a hypothetical foreign structure is reconstructed
entirely from host fields, it uses the host coupling.  If it retains a genuinely
different aligned coupling, then either its same-sector assignment or a parent
premise is false.

This theorem does **not** say that such an enlarged object must disintegrate.
Its stability and interactions are model-dependent.

## 6. Combined theorem — free across models, fixed within a sector

The already proved nonselection witnesses and SAI together establish the
compatible pair

\[
 \left|\left\{\alpha_\chi(M):M\in
 \mathfrak M_\alpha({\cal A}_{\rm RFT};\chi)
 \right\}\right|\ge2,                                   \tag{9}
\]

for the witnessed model extensions.  For each fixed \(M\) satisfying SAI1–SAI8
and having nonempty \({\cal R}^{\rm same}_\alpha(M;\chi)\),

\[
 \left\{\alpha_R(\chi):
 R\in{\cal R}^{\rm same}_\alpha(M;\chi)\right\}
 =\{\alpha_M(\chi)\}.                                   \tag{10}
\]

Equation (9) is a model-class nonselection result.  It does not by itself prove
that every member is a complete universe containing atoms, chemistry, and
electromagnetic records.  Equation (10) is a within-parent inheritance theorem.
It does not select the parent.

### Active-EM strengthening

The companion canonical-\(U(1)\) cavity theorem goes beyond the spectator
witness.  At one fixed physical context, one conserved current pulse drives a
canonically normalized photon mode through \(eJ\cdot A\), with
\(e^2=4\pi\alpha\).  The same cavity, current, write duration, positive hold,
query, graph, and relevance floor are used in every parent.  It proves

\[
 D_{\rm TV}(\alpha)=1-e^{-4\pi\alpha B^2}.              \tag{11}
\]

The contrast is zero when the electromagnetic interaction is removed and
positive when it is present. Every event-to-query route crosses that vertex,
and the completed episode satisfies both `REC` and `FCLPD_W`. Provided the
frozen ideal-parent validity domain \(I_\chi\) contains them, two fixed parents
with

\[
 \alpha_1={9\over400\pi},
 \qquad
 \alpha_2={9\over100\pi}                               \tag{12}
\]

give exact contrasts \(1/2\) and \(15/16\) under identical controls.  Hence the
finite active-EM record-forming set obeys

\[
 \boxed{|{\cal A}_{\rm RF}^{\rm EM,toy}|\ge2.}          \tag{13}
\]

Thus present recordhood does not numerically select alpha even when the
alpha-bearing electromagnetic interaction is load-bearing rather than a
spectator.  The full proof and its ideal-current/detector ceilings are in
`ALTERNATIVE_RECORD_WORLD.md`.

Equation (13) is not a theorem that two complete universes with stable atoms,
chemistry, gravity, and self-consistent dynamical charged matter exist.  That
stronger set is separately denoted \({\cal A}_{\rm RF}^{\rm full}\).

## 7. AWAI — Actual-World Alpha Inheritance Theorem

This is an empirically anchored theorem about our universe, not a theorem about
all possible universes.

Let \(W_{\rm obs}\) denote the actual world and let \(\chi_0\) be a complete
reference context in which independent measurements anchor the visible
electromagnetic coupling and its tested universality to declared uncertainties.
Taking the ideal parent parameter underlying those measurements as an empirical
premise, freeze

\[
 \alpha_{W_{\rm obs}}(\chi_0)=\alpha_{\rm obs}(\chi_0),  \tag{14}
\]

where the equality can be replaced by the prospectively reported measurement
interval when experimental uncertainty is being scored.

Define \({\rm SAME\_VISIBLE\_U1}(r,W_{\rm obs};\chi)\) independently of alpha
agreement to mean that the alpha-sensitive dynamics of \(r\) are physically
traced to the same visible photon eigenmode and parent action as (14), with
mixing, background, thresholds, medium, representation, and matching aligned.
Define \({\rm ACTVIS}(r,W_{\rm obs})\) to mean all of the following:

1. \(r\) is an independently admitted actual bona-fide record in
   \(W_{\rm obs}\);
2. \({\rm SAME\_VISIBLE\_U1}(r,W_{\rm obs};\chi)\) holds; and
3. its subsystem dynamics satisfies SAI1–SAI8 rather than adding an undeclared
   sector, phase, modulus, or interface.

Let \({\cal T}_{\chi\leftarrow\chi_0}\) be the unique declared RG, threshold,
scheme, vacuum, and context transport.  Then

\[
 \boxed{
 {\rm ACTVIS}(r,W_{\rm obs})
 \Longrightarrow
 \alpha_r(\chi)=
 {\cal T}_{\chi\leftarrow\chi_0}
 \!\left[\alpha_{\rm obs}(\chi_0)\right].}              \tag{15}
\]

### Proof

By ACTVIS2–3, \(r\) is a same-sector subsystem of the measured visible parent.
SAI-1 gives equality at a common reference context; SAI-2 gives the unique
transport to \(\chi\).  Substitution of the empirical premise (14) gives (15).
QED.

Equation (15) is not circular because SAME_VISIBLE_U1 ancestry is established
without selecting records by their agreement with alpha.  Measurement identifies
the parent coefficient; the theorem extends it conditionally to every unmeasured
SAI-scoped same-sector record and forbids an independent record-level choice.
It is an empirically anchored physical proof of the observed alpha's
**actual-world necessity and universality for that visible-sector domain**.
Extrapolation beyond measured contexts remains conditional on sector
completeness and the SAI4–SAI6 transport premises.

Possible worlds with other alpha sectors do not weaken AWAI.  They show only
that the numerical value is contingent on, or selected by, the parent rather
than derived from the portable record axioms.

### Relation to URFT and `DCL_phys`

The adopted URFT domain quantifies over actual bona-fide finite-mission records,
and `DCL_phys(r)` requires an ontic physical realization.  On the visible-EM
branch, AWAI is therefore a necessary physical-instantiation compatibility law:
a packet requiring an aligned \(\alpha'\ne\alpha_{\rm obs}\) cannot be the
claimed ordinary visible-sector realization.

Bare `REC(r)` or `DCL_phys(r)` does not say that every record uses
electromagnetism or identify its gauge ancestry.  SAME_VISIBLE_U1 remains
necessary: a non-electromagnetic record has no visible alpha to infer, while a
dark-sector record has another gauge observable and a portal rather than a
different value of the same visible observable.

This establishes **that** nature's measured sector governs its same-sector
records.  What remains open is why this phase/value is realized rather than
another admissible one, and whether a cosmological or record-related process
established it.

## 8. The alpha `ALLOW/REQUIRE/SELECT` ladder

The results are most compactly typed at three different physical levels.

### Level A — finite-record allowance

Define

\[
 \operatorname{ALLOW}_{\rm RF}^{\rm EM,toy}
 (a;\Pi,\chi)
 \quad\Longleftrightarrow\quad
 a\in{\cal A}_{\rm RF}^{\rm EM,toy}(\Pi;\chi).          \tag{16}
\]

The active-EM cavity theorem proves that this predicate is true for at least two
distinct values, and for the exact validity-domain intersection stated in the
companion theorem. Thus present finite recordhood **allows** a set; it does not
**require** one theory-absolute numeral.

### Level B — host-sector requirement

Once a world \(W\) realizes one connected parent trajectory
\(\alpha_W(\chi)\), define

\[
 \operatorname{REQUIRE}_{W}^{\rm sameU1}(r,a;\chi)
 \quad\Longleftrightarrow\quad
 \text{every SAI-scoped realization of }r\text{ in }W
 \text{ has base coupling }a\text{ at }\chi.            \tag{17}
\]

SAI and AWAI prove

\[
 \boxed{
 \operatorname{ACTVIS}(r,W)
 \Longrightarrow
 \operatorname{REQUIRE}_{W}^{\rm sameU1}
 \bigl(r,\alpha_W(\chi);\chi\bigr).}                   \tag{18}
\]

This is the second-level structure: RFT can allow several parent settings while
ordinary embedding in one realized parent requires its one common setting. The
requirement is conditional on host ancestry; it is not a derivation of the
host's boundary value.

### Level C — complete-universe requirement

Define

\[
 \operatorname{ALLOW}_{\rm full}(a)
 \quad\Longleftrightarrow\quad
 a\in{\cal A}_{\rm RF}^{\rm full}.                     \tag{19}
\]

Then a theory-absolute numerical requirement would be the singleton theorem

\[
 \operatorname{REQUIRE}_{\rm full}(a_*)
 \quad\Longleftrightarrow\quad
 {\cal A}_{\rm RF}^{\rm full}=\{a_*\}.                 \tag{20}
\]

Equation (20) is open. A non-singleton full set would instead establish
contingency.

Finally, \(\operatorname{SELECT}_{\rm parent}(a_*)\)—a law describing why one
allowed parent setting became actual—is a fourth and independent proposition.
Neither finite-record allowance nor host-sector inheritance proves that such a
selector exists. The realized RG boundary condition may simply be a contingent
fact of the parent phase.

## 9. Three different alpha questions

1. **What value governs records in our realized visible sector? — closed by
   AWAI.**  Measurement anchors the parent value; every SAI-scoped actual record
   inherits it.  The true parent parameter need not be known to infinitely many
   decimal places for the symbolic equality and experimental interval statement
   to hold.
2. **Does present RFT recordhood select one alpha? — refuted.**  The canonical
   cavity construction proves \(|{\cal A}_{\rm RF}^{\rm EM,toy}|\ge2\) even with
   fixed controls and load-bearing electromagnetism.
3. **Can more than one complete matter/gravity universe support records? —
   open.**  The alternatives are
   \(|{\cal A}_{\rm RF}^{\rm full}|>1\), in which case alpha is contingent, or
   \({\cal A}_{\rm RF}^{\rm full}=\{\alpha_{\rm obs}\}\), in which case an
   exhaustive full-universe consistency theorem selects it.

A parameter-free numerical derivation is therefore a possible deeper outcome,
not a presumed missing mechanism.  Nature may simply realize one member of a
non-singleton full set.  Scientific laws routinely contain empirically fixed
constants; that possibility does not undo AWAI.

## 10. Exact proof ceiling

SAI closes:

- common same-sector inheritance after a parent sector is fixed;
- common RG ancestry rather than one scale-independent numeral;
- impossibility of a seamless, genuinely different-alpha same-sector
  subsystem; and
- compatibility of that result with RFT alpha nonselection; and
- actual-world necessity and universality of the empirically anchored visible
  alpha for ACTVIS records; and
- nonselection by finite recordhood even when a canonical electromagnetic write
  is load-bearing.

SAI does not close:

- emergence of compact \(U(1)\) from RFT;
- completeness of the actual gauge sector;
- viability of ordinary records at every counterfactual alpha;
- cardinality of the complete-universe set
  \({\cal A}_{\rm RF}^{\rm full}\);
- first-principles prediction of \(\alpha_{\rm obs}\);
- the physical process that established the actual phase; or
- the dynamics of a portal, phase wall, or foreign-sector object.

Nor does SAI prove `DCL_phys(r)`, `REC -> DCL_phys`, or U-DCL.  Conversely, the
adopted U-DCL postulate supplies neither \(U(1)\) nor alpha.  U-DCL ranges over all
actual bona-fide finite-mission records; SAI ranges only over the nonempty typed
same-visible-EM subclass.  Their proof obligations are compatible and neither
closes the other's physical antecedent.
