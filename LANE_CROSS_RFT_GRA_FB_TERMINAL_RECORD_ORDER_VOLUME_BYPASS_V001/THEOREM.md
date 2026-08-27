# Terminal-record operational order-and-volume bypass theorem

**Lane ID:** `CROSS-RFT-GRA-FB-TERMINAL-RECORD-ORDER-VOLUME-BYPASS-V001`

**Official short name:** `TROV`

**Date:** 2026-08-27

**Status:** `BOUNDED_PROOF_DRAFT__INDEPENDENT_AUDIT_PENDING`

**Claim class:** exact URFT/TRL-to-directed-influence lemma; exact
cross-mission union/preorder theorem and countermodels; sharply conditional
Hawking--King--McCarthy/Malament plus four-volume reconstruction; exact premise
comparison against the Fisher/FERS route

**Not claimed:** that the union of recorded influences is already a partial
order; that absence of a positive signal proves spacelike separation; that an
arbitrary order is manifold-like; that event count is metric four-volume; that
causal reconstruction supplies Einstein dynamics; or that records cause a
metric merely because they reveal causal influence

## 1. Proposed bypass

RMLB shows that AURFT/U-DCL plus terminal-record localization does not force
one physically readable classical query.  This lane tests a route that does
not ask for one.  Allow incompatible instruments and define geometry from the
union of all prospectively qualified positive terminal-record influence
missions:

\[
 \text{qualified record influence}
 \longrightarrow \text{operational order}
 \longrightarrow \text{causal conformal class}
 \xrightarrow{\text{physical four-volume}}g.       \tag{FB01}
\]

This can bypass Fisher-query soldering only if the first arrow supplies a
complete global chronology rather than a context-dependent collection of
finite arrows.

## 2. Nonmetric event and mission definitions

Let `E` be a set of event identities fixed by apparatus/lineage custody, not by
coordinates, distances, light cones, or a metric.  A mission
`mu=(e,f,a,a',Q_f,z)` is **positive and qualified** when:

1. `e,f in E`, with `e != f`, are respectively the authenticated
   intervention-setting and terminal-query events under one frozen
   cross-mission identity rule;
2. `do_e(a)` versus `do_e(a')` is randomized or otherwise causally identified;
3. `Q_f` is the complete registered query of an independently bona-fide
   terminal record;
4. one actual `DCL_phys` realization types `e` as a pre-frontier
   source/controller port and `f` as its registered post-frontier query, retains
   the complete source, controller, environment, boundary, failure, and bypass
   census, and forbids a direct source-label side input to that query; and
5. the frozen lower margin is positive:

\[
 \Delta_\mu
 :=D_{\rm TV}\!\left(
 P(Q_f\mid do_e(a),z),P(Q_f\mid do_e(a'),z)
 \right)-\epsilon_\mu>0.                            \tag{FB02}
\]

Use trace/diamond process distinction in a quantum packet; only strict
nonidentity is needed.  Define

\[
 e\leadsto_R f
 \quad\Longleftrightarrow\quad
 \text{at least one positive qualified mission }\mu:e\to f, 
                                                               \tag{FB03}
\]

and let `preceq_R` be its reflexive-transitive closure.  When this closure is
antisymmetric, write

\[
 e\prec_R f\quad\Longleftrightarrow\quad
 e\preceq_R f\ \text{ and }e\ne f.                 \tag{FB03a}
\]

This typing matters: `preceq_R` is reflexive, whereas a chronological relation
`I+` is irreflexive.

The definition uses only event custody, interventions, record queries, and
their DCL incidence.  It does not assume a continuum or metric.

## 3. Theorem TROV-1 -- exact record-influence consequence

Assume AURFT/U-DCL and the terminal-record premise for every mission admitted
to (FB03).  Then:

1. every arrow `e leadsto_R f` is a same-parent directed operational influence
   witnessed by a bona-fide record, not a symmetric correlation;
2. incompatible instruments can contribute different arrows without a joint
   outcome distribution or one common classical POVM;
3. `preceq_R` is the unique minimal preorder containing all witnessed arrows;
   and
4. `preceq_R` is a partial order exactly when the union graph has no directed
   cycle between distinct event identities.

### Proof

For a positive mission, the intervention changes the terminal query law by
(FB02).  D1 places the source/intervention before the registered query in one
finite acyclic external incidence with no return or bypass.  D2 supplies the
joint apparatus-inclusive state, D3 supplies the exact complete instrument,
and D4 prevents the query kernel from receiving the source label except through
the declared physical wires.  Under causal factorization, if no directed
physical route from the authenticated intervention port to the query existed,
the `do`-change could not change the complete query law.  Hence the positive
difference is directed influence along that mission's physical incidence.
This proves item 1.

Each mission has its own complete instrument inside the allowed U-DCL family;
no joint distribution over counterfactual outcomes is used, proving item 2.
The reflexive-transitive closure of a relation is, by definition, the
intersection of all preorders containing it, proving item 3.  Such a closure
is antisymmetric exactly when no two distinct vertices lie in one directed
cycle, proving item 4. QED.

This is the exact URFT/TRL-to-operational-order lemma.  It earns the arrows and
the minimal preorder.  It does **not** earn equality to complete physical
chronology.

Per-mission DCL custody is already part of the word **qualified** in (FB03).
U-DCL universalizes the availability of that witness over admitted actual
records; it does not create a positive mission, choose the cross-mission event
identity, or prove causal completeness.

## 4. Three decisive non-implications

### 4.1 Per-record acyclicity does not give cross-mission acyclicity

Mission `mu_1` may have the finite DCL graph `A -> B` in context `z_1`, while
mission `mu_2` has `B -> A` in context `z_2`.  Each record separately has a
valid acyclic DCL witness.  If `A,B` are nevertheless glued as the same event
identities, their union cycles.  Often the correct repair is to distinguish
the context-indexed events `A_1,A_2,B_1,B_2`; U-DCL does not itself supply that
cross-mission identity decision.

Collapsing a strongly connected component would produce a block order, but it
would not prove an event-level chronology and could hide the failed identity
or contextuality test.

### 4.2 Positive influence is not transitively complete

One mission can witness `A leadsto_R B` and another `B leadsto_R C` while no
mission witnesses `A leadsto_R C`: the two instruments may be incompatible,
the second context may reset the `B` carrier, or the composed mission may not
be physically admitted.  The transitive closure adds `A preceq_R C`
mathematically, but does not turn it into a positive recorded mission.

Conversely, an underlying causal route can have zero registered influence
because the available coupling is constant, depolarizing, symmetry-forbidden,
or below every admitted margin.  Therefore

\[
 \neg(e\leadsto_R f)
 \not\Longrightarrow
 \text{physical spacelike separation}.             \tag{FB04}
\]

Equality to physical chronology requires a closed intervention/read domain
and causal faithfulness or front saturation, not merely the union operation.

### 4.3 A complete material influence front need not be the light cone

The same finite record-influence process with coordinate front speed one is
compatible with

\[
 g_A=-dt^2+dx^2,
 \qquad
 g_B=-4dt^2+dx^2.                                  \tag{FB05}
\]

Its front is null for `g_A` and strictly subluminal for `g_B`.  Every terminal
record law and arrow in (FB03) can be identical in the two realizations.  Thus
even a perfectly mapped record/material front does not identify the causal
cone unless the admitted probe domain is physically complete and the front is
proved to be the universal maximal front.  Two probes inside an unclosed
domain do not exclude an omitted faster sector.

These countermodels reproduce the exact boundary of the existing finite
operational-order theorem; terminal recordhood repairs arrow custody, not
global causal completeness.

## 5. What causal order plus volume actually yields

The Hawking--King--McCarthy/Malament result is a **uniqueness theorem on an
already admitted Lorentzian manifold**.  In the form needed here, assume:

1. one connected, time-oriented, smooth four-manifold `M` with a smooth
   Lorentzian metric exists (and therefore has the `C2` regularity needed
   downstream);
2. it is future- and past-distinguishing (strong causality is sufficient);
3. the strict cycle-free completion on the limiting event set `E_infty` has
   inner/outer convergence controlling the null boundary; and
4. event custody supplies a bijection `iota:E_infty -> M` that is a
   chronological isomorphism in the arrow orientation used here:

   \[
    e\prec_R f
    \quad\Longleftrightarrow\quad
    \iota(f)\in I^+(\iota(e)).
   \]

Then the complete chronology determines the topology, differentiable
structure, time orientation, and conformal Lorentzian class `[g]`, up to the
appropriate diffeomorphism.  It does not determine the conformal factor and it
does not prove that an arbitrary preorder has such a realization.

This clause uses strict chronology deliberately.  If an operational
construction instead retains exact null-endpoint arrows and proposes the full
causal relation `J+`, it must invoke the corresponding causal-isomorphism
theorem or independently recover `I+`; one may not silently identify a generic
positive-influence relation with chronology.

Now additionally supply a smooth positive measure `mu` independently
identified and absolutely calibrated as physical spacetime four-volume.  For
one representative `g_0 in [g]`, write

\[
 d\mu=f\,dV_{g_0},\qquad f>0.                       \tag{FB06}
\]

If `g=Omega^2 g_0`, then in four dimensions
`dV_g=Omega^4 dV_g0`; hence

\[
 \boxed{
 \Omega=f^{1/4},\qquad
 g=f^{1/2}g_0,qquad
 dV_g=d\mu.}                                       \tag{FB07}
\]

Thus complete chronology plus physical four-volume determines one Lorentzian
metric up to diffeomorphism.  If volume is calibrated only up to a constant,
one global metric scale remains.

The words **physical four-volume** are load-bearing.  An event count, record
count, gamma exponent, Hilbert dimension, spatial volume, or arbitrary positive
cell weight is not four-volume without a constitutive calibration.  Likewise,
an arbitrary causal set need not be manifold-like.  Existing M1--M4 already
prove the projective causal-measure and sharply conditional reconstruction
steps once these premises are supplied.

## 6. Minimal residual law for this route

The bypass replaces RLS with the following two-clause physical law.

> **RCV -- record causal-volume realization law.** On one same-parent domain:
> (A) nonmetric cross-mission event custody, a closed complete intervention/read
> census, support-faithful scale maps, and causal faithfulness make the strict
> part of the cycle-free completion of all qualified positive terminal-record
> influences converge to the complete distinguishing chronology of one smooth connected
> four-dimensional Lorentzian manifold, common to the maximal fronts of clocks,
> matter, EM, and independent probes; the compatible event measure is
> independently calibrated as its absolute physical four-volume.  (B) a
> prospective KEEP versus whole-lineage BREAK/reprepare intervention is run
> with clause A valid in **both** arms under one common event identification,
> volume calibration, and intact independent reconstruction-probe/read domain.
> With complete stress, work, controller, boundary, EM, probe, failure, and
> quarantine ports matched, the intervention changes the completed order,
> four-volume, or reconstructed metric through the authenticated record
> lineage.  Merely deleting a record, query, or reconstruction probe does not
> satisfy this clause.

Clause A contains two irreducible existence statements that Malament does not
provide: **manifold/cone realization** and **metric-volume identity**.  Clause B
is separately necessary because terminal records can merely detect a causal
geometry created by another sector.  Same parenthood excludes a late
cross-model join but not a common-cause or spectator explanation.

One passing clause-B intervention proves one constitutive effect of that
authenticated record lineage on the reconstructed geometry.  It does not by
itself prove that every geometric degree of freedom or all gravity originates
in records, nor does it supply a universal response law, curvature equation,
or Einstein dynamics.

Once clause A holds, (FB07) reconstructs the kinematic metric without a
Fisher/QFI query, categorical fidelity, or Regge edge assignment.  Its `C2`
metric has a mathematical Levi--Civita connection.  Physical transport beyond
the maximal cone still requires the common matter/EM/probe coupling and
torsion/nonmetricity/dispersion controls; Einstein dynamics still requires
SDCP/RIEHB or another independently proved response law.

## 7. Honest premise comparison

| Burden | Fisher/FERS route | Operational order/volume route |
|---|---|---|
| incompatible reads | must be dominated by one sufficient classical query | allowed as separate positive missions |
| local metric | query/state fidelity and QFI saturation; Chentsov naturality; scale/time solder | not used |
| causal cone | obtained only after physical localization/clock/probe soldering | requires all-probe maximal-front completeness and causal faithfulness |
| continuum | shape-regular Regge/common-query refinement and a specified smooth reconstruction | manifold-like order convergence, nonmetric event gluing, and support-faithful causal-measure refinement |
| absolute scale | prospective clock/information-to-length bind | independently calibrated physical four-volume |
| metric variation/dynamics | RW-R/SDCP plus RIEHB | still RW-R/SDCP plus RIEHB; order uniqueness does not supply variations or an action |
| record origin | matched lineage intervention RW-A | the same matched lineage intervention, clause RCV-B |

The order/volume route is a real bypass of the **one-query Fisher soldering
gate**.  It is not a lower-premise proof.  It exchanges a local statistical
sufficiency law for a global causal-completeness, manifold-likeness, and
four-volume law.  Its most difficult new assertion is that the union of all
recorded positive influences is complete enough to recover the universal
maximal chronology, including causal relations whose available material
channels may carry no detectable signal.

Relative to the current program, RCV-A can replace RLS plus the Fisher/Regge
kinematic edge-to-metric construction.  It cannot replace multiscale
convergence: it replaces that route's Regge refinement with RCV-A's own
support-faithful order/measure convergence, manifold realization, and
chronological-isomorphism premises.  It cannot replace the six-spatial/
constraint or dense-deformation closure, the positive induced Ricci
coefficient, complete stress, or record ancestry.  The two routes are
therefore complementary experimental strategies rather than one being a proof
of the other.

## 8. Exact disposition

`AURFT_TRL_AND_POSITIVE_COMPLETE_RECORD_MISSIONS_PROVE_SOUND_DIRECTED_OPERATIONAL_ARROWS_AND_THEIR_MINIMAL_PREORDER__CROSS_MISSION_UNION_IS_NOT_AUTOMATICALLY_ACYCLIC_TRANSITIVE_COMPLETE_OR_MAXIMAL__MALAMENT_HKMM_FIX_THE_CONFORMAL_CLASS_ONLY_AFTER_MANIFOLD_AND_COMPLETE_CHRONOLOGY_ARE_SUPPLIED__ABSOLUTE_PHYSICAL_FOUR_VOLUME_THEN_FIXES_THE_METRIC__FISHER_RLS_CAN_BE_BYPASSED_BUT_IS_REPLACED_BY_GLOBAL_RCV_MANIFOLD_CONE_VOLUME_AND_RECORD_ANCESTRY_LAW__EINSTEIN_DYNAMICS_REMAINS_DOWNSTREAM`

## Primary theorem anchors

- Hawking, King, and McCarthy, [“A new topology for curved space-time which
  incorporates the causal, differential, and conformal
  structures”](https://doi.org/10.1063/1.522874).
- Malament, [“The class of continuous timelike curves determines the topology
  of spacetime”](https://doi.org/10.1063/1.523436).
- Surya, [“The causal set approach to quantum
  gravity”](https://arxiv.org/abs/1903.11544), especially the explicit warning
  that not all causal sets are manifold-like and that order/number-to-continuum
  correspondence is an additional approximation problem.
