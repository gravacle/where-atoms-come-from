# THE UNIVERSAL RECORD MODEL — the program's executable representation

`model/project_model.py`

```python
from project_model import URM

m = URM()
s = m.surface(name, mechanism, dE, E_b, T, f0,
              provenance="<pinned source>")
```

The URM is where this program adds observations and laws. `model/record_model.py` remains its
first-principles corner engine: from `(H, {L_k})` it constructs every record admitted by the five
clauses. `model/project_model.py` carries that engine into the world tier with declared surfaces,
provenance refusals, layer methods, and validator gates.

Exactly one ledger row is `PROVED`: `C-71`. The URM does not upgrade a row's status; it makes the
row's computation, entry conditions, scope, and tests inspectable.

## THE THREE ENTRY DOORS

Every landed feature must use one of these paths:

1. **A new surface or venue** enters through a URM provenance gate. World-tier inputs require a
   pinned source; exact corners must self-declare `provenance="DEF-A"`. The consuming layer rechecks
   the declaration where bypass would otherwise be possible.
2. **A new law** enters as a layer method and a validator gate with a failing branch and a positive
   control. A sealed lane by itself is evidence, not a URM feature.
3. **An external number** enters as a pinned comparison beside computed output, with its units,
   tolerance, extraction semantics, and a power control stated. A source substring is provenance,
   not empirical confirmation.

Python metadata is enforceable at these program boundaries; it is not a cryptographic custody
mechanism. Sealed historical lanes may still contain raw constructors, but no new observation can
silently use them as the public path.

## THE LAYERS

| layer | public URM surface | what it carries |
|---|---|---|
| definition/laws | `clauses`, `lifetime`, `steady_value` | the five clauses, rates, and values |
| formation | `configuration`, `formation_occupancy`, `formation_orientation` | how a record is written and read |
| corner | `corner` | the DEF-A exact idealisation, with explicit self-declaration on `URM` |
| geometry/roles | geometry and role delegates | located record geometry and the three-role ledger |
| arrow | `arrow_threshold`, `arrow_ledger`, `arrow_invariance`, `arrow_history`, `arrow_redundancy`, `arrow_observation` | record-copy threshold, history, and fragments |
| count law | `census`, `count_widths` | the surviving-record staircase and both durability widths |
| classes | `coupling_venue`, `reachable_class`, `critical_kernel`, and related delegates | subcritical, critical, and supercritical coupling classes |
| writing | `writing_kernel_verdict`, `writing_uniformity`, `writing_transport`, `writing_trail_*`, `writing_gap` | conservation, criticality, transport, and trail diagnostics |
| physics-role status | `roles` | claim-scoped EM, gamma, Alpha-inheritance, and gravity status, including the exact first-record gamma-information seed, the SAI/AWAI actual-visible-sector theorem, the adopted off-shell/on-shell RGRL distinction, and the conditional record-front cone refinement; descriptive metadata, not a gravity solver |
| U-DCL working postulate | `udcl_postulate`, `udcl_postulate_certificate` | the checksum-pinned adoption decision, typed conditional theorem, and full transitive axiomatic closure `U-DCL -> universal Coverage-U`; zero caller inputs; natural validity, actualization, Born selection, A5, and gravity remain open |
| formal historywise-gravity discriminant | `historywise_gravity_discriminant`, `historywise_gravity_discriminant_certificate` | a zero-input, custody-pinned finite-group nonselection/stabilizer theorem certificate; no gravity solver, outcome selector, physical GARH-D/Q decision, Born law, GR derivation, or empirical proof |
| Gravity Formation Theory | `gravity_formation_theory`, `gravity_formation_theory_certificate` | a zero-input, hash-pinned certificate for adopted RGRL, its adopted off-shell/on-shell clarification, V002 dressed-response ceilings, the conditional record-front cone theorem, the conditional pair-resolved Coulomb `DPAR` realization, FV projected rank six, fully sealed Phase A native-support/family response, the FZ/GA/GB Ward boundaries, GC finite-family geometry, sealed GD total-momentum ownership, the independently sealed GF V005 G2 observable contract, audited public searches, and the exact conditional working-theory closure; the all-`G_L` ledger, matched same-ancestry massless helicity-two pole, soft universality, Einstein self-coupling, empirical RGRL/EIR/AFR confirmation, and microscopic numerical-(G) derivation remain open |

The field-instrument family is not listed: T-51 is still independently unverified and nothing from
that lane is registered or folded into the URM.

## THE CORNER ENGINE

```python
from record_model import RecordModel

r = RecordModel(H, Ls)        # Hamiltonian and Lindblad operators are the entire corner input
records = r.records()         # every record the pair admits
family = r.independence(records)
```

No lattice, gauge group, temperature, coupling constant, code, or geometry is invented on this
path. If a value is not derivable from `(H, {L_k})`, the corner engine does not have it.

| step | what it computes | theorem |
|---|---|---|
| `star_algebra` | `A = alg{I,H,L_k,L_k†}` | C-9 |
| `commutant` | `A'` | C-9 |
| `minimal_projections` | a maximal splitting allowed by `A'` | C-10 |
| `clause_iii` | non-triviality on an eigenspace | anchor (iii) |
| `clause_iv` | trace balance on every eigenspace | C-11 / O-4 |
| `build_writer` | an admissible `U`, with `[U,H]=0` and `U†RU=-R` | C-11 |
| `commuting_family` | independent record bits | C-14 |

The exact corner count is `k = min_E v2(m_E)`: every independent record must halve every
eigenspace. The naive `floor(log2 min m_E)` control fails on `[3,3]`, `[6,6]`, and `[5,5]`.
Clause (v) remains carrier data; the model raises rather than inventing locality.

## FORMATION AND OBSERVATION

`RecordModel.evolve` shares one eigendecomposition across many readouts. `formation`, `redundancy`,
and the arrow layer then score what the environment and its fragments hold. A coupling may be a
product operator, a distributed list of system-term/bath-site pairs, or a full interaction operator;
those are physically different inputs and are not silently interchanged.

`channel()` uses the corrected G-16 criterion: the coupling's compression onto the code space must
have a non-zero component along the record. Anticommuting with the writer is necessary on the gated
venues and is not sufficient in general.

## U-DCL WORKING POSTULATE

The two zero-argument U-DCL methods verify the adopted decision, the hostile-audited decision
basis, the sealed standard-causal dependency, the typed conditional-theorem lane, and the full
axiomatic closure lane. They freshly reproduce both the 72/72 typed checks and the 74/74 transitive
proof-custody checks before returning an immutable certificate. The certificate distinguishes the
two authorized program statements—U-DCL is the adopted working postulate and universal Coverage-U
is proved conditional on it—from scientific claims that remain open. In particular, it assigns no
empirical weight to caller input, does not infer `REC -> DCL_phys`, and does not authorize actualization,
outcome forcing, a Born law, A5, gravity, GR, `G`, `Lambda`, U(1), or `alpha`.

The landed predicate is the typed V002 implementation: `DCL_phys` is the ontic
existence of one common K/W witness across frontier, state, history, maps,
instruments, locality, and provenance. `Cert_DCL` is the separate prospective
evidential packet. Failure to obtain a certificate is not evidence that an old
record lacks the physical structure.

One independently admitted actual bona-fide finite-mission record that fails every applicable
prospectively frozen DCL category (or the independently warranted exclusive category) falsifies
the postulate. A finite collection of successful episodes cannot prove the universal physical
antecedent.

For one exact admitted reference mission implemented as a finite acyclic causal
network, the separate composable-robustness theorem transports the RMR margin
through independently certified full-channel errors:
`D_TV(actual_1,actual_0) >= [B_RMR-E_1-E_0]_+`. Query-nonancestor nodes cancel
and each shared ancestor is charged once. This requires common complete port
types and does not create approximate `DCL_phys`, authenticate hidden ports, or
infer `REC` from contrast.

Finite recurrent apparatuses are now included by a derived rather than assumed
unrolling whenever the mission has a finite horizon and bounded scheduling
fanout and every feedback cycle has strictly positive total lower delay. The
resulting occurrence DAG gives `FDFU -> FCLPD -> DCL_phys`, with Coverage-U after
independent `REC`. Robustness charges every physical use once; reused hardware
may be charged once only under a complete joint multi-slot strategy certificate.
Instantaneous algebraic loops, Zeno accumulation, unbounded spawning/horizons,
hidden cross-use memory, and query-before-later-write packets remain outside
this subclass theorem.

## PHYSICS ROLE STATUS

`URM.roles()` now carries the proof-safe interface among electromagnetism, record
distinguishability gamma, alpha, and gravity. It is descriptive status metadata,
not a numerical gravity solver or a new theorem certificate.

The long gravity development ledger below is retained as historical custody in
`ProjectModel`. Public `URM.roles()` replaces its `GRAVITY` value with the concise
current closure status and directs exact custody/ceilings to
`gravity_formation_theory_certificate()`.

- EM is an instantiated mechanism in declared record packets. Compact \(U(1)\)
  and Maxwell response still require their independent ancestry/action/Ward
  premises.
- For the first qualified positive-margin binary record in a lineage,
  complete-query squared-fidelity \(\gamma_Q\) obeys
  \(D_{\rm TV}\ge B_{\rm rec}>0\Rightarrow
  \gamma_Q\le1-B_{\rm rec}^2<1\), so
  \(I_\gamma=-\log\gamma_Q>0\): the gamma-information seed is already
  present at record formation. This is not microscopic classical gravity.
  State-family \(\gamma_{\rm state}\) supplies an exact rank-three QFI
  coframe candidate. At the isotropic `S_4` point, GSGB's QFI convention
  `s_gamma=ell_F^2 f P` matches EO's `q_EO=(4a^2/3)P` under
  `ell_F^2 f=4a^2/3`. Equivalently, ET's Bures-line convention
  `q_gamma=(ell_B^2 f/4)P` requires `ell_B^2 f=16a^2/3`, with
  `ell_B=2 ell_F`. The symbols are kept distinct; away from isotropy an
  explicit tetrad/null binding is required. The restrictive same-parent influence-functional
  \(\gamma_{\rm IF}\)/KMS packet reconstructs its declared dressed retarded
  response. These gamma objects are not identified without a separate common
  read/dilation theorem; generic gamma is not a force and supplies neither time
  orientation nor curvature. With the separately supplied ER routing law, the
  same retained record coordinate can control exact common two-probe Lorentz
  holonomy. Because that routing is supplied, this is an intermediate witness,
  not the final same-world theorem. Its global same-sector ancestry and physical
  soldering remain open.
- Alpha is not a standalone theory in the URM. `SAI/AWAI` is the empirically
  anchored same-visible-\(U(1)\) inheritance theorem. Measurement anchors the
  actual visible-parent coupling at a reference context \(\chi_0\). For every
  record independently established to satisfy `ACTVIS` and `SAI1--SAI8`,
  \[
    \operatorname{ACTVIS}(r,W_{\rm obs})
    \Longrightarrow
    \alpha_r(\chi)=
    {\cal T}_{\chi\leftarrow\chi_0}[\alpha_{\rm obs}(\chi_0)].
  \]
  Therefore an `ACTVIS` record cannot possess a private record- or region-level
  alpha; an aligned inequivalent value falsifies its same-sector assignment or
  at least one SAI premise. Bare `REC`, `DCL_phys`, and URFT establish neither
  `ACTVIS` nor compact \(U(1)\). Finite active-EM recordhood `ALLOW`s multiple
  alphas across parent models; complete-universe numerical `REQUIRE` and any
  parent `SELECT` law remain open. AWAI proves inheritance of the empirically
  anchored value, not a parameter-free prediction of that parent value or a
  derivation of gravity or \(G\). Exact custody and the 41/41 witness are in
  `LANE_RFT_ALPHA_SECTOR_INHERITANCE_V001`.
- Classical gravity is not placed on the microscopic record surface; what is
  already present in its first qualified binary record is the nonzero
  gamma-information seed above. The active successor
  calculations now derive a record-loading collective boundary
  \(4U_d\chi_0^R(0)=1\), exact open classical degree plateaus, and—conditional
  on physically retained pair-gated serial support—an exact open one-dimensional
  TFIM phase for \(U_d>2h\). An exact finite custody-handoff model now forms the
  required derivative pair memory, composes it into arbitrary finite serial
  paths and one symmetric equal-depth branch/rejoin cycle, and implements an
  extensive reversible active/quarantine BREAK. On independently earned
  bipartite support, the degree plateau is exactly a discrete Gauss sector and
  the existing link flip generates alternating-cycle ring dynamics. On
  supplied `z=4,d_*=2` diamond ice, the exact leading F3 Hamiltonian reaches the
  published pure-kinetic \(U(1)\) liquid. Independently, the actual hard-core
  carrier makes every fixed even cubic torus a local two-switch free-energy
  basin in an open controlled regime; the same saturated force makes periodic
  diamond a saddle, so those slice results cannot be added. On supplied cubic
  eligibility, `d_*=3` is exactly the zero-charge spin-half U(1) quantum-link
  sector. Formed pair memory can lawfully gate that distinct successor field,
  but unchanged BS09 hopping follows its occupied `G_n`, not the full
  degree-six eligibility graph; the saturated carrier margin is therefore not
  additive and no direct `K_eT_e` law is installed. A joint trace has a nonzero
  carrier--incidence backresponse, and its complete `h^2Y^2` dressing is known:
  it produces carrier density/hopping response but not symmetric-detuning
  electric stiffness. The exact full thermal square has open positive and
  negative response regions and multiple tuned crossings, so `beta U_d=1` is
  not a universal threshold. Exact physical shared-edge composition preserves
  both one-cell anchor signs; its sector-correct global gluing contrast and a
  distinct same-six-site pinned-ensemble Möbius finite difference are nonzero
  and sign-changing. This remains a two-cell, one-carrier dilute result, and
  the pinning controller/work is not yet physically completed. The unchanged
  qutrit law already supplies the fixed-content hard-core two-carrier sector,
  with `(J_e^psi)^2=q_u+q_v-2q_uq_v`. Its exact 1920-state trace proves that
  finite signs depend on collision algebra and number preparation; the allowed
  global rule `N_m=m` is not a selected fixed-density phase or a full-content
  trace. The exact fixed-width transfer limit now proves that the matched
  canonical and grand intensive responses agree but remain analytic at finite
  temperature: a strip zero is a crossover, not a phase transition. On a
  supplied two-dimensional square sheet, by contrast, the unchanged incidence
  trace is exactly a record-occupation-dependent Ising interaction and one
  witness straddles its true thermodynamic critical surface. This is a
  macroscopic statistical mediation theorem. A second exact theorem restores
  positive degree lock and proves open blank-boundary-disordered/occupied-
  carrier-ordered regions at a BS13-stable witness. A typed, volume-wise
  composition now routes actual FPMH derivative records into active storage or
  retained quarantine and thereby switches those carrier responses while
  conserving total formed-relation occupation across KEEP/BREAK. This is exact
  conditional formation-to-phase causation through active occupation, not a
  history-only effect at matched active state. Pairing now has two sharper
  results. An intermediate fixed-separate-marginal construction is rejected as
  record causation because its end-to-end decoder cancels record ancestry. A
  distinct pair-rebinding circuit changes the phase by changing active joint
  correlation while complete separate fields, sharp counts, spectra, and mutual
  information remain matched. Heterogeneous accumulation is controlled by
  contour/path topology rather than scalar record density: equal-density
  placements have different exact responses. Composing FPMH with that theorem
  holds fixed 32 formation episodes, eight active authenticated records, and 24
  quarantined formed excitations while retaining the response difference. F3
  reads the active physical pattern, not provenance; the topology is supplied.
  The reciprocal audit also proves that the sealed `h=t=0` phase hold has zero
  retarded carrier-to-incidence and carrier-to-pair-memory response. Restoring
  the inherited BS06 link flip yields an exact finite carrier-conditioned Rabi
  channel, but no common-background collective limit has been earned and its
  fixed-time local response is `O(N^-2)`. The minimal incidence degree-
  defect then acquires second-order line-graph hopping but fails the isolated
  finite-residue pole gate because a flat cycle fiber touches the dispersive
  band and the degree residue vanishes; inherited flip scaling also gives
  vanishing thermodynamic bandwidth. These results do not yet supply autonomous
  support, finite-temperature phase persistence at nonzero dynamics, a common
  metric, or gravity. A separate lawful joint endpoint trace has an
  exact open conditional ALLOW-side bias without `K_eT_e`, but physical sector
  access and matched occupied-`n` switching remain open. Most directly, qualified retained occupation changes the
  virtual degree-defect denominators, producing an exact record-patterned
  electric potential and plaquette stiffness without inserting a geometry reward. If one phase carries
  the resulting electric and magnetic stiffnesses to the same infrared scale,
  their product and ratio determine an optical cone and coupling; speed Hessians
  can carry sectoral tidal curvature. Maxwell conformal invariance nevertheless
  leaves the spatial-volume factor needed for nonzero static \(G_{00}\) unresolved
  in the U(1) sector alone. A calibrated same-parent massive probe would fix
  `N=omega_0/mu` and `A=omega_0/(mu c)` and permit nonzero static `G_00`, but
  without an earned proper mass/clock/residue an exact compensator leaves that
  volume unidentifiable for one uncontrolled probe. Several fixed-mass probes
  provide an exact scale-free alternative: a common static metric exists exactly
  when their cone identities agree and all coordinate-gap ratios are spatially
  constant. Those conditions determine the `N,A` profiles and coordinate `G_00`
  up to one global scale; a varying gap ratio falsifies common fixed-mass
  geometry. Current F3 supplies only a prospective dispersion test and has not
  yet earned the required independently resolved continuum species. Its two
  bare content labels give duplicate rank-one carrier rows, and its storage
  qutrits do not propagate between vertices. The minimal bare incidence degree-
  defect candidate has failed the stable-pole gate at its first controlled
  hopping order. Exact screening of its already-present carrier-dressed one-hole
  continuation now permits static dressing or binding but proves an isolation--
  bandwidth dichotomy: exact-complement touching cannot create a uniform gap,
  while a uniformly isolated continuation has bandwidth `O(N^-1)` under
  inherited `h_N=Omega/N`. Finite-density, controlled gap-closing/double-scaling,
  and higher-sector routes remain open. The exact finite-density boundary does
  contain an `O(1)` channel, but only with macroscopic Johnson coherence:
  authenticated classical word accumulation and the conditioned dressed Gibbs
  state have `c_J=0` and `O(N^-1)` escape curvature. Separately, the earned
  `Q_uv^auth` is a partial set-valued authenticated custody field: declared
  support is symmetric, writer provenance directed, and untested pairs are
  undefined rather than zero. Explicit common-parent products realize every
  supplied finite graph, so URFT/FPMH alone entails no sparsity, dimension, or
  locality. Its exact common-generator `kappa_(v<-u)` seminorm distinguishes a
  direct generator block from two-step-only influence. The unchanged one-pair
  FPMH/PESC composition now passes it exactly:
  `kappa_KEEP(v<-u)=kappa_KEEP(u<-v)=abs(t_psi)/hbar`, while both matched BREAK
  blocks vanish under the same complete Hamiltonian. On one complete finite
  declared census, this composition now gives
  `kappa_(v<-u)=abs(t_psi) S_uv^auth/hbar` in both directions on every pair, so
  direct influence support is exactly active authenticated support. The fixed-
  number finite-time signal equals the carrier transition amplitude and obeys
  separate factorial and Lieb--Robinson-type graph-distance bounds conditional
  on a supplied maximum degree. This earns a supplied-network relational
  propagation cone without `K_eT_e`; endpoint reciprocity is not carrier-to-
  record back-reaction. In the complete existing post-formation source-off
  parent, the full authenticated support projector
  `S_e=F_e Z_e^KEEP P_e^K` commutes with the Hamiltonian. Every support-word
  probability is therefore frozen for arbitrary complete states. This proves
  passive graph retention but closes autonomous or corrective selection of the
  authenticated graph by existing source-off terms; successor incidence may
  still evolve inside fixed authenticated eligibility. The exact `K_(6,6)`,
  degree-three hostile screen proves finite component-merger/fragment access,
  but the completed known one-carrier diagonal response gives the fragmented
  `2K_(3,3)` graph strictly lower fixed-graph free energy at finite temperature.
  Physical zero-temperature band bottoms tie, and graph entropy is separate.
  The exact carrier-dressed `h^4Y` switch coefficient cancels at symmetric
  detuning; in the controlled detuned domain its correction is below `3/16`
  per labeled switch and `1/48` on the normalized merger block, so it cannot
  reverse the pure switch. The reciprocal Hermitian amplitude is not merger
  tension. The exact all-orders-in-carrier diagonal response at `O(h^2)` is a
  Sylvester/Liouville inverse. In the finite `K_(6,6)` witness its relative
  fragment/connected sign changes with temperature, while proper cold
  degenerate perturbation returns to fragment preference. It is real feedback,
  not a universal connector. The first scale-level ensemble theorem now counts
  all fixed-`d` simple bipartite successors inside dense `K_(n,n)` eligibility.
  Their rank is `exp[d n log n+O(n)]`, whereas a prospectively fixed finite-
  complexity local family has rank at most `exp[2n log n+O(n)]`. Consequently,
  for `d>=3`, fixed finite temperature, equal support fibers, and the explicit
  endpoint-extensive premise `osc(H_n)<=w n`, the local-family Gibbs weight
  vanishes as `exp[-(d-2)n log n+O(n)]` even for off-diagonal graph dynamics.
  This is conditional rather than an inherited full-parent verdict: dense
  eligibility contains `n^2` physical relation resources, fixed-amplitude
  switch kinetics has `Theta(n^2)` width, and a dense endpoint projection may
  be sparse in an enlarged record/factor/event graph. The sealed CP amplitude
  satisfies the bound for its displayed leading switch block, but the complete
  relative-parent width is not yet proved. Counting each active pair token as
  a physical factor vertex gives exact average degree `4M/(N+M)<4`, but this
  does not cure endpoint hubs: the complete endpoint projection has maximum
  degree `N-1` and radius-three saturation. In the present DU/DV carrier law
  the token prepares `n_uv` but is not a propagating intermediate, so the
  endpoint projection remains the exact direct-influence graph. Existing URFT,
  FPMH, and source-off conservation imply no system-size-uniform endpoint
  capacity. The isolated missing premise `UCAIC_sigma(B)` requires every
  simultaneously active relation to own one of at most `B` complete named
  ports at each endpoint until BREAK/deactivation, with no direct bypass.
  Conditional on it, `Delta(K)<=B`, `M<=BN/2`, every fixed-size cell type has
  only `O(N)` embeddings, and DV's cone constant is size independent. UCAIC is
  not yet derived or adopted and selects neither connectedness nor dimension.
  The strongest finite-dimensional record-capacity route is now closed
  exactly. If all `2^m` independent incident episode words live in one common
  `D`-dimensional endpoint-accessible carrier, perfect whole-word decoding
  requires `2^m<=D`; approximate whole-word decoding obeys Fano--Holevo, and
  separately readable jointly randomized bits obey
  `sum_i[1-h_2(p_i)]<=log_2 D`. A jointly stable averaged RMR margin `b_i`
  may replace `p_i` through `p_i<=(1-b_i)/2`, with an analogous weaker gamma
  bound. But ordinary per-context RMR/gamma does not supply that joint code:
  parity, central-fanout, external-edge-memory, and shared-bus countermodels
  retain perfect nominated contrasts while defeating uniform independent
  endpoint ownership. Thus information capacity becomes active degree only
  after episode-faithful ownership, persistence, a system-size-uniform
  accessible dimension, and no bypass are established separately.
  History-wise inheritance is now exact once its physical premise is stated.
  Let `K_j` be the complete typed phase state, `A_j` its authenticated support
  word, and `G_j` the simple active endpoint projection. If a prospectively
  fixed realized seed lies in `Phi_*` and one complete formation instrument
  maps every nonterminal active, reject, and quarantine outcome back into
  `Phi_*` (terminal outcomes end their histories), finite induction gives
  `K_j in Phi_*` on every realized history without selecting an outcome.
  Conditional UCAIC bounds `G_j`; DW freezes `A_j` and hence `G_j` between
  writer windows; and DV supplies its bounded cone on the controlled simple
  one-episode preparation. Existing snapshot theorems do not derive the update
  law: bounded supports can alternate, a monotone path can grow into a binary
  tree, and repeated four-cells can form a cactus. The isolated candidate is a
  recursive same-parent writer law in which new active writers descend through
  realized phase-local carriers/fronts or typed boundaries, every active
  grammar branch preserves `Phi_*`, and bypass requests are rejected or
  quarantined. That writer law and the seed phase remain physical premises,
  not consequences of current FPMH.
  The complete dense-parent width ledger now separates the exact locked
  compression from the full soft parent. Under explicit bounded edge, bulk,
  and port premises the degree-`d` compressed relative block has `O(n)`
  spectral width; its one-toggle boundary is intensive for `h_n=Omega/n`, the
  exact-`Y` second-order self-energy is `O(1)` under a uniform Sylvester gap,
  and the displayed leading switch width is `O(n^-2)`. None of these is an
  all-orders effective-width theorem. The full soft degree penalty instead has
  `Theta(n^3)` width. More decisively, on the clean classical symmetric slice,
  `Z_1/Z_d=n^2 exp(-2 beta U_d)` and every family with fixed finite maximum
  degree has Gibbs weight at most `exp[-n log n+O(n)]` at fixed finite
  temperature and fixed couplings. Thus a fixed soft degree preference inside
  dense eligibility does not thermally form finite valence on that slice.
  Hard active capacity, formation-time support-rank reduction, a growing scale,
  zero-temperature control, or nonequilibrium protection remain open escapes;
  this thermal conclusion is not exported to the nonzero carrier/current/port
  parent. The next locality theorem separates bounded cellular formation from
  genuinely subextensive collective growth. Conditional on UCAIC, a connected
  attaching lineage of bounded authenticated cells with bounded edge overlap
  yields connected endpoint and cell-dual graphs and a linear cell census
  `(N-1)/e_*<=m<=qBN/2`. This is still insufficient: replacing every edge of
  a finite trivalent tree by an authenticated four-edge diamond satisfies
  `UCAIC(6)` and the complete attaching-cell conditions while its even-radius
  balls grow exactly as `9*2^k-8`. The isolated sufficient premise is the
  physically counted subextensive authenticated frontier law
  `|partial_E B(v,r)|<=A|B(v,r)|^(1-1/D)` on every intrinsic ball. Conditional
  on that law, endpoint and cell-dual balls grow at most polynomially,
  `|B(v,r)|<=(1+Ar/D)^D`, and no uniform expander family survives. `A` and `D`
  are supplied frontier data, not a derived dimension. The live physics target
  is to derive or falsify subextensive front collision, coalescence, or sealing
  from the formation, work, port, and retention ledger.
  One exact conditional mechanism now reaches that frontier law. If two
  complete reusable record-bearing operations physically coalesce every
  future-equivalent `wAB/wBA` event front, and the resulting classes,
  transitions, and squares are separately lifted to authenticated endpoints,
  pair episodes, cells, and named bounded ports, the support is `N_0^2`. It
  obeys the uniform all-root bound
  `|partial_E B(v,r)|<=8 sqrt(2)|B(v,r)|^(1/2)`; finite count-capped missions
  obey a size-independent version and therefore compose with the cellular
  frontier theorem. Constitutive merger BREAK restores the binary tree and
  exponential growth. This is a conditional physical-role theorem, not a
  current-parent derivation: FPMH does not fuse complete fronts, present
  comparison edges do not quotient histories, current F3 has no merger
  observable, and the tuned cycle zero is not a merger-balance law.
  Keeping the append arrows yields a second exact but narrower statement. With
  `t=#A+#B` and `x=#A-#B`, operational-front reachability is the parity-lattice
  condition `Delta t>=|Delta x|`, and the two elementary steps are null for a
  flat `1+1` representative. Order-interval counting gives the exact
  combinatorial identity
  `tau_comb^2=4(|I|-ell-1)`. The object proved here is a truncated operational
  transition poset, not yet a physical spacetime-event census. A
  one-class/one-event lift, uniform physical event-volume density, the full
  conformal factor, clock/length calibration, and common-probe use remain
  open. Naively adding commuting operation types does not solve `3+1`: for
  rank at least three the raw positive orthant is polyhedral, whereas a
  Lorentz cone has continuously many null rays.
  The stochastic OPEN-front ledger further sharpens what physical regulation
  would have to do. For the exact increment `xi=G-C-2P-Q`, drift sign controls
  survival/sealing and zero drift gives only a square-root high-probability
  envelope for one prospectively fixed exploration. It does not give a
  history-wise or simultaneous all-ball frontier law. The non-tautological
  target is a branchwise seal-credit/lifetime invariant plus an explicit map
  from named OPEN tokens to intrinsic boundary edges. Separately, a phase that
  is genuinely tensor-local to a redundant full-word history factor--or is
  restricted to one fixed front block--cannot create coherent front
  propagation: partial-trace invariance and CN descent make it front-invisible.
  This theorem does not cover a `Q`-controlled phase on CN's minimal rank
  encoding, whose rank labels are reused across front sectors. A canonical
  checkerboard-like construction would use a front-active coin and conditioned
  shift before persistent which-path recording; more generally, a successful
  path recombination must become front-active or use a fully owned
  history-to-front descent break. Current F3/PESC contains no earned operator
  implementing that join.
  The exact positive complement is now also known. A complete history can be
  retained reversibly as an OPEN active disposition `(q,c)` or explicit
  terminal label together with a blind residual rank `z`. Right congruence is
  the recursive append gate, reference-stable `Z` descent prevents residual
  leakage, and failure of the further CPTP `C` descent is exactly the condition
  for bounded memory to affect the scored front. If `k=|C|` is uniform, active
  RFCD layers obey `F_n<=sum_i min(k,binom(n,i))<=2+k(n-1)`; arbitrary refined
  graph balls grow by at most the finite factor `k`, with degree at most
  `k Delta+k-1`. Last-symbol chirality gives exactly `F_n=2n` and cumulative
  count `1+R(R+1)`, whereas FULL-RESIDUAL BREAK restores `2^n`. This result
  does not restore interference between orthogonally tagged histories. No
  displayed current-parent variable supplies the required active statistic or
  coin: transported content is bounded but propagation-blind, FPMH direction
  is provenance while active `K` is undirected, and `K,n` are edge-indexed.
  A prospective extension now closes the corresponding finite propagation
  calculation on one supplied fixed authenticated support. Reusing the existing
  qutrit occupied contents as probe chirality, but newly instantiating an onsite
  coin and content-conditioned transfer, gives an exact one-carrier conditioned
  shift. A single complete route/SWAP unitary gives chirality KEEP versus a
  Q-only C-descending BREAK with the displaced content retained in garbage.
  Structural formation history is fixed and blind; the interfering alternatives
  are unrecorded probe paths, so persistent orthogonal path records recover the
  incoherent EJ limit. In the padded bulk,
  `cos Omega(k)=cos(ka)cos(theta)` and `|v_g|<=a/Delta tau`; the controlled
  scaling `a/Delta tau=c`, `theta=mc^2 Delta tau/hbar` has principal limit
  `H_eff=hbar c k sigma_z+mc^2 sigma_x`. This proves a prospective conditional
  `1+1` relativistic probe architecture, not a current-parent matter law. The
  coin, conditioned transfer, front embedding, coherent-before-terminal-REC
  schedule, physical scale, `3+1` support, common metric, and gravity remain
  unearned.
  A finite coin cannot repair that missing support. For any differentiable
  finite walk symbol, the principal differential rank is at most the number
  of independent support characters, and internal noncommutativity does not
  alter the raw positive-append orthant. On a separately supplied signed
  rank-three support, an exact Pauli split-step walk instead has an isotropic
  Weyl infrared principal cone and an explicitly anisotropic ultraviolet
  dispersion. Four channels with Gram entries `1,-1/3` give a coordinate-free
  local tetrahedral frame; exact scalar-square isotropy is equivalent to the
  Clifford anticommutators, and an exact nonzero Dirac mass requires internal
  dimension at least four. Two independently prepared probes passing the same
  scoped Clifford test establish a common-cone candidate only for those
  probes. Signed support, inverse channels, global gluing/holonomy, the coin,
  probe species, scale, universality, curvature, and gravity remain unearned.
  The corresponding conformal-scale bridge is now exact at conditional finite
  `1+1` scope. If a branch-stable proper angular frequency `mu_*`, an
  independently calibrated common null speed `c_*`, the EL-step/EI-layer clock
  lock, and an independently falsifiable one-base-`Q`-cell/one-physical-cell
  binding are physically supplied, then `Delta tau=theta/mu_*`,
  `a=c_*Delta tau`, `ds^2=-c_*^2dT^2+dX^2=-4a^2du dv`, and
  `v_square=2a^2`. Complete intervals obey
  `V_2=(a^2/2)tau_comb^2=2a^2(|I|-ell-1)`, while exact M3 refinement is the
  phase-square law `theta_m(A)^2=sum_(e in q^-1(A))theta_n(e)^2`. Without the
  independent frequency standard, an exact Weyl compensator preserves the
  RFCD order, EL phase, dimensionless momentum, and speed ratio while scaling
  metric and volume. Thus gamma, counts, and order do not secretly fix the
  physical scale. The calibration/cell binding, physical `3+1` support, continuum,
  curvature, common stress response, and gravity remain unearned.
  The sealed bounded q=4 stream witness now supplies an explicit compatible
  record-front merger with exact complete-port `S_4` covariance. Conditional on
  the declared clock/cell calibration, its count contrasts give the exact
  rank-three `A_3` tetrahedral frame and a local `3+1` Lorentz-signature simplex.
  At the isotropic point, matching this EO simplex to the QFI coframe uses the
  explicit GSGB lock `ell_F^2 f=4a^2/3`, or equivalently the ET Bures-line
  lock `ell_B^2 f=16a^2/3` with `ell_B=2 ell_F`; it is not derived, and
  anisotropic matching requires an explicit tetrad/null binding.
  The Boolean-cell theorem defines a common two-probe Lorentz connection, and
  the ER witness proves that one qualified retained record can control
  nontrivial versus identity holonomy with identical signed module inventory.
  This is common-connection curvature, not yet Levi--Civita curvature: the ER
  witness fails the natural-coframe torsion equation. GSGB therefore establishes
  a typed gamma-to-gravity connection lane, not a completed physical chain. The
  final real-world target is one same-world derivation from qualified record to
  gamma seed, physical QFI/EO soldering, selected Levi--Civita curvature under
  refinement, and universal stress/Einstein--Hilbert response, with neither
  geometry nor routing inserted as a premise.
  `GSGB-JOIN`, QFI-to-physical-coframe `G-SOLDER`, compatible torsion-free
  connection selection, controlled refinement, and the downstream gravity
  gates remain open.
  Nonlinear off-diagonal `h^4Y^2`/`h^6Y` remains open, but the active need is
  now a physical capacity/formation law rather than another coefficient
  ladder.
  Any candidate must reverse or bypass this adverse component tension without
  inserting connectivity, then
  derive stable sparse connected support, degree/dimension/distance, an isolated
  finite-residue local pole, and physical thermodynamic propagation.
  Once one common smooth Lorentzian metric is earned, a positive matched total
  Ricci coefficient, including the finite-shell contribution, conditionally
  supplies the full nonlinear leading-derivative Einstein--Hilbert action and
  complete-stress back-reaction; a composite metric also needs the explicit-force
  and dense-tangent-span gates. Autonomous support/phase selection, physical
  multi-species/common-volume realization, the common metric, symmetric rank-two
  order, global same-sector ancestry, physical soldering/refinement, coefficient
  sign, the absolute scale, and observed SI \(G\) remain the gravity-identity
  frontier.

## FORMAL HISTORYWISE-GRAVITY DISCRIMINANT

The two zero-argument historywise-gravity methods reproduce the sealed 64/64 finite witness and
return an immutable certificate bound to the closed lane manifest and independent ACCEPT audit.
The analytic result is narrow: endogenous equivariant gravity/history feedback cannot choose one
member of a fixed-point-free outcome orbit, while a transitive added-input orbit has the
mathematical capacity for a covariant deterministic selector exactly under the sealed stabilizer
criterion. The executable witnesses do not prove the general finite-group theorem, and neither the
negative result nor the positive algebraic boundary table supplies a physical selector. Every
physical, empirical, GARH-D/Q, actualization, Born-law, record-causes-gravity, and GR authorization
in the certificate is fixed false.

## GRAVITY FORMATION THEORY CERTIFICATE

The two zero-argument Gravity Formation Theory methods verify the exact closure
and audit; six pinned core theorem/audit pairs; thirty-four pinned no-laboratory
advance pairs; and the adopted clarification chain consisting of its final
adoption record and seal, the clarification source, and V002. The advance
custody now includes the q4 front/support/response sequence; the finite TT
composite-cumulant screen; the projected-ice constraint-origin screen; the
fixed-parent collective-metric origin screen; the source-rank/degree-pair
geometric-strain sequence and complete homogeneous H6 response; two public-data passes; the
finite-apparatus and NIST/BIPM \(G\)-readiness results; and the HUST-2018
processed dual-method forward, history diagnostic, nominal-kernel
reconstruction, conditional homogeneous-\(G\) quotient, and five-clock
common-potential screen. The methods return a recursively immutable
status certificate and accept neither data nor solver parameters.

The certificate now makes the adopted distinction explicit. RGRL-C supplies
full-rank or declared dense-range **off-shell constitutive ancestry** on the
local spatial-metric tangent; it does not by itself supply an on-shell force
law. The separately established or measured on-shell kernel \(H^R\) may vanish.
Under the V002 well-posed retarded quotient with no unresolved zero mode, equal
dressed source, equal remainder, and equal physical solution data give exactly
zero response. `GI21`, the physical
compatibility/type join between those derivatives, remains open and is not
supplied by RGRL-C.

The frozen SPAG V001 artifact remains historical, and its local-RGRL-C
force/common-freefall verdict labels are retired.  The no-laboratory public-data
substitute has now executed the corrected complete-source-matched Lane A.  None
of the admitted Page--Geilker, Fuchs, NIST/BIPM, or Panda packets contains the
randomized same-parent eight-cell lineage intervention; even the deliberately
generous Page--Geilker proxy has support and rank `2/8`.  Those packets therefore
do not identify the lineage interaction \(\beta_{TM}\).  The NIST/BIPM summary
uncertainties give only an optimistic planning envelope of
`0.000192869--0.000578607 nN m`, not an achieved SPAG limit.  This is an exact
data-design ceiling, not evidence that a lineage response is zero.

A bounded second pass froze 28 public searches and retained no new qualifying
same-parent lineage root. It did retain a blinded five-ensemble strontium-clock
network and a levitated femtonewton gravity-drive workbook as protocol
components. Neither contains the randomized lineage cells, complete run-level
covariance, and source ownership needed for \(\beta_{TM}\); the Panda response
holdout remained closed. The search is reproducible at the query/custody level
but non-exhaustive because returned-hit rankings were not frozen.

The independently source-calibrated \(G\) lane now has the requested finite-
apparatus forward model: exact declared extended-source Newtonian torque and
stiffness kernels, a calibrated two-mode dressed Schur response, single
derivative ownership, and the exact same-source/remainder/data zero.  A
nonsingular calibrated row identifies the product \(p=Gs\), while a free global
source scale obeys \(F(G,s)=F(Gq,s/q)\); an independent
\(s\in[s_-,s_+]\) calibration converts a product interval to
\(G\in[p_-/s_+,p_+/s_-]\).  The `15/15` execution is synthetic validation of
that calculation, not a measurement of \(G\), a record-derived numerical
constant, or a lineage source law.

The public NIST/BIPM apparatus has also been reduced to
\(\Delta N_j=GA_j+r_j\), with
\(A_j=16\Gamma_jm_sm_t/R_s\), for four configurations in free and servo modes.
The eight summary observations have a nonzero source column and nominal zero
source stiffness at the torque extrema. The `G`-only Jacobian has rank one, but
a free source scale identifies only `Gs`, and one unrestricted remainder per
configuration exactly aliases the source column. Ten run-level, covariance,
geometry, and ownership fields remain absent, so this is a real-apparatus
readiness result rather than an independent numerical \(G\) fit.

The official HUST-2018 release goes materially beyond that PDF-only boundary.
Its time-of-swing data recover one figure-level gravitational stiffness
response, while its angular-acceleration-feedback tables reproduce all three
published processed source-response forwards to within `0.2 ppm`, without
inserting an accepted value of \(G\). This independently exercises the two
ordinary-gravity ownership positions—operator/stiffness and source/forcing—and
supports the finite-apparatus non-double-counting architecture. The released
source coefficients are already processed, however, and the mass-coordinate
geometry, complete transfer, correction/remainder ownership, source-scale
calibration, conserved apparatus-stress ledger, and joint covariance remain
absent. The result is therefore a real processed-coefficient forward and
limited figure-level response extraction, not a full `GC16` execution, a new
\(G\), confirmation of RGRL/GFT, or evidence for a lineage gravitational
charge.

Two bounded continuations sharpen what that release can and cannot establish.
Equal-configuration ToS endpoint returns yield exact return/differential ranks,
reused-endpoint Gram matrices, and telescoping identities, together with a
reproducible source-present-minus-background residual. The source-present and
background panels are separate same-ordinal acquisitions, no row covariance is
published, and no matched no-excursion arm exists. The result is therefore a
history-confound diagnostic, not causal memory, record lineage, `beta_TM`,
gravity emergence, `GC16`, or a new \(G\).

The public unprocessed mass and dimension fields also determine exact
conditional homogeneous AAF/ToS Newtonian functionals and highly converged
numerical source coefficients. Their local sensitivity classes reproduce the
published geometry classes, while the nominal coefficients differ from the
authors' processed values by roughly `41--53 ppm` for AAF and `91--103 ppm`
for ToS. Four AAF pair distances plus a centroid leave two shear degrees of
freedom, so the pairwise-centred placement is a premise rather than a unique
public reconstruction. The scalar discrepancies do not identify the omitted
clamp, coating, density, CMM, shape, transfer, or complete-stress remainder.

Composing the released corrected response fields with those nominal kernels,
without importing an accepted value of \(G\), gives three exact conditional AAF
quotients at `r_norm=0` between `6.6740227e-11` and `6.6742605e-11` SI. The
seven ToS rows give zero-fibre-correction anchors between `6.6734516e-11` and
`6.6736824e-11` SI plus exact affine correction families. These are not full-
apparatus estimates: at that calculation stage the matched mass-multipole/CMM
numerator, covariance, and signed ToS correction had not been independently
reconstructed. An equally
public normalization collision shifts the AAF quotient by about `1631.54 ppm`
and the ToS anchors by `152--153 ppm`, proving that total inertia cannot replace
the missing source multipole. Thus no public physical \(G\) point or compact
interval, new \(G\), `GC16`, or GFT confirmation is promoted.

The disclosed central correction tables now recover calibrated partial source
kernels and the signed ToS anelastic corrections, so the signed central ToS
line is no longer free. The resulting zero-independent-remainder forwards span
approximately `6.6741755--6.6744092e-11` SI for AAF and
`6.6739974--6.6742407e-11` SI for ToS; processed-kernel comparator gaps narrow
to about `19--30 ppm` and `0.6--12.1 ppm`. The authors' processed kernel makes
the remaining scalar inferable but does not independently own its physical
provenance. One independently reconstructible released-row-bound harmonic
remainder—or the physical maps needed to calculate it—would unlock an
independent source-side row evaluation using the released response summary. A
hostile-audited bounded search found no such
public root on its declared surfaces and retained two dissertation leads only
for targeted acquisition. The search is non-exhaustive and advances no
numerical \(G\).

A separate four-file public mechanical archive establishes repeated-input,
same-apparatus path dependence across 101,628 samples. It does not contain the
authenticated lineage, randomized paired KEEP/BREAK intervention, matched
conventional state, common future query, covariance, or gravity observable
needed for a causal lineage or gravity estimate. Those two data blockers and
the calculations an admitted payload would unlock are now stated in
[`NO_LAB_GRAVITY_DATA_ACQUISITION_PACKET_V001.md`](NO_LAB_GRAVITY_DATA_ACQUISITION_PACKET_V001.md).

The public five-clock component has likewise been advanced without assigning
it more evidential weight than it owns. Its ten separately reanalyzed processed
pair heights form a complete K5 edge vector. The cut rank is four, the cycle
dimension six, and all 37 simple cycles were checked. Exact marginal-box
optimization gives \(\rho_*=27/82<1\), so one five-node scalar intersects every
reported marginal `1 sigma` interval. The mapping to height already assumes GR
and common local \(g\), the pair estimates share clocks and runs, and no joint
covariance is deposited. This is deterministic processed node-scalar
compatibility, not a chi-square, joint coverage statement, independent GR or
common-metric test, lineage result, or gravity-emergence evidence.

The original `AFR` theorem remains one conditional raw-front refinement route.
A distinct audited collective route now proves the exact finite identity
\(B_N^\dagger B_N=4I+A_N\), the `A3` root second moment
\(\sum_{a<b}\alpha_{ab}\alpha_{ab}^{\mathsf T}=16I/3\), the cell covolume
\(16a_*^3/(3\sqrt3)\), and a refining mathematical affine three-manifold
atlas. A supplied repeatable hop/next-slab-onsite schedule now gives an exact
child-only phase, uniform parent/child Floquet quasienergy separation, and an
exact dressed-parent function of \(B_N^\dagger B_N\). It does not supply a
static source-off stagger or an autonomous phase. A separately supplied
shift-symmetric massless action on that
support has the smooth infrared principal cone
\(-\omega^2+c_*^2|k|^2=0\). The autonomous detuned phase, physical metric solder,
massless phase, finite-speed clock scaling, common probes, physical volume,
tensor constraints, stress vertex, and gravity remain open.  A separate
Clifford stencil theorem supplies an exact prospective first-order Weyl
candidate on the same `A3` roots; it is not current F3 or gravity.

The two consecutive q4 front cosets also have a newly exact support-shape
identity: their translation-invariant append-incidence completion is the
standard diamond net, with `A3`/FCC translations, two cosets, tetrahedral
degree four, girth six, and exact local exhaustion by deep nonnegative q4
slabs.  This removes the previously independent diamond-shape premise from the
existing `d_*=2` F3 calculation.  The later finite programmed support-solder
binds the deep-interior local edges, where the leading operator has
\(J_6=63h^6/(8U_d^5)\) and \(V_6=0\).  The raw slab has no global `d_*=2`
sector, so the published `mu=0` spin-one U(1)-liquid comparator still requires
a supplied regular completion. Autonomous support selection/stability, an
unbounded or selected periodic realization, volume-uniform all-orders control,
visible electromagnetism, tensor gravity, and RGRL-B remain open.

The current-parent boundary is now sharp. On any supplied saturated q4 support,
the unchanged one-carrier F3 restriction already gives the exact scalar
incidence transfer \(H_1=C I+\epsilon_\psi I+\lambda_JD-tA\); no new hopping law
is needed for that restricted result. Current BQ4 counters and append keys do
not, however, constitute coexisting F3 sites and edges, and the source-off
parent owns no static positive child/parent detuning. The programmed Floquet
gap closes that finite controlled operation only. Moreover, the full degree-four
carrier slice and the `d_*=2` diamond-ice slice are exactly disjoint when they
use the same binary incidence field. Six lawful FPMH records can supply a
static `S4` pair-register representation, but inherited dynamics conserve every
pair projector, so all interpair retarded kernels vanish in that register
construction. On the physical programmed q4 links, however, the same six Walsh
functions now have exact finite local and shared-link response under inherited
link dynamics. This is a physical Walsh-operator-algebra realization, not a
physical PMMDC solder, automatic record qualification, or a metric field;
neither a new interaction nor a second field has been adopted.

The finite site/edge portion of that support solder is nevertheless now closed
as an exact **programmed preparation**. Given a supplied finite F3 array,
address maps, q4 edge list, cap, sources, controller schedule, and complete
ports, existing FPMH/PESC gates reversibly write and hold the q4 support word,
quarantine every nonedge and guard, and make BS09 propagate on the exact finite
q4 incidence matrix. During the claimed invariant hold, both the raw BS06 flip
and the PESC `K`-gated incidence flip must be off or continuously cancelled; a
stroboscopic return alone does not preserve the carrier evolution. This removes
a finite hardware-type obstruction and
materially advances the acoustic route. It does not make one BQ4 front into
many physical sites, autonomously select or stabilize the graph, derive the
static child detuning, or establish a scalable phase. The raw finite q4 slab
also has an extreme child of degree one, so its global `d_*=2` ice sector is
exactly empty; only deep-interior local diamond inheritance survives without a
separately supplied regular boundary completion.

On a compatible completed ice domain, the physical one-link sector is exact
`T2`, the centered pair tangent is exact rank-two `E`, and pair `A1` is fixed.
FJ nevertheless establishes an exact conditional rank-six spectral response
for the six unprojected active pair operators under its programmed `H_resp`,
including nearest-cell response and finite spreading. That is a real
six-channel response result, but it does not identify the `S4` edge module with
one continuum symmetric tensor or supply its Ward/constraint packet. The local
Fisher screen remains negative: its complement-broken rank-six covariance map
borrows a vector-mean dyad, and complement symmetrization removes that term.

The inherited finite H6 ring model has now been screened exactly on one
180-state translation-closed periodic sector. Its zero-frequency connected
two-`Q` susceptibility and four-`Q` cumulant are nonzero, and the scalar
composite-source Legendre quartic is positive in both susceptibility
eigenchannels. The essential hostile-audit repair is its type: because each
`Q` is bilinear, the four-`Q` cumulant is an eight-one-link object. Composite
Legendre amputation does not amputate four independent one-link legs or perform
a two-particle-reducible subtraction. The four selected finite composite poles
therefore establish a finite composite precursor only. The lowest is above the
finite two-one-link threshold proxy and is not energy-exclusive from one-link
response, so no below-proxy or uniquely tensor-like candidate was identified;
this is not a thermodynamic no-bound-state theorem.

The projected q4/F3 ice branch inherits exactly one local constraint species:
the scalar Abelian compact-\(U(1)\) Gauss law, with incidence rank
\(|V|-1\) on a connected finite graph and one global dependency. The closed H6
and H8 rings preserve it; the unprojected one-link tunnelling term does not.
The four apparent pair relations factor through that scalar law and are zero
operators on the ice fiber, not four independent first-class generators.
Neither the matching `S4` module count nor the equality of final polarization
counts supplies the independent rank-three vector plus scalar constraint
packet required by RGRL-B. This is a finite/current-parent result, not a
thermodynamic no-go.

The fixed-parent collective inventory now separates static and dynamical rank
exactly. The six q4 root dyads span the full static co-metric coefficient
tangent,
\[
\operatorname{rank}(\delta c\mapsto\delta B)=6,
\]
but on the maximally favorable translation-complete realization of FI's
inherited functional form \(H_F=f_F(K)\), the corresponding dressed root
sources commute with one another and with \(H_F\), so
\[
\operatorname{rank}\chi^R_{\mathcal S\mathcal S}=0.
\]
This does not erase FJ's distinct conditional rank-six response. It shows that
no presently constructed same-parent object simultaneously owns six collective
tensor configuration channels, a nondegenerate conjugate response, and the
independent vector-plus-scalar null structure. Adopted RGRL-B supplies that
architecture as a working law, not as its microscopic F3 derivation.

The first `Q4-BLOCK-STRAIN-CTP` prerequisite has now been executed rather than
left as a generic next calculation. Four tetrahedral one-edge dyads span exact
`A1+T2` rank four and have the two-dimensional diagonal-traceless `E` sector as
their linear source null. Every FQ17a additive multi-edge weight remains in
that span. Fixed source-before-Feshbach reduction preserves the null, and a
general `O(j^2)` contact may alter the `E` Hessian but cannot manufacture a
linear source-off `Q_E` or its spectral residue.

The complete branch-relative census freezes the covering-matched periodic
family `L=5*2^r` and the reduced CW/FM pure-incidence Hamiltonian. Its complete
prospectively frozen microscopic linear source has exact rank four, while every
projected word through order eight has rank at most four. The occurrence-one
degree-square tensor is one lawful prospective additive query, not a tensor
uniquely forced by the source-free Hamiltonian. Carrier, storage, formation,
feedback, boundary, controller, and port sectors are excluded by the selected
reduced parent rather than proved absent in a physical BS completion. Thus the
unreduced physical source remains underdetermined.

The constitutive boundary has now advanced beyond a named `DPAR` candidate.
Existing degree/link/ice pair operators carry the missing `E` query and have
nonzero block-local H6 matrix elements, while `H[0]` alone still does not choose
their geometric derivative. A tetrahedrally symmetric lumped capacitance is
exactly `A1`-only, but a fully grounded pair-resolved four-terminal elastance
conditionally reproduces `U_d(d-2)^2`. A complete central kernel yields
`lambda=r0 V'(r0)/(2V(r0))`; ideal fixed-coupling Coulomb gives
`lambda=-1/2` and the corresponding conditional `U_d`--alpha--length relation.
This is a local physical completion, not yet current-F3 inheritance, visible
EM, or a global gauge-invariant shared-link solder.

Under the additional explicit `S10/FV-PURE` complete-source premise, exact
Feshbach differentiation gives projected off-shell nonidentity rank six:
direct `E` rank two plus H6 ring `A1+T2` rank four, certified by six actual
matrix-element witnesses. The first finite dynamical composition then applies
the two explicit rank-closing pieces to the FO 180-state homogeneous component.
Its exact hierarchy is family rank six, component operator rank five modulo
identity, commutator rank three, and ground retarded/first-moment rank two, with
two rank-one residues. Generated `Q_diag^(2,4,6)` derivatives and folds are not
in that witness packet, so the complete homogeneous source was calculated
separately. Every such term reduces exactly to the direct pair source plus a
Hilbert identity, giving

\[
f_E(x)=1-x^2-{37\over12}x^4-{16247\over900}x^6,
\qquad \rho_E=\rho f_E(x).
\]

The complete homogeneous hierarchy therefore remains `5 -> 3 -> 2 -> 2`
with the same two poles. The finite polynomial zero near `0.5398271903` is a
through-H6 truncation stratum, not a threshold. Retaining each source
insertion's native `A/B` vertex or link-midpoint support before Fourier
resolution now proves exact `m=0` recovery and the same three `f_E`
coefficients at `m=1` over `Q(zeta_240)`. The ring remains a nonzero
independent off-diagonal source. At the two declared samples, the finite
hierarchy is `6 -> 6 -> 6 -> 6 -> 6`, TT ground-image rank is two, and four
nonzero poles respond with residue ranks `1,3,1,1`. The nonzero spatial
contraction rules out naive transversality. FZ then proves that the projected
incidence charge is identically zero and cannot supply the missing temporal
slot, while the full Ward question remains undecided without a physical
divergence, current, and contacts. GA supplies an exact fixed-charge current
lift, GB exact local ring-energy continuity, and GC the covering family
`G_L`, exact momentum lattice, rank-two TT quotient, and
`|k_min|=3 pi/(2 L a_*)`; none of those results substitutes for a massless
pole. FY's terminal custody replay passes `107/107`, making the Phase A
source/family packet fully sealed. GD closes the shortcut's `G1` total-momentum gate by an exact
equal-and-opposite recoil construction on the declared auxiliary torus while
preserving the FY source. GF V005 independently seals the prospective `G2`
observable contract: it separates canonical amplitude from ancestry, gives a
disjoint and total PASS/FAIL/INDETERMINATE classifier with lower and upper
singular-value bounds, and requires an actual Poincare little-group
helicity-`+/-2` pole rather than a rank-two scalar doublet. The direct successor
is therefore to complete the native `G_L` H6 ledger and execute the matched
same-ancestry massless-helicity-two calculation. The consolidated no-
laboratory disposition is recorded in
[`NO_LAB_GRAVITY_PROGRAM_STATUS_V001.md`](NO_LAB_GRAVITY_PROGRAM_STATUS_V001.md).

## VALIDATION

From the repository root:

```bash
python3 model/validate_urm.py          # all four landed families, then geometry/project chain
python3 model/validate_project.py      # URM surface/delegate/D-25 gates
python3 model/validate_geometry.py     # geometry gates, then project chain
python3 model/validate_formation.py    # formation engine — 17 checks
python3 model/validate_model.py        # corner existence engine — 12 checks
python3 model/count_law.py             # exact corner count — 22 checks
python3 model/validate_udcl_postulate.py  # 32 adoption/theorem custody and ceiling checks
python3 model/validate_historywise_gravity_discriminant.py  # 84 formal-only checks
python3 model/validate_gravity_formation_theory.py  # bounded zero-input custody/ceiling gate
```

The family validator is the landing gate for new URM features. Its printed totals are computed by
the runners; documentation describes the expected composition but the exit status is authoritative.
