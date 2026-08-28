# Finite programmed q4-to-F3 support-solder theorem

**Lane ID:** `GRA-FH-F3-Q4-FPSS-V001`

**Short name:** `FPSS`

**Date:** 2026-08-27

**Claim class:** exact finite q4 slab census; exact reversible programmed
site/edge preparation using already admitted BQ4/FPMH/PESC/F3 hardware types;
exact nonedge/guard quarantine under a qualified raw-flip-free,
history-blind fixed-program hold; exact finite FD support inheritance; exact
raw-slab FE boundary obstruction

**Status:**
`FINITE_Q4_SITE_EDGE_SOLDER_PROGRAMMABLE_WITH_EXISTING_TYPES__REVERSIBLE_FIXED_ORTHOGONAL_PROGRAM_AND_COMPLETE_CUSTODY_EXACT__QUALIFIED_RAW_FLIP_FREE_NONEDGE_GUARD_QUARANTINE_EXACT__FD_SUPPORT_PART_CLOSED_BUT_POSITIVE_DETUNING_OPEN__RAW_FINITE_SLAB_GLOBAL_D2_ICE_EMPTY__FE_LOCAL_INTERIOR_SUPPORT_ONLY__AUTONOMOUS_SELECTION_SCALING_AND_SAME_N_COEXISTENCE_OPEN`

**Not claimed:** derivation of many coexisting F3 sites from one active BQ4
front; autonomous choice, correction, energetic stabilization, or scalable
production of q4/diamond support; a coherent superposition of support
programs; a derived address map, edge list, mission cap, hardware allocation,
controller schedule, source, work calibration, or physical port matrices; an
FD child/parent detuning; a global FE ice phase on the raw finite slab; a
periodic completion from BQ4; simultaneous FD saturation and FE `d_*=2` on
the same incidence field; a new `K_eT_e`, second field, graph reward, phase,
visible electromagnetism, tensor gravity, or gravity closure

## 1. Exact question and inherited boundary

The bounded q4 witness supplies finite count-front labels and append keys,
but one active BQ4 factor is not a tensor product of simultaneously existing
F3 sites.  The F3 seed supplies finite adjacent layers of qutrit carrier sites
and one link qubit for every possible adjacent-layer arrow.  FPMH/PESC supplies
an exact reversible formation-and-KEEP protocol for any **supplied** finite
edge set, and the source-off support theorem supplies passive conservation of
the prepared authenticated-support word.

The narrow question is therefore:

> Once a finite q4 slab, its physical site addresses, and its edge list have
> been prospectively supplied, can existing finite hardware and gates prepare
> that support exactly, quarantine padding and nonedges, and then hold the F3
> carrier dynamics blind to formation history?

The answer is yes for a fixed orthogonal finite program.  This is a compiler
and preparation theorem, not an autonomous support-selection theorem.

The inherited ingredients are used without changing their dynamics:

1. `BQ4RSW` supplies
   \(S_d=\{m\in\mathbb N_0^4:|m|_1=d\}\), append maps
   \(m\mapsto m+e_a\), finite counters, fresh orthogonal label slots, and
   canonical unitary dilations through a supplied cap.
2. `F3-QIRN-V001` supplies two adjacent layers with the same finite number of
   qutrit carrier sites and one binary incidence factor for every possible
   adjacent-layer arrow.
3. `PESC` supplies endpoint-owned writer slots, reversible relation
   formation, the reversible KEEP route into `K_e`, the already admitted
   `K_e X_{n_e}` gate, and a complete finite custody/port census for any
   supplied finite graph.
4. `ASSC` supplies exact source-off conservation of every prepared
   authenticated-support word and explicitly proves that this retention is
   passive, not selective or corrective.

## 2. Finite q4 slab and the exact padding census

Fix `N >= 0` and a BQ4 mission cap

\[
 R\ge \max\{4,N+1\}.                                      \tag{FH01}
\]

Write

\[
 a_N:=|S_N|={N+3\choose3},\qquad
 b_N:=|S_{N+1}|={N+4\choose3},\qquad M:=b_N.          \tag{FH02}
\]

Use one literal F3 adjacent-layer slab with `M` physical sites in each layer,

\[
 V_0=\{(0,i):1\le i\le M\},\qquad
 V_1=\{(1,j):1\le j\le M\}.                          \tag{FH03}
\]

Every displayed physical site retains the **complete** F3 node allocation
`(psi,r,w,b_1,...,b_K)` and its attached support/port factors, not merely the
carrier qutrit `psi`.

Supply and freeze injective owner/address maps

\[
 \iota_0:S_N\hookrightarrow V_0,\qquad
 \iota_1:S_{N+1}\xrightarrow{\sim}V_1.              \tag{FH04}
\]

Thus the child layer uses all of its sites, while the parent layer contains

\[
 g_N=M-a_N=b_N-a_N={N+3\choose2}                     \tag{FH05}
\]

padding/guard sites.  The complete adjacent-layer link domain is

\[
 D_N=V_0\times V_1,\qquad |D_N|=M^2.                \tag{FH06}
\]

The prospectively supplied q4 edge list is

\[
 E_N=\{(\iota_0(m),\iota_1(m+e_a)):
          m\in S_N,\ a\in\{1,2,3,4\}\}.             \tag{FH07}
\]

Every pair `(m,a)` gives a distinct directed edge, so

\[
 |E_N|=4a_N,\qquad |D_N\setminus E_N|=M^2-4a_N.     \tag{FH08}
\]

Every active parent has eligible degree four.  A child `c in S_(N+1)` has

\[
 \deg_{E_N}(c)=|\{a:c_a>0\}|.                       \tag{FH09}
\]

All guard-parent edges and all cross-layer pairs not listed in (FH07) are
nonedges.

### Theorem `FPSS-1` -- exact finite site/edge census

Equations (FH02)--(FH09) give an exact embedding of the finite q4 append
incidence slab into already admitted finite F3 site/link allocation, using
only parent padding.  No q4 label is identified with another, no physical
site is counted twice, and no nonedge is deleted from the raw Hilbert space.

#### Proof

Stars-and-bars gives (FH02), and Pascal's identity gives (FH05).  If two
members of (FH07) have the same ordered endpoints, injectivity of `iota_0`
gives the same `m`, then injectivity of `iota_1` gives `m+e_a=m+e_b`, hence
`a=b`.  This proves (FH08).  A child has the parent `c-e_a` exactly when
`c_a>0`, proving (FH09).  QED.

## 3. Exact finite programmed solder

For every `e in D_N`, retain the F3/PESC link qubit and reuse it literally as
the FPMH link factor,

\[
 n_e:=a_e,                                                   \tag{FH09a}
\]

together with the authenticated pair-memory factors `(L_e,K_e,G_e)`.  For
every `e in E_N`, allocate the two
distinct endpoint-owned FPMH writer qutrit slots

\[
 A_e\equiv(u,e),\qquad B_e\equiv(v,e),\qquad e=(u,v), \tag{FH10}
\]

with owner map fixed by the displayed physical F3 endpoints.  Nonedges retain
blank `L/K/G/n` factors but receive no formation token.

PESC's full generic census may additionally contain FPMH factors for
within-layer unordered pairs.  They remain explicit blank spectators.  They
are not identified with nonexistent within-layer F3 incidence links.  The
compiler uses the tensor-local one-pair FPMH construction only on the
prospective bipartite domain `D_N=V_0 x V_1`.

The fixed program `p_N` consists of the supplied cap, address maps, edge list,
hardware allocation, source-token list, orientation, route table, controller
schedule, and port census.  It is one orthogonal controller/program basis
state.  It is not summed coherently over different graphs.

Place one declared orthogonal source token in every `A_e`, `e in E_N`, with
`B_e` blank.  Apply the already admitted FPMH product formation unitary

\[
 U_{\rm form}(E_N)=\bigotimes_{e\in E_N}U_{A_e\to B_e},       \tag{FH11}
\]

followed by the already admitted reversible parallel KEEP route

\[
 U_{\rm KEEP}:L_e\longrightarrow K_e.                       \tag{FH12}
\]

The factors in (FH11) have disjoint writer slots and targets.  PESC therefore
gives the exact support output

\[
 K_e=\mathbf1_{E_N}(e)\quad(e\in D_N),                       \tag{FH13}
\]

while preserving every writer, route, source, work, outcome, failure,
boundary, quarantine, and untouched-reference factor.

If saturated incidence is required, start every `n_e` blank and apply the
already admitted finite controlled pulse

\[
 U_{KX}=\prod_{e\in D_N}
 \exp\!\left({i\pi\over2}P^K_eX_{n_e}\right).                \tag{FH14}
\]

On a `K_e=0` block (FH14) is the identity.  On a `K_e=1` block it is
`iX_(n_e)`, so, up to an irrelevant product phase,

\[
 |K=\mathbf1_{E_N}\rangle|n=0\rangle
 \longmapsto
 |K=\mathbf1_{E_N}\rangle|n=\mathbf1_{E_N}\rangle.          \tag{FH15}
\]

During the ideal pulse (FH14), the raw ungated BS06 link flip and every
noncommuting incidence/carrier term--including diagonal incidence detuning,
degree-return, current, and incidence-gated hopping terms--are scheduled off
or exactly refocused;
otherwise their simultaneous evolution would not equal (FH14).  Only the
already admitted `K`-gated generator acts for the calibrated pulse interval.
This switching/refocusing choice and its controller, clock, drive, and work
ownership are part of supplied `p_N`, not a new Hamiltonian term.  The
declared hold Hamiltonian is restored only after the pulse ends.

Define

\[
 U_{\rm solder}(p_N)=U_{KX}\,U_{\rm KEEP}\,
                      U_{\rm form}(E_N),                    \tag{FH16}
\]

with `U_KX` omitted when only the eligibility word is to be prepared.

### Theorem `FPSS-2` -- exact reversible fixed-program compiler

For every fixed finite `N` and supplied orthogonal program `p_N`, (FH16) is a
finite reversible, reference-stable, custody-complete preparation of the q4
support word (FH13), and optionally of the saturated incidence word (FH15),
using only already admitted finite gate and register types.

#### Proof

Each FPMH formation factor has its inherited unitary dilation; distinct edge
factors commute because their writer/link targets are distinct.  The KEEP
route is reversible on retained `L/K/G` and route-history factors.  Each
factor in (FH14) is the exponential of a finite Hermitian controlled Pauli and
is unitary.  Their ordered product is therefore unitary.  Tensoring every
stage with the identity on the untouched reference proves reference
stability.  No input or output register is traced, reset, or identified, so
the inherited finite source/controller/work/failure/boundary census supplies
logical and custody port completeness.  Physical energies, port matrices, and calibration remain supplied components of `p_N`; this theorem does not
derive or measure them.  QED.

## 4. Exact padding/nonedge quarantine and history-blind hold

After (FH16), turn off source, writer, route, copy, and preparation pulses and
quarantine the old source.  Retain all their factors.  Throughout any hold in
which the explicit `D_N\setminus E_N` nonedges are required to remain blank,
the fixed program also keeps the **raw ungated** BS06 term
`-h_N sum_(e in D_N) X_e` exactly zero in the hold generator, either by
switching it off or by continuous exact cancellation.  A merely
stroboscopic echo which leaves and later returns to the blank sector does not
establish the instantaneous invariant-subspace statement below; under such an
echo only the declared sampling-time return map, not (FH17)--(FH19) as
Hamiltonian identities, would be earned.  If an incidence-changing actuator
is used in the qualified hold, it is the already admitted PESC
`-h sum_e P^K_e X_e`, whose controller/work schedule is retained.  This
qualification is necessary: the raw ungated flip would take a
`K_e=0,n_e=0` nonedge out of the blank block.

All remaining terms are already admitted PESC/F3 terms.  Every nonedge
satisfies `K_e=0,n_e=0`.  The `K_eX_(n_e)` term vanishes there; degree and
incidence energies are diagonal in `n`; and the BS09/BS11 carrier terms
contain `n_e` as a multiplicative control.  Hence

\[
 \mathcal H_{\rm quar}=
 \left(\bigotimes_{e\in D_N\setminus E_N}
 \operatorname{span}\{|0\rangle_{K_e}|0\rangle_{n_e}\}\right)
 \otimes\mathcal H_{\rm rest}                               \tag{FH17}
\]

is invariant.  Every guard carrier begins blank.  It has no occupied
incident link, the carrier onsite term preserves blankness, and its formation
couplings are off, so the all-blank guard-carrier factor is also invariant.
The degree energy of a degree-zero guard may be a supplied scalar within this
fixed admitted block; it is not a claim that guards are selected by an energy
minimum.

Let `H_hist` include every BQ4 order/provenance/scaffolding/routing factor and
every compiler formation/route history.  Let `Pi_(p_N,E_N)` fix the orthogonal
program and authenticated `K` support word, but not the subsequently active
`n` word.  On that fixed support block, the qualified hold has the exact form

\[
 [H_{\rm hold},\Pi_{p_N,E_N}]=0,
 \qquad
 \Pi_{p_N,E_N}H_{\rm hold}\Pi_{p_N,E_N}
 =\Pi_{p_N,E_N}\otimes I_{\rm hist}
  \otimes H_{\rm F3}^{(N)}(E_N),                            \tag{FH18}
\]

up to independent port Hamiltonians that act as identity on support and
formation history.  This factorization is a premise of the qualified
fixed-program hold: every history-writing coupling is off and the remaining
controller/port evolution is independent of formation history.  It is not a
claim about an arbitrary source-off F3 Hamiltonian.  ASSC additionally gives

\[
 [H_{\rm hold},P^K_e]=0,
 \qquad p_s(\tau)=p_s(0)                                   \tag{FH19}
\]

for every authenticated-support word `s`.

### Theorem `FPSS-3` -- exact quarantine and passive hold

The programmed q4 `K` support, blank nonedges, and blank padding carriers are
exact invariant sectors of the qualified source-off hold, and carrier/link
evolution is blind to the retained BQ4 and compiler histories as in (FH18).
The optional saturated `n` word is not asserted to remain fixed when the
`K`-gated actuator is on.  This is passive support retention only.  Equation
(FH19) conserves every wrong, dense, disconnected, and non-q4 `K` support
word just as exactly; it supplies no selection, correction, attraction basin,
or scalable preparation law.

The fixed-orthogonal-program premise is load-bearing.  If different support
programs were coherently superposed, the subsequent controlled F3 evolution
would generally entangle the carrier with the program/support register and
would not be the single history-blind block (FH18).

## 5. What the solder closes for the FD lane

On the saturated output (FH15), restrict to one carrier of one fixed content,
exclude the invariant blank guards, and set the optional diagonal comparator
coefficients to their already lawful zero values, in particular
`lambda_J=0`.  During this FD carrier comparator, keep both the raw ungated
BS06 flip and the PESC `K`-gated incidence flip exactly zero in the generator.
All retained incidence terms are then diagonal in `n`, so the saturated `n`
word is an invariant block rather than merely an instantaneous preparation
output.  The BS09 off-diagonal
carrier block between the two active layers is exactly

\[
 H_{\rm hop}^{(N)}=-t
 \begin{pmatrix}
 0&B_N^\dagger\\ B_N&0
 \end{pmatrix},
 \qquad
 (B_N)_{cm}=1\Longleftrightarrow c=m+e_a.                    \tag{FH20}
\]

With retained diagonal terms, the same statement holds for the off-diagonal
block, accompanied by the explicitly known degree/onsite diagonal rather
than a silently assumed constant.

### Theorem `FPSS-4` -- finite FD support part is constructible

The finite site/edge portion of `Q4-CARRIER/EDGE-LIFT` is an exact programmed
construction with existing types: BS09 propagates on the physical q4 append
incidence matrix (FH20).  The theorem does **not** supply the positive uniform
child/parent detuning `Delta` required by the FD Schur-complement/acoustic
construction, nor its work/support/maintenance ports.  It also does not prove
support energetics, an `N`-uniform scalable family in nature, a collective
phase, or the infrared acoustic action.

## 6. Exact obstruction for a global FE ice sector on the raw slab

Define the degree-two incidence set on the whole finite bipartite slab by

\[
 \Omega_2(E_N)=\{n\in\{0,1\}^{E_N}:
                  d_v(n)=2\ \text{for every active vertex }v\}.             \tag{FH21}
\]

The extreme child

\[
 c_*=(N+1,0,0,0)                                    \tag{FH22}
\]

has only the eligible parent `(N,0,0,0)`.  Therefore every subgraph of `E_N`
has `d_(c_*)(n)<=1`, and

\[
 \boxed{\Omega_2(E_N)=\varnothing.}                 \tag{FH23}
\]

### Theorem `FPSS-5` -- FE receives a local support bridge, not a raw-slab phase

The compiler closes the finite physical binding of eligible q4 edges to F3
link registers.  For slabs large enough to contain the required translated
neighborhood, the already proved deep-interior diamond-net adjacency and
local linked-hexagon operator algebra are therefore physically represented;
any locked-sector use still needs a supplied compatible boundary state or
completion.  Equation (FH23) forbids the global `d_*=2` locked ice manifold
on the raw `S_N--S_(N+1)` slab.  A periodic diamond quotient or other degree-four
boundary completion can itself be programmed with the same finite hardware
types, but its extra identifications/edges are supplied boundary physics and
are not derived from BQ4 append lineage.  Thermodynamic/all-orders U(1),
support stability, visible electromagnetism, and actual-world coexistence
remain open.

## 7. Exact FD/FE same-`n` ceiling

On every active parent, the FD slice (FH15) has

\[
 d_m(n)=4,                                            \tag{FH24}
\]

whereas the FE ice slice requires `d_m(n)=2`.  Thus even on a separately
supplied regular/periodic completion,

\[
 \{n:n_e=1\ \forall e\in E\}
 \cap\{n:d_v(n)=2\ \forall v\}=\varnothing.          \tag{FH25}
\]

The finite solder does not cure the existing same-`n` incompatibility.  FD
and FE remain distinct programmed runs or sectors.  No `K_eT_e`, second
kinetic/support field, or other new dynamics is introduced here.

## 8. Exact result and remaining antecedents

The strongest earned statement is:

\[
\boxed{
\begin{gathered}
\text{A supplied finite q4 slab can be reversibly programmed into}\
\text{the existing F3/PESC site, support, and link hardware with exact}\
\text{padding/nonedge quarantine in a qualified raw-flip-free,}\
\text{fixed-program history-blind hold.}
\end{gathered}}                                                   \tag{FH26}
\]

This removes a possible **finite hardware-type obstruction** to the site/edge
solder.  It does not remove the physical antecedents that were supplied to
the compiler:

1. the coexisting finite F3 site/link array and its owner/address map;
2. the q4 edge list, mission cap, source tokens, controller schedule, and
   complete source/work/support/boundary/failure/quarantine port realization;
3. autonomous support selection, correction, energetic stability, and a
   scalable or thermodynamic preparation law;
4. the FD positive uniform child/parent detuning and its physical ownership;
5. a non-BQ4-supplied periodic/regular FE completion and a global/all-orders
   U(1) phase; and
6. any simultaneous FD/FE parent, visible EM identification, tensor response,
   universal stress coupling, or gravity theorem.

In particular, (FH26) does not map one BQ4 count-front state into many F3
sites, and it does not implement an arbitrary coherent isometry from `Q_N`
into a one-carrier F3 space.  The many sites are supplied physical hardware,
and the q4 labels are prospectively soldered to them by the fixed-basis
address program.  This is the precise difference between a lawful finite
compiler and autonomous emergence.
