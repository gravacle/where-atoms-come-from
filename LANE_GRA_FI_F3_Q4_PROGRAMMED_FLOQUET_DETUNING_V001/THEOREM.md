# Programmed F3/q4 Floquet carrier-detuning theorem

**Lane ID:** `GRA-FI-F3-Q4-PFCD-V001`

**Short name:** `PFCD`

**Date:** 2026-08-27

**Claim class:** exact admitted-generator and schedule theorem; exact
child-only carrier onsite pulse; exact finite Floquet spectrum and uniform
quasienergy separation; exact dressed-parent q4 common-child functional
calculus with an operator-norm kernel bound; exact static-versus-programmed
detuning boundary on a qualified dual-flip-free saturated-incidence slice

**Status:**
`STATIC_SOURCE_OFF_STAGGER_REMAINS_ABSENT__EXISTING_F3_ONSITE_PLUS_NEXT_SLAB_SCHEDULE_GIVES_EXACT_CHILD_ONLY_PHASE__DUAL_RAW_AND_K_GATED_FLIPS_ZERO_ON_SATURATED_INCIDENCE_SLICE__FINITE_TWO_PULSE_FLOQUET_PARENT_BRANCH_UNIFORMLY_GAPPED__EXACT_Q4_COMMON_CHILD_FUNCTIONAL_AND_CONTROLLED_SIBLING_KERNEL__NO_NEW_BULK_CARRIER_GENERATOR_OR_REGISTER_PORT_TYPE__POSITIVE_EPSILON_EXTRA_BLANK_LAYER_REPEATABLE_OWNED_PROGRAM_CONTROL_MATRICES_AND_CALIBRATION_SUPPLIED__AUTONOMOUS_DETUNING_PHASE_RGRL_B_AND_GRAVITY_OPEN`

**Not claimed:** that the unchanged source-off F3 bulk contains an always-on
staggered onsite term; that the static Hamiltonian in `CCMAC` equation (FD05)
has been derived; that the controller schedule is autonomous; that switching
costs vanish; that the numerical values of `epsilon_psi`, `t`, or either pulse
duration are predicted; that a collective acoustic phase, physical proper
time, visible electromagnetism, a tensor response, RGRL-B, gravity, or `G` has
been derived; or that a merely stroboscopic return of the incidence word is
sufficient to preserve the carrier unitary while incidence-gated hopping is
active.

## 1. Frozen question and inherited exact results

`CCMAC` showed that q4 append incidence

\[
 B_N:\ell^2(S_N)\longrightarrow\ell^2(S_{N+1}),
 \qquad (B_N)_{cm}=1\Longleftrightarrow c=m+e_a,
 \tag{FI01}
\]

obeys

\[
 K_N:=B_N^\dagger B_N=4I+A_N,                     \tag{FI02}
\]

where `A_N` is the q4 sibling adjacency.  A static positive child offset then
gives the exact Schur/common-child branch.

`CLDNG` proved two facts that remain in force:

1. on a supplied saturated q4 support, the unchanged F3 one-carrier block
   already gives the off-diagonal scalar transfer in (FI01); and
2. the regular source-off bulk is bipartition-exchange symmetric and contains
   no positive child-only onsite term.

`FPSS` subsequently supplied an exact finite programmed realization of the
q4 site/edge support using existing BQ4/FPMH/PESC/F3 register and gate types.
It also established the legitimacy of a fixed orthogonal controller program
which switches or, where instantaneous invariance is required, continuously
cancels noncommuting terms while retaining the complete clock, work, boundary,
failure, and quarantine census.

The narrow remaining question is:

> Can the existing F3 bulk generator set and the already admitted
> fixed-program control architecture create a positive uniform child/parent
> separation without adding a new bulk carrier interaction, even though the
> static source-off bulk cannot?

The answer is **yes stroboscopically**.  The existing layer-local carrier
onsite term produces an exact child-only phase when the child layer is reused
as the first layer of the next F3 slab.  A finite hop/phase cycle then has an
exact isolated parent-connected quasienergy band and generates the same q4
common-child kernel at leading controlled order.  This is a programmed
Floquet realization, not a retroactive derivation of the static FD05 block.
Here "stroboscopic" describes the carrier quasienergy law sampled after each
complete two-pulse cycle; it does not permit the supporting incidence word to
leave its fixed block within either pulse.

## 2. Exact finite physical sector

Fix `N >= 0`.  Use the `FPSS` construction with

\[
 P\cong\ell^2(S_N),\qquad
 C\cong\ell^2(S_{N+1}),\qquad
 M=|S_{N+1}|,                                      \tag{FI03}
\]

inside two consecutive F3 layers `V_0,V_1`, each containing `M` physical
sites.  The unused sites of `V_0` are the already quarantined parent guards.
Add the next admitted F3 layer `V_2`, also of width `M`, and call its carrier
space `G`.  Every `G` carrier starts in the unique blank, every `V_1--V_2`
incidence factor starts blank, and no formation token is supplied there.

Work in one fixed carrier-content sector `x in {0,1}` with exactly one carrier
in `P \oplus C`.  The content is an inert identity multiplicity.  Freeze
the authenticated q4 support and the saturated `V_0--V_1` incidence word from
`FPSS`; end the formation, writer, route, and copy pulses.  During **both** the
hop pulse and the next-slab onsite pulse, keep the raw ungated BS06 incidence
flip and the PESC `-h sum_e P_e^K X_{n_e}` actuator exactly zero in the joint
generator, either switched off or continuously exactly cancelled.  This
condition applies to the saturated old-slab links and the blank next-slab
links.  A merely stroboscopic incidence echo which leaves and later restores
the `n` word is insufficient: while `n` is changed, BS09 implements a
different carrier Hamiltonian, so the joint return need not equal (FI11) or
(FI12).  Take the already lawful comparator value `lambda_J=0`.  Every other
retained fixed storage, incidence, and independent-port term is either a
common scalar on this block, is continuously cancelled by the supplied
program, or factors into an independent spectator unitary.  Nothing is
traced.

On `P \oplus C`, define

\[
 X_B=
 \begin{pmatrix}0&B_N^\dagger\\B_N&0\end{pmatrix},
 \qquad
 \Pi_C=\begin{pmatrix}0&0\\0&I_C\end{pmatrix}.
 \tag{FI04}
\]

During a `V_0--V_1` carrier pulse, BS09 gives, after removal of the uniform
one-carrier onsite scalar,

\[
 H_H=-tX_B.                                        \tag{FI05}
\]

The dual-flip-free premise makes the saturated old-slab `n` word an invariant
block throughout this pulse; (FI05) is therefore an exact joint-sector
restriction, not an instantaneous matrix evaluated at the prepared word.
No q4 sibling edge has been inserted: `X_B` is the physical saturated append
incidence block already earned by `FPSS`.

## 3. The existing onsite term gives an exact child-only pulse

The BS09 carrier onsite term in an active slab is

\[
 H_{\rm on}^{(\ell)}
 =\epsilon_\psi\sum_{v\in V_\ell\cup V_{\ell+1}}q_v^\psi.
 \tag{FI06}
\]

Assume the already admitted storage-domain value

\[
 \epsilon_\psi>0.                                  \tag{FI07}
\]

Now activate the next `V_1--V_2` slab while the `V_0--V_1` carrier transfer is
off.  The blank `V_1--V_2` incidence word makes its BS09 hopping identically
zero.  Both the raw incidence flip and the PESC `P^KX_n` actuator are exactly
zero in the joint generator during this pulse, as required in Section 2, and
all copy/formation terms are off.  Hence the old saturated incidence and the
next-slab blank incidence both remain fixed throughout the onsite interval.
Because `G=V_2` is carrier blank and the sole carrier remains in
`P \oplus C`, (FI06) restricts exactly to

\[
 \boxed{H_D=\epsilon_\psi\Pi_C.}                  \tag{FI08}
\]

Thus the child/parent energy contrast during this pulse is uniform and
positive.  It comes from the existing F3 carrier onsite energy, not from
BS06 link detuning, BS11 record feedback, a symbolic BS12 completion, or a
new stagger interaction.  The asymmetry is in the prospectively supplied
composition schedule and blank boundary state, not in a node-label-dependent
bulk coefficient.

Choose the finite onsite-pulse duration

\[
 \tau_D={\pi\hbar\over2\epsilon_\psi}.
 \tag{FI09}
\]

After factoring the independent spectator evolution, its exact carrier
unitary is

\[
 \boxed{U_D=e^{-i\tau_DH_D/\hbar}=\Pi_P-i\Pi_C.}   \tag{FI10}
\]

### Theorem `PFCD-1` -- exact scheduled detuning pulse

For every finite `FPSS` q4 slab with `epsilon_psi>0` and one additional blank
F3 layer, the unchanged BS09 onsite generator plus the admitted next-slab
schedule implements (FI10) exactly on the fixed one-carrier,
dual-flip-free-incidence block.  No new bulk carrier generator,
field/register type, or port type is required.  The concrete controller
couplings which isolate and time the pulse remain supplied physical matrices
within that admitted port architecture; they are not proved by the symbolic
BS12 slot.

This is not an energy-free control theorem.  Selecting the active slab,
switching or continuously cancelling the noncommuting terms, timing (FI09),
and returning to the carrier pulse exchange whatever work the physical
implementation requires
with the already retained controller/clock/work/recoil ports.  Their physical
matrices and calibration remain supplied; the energy term being applied is
owned once by (FI06).

## 4. Exact two-pulse Floquet spectrum

Let `tau_H>0` be a calibrated carrier-pulse duration and put

\[
 \eta={t\tau_H\over\hbar},
 \qquad
 U_H=e^{-i\tau_HH_H/\hbar}=e^{i\eta X_B}.           \tag{FI11}
\]

The uniform BS09 onsite energy on `V_0 \cup V_1` contributes only a common
phase during this one-carrier pulse and has been factored out.  One programmed
cycle, in physical time order hop then child phase, is

\[
 \boxed{U_F=U_DU_H.}                               \tag{FI12}
\]

The same finite orthogonal program may be repeated.  Let `T_F>0` be its
calibrated full cycle time, including any exact joint-`n`-and-carrier-identity
switching or continuous-cancellation intervals.  A support echo which is
identity on `n` only at the end is not a carrier-identity interval and cannot
be inserted into (FI12) without a fresh joint-unitary calculation.  Repetition
is a supplied controller action, not an autonomous source-off law.

Every column of `B_N` has four ones and every row has at most four ones, so

\[
 \|B_N\|\le\sqrt{\|B_N\|_1\|B_N\|_\infty}\le4.    \tag{FI13}
\]

Moreover `B_N` is injective.  If

\[
 p(z)=\sum_m p_m z^m
\]

is the homogeneous polynomial whose coefficients are a parent vector, then
`B_Np` is the coefficient vector of
`(z_1+z_2+z_3+z_4)p(z)`.  The polynomial ring is an integral domain, so
`B_Np=0` implies `p=0`.

Freeze the pulse ceiling

\[
 \boxed{0<|\eta|\le{\pi\over16}.}                 \tag{FI14}
\]

For a singular value `sigma_j \in (0,4]` of `B_N`, put

\[
 a_j=\eta\sigma_j,
 \qquad
 \omega_j=\arccos\!\left({\cos a_j\over\sqrt2}\right).
 \tag{FI15}
\]

In the corresponding parent/child singular-vector pair, (FI12) is

\[
 U_{F,j}=
 \begin{pmatrix}
  \cos a_j&i\sin a_j\\
  \sin a_j&-i\cos a_j
 \end{pmatrix},                                   \tag{FI16}
\]

and has the exact eigenvalues

\[
 \boxed{
 \lambda_{j,P}=e^{i(\omega_j-\pi/4)},
 \qquad
 \lambda_{j,C}=e^{-i(\omega_j+\pi/4)}.}            \tag{FI17}
\]

Every unpaired vector in `ker B_N^dagger \subset C` is a dark child mode with
eigenvalue `-i`.

From (FI13)--(FI14), `|a_j| \le pi/4`, and hence

\[
 {\pi\over4}\le\omega_j\le{\pi\over3}.            \tag{FI18}
\]

The parent-connected eigenphase lies in `[0,pi/12]`; the coupled child phase
lies in `[-7pi/12,-pi/2]`; and the dark child phase is `-pi/2`.  Therefore the
parent-connected band is separated from every child band on the principal
quasienergy branch by at least

\[
 \boxed{\Delta_F\ge{\pi\hbar\over2T_F}>0.}         \tag{FI19}
\]

### Theorem `PFCD-2` -- exact finite positive quasienergy separation

Under the qualified dual-flip-free isolation of Section 2 and (FI07), (FI09),
and (FI14), the programmed two-pulse unitary has one parent-connected dressed
band and only child-connected/dark bands across the uniform positive
separation (FI19).  The result is exact at every finite `N`.  It is a
quasienergy theorem for the repeated schedule, not a claim that the source-off
static Hamiltonian acquired (FF17).

## 5. Exact dressed-parent common-child dynamics

Let

\[
 \iota_P:P\longrightarrow P\oplus C,
 \qquad \iota_Pp=(p,0),                            \tag{FI19a}
\]

be the canonical parent inclusion.  Let `W_P(eta):P -> P \oplus C` be the
spectral isometry which maps each parent right singular vector continuously
to the normalized `lambda_(j,P)` Floquet eigenvector, so that
`W_P(eta)^dagger W_P(eta)=I_P`, and fix its branch and phases by
`W_P(0)=iota_P`.  Pulling the principal Floquet logarithm back through `W_P`
gives the exact parent-connected quasienergy operator

\[
 \boxed{
 H_{P}^{F}=f_F(K_N),
 \quad
 f_F(K)=-{\hbar\over T_F}\left[
 \arccos\!\left({\cos(\eta\sqrt K)\over\sqrt2}\right)
 -{\pi\over4}I\right].}                            \tag{FI20}
\]

This is an exact functional-calculus statement, not a compression which
pretends that the undressed parent subspace is invariant during a cycle.

For `z \in [0,pi^2/16]`, define

\[
 g(z)=\arccos\!\left({\cos\sqrt z\over\sqrt2}\right)-{\pi\over4}.
 \tag{FI21}
\]

Writing `x=sqrt(z)`, direct differentiation gives

\[
 g'(z)={\sin x\over2x\sqrt{1+\sin^2x}},
 \qquad
 g''(z)={x\cos x-\sin x(1+\sin^2x)
  \over4x^3(1+\sin^2x)^{3/2}}.                     \tag{FI22}
\]

The continuous limits at zero are `g'(0)=1/2` and `g''(0)=-1/3`.  For
`0 \le x \le pi/4`,

\[
 0\le\sin x-x\cos x
 =\int_0^x u\sin u\,du\le{x^3\over3},
 \qquad 0\le\sin^3x\le x^3.                       \tag{FI23}
\]

Equations (FI22)--(FI23) imply

\[
 -{1\over3}\le g''(z)\le0,
 \qquad
 {z\over2}-{z^2\over6}\le g(z)\le{z\over2}.      \tag{FI24}
\]

Spectral calculus therefore gives the exact operator inequality

\[
 \boxed{
 0\preceq H_P^F+{\hbar\eta^2\over2T_F}K_N
 \preceq {\hbar\eta^4\over6T_F}K_N^2.}             \tag{FI25}
\]

Using (FI02), the leading controlled term is

\[
 -{\hbar\eta^2\over2T_F}(4I+A_N).                 \tag{FI26}
\]

After removal of its common scalar, the q4 sibling hopping coefficient is
strictly positive in magnitude,

\[
 t_{\rm sib}^{F}={\hbar\eta^2\over2T_F}>0,          \tag{FI27}
\]

and all higher common-child corrections are retained in the exact function
(FI20) and bounded by (FI25).

### Theorem `PFCD-3` -- exact scheduled common-child branch

On the same qualified dual-flip-free slice, the finite programmed F3 cycle
produces an exact isolated dressed-parent operator which is a function of the
already earned q4 kernel
`K_N=B_N^dagger B_N=4I+A_N`.  Its first nonconstant term is the same sibling
kernel used by the static `CCMAC` Schur route, with the nonasymptotic remainder
bound (FI25).  Thus the *physical separation and virtual-return role* of the
FD detuning is constructible from existing F3 interactions and an owned
schedule even though the specific static FD05 Hamiltonian is not.

## 6. Type, energy, and port ownership

No factor or exchange is left untyped in the ideal finite construction:

| role | existing owner | exact use here |
|---|---|---|
| q4 parent/child addresses and edge list | BQ4 plus supplied `FPSS` program | labels the physical `V_0--V_1` sites and saturated append edges |
| physical sites and carrier | F3 qutrit `psi_v` factors | one content-blind carrier on `P \oplus C`; blank guards and `G` |
| support and incidence | FPMH/PESC `K_e,n_e` factors | `K` passively retained; old `n` saturated and new `n` blank under the explicitly dual-flip-free two-pulse slice |
| hopping energy | BS09 `-t n_eT_e^psi` | (FI05), counted only during the carrier pulse |
| detuning energy | BS09 `epsilon_psi q_v^psi` | (FI08), counted once during the next-slab pulse |
| slab choice and pulse order | fixed orthogonal `FPSS` controller program | selects or continuously cancels admitted terms on this exact slice; no coherent graph program |
| duration and phase reference | F3 controller/clock port | owns `tau_H`, `tau_D`, and `T_F` calibration |
| switching work, heat, and recoil | inherited work/support/reservoir ports | retains every exchange; no zero-work assertion |
| source, boundary, invalid/failure, quarantine, reset, and references | inherited `FPSS`/F3 complete port census | retained throughout; blank `G`, guards, and nonedges remain explicit |

The existing symbolic port slot is not used as a hidden staggered potential.
Actual port matrices, finite ramp errors, timing tolerance, work values, and
apparatus calibration remain supplied physical data.  The theorem proves the
ideal unitary conditional on that already declared program, exactly as the
finite support compiler does; it does not call logical scheduling free.

## 7. Exact boundary and minimal lawful antecedent

The current result separates three statements which must not be conflated:

1. **Static source-off bulk:** still no positive child-only term.  `CLDNG-3`
   remains exact.
2. **Programmed admitted-generator realization:** with `epsilon_psi>0`, one
   blank next layer, and the owned repeatable two-pulse schedule, (FI10)--(FI27)
   give a uniform positive Floquet separation and exact common-child branch
   without adding an interaction.
3. **Autonomous phase:** not derived.  A natural system which realizes this
   cycle without an externally supplied program would need a concrete
   autonomous clock/controller/work dilation and a stability theorem.

Thus the minimal lawful antecedent is narrower than the staggered interaction
(FF17):

**`PROGRAMMED-NEXT-SLAB-DETUNING`.**  Choose the already admitted
`epsilon_psi>0`, `t!=0` sector; retain one blank next F3 layer; and supply a
repeatable, exactly isolated hop/next-slab-onsite schedule in which both raw
and PESC `K`-gated incidence flips are zero in the joint generator throughout
both pulses, with its existing complete controller, clock, work, support,
boundary, failure, quarantine, reset, and reference ownership.

If `epsilon_psi=0`, if `t=0`, if the next layer is not blank, or if the
controller cannot isolate and repeat the two admitted generators, this
construction does not apply.  In that case a concrete owned stagger such as
(FF17), or another separately proved physical antecedent, remains necessary.

The result closes the finite programmed detuning obstruction only.  It does
not derive the massless collective coordinate, `chi` or `kappa`, the
refinement/time scale binding, pair-field dynamics, common probes, tensor
constraints, universal stress coupling, RGRL-B, gravity, or `G`.
