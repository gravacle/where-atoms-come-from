# Reciprocal carrier-to-incidence back-propagation screen

**Lane ID:** `GRA-DQ-F3-RBPS-V001`

**Date:** 2026-08-27

**Claim class:** exact negative real-time response theorem on the sealed
GRA-DJ/GRA-DN hold; exact static conditional-link reciprocity; exact
finite-volume restored-BS06 one-link Rabi observable

**Disposition:**
`ZERO_RETARDED_CARRIER_TO_INCIDENCE_AND_K_RESPONSE_IN_THE_SEALED_HOLD__STATIC_GIBBS_RECIPROCITY_ONLY__RESTORED_BS06_FINITE_LINK_RESPONSE_IS_THE_NEXT_LAWFUL_TEST__NO_FUTURE_RECORD_FORMATION_OR_GRAVITY_CLAIM`

**Not claimed:** an autonomous post-GRA-DN feedback protocol; a change in
future FPMH formation probability; a common-background comparison of the
GRA-DN ordered and disordered arms; a nonzero fixed-time local thermodynamic
response; a metric, stress tensor, universal coupling, Newton `G`, or gravity

## 1. Exact sealed slice and question

On a supplied square support let

\[
 n_e={1-Z_e\over2},\qquad
 C_e:=(J^\psi_e)^2
 =q_u+q_v-2q_uq_v={1-s_us_v\over2}.                  \tag{DQ01}
\]

Here `C_e` is the fixed-content hard-core carrier-current square and is a
projector with eigenvalues zero and one.  GRA-DJ takes

\[
 h=t=\epsilon_\psi=\eta=0                             \tag{DQ02}
\]

and anneals incidence only in a static Gibbs trace.  After the common reader,
GRA-DN takes the same `h=t=0` carrier hold, freezes the decoded incidence word,
switches off every source, imprint, route, and reader pulse, and retains only
the unchanged current-square Hamiltonian.  No FPMH writer is replayed.
The carrier trace in these phase slices is grand canonical with a zero
one-body number ledger; no fixed-total-carrier-number theorem is imported.

For homogeneous fixed storage `r in {0,1}`, the GRA-DJ hold is

\[
 H_{\rm hold}=U_d\sum_v[d_v(n)-4]^2
 +\sum_e n_e\,[a_r+\lambda_J C_e],
 \qquad a_r=\Delta-2\lambda_Rr.                      \tag{DQ03}
\]

Here `r` is the BS storage-occupation condition used by GRA-DJ.  It is not identified with the FPMH `K` factor.  Occupation alone is not an authenticated
`REC` verdict.

The GRA-DN post-reader hold is the corresponding quenched-incidence block
`H_psi=lambda_J sum_e n_e C_e`.  Both Hamiltonians act as the identity on the
retained FPMH `K` registers.

The question is deliberately dynamical: after a record-supported carrier
collective state has formed, can this unchanged hold itself change incidence
or the support of an authenticated `K` record?

## 2. Exact obstruction: the retarded response is zero

Every term in (DQ03) is diagonal in the incidence occupations.  At `t=0` it is
also diagonal in the carrier occupations.  Therefore, for all edges `e,f` and
all FPMH occupation projectors `P_i^K`,

\[
 [H_{\rm hold},n_e]=[H_{\rm hold},C_f]
 =[H_{\rm hold},P_i^K]=0.                             \tag{DQ04}
\]

In particular,

\[
 n_e(\tau)=n_e,\qquad P_i^K(\tau)=P_i^K.              \tag{DQ05}
\]

Let `B_psi` be any carrier-only perturbation observable, not necessarily
`C_f`.  Since carrier and incidence/`K` act on distinct tensor factors,

\[
 \boxed{
 \chi^R_{n_e,B_\psi}(\tau)
 ={i\over\hbar}\theta(\tau)
 \langle[n_e(\tau),B_\psi]\rangle=0,
 \qquad
 \chi^R_{P_i^K,B_\psi}(\tau)=0.}                    \tag{DQ06}
\]

The conclusion is state-independent, volume-independent, and exact at every
time.  Hence every thermodynamic subsequential limit is also identically
zero.  BS11 supplies a joint energy, but at `h=0` it has no operator that can
change `n`; it acts as the identity on `K` and cannot change a formation
verdict.

## 3. What survives: static equilibrium conditional reciprocity

Fix every incidence except one edge `e=uv`.  Let `d_u,d_v in {0,1,2,3}` be the
endpoint degrees excluding `e`.  Adding that edge costs

\[
\begin{split}
 \Delta_eE(c)
 &=a_r+\lambda_Jc\\
 &\quad+U_d[(d_u+1-4)^2-(d_u-4)^2
            +(d_v+1-4)^2-(d_v-4)^2]\\
 &=A_e+\lambda_Jc,\\
 A_e&:=a_r+2U_d(d_u+d_v-7),\qquad c\in\{0,1\}.       \tag{DQ07}
\end{split}
\]

Consequently the exact GRA-DJ one-edge Gibbs conditional is

\[
 \boxed{
 \Pr(n_e=1\mid n_{f\ne e},s,r)
 ={1\over1+\exp[\beta(A_e+\lambda_JC_e)]}.}          \tag{DQ08}
\]

For `lambda_J>0`, an aligned endpoint pair (`C_e=0`) is statically more
favorable to incidence than an anti-aligned pair (`C_e=1`) on the same
background.  If a prescribed carrier state `sigma` has

\[
 p_\sigma=\operatorname{Tr}(\sigma C_e)
 ={1-g_\sigma\over2},\qquad
 g_\sigma=\langle s_us_v\rangle_\sigma,              \tag{DQ09}
\]

then its mean link-addition energy is

\[
 \overline{\Delta_eE}_\sigma=A_e+\lambda_Jp_\sigma,
 \qquad
 \overline{\Delta_eE}_{D}-\overline{\Delta_eE}_{O}
 ={\lambda_J\over2}(g_O-g_D).                        \tag{DQ10}
\]

The labels "ordered" and "disordered" do not determine this sign by
themselves.  The local nearest-neighbor correlator, not magnetization alone,
fixes `p_sigma`; a symmetric mixture of plus and minus ordered phases can have
zero magnetization and the same `C_e` statistics as either ordered phase.

For a finite equilibrium partition function with commuting sources
`xi n_e+zeta C_f`,

\[
 {\partial^2F\over\partial\xi\,\partial\zeta}
 =-\beta\,\operatorname{Cov}(n_e,C_f).                \tag{DQ11}
\]

Equation (DQ11) is symmetric static susceptibility of one Gibbs endpoint.
It supplies neither a retarded time arrow nor a transition law.  In GRA-DN,
incidence is quenched after the reader, so (DQ10) changes the energy assigned
to a fixed word and cannot change that word's probability dynamically.

## 4. No future FPMH formation feedback without replay

The sealed hold contains no `K` operator.  Its tensor factorization is

\[
 H_{\rm hold}=I_K\otimes H_{n\psi}.                   \tag{DQ12}
\]

The inherited controlled-link reader, when used, preserves `K`; in the hold
it is off.  The FPMH writer that creates the qualified formation/sham
distinction is also off and is never replayed.  It follows from (DQ05)--(DQ06)
that the carrier state cannot change `Pr(K_i=1)`, a final qualified query law,
or the formation-versus-sham total-variation distance.

Turning on BS10/FPMH source or writer couplings later would define another
mission with new ports and histories.  No such mission is sealed here, and no
claim about its formation probability follows from the static joint energy.

## 5. Next lawful observable: restored-BS06 conditional-link Rabi law

There is one exact finite calculation that uses no new interaction.  Restore
the already-declared BS06 flip `-h_N X_e`, set `t=0`, and use the same
conditional one-link protocol as the F3 parent: freeze storage, all other
incidences, endpoint degree contributions, and all source/writer ports.  This
is a conditional finite-link calculation, not autonomous full-sheet dynamics.

On a carrier-current sector `C_e=c`, the Hamiltonian is

\[
 H_e(c)=-h_NX_e+F_c n_e,
 \qquad F_c=A_e+\lambda_Jc,
 \qquad h_N={\Omega\over N}.                           \tag{DQ13}
\]

For an initially blank link `n_e=0`, exact exponentiation gives

\[
 \boxed{
 R(F_c,\tau):=\Pr_c(0\to1;\tau)
 ={4h_N^2\over F_c^2+4h_N^2}
 \sin^2\!\left({\tau\over2\hbar}
 \sqrt{F_c^2+4h_N^2}\right).}                        \tag{DQ14}
\]

Because `C_e` is a projector conserved at `t=0`, an arbitrary prescribed
endpoint carrier density operator `sigma` gives

\[
\boxed{
 P_\sigma(\tau)
 =(1-p_\sigma)R(A_e,\tau)
 +p_\sigma R(A_e+\lambda_J,\tau).}                   \tag{DQ15}
\]

This finite conditional formula likewise requires no fixed-number ensemble;
its only carrier input is the complete prescribed endpoint density operator.

Thus two carrier preparations on the **same** fixed storage/incidence/degree
background obey

\[
 P_O-P_D=(p_O-p_D)
 [R(A_e+\lambda_J,\tau)-R(A_e,\tau)].                \tag{DQ16}
\]

This is generically nonzero at finite `N`.  The route contains an explicit
incidence actuator, namely the restored BS06 `X_e`; BS11 only changes its
detuning.  BS11 alone still cannot flip the link.

For a `Z`-eigenstate link, the carrier dependence first enters at fourth order:

\[
 R(F,\tau)={h_N^2\tau^2\over\hbar^2}
 -{h_N^2(F^2+4h_N^2)\tau^4\over12\hbar^4}
 +O(\tau^6).                                         \tag{DQ17}
\]

At `h_N=0`, (DQ14)--(DQ17) reduce exactly to the obstruction in Section 2.

## 6. Why this is not yet a GRA-DN reciprocal theorem

GRA-DN forms its ordered carrier state on the full `KEEP` incidence sheet and
its uniform product carrier state on the blank `REBIND` sheet.  Those arms do
not hold `n`, the degree background, or `A_e` fixed.  Feeding their carrier
states into (DQ16) would therefore compare different backgrounds and would
reproduce the already-known joint endpoint bias rather than isolate
carrier-to-incidence response.

A matched test must prepare two carrier inputs with different `p_sigma` on the
same fixed `n`, `K`, storage, and degree background, then apply the same BS06
pulse.  Such reset, transport, or re-preparation ports and their work/heat
ledger are not supplied by GRA-DN.  The common pulse is exact, but equality of
carrier-preparation work or of complete laboratory work distributions is not
claimed.

Finally, `h_N=Omega/N` makes the fixed-time local response vanish.  From
`sin^2 x <= x^2`,

\[
 0\le R(F,\tau)\le{\Omega^2\tau^2\over N^2\hbar^2}. \tag{DQ18}
\]

No sealed collective summation, time scaling, relaxation law, or nonzero
thermodynamic response has been proved.  The exact hold result remains zero in
every volume; the restored-BS06 observable is a finite conditional next test.

## 7. Scientific ceiling

This lane proves a sharp obstruction, not backreaction: the actual GRA-DJ and
GRA-DN `h=t=0` holds have zero retarded carrier-to-incidence and carrier-to-`K`
response.  Their equilibrium carrier/incidence covariance is static only.
The smallest lawful nonzero observable is the finite conditional Rabi contrast
(DQ16), which restores an existing F3 actuator and requires a common
background not supplied by the GRA-DN phase comparison.

Nothing here proves future record formation feedback, semantic `REC`
necessity, unique microscopic causation, autonomous support, emergent metric,
physical stress response, universal stress coupling, nonlinear gravitational
closure, Newton `G`, or gravity.
