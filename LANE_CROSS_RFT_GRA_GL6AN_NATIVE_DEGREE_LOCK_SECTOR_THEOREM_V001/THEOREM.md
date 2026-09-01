# Native F3 degree-lock sector theorem and finite-h linear-Ward no-go

**Short name:** `GL6AN V001`  
**Date:** 2026-08-31  
**Status:** author frozen after hostile-audit scope repairs and exact replay;
fresh independent hostile audit required before promotion  
**Claim class:** exact algebraic result on an explicitly admissible inherited
homogeneous F3/A3 parameter subfamily; exact local pair-sector reduction;
exact linear-charge no-go; controlled finite-order strong-lock matrix element

**Not claimed:** that the GL6AL finite point `(-6,1)` is critical or belongs
to any named phase; a nonlinear conservation law; a finite-`h` Ward identity;
a pole, gaplessness, hydrodynamics, physical momentum, length, or speed; a
Lorentz/common physical cone; Ricci/Einstein form; gravity; or `G`.

## 1. Native parent and question

Use the exact infinite authenticated incidence of `GL6AK`.  Its active links
are

\[
 e=(x,a)\in\mathbb L=A_3\times\{1,2,3,4\},
 \qquad A_3=\{x\in\mathbb Z^4:\mathbf1^Tx=0\}.              \tag{AN01}
\]

Regard each active link as joining an original parent node `P_x` to the
literal shared-child node `C_{x+e_a}`.  Every original node has four incident
active links.  Two active links interact exactly when they have one original
node in common; this is precisely the within-parent/shared-child graph of
`GL6AK` equation (AK04).  Each pair has one owner and every active link has
six pair neighbors.

After removal of the inherited common scalar and with `h>0`, the selected
all-formed Hamiltonian is

\[
 H=-h\sum_eX_e+\varepsilon_\star\sum_en_e
      +2U_d\sum_{\{e,f\}\in E_\infty}n_en_f+C.              \tag{AN02}
\]

GL6AL observed an equal-time redistribution toward its local pair `E` sector
at `r_epsilon=-6r_U`, where

\[
 r_\varepsilon={\varepsilon_\star\over h},\qquad
 r_U={U_d\over h}.                                          \tag{AN03}
\]

The onsite coefficient is inherited rather than freely tunable:

\[
 \varepsilon_\star=\Delta+2U_d(1-2d_\star),
 \qquad h,\Delta>0,\quad U_d\ge0.                            \tag{AN03a}
\]

Consequently, for `U_d>0`, the proposed lock line belongs to the strict
inherited domain exactly when

\[
 \varepsilon_\star=-6U_d
 \quad\Longleftrightarrow\quad
 \Delta=4U_d(d_\star-2),                                    \tag{AN03b}
\]

which requires `d_star>2`.  The domain is nonempty; for example

\[
 d_\star=3,\qquad \Delta=4U_d>0                             \tag{AN03c}
\]

is an exact inherited witness.  At the commonly explored `d_star=2`, the
same algebraic line would require `Delta=0` and is only a boundary extension,
not a point in the strict inherited `Delta>0` domain.  Every positive
strong-lock conclusion below is asserted on the admissible `d_star>2`
subfamily (or explicitly as algebra at the boundary), never by treating
`epsilon_star` as an independent coupling.

That mutable finite-size observation is only the question generator.  None
of its sampled values enter the proof.

## 2. Exact degree-lock identity

For each original incidence node `v`, put

\[
 k_v=\sum_{e\ni v}n_e,\qquad q_v=k_v-2.                    \tag{AN04}
\]

At one degree-four node, with occupation number `k`,

\[
 2{\binom{k}{2}}-3k=k^2-4k=(k-2)^2-4.                     \tag{AN05}
\]

Every pair in (AN02) has one original-node owner and every active link has
two original endpoints.  Therefore, on the ray

\[
 \boxed{\varepsilon_\star=-6U_d},                           \tag{AN06}
\]

the interaction identity is exactly

\[
 2U_d\sum_{\{e,f\}}n_en_f-6U_d\sum_en_e
   =U_d\sum_v\big[(k_v-2)^2-4\big].                        \tag{AN07}
\]

Thus, modulo the extensive common scalar,

\[
 \boxed{H_{\rm lock}=-h\sum_eX_e+U_d\sum_vq_v^2+C'.}       \tag{AN08}
\]

Equivalently, with `Z_e=1-2n_e`, all one-link `Z` terms cancel and

\[
 H_{\rm lock}=-h\sum_eX_e+{U_d\over2}
       \sum_{\{e,f\}\in E_\infty}Z_eZ_f+C''.               \tag{AN09}
\]

Equations (AN07)--(AN09) are native algebraic identities, not an imported
model analogy.  For `U_d>0`, the diagonal minimum is the local degree-two
sector `q_v=0`.  For `U_d<0`, it is not a positive lock and the strong-lock
conclusions below do not apply.

## 3. Exact finite-h conservation result

The transverse term does not preserve the local degree.  Since
`[X_e,n_e]=iY_e`,

\[
 [H_{\rm lock},q_v]=-ih\sum_{e\ni v}Y_e.                   \tag{AN10}
\]

More generally, for a linear degree candidate

\[
 Q_c=\sum_vc_vq_v,
\]

Pauli-string independence gives

\[
 [H_{\rm lock},Q_c]
   =-ih\sum_{e=\{v,w\}}(c_v+c_w)Y_e.                       \tag{AN11}
\]

It vanishes only if `c_v+c_w=0` on every original incidence edge.  On the
connected bipartite incidence graph this fixes `c` up to the parent/child
sign.  A finite-support solution on the infinite graph is therefore zero.
On a connected balanced periodic quotient, the sole formal solution gives

\[
 \sum_{v\in P}q_v-\sum_{v\in C}q_v=0                       \tag{AN12}
\]

identically, because every active link has one parent and one child endpoint.
Hence there is no nontrivial conserved local or global *linear degree
charge* at finite `h`.  This rules out a linear degree-continuity/Ward
argument from (AN08) alone.  It does not prove the absence of every possible
nonlinear or emergent conservation law.

There is an exact discrete symmetry of the infinite homogeneous degree-four
interaction.  Let `mathcal C` be the global on-site product automorphism
defined locally by `mathcal C(Z_e)=-Z_e` and
`mathcal C(X_e)=X_e`.  Since every original incidence node has all four
links, on every local degree operator and local Hamiltonian term it gives

\[
 \mathcal C(q_v)=-q_v,
 \qquad \mathcal C(q_v^2)=q_v^2,
 \qquad \mathcal C(X_e)=X_e.                                \tag{AN13}
\]

Thus the infinite homogeneous interaction is invariant term by term.  This
is a global two-valued symmetry, not a continuous Ward identity.
The same statement holds on a degree-four periodic quotient, and for a fully
collared interior term when all four incident links are flipped.  It does
**not** hold for a generic finite open product over `F`: at a boundary node,
unflipped outside links remain and `k_v` maps to neither `4-k_v` nor `-q_v`.

## 4. Exact local `A1/E/T2` filter in the locked sector

At one original node write `Z_a=1-2n_a` and define its six pair observables

\[
 M_{ab}=Z_aZ_b,\qquad 1\le a<b\le4.                         \tag{AN14}
\]

The lock `k=2` is equivalent to `sum_a Z_a=0`.  It follows exactly that

\[
 \sum_{a<b}M_{ab}=-2,
 \qquad
 \sum_{b\ne a}M_{ab}=-1\quad(a=1,\ldots,4).                \tag{AN15}
\]

Let `R` be the unsigned `4 by 6` vertex/pair incidence matrix of `K_4`.
Then variations of a locked pair vector obey

\[
 R\,\delta M=0.                                             \tag{AN16}
\]

Exact row reduction gives `rank(R)=4` and `dim ker(R)=2`.  In pair order
`(12,13,14,23,24,34)`, a basis is

\[
 (1,0,-1,-1,0,1),\qquad(0,1,-1,-1,1,0).                   \tag{AN17}
\]

Both vectors have line-graph eigenvalue `-2` and opposite-pair eigenvalue
`+1`; this identifies the kernel exactly as the `E` summand in the native
`A1 \oplus E \oplus T2` six-pair decomposition used by GL6AL.  Thus

\[
 \boxed{\text{within }q_v=0:\quad A1\text{ is fixed},\quad
 E\text{ can fluctuate},\quad T2\text{ is zero}.}           \tag{AN18}
\]

For the exactly uniform ensemble on the six local `k=2` states, the connected
covariance orbits are

\[
 (c_d,c_a,c_o)=\left({8\over9},-{4\over9},{8\over9}\right), \tag{AN19}
\]

and therefore

\[
 \boxed{(C_{A1},C_E,C_{T2})=\left(0,{8\over3},0\right).}    \tag{AN20}
\]

This proves the algebraic direction of the GL6AL redistribution.  It does
not say that the finite-`h/U_d=1` comparator is in the strict lock, that its
state is locally uniform, or that `E` has slow dynamics.

## 5. Incidence flat directions and the limited quadratic soft statement

Let `B` be the original node/link incidence.  Multiplying every child row by
`-1` makes it an oriented bipartite incidence without changing `B^*B`.  For
the continuous extension of the diagonal lock, a variation from any locked
configuration obeys the exact increment identity

\[
 H_0(n+\delta n)-H_0(n)=U_d\|B\,\delta n\|^2.               \tag{AN21}
\]

On a finite connected quotient,

\[
 \operatorname{rank}B=|V|-1,
 \qquad \dim\ker B=|E|-|V|+1.                               \tag{AN22}
\]

For a translation character, choose a harmless row/column phase gauge so
that the one-cell symbol is

\[
 B(\chi)=
 \begin{pmatrix}
 1&1&1&1\\
 z_1&z_2&z_3&z_4
 \end{pmatrix},
 \qquad |z_a|=1.                                            \tag{AN23}
\]

Writing `s=sum_a z_a`, the two nonzero eigenvalues of `B(chi)^*B(chi)` are

\[
 \lambda_\pm=4\pm|s|.                                      \tag{AN24}
\]

There are therefore two exact constraint-flat directions at a generic
character.  At the trivial character the rank drops to one and there are
three.  With `z_a=exp(i theta_a)` and the common phase fixed by
`sum_a theta_a=0`,

\[
 4-|s|={1\over2}\sum_a\theta_a^2+O(\theta^4).               \tag{AN25}
\]

This is an exact algebraic quadratic softening of one *constraint-Gram
eigenvalue*, equivalently one squared singular value of `B(chi)`, near the
trivial character.  The corresponding singular value itself vanishes
linearly.  The character has not been calibrated as
physical momentum; (AN25) is not a real-time dispersion, a pole, a physical
propagator, or a cone.  The continuous flat directions are not by themselves
quantum zero modes of (AN08).

## 6. Controlled native collective motion on the positive strong-lock ray

Now take `U_d>0`.  Use the explicit period-four quotient of the `A3`
incidence constructed by the verifier: parent coordinates are reduced modulo
four in three independent `A3` directions.  Call it `Q_4=(V_Q,E_Q)`.  The
verifier proves it is simple, degree four, and has no two- or four-cycles; it
has `M=|E_Q|=256` active links.  Define

\[
 H_{0,\mathcal Q}=U_d\sum_{v\in V_{\mathcal Q}}q_v^2,
 \qquad V_{\mathcal Q}=-h\sum_{e\in E_{\mathcal Q}}X_e,
 \qquad h/U_d\longrightarrow0.                              \tag{AN26}
\]

Let `P_Q` be the ordinary finite-dimensional projector onto `q_v=0` for every
`v in V_Q`.  No global projector in the infinite quasi-local algebra is
asserted.  This particular girth-at-least-six quotient is an algebraic
regulator of the homogeneous incidence, not an authenticated finite open
mission.  Flipping one active
link creates one unit defect at each endpoint, so its exact cost is `2U_d`.
Consequently, in canonical Hermitian Kato/Schrieffer--Wolff perturbation
theory,

\[
 P_{\mathcal Q}V_{\mathcal Q}P_{\mathcal Q}=0,
 \qquad H_{{\rm eff},\mathcal Q}^{(2)}
 =-{Mh^2\over2U_d}P_{\mathcal Q}.                           \tag{AN27}
\]

Equivalently, the second-order scalar shift is `-h^2/(2U_d)` per link.  No
infinite `|mathbb L|` energy is written.  A finite collared open calculation
must separately declare its flippable set and boundary convention before
using an extensive count.

The fourth-order term is also exactly scalar on this finite quotient.  Put
`Q_Q=1-P_Q` and `R_Q=-Q_Q/H_(0,Q)`.  With
`P_Q V_Q P_Q=0`, the canonical expression is

\[
 H_{{\rm eff},\mathcal Q}^{(4)}
 =P_{\mathcal Q}V_{\mathcal Q}R_{\mathcal Q}V_{\mathcal Q}
  R_{\mathcal Q}V_{\mathcal Q}R_{\mathcal Q}V_{\mathcal Q}P_{\mathcal Q}
 -{1\over2}\{P_{\mathcal Q}V_{\mathcal Q}R_{\mathcal Q}^2
  V_{\mathcal Q}P_{\mathcal Q},
  P_{\mathcal Q}V_{\mathcal Q}R_{\mathcal Q}V_{\mathcal Q}P_{\mathcal Q}\}.
                                                                    \tag{AN27a}
\]

For two distinct flipped links, the intermediate two-flip energy in units of
`U_d` is `2` for adjacent opposite occupations, `6` for adjacent equal
occupations, and `4` for disjoint links.  A locked degree-two configuration
on the degree-four incidence has respectively `2M`, `M`, and
`binom(M,2)-3M` such unordered pairs.  Exact Q-only path counting gives

\[
 P_{\mathcal Q}V_{\mathcal Q}R_{\mathcal Q}V_{\mathcal Q}R_{\mathcal Q}
 V_{\mathcal Q}R_{\mathcal Q}V_{\mathcal Q}P_{\mathcal Q}
 =-{(3M^2+7M)h^4\over24U_d^3}P_{\mathcal Q},
\]

while the folded term is `+M^2h^4/(8U_d^3)P_Q`.  Hence

\[
 \boxed{H_{{\rm eff},\mathcal Q}^{(4)}
 =-{7M\over24}{h^4\over U_d^3}P_{\mathcal Q}.}              \tag{AN27b}
\]

There is no off-diagonal term through fourth order on this declared `Q_4`
because its incidence girth is at least six.  This statement is not assigned
to arbitrary smaller periodic quotients, which can acquire wrapped
four-cycles.  Thus configuration dependence is not left open below order six
for `Q_4`; sixth and higher diagonal/loop terms remain unclassified.

Every finite-order perturbative word changes only finitely many links.  The
nonempty finite symmetric difference of two degree-two configurations is
therefore a finite even-degree subgraph of the original parent/child
incidence graph and decomposes into finite cycles.  The same is true directly
on the declared `Q_4`.  (An arbitrary infinite even subgraph could contain a
bi-infinite path and is not used here.)  The infinite incidence graph has no
parallel-edge two-cycle.  A four-cycle would require
`e_a-e_b=e_c-e_d` for two distinct ordered port pairs, which is impossible
because the `+1` and `-1` coordinates determine `(a,b)` uniquely.  Hence any
nonempty difference contains at least a six-cycle.

One native chordless six-cycle is, with `p=e_1-e_2` and `q=e_3-e_2`,

\[
 (0,1),(p,2),(p,3),(q,1),(q,2),(0,3).                      \tag{AN28}
\]

The exact verifier constructs on `Q_4` a locked configuration in which these
six links alternate occupied/unoccupied.  Its
periodic lift supplies an infinite locked background, and every finite collar
around the hexagon embeds into the infinite authenticated incidence by AK07.
Toggling all six links produces a distinct finite-quotient locked
configuration.  No proper subset of the six toggles is locked.

For an ordering `pi` of the six flips, let `S_{pi,j}` be its first `j` links
and let

\[
 E(S)=\sum_v q_v(S)^2.                                      \tag{AN29}
\]

Exact enumeration of all `6!=720` orderings gives

\[
 \sum_{\pi\in S_6}\prod_{j=1}^{5}{-1\over E(S_{\pi,j})}
   =-{63\over8}.                                             \tag{AN30}
\]

Because no lower-order process connects these two locked states, folded or
normalization terms cannot alter their first nonzero off-diagonal entry.
Thus the leading native collective matrix element is

\[
 \boxed{
 \langle f|H_{{\rm eff},\mathcal Q}^{(6)}|i\rangle
   =-{63\over8}{h^6\over U_d^5}.}                            \tag{AN31}
\]

Equation (AN31) is a finite linked matrix element.  Its path denominators use
only the displayed hexagon collar, so the same linked coefficient embeds in
the infinite incidence; no infinite locked-sector projector is used.  This
is a genuine multi-link collective operation enabled on the selected
all-formed lineage branch and generated by the same F3 Hamiltonian, with no
new operator inserted.  The six active-link occupations are not promoted to
six independently authenticated records.  The result is not evidence by
itself that the effective locked theory is gapless: sixth/higher diagonal
terms and the full family of connected order-six loop terms have not yet been
derived.

## 7. The theorem/no-go conclusion

The degree-lock clue does uncover native structure.  The exact ray supplies
a constrained sector, an exact local `E` filter, extensive continuous
incidence-flat directions, a quadratically soft constraint-Gram eigenvalue,
and a nonzero finite-order collective motion within the locked sector.

It simultaneously closes the tempting short argument that the square itself
is already a finite-`h` conserved density: it is not.  Neither the static
flatness nor the `E` filtering constitutes a dynamical infrared law.

The decisive next calculation is therefore bounded and native: derive the
complete connected degenerate effective Hamiltonian through sixth order,
including all sixth-order diagonal terms and every native minimal-loop move,
and then test whether its stationary `A1/E/T2` response remains soft under
increasing authenticated regions.  A finite gap would refute this ray as the
sought infrared route.  A controlled vanishing threshold would advance it,
but would still require the separate physical-cone, complete-stress, and
response-form gates.

## 8. Strict ceiling

Nothing in this theorem identifies a translation character with calibrated
physical momentum, proves a stationary pole or common cone, derives a native
stress operator, obtains Ricci/Einstein response, establishes gravity, or
calculates `G`.
