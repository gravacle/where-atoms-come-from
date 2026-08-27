# Equivariant historywise-gravity nonselection theorem

**Claim class:** exact finite-group orbit theorem; exact stabilizer criterion
for added orienting input; exact finite two-cell Poisson/metric witnesses;
conditional interface to GARH-D

**Not claimed:** physical general relativity, a derivation of a historywise
gravity map, objective actualization, a Born law, a collapse dynamics, a
physical boundary variable, relativistic covariance, or empirical evidence

## 1. Inherited GARH and MGFT objects

Fix before the realized outcome a finite physical context C of the kind used in
the Global Actual-Record History lane. It has a finite candidate completed-
history space

\[
 \mathfrak H_C
\]

and an independently defined GARH admissible core

\[
 \mathfrak A_C\subseteq\mathfrak H_C.                 \tag{1}
\]

Membership in (1) may impose the inherited dynamics, causal consistency,
conservation, bona-fide record lineage, and record--outcome identity

\[
 a_e(H)=x\quad\Longleftrightarrow\quad
 \mathsf{FirstRF}_{e,x}(H).                            \tag{2}
\]

Equation (2) is evaluated inside a candidate completed history. It does not
select that history.

Let \(\mathcal Y_C\) be a nominated raw gravity-data space. Before any orbit is
counted, quotient declared gauge, diffeomorphism, and pure-description copies
that represent the same physical record history and the same physical gravity
data:

\[
 \widehat{\mathfrak A}_C=\mathfrak A_C/\!\sim_{\rm gauge},
 \qquad
 \widehat{\mathcal Y}_C=\mathcal Y_C/\!\sim_{\rm gauge}.       \tag{2a}
\]

The physical symmetry action below is required to descend to these quotient
spaces. To avoid heavy notation, H and y henceforth denote physical equivalence
classes. Orbit cardinalities count physically distinct histories, not gauge
copies. If an alleged outcome swap is pure gauge, it does not define two
physical outcomes and is outside the fixed-point-free case.

This lane adds a nominated historywise gravity-emergence map

\[
 \mathscr E_C:\widehat{\mathfrak A}_C
 \longrightarrow\widehat{\mathcal Y}_C,
 \qquad
 \mathscr E_C(H)=Y_H.                                  \tag{3}
\]

Before a G3 reconstruction is independently earned, \(Y_H\) means generic
complete gravity-sector data: source/stress information, response kernels,
constraints, boundary data, and every required charge, dressing, radiation,
memory, and worldtube port. Only on a separately established metric domain may
one write part of it as \((\mathsf g_H,\mathsf T_H)\). Nothing here proves that
such a complete map exists or that any component is a physical spacetime metric.

Let \(\mathsf G\) be a finite group of physical symmetries of C. Use \(u\) for
a group element so it is not confused with Newton's constant or a metric. The
group acts on histories, gravity outputs, residual spaces, and any additional
physical input. The context is fixed:

\[
 uC=C.                                                  \tag{4}
\]

The quotient GARH candidate core is also symmetry covariant:

\[
 u\widehat{\mathfrak A}_C=\widehat{\mathfrak A}_C.     \tag{4a}
\]

The historywise map is equivariant when

\[
 \boxed{\mathscr E_C(uH)=u\mathscr E_C(H).}            \tag{5}
\]

Equation (5) says that physically transforming the completed history transforms
its complete gravity-sector output. It does not say that the two branch outputs
are equal or that either output already contains an earned spacetime metric.

## 2. Gravity-extended admissibility

Let \(J_C(H,y)\in\{0,1\}\) be a Boolean record/gravity admissibility predicate,
and let

\[
 R_C(H,y)\in\mathcal V_C                              \tag{6}
\]

be a possibly vector-valued field-equation, boundary, conservation, or feedback
residual. The zero vector is fixed by the group. Require

\[
 J_C(uH,uy)=J_C(H,y),                                  \tag{7}
\]

\[
 R_C(uH,uy)=uR_C(H,y).                                 \tag{8}
\]

Define the endogenous gravity-consistent set

\[
 \mathfrak A_C^{\rm grav}
 =\{H\in\widehat{\mathfrak A}_C:
 J_C(H,\mathscr E_C(H))=1,
 R_C(H,\mathscr E_C(H))=0\}.                          \tag{9}
\]

The definitions and physical referents in (3), (6), and (9) must be fixed
before scoring a realized label. Otherwise the construction is a post-outcome
recoding, not a GARH-D mechanism.

## 3. Theorem HGA1 -- orbit closure and singleton obstruction

Under (4)--(9):

1. \(\mathfrak A_C^{\rm grav}\) is a union of complete \(\mathsf G\)-orbits;
2. if \(H\in\mathfrak A_C^{\rm grav}\), then
   \(\operatorname{Orb}(H)\subseteq\mathfrak A_C^{\rm grav}\);
3. if \(\mathfrak A_C^{\rm grav}=\{H_*\}\), then
   \(uH_*=H_*\) for every \(u\in\mathsf G\); and
4. no invariant gravity-consistency law can admit exactly one member of a
   fixed-point-free definite-outcome orbit.

Equivalently,

\[
 \boxed{
 |\operatorname{Orb}(H)|>1\ \land\ H\in\mathfrak A_C^{\rm grav}
 \Longrightarrow
 |\mathfrak A_C^{\rm grav}|\ge
 |\operatorname{Orb}(H)|.}                            \tag{10}
\]

### Proof

Take \(H\in\mathfrak A_C^{\rm grav}\) and \(u\in\mathsf G\). The inherited
quotient GARH core is physically symmetry covariant, so
\(uH\in\widehat{\mathfrak A}_C\).
Equations (5) and (7) give

\[
 J_C(uH,\mathscr E_C(uH))
 =J_C(uH,u\mathscr E_C(H))
 =J_C(H,\mathscr E_C(H))=1.
\]

Equations (5) and (8) similarly give

\[
 R_C(uH,\mathscr E_C(uH))
 =uR_C(H,\mathscr E_C(H))=0.
\]

Thus \(uH\in\mathfrak A_C^{\rm grav}\), proving orbit closure and (10). A
singleton orbit has stabilizer equal to the whole group, proving item 3. QED.

### Exact interpretation

The theorem does not require \(\mathscr E_C(H)=\mathscr E_C(uH)\). The two
histories may produce visibly different branch-conditioned gravity data and,
on an independently earned metric domain, different metrics and stress
tensors. Equivariance makes those outputs transform around the same orbit, so
endogenous gravity follows the branch but does not orient the branch.

A unique symmetry-fixed history is not excluded. In an outcome-swap problem it
can be a symmetric unresolved history rather than one definite outcome.

## 4. Corollaries: nonlinear feedback and invariant scores

Let a nonlinear update rule obey

\[
 \Psi_C(uH,uy)=u\Psi_C(H,y).                           \tag{11}
\]

Its historywise self-consistency solutions

\[
 \mathfrak F_C=
 \{H\in\widehat{\mathfrak A}_C:
 \Psi_C(H,\mathscr E_C(H))=H\}                        \tag{12}
\]

form a union of complete group orbits. Indeed, applying \(u\) to (12) and using
(5) and (11) proves that \(uH\) is also a solution. Nonlinearity and
backreaction do not change the symmetry conclusion.

Likewise, let \(Q_C(H,y)\in\mathbb R\) be invariant:

\[
 Q_C(uH,uy)=Q_C(H,y).                                  \tag{13}
\]

Then the historywise score

\[
 q_C(H)=Q_C(H,\mathscr E_C(H))                         \tag{14}
\]

is constant on every orbit. Its minimizers and maximizers are unions of orbits.
An invariant action, residual norm, curvature scalar, record-quality score, or
entropy functional therefore cannot uniquely choose a nonfixed member of a
symmetric orbit.

### Corollary HGA1a -- unique symmetric deterministic flow

Let \(\mathcal Z_C\) be a space of complete initial/boundary data and let

\[
 F_C:\mathcal Z_C\longrightarrow\widehat{\mathfrak A}_C       \tag{14a}
\]

be a deterministic, single-valued, equivariant history flow. If the complete
input \(z_0\) is group fixed, then its unique history is group fixed:

\[
 uz_0=z_0\ \forall u
 \quad\Longrightarrow\quad
 uF_C(z_0)=F_C(z_0)\ \forall u.                        \tag{14b}
\]

Thus a symmetric unique-flow law cannot produce one member of a fixed-point-
free outcome orbit from complete symmetric data. This does not forbid a unique
outcome when the actual initial/boundary microstate is asymmetric; that
microstate is then the additional orienting input audited by HGA3.

### Corollary HGA1b -- unique distributional fixed law

Let \(\mathcal O\subseteq\widehat{\mathfrak A}_C\) be a finite transitive
history orbit. The group acts on probability laws by pushforward. Let

\[
 \mathcal K:\mathcal P(\mathcal O)\longrightarrow
 \mathcal P(\mathcal O)                                \tag{14c}
\]

be an equivariant update on distributions:

\[
 \mathcal K(u_*\mu)=u_*\mathcal K(\mu).                \tag{14d}
\]

If \(\mu_*\) is the unique fixed distribution of \(\mathcal K\), then

\[
 u_*\mu_*=\mu_*\quad\forall u.                        \tag{14e}
\]

On a finite transitive orbit this is the uniform law. Indeed, \(u_*\mu_*\) is
another fixed law by (14d), so uniqueness gives (14e); transitivity makes all
point weights equal. This result concerns a law over candidate histories. It
does not supply a realized history. Multiple symmetry-related fixed laws,
trajectory-level noise, or asymmetric boundary data lie outside its uniqueness
premise and require their own actualization account.

## 5. Corollary HGA2 -- reduced mean-field nonselection

Assume \(\widehat{\mathcal Y}_C\) is a finite-dimensional affine space and the
group acts affinely, with its linear part acting on differences. Let
\(\mu_C^{\rm hist}\) be a symmetry-invariant probability law over histories in
\(\widehat{\mathfrak A}_C\):

\[
 \mu_C^{\rm hist}(uH)=\mu_C^{\rm hist}(H).            \tag{15}
\]

This history probability law is not the inherited MGFT spacetime volume measure
\(\mu_H^{\rm vol}\).

Define the history-averaged gravity output

\[
 \overline Y_{\mu_C^{\rm hist}}
 =\sum_{H\in\widehat{\mathfrak A}_C}
   \mu_C^{\rm hist}(H)\mathscr E_C(H).                \tag{16}
\]

Only when metric/stress coordinates have independently been earned may this be
written as \((\overline{\mathsf g}_{\mu_C^{\rm hist}},
\overline{\mathsf T}_{\mu_C^{\rm hist}})\).

Then

\[
 u\overline Y_{\mu_C^{\rm hist}}
 =\overline Y_{\mu_C^{\rm hist}}
 \quad\text{for every }u\in\mathsf G.                \tag{17}
\]

If a deterministic selector
\(S_C:\widehat{\mathcal Y}_C\to\widehat{\mathfrak A}_C\) is equivariant,
then

\[
 S_C(\overline Y_{\mu_C^{\rm hist}})
 =uS_C(\overline Y_{\mu_C^{\rm hist}})               \tag{18}
\]

for every \(u\). It cannot return a history in a fixed-point-free outcome
orbit.

### Proof

Using (5), invariance of \(\mu_C^{\rm hist}\), and relabelling the finite sum,

\[
 u\overline Y_{\mu_C^{\rm hist}}
 =\sum_H\mu_C^{\rm hist}(H)\mathscr E_C(uH)
 =\sum_{H'}\mu_C^{\rm hist}(H')\mathscr E_C(H')
 =\overline Y_{\mu_C^{\rm hist}}.
\]

Equivariance of the selector then gives (18). QED.

This finite corollary blocks the claim that a symmetry-averaged generic gravity
output by itself deterministically selects one branch. On an independently
earned metric domain the same statement applies to the corresponding averaged
metric/stress coordinates. It is not a general no-go theorem for semiclassical
gravity, stochastic gravity, symmetry-breaking fluctuations, or branch-relative
gravitational observables.

## 6. Theorem HGA3 -- exact orienting-input criterion

Let \(\Lambda\) be a finite \(\mathsf G\)-set of prospectively fixed physical
inputs. A covariant gravity/admissibility family may depend on
\(\lambda\in\Lambda\):

\[
 J_C(uH,uy;u\lambda)=J_C(H,y;\lambda),                 \tag{19}
\]

with an analogous covariance law for every residual or score. For fixed
\(\lambda\), its unique selected history, if it exists, need only be fixed by
the stabilizer

\[
 \mathsf G_\lambda=\{u:u\lambda=\lambda\}.            \tag{20}
\]

More exactly, let the input orbit and candidate history orbit be transitive:

\[
 \Lambda\simeq\mathsf G/K,
 \qquad
 \mathcal O\simeq\mathsf G/L.                         \tag{21}
\]

A \(\mathsf G\)-equivariant deterministic selector

\[
 s:\Lambda\longrightarrow\mathcal O                  \tag{22}
\]

exists iff, after a possible conjugation of L,

\[
 \boxed{K\subseteq L.}                                \tag{23}
\]

When it exists,

\[
 |\Lambda|={|\mathsf G|\over|K|}
 \ge {|\mathsf G|\over|L|}=|\mathcal O|.             \tag{24}
\]

### Proof

Choose representatives \(\lambda_0\) and \(H_0\). If \(K\subseteq L\), set

\[
 s(u\lambda_0)=uH_0.                                  \tag{25}
\]

If \(u\lambda_0=v\lambda_0\), then \(v^{-1}u\in K\subseteq L\), so
\(uH_0=vH_0\); hence (25) is well defined and equivariant. Conversely, if an
equivariant selector sends \(\lambda_0\) to \(H_0\), every element fixing
\(\lambda_0\) must fix \(H_0\), so \(K\subseteq L\). Other choices of
\(H_0\) conjugate L. Equation (24) follows from finite orbit--stabilizer. QED.

For a free two-history swap orbit, L is trivial. The orienting input must also
have trivial stabilizer. A group-fixed scalar input has K equal to the whole
group and cannot help. If the proposed gravity map or law violates equivariance,
that non-equivariance is itself a preferred orientation and must be exposed as
additional physical structure.

## 7. Exact two-cell historywise-gravity witness

Let \(\mathsf G=\mathbb Z_2=\{1,P\}\), where

\[
 P=\begin{pmatrix}0&1\\1&0\end{pmatrix},
 \qquad
 \mathcal L=\begin{pmatrix}1&-1\\-1&1\end{pmatrix}.
                                                               \tag{26}
\]

The histories \(H_+\) and \(H_-\) are exchanged by P; \(H_s\) contains the
definite record value s and nominates the corresponding source contrast below.
The two cells have separately declared actuators and reads, so P is a physical
apparatus/outcome symmetry after gauge quotient, not a pure gauge copy.
Define source/stress contrasts and gauge-fixed potentials

\[
 \mathsf T_+=\binom{1}{-1},
 \quad \mathsf T_-=P\mathsf T_+=\binom{-1}{1},         \tag{27}
\]

\[
 \phi_+=\binom{1/2}{-1/2},
 \quad \phi_-=P\phi_+=\binom{-1/2}{1/2}.              \tag{28}
\]

They obey exactly

\[
 P\mathcal LP=\mathcal L,
 \qquad
 \mathcal L\phi_\pm=\mathsf T_\pm,
 \qquad
 \mathbf 1^T\phi_\pm=0.                              \tag{29}
\]

At cell i define a rational 1+1-dimensional weak-field metric proxy

\[
 \mathsf g_\pm(i)
 =\operatorname{diag}(-(1+\phi_{\pm,i}),
                         1-\phi_{\pm,i}).              \tag{30}
\]

Because \(\phi_{\pm,i}=\pm1/2\), each cell metric has one negative and one
positive eigenvalue and determinant \(-3/4\). The group swaps the two cells, so

\[
 \mathscr E_C(PH_+)=P\mathscr E_C(H_+)
 =\mathscr E_C(H_-).                                   \tag{31}
\]

Both histories pass the same covariant Poisson, gauge, metric-signature, and
record/source-matching rules. Both also have the same invariant field action

\[
 \mathcal S(\phi,\mathsf T)
 ={1\over2}\phi^T\mathcal L\phi-\phi^T\mathsf T
 =-{1\over2}.                                         \tag{32}
\]

Thus the branch metrics are distinct and exactly self-consistent, but the
admissible set is \(\{H_+,H_-\}\), not a singleton. For the invariant uniform
law, \(\overline\phi=0\) and \(\overline{\mathsf T}=0\); the reduced mean field
is the symmetry-fixed background and also does not select.

Equation (30) is a finite gravity-like witness. It is not a discretization proof
of Einstein's equations, a continuum metric, or evidence for physical
record-conditioned gravity.

## 8. Positive swap-odd boundary witness

Add a nominated antisymmetric boundary/source vector within the finite model,

\[
 b_+=\mathsf T_+,
 \qquad b_-=\mathsf T_-,
 \qquad Pb_+=b_-.                                      \tag{33}
\]

Let its nominated algebraic coupling score be

\[
 U(H_s,b_r)=-b_r^T\mathsf T_s.                         \tag{34}
\]

The exact table is

\[
\begin{array}{c|cc}
 &H_+&H_-\\ \hline
 b_+&-2&+2\\
 b_-&+2&-2
\end{array}.                                           \tag{35}
\]

It is jointly invariant:

\[
 U(PH,Pb)=U(H,b).                                      \tag{36}
\]

For fixed \(b_+\), the unique minimum is \(H_+\); for fixed \(b_-\), it is
\(H_-\). The input and history orbits are both free, so K=L is trivial and
HGA3 is saturated. This is a positive mathematical orienting-input witness for
a GARH-D architecture: a swap-odd orientation within the finite model plus a
covariant coupling can select.

It becomes a physical GARH-D mechanism only if the field is independently real,
fixed or generated without using the realized outcome, and accompanied by its
complete dynamics and ledgers. If b is instead merely a renamed realized label,
or is assigned after the read, (35) is circular. If b is simply postulated to
take one value with Born weights and no deeper physical account, the primitive
actual draw has been moved to b rather than derived.

## 9. Operational ensemble and ancillary-composition gate

Suppose a GARH-D proposal claims standard unconditioned quantum operational
equivalence on a declared finite-dimensional input Hilbert space
\(\mathcal H_{\rm in}\), with \(d=\dim\mathcal H_{\rm in}\). Let
\(\bar\Phi\) be its unconditioned local state transformation. The inherited
actualization result requires the packet to establish all of the following on
that domain:

1. two operational preparations or remotely steerable decompositions of the
   same density operator produce the same unconditioned registered law;
2. a remote choice between such uncommunicated decompositions cannot change the
   local marginal;
3. ordinary randomized preparation obeys affine mixing,
   \[
    \bar\Phi\!\left(\sum_i p_i\rho_i\right)
    =\sum_i p_i\bar\Phi(\rho_i);                       \tag{37}
   \]
4. for every finite ancilla dimension n,
   \[
    (\bar\Phi\otimes\operatorname{id}_n)(\Omega)\ge0
    \quad\text{for every }\Omega\ge0,                 \tag{38}
   \]
   so the linear extension is completely positive. In finite input dimension it
   is equivalently sufficient to establish this at \(n=d\).

Nonlinear or stochastic conditioned histories are not excluded by (37)--(38);
the gate concerns their unconditioned operational law. If a proposal abandons
preparation/decomposition independence, remote-ensemble no-signalling, affine
mixing, or ancillary complete positivity, it must name the abandoned premise,
freeze the affected domain and scale, and state a prospective preparation-
context, signalling, or non-CP composition prediction and falsifier. A proposal
that permits only a restricted ancillary family may claim positivity on those
compositions, not complete positivity; that restriction must use the same
named, frozen departure route. Merely remaining silent about the standard
premise is not an alternative admission route.

The two-cell boundary table defines no quantum unconditioned channel and passes
none of these operational gates by itself.

## 10. Physical GARH-D admission and GARH-Q fallback

The finite theorem permits a GARH-D claim only after one packet establishes:

1. **Prospective input:** the orienting datum, domain, clocks, boundary, and
   coupling role are fixed before the realized outcome;
2. **Physical authenticity:** the datum is a physical state, field, initial or
   boundary condition, or dynamical variable, not a label string or software
   tiebreak;
3. **Complete coupling:** the actual law couples it to the candidate history and
   gravity variables, with every relevant port and no post-cut label injection;
4. **Covariance:** the joint transformation law is stated; any preferred frame,
   foliation, orientation, or boundary is exposed rather than hidden;
5. **Conservation:** outcome-resolved energy, momentum, angular momentum, charge,
   apparatus, reservoir, gravitational, and boundary exchanges close;
6. **Selection:** the complete equations yield one admissible history for the
   realized input, rather than merely assigning equal branch scores;
7. **Outcome law:** the distribution or typicality law for the orienting input
   and resulting histories is supplied and tested against Born statistics or a
   frozen deviation;
8. **Operational ensemble/composition:** preparation/decomposition
   independence, unconditioned remote-ensemble no-signalling, affine mixing,
   and completely positive ancillary composition for every finite n (or
   equivalently \(n=\dim\mathcal H_{\rm in}\)) are proved on the claimed
   standard finite-dimensional domain, or each abandoned premise or ancillary
   restriction is named with a frozen signalling, preparation-context, or
   non-CP prediction and falsifier as specified in Section 9;
9. **Record join:** the selected history satisfies the independently physical
   record predicate and any SEALED/RMR use has its own premises plus JoinCompat;
10. **Gravity scope:** the proposal states whether \(Y_H\) is generic
    gravity-sector data, a finite proxy, or contains an independently earned
    effective/physical metric, and supplies the additional G3/G4 premises
    appropriate to that claim; and
11. **Falsifier:** a prospective observation can reject the coupling, selection,
    covariance, conservation, or probability claim.

Passing only the finite orbit theorem or the two-cell example earns no physical
GARH-D admission.

If no added physical orientation or deeper law passes these gates, there are
two honest outcomes:

- leave strong actualization open; or
- explicitly adopt GARH-Q: one irreducible actual completed history and law,
  together with the physical lift from quantum outcome strings when that law is
  constructed from the inherited Born cylinders.

GARH-Q is a fallback theory choice, not a consequence of failure of one GARH-D
packet. A stochastic b with an unexplained realized draw is not automatically a
deeper derivation; without further physical content it merely relocates the
GARH-Q primitive.

## 11. Scientific ceiling

The exact advance is structural. It says what endogenous historywise gravity
cannot do under symmetry and exactly what transformation capacity an added
orienting input must have. It neither proves that nature supplies that input nor
chooses between GARH-D and GARH-Q.
