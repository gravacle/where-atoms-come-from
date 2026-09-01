# F3-QIRN microscopic action seed

**Construction ID:** `F3-QIRN-V001`

**Parent choice:** `AY-F3_COMPOSITE_EMERGENT_SPIN2`

**Weinberg--Witten premise route:** `E-EMERGENTSPACE`

**Scale route:** non-temporal microscopic-to-infrared blocking

**Claim class:** explicit candidate microscopic parent plus exact finite-size
properties and ordered falsification program

**Not claimed:** physical instantiation, a continuum limit, spacetime,
curvature, a protected spin-two pole, BRST/diffeomorphism symmetry, Einstein
dynamics, Newton \(G\), record-caused gravity, or objective actualization

## 1. Why this is the first parent

The existing scale-natural pointer/CNOT construction already supplies a finite
same-parent record and causal-composition architecture, but its incidence is
fixed. The smallest origin-bearing change is therefore to make that incidence
quantum dynamical and responsive to the same physical record carriers, while
retaining explicit source/read/reservoir/boundary factors and their custody
obligations.

The candidate is intentionally minimal. V001 has no hand-added spacetime
lattice, loop/Regge action, tensor projector, continuum dimension, or
Einstein-like vertex. If its degree-only quantum-incidence phase has no stable
local scaling window, it is rejected before additional interaction terms are
considered.

## 2. Finite microscopic state space

Fix a system size \(N\), a finite number of directed composition layers
\(\ell=0,\ldots,L\), and a fixed finite number \(K\) of reservoir-fragment
registers per bulk node. A bulk node is \(v=(\ell,i)\), \(i=1,\ldots,N\).

For every possible adjacent-layer arrow

\[
 e:(\ell,i)\longrightarrow(\ell+1,j),
\]

introduce a link qubit \(a_e\) with

\[
 n_e={1-Z^a_e\over2}\in\{0,1\}.
\]

Every node carries finite registers

\[
 (\psi_v,r_v,w_v,b_{v1},\ldots,b_{vK}),
\]

for a transported carrier, retained record, writer, and reservoir fragments.
Source, controller/clock, reader, support, and boundary registers are explicit
port factors.

### 2.1 A unique blank and no hidden record value

Each content-bearing register uses the qutrit space

\[
 \mathcal H_{\rm qutrit}=\operatorname{span}
 \{|B\rangle,|0\rangle,|1\rangle\}.
\]

Define

\[
 q=|0\rangle\langle0|+|1\rangle\langle1|,
 \qquad
 Z^r=|0\rangle\langle0|-|1\rangle\langle1|,
 \qquad Z^r|B\rangle=0.                              \tag{BS01}
\]

Thus \(q=1\) is physical occupation of a record-bearing register; content
exists only on that occupied subspace. There is no state
`blank x unobserved-content` that could carry hidden lineage.

The content exchange

\[
 C=|B\rangle\langle B|+|0\rangle\langle1|
   +|1\rangle\langle0|                              \tag{BS02}
\]

fixes the blank and swaps the two contents. The complete bulk action below is
covariant under simultaneous application of \(C\) to every content register.

Importantly, \(q=1\) is **not** `REC` or `SEALED` by definition. Those remain
history- and mission-level properties requiring independent formation,
custody, persistence/readability, source-off, and no-replay gates.

### 2.2 Relabeling status

The bulk Hamiltonian is invariant under independent within-layer node
permutations

\[
 \mathcal G_N=\prod_{\ell=0}^{L}S_N^{(\ell)}.         \tag{BS03}
\]

No gauge quotient follows from (BS03): states and observables have not been
group-averaged or quotiented. It is an exact finite relabeling
covariance/global symmetry. Comparisons must be permutation-invariant or
prospectively anchored by physical ports. The layer order is a primitive
composition schedule. It is not a spatial coordinate, a metric, a proved
causal order, a Lorentzian time coordinate, or a dated cosmological history.

### 2.3 Initial-state family

The action does not by itself fix a mission state. Each calculation must freeze
one normalized \(\rho_N(x,z)\) containing the link state, unique blanks for all
declared blank record/reservoir targets, prospective writer/source content,
controller and boundary states, and every crossing reference/correlation. Bulk
preparations must transform covariantly under (BS03); content preparations and
reads must transform under (BS02). An unlisted correlation is not traced away.

V001 therefore defines a parametric action family plus an initial/port
completion contract. It is not yet one physically instantiated complete parent.

## 3. Bounded step action

Let \(V_\ell=\{(\ell,i):1\le i\le N\}\),
\(A_\ell=V_\ell\cup V_{\ell+1}\), and let \(P_\ell\) be the prospectively
assigned ports active in that slab. All displayed microscopic coefficients are
real. For each layer use the finite Hermitian generator

\[
 H_{N,\ell}=H_{\rm inc}+H_{\rm car}+H_{\rm form}
             +H_{\rm fb}+H_{\rm port},                \tag{BS04}
\]

and the ordered mission evolution

\[
 U_N=\overleftarrow{\prod_{\ell=0}^{L-1}}
       \exp[-i\delta\tau_0 H_{N,\ell}/\hbar].          \tag{BS05}
\]

The boundary clock interval \(\delta\tau_0\) is a microscopic operational
calibration, not an assumed spacetime coordinate. If the switching is not made
autonomous, its work and timing uncertainty belong to the explicit
controller/clock port in \(H_{\rm port}\).

### 3.1 Incidence sector

Let

\[
 d^{\rm out}_{\ell i}=\sum_j n_{\ell ij},
 \qquad
 d^{\rm in}_{\ell+1,j}=\sum_i n_{\ell ij}.
\]

With \(\Omega>0\) and \(h_N=\Omega/N\), define

\[
\begin{split}
 H_{\rm inc}={}&-h_N\sum_{e\in E_\ell}X^a_e
 +\Delta\sum_{e\in E_\ell}n_e\\
&+U_d\sum_i(d^{\rm out}_{\ell i}-d_*)^2
+U_d\sum_j(d^{\rm in}_{\ell+1,j}-d_*)^2 .             \tag{BS06}
\end{split}
\]

No distance or dimension occurs. The \(1/N\) flip scaling and the degree
penalty make a sparse extensive phase possible without choosing a lattice.
Whether this actually produces a locally geometric phase is a calculation,
not an assumption.

### 3.2 Content-symmetric transfer and copying

For two qutrit registers define the Hermitian content-preserving transfer and
current operators

\[
\begin{split}
 T_{uv}&=\sum_{x=0}^1
 (|B,x\rangle\langle x,B|+|x,B\rangle\langle B,x|),\\
 J_{uv}&=i\sum_{x=0}^1
 (|B,x\rangle\langle x,B|-|x,B\rangle\langle B,x|),   \tag{BS07}
\end{split}
\]

and the blank-target copy generator

\[
 K_{a\to b}=\sum_{x=0}^1 |x\rangle\langle x|_a\otimes
 (|x\rangle\langle B|_b+|B\rangle\langle x|_b).       \tag{BS08}
\]

These operators are invariant under simultaneous content exchange. The copy
operator copies only the declared orthogonal record alphabet; it does not clone
an arbitrary quantum state.

The carrier and formation terms are

\[
 H_{\rm car}=\epsilon_\psi\sum_{v\in A_\ell}q^\psi_v
 -t\sum_{e:u\to v}n_eT^{\psi}_{uv},                   \tag{BS09}
\]

\[
\begin{split}
 H_{\rm form}=\sum_{v\in A_\ell}\bigg[&
 \epsilon_rq^r_v+\epsilon_wq^w_v
 +\sum_{k=1}^K\epsilon_{bk}q^{b_k}_v\\
 &+g_{\psi r}K_{\psi_v\to r_v}
 +g_{wr}K_{w_v\to r_v}
 +\sum_{k=1}^K g_{rbk}K_{r_v\to b_{vk}}
 \bigg].                                               \tag{BS10}
\end{split}
\]

The writer and reservoir couplings can generate content-distinguishable
environmental states while the two contents remain exactly energy-degenerate.

### 3.3 Label-blind occupied-register-to-incidence feedback

The V001 feedback term is

\[
 H_{\rm fb}=-\lambda_R\sum_{e:u\to v}n_e(q^r_u+q^r_v)
 +\lambda_J\sum_{e:u\to v}n_e(J^\psi_{uv})^2,         \tag{BS11}
\]

with \(\lambda_R,\lambda_J\ge0\). It couples incidence to the physical
presence of record-bearing storage and to a physical carrier current. It does
not couple to the record value, a logical `SEALED` flag, a hazard, or \(\gamma\).

Equation (BS11) deliberately predicts degree feedback: a record-bearing node
lowers the energy of each incident occupied link by \(\lambda_R\). This is not
silently normalized away. Runaway degree/occupation is a decisive stability
test of the seed.

No content-parity term \(n_eZ^r_uZ^r_v\), record-pair term, loop/cell term, or
hand-selected spin-two term is present in V001. Such an extension cannot be
added merely to rescue a failed scaling result.

### 3.4 Ports and energy ownership

\[
 H_{\rm port}=H_\partial+
 \sum_{v\in A_\ell,\,p\in P_\ell}g_{vp}O_vO_p        \tag{BS12}
\]

is the declared completion slot for every source, writer supply, reservoir,
clock/controller, reader, support/recoil, and boundary exchange. A physical
instantiation must freeze its matrices and assign every exchange once; the
symbolic slot alone does not prove that completion.
Logical distinguishability, `ALLOW`, success probability, and a SEALED verdict
carry no assigned energy.

Here \(H_\partial=H_\partial^\dagger\), \(g_{vp}\in\mathbb R\), and
\(O_v,O_p\) are Hermitian operators on distinct factors; any same-factor term
is explicitly symmetrized or supplied with its Hermitian conjugate. Attached
content ports transform with the bulk under (BS02), and \(H_\partial\) and the
couplings are jointly content-exchange covariant. A chosen source value is a
boundary condition, not an asymmetric Hamiltonian coefficient.

### 3.5 Finite-size and thermodynamic stability gate

Every fixed-\(N,L,K\) Hamiltonian is bounded because all factors are finite.
For an extensive lower bound, a simple sufficient V001 domain is

\[
 \Delta\ge 2\lambda_R+|t|+\delta_E,
 \qquad \delta_E>0,qquad U_d,\lambda_J\ge0,           \tag{BS13}
\]

with fixed \(K,d_*\) and bounded per-node/port couplings. Since
\(\|T_{uv}\|\le1\), the occupied-edge diagonal plus transfer contribution is
nonnegative up to the positive margin, while

\[
 -h_N\sum_{e\in E_\ell}X_e\ge-\Omega N.               \tag{BS14}
\]

The fixed number of local copy and port terms contributes at worst \(-cN\)
per layer. Thus the ground-energy lower bound is extensive. Parameters outside
(BS13) require a separately proved degree-stabilized bound; they are not
silently admitted.

The V001 storage domain also fixes
\(\epsilon_\psi,\epsilon_w,\epsilon_{bk}\ge0\) and a strictly positive bare
record-storage gap \(\epsilon_r\ge\delta_r>0\). These gaps do not establish
SEALED persistence; they only prevent a free negative-energy record occupancy.

For the port contribution, V001 additionally requires

\[
 \|H_\partial\|\le c_\partial N,
 \qquad \#\{(v,p):g_{vp}\ne0\}\le c_pN,
 \qquad |g_{vp}|\|O_v\|\|O_p\|\le g_{\max},           \tag{BS13a}
\]

with \(N\)-independent constants. Bounded coefficients without this count/norm
control would not imply an extensive lower bound.

An inner-layer register occurs in two consecutive active slabs because it is a
target during one scheduled interval and a source during the next. This is a
defined two-interval evolution, not an unscoped repetition of \(\sum_v\).

## 4. Exact finite-parent symmetries and causal ceiling

Equations (BS06)--(BS12) give, at fixed size:

1. exact unitarity of (BS05);
2. exact within-layer permutation covariance of the bulk;
3. exact bulk content-exchange covariance and complete-parent covariance only
   for the jointly transformed/covariant port family specified after (BS12);
4. a unique physical blank;
5. complete finite-factor boundedness and the sufficient extensive domain
   (BS13); and
6. an explicit ordered composition schedule, but no Hamiltonian
   no-backward-signalling theorem.

The two-register Hermitian gates are bidirectionally signalling on an
unrestricted intervention domain, and simultaneous target-local terms can feed
through the transfer operator. A future causal construction must either freeze
and authenticate a reduced ordered channel

\[
 \Phi_\ell(\rho)=\operatorname{Tr}_{V_\ell}\!\left[
 U_\ell(\rho\otimes\beta_{V_{\ell+1}})U_\ell^\dagger
 \right]                                               \tag{BS14a}
\]

with every target input independent of later choices, or replace the step by a
proved semicausal factorization. Calling an edge an arrow is insufficient.

The present structure is therefore only a composition schedule. The
possible-link set is all-to-all between adjacent layers, so it is not a
recovered physical spacetime cone. Suppression of apparent long links after
emergent localization is an open quantitative obligation.

## 5. Gamma is derived from this action, not inserted into it

For two content preparations \(x=0,1\), contract the complete declared
environmental output of (BS05) to branch states \(\rho_E^{(x)}\) and define the
record-read fidelity diagnostic

\[
 \gamma_{\rm rec}=F(\rho_E^{(0)},\rho_E^{(1)}),        \tag{BS15}
\]

where \(F\) is the squared Uhlmann fidelity, matching the frozen program
convention.

The content-symmetric formation couplings can give \(\gamma_{\rm rec}<1\) because
they correlate different environmental states with the two values. Symmetric
energetics therefore does not mean branch-indistinguishable dynamics.

Gamma never appears on the right side of (BS04). It diagnoses distinguishability
generated by the action. The same action separately supplies the physical
occupation and currents appearing in (BS11). Any observed relation between
gamma flow and incidence response must be derived from this shared parent.

This squared Uhlmann fidelity is not generally the squared magnitude of one
fixed influence amplitude when the conditional environmental branches are
mixed. The two gamma objects are kept distinct in Section 7.

## 6. Composite symmetric-source construction

The following is a source-construction template, not yet an instantiated
composite operator. It must not equate the number of low Laplacian modes with
spacetime dimension. At a frozen blocking scale \(b\), first freeze a blocking
map \(\mathfrak B_b\), define
\(\rho_{b,0}=\mathfrak B_b(\rho_{N,0})\), construct the oriented incidence
operator \(D_b\), and set

\[
 W_b=\operatorname{diag}(w_e),\qquad
 w_e=\operatorname{Tr}(\rho_{b,0}n_e)\ge0.
\]

The positive elliptic diagnostic is then

\[
 L_b=D_b^\dagger W_bD_b\succeq0.                       \tag{BS16}
\]

Freeze an isolated retained projector

\[
 P_b={\bf1}_{[0,\Lambda_b]}(L_b).                      \tag{BS17}
\]

Its rank is merely the number of retained modes. Effective dimension must be
tested independently through local heat-kernel or Weyl scaling on a
prospectively frozen scale window, with finite-size errors and direct/sequential
blocking stability. For example,

\[
 d_{\rm spec}({\cal B};s)
 =-2{\partial\over\partial\log s}
 \log\left[{1\over|{\cal B}|}
 \sum_{v\in{\cal B}}\langle v|e^{-sL_b}|v\rangle\right]. \tag{BS18}
\]

This is an elliptic/spectral diagnostic only. It does not establish Lorentzian
signature, time orientation, or a causal cone.

Before inspecting the scaling result, the calculation must also declare
whether \(L_b\) is a slice operator testing three spatial dimensions or a
Euclideanized full-complex operator testing four effective dimensions. The
target cannot be switched after seeing which plateau appears.

Let \(\{\phi^a_b\}\) be any frozen basis of the retained source space, treated
covariantly under basis changes, and let

\[
 \Delta_e\phi^a=\phi^a(v)-\phi^a(u).                  \tag{BS19}
\]

Write the active-slab Hamiltonian as a finite sum of Hermitian terms
\(H_{N,\ell}=\sum_{\xi\in\Xi_\ell}H_\xi\). Freeze a block assignment
\(\beta_b(\xi)\), edge weights
\(m_\xi^{ab}=\Delta_e\phi^a\Delta_e\phi^b\), and all
node/port weights and contacts before response scoring. Introduce a symmetric
relational strain source by

\[
 H_{N,\ell}[j]=
 \sum_{\xi\in\Xi_\ell}\left[1-
 {1\over2}j_{ab}(\beta_b(\xi))m_\xi^{ab}\right]H_\xi
 +H_{\rm contact}[j],                                  \tag{BS20a}
\]

where \(H_{\rm contact}=O(j^2)\) is Hermitian and retains every required
identity/seagull/boundary contact. The source-conjugate composite operator is

\[
 \mathcal Q^{ab}_{\cal B}
 =-2{\partial H_{N,\ell}[j]\over\partial j_{ab}({\cal B})}
   \bigg|_{j=0}
 =\sum_{\xi:\beta_b(\xi)={\cal B}}m_\xi^{ab}H_\xi.   \tag{BS20}
\]

The term list and weights must include edge, onsite, boundary, controller, and
port contributions rather than dropping them. V001 has no microscopic cell
term; any cell contribution must be generated explicitly by the frozen
blocking map. Until \(\mathfrak B_b,D_b,W_b,\beta_b,m_\xi\), contacts, and the
source quotient are frozen, (BS20a)--(BS20) remain a template. Even when frozen,
\(a,b\) are retained relational-mode labels—not spacetime indices.

The operator, its expectation coordinate, and its response are distinct. Use
the CTP combinations \(j_c=(j_++j_-)/2\) and
\(j_\Delta=j_+-j_-\). On the physical branch \(j_\Delta=0\), define

\[
\begin{aligned}
 W[j_c,j_\Delta]=-i\hbar\log Z[j_c,j_\Delta],\qquad
 q^{ab}_{\cal B}=2{\delta W\over
             \delta j_{\Delta,ab}({\cal B})}\bigg|_{j_\Delta=0},\\
 \mathcal G_R^{ab,cd}({\cal B},{\cal B}')
 ={\delta q^{ab}_{\cal B}\over
   \delta j_{c,cd}({\cal B}')}
 =2{\delta^2W\over
 \delta j_{\Delta,ab}({\cal B})\delta j_{c,cd}({\cal B}')}
 \bigg|_{j_\Delta=0}.                                 \tag{BS21}
\end{aligned}
\]

The noise/Keldysh kernel is a different \(j_\Delta j_\Delta\) derivative. A
retarded pole is a property of \(\mathcal G_R\), not of the definition of
\(\mathcal Q\).

Before \(a,b\) may be renamed spacetime indices, the flow must derive:

- an infrared observable algebra and differential calculus;
- a locally free rank-four cotangent module with consistent gluing;
- a basis-independent localization/intertwining map from the retained source
  space;
- descent of \(\mathcal Q\) through every quotient/gauge-null direction; and
- a retarded hyperbolic principal symbol with one Lorentzian time direction.

Even that descent would not prove a graviton.

## 7. Complete complex CTP route

Use the complete operator/source vector

\[
 \mathbb J=(j_{ab},J_a,J_\psi,J_r,J_w,J_b,
             J_{\rm ports},J_{\rm contacts}).          \tag{BS22}
\]

At every finite size calculate

\[
 Z_N[\mathbb J_+,\mathbb J_-]
 =\operatorname{Tr}\!\left[
 U_N[\mathbb J_+]\rho_NU_N[\mathbb J_-]^\dagger
 \right].                                              \tag{BS23}
\]

Identity, seagull, boundary, and source-contact terms remain explicit.

For a separately declared reduced influence functional obtained by fixing the
system paths and tracing the declared environment, write

\[
 Z_{\rm IF}=e^{iW_{\rm IF}/\hbar},\qquad
 \gamma_{\rm IF}=|Z_{\rm IF}|^2,
 \qquad
 \operatorname{Im}W_{\rm IF}
 =-{\hbar\over2}\log\gamma_{\rm IF}.                 \tag{BS24}
\]

For general mixed conditional environment states there is no identity
\(\gamma_{\rm IF}=\gamma_{\rm rec}\). For example,
\(\rho_E=I/2\), \(U_0=I\), and \(U_1=Z\) give identical conditional states and
therefore \(\gamma_{\rm rec}=1\), while
\(|\operatorname{Tr}(\rho_EU_1^\dagger U_0)|^2=0\). The diagnostics coincide in
the pure conditional-copy substep proved in
`FIRST_EXACT_CALCULATION.md`. Neither determines the full complex CTP phase.

The real phase, active response, and cubic contact packet must come from the
full complex (BS23). Calculate connected response through cubic order and
Legendre-transform only on a prospectively frozen invertible quotient.

## 8. Emergent BRST is a test, never a microscopic input

The microscopic parent has finite relabeling covariance, not diffeomorphism
BRST. In the blocked large-size 1PI action one must search for continuous null
generators

\[
 \Gamma_{,I}R^I_{\alpha}=0                             \tag{BS25}
\]

whose algebra closes and satisfies Jacobi, including boundaries and anomalies.
Only after that result may one introduce ghosts, a frozen gauge-fixing fermion,
and the BRST differential and test the Slavnov--Taylor/AU05/AU06 hierarchy.

Discrete permutation invariance, a tensor-shaped two-point function, or a
closing scalar gap does not imply (BS25). Failure to develop the required null
directions kills the F3 interpretation of this seed.

## 9. Physical mechanism being tested

The candidate mechanism is now explicit:

\[
\begin{array}{c}
\text{writer/carrier transfer and value-covariant copying}\\
\Downarrow\\
\text{physical occupied record carriers + redundant environmental states}\\
\Downarrow\quad(\gamma\text{ diagnoses overlap})\\
\text{label-blind lowering of incident-link energy and current backreaction}\\
\Downarrow\\
\text{possible collective reorganization of incidence under blocking}\\
\Downarrow\\
\text{candidate relational tensor response }\mathcal G_R\\
\Downarrow\\
\text{only if all gates pass: Lorentzian locality, protected spin two,}\\
\text{BRST/Ward closure, universal coupling, and Einstein infrared dynamics.}
\end{array}                                             \tag{BS26}
\]

The first two arrows are a mechanism hypothesis, not yet an instantiated join.
`FIRST_EXACT_CALCULATION.md` proves two compatible but separate sub-sector
witnesses: copying starts with \(r\) occupied, while link response pins \(q^r\)
and switches formation off. The next join must start the same \(r\) register
blank, form it in an authenticated ordered `REC` mission, and then let that
same unrepinned register change incidence. Every collective line remains a
calculation or falsifier.
