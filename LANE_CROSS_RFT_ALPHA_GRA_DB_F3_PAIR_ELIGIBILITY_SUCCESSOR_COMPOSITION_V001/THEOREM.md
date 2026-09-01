# Pair-eligibility / successor-incidence composition and carrier-force obstruction

**Lane ID:** `CROSS-DB-F3-PESC-V001`

**Date:** 2026-08-26

**Claim class:** exact finite composition of an already-formed pair-memory
eligibility field with the existing F3 successor-incidence field; exact
BS09 carrier-graph typing; exact finite formation/ownership and Hilbert census;
exact obstruction to importing the CTSSF degree-six carrier margin

**Disposition:**
`PAIR_MEMORY_TO_SUCCESSOR_ELIGIBILITY_COMPOSITION_EXACT__CUBIC_DEGREE_THREE_CZ_TYPING_LAWFUL__UNCHANGED_BS09_CARRIER_PROPAGATES_ON_G_N_NOT_G_ELL__CTSSF_DEGREE_SIX_MARGIN_NOT_IMPORTABLE_WITHOUT_NEW_KINETIC_LAW__FINITE_PROGRAMMED_FORMATION_ONLY__NO_PHASE_SELECTION_OR_GRAVITY`

**Not claimed:** a direct retained-memory-to-carrier hopping law; a second
saturated support field; autonomous production or thermodynamic selection of
cubic support; a physical two-switch after the successor field is occupied;
an `N`-uniform support phase; a cubic quantum-link phase; visible
electromagnetism; alpha; a metric, spin two, gravity, or `G`

## 1. Freeze the unchanged factors and ask the exact typing question

Let `V` contain `N` physical vertex identities and

\[
 \mathcal E=\{\{u,v\}:u<v\},\qquad
 P=|\mathcal E|={N(N-1)\over2}.                     \tag{DB01}
\]

For every potential pair `e`, retain the exact FPMH binary factors `L_e,K_e,G_e`
and its link qubit `a_e`.  Reuse that link factor literally as the successor
incidence field:

\[
 n_e:=|1\rangle\!\langle1|_{a_e},\qquad X_e^n:=X_{a_e}.          \tag{DB02}
\]

The exact FPMH/ALLOW gate is therefore

\[
 H_{KX}=-h\sum_{e\in\mathcal E}|1\rangle\!\langle1|_{K_e}
              \otimes X_e^n.                       \tag{DB03}
\]

On a retained pair-memory eigenpath `K_e=ell_e`, define

\[
 G_\ell=(V,E_\ell),\qquad
 E_\ell=\{e:\ell_e=1\},                            \tag{DB04}
\]

and the distinct instantaneous incidence graph

\[
 G_n(\ell)=\bigl(V,\{e\in E_\ell:n_e=1\}\bigr).     \tag{DB05}
\]

The question is whether `G_ell` can both be the degree-six cubic eligibility
skeleton for a degree-three `n` field and retain the carrier force calculated
by CTSSF on the full degree-six graph, without adding an interaction.

## 2. Exact admitted source-off sector and unchanged parent

Every noneligible `n_e` begins blank.  Since (DB03) has no flip when
`ell_e=0`, the prospectively declared blank-nonedge sector is invariant:

\[
 \mathcal H_n^{\rm adm}(\ell)
 :=\operatorname{span}\{|n\rangle:n_e=0
       \text{ whenever }\ell_e=0\},
 \qquad \dim\mathcal H_n^{\rm adm}(\ell)=2^{|E_\ell|}.          \tag{DB06}
\]

All conditional propagators, traces, and free energies below are restricted to
the complete admitted hold sector `H_adm(ell)` with `K=ell`, (DB06), and the
matched memory/port block.  A trace over the raw `P` successor qubits would
include a spectator factor `2^{P-|E_ell|}` and is not the declared mission.

With `U_d>0` and `|E_R|<2U_d`, the already-declared lineage-gated F3 link
Hamiltonian is

\[
 \begin{split}
 d_v^n(\ell)&=\sum_{e\in E_\ell:e\ni v}n_e,\\
 H_n[\ell]&=E_R\sum_{e\in E_\ell}n_e
 +U_d\sum_v\bigl(d_v^n(\ell)-3\bigr)^2
 -h\sum_{e\in E_\ell}X_e^n.                       \tag{DB07}
 \end{split}
\]

Now retain the unchanged BS09 qutrit carrier `psi_v`.  Its hopping term is not
controlled directly by `K`.  It is controlled by the instantaneous incidence
occupation:

\[
 \boxed{
 H_{\rm car}[n]
 =\epsilon_\psi\sum_vq_v^\psi
 -t\sum_{e\in E_\ell}n_eT_e^\psi.}                 \tag{DB08}
\]

Consequently the exact active carrier graph is

\[
 \boxed{G_\psi(\ell,n)=G_n(\ell),\quad\text{not }G_\ell.}       \tag{DB09}
\]

The BS11 current-square term, when retained, is also proportional to `n_e`;
it does not create a direct `K_eT_e^psi` channel.  Source/writer/route factors
may decouple during the hold, but `n` and `psi` do not form a tensor-sum
Hamiltonian because (DB08) is an operator-valued incidence--carrier coupling.
The strongest exact partition factorization is only

\[
 Z_{\rm DB}^{\rm adm}(\ell)
 =Z_{\rm mp}(\ell)Z_{n\psi}(\ell),                 \tag{DB10}
\]

not `Z_mp Z_n Z_psi`.

On the lawful `h=0`, `lambda_J=0` comparator inside the degree-three locked
sector, the exact joint factor is

\[
 Z_{n\psi}(\ell)
 =\sum_{n\in\Omega_3(G_\ell)}
 e^{-\beta E_{\rm diag}(n)}\,
 \Xi_\psi(G_n),                                    \tag{DB11}
\]

where `Xi_psi(G_n)=Tr_psi exp[-beta H_car(n)]`.  Even if homogeneous detuning
makes `E_diag` scalar on `Omega_3`, the carrier partition function remains
inside the sum.  At nonzero `h`, incidence transitions and carrier propagation
remain one joint trace.

### Theorem PESC-1 -- exact unchanged-parent typing

An already-formed retained pair-memory eigenpath lawfully supplies eligibility
for the distinct FPMH/F3 successor field `n=a`.  In the unchanged parent, the
BS09 carrier then propagates on `G_n`, not on the whole eligibility graph
`G_ell`.  The pair memory gates `X_n`; it does not gate carrier transfer.

## 3. Cubic degree-three successor typing is lawful

Condition on `G_ell=T_L^3` with even `L>=8`.  Then

\[
 N=L^3,\qquad |E_\ell|=m=3N,\qquad d_v(\ell)=6.     \tag{DB12}
\]

Let

\[
 \Omega_3(G_\ell)
 =\{n\in\{0,1\}^{E_\ell}:d_v^n(\ell)=3\ \forall v\}.            \tag{DB13}
\]

Its exact global basis count is

\[
 |\Omega_3(G_\ell)|
 =\left[\prod_{v\in V}x_v^3\right]
   \prod_{e=\{u,v\}\in E_\ell}(1+x_ux_v).          \tag{DB14}
\]

The local count is `binom(6,3)=20`, but the global dimension is (DB14), not
`20^N`.  Every finite six-regular bipartite graph decomposes into six perfect
matchings; the union of any three is a member of (DB13).  A simple bipartite
degree-preserving two-switch remains six-regular and therefore retains a
nonempty degree-three sector.

On cubic support this construction is explicit: orient from the even to odd
sublattice.  The six displacements `+/- e_i` are six disjoint perfect
matchings, and any three give a locked incidence state.

At each fixed finite `L`, in the controlled `CROSS-CZ` perturbative domain,

\[
 G_v=\eta_v[d_v^n(\ell)-3]=0                       \tag{DB15}
\]

and the leading projected square move is the exact spin-`1/2` quantum-link
plaquette Hamiltonian.  This imports no volume-uniform remainder or phase.

### Theorem PESC-2 -- finite cubic eligibility-to-CZ bridge

The retained degree-six cubic pair graph can lawfully be the eligibility
skeleton for the distinct degree-three CZ incidence field.  This closes the
static Hilbert/ALLOW interface only; it does not retain the CTSSF carrier
margin, select the support, or prove a cubic phase.

## 4. Exact finite pair-memory formation and ownership

There is an exact finite FPMH witness for any supplied finite edge set `E_*`.
For every `e={u,v}`, choose an arbitrary prospectively fixed orientation
`a_e -> b_e`; on bipartite/cubic support specialize it even-to-odd.  Allocate
two distinct writer qutrit incidence slots

\[
 A_e\equiv(a_e,e),\qquad B_e\equiv(b_e,e),           \tag{DB16}
\]

physically owned by the displayed endpoint vertices.  The owner map
`o(A_e)=a_e`, `o(B_e)=b_e` is frozen in the source/transducer port census.
Thus the stored pair relation is between the common vertex identities, not
between hidden slot vertices.

Let `rho_*` be an arbitrary joint state on the `m` two-dimensional writer
content subspaces; it may be entangled between edges.  Formation places each
token in `A_e` with `B_e` blank; the matched sham places the same joint state
in the `B_e` slots with `A_e` blank.  Apply

\[
 U_*=\bigotimes_{e\in E_*}U_{A_e\to B_e}.           \tag{DB17}
\]

The factors act on distinct slots and targets.  Exact FPMH linearity gives

\[
 \begin{aligned}
 F:&\quad \rho_*^A\longmapsto
 |B\cdots B\rangle\!\langle B\cdots B|_A
 \otimes\rho_*^B\otimes|1\cdots1\rangle_L,\\
 S:&\quad \rho_*^B\longmapsto
 |B\cdots B\rangle\!\langle B\cdots B|_A
 \otimes\rho_*^B\otimes|0\cdots0\rangle_L.        \tag{DB18}
 \end{aligned}
\]

Parallel KEEP routes `L_e` into `K_e` and gives `ell=1_{E_*}`; nonedges stay
blank.  The two-stage finite serial schedule is

\[
 U_{\rm serial}(\tau)
 =e^{-iH_{\rm DB}^{\rm adm}\tau/\hbar}U_{\rm KEEP}U_*.          \tag{DB19}
\]

Stage I writes and routes `ell`; Stage II turns every source/writer/route pulse
off, conserves `K`, and evolves the unchanged joint `n+psi` parent.  No token
is cloned and no Stage-I output is copied when it becomes a Stage-II control.
For cubic `E_*`, this programs the target graph in the custody schedule.  It
is finite physical-model instantiation, not autonomous graph selection.

## 5. Exact factor and port census

For cubic `m=3N`, the raw finite device contains:

- `2m=6N` endpoint-owned writer qutrit incidence slots;
- `N` BS09 carrier qutrits `psi_v`;
- `3P` pair-memory bits `(L,K,G)`;
- `P` successor/link bits `n=a`;
- `m` prospectively fixed route coordinates; and
- a finite `H_aux` containing owner-preserving source transducers, controllers,
  clocks, work/reference, complete outcome/failure history, and boundary ports.

If `D_aux=dim H_aux` excludes the displayed route coordinates,

\[
 \boxed{D_{\rm raw}=D_{\rm aux}\,3^{N+2m}2^{4P+m}
 =D_{\rm aux}\,3^{7N}2^{2N(N-1)+3N}.}              \tag{DB20}
\]

The fixed full-KEEP cubic admitted `n+psi` factor has dimension

\[
 3^N2^m,                                             \tag{DB21}
\]

and its degree-three locked subspace has dimension

\[
 3^N|\Omega_3(G_\ell)|.                             \tag{DB22}
\]

The successor `n` is the FPMH link `a`; quarantine is retained; and writer
slots are distinct from the BS09 carrier.  No factor is counted twice.

## 6. Exact CTSSF margin obstruction

The CTSSF carrier calculation uses the saturated carrier graph

\[
 H_\psi^{\rm sat}[\ell]
 =\epsilon_\psi\sum_vq_v^\psi
 -t\sum_{e\in E_\ell}T_e^\psi,                    \tag{DB23}
\]

which is the BS09 slice `n_e=1` on every support edge.  On cubic support this
has active degree six.  Its strict finite-volume switch margin is a theorem
about (DB23).

The CZ locked sector instead has exactly three occupied incidence edges at
every vertex and therefore uses (DB08).  For every `n in Omega_3(G_ell)`,

\[
 H_\psi^{\rm sat}[\ell]-H_{\rm car}[n]
 =-t\sum_{e\in E_\ell}(1-n_e)T_e^\psi\ne0.          \tag{DB24}
\]

The mismatch is spectral, not a constant ledger term.  With adjacency
matrices,

\[
 \operatorname{Tr}A_\ell^2=2|E_\ell|=6N,
 \qquad
 \operatorname{Tr}A_n^2=2|E(G_n)|=3N.              \tag{DB25}
\]

Thus even the one-carrier high-temperature expansion differs at second order.
The full-support closed-return counts that produce the CTSSF restoring margin
are not the closed-return counts of the degree-three carrier graph.

FPMH supplies `K_eX_e^n`; it does not supply `K_eT_e^psi`.  `GRA-CH` labels
`-tA_ell` a saturated one-carrier **diagnostic** and expressly withholds the
claim that the full action saturates every allowed link.  `GRA-DA` independently
proves that the physical BS09 carrier graph is `G_n` and that saturated
degree-six and CZ degree-three slices cannot be added.

It follows that the tempting ledger

\[
 F_{\rm CTSSF}(\ell)+F_{\rm CZ}(\ell)               \tag{DB26}
\]

is not an unchanged-parent free energy.  It mixes incompatible `n` sectors or
silently adds another carrier channel.  The CTSSF margin cannot be imported
into the lawful bridge proved in Sections 2--3.

### Theorem PESC-3 -- exact no-import theorem

In the unchanged FPMH+F3 parent, retained `ell` can gate the CZ successor
field, but the CTSSF degree-six carrier margin does not survive as a separate
additive force.  The same-parent support endpoint observable is instead

\[
 F_{n\psi}^{\rm adm}(\ell)
 =-\beta^{-1}\log
 \operatorname{Tr}_{\mathcal H_{\rm adm}(\ell)}
 e^{-\beta[H_n[\ell]+H_{\rm car}[n]+H_{\rm fb}+H_{\rm ports}]}, \tag{DB27}
\]

and must be evaluated anew on matched two-switch endpoints.

## 7. The genuinely new kinetic ingredient, not installed

To retain a carrier on every `ell` edge while `n` remains degree three would
require at least one new physical ingredient, for example

\[
 H_{KT}^{\rm new}=-t_s\sum_eK_eT_e^\chi             \tag{DB28}
\]

for a separately declared carrier `chi`, or a separately formed saturated
support-link field `s_e` with its own custody and `-t_s sum_e s_eT_e`.  Either
choice changes the parent, adds a physical factor or coupling, and requires a
fresh no-double-count and port audit.  Neither is implied by FPMH, BS09, or the
CH diagnostic, and neither is adopted here.

The shortest unchanged-parent continuation is therefore physics, not
machinery: calculate (DB27) for the same cubic-support two-switches, including
the degree-three incidence sum, BS09 carrier return, occupied-`n` custody, and
matched port/work ledger.  That calculation may reveal a new support force;
it may not borrow the CTSSF sign.

## 8. Relabelling covariance and final ceiling

Under a vertex permutation `pi`, simultaneously permute `K,L,G,n,psi` and
each endpoint-owned writer slot `(v,e)->(pi v,pi e)`, including its owner map.
Then the admitted parent is unitarily covariant:

\[
 U_\pi H_{\rm DB}^{\rm adm}[\ell]U_\pi^\dagger
 =H_{\rm DB}^{\rm adm}[\pi\ell].                   \tag{DB29}
\]

No distance, cubic reward, or target-grid energy has been inserted.  The
cubic graph appears only as a supplied retained-memory state and programmed
finite formation schedule.  The result closes the finite eligibility/CZ
typing bridge and sharply obstructs the attempted carrier-force join.  It does
not prove support selection, a phase, electromagnetism, geometry, or gravity.

