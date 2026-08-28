# Independent hostile audit: inherited F3/q4 order-eight TT-kernel boundary

**Lane:** `GRA-FM-F3-Q4-IEKB-V001`

**Audit date:** 2026-08-27

**Disposition:** `PASS_WITH_SCOPED_CONTINUUM_CORRECTION`

## 1. Independence, custody, and scope

I audited the frozen `d_*=2`, `E_R=0` F3 single-link-flip parent rather than
assuming the packet's conclusion. I checked the bytes of every dependency
listed in `THEOREM.md`; all seven SHA-256 values match. I then reproduced the
order-eight endpoint census, the full alternating-octagon path sum, the
diagonal classification, and the strict response expansion. No shared model
or register file was used as evidence.

The fixed-support lattice result survives. One continuum sentence did not:
the original draft attributed low-momentum analyticity to the dressed 1PI
kernel and could be read as a universal statement that a finite-order term
cannot create a pole. I required a scoped correction. Analyticity is earned
for the bare finite-range Fourier vertex. A massless dressed 1PI function may
have threshold logarithms or other nonanalyticities, and a Hamiltonian
truncated at finite order can in principle bind or shift a pole when it is
solved or resummed nonperturbatively. The corrected theorem now makes the
narrow exact statement: a strict single-insertion susceptibility does not
itself contain a new Bethe--Salpeter denominator and therefore does not
establish a new isolated tensor pole.

## 2. Independent order-eight path classification

For an eight-letter flip word, the set of links with odd multiplicity is the
symmetric difference between the initial and final ice configurations. Its
degree is even at every vertex. On a simple bipartite graph of girth at least
six, a nonempty even-degree subgraph supported by at most eight edges is:

- one alternating hexagon when the odd support has six edges; or
- one alternating simple octagon when the odd support has eight edges.

Two cycles require at least twelve edges. A non-cycle connected Eulerian
graph with a degree-four vertex also requires more than eight edges under the
same girth bound. Thus the off-diagonal list is complete. The six-edge case
uses two extra letters: either one hexagon link has multiplicity three, or
one link outside the odd support has multiplicity two. Both are dressed
hexagon transitions. The eight-edge case has eight distinct links. The
diagonal multiplicity partitions are exactly

`(8)`, `(6,2)`, `(4,4)`, `(4,2,2)`, and `(2,2,2,2)`.

The `(8)` class returns to `P_2` after its first pair and is folded, as stated.
This reproduces the theorem's exhaustive scalar / dressed-hexagon / new-
octagon operator topology.

## 3. Independent `8!=40320` octagon coefficient and sign

For an alternating cycle, every nonempty proper toggled edge subset `S` is a
union of path components. Internal degree changes cancel along each path;
only its two endpoints remain. Hence, in units of `U_d`,

\[
 \Delta(S)/U_d=2c(S),
\]

where `c(S)` is the number of connected components of `S` on the octagon.
This provides an independent denominator rule not based on the draft's
thirteen-class table. A `2^8` subset recursion using that rule gives

\[
 \sum_{\pi\in S_8}\prod_{r=1}^{7}{1\over 2c(S_r(\pi))}
 ={429\over16}.
\]

Direct enumeration of all `8!=40320` orders separately reproduces the
thirteen denominator classes and every multiplicity printed in the theorem.
Every proper subset is a forest and has positive gap; the full set alone
returns to ice. Eight factors of `-h` and seven negative resolvents give the
matrix-element sign `-`, so

\[
 \langle n\triangle C_8|H_{\rm eff}^{(8)}|n\rangle
 =-{429h^8\over16U_d^7}.
\]

Because no lower kernel has an octagon endpoint, there is no lower-order fold
with which this first-appearance coefficient can mix.

## 4. Diagonal scalar and `V_8=0`

An irreducible diagonal eight-word uses at most four distinct links. Its
denominators see only the occupation colors of those links and their shared-
vertex incidence. Girth at least six makes every such incidence graph a
forest. A color-preserving connected-tree embedding count is fixed by `z=4`,
`d=2`, and the colors of the already used incident edges; collision would
create a cycle of length at most four. Counts of disconnected types then
follow by inclusion--exclusion from those fixed connected counts. The
colored four-edge census lemma is therefore sufficient and no external link
can enter the gap formula.

As an executable hostile witness, I reproduced the two degree-two states on
the `PG(2,3)` incidence graph, the equality of their colored two-, three-, and
four-edge censuses, and the exact independent irreducible diagonal checksum

`2526594309109/13608000`

on both sides of the hexagon flip. Through order four the lower kernels and
their energy derivatives are scalar. At order six the only non-scalar lower
operator is the hexagon transition, so its order-eight self-consistency fold
dresses that transition rather than generating a diagonal potential. No
square of the order-six operator enters as early as order eight. I therefore
accept the fixed-Feshbach-convention result

\[
 H_{\rm diag}^{(8)}=\epsilon_8P_2,\qquad V_8=0.
\]

This is a fixed finite-support statement, not a volume-uniform all-orders
claim.

## 5. Why order six is already interacting

The leading term is a compact, hard-core, flippability-projected ring
Hamiltonian. Its action depends on the constrained configuration and is not
an exact quadratic oscillator Hamiltonian. It is therefore already an
interacting/non-Gaussian microscopic model at order `h^6`. Gaussian Maxwell
describes the asserted infrared fixed point. The octagon is the first new
loop topology, not the first interaction, and the packet now states this
correctly.

## 6. Correct pole boundary and remaining no-lab calculation

The corrected boundary is precise. A bare finite-range octagon or dressed-
hexagon insertion has an analytic lattice Fourier vertex. The strict
single-insertion expression

\[
 \Pi_2+\Pi_2{\cal K}_{TT}\Pi_2
\]

does not contain the resummed denominator
`(1-K_TT Pi_2)^{-1}` and cannot by itself demonstrate a new isolated pole.
That fact is not a nonperturbative no-pole theorem. Repeated insertions or an
exact spectral solution of the finite-order Hamiltonian remain capable in
principle of producing binding, and Dyson resummation may move an existing
pole.

The actual next target is now stated in the right order:

1. calculate the commonly normalized connected TT four-point function of the
   fixed pure-ice Hamiltonian, with q4 cell volume, source normalization,
   external-state convention, and field residues fixed;
2. amputate the external photon legs and separate pieces reducible in the
   selected two-photon channel to obtain the two-particle-irreducible kernel;
3. solve the momentum-dependent Bethe--Salpeter equation or the corresponding
   finite-volume spectrum;
4. test isolated-level survival and nonzero thermodynamic TT residue, a
   common linear cone, the helicity-two representation, and a protecting
   rank-two Ward/constraint identity.

A nonzero scalar matching coefficient is only the first discriminator. It is
not a pole, and a pole without the Ward, cone, and residue tests is not the
claimed gravity carrier.

## 7. Final disposition

After the scoped correction, the packet earns the following and no more:

- complete order-eight endpoint/operator classification on the declared
  finite supports;
- exact alternating-octagon coefficient `429/16` and negative matrix-element
  sign;
- dressed hex transitions plus no configuration-dependent diagonal,
  `V_8=0`;
- the recognition that the `h^6` compact ring model is already interacting;
  and
- a strict single-insertion boundary that identifies, but does not prejudge,
  the connected-four-point, resummation, Ward, and residue calculation.

**Final audit result:** `PASS_WITH_SCOPED_CONTINUUM_CORRECTION`.
