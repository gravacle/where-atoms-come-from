# GL6AT admitted primary research sources

## Source policy

Only original research papers are admitted below.  Reviews, theses, lecture
notes, news, repository summaries, and search-result snippets were not used as
claim authorities.  The arXiv version is supplied for direct public access;
the journal DOI is the publication identity.  Equation numbers refer to the
paper versions named here.

## P1 — exact model and graph/operator mapping

Michael Hermele, Matthew P. A. Fisher, and Leon Balents, “Pyrochlore photons:
The U(1) spin liquid in a S=1/2 three-dimensional frustrated magnet,”
*Physical Review B* **69**, 064404 (2004).

- [Journal DOI](https://doi.org/10.1103/PhysRevB.69.064404)
- [arXiv:cond-mat/0305401v3](https://arxiv.org/abs/cond-mat/0305401v3)
- DOI: `10.1103/PhysRevB.69.064404`

Admitted evidence:

- Eqs. (6)--(9): easy-axis perturbation theory produces the hexagonal ring
  operator, and an explicit unitary changes it to the negative, sign-free
  kinetic convention.
- Text surrounding Eqs. (9)--(11): pyrochlore sites are diamond links;
  `S^z=+1/2` is an occupied dimer; `S_t^z=0` is exactly two occupied links at
  every diamond site; only alternating hexagons flip, while nonflippable
  hexagons are annihilated.
- Eqs. (10)--(11): adding the flippability potential gives the exact
  equal-weight RK ground state at `v/g=1`.
- Eqs. (72)--(75): within the Gaussian Coulomb description, the plaquette
  kinetic-energy density is a composite magnetic-field observable with
  connected equal-time decay `1/R^8`.

Use ceiling:

- The exact graph/Hilbert/operator map and the RK state are admitted.
- The paper establishes an analytic effective-theory Coulomb phase adjacent
  to the RK point, but explicitly treats the pure easy-axis/pure ring point as
  a speculation in its concluding discussion.  It is not authority for a
  rigorous or exact `v/g=0` phase theorem.
- Its third-order XXZ coefficient is not substituted for the independently
  sealed GL6AO coefficient.

## P2 — direct zero-potential numerical phase study

Nic Shannon, Olga Sikora, Frank Pollmann, Karlo Penc, and Peter Fulde,
“Quantum Ice: a quantum Monte Carlo study,” *Physical Review Letters* **108**,
067204 (2012).

- [Journal DOI](https://doi.org/10.1103/PhysRevLett.108.067204)
- [arXiv:1105.4196v3](https://arxiv.org/abs/1105.4196v3)
- DOI: `10.1103/PhysRevLett.108.067204`

Admitted evidence:

- Eq. (2): the dimensionless Hamiltonian is the alternating-hexagon kinetic
  term plus a diagonal flippability counter `mu`, acting on all ice
  configurations.  Restoring GL6AO's coefficient gives `mu=v/g`.  The paper
  identifies zero potential as its minimal three-dimensional quantum-ice
  model.
- The text immediately after Eq. (2) distinguishes this fully-packed-loop
  Hilbert space from the degree-one diamond quantum dimer model.
- Figs. 4--5: zero-temperature GFMC finds a flux-sector level crossing at
  `v/g=-0.50+-0.03`; at `v/g=0`, flux-energy data from 320-, 640-, and
  1280-site QMC clusters and an 80-site exact diagonalization collapse as
  expected for a `U(1)` liquid.
- The reported phase summary calls the evidence unambiguous for `v/g>-0.3`
  and places zero deep in the liquid region; the interval between roughly
  `-0.5` and `-0.3` remains cluster-sensitive.

Use ceiling:

- This is the most direct admitted phase evidence for the exact comparison
  point.  GFMC/ED finite-size evidence is numerical, not a mathematical proof
  of the thermodynamic phase, a pole, or a spectral gap.
- Magnetic language, material interpretations, and degree-violating defects
  are not inherited by the abstract GL6AO link system.

## P3 — effective dispersion and single-link spectral visibility

Owen Benton, Olga Sikora, and Nic Shannon, “Seeing the light: experimental
signatures of emergent electromagnetism in a quantum spin ice,” *Physical
Review B* **86**, 075154 (2012).

- [Journal DOI](https://doi.org/10.1103/PhysRevB.86.075154)
- [arXiv:1204.1325v2](https://arxiv.org/abs/1204.1325v2)
- DOI: `10.1103/PhysRevB.86.075154`

Admitted evidence:

- Eqs. (3)--(5): the same `-g` alternating-hexagon flip plus potential
  Hamiltonian, with the RK point at potential `v=g`.
- Eqs. (61), (67), and (69)--(71): lattice Gaussian theory has two transverse
  physical modes; for positive stiffness the long-wavelength dispersion is
  linear, while the RK limit is quadratic.
- Eqs. (93)--(95): the continuum crossover form is
  `omega(k)=c|k| sqrt(1+(lambda_c|k|/(2pi))^2)`.
- Eq. (105): the low-energy microscopic spin/link structure factor contains a
  one-photon delta function with weight proportional to `omega(k)`.
- Eqs. (112)--(117): finite-size ground-state-energy scaling on 432-, 1024-,
  and 2000-site clusters, combined with the lattice field theory, gives
  `c=(0.6+-0.1)g a_0/hbar` at zero potential; a single-mode calculation gives
  a consistent upper bound.

Use ceiling:

- QMC directly tests equal-time structure factors and finite-size energies.
  The real-frequency delta function and dispersion are field-theory
  predictions calibrated against those data, not direct analytic continuation
  or real-time simulation of the microscopic Hamiltonian.
- The paper expressly retains only low-energy photon contributions in these
  structure factors.  It does not compute the sealed local pair-`E` response.
- Its `k`, `a_0`, and speed require a physical diamond embedding and cannot be
  assigned to an uncalibrated `A3` character.

## P4 — pair/tensor mechanism, with a representation mismatch

Jianlong Fu, Jeffrey G. Rau, Michel J. P. Gingras, and Natalia B. Perkins,
“Fingerprints of quantum spin ice in Raman scattering,” *Physical Review B*
**96**, 035136 (2017).

- [Journal DOI](https://doi.org/10.1103/PhysRevB.96.035136)
- [arXiv:1703.03836v1](https://arxiv.org/abs/1703.03836v1)
- DOI: `10.1103/PhysRevB.96.035136`

Admitted evidence:

- Eqs. (45a)--(45c): the active cubic `T2g` Raman components are complementary
  bond-pair differences `R_01-R_23`, `R_02-R_13`, and `R_03-R_12`.  At the
  paper's approximation order, the `A1g` vertex is proportional to the
  Hamiltonian and has no inelastic response.
- Eq. (57c) in the explicitly linked arXiv v1, renumbered Eq. (55c) in the
  published paper: the Ising part of the Raman operator is a local spin pair
  `S^z_mu S^z_nu`, represented as a bilinear of emergent electric fields.
- Published Eqs. (69)--(72), corresponding to the shifted formulas through
  arXiv-v1 Eq. (74): the gauge-only response creates two opposite-momentum
  photons.  No infrared power is imported: the published prose following
  Eq. (73) reports `Omega^2`, whereas direct three-dimensional power counting
  of published Eq. (70c) / arXiv-v1 Eq. (72c), using a linear `epsilon_p` and
  its displayed nonzero small-`p` form factor, gives `Omega^4`.  This packet
  admits the channel mechanism but treats the exponent as internally
  unresolved in the cited source.
- Section VI in arXiv v1 / Section VII in the published paper: the authors
  enumerate the enlarged-Hilbert-space,
  slave-particle, mean-field gauge, and microscopic Raman-vertex
  approximations.

Use ceiling:

- The result is primary mechanism-level evidence for a pair/tensor composite
  coupling to a two-photon continuum.
- It is not a calculation of the pure projected zero-potential Hamiltonian's
  exact spectral measure.  Its three-dimensional `T2g` complementary-pair
  differences are not the two-dimensional centered complementary-pair sums
  defining the sealed local pair `E` irrep.
- The paper's spinon spectra live in an enlarged XXZ description and are not
  imported into the strict locked GL6AO Hilbert space.

## Search disposition

No admitted primary paper was found that evaluates the exact retarded or
positive-frequency spectral measure of the sealed local pair-`E` operator for
the pure `v/g=0`, degree-two, alternating-hexagon Hamiltonian.  The absence is
a screen result, not a proof that no such paper exists.
