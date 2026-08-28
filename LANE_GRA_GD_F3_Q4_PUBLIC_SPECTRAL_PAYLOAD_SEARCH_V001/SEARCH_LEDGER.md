# Reproducible bounded search ledger

**Search date:** 2026-08-28  
**Scope:** source-bound ED/QMC arrays or code for the exact diamond quantum-ice /
quantum-dimer H6 model, followed by the closest dynamic pyrochlore sources.  No
generic gravity literature was searched.

## 1. Exact-model roots

1. [Shannon et al., “Quantum Ice: a quantum Monte Carlo study”](https://arxiv.org/abs/1105.4196),
   DOI `10.1103/PhysRevLett.108.067204`.  Queries used the exact title, DOI,
   `data`, `code`, `repository`, `supplementary`, GitHub, and Zenodo.  The
   downloaded arXiv v3 source archive has SHA-256
   `e0ed17fb05a443c8adc04a04cdf122f0bdd9da3631b0a1f4ab5be180d05daad6`.
   Inventory: one TeX file, six EPS figures, and one figure PDF.  It contains
   no `.dat`, CSV, HDF5, notebook, or program source.  Physics use: `mu=0` flux
   scaling and static small-wave-vector structure factor; not GC16.
2. [Sikora et al., “Extended quantum U(1)-liquid phase in a
   three-dimensional quantum dimer model”](https://arxiv.org/abs/1105.1322),
   DOI `10.1103/PhysRevB.84.115129`.  The downloaded arXiv source archive has
   SHA-256 `c69f465001be8bbcf77446d5f8258781136bda070e8d8d058158d42e73817521`.
   Inventory: one TeX file and 25 EPS figures; no raw arrays or code.  Physics
   use: exact diamond-QDM phase/flux evidence; not a response packet.

## 2. Closest published dynamic calculation

3. [Huang et al., “Dynamics of Topological Excitations in a Model Quantum
   Spin Ice”](https://arxiv.org/abs/1707.00099), DOI
   `10.1103/PhysRevLett.120.167202`.  The downloaded arXiv source archive has
   SHA-256 `da3dce95864b8ae4551a38b8e9f64a571fcefc24c8f23924bce85162c0e670f4`.
   Inventory: `main.tex` and five PDF figures only.  The paper reports QMC-SAC
   spin spectra for the full XXZ model on `L=8`, including a photon-like
   `Szz(q,omega)` branch, but releases neither the raw `G(q,tau)`/covariance nor
   code in the source archive.  Its spin operator and `J_z` clock do not bind to
   the complete FY/GC tensor source and `J_6` clock.

## 3. New machine-readable and executable root

4. [Zhou et al., “Quantum Fisher information as a thermal probe in frustrated
   magnets through insights from quantum spin ice”](https://www.nature.com/articles/s41467-026-74589-6),
   DOI `10.1038/s41467-026-74589-6`.  The publisher explicitly points to:
   - [data DOI `10.25442/hku.32404548.v1`](https://doi.org/10.25442/hku.32404548);
   - [code DOI `10.25442/hku.32412273.v1`](https://doi.org/10.25442/hku.32412273).
5. The data API metadata SHA-256 is
   `7e2e4bb8db3bb529970ed63794654e6666d7cfb902c3337c4baa23d1c6a44f46`.
   Ten downloaded text tables total 255,599 bytes.  Their headers are
   `x_index x_value T value` or `source_id T fQ fQ_err`; they are derived QFI
   figure data, not spectra or imaginary-time arrays.  Some rows carry a scalar
   `fQ_err` and others carry zero, but no file supplies the time-by-time,
   channel-by-channel covariance needed to refit GC16.
6. The code API metadata SHA-256 is
   `cc41ffa163c01a63d47a9bd493ad74638434efe39639b7bad5f91ffad6ebd568`.
   Static inspection of the downloaded files found:
   - `Pyrochlore.f90` and `Pyrochlore_Con.f90`: SSE/MDL-QMC with real- and
     imaginary-time sublattice spin correlator arrays and averaged output;
   - `initmulflux.f90`: multi-flux initialization/update support;
   - `QED-main.zip`: ED/Lanczos/DSSF framework, archive comment commit
     `f672e0c9e6648c90c9c9564bbd8d93d0061d7899`.
   The QED `HamiltonianBuilder` owns only one-/two-/three-body terms and its
   `ring_exchange` implementation raises a runtime error because even four-body
   exchange is unsupported.  No six-body H6 term or GC/FY source is present.

## 4. Repository-native bounded searches

7. Figshare API `articles/search`, query `quantum spin ice`, limit 100.
   Response SHA-256
   `5c1daf429dbba5c64d1b109e2a0f45a38c521ddfe83a73e3561954343dcafc28`.
   The only directly relevant numerical/code pair was the Zhou/HKU release.
   Record 9409685 was a 2D Rydberg journal PDF only.
8. Zenodo API query `quantum spin ice`, size 25.  Response SHA-256
   `d30193745585daf5618c83b24346ebd127275c8cb94934eb745d4b98acac897f`.
   Relevant-looking hits concerned 2D Rydberg, classical/artificial spin ice,
   or unrelated geometries; none bound to the 3D diamond H6 tensor response.
9. Zenodo API query `quantum dimer diamond`, size 25.  Response SHA-256
   `6188c261ddde2be553af6fcd25fd39c8bfec2b6489a996eac77a83eef0cab0a2`.
   No exact diamond-lattice pure-H6 spectral packet was returned.

## 5. Search ceiling

This was a bounded search of exact titles/DOIs, publisher supplements, arXiv
source archives, GitHub-indexed results, Figshare, and Zenodo.  It does not
establish absence from every server, private author archive, or future release.
Its positive result is exact: the inspected public objects and their released
fields do not instantiate GC16--GC19, while the 2026 code release is a credible
starting substrate for a new matched calculation.
