# Exact payload-binding matrix

`GC16--GC19` require all rows below to bind to one calculation.  A check means
the public object owns that item as released; “method only” means code could be
modified to produce it but the released result does not contain it.

| Required object | Shannon/Sikora exact diamond model | Huang 2018 dynamic QSI | Zhou 2026 data | Zhou 2026 code |
|---|---|---|---|---|
| Pure projected `H=-J_6 sum_C B_C`, `mu=0` | yes at the admitted effective model point | no: full XXZ | no: full XXZ-derived tables | no: full XXZ QMC; ED has no native H6 |
| `G_5,G_10,G_20` FS family | published static/flux sizes differ | no: `L=8` only | no: QMC `L=3,4`, ED 16 sites | programmable lattice substrate only; no released matched run |
| fixed zero-flux and complement sector | phase evidence, not event/run sector ledger | no released fixed-sector ledger | no released fixed-sector ledger | must be added and verified |
| complete native pair+H2+H4+H6+ring tensor source | absent | absent (`S^z,S^+/-`) | absent (spin/QFI channels) | absent; custom estimator required |
| GC19 rays plus conjugates | absent as arrays | plotted high-symmetry paths, not GC packet | Gamma, Gamma-prime, X aggregates | arbitrary-momentum infrastructure is suggestive, not a released GC run |
| prospective `1/sqrt(2L^3)` count convention and `Z_Q` | absent | spin normalization differs | QFI density normalization differs | must be declared before measurement |
| full connected `2x2` TT `C_AB(tau,k)` | absent | no; dynamic spin trace/projections only | absent | method only after custom source insertion |
| raw time/channel covariance or independent bins | absent | absent from arXiv source | absent; published `fQ_err` columns are zero in the inspected figure tables | averaged outputs exist; bin/covariance retention must be added |
| one-link photon comparator in the same `J_6` clock | static Maxwell evidence only | photon-like `Szz`, clock `J_z` | integrated QFI, clock `J_z` | must be added and calibrated to F3 `J_6` |
| fitted gaps, full-projector residues, continuum weight | absent | figure-level SAC spectra only | absent | calculable only after the matched run |

## Candidate typing

- **Shannon et al. 2012 / Sikora et al. 2011:** admissible for the existence and
  static properties of the `mu=0` diamond-lattice U(1) phase; not a dynamic
  native-source packet.
- **Huang et al. 2018:** useful proof-of-method that pyrochlore QMC plus SAC can
  expose a photon in a spin channel; wrong Hamiltonian/source and figure-only
  public numerical content.
- **Zhou et al. 2026 data:** genuinely machine-readable, but the released arrays
  are QFI aggregates, not raw correlators or Lehmann data.
- **Zhou et al. 2026 code:** the strongest implementation bridge.  It reduces
  software work but leaves the decisive physics/operator bindings undone.

Therefore no candidate has all rows.  Combining rows from different candidates
would silently mix Hamiltonians, operators, momenta, sectors, and clocks and is
not an admissible GC/FY inference.

