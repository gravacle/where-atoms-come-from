# T-53A — first raw world-surface write observation

This lane exercises a standalone measurement adapter on raw instrument data from a
physical magnetic record surface.  The adapter is not yet a public URM input door.  This
lane does **not** complete T-53, prove a general theory, or license a universal claim.

## Actual surface and source

- Surface class: magnetic remanence in five discrete sediment samples.
- Instrument: Lake Shore 8600 vibrating-sample magnetometer at the Institute for Rock
  Magnetism, University of Minnesota.
- Source: Reilly et al., *IODP Site U1537 1.7–3.3 Ma Rock Magnetic Data*.
- Dataset DOI: `10.5281/zenodo.14564186`.
- Related paper DOI: `10.1029/2025PA005360`.
- License: CC BY 4.0.

The retained raw inputs are the five `Step 1 Hysteresis Measurement.csv` files and their
five matching `Step 2 Remanence Curves.csv` files.  `acquire.py` extracts only those ten
files from a source archive and verifies each source checksum against the retrieved Zenodo
record metadata.

## Physical question measured

For each sample, the hysteresis loop applies both signs of a saturating field and reads the
magnetic moment again near zero applied field.  The DCD/remanence sequence starts after
positive saturation, applies progressively stronger reverse-field pulses, and measures the
remanent moment at zero field after every pulse.  Consequently the raw data can report:

1. whether opposite writer directions leave distinguishable remanent states after the
   writer field is removed;
2. how much remanence survives small reverse-field perturbations; and
3. the reverse-field scale at which the remanent sign changes.

Those are observations of three candidate formation terms—writer intervention,
post-writer persistence, and a protection threshold—on a real surface.  The files do not
contain randomized blinded trials, an independent no-write cohort, long-time retention, or
multiple observer fragments.  They therefore cannot establish the full formation process,
the necessity of the amended clauses, or universality.

## Custody and scoring

This source was selected and inspected after its outcomes were public.  All numerical
results are **retrospective/exploratory**.  No tolerance learned from these files may be
used to score a confirmatory claim.  A later structurally different surface must use a
frozen rule and blind holdout.

The lane must emit measured quantities and scope only.  It must never emit `PROVED`, an
ALL-PASS theory verdict, or a T-53 status change.
