# WAC WORLD OBSERVATION PROTOCOL V001

## T-53A observation-normalization protocol

Status: retrospective, measurement-only, no scientific scoring rule.

For each of the five physical sediment samples, retain two kinds of Lake Shore 8600 VSM
observations:

1. Interpolate the descending and ascending hysteresis branches linearly to zero applied
   field.  Record the resulting magnetic moment after the preceding positive or negative
   saturation history.
2. For every DCD/remanence row, record the applied reverse-field pulse and the magnetic
   moment measured after the field returns to zero.  The source header must state
   `Measure moment at applied fields: False`.

The normalized table uses an ordinal within-sample time because the two source protocols
have separate clocks.  `source_time_s` preserves each protocol's reported relative time.
The readout coordinate is zero tesla; `writer_field_T` carries the preceding writer or
reverse-pulse field.

This protocol classifies the bundle as `WRITE_POST_ONLY`.  It measures writer-conditioned
post-write configurations.  It lacks a pre-write/no-record state, randomized intended
messages, a no-write cohort, a common-hold survival series, and an independently frozen
formation predicate.  Therefore neither this normalization nor a successful ingestion
certificate is a record-formation proof.
