# Dependencies and custody

GL6AS uses exactly seventeen frozen objects:

- the GL6AO theorem, author manifest/seal, and independent audit
  theorem/manifest/seal;
- the GL6AP theorem, author manifest/seal, and independent audit
  theorem/manifest/seal; and
- the GL6AQ theorem and author manifest, plus its independent audit
  theorem/manifest/seal.

GL6AQ deliberately has no author `SEAL.sha256`; its author bytes are frozen
by its manifest and pinned by the separately sealed hostile audit.  GL6AS
does not invent a missing author seal.

Every row is pinned in `DEPENDENCIES.sha256`.  No GL6AN object is imported
directly: its content reaches GL6AS only through the audited transitive
custody of AO/AP/AQ.  No mutable reconnaissance lane is a premise.

The exact verifier imports no upstream executable code.
