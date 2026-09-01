# GL6AH reviewed dependency provenance

GL6AH uses the conventions and already audited conclusions of these frozen
packets; it does not edit them:

- `LANE_CROSS_RFT_GRA_GL6AA_RECORD_AUTHENTICATED_SHARED_CHILD_ATLAS_V001/`
  with hostile audit `AUDIT_G_GL6AA_RECORD_AUTHENTICATED_SHARED_CHILD_ATLAS_V001/`;
- `LANE_CROSS_RFT_GRA_GL6AB_E2_MULTI_CONNECTOR_TRIANGLE_V001/`
  with hostile audit `AUDIT_G_GL6AB_E2_MULTI_CONNECTOR_TRIANGLE_V001/`;
- `LANE_CROSS_RFT_GRA_GL6AF_FORMATION_PATTERN_E2_SOURCE_V001/`
  with hostile audit `AUDIT_G_GL6AF_FORMATION_PATTERN_E2_SOURCE_V001/`;
- `LANE_CROSS_RFT_GRA_GL6AG_N1_FORMATION_NEIGHBOR_PROPAGATION_V001/`
  with hostile audit `AUDIT_G_GL6AG_N1_FORMATION_NEIGHBOR_PROPAGATION_V001/`;
- `LANE_CROSS_RFT_GRA_GL6Z_MULTICELL_CTP_K2_JET_V001/`
  with hostile audit `AUDIT_G_GL6Z_MULTICELL_CTP_K2_JET_V001/`.

Imported conventions are the physical-`K` branch typing, GL6AB fixed `E`
matrix and `w_ab` columns, GL6Z pair order and incidence representation,
GL6AA authenticated shared-child indexing, and the exact GL6AG `N=1`
coefficients used only as a reconciliation target.  Both supplied GL6AH
replays reconstruct their coefficients from Hamiltonians and do not import
numeric response values from those packets.

This is a readable provenance note, not a hash freeze.  Dependency hashes and
a packet manifest must be generated only after hostile review authorizes a
freeze.
