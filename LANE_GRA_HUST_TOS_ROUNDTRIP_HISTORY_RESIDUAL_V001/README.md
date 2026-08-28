# HUST ToS round-trip history residual lane

This isolated no-laboratory lane computes equal-configuration endpoint-return
residuals from the pinned official HUST Figure-2 workbook. It does not modify
the source workbook or any canonical model/register file.

Run from the repository root:

```bash
python3 LANE_GRA_HUST_TOS_ROUNDTRIP_HISTORY_RESIDUAL_V001/analyze_hust_tos_roundtrip_history.py
python3 LANE_GRA_HUST_TOS_ROUNDTRIP_HISTORY_RESIDUAL_V001/verify_hust_tos_roundtrip_history.py
python3 LANE_GRA_HUST_TOS_ROUNDTRIP_HISTORY_RESIDUAL_V001/verify_hust_tos_roundtrip_history_hostile_audit.py
```

The observable is a descriptive apparatus/source-history diagnostic on
three-day period summaries. Because there is no matched no-excursion trajectory,
it is specifically a history-confound diagnostic rather than an identified
hysteresis or memory effect. It is not `beta_TM`, a record-lineage intervention,
statistical evidence, gravity emergence, or a new value of \(G\).

The source-present and source-absent loops are paired only by sequence ordinal;
they are not simultaneous, randomized, or authenticated matched trials. Decimal
results reproduce the printed workbook cells, and extra computational digits do
not represent measurement precision.
