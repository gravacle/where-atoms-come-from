# GL6AH — frozen bulk/boundary irrep separation

This author-frozen packet separates the earliest direct-edge pair-mean
signal from the later fixed-`E` terms previously seen in GL6AG.  The direct
order-six signal is an equal `A1+T2` vector on every general-`N` edge.  The
displayed order-twelve and order-sixteen `E` terms are sums of
receiver-helper chains and cancel at a complete homogeneous endpoint.

Run the local exact replay:

```sh
python3 -B verify_local_bulk_boundary.py
```

Run the full `N=1` `2^6` connector-support census:

```sh
c++ -O3 -std=c++17 -I/opt/homebrew/include -L/opt/homebrew/lib \
  verify_n1_connector_supports.cpp -lgmpxx -lgmp \
  -o /tmp/verify_gl6ah_supports
/tmp/verify_gl6ah_supports
```

Or run both plus packet/scope checks with:

```sh
python3 -B verify_mutable_packet.py
```

Status: independently hostile-reviewed clean and author-frozen.  The manifest
and seal pin these bytes.  A distinct post-freeze custody/replay audit remains
required before promotion.
