#!/usr/bin/env python3
"""Validate a proof-frontier manifest and emit its non-authoritative certificate."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from proof_frontier import ProofFrontierRefusal, certificate_json, load_proof_frontier


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    arguments = parser.parse_args()
    try:
        frontier = load_proof_frontier(arguments.manifest)
    except ProofFrontierRefusal as exc:
        print(str(exc), file=sys.stderr)
        return 2
    print(certificate_json(frontier), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

