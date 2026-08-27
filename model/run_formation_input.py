#!/usr/bin/env python3
"""Validate a generic V002 formation input and optional execution envelope."""

from __future__ import annotations

import argparse
import sys

from formation_input import (
    FormationRefusal,
    attach_formation_execution,
    certificate_json,
    load_formation_input,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Validate custody and frozen generic predicates without authorizing "
            "a scientific formation, universal, gravity, or completion result."
        )
    )
    parser.add_argument("manifest", help="path to WAC_FORMATION_BUNDLE_V002 input")
    parser.add_argument(
        "--execution",
        help="optional WAC_FORMATION_EXECUTION_V002 envelope to attach after input lock",
    )
    arguments = parser.parse_args()
    try:
        formation = load_formation_input(arguments.manifest)
        result = (
            attach_formation_execution(formation, arguments.execution)
            if arguments.execution
            else formation
        )
    except FormationRefusal as exc:
        print(str(exc), file=sys.stderr)
        return 2
    sys.stdout.write(certificate_json(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
