#!/usr/bin/env python3
"""Evaluate one WAC_GAMMA_FLOW_V001 envelope and print its certificate."""

from __future__ import annotations

import argparse
import sys

from gamma_flow import GammaFlowRefusal, certificate_json, load_gamma_flow


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate a Repair3-bounded, origin-neutral gamma-flow envelope. "
            "Synthetic input can never authorize scientific proof output."
        )
    )
    parser.add_argument("manifest", help="path to WAC_GAMMA_FLOW_V001 JSON")
    arguments = parser.parse_args()
    try:
        result = load_gamma_flow(arguments.manifest)
    except GammaFlowRefusal as exc:
        print(str(exc), file=sys.stderr)
        return 2
    sys.stdout.write(certificate_json(result))
    return 0 if result.certificate()["accepted"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
