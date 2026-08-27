#!/usr/bin/env python3
"""Validate a world-observation bundle and emit its measurement-only certificate."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from world_observation import ObservationRefusal, certificate_json, load_world_observation


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        rendered = certificate_json(load_world_observation(args.manifest))
    except ObservationRefusal as exc:
        print(str(exc), file=sys.stderr)
        return 2
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
