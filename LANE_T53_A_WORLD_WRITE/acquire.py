#!/usr/bin/env python3
"""Extract the ten selected U1537 VSM files from a locally supplied Zenodo archive.

The command is deliberately offline.  A registrar retrieves the public record metadata and
archive, then supplies both paths.  This program verifies the archive hash, the record DOI,
the CC-BY-4.0 licence, exact membership, and every selected file's repository MD5 before
writing raw/.  It never silently updates a source.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import zipfile


RECORD_ID = "14564186"
DOI = "10.5281/zenodo.14564186"
LICENSE = "cc-by-4.0"
ARCHIVE_SHA256 = "202a932bcf30af20087c21277c9d19f2b7217a69d849a198db9016f170f330ab"
SUFFIXES = (
    "Step 1 Hysteresis Measurement.csv",
    "Step 2 Remanence Curves.csv",
)


def digest(path: Path, algorithm: str = "sha256") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("raw"))
    args = parser.parse_args()

    if digest(args.archive) != ARCHIVE_SHA256:
        raise SystemExit("REFUSED: source archive SHA-256 differs from the frozen acquisition")

    metadata_bytes = args.metadata.read_bytes()
    metadata = json.loads(metadata_bytes)
    if str(metadata.get("id")) != RECORD_ID or metadata.get("metadata", {}).get("doi") != DOI:
        raise SystemExit("REFUSED: metadata is not the frozen Zenodo record")
    licence = metadata.get("metadata", {}).get("license", {}).get("id")
    if licence != LICENSE:
        raise SystemExit("REFUSED: source licence differs from CC-BY-4.0")

    repository_files = {item["key"]: item for item in metadata.get("files", [])}
    selected = sorted(name for name in repository_files if name.endswith(SUFFIXES))
    if len(selected) != 10:
        raise SystemExit(f"REFUSED: expected exactly 10 selected VSM files, found {len(selected)}")

    args.output.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.archive) as archive:
        if sorted(archive.namelist()) != sorted(repository_files):
            raise SystemExit("REFUSED: archive membership differs from record metadata")
        inventory = []
        for name in selected:
            payload = archive.read(name)
            expected = repository_files[name]["checksum"]
            if not expected.startswith("md5:"):
                raise SystemExit(f"REFUSED: unsupported repository checksum for {name}")
            got_md5 = hashlib.md5(payload).hexdigest()  # repository-supplied custody checksum
            if got_md5 != expected.removeprefix("md5:"):
                raise SystemExit(f"REFUSED: repository checksum mismatch for {name}")
            destination = args.output / name
            destination.write_bytes(payload)
            inventory.append(
                {
                    "file": name,
                    "bytes": len(payload),
                    "md5_repository": got_md5,
                    "sha256": hashlib.sha256(payload).hexdigest(),
                }
            )

    acquisition = {
        "record_id": RECORD_ID,
        "doi": DOI,
        "related_publication_doi": "10.1029/2025PA005360",
        "license": LICENSE,
        "archive_sha256": ARCHIVE_SHA256,
        "metadata_sha256": hashlib.sha256(metadata_bytes).hexdigest(),
        "selected_files": inventory,
    }
    (args.output / "SOURCE.json").write_text(
        json.dumps(acquisition, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(acquisition, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
