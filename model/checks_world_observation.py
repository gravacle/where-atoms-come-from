#!/usr/bin/env python3
"""Fast synthetic checks for the closed world-observation input contract."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import tempfile

import world_observation as contract
from world_observation import (
    CORE_COLUMNS,
    ObservationRefusal,
    certificate_json,
    load_world_observation,
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_bundle(root: Path) -> Path:
    protocol = root / "PROTOCOL.txt"
    source = root / "raw.dat"
    data = root / "observations.csv"
    protocol.write_text(
        "# WAC WORLD OBSERVATION PROTOCOL V001\n\nsynthetic protocol fixture\n",
        encoding="utf-8",
    )
    columns = list(CORE_COLUMNS) + ["writer_level"]
    with data.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(columns)
        writer.writerow(["r1", "sample-a", "write-plus", "POST_WRITE",
                         "RECORD_READOUT", "0", "0", "2", "1"])
        writer.writerow(["r2", "sample-a", "write-minus", "POST_WRITE",
                         "RECORD_READOUT", "1", "0", "-2", "-1"])
    source.write_bytes(data.read_bytes())
    manifest = {
        "schema": "WAC_WORLD_OBSERVATION_V001",
        "surface_id": "synthetic-surface",
        "sample_ids": ["sample-a"],
        "run_id": "synthetic-run",
        "instrument_id": "synthetic-instrument",
        "source_uri": "urn:wac:synthetic-fixture",
        "evidence_class": "SYNTHETIC_TEST_ONLY",
        "observation_scope": "WRITE_POST_ONLY",
        "protocol_timing": "PROSPECTIVE_PRE_OUTCOME",
        "protocol_frozen_at_utc": "2026-08-22T00:00:00Z",
        "protocol_path": protocol.name,
        "protocol_sha256": digest(protocol),
        "data_file": data.name,
        "data_sha256": digest(data),
        "columns": columns,
        "units": {"time": "s", "coordinate": "1", "value": "1"},
        "extra_columns": [{"name": "writer_level", "type": "float", "unit": "1"}],
        "controls": ["OPPOSITE_WRITE"],
        "normalization": {
            "adapter_id": "DIRECT_CSV_COPY_TEST_ONLY_V001",
            "adapter_sha256": digest(Path(contract.__file__).resolve()),
            "source_paths": [source.name],
        },
        "source_artifacts": [
            {"path": source.name, "sha256": digest(source), "media_type": "text/csv"}
        ],
    }
    path = root / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def mutate_manifest(path: Path, mutation) -> None:
    value = json.loads(path.read_text(encoding="utf-8"))
    mutation(value)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sync_normalized_source(path: Path) -> None:
    value = json.loads(path.read_text(encoding="utf-8"))
    root = path.parent
    data = root / value["data_file"]
    source = root / value["normalization"]["source_paths"][0]
    source.write_bytes(data.read_bytes())
    value["data_sha256"] = digest(data)
    for artifact in value["source_artifacts"]:
        if artifact["path"] == source.name:
            artifact["sha256"] = digest(source)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def refuses(path: Path, contains: str) -> None:
    try:
        load_world_observation(path)
    except ObservationRefusal as exc:
        assert contains in str(exc), (contains, str(exc))
    else:
        raise AssertionError(f"expected refusal containing {contains!r}")


def main() -> int:
    checks = 0
    with tempfile.TemporaryDirectory(prefix="wac-world-observation-") as directory:
        root = Path(directory)

        manifest = write_bundle(root)
        observation = load_world_observation(manifest)
        certificate = observation.certificate()
        assert certificate["coverage"]["row_count"] == 2
        assert certificate["scope_classification"] == "CONFIGURATION_EVIDENCE_ONLY"
        assert certificate["scientific_verdict"] == "NONE_NOT_SCORED"
        assert not certificate["record_formation_proof_authorized"]
        assert certificate_json(observation) == certificate_json(load_world_observation(manifest))
        checks += 1

        manifest = write_bundle(root)
        mutate_manifest(manifest, lambda value: value.__setitem__("unexpected", True))
        refuses(manifest, "manifest key closure")
        checks += 1

        manifest = write_bundle(root)
        mutate_manifest(manifest, lambda value: value.__setitem__("data_sha256", "0" * 64))
        refuses(manifest, "data_sha256 mismatch")
        checks += 1

        manifest = write_bundle(root)
        mutate_manifest(manifest, lambda value: value.__setitem__("protocol_path", "../escape"))
        refuses(manifest, "safe relative path")
        checks += 1

        manifest = write_bundle(root)
        mutate_manifest(
            manifest,
            lambda value: value["source_artifacts"][0].__setitem__("sha256", "0" * 64),
        )
        refuses(manifest, "source artifact hash mismatch")
        checks += 1

        manifest = write_bundle(root)
        data = root / "observations.csv"
        payload = data.read_text(encoding="utf-8").replace(",2,1\n", ",nan,1\n")
        data.write_text(payload, encoding="utf-8")
        sync_normalized_source(manifest)
        refuses(manifest, "value is non-finite")
        checks += 1

        manifest = write_bundle(root)
        data = root / "observations.csv"
        payload = data.read_text(encoding="utf-8").replace("writer_level", "wrong_header")
        data.write_text(payload, encoding="utf-8")
        sync_normalized_source(manifest)
        refuses(manifest, "CSV header differs")
        checks += 1

        manifest = write_bundle(root)
        mutate_manifest(manifest, lambda value: value.__setitem__("observation_scope", "ALL_PROVED"))
        refuses(manifest, "observation_scope is not registered")
        checks += 1

        manifest = write_bundle(root)
        mutate_manifest(manifest, lambda value: value.__setitem__("observation_scope", "FULL_WRITE_HOLD_READ"))
        refuses(manifest, "complete joined trajectory")
        checks += 1

        manifest = write_bundle(root)
        data = root / "observations.csv"
        lines = data.read_text(encoding="utf-8").splitlines()
        data.write_text("\n".join(lines + [lines[1]]) + "\n", encoding="utf-8")
        sync_normalized_source(manifest)
        refuses(manifest, "duplicate row_id")
        checks += 1

        manifest = write_bundle(root)
        mutate_manifest(manifest, lambda value: value.__setitem__("controls", ["NONE", "OPPOSITE_WRITE"]))
        refuses(manifest, "NONE cannot be combined")
        checks += 1

        manifest = write_bundle(root)
        mutate_manifest(manifest, lambda value: value.__setitem__("sample_ids", []))
        refuses(manifest, "sample_ids")
        checks += 1

        manifest = write_bundle(root)
        link = root / "outside-link"
        link.symlink_to(Path("/etc/hosts"))
        mutate_manifest(
            manifest,
            lambda value: value["source_artifacts"].__setitem__(
                0,
                {"path": link.name, "sha256": digest(Path("/etc/hosts")),
                 "media_type": "text/plain"},
            ),
        )
        refuses(manifest, "escapes the observation bundle through a symlink")
        checks += 1

        manifest = write_bundle(root)
        mutate_manifest(
            manifest,
            lambda value: value.__setitem__("protocol_frozen_at_utc", "2026-99-99T00:00:00Z"),
        )
        refuses(manifest, "not a real UTC date")
        checks += 1

        manifest = write_bundle(root)
        mutate_manifest(manifest, lambda value: value.__setitem__("source_uri", "not-a-uri"))
        refuses(manifest, "source_uri must use")
        checks += 1

        manifest = write_bundle(root)
        data = root / "observations.csv"
        lines = data.read_text(encoding="utf-8").splitlines()
        later_same_event = lines[1].replace("r1,", "r3,").replace(",0,0,2,", ",-1,0,2,")
        data.write_text("\n".join(lines + [later_same_event]) + "\n", encoding="utf-8")
        sync_normalized_source(manifest)
        refuses(manifest, "time decreases within record/event")
        checks += 1

        manifest = write_bundle(root)
        text = manifest.read_text(encoding="utf-8")
        manifest.write_text(text.replace("{\n", "{\n  \"schema\": \"DUPLICATE\",\n", 1),
                            encoding="utf-8")
        refuses(manifest, "duplicate JSON member")
        checks += 1

        manifest = write_bundle(root)
        data = root / "observations.csv"
        data.write_text(data.read_text(encoding="utf-8") + '"unterminated\n', encoding="utf-8")
        sync_normalized_source(manifest)
        refuses(manifest, "CSV is unreadable")
        checks += 1

        manifest = write_bundle(root)
        data = root / "observations.csv"
        header = data.read_text(encoding="utf-8").splitlines()[0]
        fragmented = [
            header,
            "f1,sample-a,trial-before,BEFORE,RECORD_READOUT,0,0,0,0",
            "f2,sample-a,trial-write,WRITE,WRITER_CONTROL,1,0,1,1",
            "f3,sample-a,trial-hold,HOLD,ENVIRONMENT,2,0,1,1",
            "f4,sample-a,trial-read,READ,RECORD_READOUT,3,0,1,1",
        ]
        data.write_text("\n".join(fragmented) + "\n", encoding="utf-8")
        sync_normalized_source(manifest)
        mutate_manifest(
            manifest, lambda value: value.__setitem__("observation_scope", "FULL_WRITE_HOLD_READ")
        )
        refuses(manifest, "complete joined trajectory")
        checks += 1

        manifest = write_bundle(root)
        data = root / "observations.csv"
        data.write_text(data.read_text(encoding="utf-8").replace(",2,1\n", ",3,1\n"),
                        encoding="utf-8")
        mutate_manifest(manifest, lambda value: value.__setitem__("data_sha256", digest(data)))
        refuses(manifest, "not the registered adapter output")
        checks += 1

        manifest = write_bundle(root)
        linked_manifest = root / "manifest-link.json"
        linked_manifest.symlink_to(manifest.name)
        refuses(linked_manifest, "manifest must not be a symlink")
        checks += 1

        manifest = write_bundle(root)
        source = root / "raw.dat"
        alias = root / "raw-alias.dat"
        os.link(source, alias)
        mutate_manifest(
            manifest,
            lambda value: value["source_artifacts"].append(
                {"path": alias.name, "sha256": digest(alias), "media_type": "text/csv"}
            ),
        )
        refuses(manifest, "duplicate physical source artifact")
        checks += 1

        manifest = write_bundle(root)
        protocol = root / "PROTOCOL.txt"
        mutate_manifest(
            manifest,
            lambda value: (
                value.__setitem__("data_file", protocol.name),
                value.__setitem__("data_sha256", digest(protocol)),
            ),
        )
        refuses(manifest, "must be distinct physical files")
        checks += 1

        manifest = write_bundle(root)
        protocol = root / "PROTOCOL.txt"
        protocol.write_text("# unregistered protocol\n", encoding="utf-8")
        mutate_manifest(manifest, lambda value: value.__setitem__("protocol_sha256", digest(protocol)))
        refuses(manifest, "protocol document marker")
        checks += 1

        real_bundle = root / "real-bundle"
        real_bundle.mkdir()
        real_manifest = write_bundle(real_bundle)
        alias_bundle = root / "alias-bundle"
        alias_bundle.symlink_to(real_bundle, target_is_directory=True)
        refuses(alias_bundle / real_manifest.name, "bundle directory must not be a symlink")
        checks += 1

        manifest = write_bundle(root)
        mutate_manifest(
            manifest,
            lambda value: value.__setitem__("sample_ids", ["sample-a", "ghost-sample"]),
        )
        refuses(manifest, "every declared sample")
        checks += 1

        manifest = write_bundle(root)
        data = root / "observations.csv"
        data.write_text(data.read_text(encoding="utf-8").replace("POST_WRITE", "RATE"),
                        encoding="utf-8")
        sync_normalized_source(manifest)
        mutate_manifest(
            manifest,
            lambda value: (
                value.__setitem__("sample_ids", ["sample-a", "sample-b"]),
                value.__setitem__("observation_scope", "RATE_STABILITY_ONLY"),
            ),
        )
        refuses(manifest, "every declared sample")
        checks += 1

        manifest = write_bundle(root)
        data = root / "observations.csv"
        lines = data.read_text(encoding="utf-8").splitlines()
        rate_a = lines[1].replace("POST_WRITE", "RATE")
        post_b = lines[2].replace("sample-a", "sample-b")
        data.write_text("\n".join([lines[0], rate_a, post_b]) + "\n", encoding="utf-8")
        sync_normalized_source(manifest)
        mutate_manifest(
            manifest,
            lambda value: (
                value.__setitem__("sample_ids", ["sample-a", "sample-b"]),
                value.__setitem__("observation_scope", "RATE_STABILITY_ONLY"),
            ),
        )
        refuses(manifest, "RATE row for every sample")
        checks += 1

    print(f"WORLD_OBSERVATION_CHECKS: {checks}/28 PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
