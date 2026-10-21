#!/usr/bin/env python3
"""Record coverage producers and enforce complete current-commit collection.

The compilation-inventory check preserves zero-hit template/header records from
old all-target builds. The final line union contains production sources only;
execution of tests contributes hits, never test-source denominator entries.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path


def record(args):
    build = args.build
    definition = json.loads(os.environ["CI_TEST_DEFINITION"])
    compiled = json.loads((build / "ci-built-targets.json").read_text())
    manifests = list(build.glob("RodinTestTargets-*.tsv"))
    with manifests[0].open() as stream:
        configured = list(csv.DictReader(stream, delimiter="\t"))
    report = build / "RodinCoverage.info"
    document = {"definition": definition, "commit": os.environ["GITHUB_SHA"],
                "configured": configured, "compiled": compiled,
                "coverage_expected": os.environ.get("RODIN_COVERAGE") == "ON" and bool(compiled),
                "uploaded": os.environ.get("COVERAGE_UPLOADED") == "success",
                "build": os.environ.get("BUILD_OUTCOME"),
                "tests": os.environ.get("TEST_OUTCOME"),
                "capture": os.environ.get("CAPTURE_OUTCOME"),
                "report_exists": report.is_file(),
                "report_sha256": hashlib.sha256(report.read_bytes()).hexdigest() if report.is_file() else None}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "receipt.json").write_text(json.dumps(document, indent=2) + "\n")
    for path in [build / "ci-selection.json", build / "ci-invocations.json", build / "ci-built-targets.json", report,
                 build / "coverage-audit.json"]:
        if path.is_file():
            (args.output / path.name).write_bytes(path.read_bytes())


def compilation_context(definition):
    profile = "baseline" if definition["profile"] in ["slow", "coverage"] else definition["profile"]
    return profile, definition["build_type"], definition["multithreaded"], definition["asan"], definition["ubsan"]


def verify_receipts(plan, folder, commit):
    expected = {entry["id"]: entry for entries in plan.values() for entry in entries}
    received = {}
    inventory = {}
    built = {}
    for path in folder.rglob("receipt.json"):
        receipt = json.loads(path.read_text())
        definition = receipt["definition"]
        identity = definition["id"]
        if identity not in expected or expected[identity] != definition or identity in received:
            raise ValueError("Unexpected or duplicate producer: " + identity)
        if receipt["commit"] != commit:
            raise ValueError("Stale producer: " + identity)
        expected_tests = "skipped" if definition["baseline"] else "success"
        if receipt["build"] != "success" or receipt["tests"] != expected_tests:
            raise ValueError("Producer failed: " + identity)
        if not definition["baseline"]:
            selection_path = path.parent / "ci-selection.json"
            if not selection_path.is_file():
                raise ValueError("Missing invocation manifest: " + identity)
            selection = json.loads(selection_path.read_text())
            if selection["group"] != definition["group"] or selection["profile"] != definition["profile"]:
                raise ValueError("Wrong invocation manifest: " + identity)
        if receipt["coverage_expected"]:
            report = path.parent / "RodinCoverage.info"
            if not (receipt["uploaded"] and receipt["report_exists"] and receipt["capture"] == "success"
                    and report.is_file() and hashlib.sha256(report.read_bytes()).hexdigest() == receipt["report_sha256"]):
                raise ValueError("Coverage upload incomplete: " + identity)
            context = compilation_context(definition)
            built.setdefault(context, set()).update(t["target"] for t in receipt["compiled"])
            if definition["baseline"] and definition["group"] == "src/Rodin":
                inventory[context] = {t["target"] for t in receipt["configured"]}
        received[identity] = receipt
    missing = expected.keys() - received.keys()
    if missing:
        raise ValueError("Missing producer records: " + ", ".join(sorted(missing)))
    if len(inventory) != 6:
        raise ValueError("Expected four optional-backend and two baseline compilation inventories")
    for context, required in inventory.items():
        absent = required - built.get(context, set())
        if absent:
            raise ValueError("Missing initial-coverage compilation targets for " + str(context) + ": " + str(sorted(absent)))
    return received


def production_lines(paths):
    lines = {}
    for path in paths:
        source = None
        for text in path.read_text().splitlines():
            if text.startswith("SF:"):
                raw = text[3:].replace("\\", "/")
                source = "src/Rodin/" + raw.split("src/Rodin/", 1)[1] if "src/Rodin/" in raw else None
                if source is None or source.startswith("src/Rodin/Test/"):
                    raise ValueError("Non-production source in coverage report: " + raw)
            elif text == "end_of_record":
                source = None
            elif text.startswith("DA:"):
                if source is None:
                    raise ValueError("Line record without a production source")
                line, count, *_ = text[3:].split(",")
                key = (source, int(line))
                lines[key] = lines.get(key, False) or int(count) > 0
    if not lines:
        raise ValueError("No production line coverage received")
    return lines


def finalize(args):
    needs = json.loads(os.environ["CI_JOB_RESULTS"])
    failures = {name: job["result"] for name, job in needs.items() if job["result"] != "success"}
    if failures:
        raise ValueError("Required test/coverage jobs did not succeed: " + str(failures))
    plan = json.loads(args.plan.read_text())
    verify_receipts(plan, args.artifacts, os.environ["GITHUB_SHA"])
    reports = list(args.artifacts.rglob("RodinCoverage.info"))
    # Legacy targeted producers retain their original artifacts and uploads.
    # Successful job results alone cannot substitute for a missing LCOV artifact.
    for prefix, count in [("coverage-audit-complex-h1-", 2),
                          ("coverage-audit-petsc-", 4)]:
        found = [p for p in args.artifacts.iterdir() if p.name.startswith(prefix)]
        if len(found) != count or any(not (p / "RodinCoverage.info").is_file() for p in found):
            raise ValueError("Missing targeted producer artifacts: " + prefix)
    lines = production_lines(reports)
    hit = sum(lines.values())
    rate = 100 * hit / len(lines)
    print(f"Complete production coverage: {hit}/{len(lines)} lines ({rate:.4f}%)")
    args.output.write_text(json.dumps({"hit": hit, "lines": len(lines), "percent": rate,
                                      "reports": len(reports)}, indent=2) + "\n")



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=["record", "finalize", "threshold"])
    parser.add_argument("--build", type=Path, default=Path("build"))
    parser.add_argument("--output", type=Path, default=Path("coverage-summary.json"))
    parser.add_argument("--plan", type=Path, default=Path("ci-test-plan.json"))
    parser.add_argument("--artifacts", type=Path, default=Path("ci-artifacts"))
    args = parser.parse_args()
    if args.operation == "record":
        record(args)
    elif args.operation == "finalize":
        finalize(args)
    else:
        summary = json.loads(args.output.read_text())
        if 100 * summary["hit"] < 95 * summary["lines"]:
            raise ValueError("Production coverage is below the required 95%")


if __name__ == "__main__":
    try:
        main()
    except ValueError as error:
        raise SystemExit(str(error)) from error
