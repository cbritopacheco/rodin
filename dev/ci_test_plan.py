#!/usr/bin/env python3
"""Plan hierarchical jobs from current test registration directories.

This discovers scheduling groups only. The configured CMake target manifest is
always authoritative for executable names and enabled backend conditions.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import itertools
import json
import os
from pathlib import Path
import re

from ci_tests import COMPONENTS, ROOT, display_group, execution_group


@functools.lru_cache(maxsize=8)
def groups(suite, root=ROOT):
    result = {}
    for cmake in sorted((root / "tests" / suite).rglob("CMakeLists.txt")):
        source = cmake.read_text()
        if not re.search(r"\b(add_executable|add_test)\s*\(", source):
            continue
        if re.search(r"(?<!rodin_)\bgtest_discover_tests\s*\(", source):
            raise ValueError("Use rodin_discover_tests so target profiles are recorded: " + str(cmake))
        owner = str(cmake.parent.relative_to(root / "tests"))
        scheduled = execution_group(owner)
        slow = bool(re.search(r"\bslow\b", source, re.IGNORECASE))
        result[scheduled] = result.get(scheduled, False) or slow
        for target, component in COMPONENTS.items():
            if target in source:
                scheduled = execution_group(owner + "/" + component)
                result[scheduled] = result.get(scheduled, False) or slow
    return result



def native_convergence_groups(root=ROOT):
    dedicated = {"Diffusion", "GeometryApproximation", "Helmholtz", "LinearElasticity",
                 "NonlinearPoisson", "P0Projection", "ReactionDiffusion", "Stokes"}
    return [group for group in groups("convergence", root)
            if "/PETSc" not in group and group not in {
                "convergence/isoparametric/" + component for component in dedicated}]


def entry(profile, group, build_type="Release", multithreaded="ON",
          asan="OFF", ubsan="OFF", baseline=False, backend=None):
    backend = backend or {"unit": "optional", "manufactured": "default"}.get(profile, "baseline")
    coordinates = [profile, group, build_type, multithreaded, asan, ubsan, str(baseline), backend]
    identity = hashlib.sha256(json.dumps(coordinates).encode()).hexdigest()[:20]
    path = display_group(group)
    if baseline:
        path = "Coverage / Baseline / " + ("Optional / " if profile == "unit" else "") + path
    elif profile == "slow":
        path += " / Slow"
    elif profile == "coverage":
        path = "Coverage / " + path
    if backend == "baseline" and profile in ["unit", "manufactured"] and not baseline:
        path += " / Backends / Baseline"
    return {"profile": profile, "backend": backend, "group": group, "build_type": build_type,
            "multithreaded": multithreaded, "asan": asan, "ubsan": ubsan,
            "baseline": baseline, "name": path, "id": identity}


def plan(root=ROOT):
    unit = groups("unit", root)
    unit["cmake/ThirdPartyIsolation"] = False
    manufactured = groups("manufactured", root)
    matrices = {name: [] for name in ["unit_1", "unit_2", "unit_3", "manufactured_1",
                                     "manufactured_2", "slow", "baseline"]}
    for suite, candidates, count in [("unit", unit, 3), ("manufactured", manufactured, 2)]:
        for index, group in enumerate(candidates):
            key = suite + "_" + str(index % count + 1)
            for build_type, mt, asan, ubsan in itertools.product(
                    ["Debug", "Release"], ["ON", "OFF"], ["ON", "OFF"], ["ON", "OFF"]):
                matrices[key].append(entry(suite, group, build_type, mt, asan, ubsan))
    # Preserve the former Coverage backend configuration in ordinary lanes.
    # Slow tests already execute in this configuration and are not rerun here.
    for suite, candidates, count in [("unit", unit, 3), ("manufactured", manufactured, 2)]:
        for index, group in enumerate(candidates):
            if group.startswith("cmake/"):
                continue
            for mt in ["ON", "OFF"]:
                matrices[suite + "_" + str(index % count + 1)].append(
                    entry(suite, group, multithreaded=mt, backend="baseline"))
    for group, slow in {**unit, **manufactured}.items():
        for mt in ["ON", "OFF"]:
            if slow:
                matrices["slow"].append(entry("slow", group, multithreaded=mt))
    # Optional-backend Unit used to compile manufactured TUs without running
    # them. Retain their zero-hit instantiations once per instrumented profile.
    for group in manufactured:
        for build_type, mt in itertools.product(["Debug", "Release"], ["ON", "OFF"]):
            matrices["baseline"].append(entry("unit", group, build_type, mt, baseline=True))
    # The default all build also compiled production libraries even when no test
    # links one. Preserve that compilation inventory in bounded baseline jobs.
    for profile, build_types in [("unit", ["Debug", "Release"]), ("coverage", ["Release"]),
                                 ("manufactured", ["Debug", "Release"])]:
        for build_type, mt in itertools.product(build_types, ["ON", "OFF"]):
            matrices["baseline"].append(entry(profile, "src/Rodin", build_type, mt, baseline=True))
    for key, rows in matrices.items():
        if not rows or len(rows) > 256:
            raise ValueError(f"{key} has {len(rows)} jobs; adjust the named family bands (limit 256)")
    return matrices


def targeted_plan(native):
    """Read concrete execution matrices so added partitions cannot escape the gate."""
    import yaml
    workflow = yaml.safe_load((ROOT / ".github/workflows/Tests.yml").read_text())
    result = []
    for job, definition in workflow["jobs"].items():
        if not any(step.get("id") == "coverage_record" for step in definition.get("steps", [])):
            continue
        matrix = definition["strategy"]["matrix"]
        rows = native if job == "ConvergenceBaseline" else matrix["include"]
        for row in rows:
            result.append({"job": job, "matrix": row})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("ci-test-plan.json"))
    args = parser.parse_args()
    matrices = plan()
    from ci_builds import plan as build_plan
    import yaml
    builds = build_plan(yaml.safe_load((ROOT / ".github/workflows/Tests.yml").read_text()))
    args.output.write_text(json.dumps(matrices, indent=2) + "\n")
    args.output.with_name("ci-build-plan.json").write_text(json.dumps(builds, indent=2) + "\n")
    native = [{"multithreaded": mt, "suite": "baseline", "group": group,
               "path": display_group(group)} for group in native_convergence_groups() for mt in ["ON", "OFF"]]
    args.output.with_name("ci-convergence-plan.json").write_text(json.dumps(native, indent=2) + "\n")
    args.output.with_name("ci-coverage-targets.json").write_text(
        json.dumps(targeted_plan(native), indent=2) + "\n")
    if len(native) > 256:
        raise ValueError("Native convergence exceeds the matrix limit; split by refinement family")
    if os.environ.get("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a") as stream:
            for name, rows in matrices.items():
                stream.write(name + "=" + json.dumps({"include": rows}, separators=(",", ":")) + "\n")
            stream.write("native_convergence=" + json.dumps({"include": native}, separators=(",", ":")) + "\n")
    print("Dedicated compile contracts:", {name: len(rows) for name, rows in builds.items()})
    print("Planned jobs:", {name: len(rows) for name, rows in matrices.items()})


if __name__ == "__main__":
    main()
