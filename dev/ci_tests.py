#!/usr/bin/env python3
"""Build and run a configured test group without discovering unbuilt siblings.

CMake exports target ownership, including optional and newly registered targets.
The runner retains CMake's generated commands and properties in a scoped CTest
view. GoogleTest discovery runs only for built targets, outside coverage. Tests
still execute through CTest, with their original working directories and ranks.
"""

from __future__ import annotations

import argparse
import csv
import functools
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
# Keep the fallback expression below Linux's 128 KiB single-argument limit.
MAX_SELECTION_EXPRESSION_BYTES = 100_000
# Component files within large registration directories receive their own jobs.
COMPONENTS = {
    "RodinVariationalQuadratureRuleSpecializationTest": "QuadratureRuleSpecialization",
    "RodinVariationalMatrixQuadratureTest": "MatrixQuadrature",
    "RodinVariationalQuadratureExpressionTest": "QuadratureExpression",
    "RodinVariationalH1MassQuadratureTest": "MassQuadrature",
    "RodinVariationalH1ShapeFunctionTest": "ShapeFunction",
    "RodinManufacturedLinearElasticity3D": "LinearElasticity3D",
    "RodinManufacturedMPIH1PoissonTest": "Poisson",
    "RodinManufacturedMPIH1LinearElasticityTest": "LinearElasticity",
}
NATIVE_RESTRICTION = r"^RodinMPIH1SubMeshRestrictionTest_.*_np(1|2|3|4|8)$"
UNIT_SLOW_GEOMETRY = (
    r"^RodinMPIAssemblyGeometrySlowTest_np(1|2|3|4|8)_"
    r"(Segment|Triangle|Quadrilateral|Tetrahedron|Pyramid|Hexahedron|Wedge)$"
)


def group_for(target):
    """Return the exact owner group; a parent does not absorb child groups."""
    owner = target["owner"].removeprefix("tests/")
    component = COMPONENTS.get(target["target"])
    return owner + "/" + component if component else owner


def execution_group(owner):
    """Batch short sibling suites at their common parent; isolate long branches.

    Source ownership remains unchanged and supplies each test's detailed labels.
    This avoids repeating runner provisioning for every tiny directory and every
    sanitizer configuration while keeping expensive components independent.
    """
    if owner.startswith("unit/Rodin/"):
        branches = ["Adaptation", "Geometry", "Models/Solid", "MPI/Assembly", "MPI/Variational",
                    "MPI", "PETSc", "Variational/H1/MassQuadrature", "Variational/H1/ShapeFunction",
                    "Variational/H1", "Variational/QuadratureRuleSpecialization",
                    "Variational/MatrixQuadrature", "Variational/QuadratureExpression", "Variational"]
        for branch in branches:
            prefix = "unit/Rodin/" + branch
            if owner == prefix or owner.startswith(prefix + "/"):
                return prefix
        return "unit/Rodin"
    if owner.startswith("manufactured/Rodin/"):
        branches = ["MPI/H1/Poisson", "MPI/H1/LinearElasticity", "H1", "P1/LinearElasticity3D", "P1"]
        for branch in branches:
            prefix = "manufactured/Rodin/" + branch
            if owner == prefix or owner.startswith(prefix + "/"):
                return prefix
        return "manufactured/Rodin"
    return owner


def display_group(group):
    parts = group.split("/")
    parts[0] = {"unit": "Unit", "manufactured": "Manufactured",
                "convergence": "Convergence", "cmake": "CMake",
                "src": "Libraries"}.get(parts[0], parts[0])
    names = {"h": "H", "p": "P", "hp": "HP", "isoparametric": "Isoparametric",
             "petsc": "PETSc", "mpi": "MPI"}
    return " / ".join(names.get(part, part[:1].upper() + part[1:]) for part in parts)


def targets(build):
    manifests = sorted(build.glob("RodinTestTargets-*.tsv"))
    if len(manifests) != 1:
        raise ValueError("Use a single-configuration build with one RodinTestTargets manifest")
    with manifests[0].open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


@functools.lru_cache(maxsize=4)
def configured_profiles(build):
    configured = targets(build)
    path_index = {row["executable"]: row for row in configured if row["executable"]}
    requirements = {row["target"]: set(row.get("profiles", "").split(";")) for row in configured}
    # Explicit MPI/convergence registrations are already generated at configure
    # time, so their commands and labels can be inspected without executing them.
    for source in (build / "tests").rglob("CTestTestfile.cmake"):
        if not (ROOT / source.parent.relative_to(build)).is_dir():
            continue
        parsed = list(commands(source.read_text()))
        labels = {}
        for name, args, _ in parsed:
            if name == "set_tests_properties":
                boundary = args.index("PROPERTIES")
                props = dict(zip(args[boundary + 1::2], args[boundary + 2::2]))
                for test in args[:boundary]:
                    labels[test] = props.get("LABELS", "").split(";")
        for name, args, _ in parsed:
            if name == "add_test":
                target = executable_target(args[1:], path_index)
                if not target and not args[0].startswith("CMake.ThirdPartyIsolation."):
                    raise ValueError("Unmapped explicit CTest invocation: " + args[0])
                if target:
                    test = {"name": args[0], "properties": [
                        {"name": "LABELS", "value": labels.get(args[0], [])}]}
                    requirements[target["target"]].update(
                        p for p in ["unit", "manufactured", "slow", "coverage", "convergence-baseline"] if selected(p, test))
    if (ROOT / "dev/ci_test_plan.py").is_file():
        # Planning reads source directories, while target names/conditions come
        # from CMake. Reject a future registration that would escape either plan.
        from ci_test_plan import groups, native_convergence_groups
        planned = {**groups("unit"), **groups("manufactured")}
        for target in configured:
            if (target["kind"] == "test" and "convergence-baseline" in requirements[target["target"]]
                    and "/PETSc" not in group_for(target)
                    and group_for(target) not in native_convergence_groups()):
                raise ValueError("Convergence group is missing from the execution plan: " + group_for(target))
            if target["kind"] == "test" and target["owner"].startswith(("tests/unit/", "tests/manufactured/")):
                group = execution_group(group_for(target))
                if group not in planned or ("slow" in requirements[target["target"]] and not planned[group]):
                    raise ValueError("Configured test group is missing from the execution plan: " + group)
    return requirements


def selected_targets(build, group, exclude="", profile="", baseline=False):
    configured = targets(build)
    requirements = configured_profiles(build)
    result = []
    for target in configured:
        if group == "src/Rodin":
            keep = target["kind"] == "library"
            if (baseline and target["kind"] == "test" and
                    (profile != "unit" or target["owner"].startswith("tests/unit/"))):
                # Retain test TUs from the old all build which no execution lane
                # compiles. Optional Unit does not have the baseline Slow lane.
                profiles = (["unit"] if profile == "unit" else
                            ["manufactured"] if profile == "manufactured" else ["unit", "manufactured", "slow"])
                keep = not any(p in requirements[target["target"]] for p in profiles)
        else:
            keep = target["kind"] == "test" and (execution_group(group_for(target)) == group or group_for(target) == group or
                    (baseline and group == "convergence" and target["owner"].startswith("tests/convergence/")))
            if profile and not baseline:
                keep = keep and profile in requirements[target["target"]]
        if keep and not (exclude and re.search(exclude, target["target"])):
            result.append(target)
    return result


def quote(value):
    """Quote a literal CMake argument without changing semicolons or newlines."""
    marker = "="
    while "]" + marker + "]" in str(value):
        marker += "="
    return "[" + marker + "[" + str(value) + "]" + marker + "]"


def commands(text):
    """Read generated CMake commands, preserving their complete original text.

    Only token identification is needed; copied commands are never reserialized.
    Quoted/bracket arguments, escaped quotes and comments may contain parentheses.
    """
    index = 0
    length = len(text)
    bracket_pattern = re.compile(r"\[(=*)\[")
    comment_pattern = re.compile(r"#\[(=*)\[")
    command_pattern = re.compile(r"[A-Za-z_][A-Za-z_0-9]*\s*\(")

    def token():
        nonlocal index
        char = text[index]
        bracket = bracket_pattern.match(text, index)
        if bracket:
            end = text.find("]" + bracket[1] + "]", index + len(bracket[0]))
            if end < 0:
                raise ValueError("Unterminated CMake bracket argument")
            value = text[index + len(bracket[0]):end]
            index = end + len(bracket[1]) + 2
            return value
        if char == '"':
            index += 1
            value = ""
            while index < length:
                char = text[index]
                index += 1
                if char == '"':
                    return value
                if char == "\\":
                    if index >= length:
                        break
                    char = text[index]
                    index += 1
                    char = {"n": "\n", "r": "\r", "t": "\t"}.get(char, char)
                value += char
            raise ValueError("Unterminated CMake quoted argument")
        start = index
        while index < length and not text[index].isspace() and text[index] not in "()":
            if text[index] == "\\":
                index += 1
            index += 1
        return text[start:index]

    while index < length:
        if text[index].isspace():
            index += 1
            continue
        if text[index] == "#":
            bracket = comment_pattern.match(text, index)
            if bracket:
                end = text.find("]" + bracket[1] + "]", index + len(bracket[0]))
                if end < 0:
                    raise ValueError("Unterminated CMake bracket comment")
                index = end + len(bracket[1]) + 2
            else:
                index = text.find("\n", index)
                if index < 0:
                    return
            continue
        start = index
        match = command_pattern.match(text, index)
        if not match:
            raise ValueError("Unexpected generated CMake syntax: " + text[index:index + 80])
        name = match[0].split("(")[0].strip().lower()
        index += len(match[0])
        depth = 1
        arguments = []
        while index < length and depth:
            if text[index].isspace():
                index += 1
            elif text[index] == ")":
                depth -= 1
                index += 1
            elif text[index] == "(":
                depth += 1
                arguments.append("(")
                index += 1
            else:
                arguments.append(token())
        if depth:
            raise ValueError("Unterminated CMake command")
        yield name, arguments, text[start:index]


def executable_target(command, configured):
    paths = configured if isinstance(configured, dict) else {
        row["executable"]: row for row in configured if row["executable"]}
    for argument in command:
        if argument in paths:
            return paths[argument]
    return None


def scoped_tree(build, group, configured, profile=""):
    """Retain generated registrations only for this exact group."""
    selected = selected_targets(build, group, profile=profile)
    chosen = {row["target"] for row in selected}
    stage = build / "ci-tests"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir()
    directories = []
    for source in sorted((build / "tests").rglob("CTestTestfile.cmake")):
        # Only the configured source-tree directories belong to this build.
        relative = source.parent.relative_to(build)
        if not (ROOT / relative).is_dir():
            continue
        parsed = list(commands(source.read_text()))
        retained = {}
        includes = []
        for name, args, raw in parsed:
            if name == "add_test":
                target = executable_target(args[1:], configured)
                if target and target["target"] in chosen:
                    retained[args[0]] = raw
                elif (not target and group == "cmake/ThirdPartyIsolation"
                      and args[0].startswith("CMake.ThirdPartyIsolation.")):
                    retained[args[0]] = raw
            elif name == "include":
                path = Path(args[0])
                if any(re.fullmatch(re.escape(t) + r"\[\d+\]_include\.cmake", path.name)
                       for t in chosen):
                    includes.append(raw)
        properties = []
        working = set()
        for name, args, raw in parsed:
            if name == "set_tests_properties":
                boundary = args.index("PROPERTIES")
                names = args[:boundary]
                wanted = [n for n in names if n in retained]
                if wanted:
                    if wanted != names:
                        raise ValueError("Shared property block crosses test groups: " + str(source))
                    properties.append(raw)
                    if "WORKING_DIRECTORY" in args[boundary + 1::2]:
                        working.update(wanted)
        if retained or includes:
            destination = stage / relative
            destination.mkdir(parents=True, exist_ok=True)
            body = list(retained.values()) + properties + includes
            for test in retained.keys() - working:
                body.append("set_tests_properties(" + quote(test)
                            + " PROPERTIES WORKING_DIRECTORY " + quote(source.parent) + ")")
            (destination / "CTestTestfile.cmake").write_text("\n".join(body) + "\n")
            directories.append(relative)
    (stage / "CTestTestfile.cmake").write_text("\n".join(
        "subdirs(" + quote(d) + ")" for d in directories) + "\n")
    return stage


def discover(stage, build):
    environment = dict(os.environ, GCOV_PREFIX=str(build.parent / (build.name + "-coverage-discovery-excluded")))
    result = subprocess.run(["ctest", "--test-dir", str(stage), "--show-only=json-v1"],
                            check=True, capture_output=True, text=True, env=environment)
    report = json.loads(result.stdout)
    missing = [test["name"] for test in report["tests"] if "NOT_BUILT" in test["name"]]
    if missing:
        raise ValueError("Selected executable was not built: " + str(missing))
    return report


def properties(test):
    return {prop["name"]: prop["value"] for prop in test.get("properties", [])}


def prepare(build, group, profile=""):
    configured = targets(build)
    stage = scoped_tree(build, group, configured, profile)
    before = discover(stage, build)
    for test in before["tests"]:
        target = executable_target(test.get("command", []), configured)
        owner = Path(target["owner"] if target else "tests")
        path = stage / owner / "CTestTestfile.cmake"
        label = display_group(group_for(target) if target else group).replace(" / ", "/")
        hierarchy = ["/".join(label.split("/")[:i]) for i in range(1, len(label.split("/")) + 1)]
        labels = list(dict.fromkeys(properties(test).get("LABELS", []) + hierarchy))
        with path.open("a") as stream:
            stream.write("set_tests_properties(" + quote(test["name"])
                         + " PROPERTIES LABELS " + quote(";".join(labels)) + ")\n")
    after = discover(stage, build)
    # The scoped view must not alter commands, timeouts, MPI arguments or resources.
    def invariant(report):
        return sorted(json.dumps({"name": t["name"], "command": t.get("command"),
                                  "properties": {k: v for k, v in properties(t).items() if k != "LABELS"}},
                                 sort_keys=True) for t in report["tests"])
    if invariant(before) != invariant(after):
        raise ValueError("Hierarchy annotation changed an invocation")
    (build / "ci-invocations.json").write_text(json.dumps(after, indent=2))
    return stage, after


def selected(profile, test):
    labels = set(properties(test).get("LABELS", []))
    name = test["name"]
    if profile == "convergence-baseline":
        excluded = (r"^RodinConvergenceIsoparametric(P0Projection|GeometryApproximation|.*Diffusion|NonlinearPoisson|LinearElasticity|Stokes|Helmholtz)_"
                    r"|^RodinConvergenceHStokesStability_")
        return "convergence" in labels and not re.search(excluded, name)
    if profile == "unit":
        return ("unit" in labels and "slow" not in labels) or bool(re.fullmatch(UNIT_SLOW_GEOMETRY, name))
    if profile == "manufactured":
        return "manufactured" in labels and "slow" not in labels
    if profile == "slow":
        return "slow" in labels and not re.fullmatch(NATIVE_RESTRICTION, name)
    if profile == "coverage":
        return (("unit" in labels and "distributed" not in labels)
                or ("distributed" in labels and "slow" not in labels)
                or ("manufactured" in labels and "slow" not in labels and "distributed" not in labels))
    raise ValueError("Unknown execution profile: " + profile)


def run_group(build, group, profile, backend=None):
    stage, report = prepare(build, group, profile)
    chosen = [test for test in report["tests"] if selected(profile, test) and
              not (backend == "baseline" and profile == "unit" and
                   "slow" in properties(test).get("LABELS", []))]
    inventory = {"group": group, "profile": profile, "registered": len(report["tests"]),
                 "selected": chosen, "targets": selected_targets(build, group, profile=profile)}
    (build / "ci-selection.json").write_text(json.dumps(inventory, indent=2))
    if not chosen:
        # Backend conditions and slow-only binaries can legitimately leave an
        # empty runtime selection. Their compilation still supplies the baseline.
        print("No invocations for this group/profile; recorded explicit empty selection", flush=True)
        return 0
    version = subprocess.check_output(["ctest", "--version"], text=True)
    match = re.search(r"version (\d+)\.(\d+)", version)
    if match and tuple(map(int, match.groups())) >= (3, 29):
        selection_file = build / "ci-selected-tests.txt"
        selection_file.write_text("\n".join(t["name"] for t in chosen) + "\n")
        selection_args = ["--tests-from-file", str(selection_file)]
    else:
        expression = "^(" + "|".join(re.escape(t["name"]) for t in chosen) + ")$"
        if len(expression.encode()) > MAX_SELECTION_EXPRESSION_BYTES:
            raise ValueError("Large scoped selections require CTest 3.29 or newer")
        selection_args = ["-R", expression]
    # Keep CTest serial: distributed invocations retain their registered ranks.
    args = ["ctest", "--test-dir", str(stage)] + selection_args + [
        "--output-on-failure", "--no-tests=error", "-j", "1"]
    # Verify CTest selects the same invocations before running any tests.
    env = dict(os.environ, GCOV_PREFIX=str(build.parent / (build.name + "-coverage-discovery-excluded")))
    verification = json.loads(subprocess.check_output(
        args + ["--show-only=json-v1"], text=True, env=env))
    if sorted(t["name"] for t in verification["tests"]) != sorted(t["name"] for t in chosen):
        raise ValueError("CTest selection differs from the recorded invocation manifest")
    outcome = subprocess.run(args).returncode
    if outcome:
        subprocess.run(["ctest", "--test-dir", str(stage), "--rerun-failed", "-V", "-j", "1"])
    return outcome


def configuration(group, profile, build_type, multithreaded, asan, ubsan, backend=None):
    backend = backend or {"unit": "optional", "manufactured": "default"}.get(profile, "baseline")
    coverage = asan == "OFF" and ubsan == "OFF"
    options = {
        "CMAKE_BUILD_TYPE": build_type,
        "RODIN_USE_MPI": "ON", "RODIN_USE_PETSC": "ON",
        "RODIN_USE_ASAN": asan, "RODIN_USE_UBSAN": ubsan,
        "RODIN_BUILD_DOC": "OFF", "RODIN_INSTALL_RESOURCES": "OFF",
        "RODIN_USE_MCSS": "OFF", "RODIN_BUILD_EXAMPLES": "OFF",
        "RODIN_BUILD_BENCHMARKS": "OFF", "RODIN_BUILD_CONVERGENCE_TESTS": "OFF",
        "RODIN_BUILD_SRC": "ON",
        "RODIN_BUILD_UNIT_TESTS": "OFF" if profile == "manufactured" else "ON",
        "RODIN_BUILD_MANUFACTURED_TESTS": "ON" if profile != "unit" or
            group.startswith("manufactured/") or group == "src/Rodin" else "OFF",
        "RODIN_CODE_COVERAGE": "ON" if coverage else "OFF",
        "RODIN_MULTITHREADED": multithreaded,
    }
    if coverage:
        options["RODIN_LCOV_FLAGS"] = "--ignore-errors gcov,gcov --ignore-errors mismatch --ignore-errors unused"
    if profile == "convergence-baseline":
        options.update({"RODIN_USE_MPI": "OFF", "RODIN_USE_PETSC": "OFF",
                        "RODIN_BUILD_UNIT_TESTS": "OFF", "RODIN_BUILD_MANUFACTURED_TESTS": "OFF",
                        "RODIN_BUILD_CONVERGENCE_TESTS": "ON", "RODIN_CODE_COVERAGE": "ON"})
    if backend == "optional":
        options.update({"CMAKE_PREFIX_PATH": os.environ.get("SUITESPARSE_PREFIX", str(ROOT / "suitesparse")),
                        "RODIN_USE_MUMPS": "OFF" if asan == "ON" else "ON",
                        "RODIN_USE_CHOLMOD": "ON", "RODIN_USE_PARU": "ON",
                        "RODIN_USE_SPQR": "ON", "RODIN_USE_UMFPACK": "ON"})
    elif backend == "baseline":
        options.update({"RODIN_USE_SCOTCH": "ON", "RODIN_USE_UMFPACK": "ON",
                        "RODIN_USE_CHOLMOD": "ON", "RODIN_USE_SPQR": "ON"})
    return options


def configure(build, group, profile, build_type, multithreaded, asan, ubsan, backend=None):
    options = configuration(group, profile, build_type, multithreaded, asan, ubsan, backend)
    subprocess.run(["cmake", "-S", str(ROOT), "-B", str(build)] +
                   ["-D" + key + "=" + value for key, value in options.items()], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=["configure", "build", "run", "list"])
    parser.add_argument("--build", type=Path, default=Path("build"))
    parser.add_argument("--group", required=True)
    parser.add_argument("--profile", choices=["unit", "manufactured", "slow", "coverage", "convergence-baseline"], default="unit")
    parser.add_argument("--backend", choices=["optional", "default", "baseline"])
    parser.add_argument("--build-type", default="Release")
    parser.add_argument("--multithreaded", choices=["ON", "OFF"], default="ON")
    parser.add_argument("--asan", choices=["ON", "OFF"], default="OFF")
    parser.add_argument("--ubsan", choices=["ON", "OFF"], default="OFF")
    parser.add_argument("--parallel", type=int, default=2)
    parser.add_argument("--exclude-target", default="")
    parser.add_argument("--baseline", action="store_true")
    args = parser.parse_args()
    build = args.build.resolve()
    if args.operation == "configure":
        configure(build, args.group, args.profile, args.build_type,
                  args.multithreaded, args.asan, args.ubsan, args.backend)
        return 0
    if args.operation == "build":
        subprocess.run(["cmake", "--build", str(build), "--target", "RodinTestManifest"], check=True)
        configured_profiles.cache_clear()
        chosen = selected_targets(build, args.group, args.exclude_target, args.profile, args.baseline)
        (build / "ci-built-targets.json").write_text(json.dumps(chosen, indent=2))
        if chosen:
            subprocess.run(["cmake", "--build", str(build), "--parallel", str(args.parallel),
                            "--target"] + [t["target"] for t in chosen], check=True,
                           env=dict(os.environ, GCOV_PREFIX=str(build.parent / (build.name + "-coverage-discovery-excluded"))))
        else:
            print("No configured executable targets for this group", flush=True)
        return 0
    if args.operation == "list":
        _, report = prepare(build, args.group, args.profile)
        print(json.dumps(report, indent=2))
        return 0
    return run_group(build, args.group, args.profile, args.backend)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, subprocess.CalledProcessError) as error:
        sys.exit(str(error))
