#!/usr/bin/env python3
"""Build a diagnostic AABB executable without changing production code.

Requires a Ninja Release build with compile_commands.json and the selected
benchmark target already built. The ordinary target measures timings; the
generated executable counts operations only. Its timings are not performance
evidence.
"""

import argparse
import json
from pathlib import Path
import re
import shlex
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", choices=["RodinAABBWorkload", "RodinMPIAABBBenchmarks"],
                        default="RodinAABBWorkload")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    build = args.build.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = (root / "src/Rodin/Location/AABB.h").read_text()

    def replace_once(before, after):
        nonlocal source
        if source.count(before) != 1:
            raise RuntimeError(f"Instrumentation anchor changed: {before!r}")
        source = source.replace(before, after, 1)

    replace_once("namespace Rodin::Location\n{", """namespace Rodin::Location
{
  struct AABBWorkCounters
  {
    size_t candidates = 0, transforms = 0, jacobians = 0, iterations = 0;
    size_t retries = 0, boxCandidates = 0, projectionRejected = 0, indexBytes = 0;
  };
  inline thread_local AABBWorkCounters aabbWork;
""")
    replace_once("        const Real physTol = physicalTolerance();\n\n",
                 "        ++aabbWork.candidates;\n        const Real physTol = physicalTolerance();\n\n")
    replace_once("          if (seed == 0)\n",
                 "          aabbWork.retries += seed != 0;\n          if (seed == 0)\n")
    replace_once("            if (residualNorm == Real(0))\n",
                 "            ++aabbWork.iterations;\n            if (residualNorm == Real(0))\n")
    for operation, counter, expected in [("transform", "transforms", 6),
                                          ("jacobian", "jacobians", 1)]:
        pattern = rf"^( +)(transformation\.{operation}\(.+;)$"
        source, count = re.subn(pattern,
            rf"\1++aabbWork.{counter};\n\1\2", source, flags=re.MULTILINE)
        if count != expected:
            raise RuntimeError(f"Expected {expected} {operation} anchors, got {count}")
    replace_once("              bool outsideHull = false;",
                 "              ++aabbWork.boxCandidates;\n              bool outsideHull = false;")
    replace_once("              if (outsideHull)\n                continue;",
                 "              if (outsideHull)\n              {\n                ++aabbWork.projectionRejected;\n                continue;\n              }")
    replace_once("        return index;", """        aabbWork.indexBytes = index.nodes.capacity() * sizeof(Node) +
          index.entries.capacity() * sizeof(Index) +
          (index.entryLo.capacity() + index.entryHi.capacity()) * sizeof(Bound) +
          index.projections.capacity() * sizeof(ProjectionBound) +
          index.projectionRanges.capacity() * sizeof(ProjectionRange);
        return index;""")
    header = output / "AABBDiagnostic.h"
    header.write_text(source)
    entries = json.loads((build / "compile_commands.json").read_text())
    source_name = "MPIAABB.cpp" if args.target == "RodinMPIAABBBenchmarks" else "AABBWorkload.cpp"
    entry = next(e for e in entries if e["file"].endswith("/" + source_name))
    command = shlex.split(entry["command"])
    obj = output / (args.target + "Diagnostic.o")
    command[command.index("-o") + 1] = str(obj)
    command[1:1] = ["-DRODIN_AABB_WORKLOAD_DIAGNOSTICS", "-include", str(header),
                    "-I" + str(root / "src/Rodin/Location")]
    subprocess.run(command, cwd=entry["directory"], check=True)
    lines = subprocess.check_output(["ninja", "-C", str(build), "-t", "commands",
                                     args.target], text=True).splitlines()
    link = shlex.split(lines[-1].split(" && ")[1])
    executable = output / (args.target + "Diagnostic")
    link[link.index("-o") + 1] = str(executable)
    link = [str(obj) if p.endswith("/" + source_name + ".o") else p for p in link]
    subprocess.run(link, cwd=build, check=True)
    print(executable)


if __name__ == "__main__":
    main()
