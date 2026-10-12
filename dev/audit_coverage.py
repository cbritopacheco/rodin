#!/usr/bin/env python3
"""Audit production LCOV records, including headers absent from the report.

Absence is an audit candidate, not proof that a header contains uncovered
executable code: declaration-only headers and static data are also absent.
LCOV cannot certify coverage of template specializations never instantiated.
"""

import argparse
import json
from pathlib import Path
import re


def audit(root, reports):
    root = root.resolve()
    source = root / "src" / "Rodin"
    lines = {}
    foreign = set()
    for report in reports:
        current = None
        for record in report.read_text().splitlines():
            if record.startswith("SF:"):
                path = Path(record[3:])
                if not path.is_absolute():
                    path = root / path
                path = path.resolve()
                try:
                    current = str(path.relative_to(source))
                except ValueError:
                    foreign.add(str(path))
                    current = None
                else:
                    lines.setdefault(current, {})
            elif current is not None and record.startswith("DA:"):
                number, hits, *_ = record[3:].split(",")
                values = lines[current]
                values[int(number)] = values.get(int(number), False) or int(hits) > 0
            elif record == "end_of_record":
                current = None
    files = []
    absent = []
    for path in sorted(source.rglob("*")):
        if path.suffix not in {".h", ".hpp", ".cpp"}:
            continue
        name = str(path.relative_to(source))
        text = path.read_text()
        # Preserve newlines when removing comments, so candidate locations
        # remain directly usable in the editor.
        code = re.sub(r"/\*.*?\*/|//[^\n]*",
                      lambda match: "\n" * match[0].count("\n"), text, flags=re.S)
        templates = [code.count("\n", 0, match.start()) + 1
                     for match in re.finditer(r"\btemplate\s*<", code)]
        if name not in lines or not lines[name]:
            absent.append({"file": "src/Rodin/" + name,
                           "template_declarations": templates,
                           "code_lines": sum(bool(line.strip()) for line in code.splitlines())})
            continue
        observed = lines[name]
        files.append({"file": "src/Rodin/" + name,
                      "lines": len(observed), "hits": sum(observed.values()),
                      "uncovered_lines": sorted(n for n, hit in observed.items() if not hit),
                      "template_declarations": templates})
    total = sum(file["lines"] for file in files)
    hits = sum(file["hits"] for file in files)
    return {"lines": total, "hits": hits,
            "coverage_percent": 100 * hits / total if total else None,
            "nonproduction_records": sorted(foreign),
            "files": sorted(files, key=lambda file: len(file["uncovered_lines"]), reverse=True),
            "absent_sources": sorted(absent, key=lambda file: file["code_lines"], reverse=True),
            "limitation": "Line coverage measures compiled instantiations. Absent sources and "
                          "uninstantiated specializations require a separate support-matrix audit."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", type=Path, nargs="+")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = audit(args.root.resolve(), args.reports)
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    else:
        print(encoded, end="")
    if result["nonproduction_records"]:
        parser.exit(1, "Coverage contains records outside src/Rodin.\n")
    if not result["lines"]:
        parser.exit(1, "Coverage contains no production executable lines.\n")


if __name__ == "__main__":
    main()
