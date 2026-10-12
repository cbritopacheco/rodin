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


def audit(root, reports, include_functions=False):
    root = root.resolve()
    source = root / "src" / "Rodin"
    lines = {}
    foreign = set()
    functions = {}
    for report in reports:
        current = None
        function_lines = {}
        for record in report.read_text().splitlines():
            if record.startswith("SF:"):
                function_lines = {}
                path = Path(record[3:])
                if not path.is_absolute():
                    path = root / path
                path = path.resolve()
                try:
                    relative = path.relative_to(source)
                    if relative.parts[0] == "Test":
                        raise ValueError("Library-side test helper")
                    current = str(relative)
                except ValueError:
                    foreign.add(str(path))
                    current = None
                else:
                    lines.setdefault(current, {})
            elif current is not None and include_functions and record.startswith("FNL:"):
                index, start, *_ = record[4:].split(",")
                function_lines[index] = int(start)
            elif current is not None and include_functions and record.startswith("FN:"):
                start, symbol = record[3:].split(",", 1)
                # Older LCOV records can also include an end line before the name.
                end, separator, remainder = symbol.partition(",")
                if separator and end.isdigit():
                    symbol = remainder
                function_lines[symbol] = int(start)
            elif current is not None and include_functions and record.startswith(("FNA:", "FNDA:")):
                if record.startswith("FNA:"):
                    index, hits, symbol = record[4:].split(",", 2)
                    start = function_lines.get(index)
                else:
                    hits, symbol = record[5:].split(",", 1)
                    start = function_lines.get(symbol)
                key = (current, symbol)
                observed = functions.setdefault(key, {"line": start, "hit": False})
                if observed["line"] is None:
                    observed["line"] = start
                observed["hit"] |= int(hits) > 0
            elif current is not None and record.startswith("DA:"):
                number, hits, *_ = record[3:].split(",")
                values = lines[current]
                values[int(number)] = values.get(int(number), False) or int(hits) > 0
            elif record == "end_of_record":
                current = None
    files = []
    absent = []
    for path in sorted(source.rglob("*")):
        if path.relative_to(source).parts[0] == "Test":
            continue
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
    result = {"lines": total, "hits": hits,
            "coverage_percent": 100 * hits / total if total else None,
            "nonproduction_records": sorted(foreign),
            "files": sorted(files, key=lambda file: len(file["uncovered_lines"]), reverse=True),
            "absent_sources": sorted(absent, key=lambda file: file["code_lines"], reverse=True),
            "limitation": "Line coverage measures compiled instantiations. Absent sources and "
                          "uninstantiated specializations require a separate support-matrix audit."}
    if include_functions:
        result["recorded_functions"] = {
            "observed": len(functions),
            "hit": sum(value["hit"] for value in functions.values())}
        result["zero_hit_functions"] = [
            {"file": "src/Rodin/" + file, "symbol": symbol, "line": value["line"]}
            for (file, symbol), value in sorted(functions.items()) if not value["hit"]]
        result["function_limitation"] = (
            "Zero-count recorded functions are audit candidates, including compiled template "
            "variants hidden by shared line coverage. Compiler optimization and inlining "
            "affect function counters; never-instantiated variants remain absent.")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", type=Path, nargs="+")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--functions", action="store_true",
                        help="Inventory zero-count recorded functions across compiler instantiations")
    args = parser.parse_args()
    result = audit(args.root.resolve(), args.reports, include_functions=args.functions)
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
