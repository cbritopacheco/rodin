#!/usr/bin/env python3
"""clang-tidy naming ratchet for Rodin.

Runs clang-tidy (readability-identifier-naming, per .clang-tidy) over the
aggregate translation unit dev/TidyTU.cpp with the full src/Rodin header
filter, and compares the findings against the committed baseline
dev/clang_tidy.baseline:

  * a finding NOT in the baseline is NEW -> reported with file, line,
    identifier, and suggested fix; the check fails;
  * a baseline finding that no longer occurs is FIXED -> reported so the
    baseline can be shrunk (never grown) with --update-baseline.

Findings are keyed by (file, identifier kind, identifier) — line-number
agnostic, so unrelated edits do not shift old findings into "new" status.

Under GitHub Actions each new finding is emitted as an inline annotation.

Usage:
  python3 dev/check_clang_tidy.py --build-dir build \
      [--clang-tidy BIN] [--update-baseline]
  python3 dev/check_clang_tidy.py --flags "<compiler flags>" \
      [--clang-tidy BIN] [--update-baseline]
"""

import argparse
import os
import re
import shlex
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASELINE_PATH = os.path.join(REPO, "dev", "clang_tidy.baseline")
GITHUB = os.environ.get("GITHUB_ACTIONS") == "true"

FINDING = re.compile(
    r"^(?P<file>[^\s:][^:]*):(?P<line>\d+):(?P<col>\d+):\s*"
    r"(?:error|warning):\s*invalid case style for "
    r"(?P<kind>[a-z ]+) '(?P<name>[^']+)'")


def color(code, s, enabled=True):
    return f"\033[{code}m{s}\033[0m" if enabled else s


def run(binary, flags, build_dir):
    cmd = [binary, os.path.join(REPO, "dev", "TidyTU.cpp"),
           "--header-filter=(^|/)src/Rodin/.*", "--quiet"]
    if build_dir:
        cmd.extend(["-p", os.path.abspath(build_dir)])
    else:
        cmd.extend(["--"] + shlex.split(flags))
    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    findings = {}
    fatal = []
    for ln in r.stdout.splitlines():
        m = FINDING.match(ln)
        if m:
            rel = os.path.relpath(m.group("file"), REPO) \
                if os.path.isabs(m.group("file")) else m.group("file")
            key = f"{rel}\t{m.group('kind').strip()}\t{m.group('name')}"
            findings.setdefault(key, (rel, m.group("line"), ln))
        elif "clang-diagnostic-error" in ln:
            fatal.append(ln)
    return findings, fatal


def load_baseline():
    if not os.path.exists(BASELINE_PATH):
        return None
    with open(BASELINE_PATH, encoding="utf-8") as f:
        return {ln.rstrip("\n") for ln in f
                if ln.strip() and not ln.startswith("#")}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clang-tidy", default=os.environ.get("CLANG_TIDY", "clang-tidy"))
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--build-dir",
                      help="CMake build directory containing compile_commands.json")
    mode.add_argument("--flags", help="compiler flags for the aggregate TU")
    ap.add_argument("--update-baseline", action="store_true")
    args = ap.parse_args()

    tty = sys.stdout.isatty() or GITHUB
    findings, fatal = run(args.clang_tidy, args.flags, args.build_dir)

    if fatal:
        print(color("1;31", "error:", tty),
              "the aggregate TU does not compile; fix these first:")
        for ln in fatal[:10]:
            print(" ", ln)
        return 2

    if args.update_baseline:
        with open(BASELINE_PATH, "w", encoding="utf-8") as f:
            f.write("# clang-tidy naming baseline (file<TAB>kind<TAB>identifier).\n"
                    "# Shrink by fixing names; regenerate with\n"
                    "#   python3 dev/check_clang_tidy.py --flags ... --update-baseline\n"
                    "# Never grow this file to silence a new finding.\n")
            f.writelines(k + "\n" for k in sorted(findings))
        print(f"baseline written: {len(findings)} findings -> "
              f"{os.path.relpath(BASELINE_PATH, REPO)}")
        return 0

    baseline = load_baseline()
    if baseline is None:
        print(color("1;31", "error:", tty),
              "no baseline found; run with --update-baseline first")
        return 2

    new = sorted(k for k in findings if k not in baseline)
    fixed = sorted(baseline - set(findings))

    for k in new:
        rel, line, raw = findings[k]
        print(f"{color('1', f'{rel}:{line}:', tty)} "
              f"{color('1;31', 'new naming violation:', tty)} "
              f"{raw.split('error:', 1)[-1].split('warning:', 1)[-1].strip()}")
        if GITHUB:
            _, kind, name = k.split("\t")
            print(f"::error file={rel},line={line},"
                  f"title=clang-tidy naming::invalid case style for {kind} "
                  f"'{name}'")

    if fixed:
        print(f"{color('1;32', 'fixed:', tty)} {len(fixed)} baseline "
              "finding(s) no longer occur — shrink the baseline with "
              "--update-baseline.")
    if new:
        print(f"\n{color('1;31', 'FAIL', tty)}: {len(new)} new naming "
              f"violation(s) (baseline {len(baseline)}, current "
              f"{len(findings)}).")
        return 1
    print(f"\n{color('1;32', 'OK', tty)}: no new naming violations "
          f"(baseline {len(baseline)}, current {len(findings)}, "
          f"{len(fixed)} fixed).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
