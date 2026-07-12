#!/usr/bin/env python3
"""Doxygen documentation-warning ratchet for Rodin.

Runs doxygen over the tree (HTML generation disabled — warnings only) and
compares the normalized warning set against the committed baseline
dev/doxygen_warnings.baseline:

  * a warning NOT in the baseline is NEW -> reported precisely and the
    check fails;
  * a baseline warning that no longer occurs is FIXED -> reported so the
    baseline can be shrunk (never grown) with --update-baseline.

This turns "the docs emit ten thousand warnings" into "your change may not
add a single one", without requiring the backlog to be fixed first.

The warning set depends on the doxygen version; the baseline records the
version it was generated with, and the check refuses to compare across
versions (CI pins the same version, see .github/workflows/Style.yml).

Reporting is colored and, under GitHub Actions, each new warning is also
emitted as an inline annotation.

Usage:
  python3 dev/check_doxygen_warnings.py [--update-baseline] [--doxygen BIN]
  python3 dev/check_doxygen_warnings.py --log FILE   # reuse an existing log
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASELINE_PATH = os.path.join(REPO, "dev", "doxygen_warnings.baseline")
GITHUB = os.environ.get("GITHUB_ACTIONS") == "true"


def color(code, s, enabled=True):
    return f"\033[{code}m{s}\033[0m" if enabled else s


def doxygen_version(binary):
    out = subprocess.run([binary, "--version"], capture_output=True, text=True,
                         check=True).stdout.strip()
    return out.split()[0]


def run_doxygen(binary, tmpdir):
    template = os.path.join(REPO, "doc", "Doxygen.in")
    with open(template, encoding="utf-8") as f:
        cfg = f.read()
    cfg = cfg.replace("@CMAKE_SOURCE_DIR@", REPO)
    cfg = cfg.replace("@CMAKE_BINARY_DIR@", tmpdir)
    log = os.path.join(tmpdir, "warnings.log")
    cfg += (
        "\n# --- overrides appended by dev/check_doxygen_warnings.py ---\n"
        f"OUTPUT_DIRECTORY = {tmpdir}\n"
        "GENERATE_HTML = NO\n"
        "GENERATE_LATEX = NO\n"
        "GENERATE_XML = NO\n"
        "QUIET = YES\n"
        f"WARN_LOGFILE = {log}\n"
    )
    doxyfile = os.path.join(tmpdir, "Doxyfile")
    with open(doxyfile, "w", encoding="utf-8") as f:
        f.write(cfg)
    subprocess.run([binary, doxyfile], cwd=REPO, check=False,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return log


def strip_line_number(warning):
    """Comparison key: drop the :line: so edits above a warning do not turn
    it into a "new" one. Reporting still uses the full current line."""
    return re.sub(r"^([^:]+):\d+:", r"\1:", warning)


def normalize(log_path, tmpdir=None):
    """Return {comparison_key: representative full warning line}."""
    warnings = {}
    with open(log_path, encoding="utf-8", errors="replace") as f:
        for ln in f:
            ln = ln.rstrip("\n")
            if ": warning:" not in ln and not ln.startswith("warning:"):
                continue  # drop continuation/detail lines
            ln = ln.replace(REPO + os.sep, "")
            if tmpdir:
                ln = ln.replace(tmpdir + os.sep, "")
            warnings.setdefault(strip_line_number(ln), ln)
    return warnings


def load_baseline():
    if not os.path.exists(BASELINE_PATH):
        return None, []
    version = None
    entries = []
    with open(BASELINE_PATH, encoding="utf-8") as f:
        for ln in f:
            ln = ln.rstrip("\n")
            m = re.match(r"#\s*doxygen\s+(\S+)", ln)
            if m:
                version = m.group(1)
            if ln and not ln.startswith("#"):
                entries.append(ln)
    return version, entries


def emit(new_warnings, tty):
    pat = re.compile(r"^(?P<file>[^:]+):(?P<line>\d+):\s*warning:\s*(?P<msg>.*)$")
    for w in new_warnings:
        m = pat.match(w)
        if m:
            loc = f"{m.group('file')}:{m.group('line')}:"
            print(f"{color('1', loc, tty)} {color('1;31', 'new warning:', tty)} "
                  f"{m.group('msg')}")
            if GITHUB:
                print(f"::error file={m.group('file')},line={m.group('line')},"
                      f"title=doxygen warning::{m.group('msg')}")
        else:
            print(f"{color('1;31', 'new warning:', tty)} {w}")
            if GITHUB:
                print(f"::error title=doxygen warning::{w}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--doxygen", default=os.environ.get("DOXYGEN", "doxygen"))
    ap.add_argument("--log", help="use an existing WARN_LOGFILE instead of "
                                  "running doxygen")
    ap.add_argument("--update-baseline", action="store_true")
    args = ap.parse_args()

    tty = sys.stdout.isatty() or GITHUB

    if args.log:
        version = doxygen_version(args.doxygen)
        current = normalize(args.log)
    else:
        version = doxygen_version(args.doxygen)
        with tempfile.TemporaryDirectory() as tmpdir:
            log = run_doxygen(args.doxygen, tmpdir)
            if not os.path.exists(log):
                print(color("1;31", "error:", tty),
                      "doxygen produced no warning log")
                return 2
            current = normalize(log, tmpdir)

    if args.update_baseline:
        entries = sorted(current.keys())
        with open(BASELINE_PATH, "w", encoding="utf-8") as f:
            f.write("# Doxygen warning baseline for the ratchet check.\n"
                    f"# doxygen {version}\n"
                    "# Entries are line-number-agnostic (file: warning text) so\n"
                    "# unrelated edits do not shift warnings into 'new' status.\n"
                    "# Shrink by fixing warnings, regenerate with\n"
                    "#   python3 dev/check_doxygen_warnings.py --update-baseline\n"
                    "# Never grow this file to silence a new warning.\n")
            f.writelines(w + "\n" for w in entries)
        print(f"baseline written: {len(entries)} warnings "
              f"(doxygen {version}) -> {os.path.relpath(BASELINE_PATH, REPO)}")
        return 0

    base_version, baseline = load_baseline()
    if base_version is None:
        print(color("1;31", "error:", tty),
              "no baseline found; run with --update-baseline first")
        return 2
    if base_version != version:
        print(color("1;33", "error:", tty),
              f"baseline was generated with doxygen {base_version} but this "
              f"is doxygen {version}; warning sets are not comparable.\n"
              "Install the pinned version (see .github/workflows/Style.yml) "
              "or regenerate the baseline.")
        return 2

    baseline_set = {strip_line_number(w) for w in baseline}
    current_keys = set(current.keys())
    new = sorted(current[k] for k in current_keys - baseline_set)
    fixed = sorted(baseline_set - current_keys)

    emit(new, tty)

    if fixed:
        print(f"{color('1;32', 'fixed:', tty)} {len(fixed)} baseline warning(s) "
              "no longer occur — shrink the baseline with --update-baseline.")
    if new:
        print(f"\n{color('1;31', 'FAIL', tty)}: {len(new)} new doxygen "
              f"warning(s) (baseline {len(baseline_set)}, current "
              f"{len(current_keys)}).")
        return 1
    print(f"\n{color('1;32', 'OK', tty)}: no new doxygen warnings "
          f"(baseline {len(baseline_set)}, current {len(current_keys)}, "
          f"{len(fixed)} fixed).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
