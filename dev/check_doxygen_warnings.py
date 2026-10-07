#!/usr/bin/env python3
"""Doxygen documentation-warning ratchet for Rodin.

Runs doxygen over the tree with HTML generation disabled, audits its XML
for parameter/return coverage, and compares the combined warning set against
the committed baseline dev/doxygen_warnings.baseline:

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
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

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
        "WARN_AS_ERROR = NO\n"
        f"WARN_LOGFILE = {log}\n"
    )
    doxyfile = os.path.join(tmpdir, "Doxyfile")
    with open(doxyfile, "w", encoding="utf-8") as f:
        f.write(cfg)
    # Collect all warnings for the ratchet, but never accept a failed run.
    subprocess.run([binary, doxyfile], cwd=REPO, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # Extract every function for the contract audit, including private/static
    # helpers. EXTRACT_ALL disables undocumented-entity warnings, so collect
    # the normal warning log above before doing this separate XML pass.
    # Audit the XML commands used by the published m.css documentation too.
    # The native configuration intentionally expands these aliases to nothing.
    with open(os.path.join(REPO, "doc", "Doxygen.mcss.in"), encoding="utf-8") as f:
        mcss = "\n".join(line for line in f.read().splitlines()
                         if not line.startswith("@INCLUDE"))
    cfg += "\n" + mcss + "\n"
    cfg += (
        "\n# --- full extraction for the parameter/return XML audit ---\n"
        "EXTRACT_ALL = YES\n"
        "EXTRACT_PRIVATE = YES\n"
        "EXTRACT_STATIC = YES\n"
        "GENERATE_XML = YES\n"
        f"WARN_LOGFILE = {tmpdir}/extraction-warnings.log\n"
    )
    with open(doxyfile, "w", encoding="utf-8") as f:
        f.write(cfg)
    subprocess.run([binary, doxyfile], cwd=REPO, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return log


def audit_xml(xml_directory):
    """Check parameter/return coverage in Doxygen's extracted C++ functions.

    Doxygen's warnings omit unnamed parameters and some internal helpers.
    XML retains those declarations, so check their contracts directly too.
    """
    def text(element):
        return "".join(element.itertext()).strip() if element is not None else ""

    findings = {}
    function_count = 0
    files = sorted(Path(xml_directory).glob("*.xml"))
    if not (Path(xml_directory) / "index.xml").is_file():
        raise ValueError("Doxygen produced no XML index")
    for xml_file in files:
        if xml_file.name in ("index.xml", "Doxyfile.xml"):
            continue
        for member in ET.parse(xml_file).getroot().iter("memberdef"):
            if member.get("kind") != "function":
                continue
            location = member.find("location")
            if location is None:
                continue
            filename = location.get("file", "")
            if os.path.isabs(filename):
                filename = os.path.relpath(filename, REPO)
            elif not filename.startswith("src/") and Path(REPO, "src", filename).is_file():
                filename = "src/" + filename
            if not filename.startswith("src/") or Path(filename).suffix not in (".h", ".hpp", ".cpp"):
                continue
            name = text(member.find("name"))
            # A Boost export macro is parsed as a function by Doxygen.
            if name == "BOOST_CLASS_EXPORT":
                continue
            function_count += 1
            args = text(member.find("argsstring"))
            documented = set()
            for item in member.findall('.//parameterlist[@kind="param"]/parameteritem'):
                if text(item.find("parameterdescription")):
                    documented.update(text(parameter) for parameter in
                                      item.findall("parameternamelist/parametername"))
            parameters = []
            for parameter in member.findall("param"):
                param_type = text(parameter.find("type"))
                param_name = text(parameter.find("declname")) or text(parameter.find("defname"))
                # Doxygen 1.14 can split the reference after a nested decltype
                # template argument into a separate XML param node.
                if param_type in ("&", "&&") and parameters and not parameters[-1][1]:
                    previous_type, _ = parameters.pop()
                    parameters.append((previous_type + " " + param_type, param_name))
                else:
                    parameters.append((param_type, param_name))
            messages = []
            for index, (param_type, param_name) in enumerate(parameters, 1):
                if param_type in ("void", "..."):
                    continue
                if not param_name or param_name == "...":
                    messages.append(f"parameter #{index} of member {name} is unnamed "
                                    "and lacks parameter documentation")
                elif param_name not in documented:
                    messages.append(f"parameter '{param_name}' of member {name} "
                                    "lacks parameter documentation")
            return_type = re.sub(
                r"\b(constexpr|consteval|inline|virtual|static|friend|explicit)\b",
                "", text(member.find("type"))).strip()
            # Conversion operators have no <type>; their result type is in
            # the function name (for example, "operator bool").
            if not return_type and name.startswith("operator "):
                return_type = name[len("operator "):].strip()
            has_return = any(text(section) for section in
                             member.findall('.//simplesect[@kind="return"]'))
            has_retval = any(text(item.find("parameterdescription")) for item in
                             member.findall('.//parameterlist[@kind="retval"]/parameteritem'))
            if (return_type and return_type != "void" and "=delete" not in args
                    and not has_return and not has_retval):
                messages.append(f"return type of member {name} lacks return documentation")
            for message in messages:
                warning = f"{filename}:{location.get('line', '1')}: warning: {message}"
                findings.setdefault(strip_line_number(warning), warning)
    if not function_count:
        raise ValueError("Doxygen XML contains no extracted C++ functions")
    return findings


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
            try:
                log = run_doxygen(args.doxygen, tmpdir)
            except subprocess.CalledProcessError as error:
                print(color("1;31", "error:", tty),
                      f"doxygen failed with exit status {error.returncode}")
                return 2
            if not os.path.exists(log):
                print(color("1;31", "error:", tty),
                      "doxygen produced no warning log")
                return 2
            current = normalize(log, tmpdir)
            try:
                current.update(normalize(os.path.join(tmpdir, "extraction-warnings.log"),
                                         tmpdir))
                xml_findings = audit_xml(os.path.join(tmpdir, "xml"))
                current.update(xml_findings)
                print(f"Doxygen full-extraction XML audit: {len(xml_findings)} "
                      "missing parameter/return descriptions.")
            except (ValueError, ET.ParseError, OSError) as error:
                print(color("1;31", "error:", tty), f"Doxygen XML audit failed: {error}")
                return 2

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
