#!/usr/bin/env python3
"""Rodin house-style linter.

Checks the rules that clang-format and clang-tidy cannot express:

  guard      include guard must be RODIN_<PATH>_H derived from the file path
  license    every source file starts with the Boost Software License block
  filedoc    every header under src/Rodin carries Doxygen @brief documentation
  docstyle   multiline Doxygen uses /** ... */; one-line Doxygen uses ///
  forbraces  for-loop bodies spanning multiple physical lines require braces
  pragma     include guards, never #pragma once
  petsc      PETSc headers are only included under src/Rodin/PETSc/

Reporting: colored, clang-style diagnostics with the offending line and a
caret; under GitHub Actions (GITHUB_ACTIONS=true) each finding is also
emitted as a workflow annotation so it appears inline in the pull request.

Existing violations are recorded in dev/style_lint.baseline and reported
dimmed without failing the build; only NEW violations fail. Regenerate the
baseline with --update-baseline after intentional fixes (shrinking it) —
never grow it to silence a new finding.

Usage:
  python3 dev/style_lint.py [--update-baseline] [--no-color] [paths...]
"""

import argparse
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(REPO, "src", "Rodin")
BASELINE_PATH = os.path.join(REPO, "dev", "style_lint.baseline")

GITHUB = os.environ.get("GITHUB_ACTIONS") == "true"


class Palette:
    def __init__(self, enabled):
        def c(code):
            return (lambda s: f"\033[{code}m{s}\033[0m") if enabled else (lambda s: s)
        self.bold = c("1")
        self.red = c("1;31")
        self.green = c("1;32")
        self.yellow = c("1;33")
        self.blue = c("1;34")
        self.magenta = c("1;35")
        self.cyan = c("36")
        self.dim = c("2")


class Finding:
    def __init__(self, check, path, line, message, source=None, caret_span=None,
                 suggestion=None):
        self.check = check
        self.path = path            # repo-relative
        self.line = line            # 1-based, or None
        self.message = message
        self.source = source        # offending source line, if any
        self.caret_span = caret_span  # (start_col, length), 0-based
        self.suggestion = suggestion

    def key(self):
        return f"{self.check}:{self.path}"


def rel(path):
    return os.path.relpath(path, REPO)


def read_lines(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        return f.read().splitlines()


def expected_guard(relpath):
    # src/Rodin/Variational/H1/H1.h   -> RODIN_VARIATIONAL_H1_H1_H
    # src/Rodin/Variational/H1/H1.hpp -> RODIN_VARIATIONAL_H1_H1_HPP
    # (.h and .hpp must NOT share a guard: the .hpp companions hold the
    # out-of-line inline definitions and would be guarded out.)
    parts = relpath.split(os.sep)[2:]  # strip src/Rodin
    stem = parts[-1]
    suffix = "_HPP" if stem.endswith(".hpp") else "_H"
    stem = re.sub(r"\.h(pp)?$", "", stem)
    parts = parts[:-1] + [stem]
    return "RODIN_" + "_".join(p.upper() for p in parts) + suffix


CPP_COMMENTS = re.compile(
    r'R"(?P<delimiter>[^ ()\\\t\r\n]{0,16})\(.*?\)(?P=delimiter)"'
    r'|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\''
    r'|(?P<comment>//[^\n]*|/\*.*?\*/)', re.DOTALL)

CPP_TOKENS = re.compile(
    r"\b[0-9][0-9A-Za-z_'.]*|" + CPP_COMMENTS.pattern
    + r'|(?P<directive>^[ \t]*#(?:[^\n]*\\\n)*[^\n]*)'
    + r'|(?P<token>[A-Za-z_]\w*|[^\s])', re.DOTALL | re.MULTILINE)


def multiline_for_bodies(text):
    """Yield unbraced for bodies whose code spans multiple physical lines.

    Match complete statements, including nested controls, rather than counting
    semicolons inside loop headers, calls, initializers, or lambda bodies.
    Comments and preprocessor directives do not supply code tokens; literals
    remain single tokens with their original source extent.
    """
    tokens = [(m.group("token") or m.group(), m.start(), m.end())
              for m in CPP_TOKENS.finditer(text)
              if m.group("comment") is None and m.group("directive") is None]
    pairs = {}
    stack = []
    for i, (value, _, _) in enumerate(tokens):
        if value in ("(", "[", "{"):
            stack.append(i)
        elif value in (")", "]", "}") and stack:
            opening = stack.pop()
            if {"(": ")", "[": "]", "{": "}"}[tokens[opening][0]] == value:
                pairs[opening] = i

    def value(i):
        return tokens[i][0] if i < len(tokens) else None

    def control_body(i):
        i += 1
        if value(i) in ("constexpr", "co_await"):
            i += 1
        if value(i) != "(" or i not in pairs:
            return None
        return pairs[i] + 1

    def skip_attributes(i):
        while (i is not None and value(i) == "[" and value(i + 1) == "["
               and i in pairs):
            i = pairs[i] + 1
        return i

    def statement_end(i):
        i = skip_attributes(i)
        if i is None or i >= len(tokens):
            return None
        first = value(i)
        if first == "{":
            return pairs.get(i)
        if first in ("for", "if", "while", "switch"):
            end = statement_end(control_body(i))
            if first == "if" and end is not None and value(end + 1) == "else":
                return statement_end(end + 2)
            return end
        if first == "do":
            end = statement_end(i + 1)
            if end is None or value(end + 1) != "while":
                return None
            tail = control_body(end + 1)
            return tail if tail is not None and value(tail) == ";" else None
        if first == "try":
            end = statement_end(i + 1)
            while end is not None and value(end + 1) == "catch":
                end = statement_end(control_body(end + 1))
            return end
        if value(i + 1) == ":" and value(i + 2) != ":":
            return statement_end(i + 2)
        while i < len(tokens):
            if value(i) == ";":
                return i
            if value(i) == "}":
                return None
            if value(i) in ("(", "[", "{"):
                if i not in pairs:
                    return None
                i = pairs[i]
            i += 1
        return None

    for i, (keyword, start, _) in enumerate(tokens):
        if keyword != "for":
            continue
        body = control_body(i)
        if body is None or value(skip_attributes(body)) == "{":
            continue
        end = statement_end(body)
        if end is not None:
            body_start, body_end = tokens[body][1], tokens[end][2]
            if "\n" in text[body_start:body_end]:
                yield start, body_start, body_end


def doxygen_style_edits(text):
    """Yield standalone comment replacements, ignoring strings and trailing docs."""
    comments = []
    for match in CPP_COMMENTS.finditer(text):
        comment = match.group("comment")
        if comment is None:
            continue
        start, end = match.span()
        line_start = text.rfind("\n", 0, start) + 1
        line_end = text.find("\n", end)
        if line_end < 0:
            line_end = len(text)
        indent = text[line_start:start]
        if indent.strip() or text[end:line_end].strip():
            continue
        if comment.startswith("///") and not comment.startswith(("///<", "////")):
            content = comment[3:]
            if content.startswith(" "):
                content = content[1:]
            comments.append((line_start, end, indent, content))
        elif comment.startswith("/**") and not comment.startswith(("/**<", "/***")):
            content = [re.sub(r"^\s*\* ?", "", line).strip()
                       for line in comment[3:-2].splitlines()]
            content = [line for line in content if line]
            if len(content) == 1:
                yield line_start, end, indent + "/// " + content[0]
    group = []
    for entry in comments:
        if group and (text[group[-1][1]:entry[0]] != "\n"
                      or group[-1][2] != entry[2]):
            if len(group) > 1:
                yield _block_edit(group)
            group = []
        group.append(entry)
    if len(group) > 1:
        yield _block_edit(group)


def _block_edit(group):
    indent = group[0][2]
    body = "\n".join(indent + " *" + (" " + entry[3] if entry[3] else "")
                     for entry in group)
    return group[0][0], group[-1][1], indent + "/**\n" + body + "\n" + indent + " */"


def for_brace_findings(relpath, lines):
    """Check loop bodies without applying library-only header conventions."""
    text = "\n".join(lines)
    for start, _, _ in multiline_for_bodies(text):
        line = text.count("\n", 0, start) + 1
        yield Finding(
            "forbraces", relpath, line,
            "for-loop bodies spanning multiple physical lines require braces",
            source=lines[line - 1],
            suggestion="enclose the loop body in Allman braces; only a body "
                       "on one physical line may omit them")


def check_file(relpath, lines):
    findings = list(for_brace_findings(relpath, lines))
    text = "\n".join(lines)
    is_header = relpath.endswith((".h", ".hpp"))

    for start, end, replacement in doxygen_style_edits(text):
        line = text.count("\n", 0, start) + 1
        findings.append(Finding(
            "docstyle", relpath, line,
            "use /** ... */ for multiline Doxygen documentation and /// for "
            "single-line documentation",
            source=lines[line - 1],
            suggestion="replace this comment with " + replacement.splitlines()[0].strip()))

    # license: the Boost license block must appear near the top.
    head = "\n".join(lines[:10])
    if "Boost Software License" not in head:
        findings.append(Finding(
            "license", relpath, 1,
            "missing Boost Software License header block",
            source=lines[0] if lines else "",
            suggestion="start the file with the standard license comment "
                       "(copy it from any neighbouring file)"))

    if is_header:
        # pragma: no #pragma once.
        for i, ln in enumerate(lines):
            if re.match(r"\s*#\s*pragma\s+once", ln):
                findings.append(Finding(
                    "pragma", relpath, i + 1,
                    "#pragma once is not used in Rodin; use a RODIN_*_H "
                    "include guard",
                    source=ln, caret_span=(ln.find("#"), len(ln.strip())),
                    suggestion=f"replace with #ifndef {expected_guard(relpath)}"))

        # guard: first #ifndef must equal the path-derived guard.
        guard = expected_guard(relpath)
        m = None
        for i, ln in enumerate(lines):
            m = re.match(r"\s*#\s*ifndef\s+(\S+)", ln)
            if m:
                actual = m.group(1)
                if actual != guard:
                    col = ln.find(actual)
                    findings.append(Finding(
                        "guard", relpath, i + 1,
                        f"include guard '{actual}' does not match the "
                        f"path-derived name",
                        source=ln, caret_span=(col, len(actual)),
                        suggestion=f"rename the guard to {guard} "
                                   "(#ifndef, #define and the #endif comment)"))
                break
        if m is None:
            findings.append(Finding(
                "guard", relpath, 1,
                f"no include guard found; expected #ifndef {guard}",
                source=lines[0] if lines else ""))

        # filedoc: headers must carry Doxygen documentation.
        if "@brief" not in text:
            findings.append(Finding(
                "filedoc", relpath, 1,
                "header has no Doxygen documentation (no @brief anywhere)",
                suggestion="add a file-level /** @file ... @brief ... */ "
                           "block before the include guard"))

    # petsc: PETSc includes only under src/Rodin/PETSc/.
    inside_petsc = relpath.startswith(os.path.join("src", "Rodin", "PETSc"))
    if not inside_petsc:
        for i, ln in enumerate(lines):
            m = re.match(r'\s*#\s*include\s*[<"](petsc[^">]*)[">]', ln)
            if m:
                findings.append(Finding(
                    "petsc", relpath, i + 1,
                    f"PETSc header '{m.group(1)}' included outside "
                    "src/Rodin/PETSc/ — backend-independent code must not "
                    "touch PETSc",
                    source=ln, caret_span=(ln.find("#"), len(ln.strip())),
                    suggestion="move the PETSc-facing code under "
                               "src/Rodin/PETSc/ or route through the "
                               "existing wrappers"))
    return findings


def collect_files(paths):
    exts = (".h", ".hpp", ".cpp")
    files = []
    roots = paths if paths else [SRC]
    for root in roots:
        root = os.path.abspath(root)
        if os.path.isfile(root):
            if root.endswith(exts):
                files.append(root)
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = [d for d in dirnames if d not in ("build", "third-party")]
            for fn in filenames:
                if fn.endswith(exts):
                    files.append(os.path.join(dirpath, fn))
    return sorted(files)


def load_baseline():
    if not os.path.exists(BASELINE_PATH):
        return set()
    with open(BASELINE_PATH) as f:
        return {ln.strip() for ln in f
                if ln.strip() and not ln.startswith("#")}


def print_finding(p, f, baselined):
    loc = f"{f.path}:{f.line}:" if f.line else f"{f.path}:"
    tag = f"[{f.check}]"
    if baselined:
        print(p.dim(f"{loc} baselined: {f.message} {tag}"))
        return
    print(f"{p.bold(loc)} {p.red('error:')} {f.message} {p.magenta(tag)}")
    if f.source is not None and f.source.strip():
        print(f"  {p.cyan(f.source)}")
        if f.caret_span:
            start, length = f.caret_span
            print("  " + " " * max(start, 0) + p.green("^" + "~" * max(length - 1, 0)))
    if f.suggestion:
        print(f"  {p.yellow('note:')} {f.suggestion}")
    if GITHUB:
        line = f.line or 1
        print(f"::error file={f.path},line={line},"
              f"title=style_lint {f.check}::{f.message}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="*",
                    help="files or directories to lint (default: src/Rodin)")
    ap.add_argument("--update-baseline", action="store_true",
                    help="rewrite dev/style_lint.baseline from a full-tree run")
    ap.add_argument("--no-color", action="store_true")
    args = ap.parse_args()

    color = not args.no_color and (sys.stdout.isatty() or GITHUB)
    p = Palette(color)

    files = collect_files([] if args.update_baseline else args.paths)
    all_findings = []
    for path in files:
        all_findings.extend(check_file(rel(path), read_lines(path)))

    # Tests, examples, bindings, and development utilities follow the same loop
    # rule, without inheriting the library's include-guard and header-doc rules.
    if args.update_baseline or not args.paths:
        extra_files = collect_files([os.path.join(REPO, root)
                                     for root in ("tests", "examples", "py", "dev",
                                                  "src/Rodin.h")])
        for path in extra_files:
            all_findings.extend(for_brace_findings(rel(path), read_lines(path)))
        files.extend(extra_files)

    if args.update_baseline:
        keys = sorted({f.key() for f in all_findings})
        with open(BASELINE_PATH, "w") as f:
            f.write("# Known pre-existing style_lint violations (check:path).\n"
                    "# Shrink by fixing; never grow to silence a new finding.\n"
                    "# Regenerate with: python3 dev/style_lint.py --update-baseline\n")
            f.writelines(k + "\n" for k in keys)
        print(f"baseline written: {len(keys)} entries -> {rel(BASELINE_PATH)}")
        return 0

    baseline = load_baseline()
    new = [f for f in all_findings if f.key() not in baseline]
    old = [f for f in all_findings if f.key() in baseline]

    for f in old:
        print_finding(p, f, baselined=True)
    for f in new:
        print_finding(p, f, baselined=False)

    n_files = len(files)
    if new:
        print(f"\n{p.red('FAIL')}: {len(new)} new style violation(s) "
              f"({len(old)} baselined) in {n_files} file(s) checked.")
        return 1
    print(f"\n{p.green('OK')}: no new style violations "
          f"({len(old)} baselined) in {n_files} file(s) checked.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
