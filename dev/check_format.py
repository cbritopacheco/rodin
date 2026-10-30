#!/usr/bin/env python3
"""Changed-lines clang-format check for Rodin.

Only lines touched by a change must satisfy the canonical format. Complete
files are formatted through format_cpp.py to retain scope context; proposed
edits are intersected with changed lines from Git. This also enforces the
plain-struct indentation that clang-format alone cannot express. The check
does not write the index or temporary Git trees. Workspace files are unchanged
unless --fix is explicitly requested.

Fix locally with:

  python3 dev/check_format.py --base REF --fix [paths...]
  python3 dev/format_cpp.py -i <file>  # formats the selected complete file

Usage:
  python3 dev/check_format.py [--base REF] [--head REF]
                              [--binary clang-format] [--no-color] [--fix]
                              [paths...]

Defaults: base = merge-base with origin/master (or HEAD~1), head = the
working tree (so uncommitted changes are checked too).
"""

import argparse
import difflib
import os
import re
import subprocess
import sys

from format_cpp import format_source
from style_lint import CPP_TOKENS

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GITHUB = os.environ.get("GITHUB_ACTIONS") == "true"

# Only these trees are subject to formatting.
CODE_PATHS = ["src", "examples", "tests", "dev"]

HUNK_RE = re.compile(r"^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@")


def color(code, s, enabled=True):
    return f"\033[{code}m{s}\033[0m" if enabled else s


def git(*args, check=True):
    return subprocess.run(["git", *args], cwd=REPO, capture_output=True,
                          text=True, check=check)


def resolve_base(base):
    if base:
        if git("rev-parse", "--verify", "--quiet", base + "^{commit}",
               check=False).returncode == 0:
            return base
        print(f"note: base '{base}' is not a valid commit; falling back")
    for candidate in ("origin/master", "origin/develop"):
        r = git("merge-base", candidate, "HEAD", check=False)
        if r.returncode == 0 and r.stdout.strip():
            return r.stdout.strip()
    return "HEAD~1"


def formatting_hunks(content, formatted):
    """Yield zero-context edits and their original inclusive line ranges.

    Disable the sequence matcher's popularity heuristic: repeated braces and
    common declaration lines must still anchor edits in large C++ files.
    Insertions are associated with both neighbours, including the first line.
    """
    original = content.splitlines(keepends=True)
    canonical = formatted.splitlines(keepends=True)
    # Anchor equal source lines without their indentation. Otherwise a long
    # reindent of repetitive code can be represented as a distant deletion and
    # insertion, which is unsafe to apply only in part.
    matcher = difflib.SequenceMatcher(
        None, [line.lstrip(" \t") for line in original],
        [line.lstrip(" \t") for line in canonical], autojunk=False)
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            edits = [(i, i + 1, j1 + i - i1, j1 + i - i1 + 1)
                     for i in range(i1, i2) if original[i] != canonical[j1 + i - i1]]
        else:
            edits = [(i1, i2, j1, j2)]
        for a, b, c, d in edits:
            old_count, new_count = b - a, d - c
            old_start = a + 1 if old_count else a
            new_start = c + 1 if new_count else c
            lines = [f"@@ -{old_start},{old_count} +{new_start},{new_count} @@"]
            lines += ["-" + line.rstrip("\r\n") for line in original[a:b]]
            lines += ["+" + line.rstrip("\r\n") for line in canonical[c:d]]
            yield (a + 1, b, lines) if old_count else (a, a + 1, lines)


def source_tokens(content):
    """Fingerprint code and literals, conservatively retaining directives verbatim."""
    return tuple(m.group("token") or m.group() for m in CPP_TOKENS.finditer(content)
                 if m.group("comment") is None)


def format_changed_source(content, formatted, touched):
    """Apply canonical edits intersecting touched original lines, and no others.

    The edits come from the same context-complete diff used by the checker.
    Replacement text is taken from the canonical source, retaining exact line
    endings and pure insertions rather than reconstructing it from diagnostics.
    """
    original = content.splitlines(keepends=True)
    canonical = formatted.splitlines(keepends=True)
    result = []
    cursor = 0
    for lo, hi, text in formatting_hunks(content, formatted):
        if not touched.intersection(range(lo, hi + 1)):
            continue
        match = HUNK_RE.match(text[0])
        start, count, new_start, new_count = map(int, match.groups())
        i = start - 1 if count else start
        j = new_start - 1 if new_count else new_start
        result.extend(original[cursor:i])
        result.extend(canonical[j:j + new_count])
        cursor = i + count
    result.extend(original[cursor:])
    fixed = "".join(result)
    if source_tokens(fixed) != source_tokens(content):
        raise ValueError("Partial formatting would change C++ tokens, literals or "
                         "directives; this file was not written")
    return fixed


def annotate(path, lo, hi):
    print(f"::error file={path},line={lo},endLine={max(hi, lo)},"
          f"title=clang-format::lines {lo}-{max(hi, lo)} of {path} are not "
          "canonically formatted (run: python3 dev/format_cpp.py -i <file>)")


def print_diff_lines(lines, tty):
    for ln in lines:
        if ln.startswith("@@"):
            print(color("36", ln, tty))
        elif ln.startswith("+"):
            print(color("32", ln, tty))
        elif ln.startswith("-"):
            print(color("31", ln, tty))
        else:
            print(ln)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=os.environ.get("FORMAT_BASE", ""))
    ap.add_argument("--head", default="",
                    help="check base..HEAD instead of base..worktree")
    ap.add_argument("--binary", default="clang-format")
    ap.add_argument("--git-clang-format", help=argparse.SUPPRESS)
    ap.add_argument("--no-color", action="store_true")
    ap.add_argument("--fix", action="store_true",
                    help="apply canonical edits touching changed working-tree lines")
    ap.add_argument("paths", nargs="*", help="limit checks/fixes to these paths")
    args = ap.parse_args()
    if args.fix and args.head:
        ap.error("--fix cannot be combined with --head; historical files are read-only")

    tty = not args.no_color and (sys.stdout.isatty() or GITHUB)
    base = resolve_base(args.base)
    while True:
        # Read changed ranges without temporary Git trees or index writes. Format
        # with complete scope context, then intersect canonical edits with those
        # ranges. The same canonical formatter is used locally and in CI.
        revision = [base, args.head] if args.head else [base]
        changed = git("diff", "--no-ext-diff", "--no-color", "--unified=0", *revision,
                      "--", *(args.paths or CODE_PATHS))
        ranges = {}
        path = None
        for line in changed.stdout.splitlines():
            if line.startswith("+++ "):
                path = line[6:] if line.startswith("+++ b/") else None
            match = HUNK_RE.match(line)
            if match and path and path.endswith((".h", ".hpp", ".cpp", ".cc", ".cxx")):
                start = int(match.group(3))
                count = int(match.group(4) or 1)
                ranges.setdefault(path, set()).update(
                    range(start, start + count) if count else (start, start + 1))

        real = []  # (path, lo, hi, text_lines)
        fixed_files = 0
        for path, touched in ranges.items():
            if args.head:
                content = git("show", f"{args.head}:{path}").stdout
            else:
                with open(os.path.join(REPO, path), encoding="utf-8") as fh:
                    content = fh.read()
            formatted = format_source(args.binary, os.path.join(REPO, path), content)
            if args.fix:
                try:
                    fixed = format_changed_source(content, formatted, touched)
                except ValueError as error:
                    print(f"error: {path}: {error}")
                    return 2
                if fixed != content:
                    with open(os.path.join(REPO, path), "w", encoding="utf-8") as fh:
                        fh.write(fixed)
                    fixed_files += 1
                    print(f"formatted changed regions: {path}")
                continue
            for lo, hi, text in formatting_hunks(content, formatted):
                if touched & set(range(lo, hi + 1)):
                    real.append((path, lo, hi, text))

        if not args.fix or fixed_files == 0:
            break

    if not real:
        if args.fix:
            print(f"Applied and verified changed-region formatting against {base[:12]}.")
            return 0
        print(f"{color('1;32', 'OK', tty)}: changed lines against "
              f"{base[:12]} are clang-formatted.")
        return 0

    print(f"{color('1;31', 'FAIL', tty)}: changed lines are not "
          "clang-formatted. Required changes:\n")
    last = None
    for path, lo, hi, text in real:
        if path != last:
            print(color("1", f"--- a/{path}", tty))
            print(color("1", f"+++ b/{path}", tty))
            last = path
        print_diff_lines(text, tty)
        if GITHUB:
            annotate(path, lo, hi)
    print(f"\n{color('1;33', 'fix:', tty)} run "
          f"`python3 dev/check_format.py --base {base[:12]} --fix` "
          "and commit the result.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
