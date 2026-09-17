#!/usr/bin/env python3
"""Changed-lines clang-format check for Rodin.

Only the lines a change touches must satisfy .clang-format — the existing
tree is never reformatted wholesale (see the "Known deviations" note in
.clang-format). That policy is defined by ``git clang-format`` against a
base revision.

Why the extra full-file cross-check:

  ``git clang-format`` reformats only the changed line *ranges*. With a
  partial range clang-format lacks the surrounding scope context, so with
  options like ``IndentAccessModifiers: true`` it mis-computes indentation
  and proposes changes to lines that are in fact already correct. The result
  is base-dependent false positives — a file "fails" against one base and
  passes against another with byte-identical content.

  So every hunk ``git clang-format`` proposes is validated against the
  canonical, context-complete result of formatting the *whole* file: a hunk
  is a real violation only if full-file clang-format would touch the same
  original lines. Hunks the full-file format leaves alone are partial-range
  artifacts and are dropped. This keeps the changed-lines-only policy while
  making the check base-independent.

Fix locally with:

  git clang-format <base>          # rewrites the changed lines in place
  clang-format -i <file>           # or format the whole file, base-independent

Usage:
  python3 dev/check_format.py [--base REF] [--head REF]
                              [--binary clang-format] [--no-color]

Defaults: base = merge-base with origin/master (or HEAD~1), head = the
working tree (so uncommitted changes are checked too).
"""

import argparse
import difflib
import os
import re
import shutil
import subprocess
import sys

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


def find_git_clang_format(binary, override=None):
    """Locate git-clang-format matching the clang-format binary."""
    if override:
        path = shutil.which(override) or override
        if os.path.exists(path):
            return path
        print(f"error: git-clang-format '{override}' not found")
        sys.exit(2)
    suffix = ""
    m = re.search(r"clang-format(.*)$", os.path.basename(binary))
    if m:
        suffix = m.group(1)
    for name in (f"git-clang-format{suffix}", "git-clang-format"):
        path = shutil.which(name)
        if path:
            return path
    print("error: git-clang-format not found in PATH")
    sys.exit(2)


def parse_hunks(diff):
    """Split a unified diff into {path: [(orig_start, orig_end, text)]}.

    orig_start/orig_end are 1-based inclusive line numbers on the '-' side
    (the current file); a pure insertion (orig_count == 0) anchors on the
    line it follows.
    """
    files = {}
    path = None
    cur = None  # (orig_start, orig_end, [lines])
    for ln in diff.splitlines():
        m = re.match(r"^\+\+\+ b/(.*)$", ln)
        if m:
            path = m.group(1)
            continue
        if ln.startswith("--- "):
            continue
        m = HUNK_RE.match(ln)
        if m and path is not None:
            ostart = int(m.group(1))
            ocount = int(m.group(2)) if m.group(2) is not None else 1
            if ocount == 0:
                lo, hi = ostart, ostart  # insertion after line `ostart`
            else:
                lo, hi = ostart, ostart + ocount - 1
            cur = (lo, hi, [ln])
            files.setdefault(path, []).append(cur)
        elif cur is not None and (ln.startswith(("+", "-", " ")) or ln == ""):
            cur[2].append(ln)
        else:
            cur = None
    return files


def full_file_bad_lines(binary, path, content):
    """Original 1-based line numbers that full-file clang-format would change."""
    r = subprocess.run([binary, f"--assume-filename={path}"],
                       input=content, capture_output=True, text=True)
    if r.returncode != 0:
        return None
    if r.stdout == content:
        return set()
    o = content.splitlines(keepends=True)
    f = r.stdout.splitlines(keepends=True)
    bad = set()
    for tag, i1, i2, _j1, _j2 in difflib.SequenceMatcher(
            None, o, f, autojunk=False).get_opcodes():
        if tag == "equal":
            continue
        if tag == "insert":
            bad.update({i1, i1 + 1})  # 1-based neighbours of the anchor
        else:
            bad.update(range(i1 + 1, i2 + 1))
    return bad


def annotate(path, lo, hi):
    print(f"::error file={path},line={lo},endLine={max(hi, lo)},"
          f"title=clang-format::lines {lo}-{max(hi, lo)} of {path} are not "
          "clang-formatted (run: git clang-format)")


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
    ap.add_argument("--git-clang-format", dest="gcf", default=None,
                    help="explicit git-clang-format script (e.g. when the "
                         "clang-format binary comes from a pip wheel)")
    ap.add_argument("--no-color", action="store_true")
    args = ap.parse_args()

    tty = not args.no_color and (sys.stdout.isatty() or GITHUB)
    base = resolve_base(args.base)
    gcf = find_git_clang_format(args.binary, args.gcf)

    # Candidate set: the changed-line hunks git-clang-format wants to rewrite.
    cmd = [gcf, "--binary", args.binary, "--diff", base]
    if args.head:
        cmd.append(args.head)
    cmd += ["--", *CODE_PATHS]
    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    out = r.stdout.strip()

    if r.returncode != 0 and not out:
        print(r.stderr.strip() or "error: git-clang-format failed")
        return 2

    clean_markers = ("no modified files to format",
                     "clang-format did not modify any files")
    if any(m in out for m in clean_markers) or not out:
        print(f"{color('1;32', 'OK', tty)}: changed lines against "
              f"{base[:12]} are clang-formatted.")
        return 0

    if not out.startswith("diff") and "---" not in out:
        print(out)
        print(r.stderr)
        return 2

    # Validate each proposed hunk against the full-file (context-complete)
    # result, dropping partial-range false positives.
    hunks = parse_hunks(out)
    bad_cache = {}
    real = []  # (path, lo, hi, text_lines)
    for path, file_hunks in hunks.items():
        if path not in bad_cache:
            if args.head:
                cf = git("show", f"{args.head}:{path}", check=False)
                content = cf.stdout if cf.returncode == 0 else None
            else:
                try:
                    with open(os.path.join(REPO, path), encoding="utf-8") as fh:
                        content = fh.read()
                except (OSError, UnicodeDecodeError):
                    content = None
            bad_cache[path] = (full_file_bad_lines(args.binary, path, content)
                               if content is not None else None)
        bad = bad_cache[path]
        for lo, hi, text in file_hunks:
            # None => could not format the file: keep the hunk to be safe.
            if bad is None or bad & set(range(lo, hi + 1)):
                real.append((path, lo, hi, text))

    if not real:
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
          f"`git clang-format {base[:12]}` (or `clang-format -i <file>`) "
          "and commit the result.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
