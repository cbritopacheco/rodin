#!/usr/bin/env python3
"""Changed-lines clang-format check for Rodin.

Runs git-clang-format against a base revision so that ONLY the lines a
change touches must satisfy .clang-format — the existing tree is never
reformatted wholesale (see the "Known deviations" note in .clang-format).

On failure the offending hunks are printed as a colored unified diff and,
under GitHub Actions, each file/line range is emitted as an inline
annotation.
Fix locally with:

  git clang-format <base>          # rewrites the changed lines in place

Usage:
  python3 dev/check_format.py [--base REF] [--head REF]
                              [--binary clang-format] [--no-color]

Defaults: base = merge-base with origin/master (or HEAD~1), head = the
working tree (so uncommitted changes are checked too).
"""

import argparse
import os
import re
import shutil
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GITHUB = os.environ.get("GITHUB_ACTIONS") == "true"

# Only these trees are subject to formatting.
CODE_PATHS = ["src", "examples", "tests", "dev"]


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


def annotate(diff):
    """Emit GitHub annotations from the unified diff."""
    current = None
    for ln in diff.splitlines():
        m = re.match(r"\+\+\+ b/(.*)$", ln)
        if m:
            current = m.group(1)
            continue
        m = re.match(r"@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@", ln)
        if m and current:
            start = int(m.group(1))
            count = int(m.group(2) or "1")
            end = start + max(count - 1, 0)
            print(f"::error file={current},line={start},endLine={end},"
                  f"title=clang-format::lines {start}-{end} of {current} are "
                  "not clang-formatted (run: git clang-format)")


def print_diff(diff, tty):
    for ln in diff.splitlines():
        if ln.startswith(("+++", "---")):
            print(color("1", ln, tty))
        elif ln.startswith("@@"):
            print(color("36", ln, tty))
        elif ln.startswith("+"):
            print(color("32", ln, tty))
        elif ln.startswith("-"):
            print(color("31", ln, tty))
        else:
            print(ln)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
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

    # Execute git-clang-format directly: depending on the platform it may be
    # a Python script or a shell wrapper (e.g. MacPorts), so its own shebang
    # must pick the interpreter.
    cmd = [gcf, "--binary", args.binary, "--diff", base]
    if args.head:
        cmd.append(args.head)
    cmd += ["--", *CODE_PATHS]

    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    out = r.stdout.strip()

    if r.returncode != 0 and not out:
        # git-clang-format itself failed (bad ref, missing binary, ...).
        print(r.stderr.strip() or "error: git-clang-format failed")
        return 2

    clean_markers = ("no modified files to format",
                     "clang-format did not modify any files")
    if any(m in out for m in clean_markers) or not out:
        print(f"{color('1;32', 'OK', tty)}: changed lines against "
              f"{base[:12]} are clang-formatted.")
        return 0

    if not out.startswith("diff") and "---" not in out:
        # git-clang-format failed for another reason (bad ref, etc.).
        print(out)
        print(r.stderr)
        return 2

    print(f"{color('1;31', 'FAIL', tty)}: changed lines are not "
          "clang-formatted. Required changes:\n")
    print_diff(out, tty)
    if GITHUB:
        annotate(out)
    print(f"\n{color('1;33', 'fix:', tty)} run "
          f"`git clang-format {base[:12]}` and commit the result.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
