#!/usr/bin/env python3
"""Draft a Doxygen @brief for each GoogleTest / Benchmark macro in a file.

For every TEST/TEST_F/TEST_P/TYPED_TEST/BENCHMARK*/... macro that lacks a
preceding doc comment, insert a `/// @brief <text>` line, where <text> is:
  1. the test's own leading `// ...` body comment, if present (promoted), or
  2. a sentence derived from the test-case name (2nd macro argument).

This is a DRAFTING aid: the derived briefs are meant to be reviewed and
refined, not committed blind. Idempotent — never touches a macro that
already has a /// or /** above it.

Usage: python3 dev/document_tests.py <file...>   # edits in place
       python3 dev/document_tests.py --dry <file> # print, don't write
"""
import re
import sys

MACRO = re.compile(r'^(\s*)((?:TEST|TEST_F|TEST_P|TYPED_TEST|TYPED_TEST_P|'
                   r'BENCHMARK_F|BENCHMARK_DEFINE_F|BENCHMARK_TEMPLATE_F)'
                   r'|BENCHMARK)\s*\(([^,]+),\s*([^)]+)\)')


def split_words(name):
    name = name.strip()
    # split snake_case and CamelCase / digit boundaries
    parts = re.split(r'[_]+', name)
    words = []
    for p in parts:
        words += re.findall(r'[A-Z]+(?=[A-Z][a-z])|[A-Z]?[a-z0-9]+|[A-Z]+|\d+', p)
    return [w for w in words if w]


def derive_from_name(testname):
    words = split_words(testname)
    if not words:
        return None
    sentence = " ".join(words)
    return sentence[0].upper() + sentence[1:] + "."


def leading_comment(lines, idx):
    """If the first non-blank line after the macro's opening brace is a //
    comment, return its text."""
    for j in range(idx + 1, min(idx + 4, len(lines))):
        s = lines[j].strip()
        if s in ("{", ""):
            continue
        m = re.match(r'//\s*(.*)', s)
        if m and m.group(1):
            t = m.group(1).strip()
            return t[0].upper() + t[1:] + ("" if t.endswith(".") else ".")
        return None
    return None


def already_documented(lines, idx):
    for j in range(idx - 1, max(idx - 3, -1), -1):
        s = lines[j].strip()
        if not s:
            continue
        return s.startswith("///") or s.endswith("*/") or s.startswith("//!")
    return False


def process(path, dry=False):
    lines = open(path, encoding="utf-8").read().splitlines()
    out, inserted = [], 0
    for i, ln in enumerate(lines):
        m = MACRO.match(ln)
        if m and not already_documented(lines, i):
            indent, _, _suite, name = m.groups()
            brief = leading_comment(lines, i) or derive_from_name(name)
            if brief:
                out.append(f"{indent}/// @brief {brief}")
                inserted += 1
        out.append(ln)
    if dry:
        for ln in out:
            if "@brief" in ln:
                print(ln)
    elif inserted:
        open(path, "w", encoding="utf-8").write("\n".join(out) + "\n")
    return inserted


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--dry"]
    dry = "--dry" in sys.argv
    total = 0
    for p in args:
        n = process(p, dry)
        total += n
        if not dry:
            print(f"{n:4d}  {p}")
    if not dry:
        print(f"total drafted: {total}")
