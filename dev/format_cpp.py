#!/usr/bin/env python3
"""Canonical C++ formatting: clang-format plus plain-struct indentation.

clang-format 18's IndentAccessModifiers indents access labels by two spaces
and members by four, but also indents plain struct members by four. The house
style requires two spaces in structs without explicit access sections. This
pass removes that single extra level, including in nested declarations; it
never changes tokens or string-literal contents. Use this entry point rather
than raw clang-format when formatting Rodin C++ files.
"""

import argparse
from pathlib import Path
import re
import subprocess
import sys

from style_lint import CPP_TOKENS


def plain_struct_indentation(text):
    """Remove one clang-format indentation level inside each plain struct."""
    matches = list(CPP_TOKENS.finditer(text))
    tokens = [(m.group("token") or m.group(), m.start(), m.end())
              for m in matches
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

    lines = text.splitlines(keepends=True)
    # A multiline literal's interior whitespace is data, not indentation.
    protected = set()
    for m in matches:
        if m.group("comment") is not None and m.group().startswith("/*"):
            start = text.count("\n", 0, m.start())
            end = text.count("\n", 0, m.end())
            protected.update(range(start + 1, end + 1))
        if (m.group("token") is None and m.group("comment") is None
                and m.group("directive") is None):
            start = text.count("\n", 0, m.start())
            end = text.count("\n", 0, m.end())
            protected.update(range(start + 1, end + 1))

    shifts = [0] * len(lines)
    for i, (value, _, _) in enumerate(tokens):
        if value != "struct" or (i and tokens[i - 1][0] == "enum"):
            continue
        j = i + 1
        angle = 0
        named = False
        inherited = False
        while j < len(tokens):
            value = tokens[j][0]
            if value == "[" and j in pairs:  # C++ attributes.
                j = pairs[j] + 1
                continue
            if value == "(" and angle and j in pairs:
                j = pairs[j] + 1
                continue
            if angle == 0 and value == ":":
                inherited = True
            if angle == 0 and not inherited and re.fullmatch(r"[A-Za-z_]\w*", value):
                if named and value != "final":
                    break  # `struct Type variable{...}` is not a definition.
                named = True
            if value == "<":
                angle += 1
            elif value == ">":
                if angle == 0:
                    break  # Elaborated type in a template argument.
                angle -= 1
            elif value in (";", "=", "}") or (
                    value in ("(", ")", "]", ",", "*", "&") and angle == 0):
                break  # Forward declaration, variable, or function return type.
            elif value == "{" and angle == 0:
                end = pairs.get(j)
                if end is None:
                    break
                k = j + 1
                explicit = False
                while k < end:
                    if (tokens[k][0] in ("public", "protected", "private")
                            and tokens[k + 1][0] == ":"):
                        explicit = True
                        break
                    if tokens[k][0] in ("(", "[", "{") and k in pairs:
                        k = pairs[k] + 1
                    else:
                        k += 1
                if not explicit:
                    first = text.count("\n", 0, tokens[j][1]) + 1
                    last = text.count("\n", 0, tokens[end][1])
                    for line in range(first, last):
                        shifts[line] += 2
                break
            j += 1

    result = []
    for i, line in enumerate(lines):
        if i not in protected and shifts[i]:
            leading = len(line) - len(line.lstrip(" "))
            line = line[min(leading, shifts[i]):]
        result.append(line)
    # ReflowComments:false deliberately leaves comment indentation unchanged
    # in clang-format. Align standalone comments with the following code after
    # removing record levels. Normalize star prefixes independently: clang-format
    # may move the opening line without moving the comment body. Preserve text
    # after each star and verbatim lines without a star prefix.
    for m in matches:
        if m.group("comment") is None:
            continue
        first = text.count("\n", 0, m.start())
        last = text.count("\n", 0, m.end())
        prefix = text[text.rfind("\n", 0, m.start()) + 1:m.start()]
        if prefix.strip():
            continue  # Trailing comments stay attached to their code.
        following = next((token for token in tokens if token[1] >= m.end()), None)
        if following is None:
            continue
        next_line = text.count("\n", 0, following[1])
        target = len(result[next_line]) - len(result[next_line].lstrip(" "))
        if following[0] == "}" and not shifts[first]:
            continue  # No declaration follows this comment outside a plain struct.
        if following[0] == "}":
            target += 2
        result[first] = " " * target + result[first].lstrip(" \t")
        for line in range(first + 1, last + 1):
            stripped = result[line].lstrip(" \t")
            if stripped.startswith("*"):
                result[line] = " " * (target + 1) + stripped
    return "".join(result)


def format_source(binary, path, content):
    """Format complete source with scope context and the repository style."""
    result = subprocess.run(
        [binary, "--style=file", f"--assume-filename={path}"],
        input=content, capture_output=True, text=True, check=True)
    return plain_struct_indentation(result.stdout)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="clang-format")
    parser.add_argument("-i", action="store_true", help="format files in place")
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()
    if not args.i and len(args.files) != 1:
        parser.error("stdout formatting accepts exactly one file")
    for filename in args.files:
        path = Path(filename)
        formatted = format_source(args.binary, str(path), path.read_text())
        if args.i:
            path.write_text(formatted)
        else:
            sys.stdout.write(formatted)


if __name__ == "__main__":
    main()
