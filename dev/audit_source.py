#!/usr/bin/env python3
"""Inventory production definitions independently of C++ instantiation.

Requires tree-sitter and tree-sitter-cpp. Parsing deliberately does not run the
preprocessor: inactive configurations and uninstantiated templates stay visible.
Parser recovery locations are reported, never silently treated as audited code.
An absent compiler line is an investigation candidate, not proof of nonexecution
(inlining and constant evaluation can remove runtime attribution).
"""

import argparse
import hashlib
import json
from pathlib import Path
import re


def descendants(node):
    yield node
    for child in node.children:
        yield from descendants(child)


def ancestors(node):
    while node is not None and node.parent is not None:
        node = node.parent
        yield node


def audit(root, line_inventory):
    from tree_sitter import Language, Parser
    import tree_sitter_cpp

    parser = Parser(Language(tree_sitter_cpp.language()))
    files, definitions, exceptions, recoveries = [], [], [], []
    for path in sorted((root / "src").rglob("*")):
        relative = path.relative_to(root)
        if path.suffix not in {".h", ".hpp", ".cpp", ".cc", ".cxx"}:
            continue
        if "Test" in relative.parts:
            continue
        source = path.read_bytes()
        # This declaration-generating macro lacks a trailing semicolon in Rodin.
        # Blank only its spelling, preserving every source byte offset and line.
        parsed = re.sub(rb"BOOST_SERIALIZATION_SPLIT_MEMBER\(\)",
                        lambda match: b" " * len(match[0]), source)
        tree = parser.parse(parsed)
        name = relative.as_posix()
        nodes = list(descendants(tree.root_node))

        def text(node):
            return source[node.start_byte:node.end_byte].decode()

        observed = {int(line): hit for line, hit in line_inventory.get(name, {}).items()}
        function_count = 0
        for node in nodes:
            if node.type == "function_definition":
                body = node.child_by_field_name("body")
                if body is not None:
                    function_count += 1
                    start, end = body.start_point.row + 1, body.end_point.row + 1
                    lines = {line: hit for line, hit in observed.items()
                             if start <= line <= end}
                    definition_lines = {line: hit for line, hit in observed.items()
                                        if node.start_point.row + 1 <= line <= end}
                    definitions.append({
                        "file": name, "line": node.start_point.row + 1,
                        "end_line": node.end_point.row + 1, "body_line": start,
                        "signature": source[node.start_byte:body.start_byte].decode().strip(),
                        "template": any(p.type == "template_declaration" for p in ancestors(node)),
                        "empty_body": not any(c.is_named and c.type != "comment" for c in body.children),
                        "has_initializers": any(c.type == "field_initializer_list" for c in node.children),
                        "measured_body_lines": len(lines), "hit_body_lines": sum(lines.values()),
                        "measured_definition_lines": len(definition_lines),
                        "hit_definition_lines": sum(definition_lines.values()),
                        "uncovered_lines": sorted(line for line, hit in lines.items() if not hit),
                        "uncovered_definition_lines": sorted(line for line, hit in definition_lines.items() if not hit),
                        "parse_error": node.has_error,
                    })
            if node.type == "ERROR" or node.is_missing:
                recoveries.append({"file": name, "line": node.start_point.row + 1,
                                   "end_line": node.end_point.row + 1, "text": text(node)[:180]})
            if node.type == "throw_statement":
                exceptions.append({"file": name, "line": node.start_point.row + 1,
                                   "kind": "throw", "text": text(node)})
            if node.type != "qualified_identifier" or not re.fullmatch(
                    r"(?:Rodin::)?Alert::(?:Exception|MemberFunctionException|NamespacedException|ClassException)",
                    text(node)):
                continue
            function = next((p for p in ancestors(node) if p.type == "function_definition"), None)
            context, signature = "type_reference", None
            if function is not None:
                body = function.child_by_field_name("body")
                signature = source[function.start_byte:body.start_byte].decode().strip()
                declarator = function.child_by_field_name("declarator")
                owner = any(p.type in {"class_specifier", "struct_specifier"}
                            for p in ancestors(function))
                qualified = declarator is not None and "::" in text(declarator).split("(", 1)[0]
                context = ("static" if re.search(r"\bstatic\b", signature)
                           else "instance_member" if owner or qualified else "free")
            exceptions.append({"file": name, "line": node.start_point.row + 1,
                               "kind": text(node), "context": context, "signature": signature})
        files.append({"file": name, "sha256": hashlib.sha256(source).hexdigest(),
                      "definitions": function_count, "parse_error": tree.root_node.has_error})
    return {"files": files, "definitions": definitions, "exceptions": exceptions,
            "parse_errors": recoveries,
            "limitations": ["Syntactic inventory, independent of preprocessing and instantiation.",
                            "Parser recovery locations require review; generated macro bodies are not expanded.",
                            "Line attribution does not prove execution of each template specialization.",
                            "No source coverage percentage is inferred from this inventory."]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--line-inventory", type=Path,
                        help="Optional JSON mapping repository-relative files to line/hit dictionaries")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    lines = json.loads(args.line_inventory.read_text()) if args.line_inventory else {}
    result = audit(args.root.resolve(), lines)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: len(result[key]) for key in
                      ("files", "definitions", "exceptions", "parse_errors")}))


if __name__ == "__main__":
    main()
