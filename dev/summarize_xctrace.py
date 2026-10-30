#!/usr/bin/env python3
"""Summarize exported Instruments time-profile samples with explicit selectors.

Export the time-profile table with xctrace, then select both the workload's
owning frame and the operation's kernel through repeated --require expressions.
This excludes fixture construction that can otherwise dominate a capture.
Weights are sampled CPU time, not elapsed time. Inclusive frames and categories
can overlap; their percentages must not be added. Referenced values are retained
even in rows without stacks. XML rows are discarded after processing rather
than retaining the complete trace in memory.
"""

import argparse
from collections import Counter
import json
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET


class XctraceProfile:
  """Attribute selected samples without interpreting them as benchmark timings."""

  def __init__(self, source, required=(), process=None, categories=(), limit=20,
               excluded=()):
    self.source = source
    self.required = tuple(re.compile(pattern) for pattern in required)
    self.excluded = tuple(re.compile(pattern) for pattern in excluded)
    self.process = re.compile(process) if process else None
    self.categories = tuple((name, re.compile(pattern))
                            for name, pattern in categories)
    if len({name for name, _ in self.categories}) != len(self.categories):
      raise ValueError("Category names must be unique")
    self.limit = limit

  def summarize(self):
    """Return JSON-compatible attribution; reject unresolved references."""
    frames, stacks, weights, processes = {}, {}, {}, {}
    leaf, inclusive, categories = Counter(), Counter(), Counter()
    total_weight = selected_weight = samples = selected_samples = 0
    schema_seen = False
    parents = []
    for event, element in ET.iterparse(self.source, events=("start", "end")):
      if event == "start":
        parents.append(element)
        if element.tag == "schema" and element.get("name") == "time-profile":
          schema_seen = True
        continue
      if element.tag == "row":
        for process in element.iter("process"):
          if process.get("id"):
            processes[process.get("id")] = process.get("fmt", "")
        weight_node = element.find("weight")
        weight = None
        if weight_node is not None:
          if weight_node.get("ref"):
            weight = weights[weight_node.get("ref")]
          else:
            weight = int(weight_node.text)
            if weight_node.get("id"):
              weights[weight_node.get("id")] = weight
          if weight < 0:
            raise ValueError("A time-profile weight cannot be negative")
          samples += 1
          total_weight += weight
        stack_node = element.find("tagged-backtrace")
        names = ()
        if stack_node is not None:
          if stack_node.get("ref"):
            names = stacks[stack_node.get("ref")]
          else:
            names_list = []
            for frame in stack_node.findall("frame"):
              if frame.get("ref"):
                name = frames[frame.get("ref")]
              else:
                name = frame.get("name", "<unnamed>")
                if frame.get("id"):
                  frames[frame.get("id")] = name
              names_list.append(name)
            names = tuple(names_list)
            if stack_node.get("id"):
              stacks[stack_node.get("id")] = names
        process_node = element.find("process")
        process_name = ""
        if process_node is not None:
          identity = process_node.get("ref") or process_node.get("id")
          process_name = processes[identity]
        selected = names and weight is not None and (
          self.process is None or self.process.search(process_name)) and all(
          any(pattern.search(name) for name in names) for pattern in self.required) and not any(
          pattern.search(name) for pattern in self.excluded for name in names)
        if selected:
          selected_samples += 1
          selected_weight += weight
          leaf[names[0]] += weight
          for name in set(names):
            inclusive[name] += weight
          for category, pattern in self.categories:
            if any(pattern.search(name) for name in names):
              categories[category] += weight
        if len(parents) > 1:
          parents[-2].remove(element)
        element.clear()
      parents.pop()
    if not schema_seen:
      raise ValueError("Expected an exported time-profile schema")
    if selected_samples == 0:
      raise ValueError("No samples match the explicit process/stack selection")
    if selected_weight == 0:
      raise ValueError("Selected samples have no positive time weight")

    def entries(counter):
      return [{"name": name, "sampled_ms": weight / 1_000_000,
               "percent_selected": 100 * weight / selected_weight}
              for name, weight in sorted(counter.items(),
                key=lambda item: (-item[1], item[0]))[:self.limit]]

    return {
      "source": str(self.source),
      "required_stack_patterns": [pattern.pattern for pattern in self.required],
      "excluded_stack_patterns": [pattern.pattern for pattern in self.excluded],
      "process_pattern": self.process.pattern if self.process else None,
      "category_patterns": {name: pattern.pattern for name, pattern in self.categories},
      "samples_with_weight": samples,
      "selected_samples": selected_samples,
      "total_sampled_ms": total_weight / 1_000_000,
      "selected_sampled_ms": selected_weight / 1_000_000,
      "selected_percent_total": 100 * selected_weight / total_weight,
      "leaf": entries(leaf),
      "inclusive": entries(inclusive),
      "categories": entries(categories),
      "interpretation": "Sample weights, not wall time; inclusive categories overlap",
    }


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("source", type=Path)
  parser.add_argument("--require", action="append", default=[],
                      help="Regex that must match a frame; repeated selectors are ANDed")
  parser.add_argument("--process", help="Regex selecting the recorded process name")
  parser.add_argument("--exclude", action="append", default=[],
                      help="Reject a stack if any frame matches; repeat to add exclusions")
  parser.add_argument("--category", action="append", default=[], metavar="NAME=REGEX",
                      help="Inclusive attribution; categories may overlap")
  parser.add_argument("--limit", type=int, default=20)
  args = parser.parse_args()
  if args.limit <= 0:
    parser.error("--limit must be positive")
  categories = []
  for category in args.category:
    name, separator, pattern = category.partition("=")
    if not separator or not name or not pattern:
      parser.error("--category expects a nonempty NAME=REGEX pair")
    categories.append((name, pattern))
  try:
    result = XctraceProfile(args.source, args.require, args.process,
                           categories, args.limit, args.exclude).summarize()
  except (OSError, ET.ParseError, KeyError, ValueError, re.error) as error:
    parser.exit(2, f"Cannot summarize trace: {error}\n")
  json.dump(result, sys.stdout, indent=2)
  sys.stdout.write("\n")


if __name__ == "__main__":
  main()
