# Developer Tooling

The `dev/` directory contains repository-maintenance tools: ratchets,
linters, aggregate translation units, baselines, and helper scripts used by
developers and CI. These files do not define Rodin's numerical API; they keep
the tree reviewable while large historical cleanups are handled incrementally.

Use this directory when you need to answer one of these questions:

- Did my change introduce a new documentation warning?
- Did my change add a new naming or house-style violation?
- Are the lines I touched formatted?
- What baseline should shrink after I fixed old warnings?
- What aggregate translation unit should clang-tidy analyze?

## Ratchets and Baselines

Several checks use a ratchet pattern. A ratchet records the current known
backlog in a committed baseline file, compares the current run against that
baseline, and fails only on new findings. When old findings are fixed, rerun
the check with `--update-baseline` so the baseline shrinks.

Never grow a baseline to silence a new finding. Fix the finding instead, or
discuss why the rule should change.

### Doxygen Warnings

`check_doxygen_warnings.py` runs Doxygen in warning-only mode and compares the
normalized warning set against `doxygen_warnings.baseline`.

```sh
python3 dev/check_doxygen_warnings.py
python3 dev/check_doxygen_warnings.py --update-baseline
python3 dev/check_doxygen_warnings.py --doxygen /path/to/doxygen
python3 dev/check_doxygen_warnings.py --log /path/to/warnings.log
```

The baseline records the Doxygen version that produced it. Warning sets are
not compared across different versions because Doxygen changes wording and
diagnostics between releases.

### clang-tidy Naming

`check_clang_tidy.py` runs the naming rules from `.clang-tidy` over
`TidyTU.cpp` and compares findings against `clang_tidy.baseline`.

```sh
cmake -S . -B build -DRODIN_BUILD_TIDY=ON
python3 dev/check_clang_tidy.py --build-dir build
python3 dev/check_clang_tidy.py --flags "<compiler flags>"
python3 dev/check_clang_tidy.py --flags "<compiler flags>" --update-baseline
```

`TidyTU.cpp` is an aggregate translation unit that includes Rodin headers so
clang-tidy can inspect the library consistently. The CMake-backed invocation
also checks optional modules enabled in that build with their actual dependency
include paths and feature definitions.

### House Style

`style_lint.py` checks conventions that clang-format and clang-tidy do not
express:

- Boost Software License header blocks
- path-derived include guards
- no `#pragma once`
- Doxygen `@brief` coverage in headers
- PETSc includes only under `src/Rodin/PETSc/`

```sh
python3 dev/style_lint.py
python3 dev/style_lint.py src/Rodin/Variational
python3 dev/style_lint.py --update-baseline
```

Existing findings live in `style_lint.baseline`. As with every ratchet,
update it only to shrink the backlog.

### Changed-Line Formatting

`check_format.py` runs `git-clang-format` against the merge base and only
checks lines touched by the current change.

```sh
python3 dev/check_format.py
python3 dev/check_format.py --base origin/develop
python3 dev/check_format.py --binary clang-format
```

To fix reported hunks locally:

```sh
git clang-format <base>
```

## Dependency Graphs

`include_dependency_graph.py` builds an include graph for the codebase. It is
useful when checking whether a header-level change introduces an unwanted
dependency or when planning a refactor across modules.

Install the Python packages in `Requirements.txt` before running graph tools:

```sh
python3 -m pip install -r dev/Requirements.txt
python3 dev/include_dependency_graph.py --help
```

## CI Status

Static-analysis reports are published for the main branches:

| Branch  | Static Analysis |
|:-------:|:---------------:|
| master  | [![Static Analysis](https://github.com/cbritopacheco/rodin/actions/workflows/StaticAnalysis.yml/badge.svg?branch=master)](https://cbritopacheco.github.io/rodin/static/dev/refs/heads/master/) |
| develop | [![Static Analysis](https://github.com/cbritopacheco/rodin/actions/workflows/StaticAnalysis.yml/badge.svg?branch=develop)](https://cbritopacheco.github.io/rodin/static/dev/refs/heads/develop/) |

## What Belongs Here

Keep repository-wide developer utilities in `dev/` when they are meant to be
run directly by maintainers or CI. Prefer small scripts with explicit command
line interfaces, deterministic output, and committed baselines when the check
is intentionally incremental.

Do not place generated reports, experiment outputs, meshes, solver logs, or
one-off scratch files here. Those belong in a scratch directory and should not
be committed.
