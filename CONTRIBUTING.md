# Contributing to Rodin

This document describes the code style, how it is enforced, and how to run
the checks locally. The style is not aspirational — it is derived from the
existing code, and the tooling is configured to match it.

## The style in five layers

| Layer | What | Enforced by |
|-------|------|-------------|
| 1. Formatting | Indentation, braces, wrapping | [`.clang-format`](.clang-format) — CI checks **changed lines only** |
| 2. Naming | Identifier conventions | [`.clang-tidy`](.clang-tidy) + [`dev/check_clang_tidy.py`](dev/check_clang_tidy.py) — ratchet against a findings baseline |
| 3. House rules | Include guards, license blocks, header docs, PETSc containment | [`dev/style_lint.py`](dev/style_lint.py) |
| 4. Documentation | No new Doxygen warnings | [`dev/check_doxygen_warnings.py`](dev/check_doxygen_warnings.py) |
| 5. Design | The semantic style | Review, guided by [`doc/agents/philosophy.md`](doc/agents/philosophy.md) |

All four mechanical layers run in CI ([Style workflow](.github/workflows/Style.yml))
and report failures as inline pull-request annotations with the exact file,
line, offending source, and suggested fix.

## Quick reference: the style itself

- **Formatting**: 2-space indent, Allman braces, `public:` indented inside
  the class with members one level deeper, constructor initializers one per
  line starting with `: ` on their own line, empty bodies as `{}` on their
  own line, pointers/references bind left (`Base* copy()`), includes are
  **never reordered** (order is load-bearing).
- **Naming**: types/namespaces/template parameters `CamelCase`;
  functions/methods `camelBack` (accessors `getX`/`setX`, chainable setters
  return `*this`); private/protected members `m_camelBack`; parameters,
  locals, public members `camelBack`; macros `RODIN_UPPER_CASE`.
  **Mathematical notation is welcome** where an identifier denotes a
  mathematical object: `F`, `Jinv`, `I1`, `m_I4` are all conformant.
- **Files**: one header = one concept; include guard `RODIN_<PATH>_H`
  derived from the file path (never `#pragma once`); Boost license block at
  the top; a `/** @file ... @brief ... */` documentation block; a class's
  `FormLanguage::Traits` specialization, its definition, its
  specializations, and its deduction guides live together in that header.
- **Boundaries**: PETSc calls only under `src/Rodin/PETSc/`
  (`assert(ierr == PETSC_SUCCESS)` after each call — no checking macros);
  third-party integrations stay in their own directories; the core never
  includes them.
- **Design** (layer 5, reviewed not linted): the program text reads as the
  mathematics; expression nodes are cloned, not shared; behavior varies by
  template parameter, not by subclass or flag; constraints inside
  optimization solves are smooth penalties, never hard projections; changes
  are minimal and behavior-preserving — no speculative fast paths. The full
  statement lives in [`doc/agents/philosophy.md`](doc/agents/philosophy.md)
  and [`doc/agents/conventions.md`](doc/agents/conventions.md).

## Running the checks locally

```sh
# 1. Formatting — check, then fix, the lines you changed:
python3 dev/check_format.py --base origin/master
git clang-format origin/master        # applies the fixes in place

# 2. Naming — full-tree ratchet over the aggregate TU:
python3 dev/check_clang_tidy.py --flags "-std=c++20 -Isrc -Ibuild/src \
    -Ithird-party/termcolor/include -isystem /usr/include/eigen3"

# 3. House rules:
python3 dev/style_lint.py                       # full tree
python3 dev/style_lint.py src/Rodin/Geometry    # one subtree

# 4. Documentation ratchet (requires doxygen 1.14.0, the pinned version):
python3 dev/check_doxygen_warnings.py
```

Requirements: clang-format ≥ 16 and clang-tidy ≥ 16 (CI uses 18),
Python ≥ 3.8, doxygen 1.14.0 for the ratchet.

## The ratchet principle

The existing tree is not fully conformant, and wholesale cleanup commits
are deliberately avoided (they poison `git blame` and conflict with every
open branch). Instead every check is a **ratchet**:

- clang-format applies to **changed lines only**;
- clang-tidy compares against `dev/clang_tidy.baseline` — new naming
  violations fail, old ones are tracked until fixed;
- `dev/style_lint.baseline` records pre-existing house-rule violations —
  new ones fail, old ones are reported dimmed. Fixing files shrinks the
  baseline (`--update-baseline`); growing it to silence a finding is not
  acceptable;
- `dev/doxygen_warnings.baseline` records pre-existing documentation
  warnings under the same contract.

Touch a file, leave it cleaner than you found it — the baselines only move
in one direction.

## Tests

A change is not done until the unit test executable for the touched module
passes (`build/tests/unit/Rodin/<Module>/Rodin<Module><Component>Test` —
ctest is not wired up), and, for assembly/solver changes, the relevant
`tests/manufactured` suites pass too. Numerical claims are made with
measurements (see `doc/agents/testing.md` for the house test patterns).

## Housekeeping

- Run examples from a scratch directory — they write `*.h5`/`*.xdmf`/logs
  into the CWD; never commit these artifacts.
- Branches: `master` is the default; development happens on `module/*`,
  `model/*` topic branches off `develop`.
