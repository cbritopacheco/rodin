# Developer tooling

`dev/` contains repository-maintenance tools used by developers and CI.
It is separate from the numerical library in `src/Rodin/` and the numerical
suites in `tests/`.

## Structure

| Files | Responsibility |
| --- | --- |
| `ci_test_plan.py` | Test directory scheduling, dedicated coverage partition inventory and shared-build planning. |
| `ci_tests.py` | Configured target selection, scoped GoogleTest discovery and serial CTest execution. |
| `ci_builds.py` | Compatible build definitions, immutable artifact publication and validated restoration. |
| `ci_coverage.py` | Execution receipts, coverage contribution validation, pooling and the final production coverage gate. |
| `audit_coverage.py`, `audit_source.py` | Compiled coverage inventories and source-level function/exception inventories, respectively. |
| `format_cpp.py`, `check_format.py` | Canonical C++ formatting and the changed-line formatting check. |
| `style_lint.py` | House rules beyond formatter and compiler checks. |
| `check_clang_tidy.py`, `TidyTU.cpp` | Naming checks over an aggregate translation unit with configured dependencies. |
| `check_doxygen_warnings.py` | Doxygen diagnostics and generated XML documentation checks. |
| `*.baseline` | Committed reference findings for incremental checks. |
| `include_dependency_graph.py` | Include-dependency graph generation. |
| `profile_aabb.py`, `summarize_xctrace.py` | Profiling support and trace summarization. |
| `test_*.py` | Regression tests for the developer tools. |
| `Requirements.txt` | Dependencies for the developer utilities that use them. |

## Relationships

The CI planner reads test registration directories; configured CMake manifests
supply authoritative target information. Build artifacts contain compiled code,
while execution jobs own their runtime counters. Coverage collection validates
those independent contributions before taking their production-line union.
The [workflow overview](../.github/workflows/README.md) describes the job graph.

Incremental style and documentation checks compare current findings with their
committed baselines. Formatting uses the canonical formatter rather than a
separate CI-only interpretation of indentation.

## Supporting documents

- [Tooling procedures](tooling.md): invocation and baseline maintenance.
- [CI test registration guide](ci-testing.md): adding tests, modules and suites.
- [Coverage and source review history](coverage-history.md): investigation
  records and historical measurements, separate from the directory architecture.

Generated profiles, meshes, solver logs and scratch outputs belong outside this
directory. Shared tracked resources live in `resources/` under Git LFS.
