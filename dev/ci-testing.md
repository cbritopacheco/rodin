# CI test registration and extension

The [test suite README](../tests/README.md) describes the categories and local
runner. The [workflow overview](../.github/workflows/README.md) describes the
shared-build graph and coverage artifacts. This guide identifies what must
change when a test or lane is extended.

Linux dependencies are provisioned through
[`apt-dependencies`](../.github/actions/apt-dependencies/action.yml). Keep package
lists in the shared environment action or dedicated lane source, rather than
adding bare `apt-get` calls to generated workflows. Download caches are keyed by
runner image, architecture, requested packages and resolved versions, so test
selection and MPI ranks do not create new package caches. Every lane still runs
and captures fresh coverage. See [dependency and disk policies](../.github/workflows/README.md#linux-dependencies-and-disk-headroom)
for mirror deadlines, exact-version fallback and conditional cleanup. Validate
changes with `test_ci_environment.py`, the full serial CI guard suite and the
generated-workflow consistency check.

When diagnosing a timeout, distinguish the CTest process watchdog from the
workflow job allowance. Compare matched configurations and complete histories;
a terminated test supplies only a runtime lower bound. Preserve every level,
control, quadrature order and rank when splitting execution. Keep timing
measurements and investigation reports outside the repository.

## Choose the scope of the change

| Change | Owning files and checks |
| --- | --- |
| Another case or source in an existing executable | Test source and its local `CMakeLists.txt`; inspect registration and execution filters. |
| Another ordinary component or module | Local and parent CMake registrations; ordinary target discovery is automatic. Review grouping in `ci_tests.py` only if the new component needs separate execution scheduling. |
| A Slow selection | CMake labels and complementary filters; confirm ordinary and Slow selections cover the intended cases without overlap. |
| A dedicated physics, geometry, degree or rank partition | The relevant matrix and execution selector in `.github/ci/targeted-lanes.yml`; `ci_builds.py` if the required targets or compile configuration change. Preserve MPI commands and ranks. |
| A new backend or compile configuration | `ci_tests.py` configuration or `ci_builds.py` definition, `ci_build_graph.py` configuration routing, and the test-environment action. Each configuration automatically gets its own hierarchy DAG. |
| A new execution/publication family | Planner inventory, workflow dependencies, report artifact family, publisher scheduling and finalizer dependencies; verify every eligible contribution is published exactly once. |
| Oversized compilation or execution | Measure the hierarchy group; review its cut level in `ci_build_graph.py` or execution grouping in `ci_tests.py` separately. Do not remove configurations or tests. |

The Python files above live in `dev/`; workflow files live in
`.github/workflows/`. The environment action is
[test-environment/action.yml](../.github/actions/test-environment/action.yml).
The existing infrastructure regressions in `dev/test_ci_*.py` guard registration,
selection, build compatibility and coverage completeness.

## Build ownership, execution lanes and publication families

A target's physical CMake registration directory remains its source owner.
Compilation stops at explicit hierarchy levels, defined by `build_group()` in
`ci_build_graph.py`:

| Category | Build group boundary | Included descendants |
| --- | --- | --- |
| Unit | `unit/Rodin/<Module>` | The module's entire subtree, including MPI and PETSc subdirectories. |
| Manufactured | `manufactured/Rodin/<Module>` | The module's entire subtree, including Solver/NewtonSolver and MPI/H1. |
| Convergence | `convergence/{h,p,hp,isoparametric}` | Every selected physics subtree under that study type. |

Targets registered directly at a category or Rodin root remain in that parent's
own group; they do not absorb separately grouped descendants. Each configured
executable has exactly one producer. Compatible regular, Slow, geometry and
rank partitions share that producer. Different compiler/backend/build-type/
threading/sanitizer configurations retain separate producers.

There is no automatic twelve-target limit, numbered build subdivision or
component-specific build exception. Execution components remain independent
lanes. Descending below an agreed build boundary requires an explicit policy
change justified by measured build time or memory, with matching dependency
and completeness guards. Existing compiler concurrency safeguards remain,
including serial compilation for the consolidated Unit/MPI group.

Coverage publication groups completed executions by family. The current graph
has seven ordinary and ten dedicated publisher families. Adding a case to an
existing lane does not add another Codecov upload. Adding a lane to an existing
family extends that family's mandatory report inventory. Adding a family must
also extend publication and finalization dependencies; a successful test job
alone does not establish a complete coverage campaign.

Configuration-local family publishers permit at most 64 submissions in the
present plan. They release independently; a slow configuration does not delay
another configuration's completed contributions.

### Regenerate the dependency graph

GitHub dependencies cannot select one child of a matrix. The generated
`TestConfiguration_*.yml` variants declare a producer for each applicable
hierarchy group inside each compatible configuration. Execution jobs wait only
for their required groups. Different groups can compile independently, and the
original Slow/geometry/rank partitions remain separate execution jobs.

After changing registration directories, execution grouping, dedicated lane
matrices or selectors, run from the repository root:

```sh
python3 dev/ci_schedule.py --write
python3 dev/ci_schedule.py
```

The first command regenerates the eleven `TestConfiguration_*.yml` variants
and ten dedicated `Test<Family>.yml` execution wrappers;
the second checks that committed output is current. Edit their source recipes,
not the generated files. The common-library recipe lives in `TestCommon.yml`;
dedicated numerical recipes and their full matrices live in
[targeted-lanes.yml](../.github/ci/targeted-lanes.yml).

Backend/instrumentation variants contain only the execution groups required by
their plan. Sanitizer variants have no coverage publisher jobs, and each
dedicated-family wrapper has just its own job. Avoid adding universal workflow
branches guarded by backend/family `if` conditions: those create skipped check
rows even when no work was ever scheduled for them. Configured ownership checks
still reject any required producer missing from a generated variant.

Source registration references establish scheduling bounds, not executable
eligibility. Configured CMake manifests still choose all build targets and must
confirm that every required producer is present in the generated dependency
edges and that every declared group has configured executables. Disabled optional
directories require matching scheduling bounds in `ci_schedule.py`; do not create
empty producers or skip conditions for them. Dynamically constructed dedicated target names need an explicit source
edge in `ci_schedule.py`; the two diffusion boundary aliases document this case.
A missing or ambiguous edge fails common planning before any test can be omitted.

## Registration procedures

### Add a test or GoogleTest suite to an existing executable

Add `TEST`, `TEST_F` or parameterized tests to an already compiled source file.
If adding a source file, also add it to the executable's CMake source list.
GoogleTest discovery registers the new cases after rebuilding; ordinary grouped
CI needs no new workflow entry. A new fixture or GoogleTest suite name does not
by itself require a new executable or GitHub job.

Check existing `TEST_FILTER` expressions and dedicated workflow selectors:
discovery cannot include a case excluded by its registration or execution
filter. Extend the appropriate partition when adding a case outside those
selectors, and verify that it is selected exactly once where intended.

### Add a component executable in an existing module

Follow the nearest module registration. For example, a new Math component can
be registered in [the Math CMake file](../tests/unit/Rodin/Math/CMakeLists.txt):

```cmake
add_executable(RodinMathNewComponentTest NewComponentTest.cpp)
target_link_libraries(RodinMathNewComponentTest
  PRIVATE GTest::gtest GTest::gtest_main Rodin::Math)
rodin_discover_tests(RodinMathNewComponentTest PROPERTIES LABELS unit)
```

The configured target manifest includes the executable automatically, and the
existing module group builds and discovers it. Do not add an executable-name
list to the workflow. If the component needs its own scheduling group because
of cost, add an explicit component entry and, when needed, an execution-group
mapping in `ci_tests.py`; update the CI guard expectations too.

MPI/PETSc executables with their own initialization must follow their backend's
existing `main` and link pattern rather than blindly adding `gtest_main`.
Use the neighboring explicit MPI `add_test` registrations, preserving launcher
flags, rank counts, timeouts, labels and external-library sanitizer handling.
See [MPI Geometry registration](../tests/unit/Rodin/MPI/Geometry/CMakeLists.txt).

### Add tests for another module

1. Mirror the production module under `tests/unit/Rodin/<Module>/`; add its
   component sources and `CMakeLists.txt` as above.
2. Include the directory with `add_subdirectory(<Module>)` in its parent test
   CMake file. Apply the same backend/feature condition as the production module.
   Merely creating a directory does not make CMake configure it.
3. Add manufactured regressions under the appropriate space, physics or backend
   directory when assembly/solver behavior needs an end-to-end check. Include
   those directories in their parents as well.
4. Regenerate the plan and inspect the configured target manifest and scoped
   selection. The new module may join a parent group rather than create a job.

Adding a new production library also requires its normal source CMake/export
registration; test discovery does not create library targets. See the
[module architecture](../doc/agents/architecture.md) and
[extension procedures](../doc/agents/workflows.md).

### Add manufactured, convergence or slow suites

Use `manufactured` for fixed-mesh numerical checks and `convergence` for
multi-resolution rate studies; preserve the full refinement history in each
study. Follow the existing
[manufactured registration](../tests/manufactured/Rodin/H1/CMakeLists.txt) and
[convergence guide](../tests/convergence/README.md).

For ordinary native convergence groups, directory discovery supplies baseline
scheduling. Dedicated native isoparametric and PETSc convergence families use
explicit matrices in `.github/ci/targeted-lanes.yml`: adding a new executable/physics family requires
updating its execution selections and `ci_builds.py` target recipe. Adding only
a rank, geometry or study partition can reuse an existing recipe if its compiled
target/configuration is unchanged. Confirm timeout budgets and all required
coverage receipts; do not assume every `convergence` label creates a dedicated
PETSc lane automatically.

For slow tests, label the registration `slow` together with its category. When
splitting one executable's cases between ordinary and slow registrations, use
complementary `TEST_FILTER` expressions so cases are neither repeated nor lost.
The Slow variant retains the directory hierarchy and reuses the regular lane's
compatible executable. It changes execution selection and watchdog policy,
not compiler flags or build ownership.

All execution families consume directory builds from the configured target
manifest. Common objects are compiled once per compatible configuration;
executables have one producer at the agreed module/study-type boundary,
including all selected descendants. Compiler, backend, optimization, threading and
sanitizer variants remain distinct. Rank and geometry selections share builds.
No execution lane should gain an inline configure/build command.

Every eligible lane saves its own production LCOV report and receipt. A new
partition must be included in its family's expected publication inventory.
Publisher jobs validate complete families and combine their contributions into
bounded incremental Codecov submissions. Finalization checks that every raw
report was published exactly once and that the published line union equals
the complete raw union before enforcing 95% coverage.

### Add an entirely new top-level test category

This is a CI architecture change, unlike adding a GoogleTest suite. Wire the new
category and enable option into CMake, define its registration labels/profile,
extend profile recording and selection in `rodin_discover_tests`/`ci_tests.py`,
and add planner/workflow scheduling. Define its backend configurations, timeout
policy and coverage eligibility. If it contributes coverage, extend the expected
receipt and compilation-inventory contracts and finalizer dependencies. Update
the registration guards and this guide. Prefer the existing categories when
they already express the test's purpose.

## Verify a registration change locally

From the repository root, inspect the generated plan outside the source tree
(Python requires the same pinned PyYAML as the workflow preflight):

```sh
python3 dev/ci_environment.py python
mkdir -p /tmp/rodin-ci-plan
python3 dev/ci_test_plan.py --output /tmp/rodin-ci-plan/ci-test-plan.json
```

This also writes `ci-build-plan.json`, `ci-convergence-plan.json` and
`ci-coverage-targets.json` beside the requested output. Compare the generated
rows before and after an edit, including compile flags, threading, sanitizers,
geometry/study selectors and MPI ranks. The configured ownership index is
generated later from actual CMake manifests; a source-only plan cannot prove
that an optional target exists or has a compatible build producer.

For an already configured build with the required backend options, build and
inspect a group, then run it. The build may compile in parallel; tests run one
invocation at a time:

```sh
python3 dev/ci_tests.py build --build build --group unit/Rodin/Geometry --profile unit
python3 dev/ci_tests.py list --build build --group unit/Rodin/Geometry --profile unit
python3 dev/ci_tests.py run --build build --group unit/Rodin/Geometry --profile unit
```

Choose the actual group from the generated plan, rather than assuming every
directory is a standalone group. Inspect `build/RodinTestTargets-*.tsv` to confirm
that the expected target is enabled. An optional-backend group can legitimately
have an empty selection in a configuration where its targets are disabled.

For changes to planning, shared compilation or coverage collection, run the CI
infrastructure guards serially:

```sh
python3 -m unittest discover -s dev -p 'test_ci_*.py' -v
```

Some artifact/coverage integration checks require Linux, CMake, a C++ compiler,
LCOV and matching gcov. A skipped tool-dependent check does not validate that
integration. Intentional changes to compilation contracts require reviewing and
updating the pinned guard expectations, rather than weakening them to accept
arbitrary output. In particular, `test_ci_build_graph.py` contains snapshots of
execution definitions and contracts: inspect the concrete inventory difference
before updating a digest. Preserve every previously required configuration,
test and rank unless an explicit change to that contract is authorized.
Matrix families must remain within GitHub's 256-entry limit;
the planner rejects an oversized family so it can be split explicitly.

For workflow edits, run actionlint as well. For build-graph edits, use the
Linux artifact/coverage integration guards: they check that common objects are
not rebuilt, restored selections contain their executables, and independent
executions start with clean counters. For publication changes, retain the
raw-to-published union equality check and the complete 95% final gate.

See also the [testing rules](../doc/agents/testing.md),
[benchmark guide](../tests/benchmarks/README.md) and
[installation checks](../tests/installation/README.md).

### Artifact retrieval and coverage collection

Build and execution lanes use exact names from their shared-build index.
Publishers and finalization use `dev/ci_artifacts.py` to derive their complete
input list from the saved plans; do not introduce workflow-wide artifact globs
or `gh run download` fallback scans. Large runs can exceed 1,000 artifacts, and
repeated REST enumeration can exhaust the installation token quota.

For a new producer, retain the naming contract in `ci_coverage.py`,
`TestGroup.yml` and `ci_artifacts.py`. Run the serial `test_ci_artifacts.py`
regressions: they compare the full configuration plan with publisher inputs and
exercise more than 1,000 artifacts, missing/empty downloads and clean transfer
retries. Their Node tests run through the Python CI guard suite, one at a time.
The final receipt/hash/inventory checks and 95% gate remain authoritative.
See the [workflow architecture](../.github/workflows/README.md) for the complete
hierarchy and coverage flow.
