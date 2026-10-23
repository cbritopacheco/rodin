# CI test registration and extension


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
explicit matrices in `Tests.yml`: adding a new executable/physics family requires
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
executables have one physical-directory producer, subdivided into numbered
groups above twelve targets. Compiler, backend, optimization, threading and
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
python3 -m pip install PyYAML==6.0.2
mkdir -p /tmp/rodin-ci-plan
python3 dev/ci_test_plan.py --output /tmp/rodin-ci-plan/ci-test-plan.json
```

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
arbitrary output. Matrix families must remain within GitHub's 256-entry limit;
the planner rejects an oversized family so it can be split explicitly.

See also the [testing rules](../doc/agents/testing.md),
[benchmark guide](../tests/benchmarks/README.md) and
[installation checks](../tests/installation/README.md).
