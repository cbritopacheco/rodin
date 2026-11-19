# Rodin 1.0.0 release checklist

This checklist tracks preparation of the 1.0.0 release. An unchecked item is
pending verification. The historical 0.0.1 installation report is not evidence
for this release.

## Release contents

- [x] Set the CMake, Python package and citation versions to 1.0.0.
- [x] Prepare [release notes](RELEASE_NOTES.md) and a changelog entry.
- [x] Document current dependencies, defaults and installation limitations.
- [x] Retain the Boost Software License in [LICENSE](LICENSE).
- [ ] Confirm the release date in `CITATION.cff` and replace `Unreleased` in
  `CHANGELOG.md` with the publication date.

## Local verification

Run builds outside the source tree. The minimal configuration below disables
PETSc, OpenMP, examples and resource installation. Resource-dependent tests
still require hydrated Git LFS files.

```bash
git submodule update --init --recursive
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release \
  -DRODIN_USE_PETSC=OFF -DRODIN_MULTITHREADED=OFF \
  -DRODIN_BUILD_EXAMPLES=OFF -DRODIN_INSTALL_RESOURCES=OFF \
  -DRODIN_BUILD_UNIT_TESTS=ON -DRODIN_BUILD_MANUFACTURED_TESTS=ON
cmake --build build/release --parallel 4
ctest --test-dir build/release/tests -L unit -LE 'slow|distributed' --output-on-failure
ctest --test-dir build/release/tests -L manufactured -LE 'slow|distributed' --output-on-failure
```

The existing [scoped runner](dev/ci-testing.md) can build and execute selected
directories without discovering unbuilt siblings. Record exactly which tests
were run; a selected suite does not certify the complete matrix.

- [x] Build and install the library, then compile and run the installed-package
  Poisson test, including its package/header version checks.
- [x] Run representative unit and manufactured-solution tests.
- [x] Verify rejection of the previous octal version encoding.
- [x] Review the final diff and run `git diff --check`.

Local checks on 2026-11-19 used macOS arm64, Clang 19, Boost 1.88, Eigen3
and HDF5, with PETSc and Rodin MPI support disabled and OpenMP enabled.
MPI compiler wrappers supplied the parallel HDF5 installation's include/link
requirements. The following selected executables passed:

- `RodinUtilityOptionalReferenceTest`: 4 tests.
- `RodinMathVectorTest`: 37 tests.
- `RodinManufacturedPoisson`: 72 tests.

The installed-package Poisson consumer compiled and ran. Reintroducing the old
padded version code in a scratch installation made its static assertion fail;
restoring the decimal code restored a passing build and run. CMake package
checks accepted requests for 1.0 and 1.0.0 and rejected 0.0.1, 1.0.1 and 2.0.0.
Changed-line formatting, the installation-test style check and whitespace
checks passed. These selected local results leave the complete CI gates below
pending.

The MMG wrapper also built after setting the local toolchain option
`CMAKE_CXX_FLAGS=-include cmath`. Without that option, MMG's C-header include
order triggered missing `NAN`/`INFINITY` diagnostics in this host's Clang 19
standard library. This workaround changed no tracked source or submodule.

The automated installation matrix uses:

```bash
RODIN_INSTALL_RESOURCES=OFF bash tests/installation/test_installation.sh
```

It covers native, MPI and PETSc installs on Linux and macOS through
[Tests.yml](.github/workflows/Tests.yml). Local toolchain paths and dependency
prefixes may need to be supplied for the host system.

## Required current-commit CI evidence

- [ ] All build and test configurations pass, including MPI/PETSc,
  manufactured solutions, convergence, slow tests and sanitizers.
- [ ] Installation matrix passes, including downstream compilation and execution.
- [ ] Style, formatting, clang-tidy and Doxygen checks pass.
- [ ] Documentation builds and its generated pages render correctly.
- [ ] Benchmarks complete; inspect the published results for regressions.
- [ ] `Coverage / Finalize` succeeds for the release commit and the complete
  current-commit report satisfies the 95% production-line floor.
- [ ] Verify stable finalization checks are required by branch protection.

See [the workflow guide](.github/workflows/README.md) for required finalization
receipts. Partial uploads and previous-commit results cannot satisfy these gates.

## Publication

Publish only after the PR is reviewed and all release gates pass.

- [ ] Merge the release PR into `master` and verify the final commit's checks.
- [ ] Tag that reviewed commit `v1.0.0` and push the tag.
- [ ] Prepare a source distribution containing initialized submodules. If it
  includes resources, hydrate Git LFS objects and verify there are no pointers.
  GitHub's automatic source archives do not include submodule contents.
- [ ] Test a fresh build and downstream install from the extracted distribution.
- [ ] Publish the GitHub release with [RELEASE_NOTES.md](RELEASE_NOTES.md),
  the tested distribution and its checksum.
- [ ] Verify the release download, documentation links and citation date.

Python bindings remain experimental; this checklist does not authorize a
PyPI upload or certify Python wheels. Native binary packages are not part of
this source-release preparation.
