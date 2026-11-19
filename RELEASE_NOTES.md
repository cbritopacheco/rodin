# Rodin v1.0.0

Rodin is a modular C++20 finite element framework for solving partial
differential equations and developing shape and topology optimization
algorithms. Its variational form language expresses weak formulations through
typed C++ objects.

## Release scope

This release brings the development line beyond the original P0/P1 release:

- High-order H1 spaces and curved/isoparametric geometry, with h-, p- and
  hp-refinement verification suites.
- Volume and facet formulations, boundary and identification constraints,
  mixed problems, and linear and nonlinear solvers.
- Local Eigen, OpenMP, MPI and PETSc execution paths, with optional direct
  solver integrations.
- Mesh refinement, point location, MMG-based remeshing and adaptation,
  level-set evolution, and solid mechanics and heart models.
- Mesh and field I/O, including HDF5/XDMF output.
- Unit, manufactured-solution and convergence tests, published benchmarks,
  production coverage gates, and expanded tutorials and API documentation.

Backend support depends on the selected operation and finite element space;
consult [the backend support guide](doc/agents/backend-support.md) for the
supported combinations. The 1.0.0 version number does not make every optional
integration available in an installed package.

## Installation

A recursive clone is required because the build uses bundled submodules.
The following library and headers configuration uses Eigen solvers and does
not require PETSc, OpenMP or the resource tree:

```bash
git clone --branch v1.0.0 --recursive https://github.com/cbritopacheco/rodin.git
cd rodin
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DRODIN_USE_PETSC=OFF -DRODIN_MULTITHREADED=OFF \
  -DRODIN_BUILD_EXAMPLES=OFF -DRODIN_INSTALL_RESOURCES=OFF \
  -DCMAKE_INSTALL_PREFIX="$HOME/.local"
cmake --build build --parallel 4
cmake --install build
```

Required dependencies are CMake 3.16 or later, a C++20 compiler, Boost 1.74
or later, Eigen3 and HDF5's C library. The default configuration additionally
requires PETSc and OpenMP. To use or install the example and test resources,
install Git LFS and run `git lfs pull` before enabling resource installation.
Detailed instructions are in [INSTALL.md](INSTALL.md).

Downstream CMake projects should enable C and CXX, use
`find_package(Rodin 1.0.0 REQUIRED)`, and link the required `Rodin::` targets.
Rebuild applications against this release; source compatibility with 0.0.1
has not been established across the intervening API changes. Current tutorials
and examples provide the reference usage.

CMake version matching now requires the same major version. A request for
0.0.1 will not silently accept 1.0.0; update the requested version after
checking the application against this release.

## Installation limits

- `Rodin::MMG` is available when building in the source tree but is not
  exported by the installed CMake package.
- Python bindings are experimental and require an explicit source build;
  this release preparation does not include Python wheels.
- Optional solvers require their external dependencies and matching CMake
  feature options. A backend configuration is not a guarantee that every
  finite element/operator combination is implemented.

See [CHANGELOG.md](CHANGELOG.md) for release history and
[RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md) for the outstanding publication
gates. Rodin is distributed under the [Boost Software License 1.0](LICENSE).
