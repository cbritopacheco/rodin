# Rodin

[![License](https://img.shields.io/badge/license-BSL--1.0-green)](https://github.com/cbritopacheco/rodin/blob/master/LICENSE)

| Branch      |  Matrix  | Documentation |
|:-----------:|:--------:|:-------------:|
| master      | [![Build](https://github.com/cbritopacheco/rodin/actions/workflows/Build.yml/badge.svg?branch=master)](https://github.com/cbritopacheco/rodin/actions/workflows/Build.yml) | [![Documentation](https://github.com/cbritopacheco/rodin/actions/workflows/Documentation.yml/badge.svg)](https://cbritopacheco.github.io/rodin/) |
| develop     | [![Build](https://github.com/cbritopacheco/rodin/actions/workflows/Build.yml/badge.svg?branch=develop)](https://github.com/cbritopacheco/rodin/actions/workflows/Build.yml) | |

## Table of Contents

1. [Building the project](#building-the-project)
2. [Features](#features)
3. [Requirements](#requirements)
4. [CMake Options](cmake-options)
5. [Building the documentation](#building-the-documentation)


## Building the project

```
git clone --recursive https://github.com/carlos-brito-pacheco/rodin
cd rodin
mkdir build && cd build
cmake ..
make -j4
```

## Features

### Embedded DSL for FEM modelling

Rodin comes with a native C++17 domain specific language (DSL) for assembling
and solving variational formulations.

```c++
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
int main(int, char**)
{
  const char* meshFile = "../resources/mfem/poisson-example.mesh";
  int Gamma = 1;
  Mesh Omega;
  Omega.load(meshFile);
  H1 Vh(Omega);
  TrialFunction u(Vh);
  TestFunction  v(Vh);
  ScalarFunction f(1.0);
  ScalarFunction g(0.0);
  Solver::UMFPack solver;
  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, g).on(Gamma);
  poisson.solve(solver);
  return 0;
}
```

### Third-Party integrations

#### MMG

```c++
MMG::Mesh Omega;
Omega.load(meshFile, IO::FileFormat::MEDIT);

MMG::MeshOptimizer().setHMax(hmax) // maximal edge size
                    .setHMin(hmin) // minimal edge size
                    .setGradation(hgrad) // ratio between two edges
                    .setHausdorff(hausd) // curvature refinement
                    .optimize(Omega);
```

## Roadmap

List of features and modules that are in the works:
  - Discontinuous Galerkin methods
  - `Rodin::Plot` module

## Requirements

- [CMake 3.12.0+](https://cmake.org/)
- [Boost 1.65+](https://www.boost.org/)
- [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html)

Any of these should be available for quick install from your standard package
manager.

## CMake Options

| Option                 | Description                                       |
|------------------------|---------------------------------------------------|
| RODIN_BUILD_EXAMPLES   | Builds the examples in the `examples/` directory. |
| RODIN_BUILD_DOC        | Builds the documentation using Doxygen            |
| RODIN_USE_MCSS         | Builds the documentation using Doxygen and m.css  |
| RODIN_BUILD_SRC        | Build the Rodin source code                       |
| RODIN_BUILD_EXAMPLES   | Build the Rodin examples                          |
| RODIN_BUILD_DOC        | Build the Rodin documentation                     |
| RODIN_USE_MCSS         | Use m.css style documentation                     |
| RODIN_WITH_PLOT        | Build the Rodin::Plot module                      |
| RODIN_USE_MPI          | Build with MPI support                            |
| RODIN_USE_OPENMP       | Build with OpenMP support                         |
| RODIN_USE_SUITESPARSE  | Build with SuiteSparse support                    |
| RODIN_SILENCE_WARNINGS | Silence warnings outputted by Rodin               |
| RODIN_BUILD_PY         | Build Python bindings                             |

## Building the documentation

See [this page](doc/README.md) to see how to build the documentation.
