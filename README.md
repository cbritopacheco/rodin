# Rodin [![License](https://img.shields.io/badge/license-BSL--1.0-green)](https://github.com/cbritopacheco/rodin/blob/master/LICENSE)

Rodin is a lightweight and modular finite element framework which provides many of the associated functionalities that are needed when implementing shape and topology optimization algorithms. These functionalities range from refining and remeshing the underlying shape, to providing elegant mechanisms to specify and solve variational problems.

It is named after the French sculptor Auguste Rodin, considered the founder of modern sculpture.

The library is still in development. It is primarily maintained by [Carlos Brito-Pacheco](https://edp-ljk.imag.fr/author/carlos-brito-pacheco/) and was developed to generate examples for his ongoing PhD.

Any contributors are warmly encouraged and any help or comments are always appreciated!

## Status

| Branch      |  Matrix  | Tests | Code Coverage | Benchmarks | Documentation |
|:-----------:|:--------:|:-----:|:-------------:|:----------:|:-------------:|
| master      | [![Build](https://github.com/cbritopacheco/rodin/actions/workflows/Build.yml/badge.svg?branch=master)](https://github.com/cbritopacheco/rodin/actions/workflows/Build.yml?query=branch%3Amaster) | [![Tests](https://github.com/cbritopacheco/rodin/actions/workflows/Tests.yml/badge.svg?branch=master)](https://github.com/cbritopacheco/rodin/actions/workflows/Tests.yml?query=branch%3Amaster) | [![codecov](https://codecov.io/gh/cbritopacheco/rodin/branch/master/graph/badge.svg?token=gwEZOnQje1)](https://app.codecov.io/gh/cbritopacheco/rodin/tree/master)  | [![Benchmarks](https://github.com/cbritopacheco/rodin/actions/workflows/Benchmarks.yml/badge.svg?branch=master)](https://cbritopacheco.github.io/rodin/benchmarks/refs/heads/master/) | [![Documentation](https://github.com/cbritopacheco/rodin/actions/workflows/Documentation.yml/badge.svg?branch=master)](https://cbritopacheco.github.io/rodin/docs/refs/heads/master) |
| develop     | [![Build](https://github.com/cbritopacheco/rodin/actions/workflows/Build.yml/badge.svg?branch=develop)](https://github.com/cbritopacheco/rodin/actions/workflows/Build.yml?query=branch%3Adevelop) | [![Tests](https://github.com/cbritopacheco/rodin/actions/workflows/Tests.yml/badge.svg?branch=develop)](https://github.com/cbritopacheco/rodin/actions/workflows/Tests.yml?query=branch%3Adevelop) | [![codecov](https://codecov.io/gh/cbritopacheco/rodin/branch/develop/graph/badge.svg?token=gwEZOnQje1)](https://app.codecov.io/gh/cbritopacheco/rodin/tree/develop) | [![Benchmarks](https://github.com/cbritopacheco/rodin/actions/workflows/Benchmarks.yml/badge.svg?branch=develop)](https://cbritopacheco.github.io/rodin/benchmarks/refs/heads/develop/) | [![Documentation](https://github.com/cbritopacheco/rodin/actions/workflows/Documentation.yml/badge.svg?branch=develop)](https://cbritopacheco.github.io/rodin/docs/refs/heads/develop) |

## Table of Contents

1. [Building the project](#building-the-project)
2. [Features](#features)
3. [Third-Party integrations](#third-party-integrations)
5. [Requirements](#requirements)
6. [CMake options](#cmake-options)
7. [Building the documentation](#building-the-documentation)


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

For example, given a domain $\Omega$ with boundary $\Gamma := \partial \Omega$, the Poisson problem:
```math
\left\{
\begin{aligned}
 -\Delta u &= f && \text{in } \Omega\\
 u &= 0 && \text{on } \Gamma \ ,
\end{aligned}
\right.
```
has the associated weak formulation:
```math
\text{Find} \ u \in H^1(\Omega) \quad \text{s.t.} \quad \forall v \in H^1_0(\Omega), \quad \int_\Omega \nabla u \cdot \nabla v \ dx = \int_\Omega f v \ dx, \quad \text{with } \quad H^1_0(\Omega) := \{ v \in H^1(\Omega) \mid v = 0 \text{ on } \Gamma \}
```

which can be quickly implemented via the following lines of code:

```c++
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
const Geometry::Attribute Gamma = 1;
int main(int, char**)
{
  Mesh Omega;
  Omega = Omega.UniformGrid(Polytope::Type::Triangle, 16, 16);
  mesh.getConnectivity().compute(1, 2);
  P1 Vh(Omega);
  TrialFunction u(Vh); TestFunction v(Vh);
  ScalarFunction f(1.0); ScalarFunction g(0.0);
  Solver::SparseLU solver;
  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, g).on(Gamma);
  poisson.solve(solver);
  return 0;
}
```

## Third-Party integrations

### MMG

[MMG](https://github.com/MmgTools/mmg) is an open source software for bidimensional and tridimensional surface and volume remeshing.

- Loading the mesh:
  ```c++
  MMG::Mesh Omega;
  Omega.load(meshFile, IO::FileFormat::MEDIT);
  ```

- Optimizing the mesh:
  ```c++
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

- [CMake 3.16.0+](https://cmake.org/)
- [Boost 1.65+](https://www.boost.org/)

Any of these should be available for quick install from your standard package
manager.

## CMake options

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
