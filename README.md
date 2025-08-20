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

1. [Installation](#installation)
2. [Building the project](#building-the-project)
3. [Features](#features)
4. [Examples](#examples)
5. [Testing](#testing)
6. [Python bindings](#python-bindings)
7. [Third-Party integrations](#third-party-integrations)
8. [Requirements](#requirements)
9. [CMake options](#cmake-options)
10. [Building the documentation](#building-the-documentation)
11. [Contributing](#contributing)
12. [License](#license)
13. [Citation](#citation)


## Installation

### Building from Source

#### Prerequisites

- [CMake 3.16.0+](https://cmake.org/)
- [Boost 1.74+](https://www.boost.org/)
- A C++17 compatible compiler (GCC 7+, Clang 5+, MSVC 2017+)

#### Build Steps

```bash
git clone --recursive https://github.com/cbritopacheco/rodin
cd rodin
mkdir build && cd build
cmake ..
make -j4
```

### Python Package Installation

Rodin also provides Python bindings that can be installed via pip:

```bash
# From source
pip install .

# Development install
pip install -e .
```

**Prerequisites for Python bindings:**
- Python 3.6+
- pybind11 2.10.0+
- CMake 3.12+
- scikit-build 0.15.0+

## Building the project

For detailed build instructions, see the [Installation](#installation) section above.

Quick build:

```bash
git clone --recursive https://github.com/cbritopacheco/rodin
cd rodin
mkdir build && cd build
cmake ..
make -j4
```

## Features

### Embedded form language for FEM modelling

Rodin comes with a native C++20 form language for assembling
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

int main(int, char**)
{
  Mesh Omega;
  Omega = Omega.UniformGrid(Polytope::Type::Triangle, 16, 16);
  mesh.getConnectivity().compute(1, 2);

  P1 Vh(Omega);

  TrialFunction u(Vh);
  TestFunction v(Vh);

  Solver::SparseLU solver;

  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v))
          - Integral(v)
          + DirichletBC(u, Zero());
  poisson.solve(solver);

  return 0;
}
```

<table>
  <tr>
    <td align="center">
      <img src="doc/README/Poisson.png" alt="Poisson.png" style="width:50%;">
    </td>
  </tr>
  <tr>
    <td align="center">
      Solution of the Poisson equation.
    </td>
  </tr>
</table>

### Full high level mesh access and functionalities

#### Cell, Face, Vertex Iterators

The API offers full support for iteration over _all_ polytopes of the mesh of some given dimension:

```c++
Mesh mesh;
mesh = mesh.UniformGrid(Polytope::Type::Triangle, 16, 16); // 2D Mesh

for (auto it = mesh.getCell(); it; ++it)
{
 // Access information about the cell
}

for (auto it = mesh.getFace(); it; ++it)
{
 // Access information about the face
}

for (auto it = mesh.getVertex(); it; ++it)
{
 // Access information about the vertex
}

for (auto it = mesh.getPolytope(1); it; ++it)
{
 // Access information about the face (face dimension in 2D is equal to 1)
}


```

#### Full connectivity computation

Rodin is able to compute any connectivity information on the mesh. For example, the following computes
the adjacency information from faces to cells:

```c++
Mesh mesh;
mesh = mesh.UniformGrid(Polytope::Type::Triangle, 16, 16); // 2D Mesh

mesh.getConnectivity().compute(1, 2);
```

In general, this means that given a face, we are able to obtain the incident (neighboring) cells.

However, one can also compute any connectivity information on different dimensions.
For example, for a mesh $\mathcal{T}_h \subset \mathbb{R}^d$, $d = 2$ of topological dimension $d$, we have:

```c++
// Compute connectivity between vertices and faces
// i.e. Given a vertex, give me the incident edges
mesh.getConnectivity().compute(0, 1);

// Compute connectivity between faces and cells
// i.e. Given a vertex, give me the incident cells
mesh.getConnectivity().compute(0, 2); 

// Compute connectivity between faces
// i.e. Given a face, give me the adjacent faces
mesh.getConnectivity().compute(1, 1);

// Compute connectivity between cells
// i.e. Given a cell, give me the adjacent cells
mesh.getConnectivity().compute(2, 2);

// Compute connectivity between cells and faces
// i.e. Given a cell, give me the adjacent faces
mesh.getConnectivity().compute(2, 1);

// Etc.
```

### Direct integration with Eigen solvers

### Detailed documentation

Rodin comes with comprehensive documentation that is built automatically on each merge, ensuring it's always up to date. The documentation includes:

- **API Reference**: Complete C++ API documentation with examples
- **User Guide**: Step-by-step tutorials and usage patterns  
- **Examples**: Fully documented example programs
- **Mathematical Foundations**: Theory behind the finite element methods

[The documentation can be found here](https://cbritopacheco.github.io/rodin/).

The documentation is built using Doxygen with optional m.css styling for a modern appearance.

### Support for different finite elements

### Support for different mesh and solution file formats

- MFEM
- MEDIT

### Different quadrature formulae

Rodin supports different kinds of quadrature.

- Grundmann-Moeller

[See here for the full list](https://cbritopacheco.github.io/rodin/docs/refs/heads/master/group___rodin_quadrature.html).

### SubMesh support

## Examples

Rodin comes with a comprehensive set of examples demonstrating various features and use cases. The examples are organized into the following categories:

- **PDEs**: Basic partial differential equations (Poisson, Helmholtz, Darcy, Elasticity)
- **Geometry**: Mesh operations, refinement, partitioning, and geometric computations
- **Variational**: Finite element spaces (P0, P1), projections, and form assembly
- **Shape Optimization**: Level set methods for cantilever, arch, and eigenvalue problems
- **Boundary Optimization**: Acoustic cloaking, surface cooling, and water tank problems
- **IO**: Mesh input/output in various formats (MFEM, MEDIT, Gmsh)
- **MMG**: Mesh adaptation and optimization using MMG
- **MPI**: Parallel computing examples
- **PETSc**: Integration with PETSc solvers
- **Plotting**: Visualization examples

To build and run examples:

```bash
mkdir build && cd build
cmake .. -DRODIN_BUILD_EXAMPLES=ON
make -j4
./examples/PDEs/Poisson  # Run the Poisson example
```

All examples are located in the `examples/` directory with comprehensive documentation.

## Testing

Rodin includes a comprehensive test suite covering unit tests, manufactured solutions, and benchmarks.

### Running Tests

```bash
# Build with tests enabled
cmake .. -DRODIN_BUILD_UNIT_TESTS=ON -DRODIN_BUILD_MANUFACTURED_TESTS=ON -DRODIN_BUILD_BENCHMARKS=ON
make -j4

# Run all tests
ctest

# Run specific test categories
ctest -R unit        # Run unit tests only
ctest -R manufactured # Run manufactured solution tests
ctest -R benchmarks  # Run benchmarks
```

### Test Categories

- **Unit Tests**: Component-level testing of individual modules
- **Manufactured Tests**: Verification against known analytical solutions
- **Benchmarks**: Performance testing and regression monitoring

Test results and coverage reports are automatically generated in CI/CD pipelines.

## Python bindings

Rodin provides Python bindings built with pybind11, allowing you to use Rodin's finite element capabilities from Python.

### Installation

```bash
pip install .
```

### Basic Usage

```python
import rodin

# Create a mesh
mesh = rodin.geometry.Mesh()

# Define geometry types
triangle = rodin.geometry.Type.Triangle
square = rodin.geometry.Type.Square

# Save mesh in different formats
mesh.save("output.mesh", rodin.io.FileFormat.MFEM)
```

The Python bindings expose core Rodin functionality including:
- Mesh operations and geometry handling
- File I/O in multiple formats
- Basic finite element operations

For more examples, see the `py/` directory and the Python test suite.

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
  MMG::Optimizer().setHMax(hmax) // maximal edge size
                  .setHMin(hmin) // minimal edge size
                  .setGradation(hgrad) // ratio between two edges
                  .setHausdorff(hausd) // curvature refinement
                  .optimize(Omega);
  ```

## Roadmap

List of features and modules that are in the works:
  - Discontinuous Galerkin methods
  - `Rodin::Plot` module
  - H1
  - L2
  - HDiv
  - HCurl
  - P2
  - P0
  - PETSc
  - METIS

## Requirements

### Core Dependencies

- **CMake 3.16.0+**: Build system
- **Boost 1.74+**: C++ libraries for various utilities
- **C++17 compatible compiler**: 
  - GCC 7.0+
  - Clang 5.0+
  - MSVC 2017+

### Optional Dependencies

The following dependencies are optional and enable additional features:

- **MPI**: For parallel computing support (`-DRODIN_USE_MPI=ON`)
- **OpenMP**: For shared-memory parallelization (`-DRODIN_USE_OPENMP=ON`)
- **SuiteSparse**: For additional sparse linear algebra support (`-DRODIN_USE_SUITESPARSE=ON`)
- **PETSc**: For advanced linear algebra and parallel solvers
- **MMG**: For mesh adaptation and optimization
- **Doxygen 1.8.15+**: For building documentation (`-DRODIN_BUILD_DOC=ON`)
- **Python 3.6+**: For Python bindings (`-DRODIN_BUILD_PY=ON`)
- **pybind11 2.10.0+**: For Python bindings
- **LaTeX**: For styled documentation with m.css (`-DRODIN_USE_MCSS=ON`)

### Platform Support

Rodin has been tested on:
- Linux (Ubuntu 18.04+, CentOS 7+)
- macOS (10.14+)
- Windows (MSVC 2017+)

All dependencies should be available from standard package managers (apt, yum, brew, vcpkg).

## CMake options

| Option                 | Description                                       |
|------------------------|---------------------------------------------------|
| RODIN_BUILD_SRC        | Build the Rodin source code                      |
| RODIN_BUILD_EXAMPLES   | Build the examples in the `examples/` directory  |
| RODIN_BUILD_DOC        | Build the documentation using Doxygen            |
| RODIN_USE_MCSS         | Use m.css style documentation                    |
| RODIN_BUILD_UNIT_TESTS | Build unit tests                                 |
| RODIN_BUILD_MANUFACTURED_TESTS | Build manufactured solution tests      |
| RODIN_BUILD_BENCHMARKS | Build benchmark tests                            |
| RODIN_WITH_PLOT        | Build the Rodin::Plot module                     |
| RODIN_USE_MPI          | Build with MPI support                           |
| RODIN_USE_OPENMP       | Build with OpenMP support                        |
| RODIN_USE_SUITESPARSE  | Build with SuiteSparse support                   |
| RODIN_SILENCE_WARNINGS | Silence warnings outputted by Rodin              |
| RODIN_BUILD_PY         | Build Python bindings                            |

## Building the documentation

See [this page](doc/README.md) to see how to build the documentation.

## Contributing

Contributions to Rodin are warmly welcomed! Whether you're fixing bugs, adding features, improving documentation, or providing examples, your help is appreciated.

### Getting Started

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/amazing-feature`
3. Make your changes
4. Build and test your changes: `cmake .. && make -j4 && ctest`
5. Commit your changes: `git commit -m 'Add amazing feature'`
6. Push to the branch: `git push origin feature/amazing-feature`
7. Open a Pull Request

### Development Guidelines

- Follow the existing code style and conventions
- Add tests for new functionality
- Update documentation as needed
- Ensure all tests pass before submitting
- Write clear commit messages

### Reporting Issues

Please use the [GitHub issue tracker](https://github.com/cbritopacheco/rodin/issues) to report bugs or request features.

## License

Rodin is licensed under the **Boost Software License - Version 1.0**.

This license allows for both commercial and non-commercial use, modification, and distribution of the software. The license is very permissive and compatible with most other licenses.

For the full license text, see the [LICENSE](LICENSE) file in the repository root.

## Citation

If you use Rodin in your research or projects, please cite it as:

```bibtex
@software{brito_pacheco_2023_rodin,
  author       = {Brito-Pacheco, Carlos},
  title        = {Rodin: C++17 finite element method and shape optimization framework},
  version      = {0.0.1},
  year         = {2023},
  url          = {https://github.com/cbritopacheco/rodin},
  license      = {Boost Software License - Version 1.0}
}
```

You can also find the citation information in [CITATION.cff](CITATION.cff) format in the repository.
