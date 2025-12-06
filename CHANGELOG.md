# Changelog

All notable changes to the Rodin project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.0.1] - 2024-11-11

### Overview

First official release of **Rodin**, a lightweight and modular finite element framework for shape and topology optimization.

### Features

#### Core Framework
- **C++20 Form Language**: Native embedded domain-specific language for expressing and solving variational formulations
- **Mesh Management**: Full mesh access with iterators and connectivity computation
- **Context System**: Support for both local and distributed computing contexts

#### Finite Element Support
- **P0 Elements**: Piecewise constant finite elements
- **P1 Elements**: Piecewise linear finite elements on simplices
- **Trial and Test Functions**: Intuitive API for defining variational problems

#### Geometry Module
- **Mesh Generation**: Uniform grid generation for triangles, tetrahedra, and other polytopes
- **Mesh I/O**: Support for MFEM, MEDIT, and EnSight6 file formats
- **Connectivity**: Automatic computation of mesh connectivity for any dimension
- **Attributes**: Flexible attribute management for mesh entities

#### Variational Module
- **Operators**: Gradient, divergence, curl, and other differential operators
- **Integration**: Support for volume and surface integrals
- **Boundary Conditions**: Dirichlet and Neumann boundary conditions
- **Problem Assembly**: Automatic assembly of variational problems into linear systems

#### Solvers
- **Direct Solvers**: Integration with Eigen's sparse direct solvers (SparseLU, SimplicialLLT, SimplicialLDLT)
- **Optional Solvers**: Support for UMFPACK, CHOLMOD, SPQR, KLU, SuperLU when enabled
- **PETSc Integration**: Optional PETSc support for advanced solvers (requires building from source)

#### Quadrature
- **Grundmann-Moeller**: High-order quadrature formulas for simplices
- **Gauss-Legendre**: Standard Gauss quadrature rules
- **Custom Quadrature**: Support for user-defined quadrature formulas

#### I/O Support
- **MFEM Format**: Read and write MFEM mesh files
- **MEDIT Format**: Support for MEDIT mesh and solution files
- **EnSight6 Format**: Export for visualization in EnSight

#### Third-Party Integration
- **MMG**: Mesh adaptation and remeshing (requires building from source)
- **Scotch**: Graph partitioning (optional)
- **Boost**: Serialization and filesystem utilities
- **Eigen**: Linear algebra backend

#### Threading and Performance
- **OpenMP**: Optional multithreading support
- **Link-Time Optimization**: Optional LTO for performance
- **Thread-Safe**: Thread-safety options for parallel computations

#### Testing and Validation
- **Unit Tests**: Comprehensive unit test suite using GoogleTest
- **Manufactured Solutions**: Verification tests using manufactured solutions
- **Benchmarks**: Performance benchmarking infrastructure

### Installation

#### System Requirements
- CMake 3.16.0 or later
- C++20 compatible compiler (GCC 10+, Clang 10+, or MSVC 2019+)
- Boost 1.74 or later
- Eigen3 3.4 or later

#### Quick Install
```bash
git clone --recursive https://github.com/cbritopacheco/rodin.git
cd rodin
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
make -j4
sudo make install
```

For detailed installation instructions, see [INSTALL.md](INSTALL.md).

### CMake Configuration Options

#### Build Options
- `RODIN_BUILD_SRC=ON` - Build the Rodin source code (default: ON)
- `RODIN_BUILD_EXAMPLES=ON` - Build examples (default: ON)
- `RODIN_BUILD_UNIT_TESTS=OFF` - Build unit tests (default: OFF)
- `RODIN_BUILD_MANUFACTURED_TESTS=OFF` - Build manufactured solution tests (default: OFF)
- `RODIN_BUILD_BENCHMARKS=OFF` - Build benchmarks (default: OFF)
- `RODIN_BUILD_DOC=OFF` - Build documentation (default: OFF)

#### Feature Options
- `RODIN_USE_OPENMP=OFF` - Enable OpenMP support (default: OFF, auto-enabled with MULTITHREADED)
- `RODIN_MULTITHREADED=ON` - Enable multithreading capabilities (default: ON)
- `RODIN_THREAD_SAFE=OFF` - Enforce thread-safety (default: OFF, auto-enabled with MULTITHREADED)
- `RODIN_WITH_PLOT=OFF` - Build Rodin::Plot module (default: OFF)
- `RODIN_WITH_PY=OFF` - Build Python bindings (default: OFF)

#### Solver Options
- `RODIN_USE_MPI=OFF` - Enable MPI support (default: OFF)
- `RODIN_USE_PETSC=OFF` - Enable PETSc support (default: OFF)
- `RODIN_USE_UMFPACK=OFF` - Enable UMFPACK solver (default: OFF)
- `RODIN_USE_CHOLMOD=OFF` - Enable CHOLMOD solver (default: OFF)
- `RODIN_USE_SPQR=OFF` - Enable SPQR solver (default: OFF)
- `RODIN_USE_KLU=OFF` - Enable KLU solver (default: OFF)
- `RODIN_USE_SUPERLU=OFF` - Enable SuperLU solver (default: OFF)
- `RODIN_USE_SCOTCH=OFF` - Enable Scotch support (default: OFF)

### Known Limitations

1. **RodinExternal::MMG Not Exported**: The MMG mesh adaptation library is not available in the installed version. Users who need MMG support must build Rodin from source. This is intentional as MMG is built from source as a third-party dependency.

2. **Python Bindings**: Python bindings are experimental and require building from source with `RODIN_WITH_PY=ON`.

3. **Plotting Module**: The `Rodin::Plot` module requires SDL2, GLFW, and other visualization dependencies and must be explicitly enabled with `RODIN_WITH_PLOT=ON`.

4. **Limited Element Types**: Currently only P0 and P1 elements are fully supported. P2 and higher-order elements are planned for future releases.

### Documentation

- **Installation Guide**: [INSTALL.md](INSTALL.md)
- **API Documentation**: https://cbritopacheco.github.io/rodin/docs/
- **Examples**: See the `examples/` directory for usage examples

### Examples Included

The release includes numerous examples demonstrating:
- Poisson equation solving
- Linear elasticity
- Shape optimization
- Boundary optimization
- Surface optimization
- Density optimization
- Surface evolution
- Integral equations
- And more...

### License

Rodin is distributed under the Boost Software License, Version 1.0.
See the [LICENSE](LICENSE) file for details.

### Acknowledgments

This library was developed as part of Carlos Brito-Pacheco's PhD research at Université Grenoble Alpes.

### Contributors

- Carlos Brito-Pacheco (@cbritopacheco) - Lead developer and maintainer

---

## Future Releases

### Planned for v0.1.0
- Additional finite element types (P2, Nedelec, Raviart-Thomas)
- Improved MPI support
- Additional examples and tutorials
- Enhanced documentation
- Performance optimizations
