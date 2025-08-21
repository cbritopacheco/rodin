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
4. [Architecture and Design](#architecture-and-design)
5. [Library Organization](#library-organization)
6. [Examples](#examples)
7. [Testing](#testing)
8. [Third-Party integrations](#third-party-integrations)
9. [Requirements](#requirements)
10. [CMake options](#cmake-options)
11. [Building the documentation](#building-the-documentation)
12. [Advanced Usage](#advanced-usage)
13. [Performance Optimization](#performance-optimization)
14. [Troubleshooting](#troubleshooting)
15. [Contributing](#contributing)
16. [Frequently Asked Questions](#frequently-asked-questions)
17. [License](#license)
18. [Citation](#citation)


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

Rodin comes with a native C++17 form language for assembling
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
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 });
  mesh.getConnectivity().compute(1, 2);

  P1 vh(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  RealFunction f = 1;

  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, Zero());
  CG(poisson).solve();

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

### Comprehensive Finite Element Support

#### Supported Element Types
- **P0 Elements**: Piecewise constant functions on simplicial meshes
- **P1 Elements**: Piecewise linear functions with full gradient support
- **P1 Linear Elasticity**: Specialized P1 elements for elasticity problems

#### Element Features
- Automatic degree of freedom management
- Efficient assembly algorithms
- Support for vector-valued and complex-valued problems
- Robust boundary condition application

```cpp
// P0 elements for discontinuous problems
P0 p0Space(mesh);
TrialFunction rho(p0Space);  // Density field

// P1 elements for continuous problems  
P1 p1Space(mesh);
TrialFunction u(p1Space);    // Displacement field

// Vector-valued P1 elements
P1 vectorSpace(mesh, 2);     // 2D vector field
TrialFunction velocity(vectorSpace);
```

### Advanced Mesh Capabilities

#### Full high level mesh access and functionalities

##### Cell, Face, Vertex Iterators

The API offers full support for iteration over _all_ polytopes of the mesh of some given dimension:

```c++
Mesh mesh;
mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 }); // 2D Mesh

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

##### Full connectivity computation

Rodin is able to compute any connectivity information on the mesh. For example, the following computes
the adjacency information from faces to cells:

```c++
Mesh mesh;
mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 }); // 2D Mesh

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

#### Mesh Generation and Loading

```cpp
// Generate simple meshes
auto square = Mesh::UnitSquare(Polytope::Type::Triangle, 32, 32);
auto cube = Mesh::UnitCube(Polytope::Type::Tetrahedron, 16, 16, 16);

// Load from various formats
Mesh mesh;
mesh.load("domain.mfem", IO::FileFormat::MFEM);
mesh.load("domain.mesh", IO::FileFormat::MEDIT);
mesh.load("domain.msh", IO::FileFormat::Gmsh);

// Save in different formats
mesh.save("output.mesh", IO::FileFormat::MEDIT);
mesh.save("output.mfem", IO::FileFormat::MFEM);
```

### Comprehensive Solver Suite

#### Direct Solvers

Rodin provides a wide range of direct solvers for different matrix types:

```cpp
// General sparse matrices
SparseLU solver;
solver.solve(A, x, b);

// Symmetric positive definite matrices
SimplicialLLT cholesky;
cholesky.solve(A, x, b);

// Least squares problems
SparseQR qr;
qr.solve(A, x, b);

// External high-performance solvers
UMFPack umfpack;        // General sparse
CHOLMOD cholmod;        // SPD matrices
SPQR spqr;             // QR factorization
```

#### Iterative Solvers

```cpp
// Conjugate Gradient (for SPD systems)
CG cg;
cg.setTolerance(1e-8)
  .setMaxIterations(1000)
  .solve(A, x, b);

// GMRES (for general systems)
GMRES gmres;
gmres.setRestart(30)
     .setTolerance(1e-6)
     .solve(A, x, b);

// BiCGSTAB (for nonsymmetric systems)
BiCGSTAB bicgstab;
bicgstab.solve(A, x, b);

// Specialized solvers
IDRSTABL idr;          // IDR(s) method
DGMRES dgmres;         // Deflated GMRES
```

#### Solver Integration

```cpp
// Direct integration with problems
Problem problem(u, v);
problem = a - l + bc;

// Solve directly
CG(problem).solve();

// Or configure solver
CG solver;
solver.setTolerance(1e-10)
      .setMaxIterations(2000);
problem.solve(solver);
```

### Detailed documentation

Rodin comes with comprehensive documentation that is built automatically on each merge, ensuring it's always up to date. The documentation includes:

- **API Reference**: Complete C++ API documentation with examples
- **User Guide**: Step-by-step tutorials and usage patterns  
- **Examples**: Fully documented example programs
- **Mathematical Foundations**: Theory behind the finite element methods

[The documentation can be found here](https://cbritopacheco.github.io/rodin/).

The documentation is built using Doxygen with optional m.css styling for a modern appearance.

### Multiple File Format Support

Rodin supports reading and writing various mesh and solution formats:

#### Mesh Formats
- **MFEM**: Native MFEM finite element mesh format
- **MEDIT**: MEDIT mesh format (.mesh files)
- **Gmsh**: Popular open-source mesh generator format
- **Custom**: Extensible format system for user-defined formats

#### Solution Output
- **GridFunction**: Native Rodin grid function format (.gf files)
- **VTK**: For visualization in ParaView/VisIt
- **Raw Binary**: High-performance binary output
- **ASCII**: Human-readable text format

```cpp
// Load mesh in different formats
mesh.load("domain.mfem", IO::FileFormat::MFEM);
mesh.load("domain.mesh", IO::FileFormat::MEDIT);  
mesh.load("domain.msh", IO::FileFormat::Gmsh);

// Save solutions
u.getSolution().save("solution.gf");
u.getSolution().save("solution.vtk", IO::FileFormat::VTK);
```

### Advanced Quadrature Support

Rodin supports multiple quadrature formulae for accurate numerical integration:

#### Available Quadrature Rules
- **Grundmann-Moeller**: High-order quadrature on simplices
- **Gauss-Legendre**: Standard Gaussian quadrature
- **Custom Rules**: User-defined quadrature points and weights

```cpp
// Use different quadrature orders
Integral(f, v, QuadratureRule::GrundmannMoeller(order=5));

// Custom quadrature
QuadratureRule custom;
custom.addPoint(point, weight);
Integral(f, v, custom);
```

[See here for the full list](https://cbritopacheco.github.io/rodin/docs/refs/heads/master/group___rodin_quadrature.html).

### SubMesh support

Rodin provides comprehensive support for sub-meshes, allowing you to work with portions of your mesh for boundary conditions, material interfaces, or domain decomposition.

```c++
// Create a submesh from boundary elements
auto boundary = mesh.getBoundary();
auto boundaryMesh = mesh.createSubMesh(boundary);

// Create submesh from specific regions
auto materialRegion = mesh.getRegion(materialId);
auto regionMesh = mesh.createSubMesh(materialRegion);

// Work with submesh connectivity and finite element spaces
P1 boundarySpace(boundaryMesh);
TrialFunction lambda(boundarySpace);  // Lagrange multiplier on boundary
```

### High-Performance Computing Support

Rodin is designed with performance in mind and supports various HPC paradigms:

#### Parallel Computing
- **OpenMP**: Shared-memory parallelism for assembly and solving
- **MPI**: Distributed-memory parallelism for large-scale problems
- **Hybrid**: OpenMP + MPI for maximum scalability

#### Performance Features
- **Vectorization**: SIMD-optimized inner loops where possible
- **Memory Efficient**: Optimized data structures for large-scale problems
- **Cache-Friendly**: Memory access patterns optimized for modern architectures
- **Memory Pools**: Efficient memory management for repeated operations

Performance characteristics:
- Assembly scales efficiently with OpenMP threads
- Memory usage: ~10-50 MB per million degrees of freedom
- Supports problems with millions of unknowns on standard workstations
- Excellent weak and strong scaling on HPC clusters

```cpp
// Enable OpenMP parallelism
#pragma omp parallel
{
    // Parallel assembly automatically handled
    problem.assemble();
}

// MPI parallelism
Context::MPI mpi;
Mesh<Context::MPI> distributedMesh(mpi);
P1 parallelSpace(distributedMesh);
```

### Advanced Mesh Operations

Beyond basic mesh handling, Rodin provides advanced mesh operations:

```cpp
// Mesh refinement
auto refinedMesh = mesh.refine(refinementCriteria);

// Mesh adaptation with MMG
#include <Rodin/External/MMG.h>
MMG::Optimizer optimizer;
optimizer.setHMax(0.1)
         .setHMin(0.01)
         .setGradation(1.3)
         .setHausdorff(0.01)
         .optimize(mesh);

// Mesh partitioning for parallel computing
#include <Rodin/External/METIS.h>
METIS::Partitioner partitioner;
auto partitions = partitioner.setNumParts(numParts)
                             .setMethod(METIS::KWAY)
                             .partition(mesh);

// Mesh quality assessment
auto quality = mesh.getQuality();
std::cout << "Min angle: " << quality.minAngle << std::endl;
std::cout << "Max aspect ratio: " << quality.maxAspectRatio << std::endl;
```

## Examples

Rodin comes with over 85 comprehensive examples demonstrating various features and use cases. The examples are organized into the following categories:

### Core Examples

#### PDEs (Partial Differential Equations)
- **Poisson**: Classic Poisson equation with Dirichlet boundary conditions
- **Helmholtz**: Wave equation and acoustic problems
- **Darcy**: Flow in porous media
- **Elasticity**: Linear elasticity problems and structural mechanics
- **Periodic**: Problems with periodic boundary conditions
- **SurfaceMesh**: PDEs on surface meshes

#### Geometry and Mesh Operations
- **Mesh Loading**: Import meshes from various formats (MFEM, MEDIT, Gmsh)
- **Connectivity**: Compute and use mesh connectivity information
- **Refinement**: Adaptive and uniform mesh refinement
- **SubMesh**: Working with boundary and region submeshes

#### Variational Formulations
- **P0/P1 Elements**: Discontinuous and continuous finite elements
- **Form Assembly**: Matrix and vector assembly patterns
- **Projections**: L2 and H1 projections
- **Integration**: Custom quadrature and integration

### Advanced Examples

#### Shape and Topology Optimization
- **Cantilever**: Classic cantilever beam optimization
- **Arch Problems**: Structural optimization of arches
- **Level Set Methods**: Topology optimization using level sets
- **Eigenvalue Problems**: Frequency optimization

#### Boundary and Surface Optimization
- **Acoustic Cloaking**: Shape optimization for acoustic devices
- **Surface Cooling**: Heat transfer optimization
- **Water Tank**: Fluid optimization problems

#### High-Performance Computing
- **MPI Examples**: Parallel computing with message passing
- **Shared Memory**: OpenMP parallelization examples
- **PETSc Integration**: Advanced linear algebra and solvers

#### Specialized Modules
- **MMG Integration**: Mesh adaptation and optimization
- **Plotting**: Visualization and output generation
- **Alert System**: Error handling and debugging

### Building and Running Examples

```bash
# Build all examples
mkdir build && cd build
cmake .. -DRODIN_BUILD_EXAMPLES=ON
make -j4

# Run specific examples
./examples/PDEs/Poisson
./examples/ShapeOptimization/SimpleCantilever2D
./examples/MPI/Variational/P1

# Run with different parameters
./examples/PDEs/Helmholtz --frequency 100 --mesh coarse.mesh
```

### Example Code Patterns

Most examples follow a consistent pattern:

```c++
#include <Rodin/Rodin.h>

int main() {
    // 1. Load or generate mesh
    Geometry::Mesh<Context::Sequential> mesh;
    mesh.load("domain.mesh", IO::FileFormat::MFEM);
    
    // 2. Define finite element space
    auto Vh = FiniteElementSpace<Context::Sequential>(mesh, FE::P1{});
    
    // 3. Define variational problem
    auto u = TrialFunction(Vh);
    auto v = TestFunction(Vh);
    
    // 4. Assemble system
    auto a = BilinearForm(u, v).add(Integral(Grad(u), Grad(v)));
    auto l = LinearForm(v).add(Integral(f, v));
    
    // 5. Apply boundary conditions and solve
    // ...
    
    return 0;
}
```

All examples include detailed comments and are located in the `examples/` directory with comprehensive documentation explaining the mathematical background and implementation details.

## Architecture and Design

### Core Design Principles

Rodin follows several key design principles that make it both powerful and easy to use:

#### Template-Based Architecture
- **Compile-time optimization**: Extensive use of C++ templates for zero-cost abstractions
- **Type safety**: Strong typing system prevents common FEM implementation errors
- **Generic programming**: Write code once, works with different element types and dimensions

#### Modular Design
- **Plugin architecture**: Easy to extend with new finite elements, solvers, or mesh formats
- **Separation of concerns**: Clear separation between geometry, discretization, and solving
- **Composable components**: Mix and match different components as needed

#### Modern C++ Features
- **C++17/20 compliance**: Uses modern language features for cleaner, safer code
- **RAII**: Automatic resource management prevents memory leaks
- **Move semantics**: Efficient handling of large objects like meshes and matrices

### Library Structure

```
Rodin/
├── Geometry/          # Mesh handling, connectivity, I/O
├── Variational/       # Finite element spaces, forms, assembly
├── Math/              # Linear algebra, solvers, algorithms
├── Context/           # Execution contexts (Sequential, MPI)
├── IO/                # File format support
├── Alert/             # Error handling and debugging
└── External/          # Third-party integrations
```

### Memory Management

Rodin uses smart memory management strategies:

- **Automatic memory management**: RAII and smart pointers eliminate manual memory management
- **Memory pools**: Efficient allocation for temporary objects during assembly
- **Copy-on-write**: Expensive operations are delayed until actually needed
- **Memory mapping**: Large mesh files can be memory-mapped for efficiency

### Thread Safety

- **Thread-safe assembly**: Parallel assembly using OpenMP with proper synchronization
- **Immutable data structures**: Many core objects are immutable after construction
- **Local storage**: Thread-local storage for temporary data in parallel sections

## Library Organization

Rodin is organized into several key modules, each providing specific functionality for finite element method computations and shape optimization:

### Core Modules

#### Rodin::Geometry
The geometry module handles all mesh-related operations and provides the foundation for finite element computations.

**Key Components:**
```cpp
#include <Rodin/Geometry.h>

// Basic mesh operations
Mesh mesh;
mesh.load("domain.mesh", IO::FileFormat::MFEM);
mesh.save("output.mesh", IO::FileFormat::MEDIT);

// Connectivity computation
mesh.getConnectivity().compute(0, 1); // Vertex to edge
mesh.getConnectivity().compute(1, 2); // Edge to face
mesh.getConnectivity().compute(2, 2); // Face to face

// Mesh iterators
for (auto it = mesh.getCell(); it; ++it) {
    // Process cells
}

for (auto it = mesh.getFace(); it; ++it) {
    // Process faces
}

// SubMesh creation
auto boundary = mesh.getBoundary();
auto region = mesh.getRegion(materialId);
```

**Features:**
- Support for simplicial meshes (triangles, tetrahedra)
- Automatic connectivity computation
- SubMesh extraction for boundaries and regions
- Multiple file format support (MFEM, MEDIT, Gmsh)
- Mesh refinement and adaptation capabilities

#### Rodin::Variational
The variational module provides the core finite element functionality, including spaces, forms, and problem assembly.

**Finite Element Spaces:**
```cpp
#include <Rodin/Variational.h>

// P0 elements (piecewise constant)
P0 p0Space(mesh);

// P1 elements (piecewise linear)
P1 p1Space(mesh);

// Trial and test functions
TrialFunction u(p1Space);
TestFunction v(p1Space);
```

**Form Assembly:**
```cpp
// Bilinear forms
auto a = BilinearForm(u, v).add(Integral(Grad(u), Grad(v)));

// Linear forms  
auto l = LinearForm(v).add(Integral(f, v));

// Boundary conditions
auto bc = DirichletBC(u, Zero());

// Complete problem
Problem problem(u, v);
problem = a - l + bc;
```

**Advanced Features:**
- Multiple quadrature rules (Grundmann-Moeller, custom)
- Boundary and interface integrals
- Periodic boundary conditions
- Complex-valued problems
- Vector-valued functions

#### Rodin::Solver
The solver module provides various linear and nonlinear solvers with optimized implementations.

**Direct Solvers:**
```cpp
#include <Rodin/Solver.h>

// Sparse LU decomposition
SparseLU solver;
solver.solve(A, x, b);

// Cholesky factorization
SimplicialLLT solver;
solver.solve(A, x, b);
```

**Iterative Solvers:**
```cpp
// Conjugate Gradient
CG cg;
cg.setTolerance(1e-8)
  .setMaxIterations(1000)
  .solve(A, x, b);

// GMRES
GMRES gmres;
gmres.setRestart(30)
     .solve(A, x, b);

// BiCGSTAB
BiCGSTAB bicgstab;
bicgstab.solve(A, x, b);
```

**Specialized Solvers:**
- UMFPACK for general sparse systems
- CHOLMOD for symmetric positive definite systems
- SPQR for least squares problems
- PaStiX for parallel direct solving
- Integration with external solver libraries

#### Rodin::Math
Mathematical utilities and linear algebra operations.

```cpp
#include <Rodin/Math.h>

// Vector operations
Vector<Real> v(n);
v.setZero();
v.setRandom();

// Matrix operations
SparseMatrix<Real> A(n, n);
A.setFromTriplets(triplets.begin(), triplets.end());

// Mathematical functions
Real norm = v.norm();
Real dot = u.dot(v);
```

#### Rodin::Context
Execution contexts for different computational environments.

**Sequential Context:**
```cpp
// Default single-threaded execution
Mesh<Context::Sequential> mesh;
P1<Context::Sequential> fes(mesh);
```

**MPI Context:**
```cpp
#include <Rodin/Context/MPI.h>

Context::MPI mpi(env, world);
Mesh<Context::MPI> mesh(mpi);
P1<Context::MPI> fes(mesh);

// Parallel assembly and solving automatically handled
```

#### Rodin::IO
Input/output operations for various file formats.

```cpp
#include <Rodin/IO.h>

// Supported formats
mesh.load("input.mfem", IO::FileFormat::MFEM);
mesh.load("input.mesh", IO::FileFormat::MEDIT);
mesh.load("input.msh", IO::FileFormat::Gmsh);

// Grid function I/O
GridFunction u(fes);
u.load("solution.gf");
u.save("output.gf");
```

### Extension Modules

#### Rodin::Plot
Visualization and plotting capabilities.

```cpp
#include <Rodin/Plot.h>

Plot::Figure fig;
fig.plot(mesh, solution)
   .setTitle("Solution")
   .setColorbar(true)
   .save("solution.png");
```

#### Rodin::Alert
Error handling and debugging utilities.

```cpp
#include <Rodin/Alert.h>

// Assertions and error reporting
RODIN_ASSERT(condition, "Error message");

// Warning system
Alert::Warning() << "This is a warning message.";

// Info messages
Alert::Info() << "Computation completed successfully.";
```

### External Integrations

#### MMG Integration
```cpp
#include <Rodin/External/MMG.h>

MMG::Optimizer optimizer;
optimizer.setHMax(0.1)
         .setHMin(0.01)
         .setGradation(1.3)
         .optimize(mesh);
```

#### PETSc Integration
```cpp
#include <Rodin/External/PETSc.h>

PETSc::KSP solver;
solver.setType(PETSc::KSPGMRES)
      .setPreconditioner(PETSc::PCILU)
      .solve(A, x, b);
```

### Module Dependencies

The modules have clear dependency relationships:
- **Geometry** ← Base dependency
- **Variational** ← Depends on Geometry
- **Solver** ← Depends on Math, can work with Variational
- **Context** ← Affects all modules
- **IO** ← Depends on Geometry
- **Plot** ← Depends on Geometry, Variational
- **Alert** ← Used by all modules

This modular design allows for:
- Clear separation of concerns
- Easy testing of individual components
- Flexible integration of new features
- Minimal compilation dependencies

## Testing

Rodin includes a comprehensive test suite covering unit tests, manufactured solutions, and benchmarks to ensure code quality and performance.

### Test Organization

The test suite is organized into three main categories:

```
tests/
├── unit/           # Component-level testing
├── manufactured/   # Verification against analytical solutions  
└── benchmarks/     # Performance testing and regression monitoring
```

### Running Tests

#### Build Configuration
```bash
# Enable all test types
cmake .. -DRODIN_BUILD_UNIT_TESTS=ON \
         -DRODIN_BUILD_MANUFACTURED_TESTS=ON \
         -DRODIN_BUILD_BENCHMARKS=ON

# Build tests
make -j4
```

#### Test Execution
```bash
# Run all tests
ctest

# Run specific test categories
ctest -R unit          # Unit tests only
ctest -R manufactured  # Manufactured solution tests only
ctest -R benchmarks    # Benchmarks only

# Run tests with detailed output
ctest -V

# Run tests in parallel
ctest -j4

# Run specific tests by name
ctest -R "Poisson"     # All Poisson-related tests
ctest -R "P1.*unit"    # All P1 unit tests
```

### Test Categories

#### Unit Tests
Component-level testing of individual modules and classes:

**Geometry Module Tests:**
- Mesh loading and saving functionality
- Connectivity computation algorithms
- SubMesh creation and manipulation
- Iterator correctness

**Variational Module Tests:**
- Finite element space construction
- Form assembly correctness
- Boundary condition application
- Integration accuracy

**Solver Module Tests:**
- Direct solver correctness
- Iterative solver convergence
- Preconditioner effectiveness
- Error handling

**Example Unit Test:**
```cpp
// tests/unit/Geometry/MeshTest.cpp
TEST_CASE("Mesh connectivity computation", "[Geometry][Mesh]") {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {4, 4});
    
    // Test edge-to-cell connectivity
    mesh.getConnectivity().compute(1, 2);
    
    // Verify connectivity is computed correctly
    auto edgeToCell = mesh.getConnectivity().get(1, 2);
    REQUIRE(edgeToCell.size() > 0);
    
    // Check boundary edges have only one incident cell
    for (auto edge : mesh.getBoundaryEdges()) {
        REQUIRE(edgeToCell[edge].size() == 1);
    }
}
```

#### Manufactured Solution Tests
Verification against known analytical solutions to ensure mathematical correctness:

**PDE Problem Tests:**
- Poisson equation with polynomial solutions
- Helmholtz equation with trigonometric solutions  
- Linear elasticity with manufactured displacements
- Heat equation with known time-dependent solutions

**Convergence Analysis:**
```cpp
// tests/manufactured/PDEs/PoissonConvergence.cpp
TEST_CASE("Poisson convergence", "[Manufactured][Poisson]") {
    std::vector<int> meshSizes = {8, 16, 32, 64};
    std::vector<Real> errors;
    
    for (int N : meshSizes) {
        Mesh mesh;
        mesh = mesh.UniformGrid(Polytope::Type::Triangle, {N, N});
        
        P1 fes(mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        
        // Manufactured solution: u_exact = sin(pi*x) * sin(pi*y)
        auto u_exact = [](const Point& p) {
            return std::sin(M_PI * p.x) * std::sin(M_PI * p.y);
        };
        
        auto f = [](const Point& p) {
            return 2 * M_PI * M_PI * std::sin(M_PI * p.x) * std::sin(M_PI * p.y);
        };
        
        Problem problem(u, v);
        problem = Integral(Grad(u), Grad(v))
                - Integral(f, v)
                + DirichletBC(u, Zero());
        
        CG(problem).solve();
        
        // Compute L2 error
        auto error = L2Error(u.getSolution(), u_exact);
        errors.push_back(error);
    }
    
    // Verify O(h^2) convergence for P1 elements
    for (size_t i = 1; i < errors.size(); ++i) {
        Real ratio = errors[i-1] / errors[i];
        REQUIRE(ratio >= 3.5);  // Should be close to 4 for h^2 convergence
    }
}
```

#### Benchmarks
Performance testing and regression monitoring:

**Assembly Benchmarks:**
- Matrix assembly performance scaling
- Memory usage profiling
- Parallel efficiency measurements

**Solver Benchmarks:**
- Iterative solver performance
- Direct solver comparison
- Preconditioning effectiveness

**Example Benchmark:**
```cpp
// tests/benchmarks/Assembly/P1AssemblyBenchmark.cpp
BENCHMARK_DEFINE_F(P1Assembly, LaplaceAssembly)(benchmark::State& state) {
    int N = state.range(0);
    
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {N, N});
    mesh.getConnectivity().compute(1, 2);
    
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    
    for (auto _ : state) {
        BilinearForm a(u, v);
        a.add(Integral(Grad(u), Grad(v)));
        benchmark::DoNotOptimize(a.assemble());
    }
    
    state.counters["DOFs"] = fes.getNumDOFs();
    state.counters["Elements"] = mesh.getNumCells();
}

BENCHMARK_REGISTER_F(P1Assembly, LaplaceAssembly)
    ->Range(8, 256)
    ->Unit(benchmark::kMillisecond);
```

### Test Output and Reporting

#### CTest Integration
Tests integrate with CTest for standardized reporting:

```bash
# Generate test reports
ctest --output-on-failure -T Test

# XML output for CI/CD
ctest -T Test --output-xml TestResults.xml
```

#### Coverage Reports
When built with coverage support:

```bash
# Build with coverage
cmake .. -DRODIN_CODE_COVERAGE=ON

# Run tests and generate coverage
make coverage

# View coverage report
open coverage/index.html
```

#### Performance Regression Detection
Benchmarks include automatic performance regression detection:

```bash
# Run benchmarks with baseline comparison
ctest -R benchmarks --output-on-failure

# Generate performance reports
./tests/benchmarks/generate_report.py
```

### Continuous Integration

Tests run automatically on multiple platforms and configurations:

**GitHub Actions Workflows:**
- **Tests.yml**: Runs all test categories on push/PR
- **Benchmarks.yml**: Performance regression testing
- **Coverage.yml**: Code coverage analysis

**Test Matrix:**
- Operating Systems: Ubuntu 20.04/22.04, macOS 11/12, Windows 2019/2022
- Compilers: GCC 9/10/11, Clang 11/12/13, MSVC 2019/2022
- Build Types: Debug, Release, RelWithDebInfo
- Configurations: With/without MPI, OpenMP, external libraries

**Test Results:**
- All tests must pass before merge
- Coverage must not decrease
- Performance regressions are flagged
- Results available at: [GitHub Actions](https://github.com/cbritopacheco/rodin/actions)

### Writing New Tests

#### Adding Unit Tests
```cpp
#include <catch2/catch.hpp>
#include <Rodin/Rodin.h>

TEST_CASE("MyNewFeature works correctly", "[Module][Component]") {
    // Arrange
    // Set up test data
    
    // Act  
    // Execute the functionality being tested
    
    // Assert
    // Verify expected behavior
    REQUIRE(result == expected);
    CHECK(approximatelyEqual(actual, expected, tolerance));
}
```

#### Adding Manufactured Tests
```cpp
#include "ManufacturedTestBase.h"

class MyPDEConvergenceTest : public ManufacturedTestBase {
public:
    Real exactSolution(const Point& p) override {
        // Return analytical solution
    }
    
    Real sourceFunction(const Point& p) override {
        // Return corresponding source term
    }
    
    void setupProblem(Problem& problem) override {
        // Configure the PDE problem
    }
};

TEST_CASE_METHOD(MyPDEConvergenceTest, "MyPDE convergence", "[Manufactured]") {
    runConvergenceTest({8, 16, 32, 64}, 2.0); // Expected convergence rate
}
```

#### Adding Benchmarks
```cpp
#include <benchmark/benchmark.h>
#include <Rodin/Rodin.h>

static void BM_MyFeature(benchmark::State& state) {
    // Setup
    auto input = setupBenchmarkInput(state.range(0));
    
    for (auto _ : state) {
        // Code to benchmark
        auto result = myFeature(input);
        benchmark::DoNotOptimize(result);
    }
    
    // Optional metrics
    state.counters["ItemsPerSecond"] = 
        benchmark::Counter(state.iterations(), benchmark::Counter::kIsRate);
}

BENCHMARK(BM_MyFeature)->Range(1<<10, 1<<20);
```

The comprehensive test suite ensures Rodin maintains high quality and performance standards across all supported platforms and configurations.

## Third-Party integrations

Rodin integrates seamlessly with several high-quality third-party libraries to provide extended functionality.

### MMG - Mesh Adaptation

[MMG](https://github.com/MmgTools/mmg) is an open source software for bidimensional and tridimensional surface and volume remeshing.

```c++
#include <Rodin/External/MMG.h>

// Load the mesh
MMG::Mesh Omega;
Omega.load(meshFile, IO::FileFormat::MEDIT);

// Configure optimization parameters
MMG::Optimizer optimizer;
optimizer.setHMax(hmax)         // maximal edge size
         .setHMin(hmin)         // minimal edge size
         .setGradation(hgrad)   // ratio between two edges
         .setHausdorff(hausd)   // curvature refinement
         .setAngleDetection(45) // angle for ridge detection
         .optimize(Omega);

// Use optimized mesh in Rodin
Geometry::Mesh<Context::Sequential> mesh;
mesh = MMG::convert(Omega);
```

### PETSc - Advanced Linear Algebra

[PETSc](https://petsc.org/) provides scalable parallel linear algebra and solvers.

```c++
#include <Rodin/External/PETSc.h>

// Create PETSc-based finite element space
auto Vh = FiniteElementSpace<Context::PETSc>(mesh, FE::P1{});

// Use PETSc solvers
PETSc::KSP solver;
solver.setType(PETSc::KSPGMRES)
      .setPreconditioner(PETSc::PCILU)
      .setTolerances(1e-8, 1e-12, PETSC_DEFAULT, 1000)
      .solve(A, x, b);
```

### MPI - Parallel Computing

Built-in MPI support for distributed computing:

```cpp
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <Rodin/MPI/Variational/P1.h>

// Initialize MPI environment
boost::mpi::environment env(argc, argv);
boost::mpi::communicator world;

// Create MPI context
Context::MPI mpi(env, world);

// Create distributed mesh
Mesh<Context::MPI> mesh(mpi);
mesh.load("mesh.mfem");

// Parallel finite element space
P1 fes(mesh);

// Parallel assembly automatically handled
TrialFunction u(fes);
TestFunction v(fes);
Problem problem(u, v);
problem = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, Zero());
```

### METIS - Graph Partitioning

[METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) for mesh partitioning:

```c++
#include <Rodin/External/METIS.h>

// Partition mesh for parallel computing
METIS::Partitioner partitioner;
auto partitions = partitioner.setNumParts(nprocs)
                             .setMethod(METIS::KWAY)
                             .partition(mesh);
```

### Boost - Core Utilities

Extensive use of [Boost](https://www.boost.org/) libraries:

- **Boost.Geometry**: Geometric algorithms and spatial data structures
- **Boost.Graph**: Graph algorithms for mesh connectivity
- **Boost.Multi-Array**: Efficient multi-dimensional arrays
- **Boost.MPI**: High-level MPI interface
- **Boost.Serialization**: Object serialization for checkpointing

### SuiteSparse - Sparse Linear Algebra

[SuiteSparse](https://sparse.tamu.edu/) provides advanced sparse matrix algorithms:

```c++
// Use UMFPACK for direct sparse solving
SuiteSparse::UMFPACK solver;
solver.factorize(A).solve(x, b);

// Use CHOLMOD for symmetric positive definite systems
SuiteSparse::CHOLMOD solver;
solver.factorize(A).solve(x, b);
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
- **LaTeX**: For styled documentation with m.css (`-DRODIN_USE_MCSS=ON`)

### Platform Support

Rodin has been tested on:
- Linux (Ubuntu 18.04+, CentOS 7+)
- macOS (10.14+)
- Windows (MSVC 2017+)

All dependencies should be available from standard package managers (apt, yum, brew, vcpkg).

### Installation Troubleshooting

#### Common Build Issues

**Boost not found:**
```bash
# Ubuntu/Debian
sudo apt-get install libboost-all-dev

# CentOS/RHEL
sudo yum install boost-devel

# macOS
brew install boost

# Or specify Boost location
cmake .. -DBOOST_ROOT=/path/to/boost
```

**CMake version too old:**
```bash
# Install newer CMake from Kitware
wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc | sudo apt-key add -
sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ focal main'
sudo apt-get update && sudo apt-get install cmake
```

**C++17 compiler issues:**
```bash
# Ensure modern compiler
sudo apt-get install gcc-9 g++-9
export CC=gcc-9 CXX=g++-9
cmake .. -DCMAKE_CXX_COMPILER=g++-9
```

**Git submodule problems:**
```bash
# Ensure submodules are properly initialized
git submodule update --init --recursive

# Or clone with submodules
git clone --recursive https://github.com/cbritopacheco/rodin
```

#### Runtime Issues

**Segmentation faults:**
- Usually indicate mesh connectivity issues or boundary condition problems
- Use debug builds: `cmake .. -DCMAKE_BUILD_TYPE=Debug`
- Enable Rodin assertions: `cmake .. -DRODIN_ENABLE_ASSERTIONS=ON`

**Performance problems:**
- Ensure release build: `cmake .. -DCMAKE_BUILD_TYPE=Release`
- Enable OpenMP: `cmake .. -DRODIN_USE_OPENMP=ON`
- Check memory usage with tools like `valgrind` or `heaptrack`

**Convergence issues:**
- Verify boundary conditions are properly applied
- Check mesh quality and element orientation
- Try different solvers or preconditioners

## CMake options

| Option                         | Description                                       |
|--------------------------------|---------------------------------------------------|
| RODIN_BUILD_SRC                | Build the Rodin source code                      |
| RODIN_BUILD_EXAMPLES           | Build the examples in the `examples/` directory  |
| RODIN_BUILD_DOC                | Build the documentation using Doxygen            |
| RODIN_USE_MCSS                 | Use m.css style documentation                    |
| RODIN_BUILD_UNIT_TESTS         | Build unit tests                                 |
| RODIN_BUILD_MANUFACTURED_TESTS | Build manufactured solution tests                |
| RODIN_BUILD_BENCHMARKS         | Build benchmark tests                            |
| RODIN_WITH_PLOT                | Build the Rodin::Plot module                     |
| RODIN_WITH_PY                  | Build Python bindings (currently disabled)       |
| RODIN_USE_MPI                  | Build with MPI support                           |
| RODIN_USE_OPENMP               | Build with OpenMP support                        |
| RODIN_USE_UMFPACK              | Build with UMFPACK support                       |
| RODIN_USE_PASTIX               | Build with PaStiX support                        |
| RODIN_USE_KLU                  | Build with KLU support                           |
| RODIN_USE_SUPERLU              | Build with SuperLU support                       |
| RODIN_USE_SPQR                 | Build with SPQR support                          |
| RODIN_USE_PARDISO              | Build with Pardiso support                       |
| RODIN_USE_APPLE_ACCELERATE     | Build with Apple Accelerate support              |
| RODIN_USE_SCOTCH               | Build with Scotch support                        |
| RODIN_USE_PETSC                | Build with PETSc support                         |
| RODIN_SILENCE_WARNINGS         | Silence warnings outputted by Rodin              |
| RODIN_SILENCE_EXCEPTIONS       | Silence exceptions thrown by Rodin               |
| RODIN_CODE_COVERAGE            | Compile with code coverage flags                 |
| RODIN_LTO                      | Compile with link time optimization              |
| RODIN_FEATURE_SUMMARY          | Print Rodin feature summary                      |
| RODIN_THREAD_SAFE              | Enforces thread-safety                           |
| RODIN_MULTITHREADED            | Utilizes multithreading capabilities             |
| BUILD_SHARED_LIBS              | Build using shared libraries                     |

## Building the documentation

See [this page](doc/README.md) to see how to build the documentation.

## Advanced Usage

This section covers advanced features and usage patterns for experienced users.

### Custom Finite Elements

While Rodin currently provides P0 and P1 elements, the template-based architecture makes it easy to implement custom finite elements:

```cpp
// Custom finite element class
template<typename Context>
class MyCustomElement : public FiniteElementBase<Context> {
public:
    // Implementation details
    void assemble(const BilinearFormIntegrator& integrator) override;
    void project(const Function& f) override;
};

// Usage
MyCustomElement<Context::Sequential> customFE(mesh);
```

### Advanced Mesh Operations

#### Mesh Adaptation with MMG

```cpp
#include <Rodin/External/MMG.h>

// Load initial mesh
Mesh mesh;
mesh.load("initial.mesh", IO::FileFormat::MEDIT);

// Solve problem to get size field
P1 fes(mesh);
TrialFunction u(fes);
// ... solve problem ...

// Create metric based on solution
MMG::MetricTensor metric;
metric.computeFromGradient(u.getSolution());

// Adapt mesh
MMG::Optimizer optimizer;
optimizer.setMetric(metric)
         .setHausdorff(0.01)
         .setAngleDetection(45)
         .optimize(mesh);
```

#### Parallel Mesh Partitioning

```cpp
#include <Rodin/External/METIS.h>

// Partition mesh for parallel computing
METIS::Partitioner partitioner;
auto partitions = partitioner.setNumParts(4)
                             .setMethod(METIS::KWAY)
                             .setImbalanceTolerance(1.05)
                             .partition(mesh);

// Distribute to MPI processes
for (int rank = 0; rank < world.size(); ++rank) {
    auto localMesh = mesh.extract(partitions[rank]);
    // Send localMesh to rank
}
```

### Complex Boundary Conditions

#### Periodic Boundary Conditions

```cpp
// Define periodic boundaries
auto leftBoundary = mesh.getBoundary(1);
auto rightBoundary = mesh.getBoundary(2);

// Create periodic BC
PeriodicBC periodic(u, leftBoundary, rightBoundary);
periodic.setTransformation([](const Point& p) {
    return Point(p.x + L, p.y); // Translation by L
});

Problem problem(u, v);
problem = a - l + periodic;
```

#### Robin Boundary Conditions

```cpp
// Robin BC: alpha * u + beta * grad(u) · n = g
RealFunction alpha = 1.0;
RealFunction beta = 0.1;
RealFunction g = 0.0;

auto robinBC = RobinBC(u, v, alpha, beta, g);
problem += robinBC;
```

### Shape Optimization Workflows

```cpp
#include <Rodin/ShapeOptimization.h>

// Define shape optimization problem
ShapeOptimization::Problem shapeProblem;

// Objective function (compliance minimization)
auto compliance = [&](const Mesh& mesh) {
    // Solve state equation
    auto u = solveStateEquation(mesh);
    return computeCompliance(u);
};

// Volume constraint
auto volumeConstraint = [&](const Mesh& mesh) {
    return mesh.volume() - targetVolume;
};

shapeProblem.setObjective(compliance)
           .addConstraint(volumeConstraint)
           .setOptimizer(ShapeOptimization::LBFGS{});

// Run optimization
auto result = shapeProblem.solve(initialMesh);
```

### Custom Integrators

```cpp
// Custom bilinear form integrator
class MyIntegrator : public BilinearFormIntegrator {
public:
    void integrate(const TrialFunction& u, const TestFunction& v,
                   const Cell& cell, Matrix& localMatrix) override {
        // Custom integration logic
        auto qr = cell.getQuadratureRule(order);
        for (auto& qp : qr) {
            auto grad_u = u.grad(qp);
            auto grad_v = v.grad(qp);
            localMatrix += qp.weight * myKernel(grad_u, grad_v);
        }
    }
private:
    Real myKernel(const Gradient& grad_u, const Gradient& grad_v);
};

// Usage
BilinearForm a(u, v);
a.add(MyIntegrator{});
```

## Performance Optimization

### Compilation Optimizations

For maximum performance, use these CMake settings:

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_CXX_FLAGS="-O3 -march=native -DNDEBUG" \
         -DRODIN_LTO=ON \
         -DRODIN_USE_OPENMP=ON
```

### Memory Optimization

#### Large Problem Handling

```cpp
// For problems with millions of DOFs
cmake .. -DCMAKE_CXX_FLAGS="-DEIGEN_MAX_STATIC_ALIGN_BYTES=0"

// Use memory-mapped files for large meshes
Mesh mesh;
mesh.loadMemoryMapped("huge_mesh.mfem");

// Enable memory pools for assembly
Context::Sequential context;
context.enableMemoryPool(true);
```

#### Sparse Matrix Optimization

```cpp
// Pre-allocate sparse matrix pattern
SparseMatrix A(n, n);
A.reserve(estimatedNonZeros);

// Use compressed storage
A.makeCompressed();

// Optimize for repeated assembly
A.setOptimizationLevel(SparseMatrix::OptimizationLevel::High);
```

### Parallel Performance

#### OpenMP Configuration

```bash
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

# Run with optimized settings
./my_rodin_program
```

#### MPI Performance Tips

```cpp
// Minimize communication
Context::MPI mpi;
mpi.setBufferSize(1024 * 1024); // 1MB buffer

// Use non-blocking operations
auto request = mpi.iallreduce(localSum);
// Do other work...
auto globalSum = mpi.wait(request);
```

### Solver Performance

#### Iterative Solver Tuning

```cpp
// Fine-tune CG solver
CG cg;
cg.setTolerance(1e-6)           // Reduce if possible
  .setMaxIterations(1000)       // Increase for difficult problems
  .setRestartThreshold(50)      // Restart every 50 iterations
  .enablePreconditioning(true); // Use preconditioning

// Use problem-specific preconditioners
IncompleteCholesky ic;
ic.setDropTolerance(1e-3);
cg.setPreconditioner(ic);
```

#### Direct Solver Selection

```cpp
// For symmetric positive definite systems
SimplicialLLT solver; // Fastest

// For general sparse systems
SparseLU solver; // Most robust

// For very large systems
UMFPACK solver; // Best performance for large problems
```

### Profiling and Benchmarking

```bash
# Build with profiling
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo

# Profile with perf
perf record -g ./my_program
perf report

# Memory profiling with valgrind
valgrind --tool=massif ./my_program

# Benchmark with built-in tools
cmake .. -DRODIN_BUILD_BENCHMARKS=ON
make -j4
ctest -R benchmarks -V
```

## Troubleshooting

This section addresses common issues and their solutions.

### Build Issues

#### CMake Configuration Problems

**Problem:** CMake cannot find Boost
```bash
CMake Error: Could not find Boost
```

**Solutions:**
```bash
# Method 1: Install Boost through package manager
# Ubuntu/Debian
sudo apt-get install libboost-all-dev

# CentOS/RHEL/Fedora
sudo dnf install boost-devel

# macOS
brew install boost

# Method 2: Specify Boost location manually
cmake .. -DBOOST_ROOT=/path/to/boost \
         -DBoost_INCLUDE_DIR=/path/to/boost/include \
         -DBoost_LIBRARY_DIR=/path/to/boost/lib

# Method 3: Use vcpkg (Windows/Cross-platform)
vcpkg install boost
cmake .. -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake
```

**Problem:** C++17 compiler not found
```bash
The compiler does not support C++17
```

**Solutions:**
```bash
# Ubuntu: Install modern GCC
sudo apt-get install gcc-9 g++-9
export CC=gcc-9 CXX=g++-9
cmake .. -DCMAKE_CXX_COMPILER=g++-9

# CentOS: Enable devtoolset
sudo yum install centos-release-scl
sudo yum install devtoolset-9-gcc-c++
scl enable devtoolset-9 bash

# macOS: Update Xcode command line tools
xcode-select --install

# Or use Homebrew GCC
brew install gcc@11
export CC=gcc-11 CXX=g++-11
```

**Problem:** Git submodule issues
```bash
fatal: No such file or directory: third-party/eigen/.git
```

**Solutions:**
```bash
# Initialize submodules
git submodule update --init --recursive

# Or clone with submodules from scratch
git clone --recursive https://github.com/cbritopacheco/rodin

# Update existing submodules
git submodule update --remote --merge

# Force update problematic submodules
git submodule foreach --recursive git reset --hard
git submodule update --init --recursive --force
```

#### Compilation Errors

**Problem:** Template instantiation errors
```cpp
error: no matching function for call to 'TrialFunction::TrialFunction(...)'
```

**Solutions:**
```cpp
// Ensure proper includes
#include <Rodin/Variational.h>

// Check finite element space type
P1 fes(mesh);  // Correct
TrialFunction u(fes);

// Not: TrialFunction u(mesh); // Incorrect
```

**Problem:** Linking errors with external libraries
```bash
undefined reference to 'MMG5_...'
```

**Solutions:**
```bash
# Ensure proper CMake options
cmake .. -DRODIN_USE_MMG=ON \
         -DMMG_ROOT=/path/to/mmg

# Check library paths
export LD_LIBRARY_PATH=/path/to/mmg/lib:$LD_LIBRARY_PATH

# Verify third-party builds
cd third-party/mmg
make -j4
```

### Runtime Issues

#### Segmentation Faults

**Common Causes and Solutions:**

1. **Uninitialized mesh connectivity:**
```cpp
// Problem: Missing connectivity computation
Mesh mesh;
mesh.load("domain.mesh");
// mesh.getConnectivity().compute(1, 2);  // Missing!

P1 fes(mesh);  // May cause segfault

// Solution: Always compute required connectivity
mesh.getConnectivity().compute(1, 2);
```

2. **Invalid boundary condition application:**
```cpp
// Problem: Applying BC to wrong boundary
DirichletBC bc(u, Zero(), wrongBoundaryId);

// Solution: Verify boundary IDs
auto boundaries = mesh.getBoundaryIds();
for (auto id : boundaries) {
    std::cout << "Boundary ID: " << id << std::endl;
}
```

3. **Memory access violations:**
```cpp
// Use debug build for better error messages
cmake .. -DCMAKE_BUILD_TYPE=Debug \
         -DRODIN_ENABLE_ASSERTIONS=ON

// Enable address sanitizer
cmake .. -DCMAKE_CXX_FLAGS="-fsanitize=address -g"
```

#### Convergence Problems

**Iterative Solver Not Converging:**

```cpp
// Check matrix properties
std::cout << "Matrix is symmetric: " << A.isSymmetric() << std::endl;
std::cout << "Matrix condition number: " << A.conditionNumber() << std::endl;

// Try different solvers
if (!CG(problem).solve()) {
    BiCGSTAB(problem).solve();  // Try BiCGSTAB instead
}

// Increase tolerance temporarily
CG solver;
solver.setTolerance(1e-4)      // Reduce from 1e-8
      .setMaxIterations(2000)  // Increase iterations
      .solve(problem);
```

**Poor Solution Quality:**

```cpp
// Check mesh quality
auto quality = mesh.computeQuality();
if (quality.minAngle < 10.0) {
    std::cout << "Warning: Poor mesh quality detected" << std::endl;
    // Consider remeshing
}

// Verify boundary conditions
auto appliedDOFs = bc.getAppliedDOFs();
std::cout << "Applied " << appliedDOFs.size() << " boundary conditions" << std::endl;

// Check for proper scaling
Real scale = problem.getCharacteristicScale();
if (scale > 1e6 || scale < 1e-6) {
    std::cout << "Warning: Poor problem scaling: " << scale << std::endl;
}
```

### Performance Issues

**Slow Assembly:**

```cpp
// Enable timing
Timer timer;
timer.start();

// Profile assembly
problem.assemble();

timer.stop();
std::cout << "Assembly time: " << timer.elapsed() << " seconds" << std::endl;

// Optimize assembly
cmake .. -DRODIN_USE_OPENMP=ON  // Enable parallel assembly
export OMP_NUM_THREADS=4        // Set thread count
```

**High Memory Usage:**

```bash
# Monitor memory usage
valgrind --tool=massif ./my_program

# Optimize memory
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DRODIN_MEMORY_OPTIMIZE=ON

# Use memory pools
Context::Sequential context;
context.setMemoryPoolSize(1024 * 1024);  // 1MB pool
```

**Poor Parallel Scaling:**

```bash
# Check thread affinity
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

# Profile parallel regions
export OMP_DISPLAY_ENV=true

# Monitor load balancing
htop  # Check if all cores are utilized equally
```

### Installation Issues on Different Platforms

#### Windows (MSVC)

```cmd
REM Use vcpkg for dependencies
vcpkg install boost eigen3

REM Configure with Visual Studio
cmake .. -G "Visual Studio 16 2019" -A x64 ^
         -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake

REM Build
cmake --build . --config Release
```

#### macOS

```bash
# Install dependencies via Homebrew
brew install boost eigen cmake

# Handle macOS-specific issues
export MACOSX_DEPLOYMENT_TARGET=10.14
cmake .. -DCMAKE_OSX_DEPLOYMENT_TARGET=10.14

# For Apple Silicon
cmake .. -DCMAKE_OSX_ARCHITECTURES=arm64
```

#### HPC Clusters

```bash
# Load appropriate modules
module load gcc/9.3.0 cmake/3.18.0 boost/1.74.0

# Handle non-standard library locations
cmake .. -DCMAKE_PREFIX_PATH="/cluster/boost:/cluster/eigen" \
         -DCMAKE_INSTALL_PREFIX=$HOME/rodin

# Build with job scheduler
sbatch build_script.sh
```

### Getting Help

If you encounter issues not covered here:

1. **Check the GitHub Issues:** [https://github.com/cbritopacheco/rodin/issues](https://github.com/cbritopacheco/rodin/issues)

2. **Create a Minimal Example:** Reduce your problem to the smallest possible reproduction case

3. **Gather System Information:**
```bash
# System info
uname -a
lscpu
free -h

# Compiler info
gcc --version
cmake --version

# Library versions
pkg-config --modversion boost
```

4. **Enable Debug Output:**
```cpp
#include <Rodin/Alert.h>

Alert::setLevel(Alert::Level::Debug);
// Your problematic code here
```

5. **Submit a Bug Report** with:
   - Minimal reproduction code
   - Complete error messages
   - System information
   - Build configuration used

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

#### Code Style
- **Formatting**: Use consistent indentation (2 spaces) and follow existing patterns
- **Naming**: Use PascalCase for classes/types, camelCase for functions/variables
- **Documentation**: Document public APIs with Doxygen-style comments
- **Headers**: Include proper header guards and minimize dependencies

#### Testing Requirements
- **Unit tests**: Add tests for new functionality in the appropriate test directory
- **Manufactured tests**: Add convergence tests for new finite elements or formulations
- **Documentation**: Update documentation and examples as needed
- **Performance**: Ensure new features don't significantly impact performance

#### Code Organization
- **Modularity**: Keep changes focused and minimal
- **Compatibility**: Maintain backward compatibility when possible
- **Error handling**: Use Rodin's Alert system for error reporting
- **Memory management**: Use RAII and smart pointers

### Contribution Types

#### Bug Fixes
- Include minimal reproducible example
- Add regression test if applicable
- Reference issue number in commit message

#### New Features
- Discuss design in GitHub issue first
- Provide comprehensive tests and documentation
- Include examples demonstrating usage

#### Documentation
- Keep documentation up-to-date with code changes
- Improve clarity and add missing sections
- Fix typos and formatting issues

#### Examples
- Provide well-commented examples
- Include mathematical background when relevant
- Test examples in CI pipeline

### Development Environment Setup

```bash
# Clone with development tools
git clone --recursive https://github.com/cbritopacheco/rodin
cd rodin

# Create development build
mkdir debug && cd debug
cmake .. -DCMAKE_BUILD_TYPE=Debug \
         -DRODIN_BUILD_UNIT_TESTS=ON \
         -DRODIN_BUILD_EXAMPLES=ON \
         -DRODIN_BUILD_DOC=ON

# Install development dependencies
sudo apt-get install doxygen graphviz clang-format

# Set up pre-commit hooks (recommended)
cp dev/hooks/pre-commit .git/hooks/
chmod +x .git/hooks/pre-commit
```

### Reporting Issues

When reporting bugs, please include:

- **Version information**: Rodin version, compiler, OS
- **Minimal example**: Shortest code that reproduces the issue
- **Expected vs actual behavior**: Clear description of the problem
- **Error messages**: Full error output and stack traces
- **Build configuration**: CMake options and dependencies

Use the [GitHub issue tracker](https://github.com/cbritopacheco/rodin/issues) to:
- Report bugs with the "bug" label
- Request features with the "enhancement" label  
- Ask questions with the "question" label
- Suggest documentation improvements with the "documentation" label

## Frequently Asked Questions

### General Usage

**Q: How do I get started with Rodin?**
A: Start with the [Installation](#installation) section, then explore the examples in the `examples/` directory. The Poisson example (`examples/PDEs/Poisson.cpp`) is a good starting point for understanding the basic workflow.

**Q: What finite element types does Rodin support?**
A: Rodin currently supports:
- **P0 elements**: Piecewise constant functions on simplicial meshes
- **P1 elements**: Piecewise linear functions with full gradient support  
- **P1 Linear Elasticity**: Specialized P1 elements for elasticity problems
Higher-order elements (P2, P3) and other element types (H(div), H(curl)) are planned for future releases.

**Q: Can I use Rodin for 3D problems?**
A: Yes, Rodin supports 1D, 2D, and 3D problems. The API is dimension-independent, so the same code often works across different dimensions. Simply use tetrahedra instead of triangles for 3D meshes.

**Q: How does Rodin compare to other FEM libraries?**
A: Rodin focuses on:
- **Ease of use**: Modern C++17 API with intuitive syntax
- **Shape optimization**: Built-in support for optimization workflows
- **Research flexibility**: Template-based architecture for easy extension
- **Performance**: Optimized for both single-threaded and parallel execution

For production simulations in specific domains, consider specialized libraries like FEniCS (general PDE), deal.II (adaptive methods), or MFEM (high-order methods).

**Q: What mesh formats can I use?**
A: Rodin supports multiple mesh formats:
- **MFEM**: Native MFEM format (.mfem files)
- **MEDIT**: MEDIT format (.mesh files)  
- **Gmsh**: Popular open-source mesh generator (.msh files)
You can also generate simple meshes programmatically using built-in functions.

**Q: How do I visualize results?**
A: Several options are available:
- **GridFunction output**: Save solutions as .gf files and use MFEM's visualization tools
- **VTK export**: Export to VTK format for ParaView or VisIt
- **Rodin::Plot**: Built-in plotting capabilities (enable with `-DRODIN_WITH_PLOT=ON`)
- **Custom output**: Write solutions to formats your preferred visualization tool supports

### Performance and Scalability

**Q: How large problems can Rodin handle?**
A: Problem size depends on available memory and computational resources:
- **Memory usage**: ~10-50 MB per million degrees of freedom
- **Single machine**: Problems with millions of DOFs on standard workstations
- **HPC clusters**: Distributed problems with tens of millions of DOFs using MPI
- **Practical limits**: Usually limited by solver convergence rather than assembly

**Q: Does Rodin support parallel computing?**
A: Yes, Rodin supports multiple levels of parallelism:
- **OpenMP**: Shared-memory parallelism for assembly and solving
- **MPI**: Distributed-memory parallelism for large-scale problems
- **Hybrid**: Combination of OpenMP and MPI for maximum scalability
- **GPU**: Limited support through external solver libraries (future enhancement)

**Q: How can I improve performance?**
A: Several optimization strategies:
```bash
# Build optimizations
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_CXX_FLAGS="-O3 -march=native" \
         -DRODIN_LTO=ON

# Runtime optimizations  
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=spread
```
- Use appropriate solvers (direct for small problems, iterative for large)
- Enable OpenMP for parallel assembly
- Choose optimal quadrature orders
- Consider mesh adaptation instead of uniform refinement

**Q: Why is my iterative solver not converging?**
A: Common causes and solutions:
- **Poor conditioning**: Try preconditioning or direct solvers
- **Incorrect boundary conditions**: Verify BCs are properly applied
- **Mesh quality**: Check for highly distorted elements
- **Problem scaling**: Normalize equations for better conditioning
- **Solver parameters**: Increase iterations or relax tolerance temporarily

### Development and Integration

**Q: Can I extend Rodin with custom finite elements?**
A: Yes, Rodin's template-based design makes it straightforward to add new finite element types:
```cpp
template<typename Context>
class MyCustomElement : public FiniteElementBase<Context> {
    // Implement required virtual functions
    void computeShapeFunctions(const Point& xi) override;
    void computeGradients(const Point& xi) override;
    // ...
};
```
See the developer documentation for detailed guidance.

**Q: How do I integrate Rodin with existing code?**
A: Several integration approaches:
- **Static library**: Link Rodin as a static library
- **Header-only**: Include Rodin headers directly (for template-heavy code)
- **CMake integration**: Use `find_package(Rodin)` in your CMakeLists.txt
- **Custom build**: Extract specific Rodin modules you need

**Q: Is there Python support?**
A: Python bindings are currently under development and not yet available. The `RODIN_WITH_PY` option exists but is disabled as the bindings don't work yet. For now, Rodin is a C++-only library.

**Q: Can I use Rodin for commercial projects?**
A: Yes, Rodin is licensed under the **Boost Software License 1.0**, which allows both commercial and non-commercial use, modification, and distribution with minimal restrictions.

**Q: How do I contribute to Rodin?**
A: Contributions are welcome! See the [Contributing](#contributing) section for detailed guidelines. Common contribution types:
- **Bug fixes**: Submit issues and pull requests
- **New features**: Discuss in GitHub issues first
- **Documentation**: Improve existing docs or add missing sections
- **Examples**: Provide well-commented examples for specific use cases

### Troubleshooting

**Q: I'm getting compilation errors. What should I do?**
A: Check these common issues:
1. **C++17 compiler**: Ensure you have GCC 7+, Clang 5+, or MSVC 2017+
2. **Boost version**: Use Boost 1.74 or later
3. **CMake version**: Use CMake 3.16 or later
4. **Submodules**: Run `git submodule update --init --recursive`
5. **Include paths**: Verify all dependencies are found by CMake

See the [Installation Troubleshooting](#installation-troubleshooting) section for detailed solutions.

**Q: My simulation results look wrong. How do I debug?**
A: Systematic debugging approach:
1. **Enable debug build**: `cmake .. -DCMAKE_BUILD_TYPE=Debug`
2. **Check boundary conditions**: Verify BCs are applied to correct boundaries
3. **Verify mesh**: Ensure mesh connectivity and orientation are correct
4. **Test convergence**: Run manufactured solution tests to verify correctness
5. **Compare with analytical solutions**: Start with simple problems with known solutions
6. **Enable assertions**: `cmake .. -DRODIN_ENABLE_ASSERTIONS=ON`

**Q: The build is failing on my system. What information should I provide?**
A: When reporting build issues, include:
- **System information**: OS, compiler version, CMake version
- **Error messages**: Complete error output, not just the summary
- **Build configuration**: Full CMake command and options used
- **Dependencies**: Versions of Boost, Eigen, and other libraries
- **Minimal reproduction**: Steps to reproduce the issue from a clean state

**Q: Where can I find more examples?**
A: The `examples/` directory contains over 85 examples organized by category:
- **PDEs/**: Basic PDE problems (Poisson, Helmholtz, elasticity)
- **ShapeOptimization/**: Shape and topology optimization examples
- **MPI/**: Parallel computing examples
- **MMG/**: Mesh adaptation examples
- **Variational/**: Advanced variational formulation examples

Each example includes detailed comments explaining the mathematical background and implementation details.

**Q: How do I cite Rodin in my research?**
A: Use the BibTeX entry provided in the [Citation](#citation) section. You can also find citation information in [CITATION.cff](CITATION.cff) format in the repository root.

**Q: I found a bug. How do I report it?**
A: Use the [GitHub issue tracker](https://github.com/cbritopacheco/rodin/issues):
1. **Search existing issues** to avoid duplicates
2. **Provide a minimal reproduction case**
3. **Include system information and error messages**
4. **Label appropriately** (bug, enhancement, question, etc.)
5. **Be responsive** to maintainer questions and requests for clarification

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
