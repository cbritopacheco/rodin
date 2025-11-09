# GitHub Copilot Instructions for Rodin

## Repository Overview

Rodin is a lightweight and modular finite element framework in C++20 that provides functionalities for shape and topology optimization algorithms. It includes:
- Native C++20 form language for assembling and solving variational formulations
- Full high-level mesh access and functionalities
- Direct integration with Eigen solvers
- Support for different finite elements and mesh formats
- Integration with third-party libraries like MMG for remeshing

## Project Structure

```
rodin/
├── src/           # Main source code (Rodin and RodinExternal)
├── examples/      # Example applications demonstrating usage
├── tests/         # Test suite
│   ├── benchmarks/        # Performance benchmarks
│   ├── manufactured/      # Manufactured solution tests
│   └── unit/             # Unit tests
├── doc/           # Documentation
├── cmake/         # CMake utility modules
├── third-party/   # Third-party dependencies (git submodules)
├── plugins/       # Plugin modules
└── py/            # Python bindings
```

## Build System

Rodin uses **CMake 3.16+** as its build system.

### Standard Build Process

```bash
# Clone with submodules
git clone --recursive https://github.com/cbritopacheco/rodin
cd rodin

# Configure
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build . -j4
```

### Important CMake Options

- `RODIN_BUILD_SRC=ON` - Build the Rodin source code
- `RODIN_BUILD_EXAMPLES=ON/OFF` - Build examples
- `RODIN_BUILD_UNIT_TESTS=ON/OFF` - Build unit tests
- `RODIN_BUILD_MANUFACTURED_TESTS=ON/OFF` - Build manufactured solution tests
- `RODIN_BUILD_BENCHMARKS=ON/OFF` - Build benchmarks
- `RODIN_BUILD_DOC=ON/OFF` - Build documentation
- `RODIN_USE_MPI=ON/OFF` - Enable MPI support
- `RODIN_USE_OPENMP=ON/OFF` - Enable OpenMP support
- `RODIN_MULTITHREADED=ON/OFF` - Enable multithreading
- `RODIN_CODE_COVERAGE=ON/OFF` - Enable code coverage

## Testing

### Running Tests

```bash
# From build directory
ctest --test-dir tests

# Verbose output on failure
ctest --test-dir tests -V --output-on-failure

# Rerun failed tests
ctest --test-dir tests --rerun-failed -V
```

### Test Types

1. **Unit Tests**: Test individual components and classes
2. **Manufactured Tests**: Verify numerical correctness using manufactured solutions for PDEs
3. **Benchmarks**: Performance measurements

## Code Standards

### Language and Compiler

- **C++ Standard**: C++20
- **Supported Compilers**: GCC 12+, GCC 14+ (for tests)
- **Code Style**: Follow existing patterns in the codebase

### Key Conventions

1. **Namespaces**: Core functionality in `Rodin::` namespace with submodules like `Rodin::Geometry::` and `Rodin::Variational::`
2. **Headers**: Use `.h` for C++ headers
3. **Includes**: Use angle brackets for Rodin headers: `#include <Rodin/Solver.h>`
4. **Naming**: 
   - PascalCase for classes and types
   - camelCase for methods and functions
   - Variables follow context-appropriate conventions

### Example Code Pattern

```cpp
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, 16, 16);
  mesh.getConnectivity().compute(1, 2);

  P1 Vh(mesh);
  TrialFunction u(Vh);
  TestFunction v(Vh);

  Problem problem(u, v);
  problem = Integral(Grad(u), Grad(v))
          - Integral(v)
          + DirichletBC(u, Zero());
  
  Solver::SparseLU solver;
  problem.solve(solver);

  return 0;
}
```

## Dependencies

### Required

- CMake 3.16.0+
- Boost 1.74+
- Eigen3
- C++20 compatible compiler

### Optional (for full functionality)

- libomp (OpenMP support)
- libsuitesparse (sparse linear algebra)
- libmetis (graph partitioning)
- MPI (parallel computing)
- MMG (mesh optimization)
- SDL2, GLFW, OpenAL, Vulkan (for Rodin::Plot module)

## Common Development Workflows

### Adding a New Feature

1. Create feature in appropriate namespace under `src/Rodin/`
2. Add corresponding examples in `examples/`
3. Add unit tests in `tests/unit/` or manufactured tests in `tests/manufactured/`
4. Build and test locally before committing
5. Ensure code coverage is maintained

### Working with Meshes

- Mesh operations are in `Rodin::Geometry::` namespace
- Compute connectivity before using: `mesh.getConnectivity().compute(d1, d2)`
- Iterate over polytopes: `mesh.getCell()`, `mesh.getFace()`, `mesh.getVertex()`

### Working with Variational Formulations

- Define function spaces (P1, P2, etc.)
- Create trial and test functions
- Build problems using `Integral`, boundary conditions, etc.
- Solve using appropriate solver from `Rodin::Solver::`

## CI/CD

The repository uses GitHub Actions with the following workflows:

- **Build.yml**: Builds on Ubuntu and macOS with various configurations
- **Tests.yml**: Runs unit tests, manufactured tests with code coverage
- **Benchmarks.yml**: Performance benchmarking
- **Documentation.yml**: Builds and deploys documentation
- **StaticAnalysis.yml**: Static code analysis

All PRs should pass CI checks before merging.

## Documentation

- Inline documentation uses Doxygen format
- Build documentation: `cmake .. -DRODIN_BUILD_DOC=ON && make`
- Published at: https://cbritopacheco.github.io/rodin/docs/

## Additional Notes

- The repository uses git submodules for third-party dependencies
- Always clone with `--recursive` or run `git submodule update --init --recursive`
- Parallel builds recommended: use `-j4` or higher
- Test environment variables: `OMP_NUM_THREADS=2` for threading tests
