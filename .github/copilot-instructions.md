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

## Installation

The Rodin library supports installation to allow downstream projects to use it via CMake's `find_package()`. This is essential for testing code changes in a realistic environment.

### Quick Installation for Testing

When making code changes, you can install to a temporary location to test:

```bash
# From the repository root
mkdir build && cd build

# Configure with a temporary install prefix
cmake .. -DCMAKE_INSTALL_PREFIX=/tmp/rodin-install -DCMAKE_BUILD_TYPE=Release

# Build and install
make -j4
make install
```

### Installation Types

#### 1. User-Local Installation (Recommended for Testing)

Install to your home directory without requiring sudo:

```bash
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local -DCMAKE_BUILD_TYPE=Release
make -j4
make install
```

#### 2. System-Wide Installation

Install globally (requires sudo):

```bash
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=Release
make -j4
sudo make install
```

#### 3. Staged Installation (For Testing)

Install to a custom directory for testing without affecting the system:

```bash
cmake .. -DCMAKE_INSTALL_PREFIX=/tmp/rodin-test -DCMAKE_BUILD_TYPE=Release
make -j4
make install
```

### What Gets Installed

After installation, the following are available:

- **Headers**: `${CMAKE_INSTALL_PREFIX}/include/Rodin/`
- **Libraries**: `${CMAKE_INSTALL_PREFIX}/lib/`
- **CMake Config**: `${CMAKE_INSTALL_PREFIX}/lib/cmake/Rodin/`
- **Resources**: `${CMAKE_INSTALL_PREFIX}/share/Rodin/resources/`

### Verifying Installation

After installing, verify it works by creating a test project:

```bash
# Create test directory
mkdir /tmp/test-rodin && cd /tmp/test-rodin

# Create CMakeLists.txt
cat > CMakeLists.txt << 'EOF'
cmake_minimum_required(VERSION 3.16)
project(TestRodin CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Rodin REQUIRED)

add_executable(test_app main.cpp)
target_link_libraries(test_app PRIVATE Rodin::Geometry Rodin::Variational)
EOF

# Create test source
cat > main.cpp << 'EOF'
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <iostream>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main() {
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, {8, 8});
  P1 Vh(mesh);
  std::cout << "Test passed! Mesh has " << mesh.getCellCount() 
            << " cells, FE space has " << Vh.getSize() << " DOFs" << std::endl;
  return 0;
}
EOF

# Build the test project
cmake . -DCMAKE_PREFIX_PATH=/tmp/rodin-install  # Adjust path as needed
make
./test_app
```

### Automated Installation Testing

The repository includes automated installation tests in `tests/installation/`:

```bash
# Run the full installation test suite
bash tests/installation/test_installation.sh

# This script will:
# 1. Build Rodin from source
# 2. Install to a temporary directory
# 3. Create a test project
# 4. Verify the test project can find and link against Rodin
# 5. Run a simple Poisson equation solver
# 6. Clean up
```

### Testing Your Code Changes

When working on code changes, use this workflow:

1. **Make your changes** to the source code
2. **Build and install** to a temporary location:
   ```bash
   cd build
   cmake .. -DCMAKE_INSTALL_PREFIX=/tmp/rodin-dev -DCMAKE_BUILD_TYPE=Debug -DRODIN_BUILD_UNIT_TESTS=ON
   make -j4 && make install
   ```
3. **Run unit tests** to verify your changes:
   ```bash
   ctest --test-dir tests/unit -V --output-on-failure
   ```
4. **Test with a downstream project** (optional but recommended):
   ```bash
   # Create a test project that uses the installed library
   mkdir /tmp/test-changes && cd /tmp/test-changes
   # ... create test project as shown above ...
   cmake . -DCMAKE_PREFIX_PATH=/tmp/rodin-dev
   make && ./test_app
   ```
5. **Run installation test** to ensure the library installs correctly:
   ```bash
   bash tests/installation/test_installation.sh
   ```

### CMake Integration for Downstream Projects

When a downstream project uses `find_package(Rodin REQUIRED)`, these targets become available:

**Core Targets:**
- `Rodin::Rodin` - Main interface
- `Rodin::Geometry` - Mesh and geometry functionality
- `Rodin::Variational` - Variational formulations
- `Rodin::Solver` - Solver interfaces
- `Rodin::Math` - Mathematical utilities
- `Rodin::IO` - Input/output functionality
- `Rodin::QF` - Quadrature formulas
- `Rodin::Assembly` - Assembly routines
- `Rodin::FormLanguage` - Form language support

**Utility Targets:**
- `Rodin::Alert` - Logging and alerts
- `Rodin::Utility` - General utilities
- `Rodin::Context` - Execution context
- `Rodin::Threads` - Threading support

**Optional Targets** (when enabled with corresponding CMake flags):
- `Rodin::Plot` - Plotting (requires `RODIN_WITH_PLOT=ON`)
- `Rodin::PETSc` - PETSc integration (requires `RODIN_USE_PETSC=ON`)
- `Rodin::MPI` - MPI support (requires `RODIN_USE_MPI=ON`)

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
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, {16, 16});
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
