# Installation Guide

This guide provides instructions for building and installing the Rodin finite element library.

## Table of Contents

1. [Requirements](#requirements)
2. [Building from Source](#building-from-source)
3. [Installation](#installation)
4. [Using Rodin in Your Project](#using-rodin-in-your-project)
5. [CMake Configuration Options](#cmake-configuration-options)

## Requirements

### Required Dependencies

- **CMake 3.16.0+** - Build system generator
- **C++20 compatible compiler** - GCC 10+, Clang 10+, or MSVC 2019+
- **Boost 1.74+** - Required components: `filesystem`, `serialization`
- **Eigen3 3.4+** - Linear algebra library

### Optional Dependencies

- **OpenMP** - For parallel computations (recommended)
- **MPI** - For distributed computing
- **PETSc** - For advanced solvers
- **Scotch** - For graph partitioning
- **SuiteSparse** - For direct solvers (UMFPACK, CHOLMOD, SPQR, KLU)
- **VTK** - For visualization

## Building from Source

### 1. Clone the Repository

```bash
git clone --recursive https://github.com/cbritopacheco/rodin.git
cd rodin
```

**Important:** Use `--recursive` to initialize git submodules (eigen, googletest, etc.).

If you've already cloned without `--recursive`, run:

```bash
git submodule update --init --recursive
```

### 2. Create Build Directory

```bash
mkdir build
cd build
```

### 3. Configure with CMake

**Basic configuration:**
```bash
cmake ..
```

**Configuration with options:**
```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/usr/local \
  -DRODIN_BUILD_EXAMPLES=ON \
  -DRODIN_MULTITHREADED=ON
```

### 4. Build

```bash
make -j$(nproc)
```

Or specify the number of parallel jobs:
```bash
make -j4
```

## Installation

### System-wide Installation (requires sudo)

```bash
sudo make install
```

This will install:
- **Libraries** to `${CMAKE_INSTALL_PREFIX}/lib`
- **Headers** to `${CMAKE_INSTALL_PREFIX}/include/Rodin`
- **CMake config files** to `${CMAKE_INSTALL_PREFIX}/lib/cmake/Rodin`
- **Resources** to `${CMAKE_INSTALL_PREFIX}/share/Rodin/resources`

### User-local Installation

```bash
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local
make -j4
make install
```

### Staged Installation (for packaging)

```bash
make install DESTDIR=/tmp/rodin-install
```

## Using Rodin in Your Project

### CMake Integration

Create a `CMakeLists.txt` file:

```cmake
cmake_minimum_required(VERSION 3.16)
project(MyProject CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the installed Rodin package
find_package(Rodin REQUIRED)

# Create your executable
add_executable(my_app main.cpp)

# Link against Rodin libraries
target_link_libraries(my_app PRIVATE
  Rodin::Geometry
  Rodin::Variational
)
```

### Example Code

```cpp
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <iostream>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main()
{
  // Create a uniform triangular mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, {16, 16});
  
  // Define finite element space
  P1 Vh(mesh);
  
  std::cout << "Mesh has " << mesh.getCellCount() << " cells" << std::endl;
  std::cout << "FE space has " << Vh.getSize() << " DOFs" << std::endl;
  
  return 0;
}
```

### Build Your Project

```bash
mkdir build && cd build
cmake ..
make
./my_app
```

If Rodin is installed in a non-standard location:

```bash
cmake .. -DCMAKE_PREFIX_PATH=/path/to/rodin/install
```

## CMake Configuration Options

### Build Options

| Option | Description | Default |
|--------|-------------|---------|
| `RODIN_BUILD_SRC` | Build the Rodin source code | `ON` |
| `RODIN_BUILD_EXAMPLES` | Build the Rodin examples | `ON` |
| `RODIN_BUILD_UNIT_TESTS` | Build Rodin unit tests | `OFF` |
| `RODIN_BUILD_MANUFACTURED_TESTS` | Build manufactured solution tests | `OFF` |
| `RODIN_BUILD_BENCHMARKS` | Build Rodin benchmarks | `OFF` |
| `RODIN_BUILD_DOC` | Build the Rodin documentation | `OFF` |

### Feature Options

| Option | Description | Default |
|--------|-------------|---------|
| `RODIN_USE_OPENMP` | Enable OpenMP support | `OFF` (auto `ON` with `RODIN_MULTITHREADED`) |
| `RODIN_MULTITHREADED` | Enable multithreading capabilities | `ON` |
| `RODIN_THREAD_SAFE` | Enforce thread-safety | `OFF` (auto `ON` with `RODIN_MULTITHREADED`) |
| `RODIN_WITH_PLOT` | Build the Rodin::Plot module | `OFF` |
| `RODIN_WITH_PY` | Build Python bindings | `OFF` |

### Solver Options

| Option | Description | Default |
|--------|-------------|---------|
| `RODIN_USE_MPI` | Enable MPI support | `OFF` |
| `RODIN_USE_PETSC` | Enable PETSc support | `OFF` |
| `RODIN_USE_UMFPACK` | Enable UMFPACK solver | `OFF` |
| `RODIN_USE_CHOLMOD` | Enable CHOLMOD solver | `OFF` |
| `RODIN_USE_SPQR` | Enable SPQR solver | `OFF` |
| `RODIN_USE_KLU` | Enable KLU solver | `OFF` |
| `RODIN_USE_SUPERLU` | Enable SuperLU solver | `OFF` |
| `RODIN_USE_PASTIX` | Enable PaStiX solver | `OFF` |
| `RODIN_USE_PARDISO` | Enable Pardiso solver | `OFF` |
| `RODIN_USE_SCOTCH` | Enable Scotch support | `OFF` |

### Other Options

| Option | Description | Default |
|--------|-------------|---------|
| `BUILD_SHARED_LIBS` | Build using shared libraries | `OFF` |
| `RODIN_LTO` | Compile with link time optimization | `OFF` |
| `RODIN_CODE_COVERAGE` | Compile with code coverage flags | `OFF` |
| `RODIN_SILENCE_WARNINGS` | Silence warnings outputted by Rodin | `OFF` |
| `RODIN_SILENCE_EXCEPTIONS` | Silence exceptions thrown by Rodin | `ON` |

## Available Rodin CMake Targets

When you `find_package(Rodin)`, the following targets become available:

### Core Targets
- `Rodin::Rodin` - Main interface library
- `Rodin::Geometry` - Mesh and geometry functionality
- `Rodin::Variational` - Variational formulations
- `Rodin::Solver` - Solver interfaces
- `Rodin::Math` - Mathematical utilities
- `Rodin::IO` - Input/output functionality
- `Rodin::QF` - Quadrature formulas
- `Rodin::Assembly` - Assembly routines
- `Rodin::FormLanguage` - Form language support

### Utility Targets
- `Rodin::Alert` - Logging and alerts
- `Rodin::Utility` - General utilities
- `Rodin::Context` - Execution context
- `Rodin::Threads` - Threading support
- `Rodin::Serialization` - Serialization support
- `Rodin::Test` - Testing utilities

### Model Targets
- `Rodin::Models::Eikonal` - Eikonal equation models
- `Rodin::Models::Advection` - Advection models

### Optional Targets (when enabled)
- `Rodin::Plot` - Plotting (requires `RODIN_WITH_PLOT=ON`)
- `Rodin::PETSc` - PETSc integration (requires `RODIN_USE_PETSC=ON`)
- `Rodin::Scotch` - Scotch integration (requires `RODIN_USE_SCOTCH=ON`)
- `Rodin::MPI` - MPI support (requires `RODIN_USE_MPI=ON`)

## Troubleshooting

### Missing Dependencies

If CMake cannot find required dependencies:

```bash
# For Eigen3
sudo apt-get install libeigen3-dev  # Ubuntu/Debian
brew install eigen                   # macOS

# For Boost
sudo apt-get install libboost-all-dev  # Ubuntu/Debian
brew install boost                      # macOS
```

### Submodules Not Initialized

If you see errors about missing third-party libraries:

```bash
git submodule update --init --recursive
```

### Custom Dependency Paths

If dependencies are installed in non-standard locations:

```bash
cmake .. \
  -DEigen3_DIR=/path/to/eigen3/share/eigen3/cmake \
  -DBoost_ROOT=/path/to/boost
```

## Getting Help

- **Documentation**: https://cbritopacheco.github.io/rodin/docs/
- **GitHub Issues**: https://github.com/cbritopacheco/rodin/issues
- **Examples**: See the `examples/` directory in the source tree

## License

Rodin is distributed under the Boost Software License, Version 1.0.
See the `LICENSE` file for details.
