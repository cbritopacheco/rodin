# Rodin

[![License](https://img.shields.io/badge/license-BSL--1.0-green)](https://github.com/cbritopacheco/rodin/blob/master/LICENSE)

| Branch      |  Matrix  | Documentation |
|:-----------:|:--------:|:-------------:|
| master      | [![Build](https://github.com/cbritopacheco/rodin/actions/workflows/Build.yml/badge.svg?branch=master)](https://github.com/cbritopacheco/rodin/actions/workflows/Build.yml) | [![Documentation](https://img.shields.io/badge/doc-master-blue)](https://cbritopacheco.github.io/rodin/) |
| develop     | [![Build](https://github.com/cbritopacheco/rodin/actions/workflows/Build.yml/badge.svg?branch=develop)](https://github.com/cbritopacheco/rodin/actions/workflows/Build.yml) | |

## Requirements

- [CMake 3.12.0+](https://cmake.org/)
- [Boost 1.65+](https://www.boost.org/)
- [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html)

Any of these should be available for quick install from your standard package
manager.

## Building the project

```
git clone --recursive https://github.com/carlos-brito-pacheco/rodin
cd rodin
mkdir build && cd build
cmake ..
make -j4
```

### CMake Options

| Option                 | Description                                       |
|------------------------|---------------------------------------------------|
| RODIN_BUILD_EXAMPLES   | Builds the examples in the `examples/` directory. |
| RODIN_BUILD_DOC        | Builds the documentation using Doxygen            |
| RODIN_USE_MCSS         | Builds the documentation using Doxygen and m.css  |

## Building the documentation

See [this page](doc/README.md) to see how to build the documentation.
