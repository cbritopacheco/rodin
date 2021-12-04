# Rodin

[![License](https://img.shields.io/github/license/cbritopacheco/rodin)](https://github.com/cbritopacheco/rodin/blob/master/LICENSE)
[![Documentation](https://img.shields.io/badge/documentation-master-blue)](https://cbritopacheco.github.io/rodin/)


| OS              | Compiler      |  Status  |
|:---------------:|:-------------:|:--------:|
| Ubuntu 18.04    | GCC 7         | [![Badge](https://byob.yarr.is/cbritopacheco/rodin/Ubuntu-18_04-gcc-7-build_badge)](https://github.com/cbritopacheco/rodin/actions/workflows/Ubuntu-18_04-gcc-7.yml) |
| Ubuntu 18.04    | Clang 7       | [![Badge](https://byob.yarr.is/cbritopacheco/rodin/Ubuntu-18_04-clang-7-build_badge)](https://github.com/cbritopacheco/rodin/actions/workflows/Ubuntu-18_04-clang-7.yml) |
| macOS Catalina  | AppleClang 12 | [![Badge](https://byob.yarr.is/cbritopacheco/rodin/macOS-10_15-clang-12-build_badge)](https://github.com/cbritopacheco/rodin/actions/workflows/macOS-10_15-clang-12.yml) |

## Requirements

- CMake 3.12.0+

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
