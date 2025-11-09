# Release Checklist for Rodin v0.0.1

This checklist outlines the steps needed to prepare and publish the first release of Rodin.

## Pre-Release Preparation

### Code Quality & Testing

- [ ] Run full test suite and ensure all tests pass
  ```bash
  cmake .. -DRODIN_BUILD_UNIT_TESTS=ON
  make -j4
  ctest --output-on-failure
  ```

- [ ] Run benchmarks to establish baseline performance
  ```bash
  cmake .. -DRODIN_BUILD_BENCHMARKS=ON
  make -j4
  # Run benchmarks
  ```

- [ ] Build documentation and verify it renders correctly
  ```bash
  cmake .. -DRODIN_BUILD_DOC=ON
  make doc
  ```

- [ ] Verify code coverage is reasonable (optional)
  ```bash
  cmake .. -DRODIN_CODE_COVERAGE=ON
  make -j4
  make RodinCoverage
  ```

### Installation Testing

- [ ] Test installation on clean system (Ubuntu/Debian)
  ```bash
  mkdir build && cd build
  cmake .. -DCMAKE_INSTALL_PREFIX=/tmp/rodin-test
  make -j4
  make install
  ```

- [ ] Verify all headers are installed correctly
  ```bash
  find /tmp/rodin-test/include -name "*.h" -o -name "*.hpp" | wc -l
  ```

- [ ] Test linking against installed library
  ```bash
  # Create test project and verify it builds
  mkdir /tmp/test-rodin && cd /tmp/test-rodin
  # Create CMakeLists.txt and main.cpp
  cmake . -DCMAKE_PREFIX_PATH=/tmp/rodin-test
  make
  ./test_app
  ```

- [ ] Verify CMake config files work correctly
  ```bash
  cmake --find-package -DNAME=Rodin -DCOMPILER_ID=GNU -DLANGUAGE=CXX -DMODE=EXIST
  ```

### Documentation

- [x] Create INSTALL.md with comprehensive installation instructions
- [ ] Update README.md with installation instructions link
- [ ] Create CHANGELOG.md documenting all features in v0.0.1
- [ ] Ensure LICENSE file is present and correct
- [ ] Update CITATION.cff with correct version and date
- [ ] Verify all examples have documentation/comments

### Version Management

- [ ] Update version number in CMakeLists.txt (currently 0.0.1)
- [ ] Update version in CITATION.cff
- [ ] Tag release in git
  ```bash
  git tag -a v0.0.1 -m "First release of Rodin"
  git push origin v0.0.1
  ```

## GitHub Release Process

### Preparing the Release

- [ ] Create release notes summarizing:
  - Main features
  - Dependencies
  - Installation instructions
  - Known limitations
  - Breaking changes (N/A for first release)

- [ ] Prepare binary artifacts (optional):
  - [ ] Ubuntu .deb package
  - [ ] macOS .pkg package
  - [ ] Windows installer
  - [ ] Source tarball

### Creating the Release

1. [ ] Go to GitHub repository → Releases → "Draft a new release"
2. [ ] Choose tag: `v0.0.1` (create new tag if not exists)
3. [ ] Release title: "Rodin v0.0.1 - First Release"
4. [ ] Add release description (use template below)
5. [ ] Attach binary artifacts if available
6. [ ] Mark as pre-release if unstable
7. [ ] Publish release

### Release Description Template

```markdown
# Rodin v0.0.1 - First Release

We are excited to announce the first release of **Rodin**, a lightweight and modular finite element framework for shape and topology optimization.

## ✨ Features

- **Embedded Form Language**: Native C++20 form language for assembling and solving variational formulations
- **Mesh Management**: Full mesh access with iterators and connectivity computation
- **Finite Elements**: Support for P0, P1 elements
- **Solvers**: Integration with Eigen solvers, optional PETSc, SuiteSparse
- **I/O Support**: MFEM, MEDIT, EnSight6 file formats
- **Quadrature**: Multiple quadrature formulas (Grundmann-Moeller, etc.)
- **Third-Party Integration**: MMG for remeshing, Scotch for partitioning

## 📦 Installation

### Quick Start

```bash
git clone --recursive https://github.com/cbritopacheco/rodin.git
cd rodin
mkdir build && cd build
cmake ..
make -j4
sudo make install
```

For detailed installation instructions, see [INSTALL.md](INSTALL.md).

### Requirements

- CMake 3.16+
- C++20 compatible compiler
- Boost 1.74+
- Eigen3 3.4+

## 📚 Documentation

- **Installation Guide**: [INSTALL.md](INSTALL.md)
- **API Documentation**: https://cbritopacheco.github.io/rodin/docs/
- **Examples**: See the `examples/` directory

## 🔧 Usage Example

```cpp
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;

int main() {
  Geometry::Mesh mesh;
  mesh = mesh.UniformGrid(Geometry::Polytope::Type::Triangle, {16, 16});
  
  Variational::P1 Vh(mesh);
  Variational::TrialFunction u(Vh);
  Variational::TestFunction v(Vh);
  
  // Define and solve your variational problem
  // ...
}
```

## ⚠️ Known Limitations

- RodinExternal::MMG is not available in installed version (must build from source for MMG support)
- Some optional features require building from source with specific CMake flags
- Documentation is work in progress

## 🐛 Bug Reports

Please report any issues on our [GitHub Issues](https://github.com/cbritopacheco/rodin/issues) page.

## 📄 License

Rodin is distributed under the Boost Software License, Version 1.0. See [LICENSE](LICENSE) for details.

## 🙏 Acknowledgments

This library was developed as part of Carlos Brito-Pacheco's PhD research.

## 📊 Project Stats

- Lines of Code: [TBD]
- Test Coverage: [TBD]
- Documentation Coverage: [TBD]
```

## Post-Release Tasks

### Verification

- [ ] Verify release appears on GitHub
- [ ] Test download and installation from release tarball
- [ ] Check that links in release notes work
- [ ] Verify documentation site is updated

### Communication

- [ ] Announce release on project website
- [ ] Post announcement on relevant forums/communities
- [ ] Update project status badges
- [ ] Share on social media (optional)

### Repository Updates

- [ ] Update main branch README with installation instructions
- [ ] Add "Installation" section to README.md
- [ ] Update badges to show latest release version
- [ ] Close milestone for v0.0.1
- [ ] Create milestone for v0.1.0 (next release)

## Optional Packaging

### Distribution Packages

- [ ] Submit to Ubuntu PPA
- [ ] Submit to Homebrew
- [ ] Submit to vcpkg
- [ ] Submit to Conan Center
- [ ] Submit to Spack

### Container Images

- [ ] Create Docker image
- [ ] Publish to Docker Hub
- [ ] Create Singularity image

## Future Considerations

### For Next Release (v0.1.0)

- [ ] Add more finite element types (P2, P0)
- [ ] Improve MPI support
- [ ] Add more examples
- [ ] Improve documentation
- [ ] Add CI/CD for automated releases
- [ ] Consider adding Python bindings to PyPI

### Long-term Goals

- [ ] Establish regular release schedule
- [ ] Set up automatic package building
- [ ] Create user mailing list
- [ ] Develop contributor guidelines
- [ ] Create roadmap document
