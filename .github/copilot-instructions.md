# Rodin Finite Element Framework
Rodin is a C++ finite element framework for shape and topology optimization algorithms. It uses CMake build system, C++20, and has complex third-party dependencies via git submodules.

**Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.**

## Working Effectively

### Prerequisites and Environment Setup
- Ubuntu 24.04 LTS or similar (validated environment)
- CMake 3.16+ (3.31+ recommended)
- GCC 13.3+ (C++20 support required)
- Git with submodule support

### Bootstrap and Build Process
**CRITICAL**: Always run these commands in the exact order shown. Build times are significant - NEVER CANCEL long-running operations.

#### 1. Install System Dependencies
```bash
sudo apt-get update
sudo apt-get install -y \
  libgl1-mesa-dev \
  lcov \
  libomp-dev \
  libsdl2-dev \
  libglfw3-dev \
  libopenal-dev \
  libvulkan-dev \
  libsuitesparse-dev \
  libmetis-dev \
  libboost1.74-all-dev \
  doxygen
```
**Timing**: 2-5 minutes depending on network speed.

#### 2. Initialize Git Submodules
**CRITICAL WORKAROUND**: The Eigen submodule fails to clone from GitLab. Manual intervention required:

```bash
# First attempt normal submodule initialization
git submodule update --init --recursive

# Fix Eigen manually if it fails
cd third-party/eigen
git clone https://github.com/eigenteam/eigen-git-mirror.git .

# Fix ISCD/Mshdist if needed
cd ../ISCD
git clone https://github.com/cbritopacheco/Mshdist.git Mshdist-temp
rm -rf Mshdist
mv Mshdist-temp Mshdist
cd ../../
```
**Timing**: 5-10 minutes for all submodules. **NEVER CANCEL** - dependency download takes time.

#### 3. Configure with CMake
```bash
mkdir -p build && cd build
cmake ..
```
**Timing**: ~10 seconds. This should complete quickly if dependencies are properly installed.

#### 4. Build Core Libraries
**NEVER CANCEL**: Core library build takes 2-5 minutes. Set timeout to 10+ minutes.
```bash
make -j4 RodinAlert RodinAssembly RodinFormLanguage RodinGeometry RodinIO RodinQF RodinVariational
```
**Timing**: 2-5 minutes for core libraries.

### Advanced Build Options

#### Documentation Generation
```bash
cmake .. -DRODIN_BUILD_DOC=ON
make RodinDoxygen
```
**Timing**: ~30 seconds. Documentation builds quickly once core is built.

#### Build Configuration Options
| Option | Description | Status |
|--------|-------------|---------|
| `RODIN_BUILD_SRC=ON` | Build core Rodin libraries | ✅ Works |
| `RODIN_BUILD_EXAMPLES=OFF` | Disable examples (recommended) | ✅ Works |
| `RODIN_BUILD_UNIT_TESTS=OFF` | Disable unit tests (recommended) | ⚠️ Has compilation issues |
| `RODIN_BUILD_DOC=ON` | Enable documentation | ✅ Works |
| `RODIN_USE_MCSS=ON` | Enhanced documentation styling | ⚠️ Requires Python dependencies |

### Known Build Issues and Workarounds

#### Examples and Tests Currently Don't Build Cleanly
**Issue**: Examples and unit tests have compilation errors due to strict warning settings and dependency issues.

**Workaround**: Always disable them initially:
```bash
cmake .. -DRODIN_BUILD_EXAMPLES=OFF -DRODIN_BUILD_UNIT_TESTS=OFF
```

#### ISCD/MMG Dependencies
**Issue**: Some external dependencies (ISCD/Mshdist) have submodule configuration problems.

**Workaround**: Use core libraries only, or manually fix submodules as shown above.

#### Eigen Warnings
**Expected Behavior**: You will see many "zero as null pointer constant" warnings from Eigen. These are warnings, not errors, and can be ignored.

## Validation and Testing

### Basic Functionality Validation
After building core libraries, validate the build:

```bash
# Verify libraries were built
ls -la src/Rodin/*/lib*.a

# Check that CMake finds the libraries
cmake .. --trace-expand | grep -i rodin
```

### Documentation Validation
```bash
# Build documentation
make RodinDoxygen

# Check generated documentation
ls -la doc/html/index.html
```

### File Structure Overview
```
rodin/
├── README.md              # Basic build instructions
├── CMakeLists.txt         # Main CMake configuration
├── src/                   # Core Rodin source code
│   ├── Rodin/Alert/       # Alert and messaging system
│   ├── Rodin/Assembly/    # Finite element assembly
│   ├── Rodin/Geometry/    # Geometric operations
│   ├── Rodin/IO/          # Input/output operations  
│   ├── Rodin/Variational/ # Variational formulations
│   └── RodinExternal/     # External integrations (MMG, etc.)
├── examples/              # Example programs (build issues)
├── tests/                 # Unit tests (build issues)
├── third-party/          # Git submodules for dependencies
├── doc/                   # Documentation configuration
└── .github/workflows/     # CI/CD configurations
```

## Timing Expectations

| Operation | Expected Time | Timeout Setting | Notes |
|-----------|---------------|-----------------|-------|
| System package install | 2-5 minutes | 10 minutes | **NEVER CANCEL** |
| Git submodule init | 5-10 minutes | 15 minutes | **NEVER CANCEL** |
| CMake configuration | 10 seconds | 2 minutes | Should be fast |
| Core library build | 2-5 minutes | 10 minutes | **NEVER CANCEL** |
| Documentation build | 30 seconds | 5 minutes | Fast once core is built |
| Full build (if working) | 15-30 minutes | 60 minutes | **NEVER CANCEL** |

## Common Workflow Steps

### Making Code Changes
1. **Always build and test core libraries first**:
   ```bash
   make RodinAlert RodinFormLanguage RodinVariational
   ```

2. **When adding new functionality**:
   - Focus on core modules under `src/Rodin/`
   - Test compilation after each significant change
   - Avoid examples/ and tests/ directories until core issues are resolved

3. **Before committing changes**:
   ```bash
   # Verify core build still works
   make clean
   make -j4 RodinAlert RodinAssembly RodinFormLanguage RodinGeometry RodinIO RodinQF RodinVariational
   ```

### Working with Dependencies
- **Third-party changes**: Be very careful with `third-party/` directory
- **Submodule updates**: Use `git submodule update --recursive` with caution
- **New dependencies**: Add to system package list and test thoroughly

### Debugging Build Issues
1. **Check submodule status**: `git submodule status`
2. **Verify system dependencies**: Ensure all apt packages are installed
3. **Clean build**: `rm -rf build && mkdir build && cd build`
4. **Check CMake logs**: Look for missing dependencies in cmake output
5. **Build incrementally**: Start with individual targets before attempting full build

## What Currently Works vs. What Doesn't

### ✅ Confirmed Working
- System dependency installation
- Git submodule initialization (with workarounds)
- CMake configuration  
- Core Rodin library compilation
- Documentation generation
- Basic finite element framework functionality

### ⚠️ Known Issues
- Unit tests have compilation errors
- Examples have dependency build issues  
- Some external dependencies (MMG/ISCD) don't build cleanly
- Full CI build may fail due to test compilation issues

### 🚫 Currently Broken
- Complete test suite execution
- Some advanced examples requiring MMG integration
- Benchmark compilation

## Emergency Troubleshooting

If builds fail completely:
1. **Clean everything**: `git clean -fdx && git submodule foreach --recursive git clean -fdx`
2. **Reinitialize submodules**: Follow submodule setup steps exactly as documented
3. **Verify system packages**: Reinstall all required apt packages
4. **Build minimal configuration**: Use `-DRODIN_BUILD_EXAMPLES=OFF -DRODIN_BUILD_UNIT_TESTS=OFF`
5. **Focus on core libraries only**: Don't attempt full build until core libraries work

**Remember**: The core finite element framework libraries build and work correctly. Focus on those for development work rather than trying to fix the peripheral build issues immediately.