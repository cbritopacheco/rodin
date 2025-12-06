# Rodin Installation Test

This directory contains tests to verify that the Rodin library can be properly installed and used by downstream projects.

## Purpose

The installation test ensures that:
1. The library builds and installs correctly
2. All required headers and libraries are installed
3. CMake configuration files are properly generated
4. Downstream projects can find the library using `find_package(Rodin REQUIRED)`
5. Projects can link against Rodin targets
6. The installed library can solve a simple Poisson equation

## Files

- `test_installation.sh` - Main test script that orchestrates the installation test
- `test_project_CMakeLists.txt` - CMakeLists.txt template for the test project
- `test_poisson.cpp` - Simple Poisson equation solver used as the test executable

## Running the Test Locally

From the repository root:

```bash
# Run the installation test
bash tests/installation/test_installation.sh
```

The script will:
1. Build Rodin from source
2. Install to a temporary directory (`build/install-test`)
3. Create a test project that uses the installed library
4. Configure, build, and run the test project
5. Verify the results
6. Clean up temporary files

To skip cleanup (for debugging):

```bash
CLEANUP=no bash tests/installation/test_installation.sh
```

## GitHub Workflow

The installation test is also run automatically via GitHub Actions on every push and pull request. See `.github/workflows/Installation.yml` for details.

## What Gets Tested

### Build and Installation
- Library compiles without errors
- Installation completes successfully
- Headers are installed to the correct location
- Libraries are installed to the correct location
- CMake config files are generated

### CMake Integration
- `find_package(Rodin REQUIRED)` works
- CMake targets are available (`Rodin::Geometry`, `Rodin::Variational`, `Rodin::Solver`)
- No CMake warnings or errors

### Compilation and Linking
- Test project compiles without errors
- All required libraries are linked correctly
- No missing symbols

### Runtime Execution
- Test executable runs successfully
- Poisson equation is solved correctly
- Results are validated

## Expected Output

When the test passes, you should see:

```
==========================================
Rodin Installation Test
==========================================
...
Step 1: Building and installing Rodin...
✓ Installation complete

Step 2: Verifying installation structure...
✓ Installation structure verified

Step 3: Creating test project...
✓ Test project created

Step 4: Configuring test project...
✓ Configuration successful

Step 5: Building test project...
✓ Build successful

Step 6: Running test executable...
========================================
✓ Poisson equation solved successfully!
========================================
Mesh cells: 450
DOFs: 256
========================================
✓ Installation test PASSED
✓ Execution successful

Step 7: Cleanup...
✓ Cleanup complete

==========================================
✓ All installation tests passed!
==========================================
```
