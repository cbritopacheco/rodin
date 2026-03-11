#!/bin/bash
# Test script to verify Rodin library installation and linking

set -e  # Exit on error

echo "=========================================="
echo "Rodin Installation Test"
echo "=========================================="

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build"
INSTALL_DIR="${BUILD_DIR}/install-test"
TEST_PROJECT_DIR="${BUILD_DIR}/test-installation-project"

echo "Repository root: ${REPO_ROOT}"
echo "Build directory: ${BUILD_DIR}"
echo "Install directory: ${INSTALL_DIR}"
echo "Test project directory: ${TEST_PROJECT_DIR}"
echo ""

# Step 1: Initialize submodules (if not already initialized)
echo "=========================================="
echo "Step 1: Checking submodules..."
echo "=========================================="

cd "${REPO_ROOT}"

# Check if submodules are already initialized
if [ -z "$(ls -A third-party/benchmark 2>/dev/null)" ]; then
  echo "Initializing submodules..."
  git submodule update --init --recursive
  echo "✓ Submodules initialized"
else
  echo "✓ Submodules already initialized"
fi

echo ""

# Step 2: Build and install Rodin
echo "=========================================="
echo "Step 2: Building and installing Rodin..."
echo "=========================================="

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Detect OS for platform-specific CMake options
if [[ "$OSTYPE" == "darwin"* ]]; then
  CMAKE_EXTRA_FLAGS="-DCMAKE_POLICY_VERSION_MINIMUM=3.5"
else
  CMAKE_EXTRA_FLAGS=""
fi

cmake .. \
  -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DRODIN_BUILD_EXAMPLES=OFF \
  -DRODIN_BUILD_UNIT_TESTS=OFF \
  -DRODIN_BUILD_MANUFACTURED_TESTS=OFF \
  -DRODIN_BUILD_BENCHMARKS=OFF \
  ${CMAKE_EXTRA_FLAGS}

# Detect number of processors (cross-platform)
if command -v nproc &> /dev/null; then
  NPROC=$(nproc)
elif command -v sysctl &> /dev/null; then
  NPROC=$(sysctl -n hw.ncpu)
else
  NPROC=2
fi

make -j${NPROC}
make install

echo "✓ Installation complete"
echo ""

# Step 3: Verify installation structure
echo "=========================================="
echo "Step 3: Verifying installation structure..."
echo "=========================================="

# Check for required directories and files
if [ ! -d "${INSTALL_DIR}/include/Rodin" ]; then
  echo "✗ FAILED: Headers not installed"
  exit 1
fi

if [ ! -d "${INSTALL_DIR}/lib/cmake/Rodin" ]; then
  echo "✗ FAILED: CMake config files not installed"
  exit 1
fi

if [ ! -f "${INSTALL_DIR}/lib/cmake/Rodin/RodinConfig.cmake" ]; then
  echo "✗ FAILED: RodinConfig.cmake not found"
  exit 1
fi

echo "✓ Installation structure verified"
echo ""

# Check if we should skip test project (e.g., for macOS with Xcode 16 issues)
if [ "${SKIP_TEST_PROJECT}" = "1" ]; then
  echo "=========================================="
  echo "Skipping test project (SKIP_TEST_PROJECT=1)"
  echo "=========================================="
  echo ""
  echo "✓ Installation verified (build and install successful)"
  echo "✓ Test project skipped due to platform limitations"
  echo ""
  echo "=========================================="
  echo "✓ Installation test completed!"
  echo "=========================================="
  exit 0
fi

# Step 4: Create test project
echo "=========================================="
echo "Step 4: Creating test project..."
echo "=========================================="

rm -rf "${TEST_PROJECT_DIR}"
mkdir -p "${TEST_PROJECT_DIR}"
cd "${TEST_PROJECT_DIR}"

# Copy test files
cp "${SCRIPT_DIR}/test_project_CMakeLists.txt" CMakeLists.txt
cp "${SCRIPT_DIR}/test_poisson.cpp" main.cpp

echo "✓ Test project created"
echo ""

# Step 5: Configure test project
echo "=========================================="
echo "Step 5: Configuring test project..."
echo "=========================================="

# Platform-specific CMake configuration for test project
if [[ "$OSTYPE" == "darwin"* ]]; then
  # On macOS with Xcode 16, we need to ensure proper SDK configuration
  # Let CMake auto-detect the SDK and set deployment target
  CMAKE_TEST_FLAGS="-DCMAKE_OSX_DEPLOYMENT_TARGET=12.0"
  
  # Also ensure we're not using any problematic inherited flags
  export CXXFLAGS=""
  export CFLAGS=""
  export LDFLAGS=""
else
  CMAKE_TEST_FLAGS=""
fi

cmake . -DCMAKE_PREFIX_PATH="${INSTALL_DIR}" ${CMAKE_TEST_FLAGS}

if [ $? -ne 0 ]; then
  echo "✗ FAILED: CMake configuration failed"
  exit 1
fi

echo "✓ Configuration successful"
echo ""

# Step 6: Build test project
echo "=========================================="
echo "Step 6: Building test project..."
echo "=========================================="

make

if [ $? -ne 0 ]; then
  echo "✗ FAILED: Build failed"
  exit 1
fi

echo "✓ Build successful"
echo ""

# Step 7: Run test executable
echo "=========================================="
echo "Step 7: Running test executable..."
echo "=========================================="

./test_app

if [ $? -ne 0 ]; then
  echo "✗ FAILED: Execution failed"
  exit 1
fi

echo "✓ Execution successful"
echo ""

# Step 8: Cleanup (optional)
echo "=========================================="
echo "Step 8: Cleanup..."
echo "=========================================="

if [ "${CLEANUP:-yes}" = "yes" ]; then
  cd "${REPO_ROOT}"
  rm -rf "${INSTALL_DIR}"
  rm -rf "${TEST_PROJECT_DIR}"
  echo "✓ Cleanup complete"
else
  echo "Skipping cleanup (CLEANUP=${CLEANUP})"
fi

echo ""
echo "=========================================="
echo "✓ All installation tests passed!"
echo "=========================================="
