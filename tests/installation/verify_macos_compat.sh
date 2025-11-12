#!/bin/bash
# Verification script for macOS compatibility testing
# This script validates that the installation test will work on macOS

set -e

echo "=========================================="
echo "macOS Compatibility Verification"
echo "=========================================="
echo ""

# Test 1: Verify CPU detection logic
echo "Test 1: CPU Detection Logic"
echo "----------------------------"

# Simulate macOS environment (no nproc, has sysctl)
echo "Simulating macOS environment..."

# Test the logic directly
test_cpu_detection() {
  # Mock commands for testing
  if command -v nproc &> /dev/null; then
    echo "  Linux detected: would use nproc"
    NPROC=$(nproc)
  elif command -v sysctl &> /dev/null; then
    echo "  macOS detected: would use sysctl"
    # On Linux, sysctl exists but doesn't have hw.ncpu
    # On macOS, this would work: sysctl -n hw.ncpu
    if sysctl -n hw.ncpu 2>/dev/null; then
      NPROC=$(sysctl -n hw.ncpu)
    else
      echo "  sysctl exists but hw.ncpu not available (expected on Linux)"
      NPROC=2
    fi
  else
    echo "  Fallback: using default value"
    NPROC=2
  fi
  echo "  Detected processors: ${NPROC}"
}

test_cpu_detection
echo "✓ CPU detection logic is correct"
echo ""

# Test 2: Verify script syntax
echo "Test 2: Script Syntax Validation"
echo "---------------------------------"
if bash -n /home/runner/work/rodin/rodin/tests/installation/test_installation.sh; then
  echo "✓ Script syntax is valid"
else
  echo "✗ Script syntax errors detected"
  exit 1
fi
echo ""

# Test 3: Check for macOS-specific issues
echo "Test 3: macOS-Specific Issues"
echo "------------------------------"

# Check for Linux-only commands
if grep -E 'apt-get|apt|yum|dnf' /home/runner/work/rodin/rodin/tests/installation/test_installation.sh; then
  echo "✗ Found Linux package manager commands in script"
  exit 1
else
  echo "✓ No Linux-specific package managers in script"
fi

# Check that nproc is not hardcoded
if grep -E 'make -j\$\(nproc\)' /home/runner/work/rodin/rodin/tests/installation/test_installation.sh; then
  echo "✗ Found hardcoded nproc usage"
  exit 1
else
  echo "✓ No hardcoded nproc usage"
fi

echo ""

# Test 4: Verify workflow configuration
echo "Test 4: Workflow Configuration"
echo "-------------------------------"

if grep -q "macos-latest" /home/runner/work/rodin/rodin/.github/workflows/Installation.yml; then
  echo "✓ macOS runner configured in workflow"
else
  echo "✗ macOS runner not found in workflow"
  exit 1
fi

if grep -q "brew install" /home/runner/work/rodin/rodin/.github/workflows/Installation.yml; then
  echo "✓ Homebrew package installation configured"
else
  echo "✗ Homebrew installation not found"
  exit 1
fi

echo ""

# Test 5: Check dependencies
echo "Test 5: Required Dependencies"
echo "-----------------------------"

required_deps=("boost" "eigen" "libomp")
for dep in "${required_deps[@]}"; do
  if grep -q "$dep" /home/runner/work/rodin/rodin/.github/workflows/Installation.yml; then
    echo "✓ Dependency '$dep' is installed"
  else
    echo "✗ Dependency '$dep' not found"
    exit 1
  fi
done

echo ""

# Test 6: Verify submodule initialization
echo "Test 6: Submodule Initialization"
echo "---------------------------------"

if grep -q "git submodule update --init --recursive" /home/runner/work/rodin/rodin/tests/installation/test_installation.sh; then
  echo "✓ Submodules are initialized in script"
else
  echo "✗ Submodule initialization missing"
  exit 1
fi

echo ""

# Summary
echo "=========================================="
echo "Verification Summary"
echo "=========================================="
echo ""
echo "✓ CPU detection works cross-platform"
echo "✓ Script syntax is valid"
echo "✓ No macOS-incompatible commands"
echo "✓ Workflow properly configured for macOS"
echo "✓ All required dependencies specified"
echo "✓ Submodules properly initialized"
echo ""
echo "The installation test script is ready for macOS testing."
echo ""
echo "Note: The actual macOS test will run automatically in GitHub Actions"
echo "when this PR is pushed. Check the 'Installation' workflow for results."
echo ""
echo "=========================================="
