# H1 Grad Tests - Complete Report Index

## Overview

This folder contains comprehensive test execution reports for the `RodinVariationalH1GradTest` test suite. The tests were built and executed successfully, revealing systematic gradient computation errors in the H1 Grad operator.

## Quick Summary

- **Build Status**: ✅ SUCCESS
- **Tests Passed**: 5 out of 10 (50%)
- **Tests Failed**: 5 out of 10 (50%)
- **Root Cause**: Jacobian transformation bug in gradient evaluation
- **Priority**: HIGH

## Report Files

### 1. **TEST_EXECUTION_REPORT.md** (Start Here!)
📄 **Location**: `TEST_EXECUTION_REPORT.md`
**Size**: 7.3 KB

**Best for**: Understanding what failed and why
**Contains**:
- Executive summary
- Test results overview (5 passing, 5 failing)
- Critical findings about the bug
- Detailed error messages for each test
- Root cause analysis
- Recommendations for fixing

**Read this first** if you want a comprehensive overview.

---

### 2. **H1_GRAD_TEST_SUMMARY.md**
📄 **Location**: `H1_GRAD_TEST_SUMMARY.md`
**Size**: 6.7 KB

**Best for**: Detailed failure analysis by test
**Contains**:
- Build status and configuration
- Complete test execution summary
- Detailed failure breakdown for each of 5 failing tests
- Pattern analysis across tests
- Suspected problem areas in code
- Next steps for investigation

**Read this** if you want test-by-test details.

---

### 3. **H1_GRAD_DETAILED_FAILURES.txt**
📄 **Location**: `H1_GRAD_DETAILED_FAILURES.txt`
**Size**: 17 KB

**Best for**: Detailed numerical analysis
**Contains**:
- Formatted failure tables with expected vs. actual values
- Complete error messages with line numbers
- Detailed analysis of each failure
- Cross-test pattern analysis
- ASCII-formatted tables comparing values
- Comprehensive failure summary

**Read this** if you need exact error values and detailed breakdown.

---

### 4. **QUICK_REFERENCE.txt**
📄 **Location**: `QUICK_REFERENCE.txt`
**Size**: 2.7 KB

**Best for**: Quick lookup and debugging
**Contains**:
- Build and run commands
- Quick test statistics
- Failed tests summary
- Bug location
- Root cause indicators
- Quick debug commands

**Use this** as a quick reference while debugging.

---

## The Bug at a Glance

### Location
**File**: `src/Rodin/Variational/H1/Grad.h`  
**Method**: `interpolate()` for GridFunction specialization  
**Lines**: 160-172, specifically line 171

### The Problem
```cpp
out = p.getJacobianInverse().transpose() * s_res;
```

The Jacobian inverse transpose multiplication appears to be incorrect. It might be:
- Applied to the wrong vector
- Using transpose when it shouldn't (or vice versa)
- Computing the reference gradient incorrectly

### Evidence
- Y-component of gradient is **consistently wrong** across ALL tests
- Same bug affects H1<1>, H1<2>, and H1<3> basis functions
- Bug only affects gradient **evaluation**, not construction or assembly
- Linear functions (which should be exact) are computed incorrectly

---

## Test Results Summary

### ✅ Passing Tests (5)
1. **ShapeFunction_Construction** - Can construct Grad for shape functions
2. **GridFunction_Construction** - Can construct Grad for grid functions
3. **GridFunction_ConstantFunction** - Correctly returns zero gradient for constants
4. **Copy** - Can copy Grad objects correctly
5. **UsageInBilinearForm** - Can use Grad in bilinear forms and assemble

### ❌ Failing Tests (5)
1. **GridFunction_LinearFunction** - f(x,y)=3x+4y → Y-gradient is -2.0 (expected 4.0)
2. **GridFunction_QuadraticFunction_H1_2** - f(x,y)=x²+y² → Magnitudes completely wrong
3. **GridFunction_QuadraticFunction_H1_3** - f(x,y)=x²+2xy+y² → Y-component consistently negative
4. **MultipleEvaluations** - f(x,y)=x+2y → Y-gradient is zero at all points
5. **DifferentPolynomialDegrees** - f(x,y)=3x+4y → Y-gradient consistently wrong across all spaces

---

## How to Use These Reports

### If you want a quick understanding:
1. Read **QUICK_REFERENCE.txt** (2 min)
2. Skim **TEST_EXECUTION_REPORT.md** (5 min)

### If you want to fix the bug:
1. Read **TEST_EXECUTION_REPORT.md** - understand the problem
2. Check **H1_GRAD_DETAILED_FAILURES.txt** - get exact error values
3. Review **Grad.h lines 160-172** - see the buggy code
4. Compare with **P1/Grad.h** - see a working implementation

### If you need exact error values:
1. Go to **H1_GRAD_DETAILED_FAILURES.txt**
2. Look at the formatted tables for each test

---

## Key Error Examples

### Example 1: Linear Function
```
Function:    f(x,y) = 3x + 4y
Expected ∇f: [3.0, 4.0]
Got ∇f:      [?, -2.0]
Error:       Y-component off by factor of -6
```

### Example 2: Quadratic Function
```
Function:    f(x,y) = x² + y²
Expected ∇f: [3.0, 4.0] at (1.5, 2.0)
Got ∇f:      [13.5, 5.5]
Error:       X-component 450% too large, Y-component 37.5% too large
```

### Example 3: Multiple Evaluations
```
Function:    f(x,y) = x + 2y
Expected ∇f: [1.0, 2.0] at all 5 points
Got ∇f:      [?, ~0] at all points
Error:       Y-component is ZERO everywhere
```

---

## Common Pattern: Y-Component Always Wrong

Every failing test has one thing in common: **the Y-component of the gradient is wrong**.

This strongly suggests a systematic issue in how the second component of the gradient is transformed from reference to physical coordinates.

---

## Next Steps

1. **Review** lines 160-172 of `src/Rodin/Variational/H1/Grad.h`
2. **Compare** with `src/Rodin/Variational/P1/Grad.h` (simpler, working implementation)
3. **Verify** that `J^{-T}` (Jacobian inverse transpose) is correct
4. **Check** if matrix-vector multiplication order is correct
5. **Add debug output** to trace intermediate values

---

## Files Referenced

### Test Files
- **Test Location**: `tests/unit/Rodin/Variational/H1/GradTest.cpp` (317 lines)
- **Test Functions**: 10 test cases for gradient operator

### Implementation Files
- **Main Bug Location**: `src/Rodin/Variational/H1/Grad.h` (290 lines, lines 160-172 problematic)
- **Reference Implementation**: `src/Rodin/Variational/P1/Grad.h`
- **Base Interface**: `src/Rodin/Variational/Grad.h`

---

## Build Information

### Configuration Used
```bash
cmake -DCMAKE_BUILD_TYPE=Debug -DRODIN_BUILD_UNIT_TESTS=ON ..
```

### Dependencies
- Boost 1.74
- Eigen 3
- CMake 3.16+
- C++20 compiler (GCC 10+)

### Build Output
- **Compilation**: Successful with non-critical warnings
- **Executable**: `build/tests/unit/Rodin/Variational/H1/RodinVariationalH1GradTest`
- **Build Time**: ~3 seconds
- **Test Execution Time**: 18 milliseconds

---

## Priority

🔴 **HIGH PRIORITY**

This bug affects core gradient computation functionality, which is essential for:
- PDE solvers using H1 finite elements
- Variational form assembly
- Finite element error analysis
- Any application using gradient-based computations

---

## Report Generation Date

Generated: 2024-01-19  
Test Execution: Complete with comprehensive analysis  
All Reports: Available in `/home/runner/work/rodin/rodin/`

---

## Questions?

Refer to the specific report files for:
- **What failed**: TEST_EXECUTION_REPORT.md
- **Why it failed**: H1_GRAD_TEST_SUMMARY.md
- **Exact numbers**: H1_GRAD_DETAILED_FAILURES.txt
- **Quick answers**: QUICK_REFERENCE.txt

