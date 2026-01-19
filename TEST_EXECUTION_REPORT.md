# H1 Grad Test Execution Report

## Executive Summary

The RodinVariationalH1GradTest successfully **built** and **executed**, but **50% of tests (5 out of 10) are failing** with systematic gradient computation errors.

### Quick Facts
- **Build Status**: ✅ SUCCESS
- **Build Time**: ~3 seconds  
- **Total Tests**: 10
- **Passed**: 5 (50%) ✅
- **Failed**: 5 (50%) ❌
- **Execution Time**: 18ms

---

## Test Results Overview

### ✅ PASSING TESTS (5)

1. **ShapeFunction_Construction** - Tests that Grad can be constructed for shape functions
2. **GridFunction_Construction** - Tests that Grad can be constructed for grid functions
3. **GridFunction_ConstantFunction** - Gradient of constant function correctly returns zero
4. **Copy** - Grad objects can be copied correctly
5. **UsageInBilinearForm** - Grad can be used in bilinear forms and assembled into sparse matrices

**Note**: The fact that UsageInBilinearForm passes while evaluation tests fail suggests the gradient operator works in assembly context but fails during point evaluation.

### ❌ FAILING TESTS (5)

| Test Name | Function Tested | Issue |
|-----------|-----------------|-------|
| GridFunction_LinearFunction | f(x,y) = 3x + 4y | Y-gradient is -2.0 instead of 4.0 |
| GridFunction_QuadraticFunction_H1_2 | f(x,y) = x² + y² | Both components completely wrong in magnitude |
| GridFunction_QuadraticFunction_H1_3 | f(x,y) = x² + 2xy + y² | Y-component systematically negative across all test points |
| MultipleEvaluations | f(x,y) = x + 2y at 5 points | Y-component is zero everywhere |
| DifferentPolynomialDegrees | f(x,y) = 3x + 4y in H1<1>, H1<2>, H1<3> | Consistent y-component errors across all polynomial degrees |

---

## Critical Findings

### 1. **Systematic Y-Component Error**
Across all failing tests, the y-component of the gradient is consistently wrong:
- Sometimes has wrong sign (positive expected, negative computed)
- Sometimes is zero (expected non-zero value)
- Magnitudes don't match in any failing case

### 2. **Not a Polynomial Order Issue**
The same function fails for H1<1> (linear), H1<2> (quadratic), and H1<3> (cubic) basis spaces. This rules out basis function approximation errors.

### 3. **Not a Projection Issue**
The constant function test passes (returns zero gradient correctly), suggesting projection is working. The issue is in gradient computation itself.

### 4. **Location of Bug**
The bug is in the `interpolate()` method of the GridFunction specialization in `/src/Rodin/Variational/H1/Grad.h` (lines 98-172).

Specifically, lines 160-172:
```cpp
static thread_local SpatialVectorType s_res;

s_res.resize(d);
s_res.setZero();

assert(d == mesh.getDimension());
const auto& gf = this->getOperand();
const auto& fes = gf.getFiniteElementSpace();
const auto& fe = fes.getFiniteElement(d, i);
const auto& rc = p.getReferenceCoordinates();
for (size_t local = 0; local < fe.getCount(); local++)
{
  const auto& basis = fe.getBasis(local);
  basis.getGradient()(rc);  // <- LINE 168: Result not used!
  s_res += gf[fes.getGlobalIndex({d, i}, local)] * basis.getGradient()(rc);
}
out = p.getJacobianInverse().transpose() * s_res;
```

**Problem**: Line 171 applies Jacobian transformation but the reference gradient (`s_res`) is computed in line 169. The Jacobian inverse transpose multiplication on line 171 may be incorrect or applied to the wrong quantity.

---

## Error Patterns

### Pattern 1: Y-Component Sign Flip
```
GridFunction_LinearFunction:
  Expected: [3.0, 4.0]
  Got:      [?, -2.0]
```

### Pattern 2: Complete Magnitude Failure
```
GridFunction_QuadraticFunction_H1_2:
  Expected: [3.0, 4.0]
  Got:      [13.5, 5.5]
  Error:    [450%, 37.5%]
```

### Pattern 3: Zero Gradient
```
MultipleEvaluations at (0,0):
  Expected: [1.0, 2.0]
  Got:      [0.0, 4.44e-16]
```

### Pattern 4: Persistent Across All Spaces
```
DifferentPolynomialDegrees:
  H1<1>: grad[1] = -2.0 (expected 4.0)
  H1<2>: grad[1] = -2.0 (expected 4.0)
  H1<3>: grad[1] = -2.0 (expected 4.0)
```

---

## Test Coverage Analysis

### What's Being Tested
- ✅ Construction of Grad objects
- ✅ Copy semantics
- ✅ Usage in bilinear forms
- ❌ **Point evaluation of gradients**
- ❌ **Correctness of computed values**
- ❌ **Multiple evaluation points**
- ❌ **Different polynomial degrees**

### What's NOT Being Tested
- Gradient on boundaries/faces (there's code for this but no tests)
- Gradient assembly (passes, but might be working by accident)
- Integration with actual finite element assembly

---

## Detailed Error Messages

### Test 1: GridFunction_LinearFunction
```
File: GradTest.cpp:75
Assertion: EXPECT_NEAR(grad_value(1), 4.0, 1e-5)
Error: |(-2.000000000000024) - 4.0| = 6.000000000000024 > 1e-5
```

### Test 2: GridFunction_QuadraticFunction_H1_2  
```
File: GradTest.cpp:173-174
At point (1.5, 2.0):
  grad[0]: |13.5 - 3.0| = 10.5 > 1e-10 ✗
  grad[1]: |5.5 - 4.0| = 1.5 > 1e-10 ✗
```

### Test 3: GridFunction_QuadraticFunction_H1_3
```
File: GradTest.cpp:216-217
At multiple points:
  Point 1: grad[1] = 0.3735 (expected 1.0)
  Point 2: grad[1] = -1.8333 (expected 2.0)
  Point 3: grad[1] = -2.0278 (expected 2.0)
```

### Test 4: MultipleEvaluations
```
File: GradTest.cpp:254-255
At 5 different points:
  All have grad[1] ≈ 0 instead of 2.0
```

### Test 5: DifferentPolynomialDegrees
```
File: GradTest.cpp:280-281, 296-297, 312-313
Across H1<1>, H1<2>, H1<3>:
  grad[1] consistently = -2.0 instead of 4.0
```

---

## Build Details

### CMake Configuration
```bash
cmake -DCMAKE_BUILD_TYPE=Debug -DRODIN_BUILD_UNIT_TESTS=ON ..
```

### Dependencies Installed
- libboost1.74-all-dev
- libeigen3-dev

### Build Warnings (Non-fatal)
- Unused variables in H1.hpp
- Sign comparison warnings in OpenMP.h
- Unrecognized compiler flags (non-critical)

---

## Recommendations

### For Immediate Investigation
1. **Review Jacobian transformation logic** - Lines 171 in Grad.h likely incorrect
2. **Check basis.getGradient() computation** - Verify reference gradient is computed correctly
3. **Verify reference-to-physical mapping** - The J^{-T} multiplication might be transposed or applied incorrectly

### For Testing
1. Add debug output to trace intermediate gradient values
2. Compare with P1/Grad.h implementation (simpler linear case)
3. Add unit tests for Jacobian transformation separately

### For Debugging
Run individual tests with verbose output:
```bash
./build/tests/unit/Rodin/Variational/H1/RodinVariationalH1GradTest \
  --gtest_filter="Rodin_Variational_H1_Grad.GridFunction_LinearFunction"
```

---

## Files Analyzed

| File | Purpose |
|------|---------|
| `/tests/unit/Rodin/Variational/H1/GradTest.cpp` | Test cases (317 lines) |
| `/src/Rodin/Variational/H1/Grad.h` | Implementation (290 lines) |
| `/src/Rodin/Variational/Grad.h` | Base interface |
| `/src/Rodin/Variational/P1/Grad.h` | Reference implementation (simpler) |

---

## Summary

The H1 Grad operator has a **systematic bug in gradient evaluation** that affects all GridFunction point evaluations. The bug appears to be in how the Jacobian inverse transformation is applied in lines 160-172 of `H1/Grad.h`.

The issue is:
- **Not** in construction or copying
- **Not** in basis function approximation
- **Not** in bilinear form assembly
- **Is** in the interpolate method's physical coordinate transformation

**Priority**: HIGH - This affects core functionality of gradient evaluation, which is essential for many PDE solvers.

