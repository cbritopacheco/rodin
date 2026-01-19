# H1 Grad Test Results Summary

## Build Status
✅ **BUILD SUCCESSFUL** - RodinVariationalH1GradTest built without errors

## Test Execution Summary
- **Total Tests**: 10
- **Passed**: 5 ✅
- **Failed**: 5 ❌
- **Total Time**: 18ms

## Passed Tests (5)
1. ✅ `Rodin_Variational_H1_Grad.ShapeFunction_Construction`
2. ✅ `Rodin_Variational_H1_Grad.GridFunction_Construction`
3. ✅ `Rodin_Variational_H1_Grad.GridFunction_ConstantFunction`
4. ✅ `Rodin_Variational_H1_Grad.Copy`
5. ✅ `Rodin_Variational_H1_Grad.UsageInBilinearForm`

## Failed Tests (5)

### 1. ❌ GridFunction_LinearFunction
**Location**: `GradTest.cpp:49-76`

**Test Description**: Tests that the gradient of a projected linear function f(x,y) = 3x + 4y is correctly computed as [3, 4]

**Failure Details**:
- Line 75 (grad_value(1) != 4.0):
  - Expected: `4.0`
  - Actual: `-2.000000000000024`
  - Difference: `6.000000000000024` (exceeds tolerance 1e-5)
  - Error Location: `GradTest.cpp:75`

**Root Cause**: The gradient component in the y-direction is being computed with the wrong sign, suggesting an issue with how gradients are computed or applied in the reference-to-physical coordinate transformation.

---

### 2. ❌ GridFunction_QuadraticFunction_H1_2
**Location**: `GradTest.cpp:144-175`

**Test Description**: Tests that the gradient of a quadratic function f(x,y) = x² + y² is correctly computed as [2x, 2y] using H1<2> (quadratic) space.

**Failure Details**:
**Point 1** (0.25, 0.25):
- Line 173 (grad_value(0)):
  - Expected: `2.0 * 1.5 = 3.0`
  - Actual: `13.500000000000107`
  - Difference: `10.500000000000107` (exceeds tolerance 1e-10)
  
- Line 174 (grad_value(1)):
  - Expected: `2.0 * 2.0 = 4.0`
  - Actual: `5.5000000000000497`
  - Difference: `1.5000000000000497` (exceeds tolerance 1e-10)

**Root Cause**: For quadratic functions, the gradient is being computed incorrectly. The error pattern suggests incorrect handling of the gradient basis functions or Jacobian inverse transformation for higher-order spaces.

---

### 3. ❌ GridFunction_QuadraticFunction_H1_3
**Location**: `GradTest.cpp:177-219`

**Test Description**: Tests that the gradient of f(x,y) = x² + 2xy + y² is correctly computed as [2x + 2y, 2x + 2y] using H1<3> (cubic) space.

**Failure Details**:
**Multiple evaluation points with significant errors**:

Point (0.25, 0.25):
- grad_value(0): Expected 1.0, Got 0.99074074074073482, Error: 0.0092592592592651846
- grad_value(1): Expected 1.0, Got 0.37345679012345617, Error: 0.62654320987654377

Point (0.5, 0.5):
- grad_value(0): Expected 2.0, Got 1.8333333333333566, Error: 0.16666666666664343
- grad_value(1): Expected 2.0, Got -1.8333333333333597, Error: 3.8333333333333597

Point (0.75, 0.25):
- grad_value(0): Expected 2.0, Got 2.0277777777777928, Error: 0.027777777777792778
- grad_value(1): Expected 2.0, Got -2.0277777777777932, Error: 4.0277777777777928

**Root Cause**: Similar to GridFunction_QuadraticFunction_H1_2, but the errors are more complex. The y-component consistently has wrong signs and magnitudes across different evaluation points, indicating systematic errors in gradient computation for higher-order spaces.

---

### 4. ❌ MultipleEvaluations
**Location**: `GradTest.cpp:221-257`

**Test Description**: Tests that the gradient of a linear function f(x,y) = x + 2y is correctly computed as [1, 2] at multiple evaluation points.

**Failure Details**:
Multiple evaluation points fail:

Point (0.0, 0.0):
- grad_value(0): Expected 1.0, Got 0.0, Error: 1.0
- grad_value(1): Expected 2.0, Got 4.4408920985006262e-16, Error: ~2.0

Point (1.0, 0.0):
- grad_value(1): Expected 2.0, Got 0.0, Error: 2.0

Points (0.0, 1.0), (0.5, 0.5), (0.25, 0.75):
- grad_value(1): Expected 2.0, Got values close to 0 or with floating-point noise

**Root Cause**: The gradient computation is yielding zero or near-zero values for the y-component across multiple points, suggesting that gradients in certain directions are not being computed at all, or the evaluation points are being handled incorrectly.

---

### 5. ❌ DifferentPolynomialDegrees
**Location**: `GradTest.cpp:259-315`

**Test Description**: Tests that the same linear function f(x,y) = 3x + 4y is correctly represented and its gradient [3, 4] is computed for different polynomial degrees (H1<1>, H1<2>, H1<3>).

**Failure Details**:
**H1<1> (linear basis)**:
- grad_value(0): Expected 3.0, Got 6.0, Error: 3.0
- grad_value(1): Expected 4.0, Got -2.0, Error: 6.0

**H1<2> (quadratic basis)**:
- grad_value(1): Expected 4.0, Got -2.000000000000024, Error: 6.000000000000024

**H1<3> (cubic basis)**:
- grad_value(1): Expected 4.0, Got -2.0000000000000946, Error: 6.0000000000000941

**Root Cause**: Across all polynomial degrees, the y-component of the gradient is consistently wrong (off by a factor of -3/2 or completely wrong), suggesting systematic issues with how the gradient direction is being handled or how reference-to-physical coordinate transformations are applied.

---

## Pattern Analysis

### Common Issues Observed:
1. **Sign Errors**: The y-component gradient frequently has the wrong sign (expected positive, got negative)
2. **Magnitude Errors**: Even when the sign is correct, magnitudes are significantly off
3. **Zero Gradients**: Some evaluation points return zero or near-zero gradients when they shouldn't
4. **Systematic Offset**: The errors follow patterns across different polynomial degrees and evaluation points

### Suspected Problem Areas (from code review):
Looking at `Grad.h` line 168-171:
```cpp
const auto& basis = fe.getBasis(local);
basis.getGradient()(rc);  // <- Line 168 result is not used!
s_res += gf[fes.getGlobalIndex({d, i}, local)] * basis.getGradient()(rc);
```

**Critical Bug Found**: Line 168 calls `basis.getGradient()(rc)` but doesn't use the result. This appears to be computing the gradient but not storing it. Then line 169 computes it again.

More importantly, the Jacobian transformation on line 171:
```cpp
out = p.getJacobianInverse().transpose() * s_res;
```

This might be incorrectly applying the Jacobian transformation. The reference gradient must be multiplied by the Jacobian inverse transpose: `J^{-T} ∘ ∇_ref u`.

---

## Test File Location
- Path: `/home/runner/work/rodin/rodin/tests/unit/Rodin/Variational/H1/GradTest.cpp`
- Implementation: `/home/runner/work/rodin/rodin/src/Rodin/Variational/H1/Grad.h`

## Next Steps
1. Review the Jacobian transformation logic in `Grad.h` lines 160-172
2. Verify that basis function gradients are being correctly computed in reference coordinates
3. Check if the physical coordinate transformation is correct
4. Validate the computation order and ensure gradients are stored before being used
5. Add debugging output to trace actual vs. expected gradient values

