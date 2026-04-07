# Rodin Math module

This directory contains core numerical types and low-level mathematical helpers used across Rodin.

## Public surface (stable headers)

- `Rodin/Math.h` (umbrella include)
- `Rodin/Math/Constants.h`
- `Rodin/Math/Common.h`
- `Rodin/Math/Vector.h`
- `Rodin/Math/Matrix.h`
- `Rodin/Math/SparseMatrix.h`
- `Rodin/Math/SpatialVector.h`
- `Rodin/Math/SpatialMatrix.h`
- `Rodin/Math/LinearSystem.h`
- `Rodin/Math/Unit.h`
- `Rodin/Math/Rad.h`
- `Rodin/Math/RungeKutta/RK2.h`
- `Rodin/Math/RungeKutta/RK4.h`
- `Rodin/Math/RootFinding/NewtonRaphson.h`

## Internal/support headers

- `ForwardDecls.h`, `Traits.h`, `Types.h`

These are primarily support headers for type traits and forward declarations.

## Notes

- Dense and sparse algebra is based on Eigen.
- `SpatialVector`/`SpatialMatrix` provide bounded-size math objects (up to `RODIN_MAXIMAL_SPACE_DIMENSION`).
- `LinearSystem` provides dense and sparse system containers with elimination/merge utilities.
