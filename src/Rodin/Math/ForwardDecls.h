/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ForwardDecls.h
 * @brief Forward declarations for the Math module.
 *
 * This file provides forward declarations of key types in the Math module,
 * allowing other headers to reference them without full definitions.
 */
#ifndef RODIN_MATH_FORWARDDECLS_H
#define RODIN_MATH_FORWARDDECLS_H

namespace Rodin::Math
{
  /**
   * @brief Represents a linear system of equations.
   *
   * Forward declaration for the LinearSystem class template, which represents
   * systems of linear equations in the form:
   * @f[
   *   \text{Operator} \cdot \text{Solution} = \text{Vector}
   * @f]
   *
   * For example, @f$ Ax = b @f$ where:
   * - @f$ A @f$ is the operator (matrix)
   * - @f$ x @f$ is the solution vector (unknown)
   * - @f$ b @f$ is the right-hand side vector
   *
   * @tparam Operator Type of the linear operator (e.g., SparseMatrix, Matrix)
   * @tparam Vector Type of the vectors (solution and right-hand side)
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref LinearSystem "LinearSystem<SparseMatrix, Vector>" | Local sparse linear system used by finite element assembly and iterative/direct sparse solvers. |
   * | @ref LinearSystem "LinearSystem<Matrix, Vector>" | Local dense linear system used by dense solvers and small algebraic systems. |
   * | @ref LinearSystem "LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>" | PETSc-backed distributed linear system assembled and solved through PETSc objects. |
   * | @ref LinearSystem "LinearSystem<Operator, Vector>" | Generic linear system interface for compatible operator/vector pairs. |
   *
   * @see <a href="class_rodin_1_1_math_1_1_linear_system.html">LinearSystem</a>
   */
  template <class Operator, class Vector>
  class LinearSystem;

  template <class Scalar>
  class SpatialVector;

  template <class Scalar>
  class SpatialMatrix;
}

#endif
