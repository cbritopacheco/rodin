/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file UMFPack.h
 * @brief UMFPACK multifrontal sparse LU solver wrapper.
 *
 * This header provides a wrapper for the UMFPACK solver from SuiteSparse,
 * a high-performance direct solver for general sparse matrices using
 * multifrontal LU factorization.
 *
 * ## Algorithm
 * UMFPACK computes:
 * @f[
 *   PAQ = LU
 * @f]
 * where @f$ P @f$ and @f$ Q @f$ are permutation matrices, @f$ L @f$ is lower
 * triangular, and @f$ U @f$ is upper triangular. The multifrontal method
 * processes the matrix using an assembly tree for efficiency.
 *
 * ## Applicability
 * - General sparse matrices (non-symmetric, non-positive definite)
 * - Large-scale problems
 * - Circuit simulation
 * - Structural analysis
 * - Computational fluid dynamics
 *
 * ## Usage Example
 * ```cpp
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
 * 
 * Solver::UMFPack solver(problem);
 * solver.solve();
 * ```
 *
 * @note This solver requires UMFPACK from SuiteSparse to be installed
 * and RODIN_USE_UMFPACK to be defined at compile time.
 *
 * @see UMFPack for the solver implementation
 */
#ifndef RODIN_SOLVER_UMFPACK_H
#define RODIN_SOLVER_UMFPACK_H

#ifdef RODIN_USE_UMFPACK

#include <optional>
#include <functional>

#include <Eigen/UmfPackSupport>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"

namespace Rodin::FormLanguage
{
  template <class LinearSystem>
  struct Traits<Solver::UMFPack<LinearSystem>>
  {
    using LinearSystemType = LinearSystem;
  };
}

namespace Rodin::Solver
{
  /**
   * @defgroup UMFPackSpecializations UMFPack Template Specializations
   * @brief Template specializations of the UMFPack class.
   * @see UMFPack
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for UMFPack
   */
  template <class LinearSystem>
  UMFPack(Variational::ProblemBase<LinearSystem>&) -> UMFPack<LinearSystem>;

  /**
   * @ingroup UMFPackSpecializations
   * @brief UMFPACK multifrontal sparse LU solver for general sparse systems.
   *
   * This solver provides access to UMFPACK, a high-performance direct solver
   * from SuiteSparse that uses multifrontal LU factorization. UMFPACK is
   * particularly effective for general sparse matrices.
   *
   * **Characteristics:**
   * - **Type**: Direct solver
   * - **Matrix requirements**: General sparse (no symmetry required)
   * - **Memory**: Efficient for sparse systems
   * - **Performance**: Excellent for large sparse systems
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class UMFPack<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
    : public LinearSolverBase<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      /// Type of scalar values in the system
      using ScalarType = Scalar;

      /// Type of solution and right-hand side vectors
      using VectorType = Math::Vector<ScalarType>;

      /// Type of system matrix (operator)
      using OperatorType = Math::SparseMatrix<ScalarType>;

      /// Type of linear system being solved
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      /// Type of problem being solved
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      /// Parent class type
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      /**
       * @brief Constructs a UMFPack solver for the given problem.
       * @param pb Reference to the problem to solve
       */
      UMFPack(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /**
       * @brief Copy constructor.
       * @param other Solver to copy from
       */
      UMFPack(const UMFPack& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Solver to move from
       */
      UMFPack(UMFPack&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Default destructor.
       */
      ~UMFPack() = default;

      /**
       * @brief Solves the linear system using UMFPACK.
       * @param[in,out] axb The linear system to solve. The solution is stored
       *                    in the system's solution vector.
       */
      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getVector());
      }

      /**
       * @brief Prints UMFPACK control parameters.
       *
       * Displays the current control settings used by UMFPACK.
       */
      void printControl()
      {
        m_solver.printUmfpackControl();
      }

      /**
       * @brief Prints UMFPACK information.
       *
       * Displays statistics about the factorization.
       */
      void printInfo()
      {
        m_solver.printUmfpackInfo();
      }

      /**
       * @brief Prints UMFPACK status.
       *
       * Displays the current status of the solver.
       */
      void printStatus()
      {
        m_solver.printUmfpackStatus();
      }

      /**
       * @brief Creates a copy of this solver.
       * @returns Pointer to a new UMFPack instance
       */
      UMFPack* copy() const noexcept override
      {
        return new UMFPack(*this);
      }

    private:
      /// Underlying Eigen UMFPACK solver
      Eigen::UmfPackLU<OperatorType> m_solver;
  };
}

#endif // #ifdef RODIN_USE_UMFPACK
#endif


